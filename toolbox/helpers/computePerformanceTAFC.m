function [predictions, theClassifierEngine, responses] = computePerformanceTAFC(nullScene, testScene, ...
    temporalSupport, nTrain, nTest, theNeuralEngine, theClassifierEngine, trainNoiseFlag, testNoiseFlag, ...
    saveResponses, visualizeAllComponents)
% Compute performance of a classifier given a null and test scene, a neural engine, and a classifier engine.
%
% Syntax:
%    [predictions, theClassifierEngine, responses] = ...
%        computePerformanceTAFC(nullScene, testScene, temporalSupport, nTrain, nTest, theNeuralEngine, ...
%       theClassifierEngine, trainNoiseFlag, testNoiseFlag, saveResponses, visualizeAllComponents)
%
% Description:
%     Train a classifier on a discrimination and report back a vector of
%     1's and 0's indicating correct and incorrect trials respectively.
%
%     This uses the ISETBioCSFGeneratorFramework and works because the uers
%     passes a set of objects with standardized API.  These describe the
%     two scenes to be discriminated, the neural pipeline that processes
%     these scenes, and the classifer.
%
% Inputs:
%     nullScene             - Null scene sequence.
%     testScene             - Test scene sequence.
%     temporalSupport       - Temporal support vector (in seconds) for
%                             scene sequences.
%     nTrain                - Number of null and test response instances
%                             used in classifer training.  The two types of
%                             instances are paired and a nTrain TAFC task is
%                             simulated.
%     nTest                 - Number of null and test response instances
%                             used in classifer training.  The two types of
%                             instances are paired and nTest TAFC
%                             trials are simulated for evaluating
%                             performance.
%     theNeuralEngine       - @neuralResponseEngine object to compute
%                             neural responses.
%     theClassifierEngine   - @responseClassifierEngine object that
%                             implements observer decision model.  This is
%                             assumed untrained if trainedNoiseFlag
%                             contains a string, and trained if
%                             trainedNoiseFlag is empty.
%     trainNoiseFlag        - String.  Type of noise to be used in training
%                             the classifier. This flag are passed to
%                             theNeuralEngine to generate the training
%                             response instances. Typically either 'none' or
%                             'random' depending on whether the desired
%                             classifier is signal known exactly ('none')
%                             or signal known statistically ('random').
%     testNoiseFlag          - String. Type of noise to be used in
%                             evaluating performance. This flag are passed to
%                             theNeuralEngine to generate the test
%                             response instances. Typically 'random'.
%     saveResponses         - Logical. Whether to return the computed
%                             response instances
%     visualAllComponentrs  - Logical. Whether to visualize or not.
%
% Outputs:
%     predictions            - Vector of 1's (correct) and 0's (incorrect)
%                              that gives trial by trial performance of the
%                              tested classifier in the TAFC task.
%                              Contains nTest entries.
%     theClassifierEngine    - Trained version of passed classifier object.
%
%     responses              - Neural responses computed
%
% Optional key/value pairs:
%     None.
%
% See also
%   t_thresholdEngine, t_spatialCsf, computeThresholdTAFC
%

% History:
%   10/23/20  dhb  Comments.

% Empty responses
responses = [];

% Train the classifier.
%
% If trainNoiseFlag is empty, then the passed classifier has already been trained
% and training is skipped.  Otherwise trainFlag is passed to the stimulus
% generation routine to indicate what type of noise (typically 'none' or
% 'random') should be used in the training.

if (~isempty(trainNoiseFlag))
    % Generate stimulus for training, NULL stimulus
    [inSampleNullStimResponses, ~] = theNeuralEngine.compute(...
        nullScene, ...
        temporalSupport, ...
        nTrain, ...
        'noiseFlags', {trainNoiseFlag});
    
    % Generate stimulus for training, TEST stimulus
    [inSampleTestStimResponses, responseTemporalSupportSeconds] = theNeuralEngine.compute(...
        testScene, ...
        temporalSupport, ...
        nTrain, ...
        'noiseFlags', {trainNoiseFlag});
    
    if (visualizeAllComponents)
        if (isfield(theNeuralEngine.neuralPipeline, 'coneMosaic'))
            diffResponse = inSampleTestStimResponses(trainNoiseFlag) - inSampleNullStimResponses(trainNoiseFlag);
            % Visualize the activation
            theNeuralEngine.neuralPipeline.coneMosaic.visualize('activation', squeeze(diffResponse), 'verticalActivationColorBarInside', true);
        
            % Also visualize the full absorptions density
            figNo = 999;
            theNeuralEngine.neuralPipeline.coneMosaic.visualizeFullAbsorptionsDensity(figNo);
        end
        
        if (isfield(theNeuralEngine.neuralPipeline, 'mRGCmosaic'))
            theNeuralEngine.neuralPipeline.mRGCmosaic.visualizeResponses(...
                responseTemporalSupportSeconds, inSampleTestStimResponses(trainNoiseFlag), ...
                'stimulusTemporalSupportSeconds', temporalSupport,...
                'stimulusSceneSequence', testScene);
        end
        
    end
    
    % Train the classifier. This shows the usage to extact information
    % from the container retrned as the first return value from the neural
    % response engine - we index the responses by the string contained in
    % the variable trainFlag (which was itself passed to the neural
    % repsonse engine above.)
    %
    % Once extracted from the container, the responses are a 3 dimensional
    % matrix, with the dimensions indexing [instancesNum x mNeuralDim x tTimeBins].
    %   instancesNum   - number of response instances
    %   mNeuralDim     - dimension of neural response at one timepoint
    %   tTimeBins      - number of time points in stimulus sequence.
    theClassifierEngine.compute('train', ...
        inSampleNullStimResponses(trainNoiseFlag), ...
        inSampleTestStimResponses(trainNoiseFlag));
    
    % Save computed response instances
    if (saveResponses)
        responses.inSampleNullStimResponses = inSampleNullStimResponses;
        responses.inSampleTestStimResponses = inSampleTestStimResponses;
    end
    
end

% Predict using trained classifier.
%
% Generate stimulus for prediction, NULL stimulus.  The variable testFlag
% indicates what type of noise is used to generate the stimuli used for
% prediction.  Typically 'random'.
[outOfSampleNullStimResponses, ~] = theNeuralEngine.compute(...
    nullScene, ...
    temporalSupport, ...
    nTest, ...
    'noiseFlags', {testNoiseFlag});

% Generate stimuli for prediction, TEST stimulus
[outOfSampleTestStimResponses, ~] = theNeuralEngine.compute(...
    testScene, ...
    temporalSupport, ...
    nTest, ...
    'noiseFlags', {testNoiseFlag});

% Do the prediction
dataOut = theClassifierEngine.compute('predict', ...
    outOfSampleNullStimResponses(testNoiseFlag), ...
    outOfSampleTestStimResponses(testNoiseFlag));

% Save computed response instances
if (saveResponses)
    responses.outOfSampleNullStimResponses = outOfSampleNullStimResponses;
    responses.outOfSampleTestStimResponses = outOfSampleTestStimResponses;
end
    
% Set return variable.  For each trial 0 means wrong and 1 means right.
% Taking mean(response) gives fraction correct.
predictions = dataOut.trialPredictions;

end