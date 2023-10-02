function [predictions, theClassifierEngine, responses] = computePerformance(task, theScenes, ...
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
%     visualAllComponents   - Logical. Whether to visualize or not.
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
nScenes = length(theScenes);
%nScenes = 2 if the input task is 'TAFC'
%nScenes = #alternatives if the input task is 'NWay_OneStimulusPerTrial'

% Train the classifier.
%
% If trainNoiseFlag is empty, then the passed classifier has already been trained
% and training is skipped.  Otherwise trainFlag is passed to the stimulus
% generation routine to indicate what type of noise (typically 'none' or
% 'random') should be used in the training.

if (~isempty(trainNoiseFlag))
    inSampleStimResponsesCell = cell(1,nScenes);
    % Generate stimuli for training
    for n = 1:nScenes
        [inSampleStimResponsesCell{n}, ~] = theNeuralEngine.compute(...
            theScenes{n}, ...
            temporalSupport, ...
            nTrain, ...
            'noiseFlags', {trainNoiseFlag});
    end
    inSampleStimResponses = combineContainers(inSampleStimResponsesCell);

    % Visualization.
    if visualizeAllComponents
        visualizeConeResps(theNeuralEngine, inSampleStimResponses,...
            trainNoiseFlag, nScenes);
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
    theClassifierEngine.compute('train', inSampleStimResponses(trainNoiseFlag),[]);

    % Save computed response instances
    if (saveResponses)
        responses.inSampleStimResponses = inSampleStimResponses;
    end
end

% Predict using trained classifier.
%
% Generate stimulus for prediction, NULL stimulus.  The variable testFlag
% indicates what type of noise is used to generate the stimuli used for
% prediction.  Typically 'random'.
switch task
    case 'TAFC'
        nTests_eachScene = nTest;
    case 'NWay_OneStimulusPerTrial'
        % Alternative implementation here ...
        nTests_eachScene = nTest/nScenes;
        assert(mod(nTest, nScenes) == 0, 'The number of test trials must be an integer multiple of the number of alternative choices');
        whichAlternatives = repmat(1:nScenes,[nTests_eachScene, 1]);
        whichAlternatives = whichAlternatives(:);
    otherwise 
        %technically this would not happen because if the following error is true, 
        % computeThreshold.m would throw an error already
        error('Invalid task name! You can either input ''TAFC'' or ''NWay_OneStimulusPerTrial''');
end

outOfSampleStimResponsesCell = cell(1, nScenes);
eStart = tic;
for n = 1:nScenes
    [outOfSampleStimResponsesCell{n}, ~] = theNeuralEngine.compute(...
        theScenes{n}, ...
        temporalSupport, ...
        nTests_eachScene, ...
        'noiseFlags', {testNoiseFlag});
end
e = toc(eStart);
fprintf('computePerformance: Took %0.1f secs to generate mean responses for all alternatives\n',e);
outOfSampleStimResponses = combineContainersMat(outOfSampleStimResponsesCell);

switch task
    case 'TAFC'
        % Do the prediction
        dataOut = theClassifierEngine.compute('predict', ...
            outOfSampleStimResponses(testNoiseFlag),[]);
    case 'NWay_OneStimulusPerTrial'
        dataOut = theClassifierEngine.compute('predict', ...
            outOfSampleStimResponses(testNoiseFlag),whichAlternatives);
end

% Save computed response instances
if (saveResponses)
    responses.outOfSampleStimResponses = outOfSampleStimResponses;
end
    
% Set return variable.  For each trial 0 means wrong and 1 means right.
% Taking mean(response) gives fraction correct.
predictions = dataOut.trialPredictions;

end

function visualizeConeResps(theNeuralEngine, inSampleTestStimResponses, ...
    inSampleNullStimResponses, trainNoiseFlag)
if (isfield(theNeuralEngine.neuralPipeline, 'coneMosaic'))
    diffResponse = inSampleTestStimResponses(trainNoiseFlag) - ...
        inSampleNullStimResponses(trainNoiseFlag);
    % Visualize the activation
    theNeuralEngine.neuralPipeline.coneMosaic.visualize('activation', ...
        squeeze(diffResponse), 'verticalActivationColorBarInside', true);

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