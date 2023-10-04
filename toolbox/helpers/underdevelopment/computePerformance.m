function [predictions, theClassifierEngine, responses] = computePerformance(theScenes, ...
    temporalSupport, nTrain, nTest, theNeuralEngine, theClassifierEngine, trainNoiseFlag, testNoiseFlag, ...
    varargin)
% Compute performance of a classifier given different scenes, a neural engine, and a classifier engine.
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
%
%
% Optional key/value pairs:
%     TAFC                  - logical. Whether this is a two-interval 
%                             forced-choice task or N-way one-stimulus-per
%                             -trial task. Default false.
%     saveResponses         - Logical. Whether to return the computed
%                             response instances
%     visualAllComponents   - Logical. Whether to visualize or not.
%
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
%   t_thresholdEngine, t_spatialCsf, computeThreshold
%

% History:
%   10/23/20  dhb  Comments.
%   04/10/24  fh   Merged computePerformanceTAFC.m and
%                   computePerformanceNWay_OneStimulusPerTrial.m by adding 
%                   a key/pair pair specifying whether the task is TAFC or 
%                   NWay_OneStimulusPerTrial.

p = inputParser;
p.addParameter('TAFC',  false, @islogical);
p.addParameter('saveResponses',false, @islogical);
p.addParameter('visualizeAllComponents', false, @islogical);

parse(p, varargin{:});
isTAFC = p.Results.TAFC;
saveResponses = p.Results.saveResponses;
visualizeAllComponents = p.Results.visualizeAllComponents;

% Empty responses
responses = [];
nScenes   = length(theScenes);
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

    %If the task is TAFC, then we need to do the following reorganization
    %of the data
    if isTAFC
        %concatenate them along the 3rd dimension (cones)
        cat1 = cat(3, inSampleStimResponsesCell{1}, inSampleStimResponsesCell{2});
        cat2 = cat(3, inSampleStimResponsesCell{2}, inSampleStimResponsesCell{1});
        %nicely put the concatenated cells back to the container
        inSampleStimResponses = containers.Map(trainNoiseFlag, {cat1, cat2});
    else %NWay_OneStimulusPerTrial
        inSampleStimResponses = combineContainers(inSampleStimResponsesCell);
    end

    % Visualization the cone excitation for selected stimulus
    if visualizeAllComponents
        %TAFC: 1st stim is the test; NWay: 1st stim is the correct stim
        theStim = 1; 
        visualizeConeResps(theNeuralEngine, inSampleStimResponses,...
            trainNoiseFlag, theStim);
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
if isTAFC
    nTests_eachScene = nTest;
    whichAlternatives = ones(nTests_eachScene, 1);
else
    % Alternative implementation here ...
    nTests_eachScene = nTest/nScenes;
    assert(mod(nTest, nScenes) == 0, ['The number of test trials must be an',...
        ' integer multiple of the number of alternative choices']);
    whichAlternatives = repmat(1:nScenes,[nTests_eachScene, 1]);
    whichAlternatives = whichAlternatives(:);
end

outSampleStimResponsesCell = cell(1, nScenes);
eStart = tic;
for n = 1:nScenes
    [outSampleStimResponsesCell{n}, ~] = theNeuralEngine.compute(...
        theScenes{n}, ...
        temporalSupport, ...
        nTests_eachScene, ...
        'noiseFlags', {testNoiseFlag});
end
e = toc(eStart);
fprintf('computePerformance: Took %0.1f secs to generate mean responses for all alternatives\n',e);

%If the task is TAFC, then we need to do the following reorganization
%of the data
if isTAFC
    %get the responses given the test stimulus
    outSampleTestStimResponses = outSampleStimResponsesCell{1}(testNoiseFlag);
    %combine the 2nd (time) and the 3rd dimensions (cones)
    outSampleTestStimResponses = outSampleTestStimResponses(:,:);
    %get the responses given the null stimulus
    outSampleNullStimResponses = outSampleStimResponsesCell{2}(testNoiseFlag);
    outSampleNullStimResponses = outSampleNullStimResponses(:,:);
    %concatenate them together along the 2nd dimension (cones)
    cat_resp = cat(2, outSampleTestStimResponses, outSampleNullStimResponses);
    %nicely put them back to the container
    outSampleStimResponses = containers.Map(testNoiseFlag, cat_resp);
else
    outSampleStimResponses = combineContainersMat(outSampleStimResponsesCell);
end

% Do the prediction
dataOut = theClassifierEngine.compute('predict', ...
    outSampleStimResponses(testNoiseFlag),whichAlternatives);

% Save computed response instances
if (saveResponses)
    responses.outSampleStimResponses = outSampleStimResponses;
end
    
% Set return variable.  For each trial 0 means wrong and 1 means right.
% Taking mean(response) gives fraction correct.
predictions = dataOut.trialPredictions;

end

function visualizeConeResps(theNeuralEngine, inSampleTestStimResponses, ...
    trainNoiseFlag, selectedScene)
if (isfield(theNeuralEngine.neuralPipeline, 'coneMosaic'))
    Responses = inSampleTestStimResponses(trainNoiseFlag);
    Responses_slc = Responses{selectedScene};
    % Visualize the activation
    theNeuralEngine.neuralPipeline.coneMosaic.visualize('activation', ...
        squeeze(Responses_slc), 'verticalActivationColorBarInside', true);

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