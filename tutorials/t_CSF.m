% Build on the threshold engine tutorial (t_thresholdEngine), here we compute
%  the threshold a static gabor stimulus of  multiple spatial frequencies
% to get the Contrast Sensitivity Function (CSF) for a Poisson 2AFC ideal observer.



% luminance stimulus
computeThreshold([0.5, 0.5, 0.5],  5);

% Compute threshold for a particular chromatic direction and spatial frequency 
% Chromatic direction is specified 
function [threshold] = computeThreshold(chromaDir, spatialFreq)

% Compute function handle for grating stimuli
sceneComputeFunction = @sceGrating;

% Retrieve the default params for the grating stimulus
gratingParams = sceGrating();

% Configure chromatic direction and and spatial frequency of the grating
% with a 90 deg orientation, and a cosine spatial phase
gratingParams.coneContrastModulation = chromaDir;
gratingParams.spatialFrequencyCyclesPerDeg = spatialFreq;
gratingParams.spatialPhaseDegs = 0;
gratingParams.orientationDegs = 90;

% Configure a disk spatial envelope
gratingParams.spatialEnvelope = 'disk';
gratingParams.minPixelsNumPerCycle = 30;
gratingParams.spatialEnvelopeRadiusDegs = 0.4;

% Configure temporal modulation: 100 ms duration for only 1 frame
gratingParams.frameDurationSeconds = 100/1000;
gratingParams.temporalModulation = 'flashed';
gratingParams.temporalModulationParams =  struct(...
    'stimOnFrameIndices', 1, 'stimDurationFramesNum', 1);

% Instantiate a sceneEngine with the above sceneComputeFunctionHandle
% and the custom grating params.
theSceneEngine = sceneEngine(sceneComputeFunction, gratingParams);

% Compute the scene sequence
visualizationContrast = 1.0;
[theSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(visualizationContrast);

% Visualize the generated scene sequence
theSceneEngine.visualizeSceneSequence(theSceneSequence, theSceneTemporalSupportSeconds);

end

function [predictions, theClassifierEngine] = computePerformance(nullScene, testScene, temporalSupport, nTrain, nTest, theNeuralEngine, theClassifierEngine, trainFlag, testFlag)

% Train the classifier.
%
% If trainFlag is empty, then the passed classifier has already been trained
% and training is skipped.  Otherwise trainFlag is passed to the stimulus
% generation routine to indicate what type of noise (typically 'none' or
% 'random') should be used in the training.
if (~isempty(trainFlag))
    % Generate stimulus for training, NULL stimulus
    [inSampleNullStimResponses, ~] = theNeuralEngine.compute(...
        nullScene, ...
        temporalSupport, ...
        nTrain, ...
        'noiseFlags', {trainFlag});
    
    % Generate stimulus for training, TEST stimulus
    [inSampleTestStimResponses, ~] = theNeuralEngine.compute(...
        testScene, ...
        temporalSupport, ...
        nTrain, ...
        'noiseFlags', {trainFlag});
    
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
        inSampleNullStimResponses(trainFlag), ...
        inSampleTestStimResponses(trainFlag));
end

% Predict using trained classifier.
%
% Generate stimulus for prediction, NULL stimulus.  The variable testFlag
% indicates what type of noise is used to generate the stimuli used for
% prediction.  Typically 'random'.
[inSampleNullStimResponses, ~] = theNeuralEngine.compute(...
    nullScene, ...
    temporalSupport, ...
    nTest, ...
    'noiseFlags', {testFlag});

% Generate stimuli for prediction, TEST stimulus
[inSampleTestStimResponses, ~] = theNeuralEngine.compute(...
    testScene, ...
    temporalSupport, ...
    nTest, ...
    'noiseFlags', {testFlag});

% Do the prediction
dataOut = theClassifierEngine.compute('predict', ...
    inSampleNullStimResponses(testFlag), ...
    inSampleTestStimResponses(testFlag));

% Set return variable.  For each trial 0 means wrong and 1 means right.
% Taking mean(response) gives fraction correct.
predictions = dataOut.trialPredictions;

end