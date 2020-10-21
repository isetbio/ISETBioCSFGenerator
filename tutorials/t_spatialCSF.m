% Compute spatial CSF in different color directions
%
% Description:
%    Use ISETBioCSFGenerator to run out CSFs in different color directions.
%    This example uses an ideal Poisson TAFC observer and circularly
%    windowed gratings of constant size.
%
% See also: t_thresholdEngine, t_modulatedGratingsSceneGeneration
%

% History:
%   10/20/20  lz   Wrote it.
%   10/21/20  dhb  More commments.

%% CSF Calculation
%
% List of spatial frequencies to be tested.
spatialFreqs = [0.5, 1, 2, 4, 8, 12, 16, 25];

% Allocate space for thresholds
logThreshold = zeros(1, length(spatialFreqs));

% Choose stimulus chromatic direction specified as a 1-by-3 vector
% of L, M, S cone contrast.  These vectors get normalized below, so only
% their direction matters in the specification.
stimType = 'luminance';
switch (stimType)
    case 'luminance'
        chromaDir = [1.0, 1.0, 1.0];        
    case 'red-green'
        chromaDir = [1.0, -1.0, 0.0];
    case 'L-isolating'
         chromaDir = [1.0, 0.0, 0.0];
end

% Set the RMS cone contrast of the stimulus. Things may go badly if you
% exceed the gamut of the monitor, so we are conservative and set this at a
% value that is within gamut of typical monitors and don't worry about it
% further for this tutorial.  A vector length contrast of 0.08 should be
% OK.
rmsContrast = 0.08;
chromaDir = chromaDir / norm(chromaDir) * rmsContrast;
assert(abs(norm(chromaDir) - rmsContrast) <= 1e-10);

%% Compute threshold for each spatial frequency
%
% Function computeThreshold is below.
thePsychometricFcnFig = figure;
for idx = 1:length(spatialFreqs)
    logThreshold(idx) = computeThreshold(chromaDir,  spatialFreqs(idx), idx, thePsychometricFcnFig);
end
set(gcf, 'Position',  [0, 0, 800, 800]);

% Convert returned log threshold to linear threshold
threshold = 10 .^ logThreshold;

%% Plot Contrast Sensitivity Function
theCsfFig = figure;
loglog(spatialFreqs, 1 ./ threshold, '-ok', 'LineWidth', 2);
xlabel('Spatial Frequency (cyc/deg)');
ylabel('Sensitivity');
set(gcf, 'Position',  [0, 0, 600, 800]);

%% Helper functions for calculating threshold and classifier performance

% Compute threshold for a particular chromatic direction and spatial frequency
% Chromatic direction is a 1-by-3 vector specifying contrast on the L, M and S Cone, respectively
% Spatial frequency is a number in the unit of cycles per degree
function [logThreshold] = computeThreshold(chromaDir, spatialFreqs, index, theFig)

% Compute function handle for grating stimuli
sceneComputeFunction = @sceGrating;

% Retrieve the default params for the grating stimulus
gratingParams = sceGrating();

% Configure chromatic direction and and spatial frequency of the grating
% with a 90 deg orientation, and a cosine spatial phase
gratingParams.coneContrastModulation = chromaDir;
gratingParams.spatialFrequencyCyclesPerDeg = spatialFreqs;
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

% Instantiate a neuralResponseEngine
neuralParams = nrePhotopigmentExcitationsWithNoEyeMovements;
neuralParams.coneMosaicParams.fovDegs = 0.25;
theNeuralEngine = neuralResponseEngine(@nrePhotopigmentExcitationsWithNoEyeMovements, neuralParams);

% The actual threshold varies enough with the different engines that we
% need to adjust the contrast range that Quest+ searches over, as well as
% the range of psychometric function slopes.
logThreshLimitLow = 2.5; logThreshLimitHigh = 0; logThreshLimitDelta = 0.02;
slopeRangeLow = 1; slopeRangeHigh = 100; slopeDelta = 2.5;

% Instantiate the PoissonTAFC responseClassifierEngine
% PoissonTAFC makes decision by performing the Poisson likelihood ratio test
classifierEngine = responseClassifierEngine(@rcePoissonTAFC);
trainFlag = 'none'; testFlag = 'random';
nTrain = 1;  nTest = 128;

% Construct a QUEST threshold estimator estimate threshold on log contrast
% Run a fixed number of trials (i.e., 10 contrast level, 1280 trials in total)
estDomain  = -logThreshLimitLow : logThreshLimitDelta : -logThreshLimitHigh;
slopeRange = slopeRangeLow: slopeDelta : slopeRangeHigh;

estimator = questThresholdEngine('minTrial', 1280, 'maxTrial', 1280, ...
    'estDomain', estDomain, 'slopeRange', slopeRange, 'numEstimator', 1);

% Generate the NULL stimulus (zero contrast)
nullContrast = 0.0;
[theNullSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(nullContrast);

% Threshold estimation with QUEST+
% Get the initial stimulus contrast from QUEST+
[logContrast, nextFlag] = estimator.nextStimulus();

% Loop over trials.
testedContrasts = [];
while (nextFlag)
    
    % Convert log contrast -> contrast
    testContrast = 10 ^ logContrast;
    
    % Have we already built the classifier for this contrast?
    testedIndex = find(testContrast == testedContrasts);
    if (isempty(testedIndex))
        % No.  Save contrast in list
        testedContrasts = [testedContrasts testContrast];
        testedIndex = find(testContrast == testedContrasts);
        
        % Generate the TEST scene sequence for the given contrast
        [theTestSceneSequences{testedIndex}, ~] = theSceneEngine.compute(testContrast);
        
        % Train classifier for this TEST contrast and get predicted
        % correct/incorrect predictions.  This function also computes the
        % neural responses needed to train and predict.
        [predictions, theTrainedClassifierEngines{testedIndex}] = computePerformance(...
            theNullSceneSequence, theTestSceneSequences{testedIndex}, ...
            theSceneTemporalSupportSeconds, nTrain, nTest, ...
            theNeuralEngine, classifierEngine, trainFlag, testFlag);
        
    else
        % Classifier is already trained, just get predictions
        predictions = computePerformance(...
            theNullSceneSequence, theTestSceneSequences{testedIndex}, ...
            theSceneTemporalSupportSeconds, nTrain, nTest, ...
            theNeuralEngine, theTrainedClassifierEngines{testedIndex}, [], testFlag);
    end
    
    % Tell QUEST+ what we ran (how many trials at the given contrast) and
    % get next stimulus contrast to run.
    [logContrast, nextFlag] = ...
        estimator.multiTrial(logContrast * ones(1, nTest), predictions);
    
    % Report what happened
    % fprintf('Current test contrast: %g, P-correct: %g \n', testContrast, mean(predictions));    
end

% Add info to correct panel of the figure
figure(theFig);
subplot(4, 4, index * 2 - 1);

% Compute the scene sequence
% Visualize the generated scene sequence
visualizationContrast = 1.0;
[theSceneSequence] = theSceneEngine.compute(visualizationContrast);
theSceneEngine.visualizeStaticFrame(theSceneSequence);

% Estimate threshold and plot/report results.  This
% does a maximumu likelihood based on the trials run, and is not subject to
% the discretization used by QUEST+.
subplot(4, 4, index * 2);
[logThreshold, para] = estimator.thresholdMLE('showPlot', true, 'pointSize', 2.5);
fprintf('Maximum likelihood fit parameters: %0.2f, %0.2f, %0.2f, %0.2f\n', ...
    para(1), para(2), para(3), para(4));

end

%% Compute performance for the classifier
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
