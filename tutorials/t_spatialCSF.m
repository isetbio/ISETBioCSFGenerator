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

%% Create neural engine
% Instantiate a neuralResponseEngine
neuralParams = nrePhotopigmentExcitationsWithNoEyeMovements;
neuralParams.coneMosaicParams.fovDegs = 0.25;
theNeuralEngine = neuralResponseEngine(@nrePhotopigmentExcitationsWithNoEyeMovements, neuralParams);

%% Instantiate the PoissonTAFC responseClassifierEngine
% PoissonTAFC makes decision by performing the Poisson likelihood ratio test
classifierEngine = responseClassifierEngine(@rcePoissonTAFC);

% Parameter associated with this classifier
classifierPara.trainFlag = 'none'; classifierPara.testFlag = 'random';
classifierPara.nTrain = 1; classifierPara.nTest = 128;

%% Parameter for threshold estimation/quest engine
% The actual threshold varies enough with the different engines that we
% need to adjust the contrast range that Quest+ searches over, as well as
% the range of psychometric function slopes.
thresholdPara.logThreshLimitLow = 2.5;
thresholdPara.logThreshLimitHigh = 0;
thresholdPara.logThreshLimitDelta = 0.02;

thresholdPara.slopeRangeLow = 1;
thresholdPara.slopeRangeHigh = 100;
thresholdPara.slopeDelta = 2.5;

questEnginePara.minTrial = 1280;
questEnginePara.maxTrial = 1280;
questEnginePara.numEstimator = 1;

%% Compute threshold for each spatial frequency
% See toolbox/helpers for the definitin of
% function 'createGratingScene' and
% function 'computeThreshold'

dataFig = figure();
for idx = 1:length(spatialFreqs)
    % create a static grating scene with a particular chromatic direction,
    % spatial frequency, and temporal duration
    gratingScene = createGratingScene('chromaDir', chromaDir, 'spatialFreq', spatialFreqs(idx));
    [logThreshold(idx), questObj] = computeThreshold(gratingScene, theNeuralEngine, classifierEngine, classifierPara, thresholdPara, questEnginePara);
    
    % Plot stimulus
    figure(dataFig);
    subplot(4, 4, idx * 2 - 1);
    
    visualizationContrast = 1.0;
    [theSceneSequence] = gratingScene.compute(visualizationContrast);
    gratingScene.visualizeStaticFrame(theSceneSequence);
    
    %  and data and psychometric curve with a marker size of 2.5
    subplot(4, 4, idx * 2);
    questObj.plotMLE(2.5);
end
set(gcf, 'Position',  [0, 0, 800, 800]);

% Convert returned log threshold to linear threshold
threshold = 10 .^ logThreshold;

%% Plot Contrast Sensitivity Function
theCsfFig = figure();
loglog(spatialFreqs, 1 ./ threshold, '-ok', 'LineWidth', 2);
xlabel('Spatial Frequency (cyc/deg)');
ylabel('Sensitivity');
set(gcf, 'Position',  [0, 0, 600, 800]);

%% Helper functions for calculating threshold and classifier performance

% Compute threshold for a particular chromatic direction and spatial frequency
% Chromatic direction is a 1-by-3 vector specifying contrast on the L, M and S Cone, respectively
% Spatial frequency is a number in the unit of cycles per degree
function [logThreshold, questObj] = computeThreshold(theSceneEngine, theNeuralEngine, classifierEngine, classifierPara, thresholdPara, questEnginePara)

% Construct a QUEST threshold estimator estimate threshold on log contrast
% Run a fixed number of trials (e.g., 10 contrast level, 1280 trials in total)
estDomain  = -thresholdPara.logThreshLimitLow : thresholdPara.logThreshLimitDelta : -thresholdPara.logThreshLimitHigh;
slopeRange = thresholdPara.slopeRangeLow: thresholdPara.slopeDelta : thresholdPara.slopeRangeHigh;

estimator = questThresholdEngine('minTrial', questEnginePara.minTrial, 'maxTrial', questEnginePara.maxTrial, ...
    'estDomain', estDomain, 'slopeRange', slopeRange, 'numEstimator', questEnginePara.numEstimator);

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
            theSceneTemporalSupportSeconds, classifierPara.nTrain, classifierPara.nTest, ...
            theNeuralEngine, classifierEngine, classifierPara.trainFlag, classifierPara.testFlag);
        
    else
        % Classifier is already trained, just get predictions
        predictions = computePerformance(...
            theNullSceneSequence, theTestSceneSequences{testedIndex}, ...
            theSceneTemporalSupportSeconds, classifierPara.nTrain, classifierPara.nTest, ...
            theNeuralEngine, theTrainedClassifierEngines{testedIndex}, [], classifierPara.testFlag);
    end
    
    % Tell QUEST+ what we ran (how many trials at the given contrast) and
    % get next stimulus contrast to run.
    [logContrast, nextFlag] = ...
        estimator.multiTrial(logContrast * ones(1, classifierPara.nTest), predictions);
    
    % Report what happened
    % fprintf('Current test contrast: %g, P-correct: %g \n', testContrast, mean(predictions));
end

[logThreshold, para] = estimator.thresholdMLE('showPlot', false);
fprintf('Maximum likelihood fit parameters: %0.2f, %0.2f, %0.2f, %0.2f\n', ...
    para(1), para(2), para(3), para(4));

questObj = estimator;

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
