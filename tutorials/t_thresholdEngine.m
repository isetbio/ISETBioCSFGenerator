% Combine everything we have to compute threshold of a simple temporal modulation
%
% Syntax:
%    t_temporalThreshold
%
% Description:
%    Demonstrates how to compute threshold using our framework (i.e.,
%    scene, neural, classifier, and threshold engine).
%
% Inputs:
%    None.
%
% Outputs:
%   None.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%   t_thresholdClassifier

% History:
%   10/05/20  dhb   Add logic to cache trained classifiers so we don't
%                   train them again.

%% Initialization
clear; close all;

%% Instantiate a sceneGenerationEngine
% Choices are:
%   'sceUniformFieldModulation
whichSceneEngine = 'sceUniformFieldModulation';
switch (whichSceneEngine)
    case 'sceUniformFieldModulation'
        % Spatially uniform field temporal modulation
        %
        % Note use of call without arguments to get default parameters
        % struct, and adjustment of size to make it small.  This speeds
        % things up for this demo.
        sceneParams = sceUniformFieldTemporalModulation;
        sceneParams.sizePixels = 5;
        theSceneEngine = sceneEngine(@sceUniformFieldTemporalModulation,sceneParams);
        
    otherwise
        error('Unknown scene engine specified');
end

%% Instantiate a neuralResponseEngine.
% Choices are:
%   'ncePhotopigmentExcitationsWithNoEyeMovements'
%   'nceScenePhotonNoise'
%
% The actual threshold varies enough with the different engines that we
% need to adjust the contrast range that Quest+ searches over, as well as
% the range of psychometric function slopes.
whichNeuralEngine = 'ncePhotopigmentExcitationsWithNoEyeMovements';
switch (whichNeuralEngine)
    case 'ncePhotopigmentExcitationsWithNoEyeMovements'
        % Basic retinal image formation and sampling by the cone mosaic.
        % Note use of neural engine to get its own default parameters and
        % adjust them.  Smaller field of view speeds things up.
        neuralParams = nrePhotopigmentExcitationsWithNoEyeMovements;
        neuralParams.coneMosaicParams.fovDegs = 0.1;
        theNeuralEngine = neuralResponseEngine(@nrePhotopigmentExcitationsWithNoEyeMovements,neuralParams);
        logThreshLimitLow = 4;
        logThreshLimitHigh = 1;
        logThreshLimitDelta = 0.05;
        slopeRangeLow = 1;
        slopeRangeHigh = 100;
        slopeDelta = 1;
        
    case 'nceScenePhotonNoise'
        % Add Poisson noise to the photon counts.
        theNeuralEngine = neuralResponseEngine(@nreScenePhotonNoise);
        logThreshLimitLow = 7;
        logThreshLimitHigh = 5;
        logThreshLimitDelta = 0.005;
        slopeRangeLow = 100;
        slopeRangeHigh = 10000;
        slopeDelta = 100;
        
    otherwise
        error('Unknown neural engine specified');
end

%% Instantiate a responseClassifierEngine
% Choices are:
%   'rcePoissonTAFC'
%   'rcePcaSVMTAFC'
whichObserver = 'rcePoissonTAFC';

% A larger nTest is usually more effective, but depending on the performance
% bottleneck of your observer, you might consider a smaller nTest.  If it
% is fast to compute the classifier for a given contrast and slow to make
% predictions for a trial, then small nTest will be better.  If
% it's slow to build a classifier for a given contrast and fast to make
% predictions, then large nTest will be better.
switch whichObserver
    case 'rcePoissonTAFC'
        % The ideal observer for a TAFC task limited by Poisson noise
        theRawClassifierEngine = responseClassifierEngine(@rcePoissonTAFC);
        
        % Noise-free instance for training, random for test
        trainFlag = 'none'; testFlag = 'random';
        nTrain = 1; nTest = 120;
        
    case 'rcePcaSVMTAFC'
        % SVM classifier wit PCA pre-processing
        theRawClassifierEngine = responseClassifierEngine(@rcePcaSVMTAFC);
        
        % Noisy instances for training and testing
        trainFlag = 'random'; testFlag = 'random';
        nTrain = 512; nTest = 120;
        
    otherwise
        error('Unknown observer specified');
end

% Generate and compute the zero contrast NULL stimulus sequence
nullContrast = 0.0;
[theNullSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(nullContrast);

%% Construct a QUEST threshold estimator estimate threshold on log contrast
estDomain  = -logThreshLimitLow : logThreshLimitDelta : -logThreshLimitHigh;
slopeRange = slopeRangeLow: slopeDelta : slopeRangeHigh;

% There are two recommended ways to setup the QUEST+ threshold engine.
%   'fixedNumber'    - run a fixed number of trials
%   'adaptiveMode'   - run until estimate reaches specified precision.
% See below for more.
questMode = 'adaptiveMode';


switch questMode
    case 'fixedNumber'
        % Run fixed number of trials.  This is done by setting 'minTrial' and
        % maxTrial values to be the same and running a single Quest+
        % object.
        estimator = questThresholdEngine('minTrial', 1e3, 'maxTrial', 1e3, ...
            'estDomain', estDomain, 'slopeRange', slopeRange, 'numEstimator', 1);
        
    case 'adaptiveMode'
        % Run 'numEstimator' interleaved Quest+ object. In this case, the
        % threshold engine calculates the running standard error (SE) among
        % those objects. The stopping criterion is triggered when 1) total number of trials >=
        % 'minTrial' AND the SE/Threshold < 'stopCriterion', OR 2) when total number
        % of trials >= 'maxTrial'
        estimator = questThresholdEngine('minTrial', 1e2, 'maxTrial', 1e4, ...
            'estDomain', estDomain, 'slopeRange', slopeRange, ...
            'numEstimator', 4, 'stopCriterion', 0.025);
        
    otherwise
        error('Unknown threshold engine mode specified');
end

%% Threshold estimation with QUEST+
[logContrast, nextFlag] = estimator.nextStimulus();
testedContrasts = [];
while (nextFlag)
    
    % log contrast -> contrast
    % Compute the TEST stimulus
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
        % responses
        [response, theTrainedClassifierEngines{testedIndex}] = computeResponse(...
            theNullSceneSequence, theTestSceneSequences{testedIndex}, ...
            theSceneTemporalSupportSeconds, nTrain, nTest, ...
            theNeuralEngine, theRawClassifierEngine, trainFlag, testFlag);
        
    else
        % Classifier is already trained, just get responses
        response = computeResponse(...
            theNullSceneSequence, theTestSceneSequences{testedIndex}, ...
            theSceneTemporalSupportSeconds, nTrain, nTest, ...
            theNeuralEngine, theTrainedClassifierEngines{testedIndex}, [], testFlag);
    end
    
    % Report what happened
    fprintf('Current test contrast: %g, P-correct: %g \n', testContrast, mean(response));
    
    % Get next stimulus contrast
    [logContrast, nextFlag] = ...
        estimator.multiTrial(logContrast * ones(1, nTest), response);
    
    % Get current threshold estimate
    [threshold, stderr] = estimator.thresholdEstimate();
    fprintf('Current threshold estimate: %g, stderr: %g \n', 10 ^ threshold, stderr);
end

%% Show results
fprintf('%d trials recorded \n', estimator.nTrial);

% Estimate threshold and plot/report results
figure();
[threshold, para] = estimator.thresholdMLE('showPlot', true, 'pointSize', 7.5);
fprintf('Maximum likelihood fit parameters: %0.2f, %0.2f, %0.2f, %0.2f\n', ...
    para(1), para(2), para(3), para(4));

%% Validation by computing the entire psychometric curve
%
% Takes a rather long time to run, change runValidation to 'true' if you
% want to run the validation
runValidation = false;
if (runValidation)
    
    logContrast = -logThreshLimitLow : logThreshLimitDelta : -logThreshLimitHigh;
    pCorrect = zeros(1, length(logContrast));
    parfor idx = 1:length(logContrast)
        testContrast = 10 ^ logContrast(idx);
        [theTestSceneSequence, ~] = theSceneEngine.compute(testContrast);
        
        pCorrect(idx) = mean(computeResponse(...
            theNullSceneSequence, theTestSceneSequence, ...
            theSceneTemporalSupportSeconds, 512, 512, ...
            theNeuralEngine, theRawClassifierEngine), trainFlag, testFlag);
    end
    
    hold on;
    plot(logContrast, pCorrect, '-ok', 'LineWidth', 1);
    ylim([0, 1]);
    
end

%% Helper function
function [response,theClassifierEngine] = computeResponse(nullScene, testScene, temporalSupport, nTrain, nTest, theNeuralEngine, theClassifierEngine, trainFlag, testFlag)

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
    
    % Train the classifier
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
response = dataOut.trialPredictions;

end
