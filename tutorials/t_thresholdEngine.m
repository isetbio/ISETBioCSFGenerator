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

%% Initialization
clear; close all;

%% Instantiate a sceneGenerationEngine with
% function handle that generates uniform field temporal modulation
sceneParams = sceUniformFieldTemporalModulation;
sceneParams.sizePixels = 5;
theSceneEngine = sceneEngine(@sceUniformFieldTemporalModulation,sceneParams);

%% Instantiate a neuralResponseEngine.
% Choices are:
%   'ncePhotopigmentExcitationsWithNoEyeMovements'
%   'nceScenePhotonNoise'
whichSceneEngine = 'ncePhotopigmentExcitationsWithNoEyeMovements';
switch (whichSceneEngine)
    case 'ncePhotopigmentExcitationsWithNoEyeMovements'
        neuralParams = nrePhotopigmentExcitationsWithNoEyeMovements;
        neuralParams.coneMosaicParams.fovDegs = 0.1;
        theNeuralEngine = neuralResponseEngine(@nrePhotopigmentExcitationsWithNoEyeMovements,neuralParams);
        logThreshLimitLow = 4;
        logThreshLimitHigh = 1;
        logThreshLimitDelta = 0.02;
        slopeRangeLow = 1;
        slopeRangeHigh = 100;
        slopeDelta = 1;
        
    case 'nceScenePhotonNoise'
        theNeuralEngine = neuralResponseEngine(@nreScenePhotonNoise);
        logThreshLimitLow = 8;
        logThreshLimitHigh = 4;
        logThreshLimitDelta = 0.02;
        slopeRangeLow = 1000;
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
% bottleneck of your observer, you might consider a smaller nTest
switch whichObserver
    case 'rcePoissonTAFC' 
        % The ideal observer for a TAFC task limited by Poisson noise
        theClassifierEngine = responseClassifierEngine(@rcePoissonTAFC);
        
        % Noise-free instance for training, random for test
        trainFlag = 'none'; testFlag = 'random';
        nTrain = 1; nTest = 120;
        
    case 'rcePcaSVMTAFC'  
        % SVM classifier wit PCA pre-processing
        theClassifierEngine = responseClassifierEngine(@rcePcaSVMTAFC);
        
        % Noisy instances for training and testing
        trainFlag = 'random'; testFlag = 'random';
        nTrain = 512; nTest = 120;
        
    otherwise
        error('Unknown observer specified');
end

% Generate and compute the zero contrast NULL stimulus (sequence)
nullContrast = 0.0;
[theNullSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(nullContrast);

% Construct a QUEST threshold estimator
% estimate threshold on log contrast
estDomain  = -logThreshLimitLow : logThreshLimitDelta : -logThreshLimitHigh;
slopeRange = slopeRangeLow: slopeDelta : slopeRangeHigh;

% Run QUEST for a total of minTrial trials. For typical usage, this is recommended

% However, when numEstimator > 1, an adpative procedure is invoked,
% for which the routine will stop if the S.E. among N parallel quest+ object
% is below the threshold specified as stopCriterion
estimator = questThresholdEngine('minTrial', 1e3, 'maxTrial', 1e4, ...
    'estDomain', estDomain, 'slopeRange', slopeRange, ...
    'numEstimator', 3, 'stopCriterion', 0.025);

%% Threshold estimation with QUEST+
[logContrast, nextFlag] = estimator.nextStimulus();
while (nextFlag)
    
    % log contrast -> contrast
    % Compute the TEST stimulus
    testContrast = 10 ^ logContrast;
    
    % Test scene generate
    [theTestSceneSequence, ~] = theSceneEngine.compute(testContrast);
    
    response = computeResponse(theNullSceneSequence, theTestSceneSequence, ...
        theSceneTemporalSupportSeconds, nTrain, nTest, ...
        theNeuralEngine, theClassifierEngine, trainFlag, testFlag);
    
    fprintf('Current test contrast: %g, P-correct: %g \n', testContrast, mean(response));
    
    [logContrast, nextFlag] = ...
        estimator.multiTrial(logContrast * ones(1, nTest), response);
    
    [threshold, stderr] = estimator.thresholdEstimate();
    fprintf('Current threshold estimate: %g, stderr: %g \n', 10 ^ threshold, stderr);
    
end

%% Show results
fprintf('%d trials recorded \n', estimator.nTrial);

% Plot results
figure();
[threshold, para] = estimator.thresholdMLE('showPlot', true, 'pointSize', 7.5);

fprintf('Maximum likelihood fit parameters: %0.2f, %0.2f, %0.2f, %0.2f\n', ...
    para(1), para(2), para(3), para(4));

%% Validation by computing the entire psychometric curve
% Takes a rather long time to run, change runValidation to 'true' if you
% want to run the validation

runValidation = false;

if (runValidation)
    
    logContrast = -4 : 0.05 : -1;
    pCorrect = zeros(1, length(logContrast));
    
    parfor idx = 1:length(logContrast)
        testContrast = 10 ^ logContrast(idx);
        [theTestSceneSequence, ~] = theSceneEngine.compute(testContrast);
        
        pCorrect(idx) = mean(computeResponse(...
            theNullSceneSequence, theTestSceneSequence, ...
            theSceneTemporalSupportSeconds, 512, 512, ...
            theNeuralEngine, theClassifierEngine), trainFlag, testFlag);
    end
    
    hold on;
    plot(logContrast, pCorrect, '-ok', 'LineWidth', 1);
    ylim([0, 1]);
    
end

%% Helper function
function response = computeResponse(nullScene, testScene, temporalSupport, nTrain, nTest, theNeuralEngine, theClassifierEngine, trianFlag, testFlag)

% Generate stimulus for training
% NULL stimulus mean response
[inSampleNullStimResponses, ~] = theNeuralEngine.compute(...
    nullScene, ...
    temporalSupport, ...
    nTrain, ...
    'noiseFlags', {trianFlag});

% TEST stimulus instances
[inSampleTestStimResponses, ~] = theNeuralEngine.compute(...
    testScene, ...
    temporalSupport, ...
    nTrain, ...
    'noiseFlags', {trianFlag});

% Run a ideal observer binary classifier on the above NULL/TEST response set
% 'Training' is basically storing the noise-free template
theClassifierEngine.compute('train', ...
    inSampleNullStimResponses(trianFlag), ...
    inSampleTestStimResponses(trianFlag));

% Generate stimulus for testing/prediction
% 'Predict' on noisy responses
% NULL stimulus mean response
[inSampleNullStimResponses, ~] = theNeuralEngine.compute(...
    nullScene, ...
    temporalSupport, ...
    nTest, ...
    'noiseFlags', {testFlag});

% TEST stimulus instances
[inSampleTestStimResponses, ~] = theNeuralEngine.compute(...
    testScene, ...
    temporalSupport, ...
    nTest, ...
    'noiseFlags', {testFlag});

dataOut = theClassifierEngine.compute('predict', ...
    inSampleNullStimResponses(testFlag), ...
    inSampleTestStimResponses(testFlag));

response = dataOut.trialPredictions;

end
