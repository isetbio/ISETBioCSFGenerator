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

% Instantiate a sceneGenerationEngine with
% function handle that generates uniform field temporal modulation
theSceneEngine = sceneEngine(@sceUniformFieldTemporalModulation);

% Instantiate a neuralResponseEngine. No responseParams passed, so we
% are using the default params specified nrePhotopigmentExcitationsWithNoEyeMovements
theNeuralEngine = neuralResponseEngine(@nrePhotopigmentExcitationsWithNoEyeMovements);

% Instantiate a responseClassifierEngine with rcePoissonTAFC which is the ideal observer for 2AFC task.
% It is rather simple to change components of the pipeline (i.e., observer), without changing other aspects of
% code (e.g., stimulus generation). Here, we can also choose to use rcePcaSVMTAFC which is a SVM based observer

%% Use rcePoissonTAFC (the ideal observer of the task)
observer_ID = 1;

%% Or use rcePcaSVMTAFC observer
observer_ID = 2;

%% Initialization, continued
% A larger nTest is usually more effective, but depending on the performance
% bottleneck of your observer, you might consider a smaller N
switch observer_ID
    case 1        
        theClassifierEngine = responseClassifierEngine(@rcePoissonTAFC);
        trianFlag = 'none'; testFlag = 'random';
        nTrian = 1; nTest = 120;
        
    case 2        
        theClassifierEngine = responseClassifierEngine(@rcePcaSVMTAFC);
        trianFlag = 'random'; testFlag = 'random';
        nTrian = 512; nTest = 120;
end

% Generate and compute the zero contrast NULL stimulus (sequence)
nullContrast = 0.0;
[theNullSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(nullContrast);

% Construct a QUEST threshold estimator
% estimate threshold on log contrast
estDomain  = -4 : 0.02 : -1;
slopeRange = 1.0 : 1.0 : 100;

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
        theSceneTemporalSupportSeconds, nTrian, nTest, ...
        theNeuralEngine, theClassifierEngine, trianFlag, testFlag);
    
    fprintf('Current test contrast: %.4f, P-correct: %.4f \n', testContrast, mean(response));
    
    [logContrast, nextFlag] = ...
        estimator.multiTrial(logContrast * ones(1, nTest), response);
    
    [threshold, stderr] = estimator.thresholdEstimate();
    fprintf('Current threshold estimate: %.4f, stderr: %.4f \n', 10 ^ threshold, stderr);
    
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
            theNeuralEngine, theClassifierEngine), trianFlag, testFlag);
    end
    
    hold on;
    plot(logContrast, pCorrect, '-ok', 'LineWidth', 1);
    ylim([0, 1]);
    
end

%% Helper function
function response = computeResponse(nullScene, testScene, temporalSupport, nTrain, nTest, theNeuralEngine, theClassifierEngine, trianFlag, testFlag)

% Code for observer (classifier or other form)
% Compute many respose instances to the NULL and TEST stimuli for training the SVM classifier

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
