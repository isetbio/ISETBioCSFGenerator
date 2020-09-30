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
theSceneEngine = sceneEngine(@uniformFieldTemporalModulation);

% Instantiate a neuralResponseEngine. No responseParams passed, so we
% are using the default params specified photopigmentExcitationsWithNoEyeMovements
theNeuralEngine = neuralResponseEngine(@photopigmentExcitationsWithNoEyeMovements);

% Instantiate a responseClassifierEngine with poissonTemplateClassifier
theClassifierEngine = responseClassifierEngine(@poissonTemplateClassifier);

% Generate and compute the zero contrast NULL stimulus (sequence)
nullContrast = 0.0;
[theNullSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(nullContrast);

% Construct a QUEST threshold estimator
% estimate threshold on log contrast
estDomain  = -3.5 : 0.01 : -1;
slopeRange = 1.0 : 1.0 : 200;

% Run QUEST for a total of minTrial trials. For typical usage, this is recommended

% However, when numEstimator > 1, an adpative procedure is invoked,
% for which the routine will stop if the S.E. among N parallel quest+ object
% is below the threshold specified as stopCriterion
estimator = QuestThresholdEngine('minTrial', 640, 'maxTrial', 640, ...
    'estDomain', estDomain, 'slopeRange', slopeRange, ...
    'numEstimator', 1, 'stopCriterion', 0.025);

%% Threshold estimation with QUEST+

[logContrast, nextFlag] = estimator.nextStimulus();

% 64 trial for each contrast level
nRepeat = 64;
while (nextFlag)
    
    % log contrast -> contrast
    % Compute the TEST stimulus
    testContrast = 10 ^ logContrast;
    [theTestSceneSequence, ~] = theSceneEngine.compute(testContrast);
    
    % Code for observer (classifier or other form)
    % Compute many respose instances to the NULL and TEST stimuli for training the SVM classifier
    
    % NULL stimulus mean response
    [inSampleNullStimResponses, ~] = theNeuralEngine.compute(...
        theNullSceneSequence, ...
        theSceneTemporalSupportSeconds, ...
        nRepeat, ...
        'noiseFlags', {'none', 'random'});
    
    % TEST stimulus instances
    [inSampleTestStimResponses, ~] = theNeuralEngine.compute(...
        theTestSceneSequence, ...
        theSceneTemporalSupportSeconds, ...
        nRepeat, ...
        'noiseFlags', {'none', 'random'});
    
    % run binary classifier on the above NULL/TEST response set (no
    % training required since it is template based)
    dataOut = theClassifierEngine.compute('predict', inSampleNullStimResponses, inSampleTestStimResponses);
    response = dataOut.response;    
    
    % Translate p-correct into binary response using binomial distribution
    % Not necessary for template-based observer
    nCorrect = binornd(nRepeat, pCorrect);
    response = [zeros(1, nRepeat - nCorrect), ones(1, nCorrect)];
    
    fprintf('Current test contrast: %.4f, P-correct: %.4f \n', testContrast, pCorrect);
    
    [logContrast, nextFlag] = ...
        estimator.multiTrial(logContrast * ones(1, nRepeat), response(randperm(length(response))));
    
    % Current threshold estimate and its standard error
    [threshold, stderr] = estimator.thresholdEstimate();
    fprintf('Current threshold estimate: %.4f, stderr: %.4f \n', threshold, stderr);
end

%% Show results
fprintf('%d trials recorded \n', estimator.nTrial);

% Plot results
figure();
[threshold, para] = estimator.thresholdMLE('showPlot', true, 'pointSize', 7.5);

fprintf('Maximum likelihood fit parameters: %0.2f, %0.2f, %0.2f, %0.2f\n', ...
    para(1), para(2), para(3), para(4));
