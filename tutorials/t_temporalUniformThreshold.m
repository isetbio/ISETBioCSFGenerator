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
%   t_responseClassifier
%   t_neuralResponseCompute
%   t_sceneGeneration

%% Initialization

% Instantiate a sceneGenerationEngine with
% function handle that generates uniform field temporal modulation
theSceneEngine = sceneEngine(@uniformFieldTemporalModulation);

% Instantiate a neuralResponseEngine. No responseParams passed, so we
% are using the default params specified photopigmentExcitationsWithNoEyeMovements
theNeuralEngine = neuralResponseEngine(@photopigmentExcitationsWithNoEyeMovements);

% Instantiate a responseClassifierEngine with pcaSVMClassifier
theClassifierEngine = responseClassifierEngine(@pcaSVMClassifier);

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
estimator = QuestThresholdEngine('minTrial', 1280, 'maxTrial', 1280, ...
    'estDomain', estDomain, 'slopeRange', slopeRange, ...
    'numEstimator', 1, 'stopCriterion', 0.05);

%% Threshold estimation with QUEST+

[logContrast, nextFlag] = estimator.nextStimulus();

% 128 trial for each contrast level
% 512 trials for training and testing
nRepeat = 128;
nTrainSample = 512;
while (nextFlag)
    
    % log contrast -> contrast
    % Compute the TEST stimulus
    testContrast = 10 ^ logContrast;
    [theTestSceneSequence, ~] = theSceneEngine.compute(testContrast);
    
    fprintf('Current test contrast: %.4f ', testContrast);
    
    % Code for observer (classifier or other form)
    % Compute many respose instances to the NULL and TEST stimuli for training the SVM classifier
    
    % NULL stimulus instances
    [inSampleNullStimResponses, ~] = theNeuralEngine.compute(...
        theNullSceneSequence, ...
        theSceneTemporalSupportSeconds, ...
        nTrainSample, ...
        'noiseFlags', {'random'});
    
    % TEST stimulus instances
    [inSampleTestStimResponses, ~] = theNeuralEngine.compute(...
        theTestSceneSequence, ...
        theSceneTemporalSupportSeconds, ...
        nTrainSample, ...
        'noiseFlags', {'random'});
    
    % Train the binary classifier on the above NULL/TEST response set
    trainingData = theClassifierEngine.compute('train',...
        inSampleNullStimResponses('random'), ...
        inSampleTestStimResponses('random'));
    
    pCorrect = trainingData.pCorrectInSample;
    fprintf('P-correct: %.3f \n', pCorrect);
    
    % Translate p-correct into binary response using binomial distribution
    % Not necessary for template-based observer
    nCorrect = binornd(nRepeat, pCorrect);
    response = [zeros(1, nRepeat - nCorrect), ones(1, nCorrect)];
    
    [logContrast, nextFlag] = ...
        estimator.multiTrial(logContrast * ones(1, nRepeat), response(randperm(length(response))));
    
end

%% Show results
fprintf('%d trials recorded \n', estimator.nTrial);

% Plot results
figure(); subplot(1, 2, 1);
[threshold, para] = estimator.thresholdMLE('showPlot', true, 'pointSize', 10);

fprintf('Maximum likelihood fit parameters: %0.2f, %0.2f, %0.2f, %0.2f\n', ...
    para(1), para(2), para(3), para(4));

%% Validation by computing the whole psychometric curve
logContrast = -3.5 : 0.125 : -1;
pCorrect = zeros(1, length(logContrast));

nTrainSample = 512;
parfor idx = 1 : length(logContrast)
    testContrast = 10 ^ logContrast(idx);
    [theTestSceneSequence, ~] = theSceneEngine.compute(testContrast);
    
    % NULL stimulus instances
    [inSampleNullStimResponses, ~] = theNeuralEngine.compute(...
        theNullSceneSequence, ...
        theSceneTemporalSupportSeconds, ...
        nTrainSample, ...
        'noiseFlags', {'random'});
    
    % TEST stimulus instances
    [inSampleTestStimResponses, ~] = theNeuralEngine.compute(...
        theTestSceneSequence, ...
        theSceneTemporalSupportSeconds, ...
        nTrainSample, ...
        'noiseFlags', {'random'});
    
    % Train the binary classifier on the above NULL/TEST response set
    trainingData = theClassifierEngine.compute('train',...
        inSampleNullStimResponses('random'), ...
        inSampleTestStimResponses('random'));
    
    pCorrect(idx) = trainingData.pCorrectInSample;
end

subplot(1, 2, 2);
plot(logContrast, pCorrect, '-ok', 'LineWidth', 2);
ylim([0, 1]);
