% Combine everything we have to compute threshold of a temporal modulation
%
% Syntax:
%    t_thresholdEngine
%
% Description:
%    Demonstrates how to compute threshold using our framework (i.e.,
%    scene, neural, classifier, and threshold engine).
%
%    In addition to demonstrating usage, this illustrates how various
%    components of the pipeline can be interchanged without affecting the
%    rest of them.
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
%   t_sceneGeneration, t_modulatedGratingsSceneGeneration,
%   t_neuralResponseCompute, t_responseClassifier
%
% History:
%   10/05/20  dhb   Add logic to cache trained classifiers so we don't
%                   train them again.

%% Initialization
clear; close all;

%% Instantiate a sceneGenerationEngine
%
% The first step is to define the spatial and temporal parameters of a
% scene.  We do this by instantiating a sceneEngine object, passing it a
% function that is responsible for defining the scene as well as a struct
% of parameters that that function accepts.
%
% The name of passed function should by convention begin with the prefix
% 'sce'.
%
% Also by convention, the passed function called with no arguments returns
% a struct of default parameters, making it easy to see what can be
% controlled.
%
% Once the sceneEngine has been instantiated, its compute method may be
% called with a standard set of arguments and return a standard set of
% values, so that the details of what it does are transparent to other
% steps in the pipeline.  See examples below for this usage.  The compute
% method needs to know how to vary the parameter that is being varied to
% find threshold. The canonical usage is to vary contrast, but nothing in
% the code actually cares about the semantics.
%
% See t_sceneGeneration and t_modulatedGratingSceneGeneration for tutorials
% that focus on how to use the sceneEngine.  See the functions
% sceUniformFieldTemporalModulation and sceGrating for example
% implementations of compute functions to use with the sceneEngine object.
%
% Choices of passed functions that can be used in this tutorial are:
%   'sceUniformFieldModulation'
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
        
        % Instantiate the sceneEngine object
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
%
% responseClassifierEngines are responsible for simulating observer
% performance on a psychophysical task, given noisy samples of the output of the
% neuralResponseEngine.  As with our other objects, they are instantiated
% with a function that does this work.  This function must be able to train
% the classifier given samples of the neural responses, and then predict
% on a simulated trial-by-trial basis whether that trial was correct or
% incorrect.  
%
% For the most part, we simulate two-alternative forced-choice trials, but
% we think the framework will work for at least some other psychophysical
% tasks.
%
% The same general conventions apply to responseClassifierEngines as to our
% other objects. Usage is illustrated below.
%
% A larger nTest is usually more effective, but depending on the performance
% bottleneck of your observer, you might consider a smaller nTest.  If it
% is fast to compute the classifier for a given contrast and slow to make
% predictions for a trial, then small nTest will be better.  If
% it's slow to build a classifier for a given contrast and fast to make
% predictions, then large nTest will be better.
%
% Choices are:
%   'rcePoissonTAFC'
%   'rcePcaSVMTAFC'
whichObserver = 'rcePoissonTAFC';
switch whichObserver
    case 'rcePoissonTAFC'
        % The ideal observer for a TAFC task limited by Poisson noise.
        % This classifier doesn't take any parameters.
        theRawClassifierEngine = responseClassifierEngine(@rcePoissonTAFC);
        
        % Below we use a wrapper routine to train the classifier and to
        % predict trial-by-trial responses.  Because that routine is
        % general, we define some parameters for it here.
        %
        % We'll illustrate a signal know exactly classifier, using the trainFlag value
        % of 'none' to indicate that the classifier should be trained with a noise free
        % version of the stimuli.  We'll test with noisy values, however.
        %
        % Similarly, since this is a template type of classifer, we can just pass one
        % noise free examplar of the NULL and TEST stimuli.  We choose here to predict
        % performance on 120 stimuli for each contrast tested.  This is because classification
        % of individual stimuli is pretty fast for this classifier, and we might as well
        % get a good estimate of performance for each contrast tested.
        trainFlag = 'none'; testFlag = 'random';
        nTrain = 1; nTest = 120;
        
    case 'rcePcaSVMTAFC'
        % SVM linear classifier with PCA pre-processing.  The usage here obtains
        % the default parameters and then passes them in. You could adjust
        % the parameters if you wanted to, for example changing the number
        % of PCA components retained for the SVM processing stage.
        rcePcaSVMTAFCParams = rcePcaSVMTAFC;
        theRawClassifierEngine = responseClassifierEngine(@rcePcaSVMTAFC,rcePcaSVMTAFCParams);
        
        % Here we use noisy exemplars to train the SVM classifier, and
        % again predict in batches of 120 trials.
        trainFlag = 'random'; testFlag = 'random';
        nTrain = 512; nTest = 120;
        
    otherwise
        error('Unknown observer specified');
end

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
        % Run 'numEstimator > 1' interleaved Quest+ objects. In this case, the
        % threshold engine calculates the running standard error (SE) among
        % those objects. The stopping criterion is triggered when 
        % 1) total number of trials >= 'minTrial' 
        % AND the 'stopCriterion'(threshold, SE) is TRUE, 
        % OR 2) when total numberof trials >= 'maxTrial'.
        
        % 'stopCriterion' could be one of two options:
        
        % 1) A single number. In this case, the criterion will simply be
        % SE < 'stopCriterion'.
        stopCriterion = 0.025;
        
        % 2) A function handle that takes the current estimate of threshold
        % and SE estimates as input arguments, and returns a boolean variable.
        % For example: we can use a relative criterion w.r.t. the magnitude of threshold        
        stopCriterion = @(threshold, se) se / abs(threshold) < 0.01;
        
        estimator = questThresholdEngine('minTrial', 2e2, 'maxTrial', 5e3, ...
            'estDomain', estDomain, 'slopeRange', slopeRange, ...
            'numEstimator', 4, 'stopCriterion', stopCriterion);
        
    otherwise
        error('Unknown threshold engine mode specified');
end

%% Generate the NULL scene 
%
% Threshold will be measured a a perturbation from this scene.  Typically
% it will correspond to zero contrast, but it doesn't have to.
%
% This call illustrates the compute method of the sceneEngine class.  We
% pass the desired contrast, and back comes a cell array as the first
% returned value, with one ISETBio scene for each time point.  The
% corresponding times are in the second returned value, an array.  Times
% are in seconds.
nullContrast = 0.0;
[theNullSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(nullContrast);

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
