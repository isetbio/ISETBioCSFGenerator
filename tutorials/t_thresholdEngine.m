% Combine everything we have to compute threshold for a modulation
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
%   t_neuralResponseCompute, t_responseClassifier, computePerformance
%
% History:
%   10/05/20  dhb   Add logic to cache trained classifiers so we don't
%                   train them again.

%% Initialization
clear; close all;

%% Instantiate a sceneEngine
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
        sceneParams.stimDurationFramesNum = 1;
        sceneParams.stimOnsetFramesIndices = 1;
        sceneParams.frameDurationSeconds = 0.1;
        
        % Instantiate the sceneEngine object
        theSceneEngine = sceneEngine(@sceUniformFieldTemporalModulation,sceneParams);
        
    otherwise
        error('Unknown scene engine specified');
end

%% Instantiate a neuralResponseEngine.
%
% The neural response engine is our model of visual processing, code that
% takes a scene sequence as input and produces the neural representation as
% output.  The scene sequence is produced by the sceneEngine object, and
% the neural responses are provided to a classifier object for making
% perceptual decisions.
%
% The structure of the neural response engine is conceptually similar to
% that of the scene engine.  The user provides a function that does the
% computations, and instantiates the neural response engine with this
% function.  The function takes a struct of parameters, which may be
% included with the instantiation.  By convention, the function will return
% a struct of its default arguments if called with no inputs.
%
% The way we have set up the neuralResponseEngine, the passed function
% takes two important key/value pairs.
%    'noiseFlags'                   - Cell array of strings containing labels
%                                     that encode the type of noise to be included
%                                     Valid values are: 
%                                        - 'none' (noise-free responses)
%                                        - 'random' (noisy response instances)
%                                     Default is {'random'}.  Only one
%                                     instance is returned when set to
%                                     'none', independent of how many are
%                                     asked for.
%   'rngSeed'                       - Integer.  Set rng seed. Empty (default) means don't touch the
%                                     seed.
%
% You can pass more than one noise flag in the cell array of strings.  The
% response engine computes for each of these, and returns the responses for
% each noise flag in a Matlab container. Containers are very convenient but
% take a little getting used to.  We provide some usage comments where the
% neural response engine is called a bit further down.  We also explain
% there the format of the returned responses.
%
% The use of two example neural response engines is illustrated below.
%
whichNeuralEngine = 'nrePhotopigmentExcitationsCmosaic';
switch (whichNeuralEngine)
    case 'nrePhotopigmentExcitationsCmosaic'
        % Basic retinal image formation and sampling by the cone mosaic.
        % Note use of neural engine to get its own default parameters and
        % adjust them.  Smaller field of view speeds things up.
        neuralParams = nrePhotopigmentExcitationsCmosaic;
        neuralParams.coneMosaicParams.sizeDegs = [0.25 0.25];
        neuralParams.coneMosaicParams.timeIntegrationSeconds = 0.1;
        theNeuralEngine = neuralResponseEngine(@nrePhotopigmentExcitationsCmosaic,neuralParams);
        
        % The actual threshold varies enough with the different engines that we
        % need to adjust the contrast range that Quest+ searches over, as well as
        % the range of psychometric function slopes.  The threshold values
        % are specified as the negative log10(threshold), for reasons
        % beyond knowing.
        logThreshLimitLow = 4;
        logThreshLimitHigh = 1;
        logThreshLimitDelta = 0.05;
        slopeRangeLow = 1/20;
        slopeRangeHigh = 10000/20;
        slopeDelta = 5/20;
        
    case 'nreScenePhotonNoise'
        % Add Poisson noise to the photon counts.  This doesn't take any
        % parameters to speak of. In the early literature on ideal
        % observers, there was interest in comparing how well the human
        % visual system performed against what was imposed by physical
        % limits and no consideration of the visual system.  This 'neural'
        % engine instantiates calculations of this sort.
        theNeuralEngine = neuralResponseEngine(@nreScenePhotonNoise);
        
        % Again, custom ranges for the code that finds threshold
        logThreshLimitLow = 7;
        logThreshLimitHigh = 5;
        logThreshLimitDelta = 0.005;
        slopeRangeLow = 100/20;
        slopeRangeHigh = 10000/20;
        slopeDelta = 100/20;
        
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
%   'rcePoisson'
%   'rceTemplateDistanceTAFC'
%   'rcePcaSVMTAFC'
whichObserver = 'rcePcaSVMTAFC';
switch whichObserver
    case 'rcePoisson'
        % The ideal observer for a TAFC task limited by Poisson noise.
        % This classifier doesn't take any parameters.
        theRawClassifierEngine = responseClassifierEngine(@rcePoisson);
        
        % Below we use a wrapper routine to train the classifier and to
        % predict trial-by-trial responses.  Because that routine is
        % general, we define some parameters for it here.
        %
        % We'll illustrate a signal known exactly classifier, using the trainFlag value
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
    
    case 'rceTemplateDistanceTAFC'
        % Template matching (nearest neighbor) classifier for TAFC
        theRawClassifierEngine = responseClassifierEngine(@rceTemplateDistanceTAFC);
        
        % Below we use a wrapper routine to train the classifier and to
        % predict trial-by-trial responses.  Because that routine is
        % general, we define some parameters for it here.
        %
        % We'll illustrate a signal known exactly classifier, using the trainFlag value
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
%
% We have found that computing out psychometric functions takes a long time
% if you use the method of constant stimuli.  Thus we use an adaptive
% psychophysical procedure, QUEST+, to make the computations efficient.
% You can learn more about QUEST+ at the gitHub site that has our QUEST+
% Matlab implementation: https://github.com/BrainardLab/mQUESTPlus.
%
% A features of QUEST+ is that you need to specify a discrete list of
% contrast values that will be tested, which in the wrapper here is also
% taken as the set of possible thresholds as QUEST+ chooses stimulus values.
% We also discretize the possible slopes of the underlying Weibull
% psychometric function.  As noted above, appropriate choices for these
% depend on the neural engine (and possibly the quality of the response
% classifier) and need to be set on the basis of experience.  There is also
% a little art to choosing the spacing of the discrete values.  Too many
% and QUEST+ will run slowly; too few and it may not be able to put stimuli
% in the right places.
estDomain  = -logThreshLimitLow : logThreshLimitDelta : -logThreshLimitHigh;
slopeRange = slopeRangeLow: slopeDelta : slopeRangeHigh;

% There are two recommended ways to setup the QUEST+ threshold engine,
% which vary in terms of how clever they try to be about how many trials to
% run before stopping.
%
% The fixed number option runs a fixed number of trials.  This is fine if
% you have a good idea of how many trials you need to obtain the accuracy
% you want.
%
% In adaptive mode, we set up multiple QUEST+ instances and run them
% interleaved. We stop when the standard error of the threshold estimates
% across them becomes small enough.  See below for more.
%
% Note the explicit setting of the PF for the questThresholdEnging.  Using
% @qpPFWeibullLog causes it all to happen in log10 units, rather than the
% default dB units.
%
% Choices:
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
            'estDomain', estDomain, 'slopeRange', slopeRange, 'numEstimator', 1, ...
            'qpPF',@qpPFWeibullLog);
        
    case 'adaptiveMode'
        % Run 'numEstimator > 1' interleaved QUEST+ objects. In this case,
        % the threshold engine calculates the running standard error (SE)
        % among those objects. Essentially, this is an internal estimate of
        % the precision of the current threshold estimate.  Given this, the
        % stopping criterion is triggered when either of the following is
        % true
        %   1) total number of trials >= 'minTrial' AND the 'stopCrition'
        %   is true.
        %   2) when total number of trials >= 'maxTrial'.
        % The total number of trials includes all the trials across the
        % multiple QUEST+ objects.
        %
        % The 'stopCriterion' itself can be defined in two ways.  
        %   1) A single number. In this case, the criterion will be
        %      true when the SE is less than the specified number.
        %   2) A function handle that takes the current estimate of threshold
        %      and SE estimates as input arguments, and returns a boolean variable.
        %      For example: we can use a relative criterion w.r.t. the magnitude of threshold  
        %
        % Comment in either of the lines below to chose which you do. As
        % with the contrast range, there is some art to choosing the
        % specific values for the stop criterion, as well as the min and
        % max number of trials.  For example, the relative stopping
        % criterion in the function handle below can terminate too early if
        % initial threshold values are large.  This can be avoided by
        % appropriate choice of minimum number of trials.
        %
        % Choices (comment in one):
        %stopCriterion = 0.025;
        stopCriterion = @(threshold, se) se / abs(threshold) < 0.02;
        
        if (nTest > 1)
            blockedVal = true;
        else
            blockedVal = false;
        end

        % Set up the estimator object.
        estimator = questThresholdEngine('blocked',blockedVal,...
            'minTrial', 2e2, 'maxTrial', 1e3, ...
            'estDomain', estDomain, 'slopeRange', slopeRange, ...
            'numEstimator', 4, 'stopCriterion', stopCriterion, ...
            'qpPF',@qpPFWeibullLog);
        
    case 'validationMode'
        % Validation mode, instead of running an adpative procedure, compute the
        % full psychometric curve at each contrast level for nTest number
        % of trials        
        estimator = questThresholdEngine('validation', true, 'nRepeat', nTest, ...
            'estDomain', estDomain, 'slopeRange', slopeRange, ...
            'qpPF',@qpPFWeibullLog);
        
    otherwise
        error('Unknown threshold engine mode specified');
end

%% Generate the NULL scene sequence
%
% Threshold will be measured as a perturbation from this scene.  Typically
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
%
% There is some logic here to cache scenes and classifiers for each
% contrast tested, so that things will run faster if a contrast is repeated
% a second time.  This can use up space, so for real problems you may need
% to think about space-time tradeoffs, caching things to disk, etc.
%
% This logic is encapsulated in function computeThreshold, but we do it
% explicitly here to illustrate the basics of that function.

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
        %
        % See helper function computePerformance for how we use the
        % components to get trial by trial predictions.
        %
        % Note that because the classifer engine is a handle class, we
        % need to use a copy method to do the caching.  Otherwise the
        % cached pointer will simply continue to point to the same
        % object, and be updated by future training.
        [predictions, tempClassifierEngine] = computePerformance(...
            {theNullSceneSequence, theTestSceneSequences{testedIndex}},...
            theSceneTemporalSupportSeconds, nTrain, nTest, ...
            theNeuralEngine, theRawClassifierEngine, trainFlag, testFlag, ...
            'saveResponses',false, 'visualizeAllComponents', false, ...
            'TAFC',true);
        
        % Copy the trained classifier
        theTrainedClassifierEngines{testedIndex} = tempClassifierEngine.copy;
    else
        % Classifier is already trained, just get predictions
        predictions = computePerformance(...
            {theNullSceneSequence, theTestSceneSequences{testedIndex}},...
            theSceneTemporalSupportSeconds, nTrain, nTest, ...
            theNeuralEngine, theTrainedClassifierEngines{testedIndex},...
            [], testFlag, 'saveResponses',false, ...
            'visualizeAllComponents', false, ...
            'TAFC',true);
    end
    
    % Report what happened
    fprintf('Current test contrast: %g, P-correct: %g \n', testContrast, mean(predictions));
    
    % Tell QUEST+ what we ran (how many trials at the given contrast) and
    % get next stimulus contrast to run.
    if (estimator.validation)
        % Method of constant stimuli
        [logContrast, nextFlag] = ...
            estimator.multiTrial(logContrast * ones(1, nTest), predictions);
    else
        % Quest
        [logContrast, nextFlag] = ...
            estimator.multiTrialQuestBlocked(logContrast * ones(1, nTest), predictions);
    end

    % Get current threshold estimate
    [threshold, stderr] = estimator.thresholdEstimate();
    fprintf('Current threshold estimate: %g, stderr: %g \n', 10 ^ threshold, stderr);
end

%% Show results
fprintf('%d trials recorded \n', estimator.nTrial);

% Estimate threshold and plot/report results.  This
% does a maximumu likelihood based on the trials run, and is not subject to
% the discretization used by QUEST+.
if(strcmp(questMode, 'validationMode'))
    plotSize = 50;
else
    plotSize = 10;
end

% Return threshold value. For the mQUESTPlus Weibull PFs, the first
% parameter of the PF fit is the 0.81606 proportion correct threshold,
% when lapse rate is 0 and guess rate is 0.5.  Better to make this an
% explicit parameter, however.
figure();
thresholdCriterion = 0.81606;
[threshold, para] = estimator.thresholdMLE('showPlot', true, 'pointSize', plotSize, ...
    'thresholdCriterion', thresholdCriterion);
fprintf('Maximum likelihood fit parameters: %0.2f, %0.2f, %0.2f, %0.2f\n', ...
    para(1), para(2), para(3), para(4));
fprintf('Threshold (criterion proportion correct %0.4f: %0.2f (log10 units)\n', ...
    thresholdCriterion,threshold);

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
        
        pCorrect(idx) = mean(computePerformanceTAFC(...
            theNullSceneSequence, theTestSceneSequence, ...
            theSceneTemporalSupportSeconds, 512, 512, ...
            theNeuralEngine, theRawClassifierEngine, trainFlag, testFlag, ...
            false, false));
    end
    
    hold on;
    plot(logContrast, pCorrect, '-ok', 'LineWidth', 1);
    ylim([0, 1]);
    
end

