function [logThreshold, logMAR, questObj, psychometricFunction, fittedPsychometricParams, thresholdPara, ...
    trialByTrialStimulusAlternatives, trialByTrialPerformance, trialByTrialWhichResponses] = ...
    t_BerkeleyAOTumblingEThreshold(options)
% Compute tumbling E threshold with AO optics
%
% This takes a number of key/value pairs that control its detailed
% conditions.  See comments under options block below.
%
% Note, in case you are tempted, that this can't be run with the meta
% contrast method because the parameter is stimulus size rather than
% contrast, and responses do not scale as a linear function of size the way
% they do with contrast.
%
% See also t_BerkeleyAOTumblingESceneEngine.

% Examples:
%{
    % Run with default parameters.  More examples are available in the
    % ISETBerkeleyAO project.  Those call into this tutorial function.
    t_BerkeleyAOTumblingEThreshold( ...
        'visualizeScene', false, ...
        'validationThresholds',[0.027]);
%}

%% Pick up optional arguments
arguments
    % Run with fast parameters overrides
    options.fastParams (1,1) logical = true;

    % Keep rng doing the same thing each time for validation.
    % 0 means don't touch the seed so we get random variation
    options.rngSeed (1,1) double = 12;

    % Print out/plot  more diagnostics, or not
    options.verbose (1,1) logical = false;
    options.visualizeEsOnMosaic (1,1) logical = false;
    options.visualizeEsWhichFrames (1,:) double = 1;
    options.visualizeEsFileBase (1,:) char = '';

    options.visualizeScene (1,1) logical = true;
    options.scenePdfFileBase (1,:) char = '';
    options.outputFiguresDir (1,:) char = '';
    options.outputResultsDir (1,:) char = '';

    options.plotPsychometric (1,1) logical = true;

    % Wavelength support
    options.wave (:,1) double = (500:5:870)';

    % Psychometric parameters
    options.letterSizesNumExamined = 9;
    options.nTest = 512;
    options.thresholdP = 0.781;
    options.minLetterSizeMinutes (1,1) double = 0.02;
    options.maxLetterSizeMinutes (1,1) double = 2;

    % Optics parameters
    options.pupilDiameterMm (1,1) double = 6;
    options.accommodatedWl  (1,1) double = 840;
    options.defocusDiopters (1,1) double =  0.05;

    % Use cone contrast
    options.useConeContrast (1,1) logical = false

     % Choose noise model
    %   Choices: 'Poisson'
    %            'Gaussian'
    options.whichNoisyInstanceNre (1,:) char = 'Poisson'
    options.gaussianSigma double = [];

     % Apply temporal filter?
    %
    % The timebase of the filter is assumed to match the frame rate, so we
    % only need to specify a list of filter values.  Since these can
    % describe gain changes as well as temporal processing per se, we do
    % not require that these sum to 1.  If you want to preserve signal
    % magnitude, you can normalize the filter values yourself, so that they
    % sum to 1.
    % 
    % This can also be set to 'photocurrentImpulseResponseBased', in which
    % case the filter values are computed on the fly
    %
    % It can also be 'watsonFilter'
    options.temporalFilterValues (1,:) = [];
    options.watsonParams_tau = 6.25;

    % Choose classifier engine
    %    rcePoisson - signal known exactly Poission max likelihood
    %    rceTemplateDistance - signal known exactly nearest L2 template
    %                 distance.
    options.whichClassifierEngine (1,:) char = 'rcePoisson'

    % AO scene parameters
    options.displayNPixels (1,1) double = 128;
    options.displayFOVDeg (1,1) double = 1.413;
    options.AOPrimaryWls (1,3) double = [840 683 543]; % [700 683 54];
    options.AOPrimaryFWHM (1,3) double = [22 27 23];
    options.AOCornealPowersUW (1,3) double = [141.4 10 10];
    options.plotDisplayCharacteristics (1,1) logical = false;
    options.chromaSpecification_type (1,:) char = 'RGBsettings';
    options.chromaSpecification_backgroundRGB (1,3) double = [1 0 0];
    options.chromaSpecification_foregroundRGB (1,3) double = [0 0 0];
    options.temporalModulationParams_frameRateHz (1,1) double = 60;
    options.temporalModulationParams_numFrame (1,1) double = 3;
    options.temporalModulationParams_xShiftPerFrameMin (1,:) double = [0 10 0];
    options.temporalModulationParams_yShiftPerFrameMin (1,:) double = [0 0 10];
    options.temporalModulationParams_backgroundRGBPerFrame (:,:) double = [0 0 0; 1 0 0; 0 0 0];
    options.temporalModulationParams_stimOnFrames (:,:) double = [0 1 0];

    % Run the validation check?  This gets overridden to empty if other
    % options change the conditions so that the validation data don't
    % apply.
    options.validationThresholds (1,:) double = []
end

%% Initialize figures
close all;

%% Stabilize rng
if (options.rngSeed ~= 0)
    rng(options.rngSeed);
end

%% Fast parameter overrides
if (options.fastParams)
    options.displayNPixels = 128;
    options.letterSizesNumExamined = 4;
    options.nTest = 64;
end

%% Map options the way we need them below

% Scene parameters default overrides
aoSceneParams = struct( ...
    'outputFiguresDir', options.outputFiguresDir;
    'visualizeScene', options.visualizeScene, ...
    'scenePdfFileBase', options.scenePdfFileBase, ...
    'pupilSizeMM', options.pupilDiameterMm, ...
    'displayNPixels', options.displayNPixels, ...
    'displayFOVDeg', options.displayFOVDeg, ...
    'AOPrimaryWls', options.AOPrimaryWls, ...
    'AOPrimaryFWHM', options.AOPrimaryFWHM, ...
    'AOCornealPowersUW', options.AOCornealPowersUW, ...
    'ambientSpd', zeros(size(options.wave)), ...
    'plotDisplayCharacteristics', options.plotDisplayCharacteristics, ...
    'chromaSpecification_type', options.chromaSpecification_type , ...
    'chromaSpecification_backgroundRGB', options.chromaSpecification_backgroundRGB , ...
    'chromaSpecification_foregroundRGB', options.chromaSpecification_foregroundRGB , ...
    'temporalModulationParams_frameRateHz', options.temporalModulationParams_frameRateHz , ...
    'temporalModulationParams_numFrame', options.temporalModulationParams_numFrame , ...
    'temporalModulationParams_xShiftPerFrame', options.temporalModulationParams_xShiftPerFrameMin/60 , ...
    'temporalModulationParams_yShiftPerFrame', options.temporalModulationParams_yShiftPerFrameMin/60 , ...
    'temporalModulationParams_backgroundRGBPerFrame', options.temporalModulationParams_backgroundRGBPerFrame , ...
    'temporalModulationParams_stimOnFrames', options.temporalModulationParams_stimOnFrames , ...
    'wave', options.wave ...
    );

% Make sure figures and results directories exist so that output writes
% don't fail
if (isempty(options.outputFiguresDir))
    outputFiguresDir = fullfile(ISETBioCSFGeneratorRootPath,'local',mfilename,'figures');
end
if (isempty(options.outputResultsDir))
    outputResultsDir = fullfile(ISETBioCSFGeneratorRootPath,'local',mfilename,'results');
end
if (~exist(outputFiguresDir,'dir'))
    mkdir(outputFiguresDir);
end
if (~exist(outputResultsDir,'dir'))
    mkdir(outputResultsDir);
end

% Define the AO scene parameters for the experiment we are modeling
%
% Get the tumbling E scene engines.
%
% At the moment cannot vary the set of orientations, but would be easy
% enough to allow this with a key/value pair passed to the called tutorial
% function.
%
% The scene engine tutorial returns its parameters, which are used below to
% try to match things up as best as possible.  One thing it doesn't return
% is its four hard coded orientaions, so we build that here for use below.
sceneOptionsCell = [fieldnames(aoSceneParams) , struct2cell(aoSceneParams)]';
[sce0,sce90,sce180,sce270,backgroundSceneEngine,sceneParams] = t_BerkeleyAOTumblingESceneEngine(sceneOptionsCell{:});
tumblingEsceneEngines = {sce0, sce90, sce180, sce270};
clear sce0 sce90 sce180 sce270

% For the background scene, we have to supply a letter size.  We don't care what
% what it is because the foreground and background colors are matched for
% the background scene.  But can't make it too big relative to the FOV.
backgroundSceneSequence = backgroundSceneEngine.compute(sceneParams.displayFOVDeg/5);

%% Set up temporal filter if we have one.
%
% Note that the nre returns the same number of frames it was passed.  So if
% you want post-stimulus responses produced by the temporal filter, you
% need to pad your input appropriately. That padding process is not
% illustrated in this tutorial script.
frameDurationSeconds = 1/backgroundSceneEngine.sceneParams.temporalModulationParams.frameRateHz;
if (~isempty(options.temporalFilterValues))
    % Photocurrent filter computed and applied in nre
    if (ischar(options.temporalFilterValues) & strcmp(options.temporalFilterValues,'photocurrentImpulseResponseBased'))
        temporalFilter.temporalSupport = '';
        temporalFilter.filterValues = options.temporalFilterValues;
        
    elseif (ischar(options.temporalFilterValues) & strcmp(options.temporalFilterValues,'watsonFilter'))
        % Watson filter, computed here
        [~,watsonParams] = WatsonFilter([],[]);
        watsonParams.tau = options.watsonParams_tau;
        temporalFilter.temporalSupport = frameDurationSeconds*(0:options.temporalModulationParams_numFrame-1);
        temporalFilter.filterValues = WatsonFilter(watsonParams,temporalFilter.temporalSupport);
        
    else
        % Filter explicitly passed
        temporalFilter.filterValues = options.temporalFilterValues;
        temporalFilter.temporalSupport = frameDurationSeconds*(0:options.temporalModulationParams_numFrame-1);
        if (length(temporalFilter.values) ~= length(temporalFilter.temporalSupport))
            error('Passed temporal filter values not matched in length to computed temporal filter support');
        end
    end
else
    temporalFilter = [];
end

%% Create neural response engine
%
% This calculates excitations in a patch of cone mosaic with Poisson
% noise, and includes optical blur.
noiseFreeResponseParams = nreNoiseFreeCMosaic([],[],[],[],'opticsType','BerkeleyAO');

% Set optics params
noiseFreeResponseParams.opticsParams.wls = options.wave;
noiseFreeResponseParams.opticsParams.pupilDiameterMM = options.pupilDiameterMm;
noiseFreeResponseParams.opticsParams.defocusAmount = options.defocusDiopters;
noiseFreeResponseParams.opticsParams.accommodatedWl = options.accommodatedWl;
noiseFreeResponseParams.opticsParams.zCoeffs = zeros(66,1);
noiseFreeResponseParams.opticsParams.defeatLCA = true;
noiseFreeResponseParams.nullStimulusSceneSequence = backgroundSceneSequence;
noiseFreeResponseParams.verbose = options.verbose;

% Cone params
noiseFreeResponseParams.coneMosaicParams.wave = options.wave;
noiseFreeResponseParams.coneMosaicParams.fovDegs = options.displayFOVDeg;
noiseFreeResponseParams.coneMosaicParams.sizeDegs = [0.5 0.5];
noiseFreeResponseParams.coneMosaicParams.timeIntegrationSeconds = ...
    frameDurationSeconds;

% Temporal filter
noiseFreeResponseParams.temporalFilter = temporalFilter;

% Handle cone contrast setting
if (options.useConeContrast)
    noiseFreeResponseParams.coneMosaicParams.outputSignalType = 'coneContrast';
else
    noiseFreeResponseParams.coneMosaicParams.outputSignalType = 'coneExcitations';
end

switch (options.whichNoisyInstanceNre)
    case 'Poisson'
        nreNoisyInstances = @nreNoisyInstancesPoisson;
        noisyInstancesParams = nreNoisyInstancesPoisson;

    case 'Gaussian'
        nreNoisyInstances = @nreNoisyInstancesGaussian;
        noisyInstancesParams = nreNoisyInstancesGaussian;
        noisyInstancesParams.sigma = options.gaussianSigma;
end

% Set up nre
theNeuralEngine = neuralResponseEngine( ...
    @nreNoiseFreeCMosaic, ...
    nreNoisyInstances, ...
    noiseFreeResponseParams, ...
    noisyInstancesParams);

%% N-way classifier.
switch (options.whichClassifierEngine)
    case {'rcePoisson'}
        classifierEngine = responseClassifierEngine(@rcePoisson);
        classifierPara = struct('trainFlag', 'none', ...
            'testFlag', 'random', ...
            'nTrain', 1, 'nTest', options.nTest);

    case {'rceTemplateDistance'}
        classifierEngine = responseClassifierEngine(@rceTemplateDistance);
        classifierPara = struct('trainFlag', 'none', ...
            'testFlag', 'random', ...
            'nTrain', 1, 'nTest', options.nTest);

    otherwise
        error('Unsupported rce specified');
end

%% Parameters for threshold estimation/quest engine
thresholdPara = struct(...
    'maxParamValue', 1, ...    
    'logThreshLimitLow', -1*log10(options.minLetterSizeMinutes/60), ...              % negative minimum log10(normalized param value) 
    'logThreshLimitHigh', -1*log10(options.maxLetterSizeMinutes/60), ...             % negative maximum log10(normalized param value)
    'logThreshLimitDelta', 0.01, ...
    'slopeRangeLow', 1/20, ...
    'slopeRangeHigh', 500/20, ...
    'slopeDelta', 2/20, ...
    'thresholdCriterion', options.thresholdP, ...
    'guessRate', 1/numel(tumblingEsceneEngines), ...
    'lapseRate', [0 0]);

% Parameters for Quest
questEnginePara = struct( ...
    'qpPF',@qpPFWeibullLog, ...
    'minTrial', options.nTest*options.letterSizesNumExamined, ...
    'maxTrial', options.nTest*options.letterSizesNumExamined, ...
    'numEstimator', 1, ...
    'stopCriterion', 0.05);

% Compute psychometric function for the 4AFC paradigm with the 4 E scenes
[logThreshold, questObj, psychometricFunction, fittedPsychometricParams, ...
    trialByTrialStimulusAlternatives,trialByTrialPerformance,trialByTrialWhichResponses] = ...
    computeThreshold(tumblingEsceneEngines, theNeuralEngine, classifierEngine, ...
    classifierPara, thresholdPara, questEnginePara, ...
    'visualizeAllComponents', ~true, ...
    'verbose', true, ...
    'TAFC', false, 'useMetaContrast', false);
logMAR = log10(10.^logThreshold*60/5);
threshold = 10.^logThreshold;

%% Plot the derived psychometric function and other things. 
if (options.plotPsychometric)
    pdfFileName = options.scenePdfFileBase;
    [stimulusLevels, pCorrect] = plotPsychometricFunction(questObj, threshold, fittedPsychometricParams, ...
        thresholdPara, pdfFileName, 'xRange', [options.minLetterSizeMinutes/60  options.maxLetterSizeMinutes/60]);

    % Print out table of stimulus levels and pCorrect
    fprintf('\nMeasured performance\n')
    for ii = 1:length(stimulusLevels)
        fprintf('%0.2f min (%0.3f deg), %0.2f pCorrect\n',60*stimulusLevels(ii),stimulusLevels(ii),pCorrect(ii));
    end
    fprintf('\n');
end
if (options.visualizeEsOnMosaic)
    % This runs but I am not sure it is actually showing the stimulus.
    % Might have to do with the fact that the stimulus is at 840 nm.
    for ff = 1:length(options.visualizeEsWhichFrames)
        if (~isempty(options.visualizeEsFileBase))
            pdfFileName = fullfile(outputFiguresDir, [options.visualizeEsFileBase '_VisualizeEsOnMosaic_Frame' num2str(options.visualizeEsWhichFrames(ff)) '.pdf']);
        else
            pdfFileName = [];
        end
        visualizeAOTumblingESimulationResults(questObj, threshold, fittedPsychometricParams, ...
            thresholdPara, options.visualizeEsWhichFrames(ff), tumblingEsceneEngines, backgroundSceneEngine, theNeuralEngine, ...
            pdfFileName);
    end
end

%% Do a check on the answer
%
% So that if we break something in the future we will have
% a chance of knowing it. The current numbers don't quite
% match the old version, but I think that is because of a change
% away from iePoisson which was not freezing the rng, and also
% other changes somewhere in stochasticity that I have not quite
% tracked down. But this validation generally passes.  Might fail
% sometimes.
validationTolerance = 0.01;
if (~isempty(options.validationThresholds))
    if (any(abs(threshold-options.validationThresholds)./options.validationThresholds > validationTolerance))
        fprintf('Threshold is %0.3f degs, %0.2f minutes; validation threshold is %0.3f degs\n',threshold,60*threshold, options.validationThresholds);
        error(sprintf('Do not replicate validation thresholds to %d%%. Check that parameters match, or for a bug.',round(100*validationTolerance)));
    else
        fprintf('Validation regression check passes\n');
        fprintf('Threshold is %0.3f degs, %0.2f minutes\n',threshold,60*threshold);
    end
else
    fprintf('No validation thresholds, validation regression check not run\n');
    fprintf('Threshold is %0.3f degs, %0.2f minutes\n',threshold,60*threshold);
end
end
