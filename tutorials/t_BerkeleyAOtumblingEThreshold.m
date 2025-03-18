function [logThreshold, logMAR, questObj, psychometricFunction, fittedPsychometricParams, ...
    trialByTrialStimulusAlternatives,trialByTrialPerformance] = ...
    t_BerkeleyAOtumblingEThreshold(options)
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
% See also t_BerkeleyAOtumblingESceneEngine.

% Examples:
%{

%}

%% Pick up optional arguments
arguments
    % Run with fast parameters overrides
    options.fastParams (1,1) logical = true;

    % Keep rng doing the same thing each time for validation
    options.rngSeed (1,1) double = 12;

    % Print out/plot  more diagnostics, or not
    options.verbose (1,1) logical = false;
    options.visualEsOnMosaic (1,1) logical = false;
    options.visualizeScene (1,1) logical = true;

    % Wavelength support
    options.wave (:,1) double = (500:5:870)';

    % Psychometric parameters
    options.letterSizesNumExamined = 9;
    options.nTest = 512;
    options.thresholdP =  0.781;

    % Optics parameters
    options.pupilDiameterMm (1,1) double = 6;
    options.accommodatedWl  (1,1) double = 840;
    options.defocusDiopters (1,1) double =  0.05;

    % Use cone contrast
    options.useConeContrast (1,1) logical = false

     % Choose noise model
    %   Choices: 'Poisson'
    %                  'Gaussian'
    options.whichNoisyInstanceNre (1,:) char = 'Poisson'
    options.gaussianSigma double = [];

     % Apply temporal filter?
    %
    % The timebase of the filter is assumed to match the frame rate, so we
    % only need to specify a list of filter values.  Since these can
    % describe gain changes as well as temporal processing per se, we do
    % not require that these sum to 1.  If you want to preserve signal
    % magnitude, you can normalize the filter values yourself, so that they
    % sum to 1. This can also be set to some string, e.g.,
    % 'photocurrentImpulseResponseBased', in which case the filter values
    % are computed on the fly
    %
    % Finally, this can be 'watsonFilter'
    options.temporalFilterValues (1,:) = []

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
    options.eHeightMin (1,1) double = 30;
    options.temporalModulationParams_frameRateHz (1,1) double = 60;
    options.temporalModulationParams_numFrame (1,1) double = 3;
    options.temporalModulationParams_xShiftPerFrame (1,:) double = [0 10/60 0];
    options.temporalModulationParams_yShiftPerFrame (1,:) double = [0 0 10/60];
    options.temporalModulationParams_backgroundRGBPerFrame (:,:) double = [0 0 0; 1 0 0; 0 0 0];

    % Run the validation check?  This gets overridden to empty if other
    % options change the conditions so that the validation data don't
    % apply.
    options.validationThresholds (1,:) double = []
end

%% Initialize figures
close all;

%% Stabilize rng
rng(options.rngSeed);

%% Fast parameter overrides
if (options.fastParams)
    options.displayNPixels = 128;
    options.letterSizesNumExamined = 4;
    options.nTest = 64;
end

%% Map options the way we need them below

% Scene parameters default overrides
aoSceneParams = struct( ...
    'visualizeScene', options.visualizeScene, ...
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
    'eHeightMin', options.eHeightMin' , ...
    'temporalModulationParams_frameRateHz', options.temporalModulationParams_frameRateHz , ...
    'temporalModulationParams_numFrame', options.temporalModulationParams_numFrame , ...
    'temporalModulationParams_xShiftPerFrame', options.temporalModulationParams_xShiftPerFrame , ...
    'temporalModulationParams_yShiftPerFrame', options.temporalModulationParams_yShiftPerFrame , ...
    'temporalModulationParams_backgroundRGBPerFrame', options.temporalModulationParams_backgroundRGBPerFrame , ...
    'wave', options.wave ...
    );

% Make sure figures and results directories exist so that output writes
% don't fail
rootPath = ISETBioCSFGeneratorRootPath;
if (~exist(fullfile(rootPath,'local',mfilename,'figures'),'dir'))
    mkdir(fullfile(rootPath,'local',mfilename,'figures'));
end
if (~exist(fullfile(rootPath,'local',mfilename,'results'),'dir'))
    mkdir(fullfile(rootPath,'local',mfilename,'results'));
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
[sce0,sce90,sce180,sce270,backgroundSceneEngine,sceneParams] = t_BerkeleyAOtumblingESceneEngine(sceneOptionsCell{:});
tumblingEsceneEngines = {sce0, sce90, sce180, sce270};
clear sce0 sce90 sce180 sce270

%% Set up temporal filter if we have one.
%
% Note that the nre returns the same number of frames it was passed.  So if
% you want post-stimulus responses produced by the temporal filter, you
% need to pad your input appropriately. That padding process is not
% illustrated in this tutorial script.
if (~isempty(options.temporalFilterValues))
    % Photocurrent filter computed and applied in nre
    if (ischar(options.temporalFilterValues) & strcmp(options.temporalFilterValues,'photocurrentImpulseResponseBased'))
        temporalFilter.temporalSupport = '';
        temporalFilter.filterValues = options.temporalFilterValues;
        
    elseif (ischar(options.temporalFilterValues) & strcmp(options.temporalFilterValues,'watsonFilter'))
        % Watson filter, computed here
        [~,watsonParams] = WatsonFilter([],[]);
        temporalFilter.temporalSupport = frameDurationSeconds*(0:framesNum-1);
        temporalFilter.filterValues = WatsonFilter(watsonParams,temporalFilter.temporalSupport);
        
    else
        % Filter explicitly passed
        temporalFilter.filterValues = options.temporalFilterValues;
        temporalFilter.temporalSupport = frameDurationSeconds*(0:framesNum-1);
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
noiseFreeResponseParams.verbose = options.verbose;

% Cone params
noiseFreeResponseParams.coneMosaicParams.wave = options.wave;
noiseFreeResponseParams.coneMosaicParams.fovDegs = options.displayFOVDeg;
noiseFreeResponseParams.coneMosaicParams.sizeDegs = [0.5 0.5];
noiseFreeResponseParams.coneMosaicParams.timeIntegrationSeconds = ...
    1/backgroundSceneEngine.sceneParams.temporalModulationParams.frameRateHz;

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
        noisyInstancesParams.sigma = gaussianSigma;
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
            'nTrain', 1, 'nTest', options.nTest)

    otherwise
        error('Unsupported rce specified')
end

%% Parameters for threshold estimation/quest engine
thresholdPara = struct(...
    'maxParamValue', 1, ...    
    'logThreshLimitLow', 2.0, ...              % minimum log10(normalized param value)
    'logThreshLimitHigh', 0.0, ...             % maximum log10(normalized param value)
    'logThreshLimitDelta', 0.01, ...
    'slopeRangeLow', 1/20, ...
    'slopeRangeHigh', 500/20, ...
    'slopeDelta', 2/20, ...
    'thresholdCriterion', options.thresholdP, ...
    'guessRate', 1/numel(tumblingEsceneEngines), ...
    'lapseRate', [0 0.02]);

% Parameters for Quest
questEnginePara = struct( ...
    'qpPF',@qpPFWeibullLog, ...
    'minTrial', options.nTest*options.letterSizesNumExamined, ...
    'maxTrial', options.nTest*options.letterSizesNumExamined, ...
    'numEstimator', 1, ...
    'stopCriterion', 0.05);

% Compute psychometric function for the 4AFC paradigm with the 4 E scenes
[logThreshold, questObj, psychometricFunction, fittedPsychometricParams, ...
    trialByTrialStimulusAlternatives,trialByTrialPerformance] = ...
    computeThreshold(tumblingEsceneEngines, theNeuralEngine, classifierEngine, ...
    classifierPara, thresholdPara, questEnginePara, ...
    'visualizeAllComponents', ~true, ...
    'verbose', true, ...
    'TAFC', false, 'useMetaContrast', false);
logMAR = log10(10.^logThreshold*60/5);
threshold = 10.^logThreshold;

%% Plot the derived psychometric function and other things. 
pdfFileName = [];
plotDerivedPsychometricFunction(questObj, threshold, fittedPsychometricParams, ...
    thresholdPara,pdfFileName, 'xRange', [0.02 0.2]);
if (options.visualEsOnMosaic)
    % This runs but I am not sure it is actually showing the stimulus.
    % Might have to do with the fact that the stimulus is at 840 nm.
    visualizeSimulationResults(questObj, threshold, fittedPsychometricParams, ...
        thresholdPara, tumblingEsceneEngines, backgroundSceneEngine, theNeuralEngine, ...
       pdfFileName);
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
if (~isempty(validationThresholds))
    if (any(abs(threshold-validationThresholds)./validationThresholds > validationTolerance))
        threshold
        validationThresholds
        error(sprintf('Do not replicate validation thresholds to %d%%. Check that parameters match, or for a bug.',round(100*validationTolerance)));
    else
        fprintf('Validation regression check passes\n');
    end
else
    threshold
    fprintf('No validation thresholds, validation regression check not run\n');
end
end
