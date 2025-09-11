function thresholdRet = t_spatialCSF(options)
% Compute spatial CSF in different color directions, with fEM
%
% Syntax:
%   thresholdRet = t_spatialCSF;
%
% Description:
%    Use ISETBioCSFGenerator to run out CSFs in different color directions.
%    This example uses an ideal Poisson observer and circularly
%    windowed gratings of constant size, with fixational eye movements.
%
%    This is set up with key/value pairs that demonstate how to mix and
%    match noise free neural response models, noisy instance models, and
%    classifier engines. Different choices are illustrated in the examples
%    in the source code.
%
% Optional key/value pairs
%    See source code arguments block for a list of key/value pairs.

% See also: t_modulatedGratingsSceneGeneration,
%           t_chromaticThresholdContour, computeThreshold, computePerformance

% History:
%   10/20/20  lqz   Wrote it.
%   10/21/20  dhb   More commments.
%   10/22/20  lqz   Restructure the code
%   10/23/20  dhb   More commments.
%   10/25/20  dhb   Change contrast vectors to column vectors.  This is PTB
%                   convention, and also the convention of my brain.
%   05/10/23  fh    Edited it to call the new functions computeThreshold.m
%                       & computePerformance.m & rcePossion.m
%   04/17/24  dhb   Remove oldWay option.  Ever forward.  Enforce sine phase.
%   12/19/24  dhb   Update for new architecture.
%   07/30/25  NPC   Numerous updates for mRGCs.

% Note additional examples that use mRGC are in t_spatialCSF_mRGCExamples.
% Split out so they get run in smaller bunches.

% Examples:
%{
     % Validate without meta contrast
     t_spatialCSF('useMetaContrast', false, ...
        'whichNoiseFreeNre', 'excitationsCmosaic', ...
        'whichNoisyInstanceNre', 'Poisson', ...
        'whichClassifierEngine', 'rcePoisson', ...
        'useConeContrast', false, ...
        'useFixationalEMs', false, ...
        'temporalFilterValues', [], ...
        'oiPadMethod', 'mean', ...
        'validationThresholds', [0.0418    0.0783    0.1540    0.6759]);

     % Validate with meta contrast
     t_spatialCSF('useMetaContrast', true, ...
        'whichNoiseFreeNre', 'excitationsCmosaic', ...
        'whichNoisyInstanceNre', 'Poisson', ...
        'whichClassifierEngine', 'rcePoisson', ...
        'useConeContrast', false, ...
        'useFixationalEMs', false, ...
        'temporalFilterValues', [], ...
        'oiPadMethod', 'mean', ...
        'validationThresholds', [0.0418    0.0783    0.1540    0.6759]);

     % Validate with meta contrast, more frames
     t_spatialCSF('useMetaContrast', true, ...
        'whichNoiseFreeNre', 'excitationsCmosaic', ...
        'whichNoisyInstanceNre', 'Poisson', ...
        'whichClassifierEngine', 'rcePoisson', ...
        'useConeContrast', false, ...
        'useFixationalEMs', false, ...
        'numberOfFrames', 4, ...
        'temporalFilterValues', [], ...
        'oiPadMethod', 'mean', ...
        'validationThresholds',  [0.0208    0.0370    0.0770    0.3300]);
     
    % Verify that Gaussian noise works, as well as template classifier
    t_spatialCSF('useMetaContrast', true, ...
        'whichNoiseFreeNre', 'excitationsCmosaic', ...
        'whichNoisyInstanceNre', 'Gaussian', ...
        'whichClassifierEngine', 'rceTemplateDistance', ...
        'useConeContrast', false, ...
        'useFixationalEMs', false, ...
        'temporalFilterValues', [], ...
        'oiPadMethod', 'mean', ...
        'validationThresholds', [0.0997    0.1703    0.3581    1.0000]);

    % Verify that rcePcaSVM works
    t_spatialCSF('useMetaContrast', true, ...
        'whichNoiseFreeNre', 'excitationsCmosaic', ...
        'whichNoisyInstanceNre', 'Poisson', ...
        'whichClassifierEngine', 'rcePcaSVM', ...
        'nTrain', 128, ...
        'useConeContrast', false, ...
        'useFixationalEMs', false, ...
        'temporalFilterValues', [], ...
        'oiPadMethod', 'mean', ...
        'validationThresholds', [0.0377    0.0796    0.1427    0.6956]);

    % Verify that sceneAsresponses sce works
    t_spatialCSF('useMetaContrast', true, ...
        'whichNoiseFreeNre', 'sceneAsResponses', ...
        'whichNoisyInstanceNre', 'Poisson', ...
        'whichClassifierEngine', 'rcePoisson', ...
        'useConeContrast', false, ...
        'useFixationalEMs', false, ...
        'temporalFilterValues', [], ...
        'oiPadMethod', 'mean', ...
        'thresholdPara', struct( ...
            'logThreshLimitLow', 2.4, ...
            'logThreshLimitHigh', 0, ...
            'logThreshLimitDelta', 0.02, ...
            'slopeRangeLow', 1/20, ...
            'slopeRangeHigh', 100/20, ...
            'slopeDelta', 1/20, ...
            'thresholdCriterion', 0.81606), ...
        'validationThresholds', [0.3520    0.1895    0.1572    0.1554]);
%}

arguments
    % Use meta contrast method to speed things up?
    options.useMetaContrast (1,1) logical = true

    % Choose noise free neural model
    %   Choices: 'excitationsCmosaic'
    %            'sceneAsResponses'
    %            'mRGCMosaic'
    options.whichNoiseFreeNre (1,:) char  = 'excitationsCmosaic'

    % If the neural model is mRGCMosaic, can specify where to pick
    % off neural responses.
    %    Choices: 'mRGCs' (default). Use mRGC signals.
    %             'cones'.  Use the cone excitations (or contrast)
    options.mRGCOutputSignalType (1,:) char = 'mRGCs'

    % If the neural model is mRGCMosaic, we need to specify the
    % eccentricity and size of the prebaked mRGCmosaic (from which we can
    % crop a submosaic, depending on the stimulus size
    options.mRGCMosaicRawEccDegs (1,2) double = [0.0 0.0];
    options.mRGCMosaicRawSizeDegs (1,2) double = [2.0 2.0];

    % If the neural model is mRGCMosaic, we can specify the Optics type to use
    options.opticsType (1,:) char  = 'loadComputeReadyRGCMosaic';

    % Choose noise model
    %   Choices: 'Poisson'
    %                  'Gaussian'
    options.whichNoisyInstanceNre (1,:) char = 'Poisson'
    options.gaussianSigma double = [];

    % Choose classifier engine
    %    rcePoisson - signal known exactly Poission max likelihood
    %    rceTemplateDistance - signal known exactly nearest L2 template
    %                 distance.
    %    rcePcaSVM  - support vector machine linear classifier after PCA.
    options.whichClassifierEngine (1,:) char = 'rcePoisson'

    % Use cone contrast rather than cone excitations
    options.useConeContrast (1,1) logical = false

    % Use fixational eye movements, and some aspects of that.
    %
    %  useFixationalEMs: Controls whether we simulate fixational EMs.
    %  testEMsMatchTrain: Cetermines whether the training and test EM paths match
    %    each other. It can only be true if the number of training and test
    %    EMs match each other below. This is basically whether EM are known
    %    exactly or not known.
    %  sameEMsEachSF: Cetermines whether the same EM paths are used for each
    %    spatial frequency. This is a science decision.  Probably makes sense
    %    for them to be different but if one is studying particular paths one
    %    might want to set this to true and also possibly take more control
    %    of the particular paths than we do in this tutorial.
    % sameEMsEachContrast: The EM paths used for each test contrast.
    %    Conceptually one might want this false, but at present the code only
    %    allows true and an error is thrown if it is false. It is added to
    %    remind us at the top level of what the underlying code is doing.
    %    By assuming it is true, we could in principle use the meta contrast
    %    machinery for the EM case. If it is false that would not be
    %    possible in principle. To make it false, we would have to push the
    %    EM control code down further into computeThreshold.
    % seedForEMs: Seed for EM path generation
    options.useFixationalEMs (1,1) logical = false
    options.testEMsMatchTrain (1,1) logical = true
    options.sameEMsEachSF (1,1) logical = false
    options.sameEMsEachContrast (1,1) = true
    options.seedForEMs (1,1) double = 1021;

    % nTrainingEMs, nTestEMs
    %
    % If EMs are turned on, nTrainEMs determines the number of
    % fixational EMs used in training, while nTestEMs determines the number
    % used at test.
    options.nTrainEMs (1,1) double = 1
    options.nTestEMs (1,1) double = 1

    % nTrain and nTest
    %
    % nTrain determines the number of training instances. When there are
    % training fixational EMs, this is the number that get used to train the
    % classiefier for each EM.
    %
    % nTest determines the number of test instances. When there are
    % fixationsal EMs, this number is spread across them equally.  Thus
    % this number must be a multiple of the number of fixational EMs
    % specified.
    options.nTrain (1,1) double = 1
    options.nTest (1,1) double = 128

    % Mosaic size and eccentricity
    options.mosaicEccDegs (1,2) double = [0 0];
    options.mosaicSizeDegs (1,2) double = [0.5 0.5];

    % Varied stimulus parameter
    options.spatialFreqs (1,:) double = [4, 8, 16, 32];
    
    % Fixed stimulus parameters
    options.employMosaicSpecificConeFundamentals (1,1) logical = false;
    options.meanLuminanceCdPerM2 (1,1) double = 100;
    options.meanChromaticityXY (1,2) double = [0.30 0.32];
    options.stimulusChroma (1,:) char = 'luminance'

    options.orientationDegs (1,1) double = 90
    options.spatialPhaseDegs (1,1) double = 90
    options.stimOnFrameIndices (1,:) double = [];
    options.numberOfFrames double = []
    options.frameDurationSeconds (1,1) double = 0.1;
    options.stimDurationTemporalCycles (1,1) double = 1
    options.temporalFrequencyHz (1,1) double = 5;
    

    % Stimulus size
    options.stimSizeDegs (1,1) double = 0.5;
    options.pixelsNum (1,1) double = 128;

    % Presentation mode
    options.presentationMode (1,:) char = 'static';

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

    % Run the validation check?  This gets overridden to empty if other
    % options change the conditions so that the validation data don't
    % apply.
    options.validationThresholds (1,:) double = []

    % Fast parameters?
    options.fastParameters (1,1) logical = true;

    % Optical image pad method.
    % Determines how image is padded for
    % convolution with optical point spread function.
    %
    % Can be 'zero' or 'mean'.
    options.oiPadMethod (1,:) char = 'zero'

    % Apply a filter to the spectra before computing responses?  See
    % t_spatialCSFFilter
    options.filter (1,1) = struct('spectralSupport',[],'transmission',[])

    % Threshold limits
    options.thresholdPara (1,1) = struct( ...
        'logThreshLimitLow', 2.4, ...
        'logThreshLimitHigh', 0, ...
        'logThreshLimitDelta', 0.02, ...
        'slopeRangeLow', 1/20, ...
        'slopeRangeHigh', 50/20, ...
        'slopeDelta', 2.5/50, ...
        'thresholdCriterion', 0.81606, ...
        'guessRate', 1/2, ...
        'lapseRate', 0);

    % Number of psychometric curve samples
    options.psychometricCurveSamplesNum (1,:) double = []

    % Verbose?
    options.verbose (1,1) logical = true

    % Visualize the output of each compute (useful for debuging code)
    % Scene and response visualization are controlled separately.
    options.visualizeEachScene (1,1) logical = false
    options.visualizeEachResponse (1,1) logical = false
    options.responseVideoFileName(1,:) char = '';
    options.responseVisualizationFunction = []
    options.maxVisualizedNoisyResponseInstances = 1
    options.maxVisualizedNoisyResponseInstanceStimuli = 1
    options.figureFileBase (1,:) char = [];
    options.resultsFileBase (1,:) char = [];
    
    % Computed thresholds filename
    options.thresholdsDataFileName (1,:) char = '';

    options.parPoolSize (1,:) char = 'default'

end % arguments


%% Set flags from key/value pairs
filter = options.filter;
useMetaContrast = options.useMetaContrast;
useConeContrast = options.useConeContrast;

useFixationalEMs = options.useFixationalEMs;
testEMsMatchTrain = options.testEMsMatchTrain;
sameEMsEachSF = options.sameEMsEachSF;
sameEMsEachContrast = options.sameEMsEachContrast;
seedForEMs = options.seedForEMs;
nTrainEMs = options.nTrainEMs;
nTestEMs = options.nTestEMs;

nTrain = options.nTrain;
nTest = options.nTest;
whichNoiseFreeNre = options.whichNoiseFreeNre;
whichNoisyInstanceNre = options.whichNoisyInstanceNre;

gaussianSigma = options.gaussianSigma;
whichClassifierEngine = options.whichClassifierEngine;
validationThresholds = options.validationThresholds;
mRGCOutputSignalType = options.mRGCOutputSignalType;

% Mosaic sizes
mosaicEccDegs = options.mosaicEccDegs;
mosaicSizeDegs = options.mosaicSizeDegs;
mRGCRawEccDegs = options.mRGCMosaicRawEccDegs;
mRGCRawSizeDegs = options.mRGCMosaicRawSizeDegs;
mRGCCropSize = mosaicSizeDegs;

% Stimulus params
employMosaicSpecificConeFundamentals = options.employMosaicSpecificConeFundamentals;
meanLuminanceCdPerM2 = options.meanLuminanceCdPerM2;
meanChromaticityXY = options.meanChromaticityXY;
orientationDegs = options.orientationDegs;
spatialPhaseDegs = options.spatialPhaseDegs;

% Temporal params
temporalFrequencyHz = options.temporalFrequencyHz;
stimOnFrameIndices = options.stimOnFrameIndices;
temporalFilterValues = options.temporalFilterValues;
numberOfFrames = options.numberOfFrames;
frameDurationSeconds = options.frameDurationSeconds;
stimDurationTemporalCycles = options.stimDurationTemporalCycles;
presentationMode = options.presentationMode;

% Size params
stimSizeDegs = options.stimSizeDegs;
pixelsNum = options.pixelsNum;

fastParameters = options.fastParameters;
oiPadMethod = options.oiPadMethod;
opticsType = options.opticsType;
thresholdPara = options.thresholdPara;

verbose = options.verbose;
visualizeEachScene = options.visualizeEachScene;
visualizeEachResponse = options.visualizeEachResponse;

responseVisualizationFunction = options.responseVisualizationFunction;
maxVisualizedNoisyResponseInstances = options.maxVisualizedNoisyResponseInstances;
maxVisualizedNoisyResponseInstanceStimuli = options.maxVisualizedNoisyResponseInstanceStimuli;

thresholdsDataFileName = options.thresholdsDataFileName;
parPoolSize = options.parPoolSize;

%% Freeze rng for replicability and validation
rng(1);

%% Close any stray figs
close all;

%% Make sure local/figures directory exists so we can write out our figures in peace
projectBaseDir = ISETBioCSFGeneratorRootPath;
if (~exist(fullfile(projectBaseDir,'local',mfilename,'figures'),'dir'))
    mkdir(fullfile(projectBaseDir,'local',mfilename,'figures'));
end
if (~exist(fullfile(projectBaseDir,'local',mfilename,'results'),'dir'))
    mkdir(fullfile(projectBaseDir,'local',mfilename,'results'));
end

%% Figure output base name
if (isempty(options.figureFileBase))
    figureFileBase = fullfile(projectBaseDir,'local',mfilename,'figures', ...
        sprintf('%s_Meta_%d_ConeContrast_%d_FEMs_%d_%s_%s_%s_%s', mfilename, ...
        useMetaContrast,useConeContrast,useFixationalEMs,whichNoiseFreeNre,whichNoisyInstanceNre,...
        whichClassifierEngine,mRGCOutputSignalType));
else
    figureFileBase = fullfile(options.figureFileBase, ...
    sprintf('%s_Meta_%d_ConeContrast_%d_FEMs_%d_%s_%s_%s_%s', mfilename, ...
    useMetaContrast,useConeContrast,useFixationalEMs,whichNoiseFreeNre,whichNoisyInstanceNre,...
    whichClassifierEngine,mRGCOutputSignalType));
end

if (~exist(figureFileBase, 'dir'))
    mkdir(figureFileBase);
    fprintf('Generated figure sub-directory at %s\n', figureFileBase);
end

if (isempty(options.resultsFileBase))
    resultsFileBase = fullfile(projectBaseDir,'local',mfilename,'results', ...
        sprintf('%s_Meta_%d_ConeContrast_%d_FEMs_%d_%s_%s_%s_%s', mfilename, ...
        useMetaContrast,useConeContrast,useFixationalEMs,whichNoiseFreeNre,whichNoisyInstanceNre,...
        whichClassifierEngine,mRGCOutputSignalType));
else
    resultsFileBase = fullfile(options.resultsFileBase, ...
        sprintf('%s_Meta_%d_ConeContrast_%d_FEMs_%d_%s_%s_%s_%s', mfilename, ...
        useMetaContrast,useConeContrast,useFixationalEMs,whichNoiseFreeNre,whichNoisyInstanceNre,...
        whichClassifierEngine,mRGCOutputSignalType));
end

if (~exist(resultsFileBase, 'dir'))
        mkdir(resultsFileBase);
        fprintf('Generated figure sub-directory at %s\n', resultsFileBase);
end

if (~isempty(options.responseVideoFileName)) && (ischar(options.responseVideoFileName))
    options.responseVideoFileName = fullfile(options.figureFileBase, options.responseVideoFileName);
end


%% Parpool
AppleSiliconParPoolManager(parPoolSize);

%% Set base values for parameters that control what we do
if (~isempty(numberOfFrames))
    framesNum = numberOfFrames;
    padFramesBefore = 0;
    padFramesAfter = 0;
elseif (~useFixationalEMs && isempty(temporalFilterValues))
    framesNum = 1;
    padFramesBefore = 0;
    padFramesAfter = 0;
else
    framesNum = 4;
    padFramesBefore = 0;
    padFramesAfter = 0;
end


%% Set up temporal filter if we have one.
%
% Note that the nre returns the same number of frames it was passed.  So if
% you want post-stimulus responses produced by the temporal filter, you
% need to pad your input appropriately. That padding process is not
% illustrated in this tutorial script.
if (~isempty(temporalFilterValues))
    % Photocurrent filter computed and applied in nre
    if (ischar(temporalFilterValues) & strcmp(temporalFilterValues,'photocurrentImpulseResponseBased'))
        temporalFilter.temporalSupport = '';
        temporalFilter.filterValues = temporalFilterValues;
        
    elseif (ischar(temporalFilterValues) & strcmp(temporalFilterValues,'watsonFilter'))
        % Watson filter, computed here
        [~,watsonParams] = WatsonFilter([],[]);
        watsonParams.tau = options.watsonParams_tau;
        temporalFilter.temporalSupport = frameDurationSeconds*(0:framesNum-1);
        temporalFilter.filterValues = WatsonFilter(watsonParams,temporalFilter.temporalSupport);

    else
        % Filter explicitly passed
        temporalFilter.filterValues = temporalFilterValues;
        temporalFilter.temporalSupport = frameDurationSeconds*(0:framesNum-1);
    end
else
    temporalFilter = [];
end


%% List of spatial frequencies to be tested.
spatialFreqs = options.spatialFreqs;

%% Chromatic direction and contrast
stimulusChroma = options.stimulusChroma;

% Choose stimulus chromatic direction specified as a 1-by-3 vector
% of L, M, S cone contrast.  These vectors get normalized below, so only
% their direction matters in the specification.
switch (stimulusChroma)
    case 'luminance'
        chromaDir = [1.0, 1.0, 1.0]';
    case 'red-green'
        chromaDir = [1.0, -1.0, 0.0]';
    case 'L-isolating'
        chromaDir = [1.0, 0.0, 0.0]';
end

% Set the RMS cone contrast of the stimulus. Things may go badly if you
% exceed the gamut of the monitor, so we are conservative and set this at a
% value that is within gamut of typical monitors and don't worry about it
% further for this tutorial.  A vector length contrast of 0.1 should be OK.
rmsContrast = 0.1;
chromaDir = chromaDir / norm(chromaDir) * rmsContrast;
assert(abs(norm(chromaDir) - rmsContrast) <= 1e-10);


%% Set grating engine parameters

if (strcmp(presentationMode, 'static'))
    % DHB: I DON'T QUITE UNDERSTAND THE TWO CASES HERE.  WHY DON'T WE JUST
    % SPECIFY STATIC OR DYNAMIC AS AN OPTION?
    stimulusDuration = framesNum*frameDurationSeconds;
    
    if (~useFixationalEMs && isempty(temporalFilter) && framesNum == 1)
        gratingSceneParams = struct( ...
            'meanLuminanceCdPerM2', meanLuminanceCdPerM2, ...
            'meanChromaticityXY', meanChromaticityXY, ...
            'fovDegs', stimSizeDegs, ...
            'presentationMode', 'flashedmultiframe', ...
            'duration', stimulusDuration, ...
            'frameDurationSeconds', stimulusDuration/framesNum, ...
            'orientation', orientationDegs, ...
            'spatialPhase', spatialPhaseDegs, ...
            'pixelsNum', pixelsNum, ...
            'filter', filter...
            );
    else
        % Dynamic stimulus parameters
        presentationMode = 'counterphasemodulated';
        gratingSceneParams = struct( ...
            'meanLuminanceCdPerM2', meanLuminanceCdPerM2, ...
            'meanChromaticityXY', meanChromaticityXY, ...
            'fovDegs', stimSizeDegs, ...
            'presentationMode', presentationMode, ...
            'duration', stimulusDuration, ...
            'frameDurationSeconds', stimulusDuration/framesNum, ...
            'temporalFrequencyHz', temporalFrequencyHz, ...
            'stimOnFrameIndices', stimOnFrameIndices, ...
            'orientation', orientationDegs, ...
            'spatialPhase', spatialPhaseDegs, ...
            'pixelsNum', pixelsNum, ...
            'filter', filter ...
            );
    end

else
    % Non-static 
    stimulusDuration = stimDurationTemporalCycles*frameDurationSeconds;

    temporalModulationParams.stimOnFrameIndices = [];
    temporalModulationParams.stimDurationFramesNum = [];
    temporalModulationParams.phaseDirection = 1;
    temporalModulationParams.stimDurationTemporalCycles = stimDurationTemporalCycles;
    temporalModulationParams.temporalFrequencyHz = temporalFrequencyHz;

    gratingSceneParams = struct( ...
        'meanLuminanceCdPerM2', meanLuminanceCdPerM2, ...
        'meanChromaticityXY', meanChromaticityXY, ...
        'spectralSupport', 400:20:750, ...
        'fovDegs', stimSizeDegs, ...
        'pixelsNum', pixelsNum, ...
        'spatialEnvelope', 'rect', ...
        'spatialEnvelopeRadiusDegs', stimSizeDegs, ...
        'orientation', orientationDegs, ...
        'spatialPhase', spatialPhaseDegs, ...
        'presentationMode', presentationMode, ...
        'duration', stimulusDuration, ...
        'frameDurationSeconds', frameDurationSeconds, ...
        'temporalModulationParams', temporalModulationParams);
end


%% Thresholds filename
if (isempty(thresholdsDataFileName))
    thresholdsDataFileName = ...
        sprintf('%sCSF_%s_Optics_%s_EccDegs_x%2.1f_%2.1f_SizeDegs_%2.1fx%2.1f_OriDegs_%2.0f_%s.mat', ...
        mRGCOutputSignalType, ...
        stimulusChroma, ...
        opticsType, ...
        mosaicEccDegs(1), mosaicEccDegs(2), ...
        mosaicSizeDegs(1),mosaicSizeDegs(2), ...
        orientationDegs, ...
        presentationMode);
end


%% Neural response engines
%
% Setup the noise-free neural response engine
switch (whichNoiseFreeNre)

    case 'mRGCMosaic'
        % Set the compute function
        nreNoiseFreeComputeFunction = @nreNoiseFreeMidgetRGCMosaic;

        % Get default params struct
        nreNoiseFreeParams = nreNoiseFreeMidgetRGCMosaic([],[],[],[], ...
            'oiPadMethod',oiPadMethod);

        % Modify certain parameters of interest
        %
        % 1. Select one of the pre-computed mRGC mosaics by specifying its
        % eccentricity, size, and type.
        if (mosaicSizeDegs(1) > mRGCRawSizeDegs(1)) || (mosaicSizeDegs(2) > mRGCRawSizeDegs(2))
            mosaicSizeDegs
            mRGCRawSizeDegs
            error('Cannot ask for mosaic larger than mRGCRawSizeDegs');
        end
        nreNoiseFreeParams.mRGCMosaicParams.eccDegs = mRGCRawEccDegs;
        nreNoiseFreeParams.mRGCMosaicParams.sizeDegs = mRGCRawSizeDegs;
        nreNoiseFreeParams.mRGCMosaicParams.rgcType = 'ONcenterMidgetRGC';

        % Adjust surround optimization substring
        if (nreNoiseFreeParams.mRGCMosaicParams.eccDegs(1) < -20)
            nreNoiseFreeParams.mRGCMosaicParams.surroundOptimizationSubString = strrep(...
                nreNoiseFreeParams.mRGCMosaicParams.surroundOptimizationSubString, ...
                'PackerDacey2002H1freeLowH1params', 'PackerDacey2002H1freeUpperH1params');
        elseif (nreNoiseFreeParams.mRGCMosaicParams.eccDegs(1) < -16)
            nreNoiseFreeParams.mRGCMosaicParams.surroundOptimizationSubString = strrep(...
                nreNoiseFreeParams.mRGCMosaicParams.surroundOptimizationSubString, ...
                'PackerDacey2002H1freeLowH1params', 'PackerDacey2002H1freeMidH1params');
        end

        % 2. We can crop the mRGCmosaic to some desired size.
        %    Passing [] for sizeDegs will not crop.
        %    Passing [] for eccentricityDegs will crop the mosaic at its center.
        nreNoiseFreeParams.mRGCMosaicParams.cropParams = struct(...
            'sizeDegs', mRGCCropSize, ...
            'eccentricityDegs', mosaicEccDegs ...
            );

        % 3. Set the input cone mosaic integration time
        nreNoiseFreeParams.mRGCMosaicParams.coneIntegrationTimeSeconds = frameDurationSeconds;

        % 4. Set the noise level.
        %
        % Post-cone summation noise is additive Gaussian noise with a desired
        % sigma. When the input is raw cone excitations, the sigma should be
        % expressed in terms of cone excitations/integration time. When the input
        % to mRGCs is cone modulations with respect to the background, which have a
        % max amplitude of 1.0, the sigma should be scaled appropriately.
        % So if your mRGC mosiac were operating directly on cone excitations, a
        % different value more like 100 or 1000 would be reasonable, depending on
        % light level.
        nreNoiseFreeParams.mRGCMosaicParams.outputSignalType = mRGCOutputSignalType;
        if (isempty(gaussianSigma))
            switch (mRGCOutputSignalType)
                case 'mRGCs'
                    gaussianSigma = 0.003;
                case 'cones'
                    gaussianSigma = 50;
                otherwise
                    error('Unknown mRGC output signal type specified');
            end
        end

        % 5. Optics type
        nreNoiseFreeParams.opticsParams.type = opticsType;

        % Handle cone contrast setting
        if (useConeContrast)
            nreNoiseFreeParams.mRGCMosaicParams.inputSignalType = 'coneContrast';
        else
            nreNoiseFreeParams.mRGCMosaicParams.inputSignalType = 'coneExcitations';
        end


        % Handle temporal filter
        nreNoiseFreeParams.temporalFilter = temporalFilter;

        % Sanity check on the amount of mRGCMosaicVMembraneGaussianNoiseSigma for
        % the specified noiseFreeParams.mRGCMosaicParams.inputSignalType
        switch (nreNoiseFreeParams.mRGCMosaicParams.inputSignalType)
            case 'coneContrast'
                % Cone contrast is in the range of [-1 1]
                if (gaussianSigma > 1)
                    error('Gaussian noise sigma (%f) seems too large when operating on ''%s''.', ...
                        gaussianSigma, ...
                        nreNoiseFreeParams.mRGCMosaicParams.inputSignalType);
                end

            case 'coneExcitations'
                % Excitations are unlikely to be of order 1.
                if (gaussianSigma < 1)
                    error('Gaussian noise sigma (%f) seems too small when operating on ''%s''.', ...
                        gaussianSigma, ...I
                        nreNoiseFreeParams.mRGCMosaicParams.inputSignalType);
                end

            otherwise
                error('Unknown mRGC signal type specified: ''%s''.', nreNoiseFreeParams.mRGCMosaicParams.inputSignalType)
        end

    case 'excitationsCmosaic'
        % Set the compute function
        nreNoiseFreeComputeFunction = @nreNoiseFreeCMosaic;

        % Get default params struct
        nreNoiseFreeParams = nreNoiseFreeCMosaic([],[],[],[], ...
            'oiPadMethod',oiPadMethod);

        % Modify certain parameters of interest
        nreNoiseFreeParams.coneMosaicParams.sizeDegs = mosaicSizeDegs;
        nreNoiseFreeParams.coneMosaicParams.eccDegs = mosaicEccDegs;

        % Set the cone mosaic integration time to match the stimulus frame duration
        nreNoiseFreeParams.coneMosaicParams.timeIntegrationSeconds = frameDurationSeconds;

        % Handle cone contrast setting
        if (useConeContrast)
            nreNoiseFreeParams.coneMosaicParams.outputSignalType = 'coneContrast';
        else
            nreNoiseFreeParams.coneMosaicParams.outputSignalType = 'coneExcitations';
        end

        % Handle temporal filter
        nreNoiseFreeParams.temporalFilter = temporalFilter;

        % Set Gaussian sigma in case we're using Gaussian noise below.
        if (isempty(gaussianSigma))
            gaussianSigma = 50;
        end

    case 'sceneAsResponses'
        % This is image photon counts, and thus provides an estimate of the
        % upper bound on performance for a photon limited system
        %
        % Override size and frame duration to put performance in reasonable range
        stimSizeDegs = 0.05;
        frameDurationSeconds = 1e-16;

        % Set up nre
        nreNoiseFreeComputeFunction = @nreNoiseFreeSceneAsResponses;
        nreNoiseFreeParams = nreNoiseFreeSceneAsResponses;
        nreNoiseFreeParams.frameDurationSeconds = frameDurationSeconds;

        % This scene engine should not use cone contrast, no matter what
        useConeContrast = false;

        % Probably wouldn't use Gaussian noise with this case, but set the
        % sigma in case someone tries it.  No idea whether this is a good
        % value.  Probably not since I just made it up.
        if (isempty(gaussianSigma))
            gaussianSigma = 50;
        end

    otherwise
        error('Unsupported noise free nre specified');
end

% Setup the noisy neural response engine
switch (whichNoisyInstanceNre)
    case 'Poisson'
        nreNoisyInstancesComputeFunction = @nreNoisyInstancesPoisson;
        nreNoisyInstancesParams = nreNoisyInstancesPoisson;

    case 'Gaussian'
        nreNoisyInstancesComputeFunction = @nreNoisyInstancesGaussian;
        nreNoisyInstancesParams = nreNoisyInstancesGaussian;
        nreNoisyInstancesParams.sigma = gaussianSigma;

    otherwise
        error('Unsupported noisy instances nre specified');
end % switch (whichNoisyInstanceNre)



if (employMosaicSpecificConeFundamentals) 
   % For this computation we need the input cone mosaic and the optics that
   % will be used for the rest of the computation

   
   switch (whichNoiseFreeNre)

       case 'mRGCMosaic'
           mosaicParams = nreNoiseFreeParams.mRGCMosaicParams;
       case 'excitationsCmosaic'
           mosaicParams = nreNoiseFreeParams.coneMosaicParams;
       otherwise
           error('noiseFreeNRE must be either ''mRGCMosaic'', or ''excitationsCmosaic''.');
   end % switch (whichNoiseFreeNre)

   [theOI,theMRGCMosaic] = generateOpticsAndMosaicFromParams(...
            nreNoiseFreeParams.opticsParams, ...
            [], ...
            mosaicParams);

    % Compute the custom cone fundamentals and store them in gratingSceneParams
    % so in the next pass, we can generate the desired scene
    maxConesNumForAveraging = 3;
    gratingSceneParams.customConeFundamentals = visualStimulusGenerator.coneFundamentalsForPositionWithinConeMosaic(...
            theMRGCMosaic.inputConeMosaic, theOI, ...
            mosaicEccDegs, gratingSceneParams.fovDegs, maxConesNumForAveraging);
end



%% If we use cone contrast, we will neeed a null scene for normalization.
%
% If we use a filter that is computed on the fly, e.g. 'photocurrentImpulseResponseBased'
% we also need a null scene to compute the photocurrent given the background
%  So, we create a nullGratingSceneEngine so we can generate the theNullStimulusScene
if (useConeContrast) || (ischar(temporalFilterValues))
    dummySpatialFrequency = 4;
    nullGratingSceneEngine = createGratingSceneEngine(chromaDir, dummySpatialFrequency, ...
        gratingSceneParams);
    nreNoiseFreeParams.nullStimulusSceneSequence = nullGratingSceneEngine.compute(0.0);
end

%% Create the neural engine
theNeuralEngine = neuralResponseEngine( ...
    nreNoiseFreeComputeFunction, ...
    nreNoisyInstancesComputeFunction, ...
    nreNoiseFreeParams, ...
    nreNoisyInstancesParams ...
    );

% Set the neuralEngine's various visualization properties
theNeuralEngine.visualizeEachCompute = visualizeEachResponse;
theNeuralEngine.responseVideoFileName = options.responseVideoFileName;
theNeuralEngine.customVisualizationFunctionHandle = responseVisualizationFunction;
theNeuralEngine.maxVisualizedNoisyResponseInstances = maxVisualizedNoisyResponseInstances;


%% If using meta contrast, set this up.
if (useMetaContrast)
    metaSceneEngineParams = sceMetaContrast;
    theMetaSceneEngine = sceneEngine(@sceMetaContrast,metaSceneEngineParams);

    % Create nreMetaContrast using the actual scene and neural engines
    metaNeuralResponseEngineNoiseFreeParams = nreNoiseFreeMetaContrast;
    metaNeuralResponseEngineNoiseFreeParams.contrast0 = 0;
    metaNeuralResponseEngineNoiseFreeParams.contrast1 = 1;
    metaNeuralResponseEngineNoiseFreeParams.neuralEngine = theNeuralEngine;

    metaNeuralResponseEngineNoisyInstanceParams = nreNoisyInstancesMetaContrast;
    metaNeuralResponseEngineNoisyInstanceParams.neuralEngine = theNeuralEngine;
end


%% Instantiate the responseClassifierEngine
%
% rcePoisson makes decision by performing the Poisson likelihood ratio
% test. This is the ideal observer for the Poisson noice cone excitations
% illustrated in this script.  But you can run other rce's as well.
%
% Set up parameters associated with use of this classifier.
switch (whichClassifierEngine)
    case {'rcePoisson'}
        whichClassifierEngine = responseClassifierEngine(@rcePoisson);
        classifierPara = struct('trainFlag', 'none', ...
            'testFlag', 'random', ...
            'nTrain', nTrain, 'nTest', nTest);

    case {'rceTemplateDistance'}
        whichClassifierEngine = responseClassifierEngine(@rceTemplateDistance);
        classifierPara = struct('trainFlag', 'none', ...
            'testFlag', 'random', ...
            'nTrain', nTrain, 'nTest', nTest);

    case 'rcePcaSVM'
        if (nTrain == 1)
            error('Does not really make sense to use rcePcaSVM with only 1 training instance');
        end

        pcaSVMParams = struct(...
            'PCAComponentsNum', 2, ...          % number of PCs used for feature set dimensionality reduction
            'crossValidationFoldsNum', 10, ...  % employ a 10-fold cross-validated linear
            'kernelFunction', 'linear', ...     % linear
            'classifierType', 'svm' ...         % binary SVM classifier
            );
        whichClassifierEngine = responseClassifierEngine(@rcePcaSVM,pcaSVMParams);
        classifierPara = struct('trainFlag', 'none', ...
            'testFlag', 'random', ...
            'nTrain', nTrain, 'nTest', nTest);

    otherwise
        error('Unsupported rce specified')
end


%% Parameters for threshold estimation/quest engine
% The actual threshold varies enough with the different engines that we
% need to adjust the contrast range that Quest+ searches over, as well as
% the range of psychometric function slopes. Threshold limits are computed
% as 10^-logThreshLimitVal.  The reason it is log units is that below we
% define the PF for the questEngine as @qpPFWeibullLog. Had we used the
% default (@qpPFWeibull), the units would have been dB.
%
% Also note explicit passing of proportion correct criterion for threshold.
% The default value of 0.81606 matches the parameterization of mQUESTPlus'
% Weibull PFs, when lapse rate is 0 and guess rate is 0.5.  But it seems
% better to pass it explicitly so we know what it is. Keeping 0.81606 for
% backward compatibilty.
%
% There are two separate structures below. The conceptual distinction
% between them is not entirely clear.  These are interpretted by
% computeThreshold.
%
% See t_spatialCSF.m for more on options of the two different mode of
% operation (fixed numer of trials vs. adaptive

if (isempty(options.psychometricCurveSamplesNum))
    questEnginePara = struct( ...
        'qpPF',@qpPFWeibullLog, ...
        'minTrial', 1280, ... %nTest*psychometricCurveSamplesNum, ...
        'maxTrial', 1280, ... %nTest*psychometricCurveSamplesNum, ...
        'numEstimator', 1, ...
        'stopCriterion', 0.05);
else
    questEnginePara = struct( ...
        'qpPF',@qpPFWeibullLog, ...
        'minTrial', nTest*psychometricCurveSamplesNum, ...
        'maxTrial', nTest*psychometricCurveSamplesNum, ...
        'numEstimator', 1, ...
        'stopCriterion', 0.05);
end
questEnginePara


%% Setup figures
% Generate a figure with a random ID
dataFig = figure(floor(sum(datevec(datetime('now'))*100))); clf;
set(dataFig, 'Position', [10 10 2000 800], 'Color', [1 1 1]);

subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', length(spatialFreqs), ...
       'heightMargin',  0.06, ...
       'widthMargin',    0.02, ...
       'leftMargin',     0.02, ...
       'rightMargin',    0.00, ...
       'bottomMargin',   0.04, ...
       'topMargin',      0.01);

for idx = 1:length(spatialFreqs)
    axTop{idx}  = subplot('Position', subplotPosVectors(1,idx).v);
    axBottom{idx} = subplot('Position', subplotPosVectors(2,idx).v);
end


%% Compute threshold for each spatial frequency
% See toolbox/helpers for functions createGratingSceneEngine, computeThreshold,
% computePeformance

logThreshold = zeros(1, length(spatialFreqs));
for idx = 1:length(spatialFreqs)
    %% Set up EM object to define fixational eye movement paths
    %
    % We can do this separately for each spatial frequency, so that each SF gets a
    % a new draw.  But if we fix the rng seed, we get the same paths for
    % each SF.
    %
    % As written, we use the same EM paths for each contrast we study for a
    % given SF.  This in principle would allow us to use the metaContrast
    % method, as long as we add code that tracks the EM paths and sets up a
    % new metaContrast calculation for each one.
    if (useFixationalEMs)
        % Check that flags are OK.  
        if (useMetaContrast)
            if (~testEMsMatchTrain) || (nTrainEMs > 1) || (nTestEMs > 1)
                error('There must be either no EMs or only one EM path total to use metaContrast method at present');
            end
        end

        % Check that we aren't trying to vary EMs across contrasts, which
        % is not currently possible.
        if (~sameEMsEachContrast)
            error('Varying EMs across contrasts in threshold determination not currently supported by computeThreshold');
        end

        % Set some train EM parameters
        emMicrosaccadeType = 'none';
        emComputeVelocitySignal = false;
        emCenterPaths = false;
        emCenterPathsAtSpecificTimeMsec = [];

        % Setting the seed to -1 means don't touch the seed or call rng()
        % when we call the fixationEM object to generate paths.  We do this
        % so we can control the rng seed here.  Either we always set it the
        % same for the calls to the EM compute method, or we just let it
        % take on its current value.
        emRandomSeed = -1;
        if (sameEMsEachSF)
            seedBeforeEMGeneration = rng(seedForEMs);
        else
            seedBeforeEMGeneration = rng;
        end

        % Set up the train EM object
        trainFixationalEMObj = fixationalEM;
        trainFixationalEMObj.microSaccadeType = emMicrosaccadeType;
        trainFixationalEMObj.randomSeed = emRandomSeed ;
        trainFixationalEMObj.compute(stimulusDuration, frameDurationSeconds, nTrainEMs, emComputeVelocitySignal, ...
            'centerPaths', emCenterPaths, 'centerPathsAtSpecificTimeMsec', emCenterPathsAtSpecificTimeMsec);

        % And this one for testing. 
        testFixationalEMObj = fixationalEM;
        testFixationalEMObj.microSaccadeType = emMicrosaccadeType;
        testFixationalEMObj.randomSeed = emRandomSeed;
        testFixationalEMObj.compute(stimulusDuration, frameDurationSeconds, nTestEMs, emComputeVelocitySignal, ...
            'centerPaths', emCenterPaths, 'centerPathsAtSpecificTimeMsec', emCenterPathsAtSpecificTimeMsec);

        % Make the EM seed what it was before the em path generation.
        rng(seedBeforeEMGeneration);
     
        % For case where test matches training, copy in training EM paths
        % to test EM object.
        if (testEMsMatchTrain)
            if (nTestEMs ~= nTrainEMs)
                error('Number of test and training EMs must match when testEMsMatchTrain is true');
            end
            testFixationalEMObj.emPosArcMin = trainFixationalEMObj.emPosArcMin;
        end
    else
        trainFixationalEMObj = [];
        testFixationalEMObj = [];
    end

    % Create a grating scene engine with a particular chromatic direction,
    % spatial frequency, and temporal duration.  Put grating in sine phase
    % becuase that keeps the spatial mean constant across spatial
    % frequencies.
    %
    % Create scene produces square scenes.  We use the min of the mosaic
    % field size to pick a reasonable size
    gratingSceneEngine = createGratingSceneEngine(chromaDir, spatialFreqs(idx),...
        gratingSceneParams);

    % Set the sceneEngine's visualizeEachCompute property
    gratingSceneEngine.visualizeEachCompute = visualizeEachScene;

    % Create the sceMetaContrast scene engine
    if (useMetaContrast)
        % Instantiate meta contrast neural engine for this spatial
        % frequency and use it to compute threshold
        metaNeuralResponseEngineNoiseFreeParams.sceneEngine = gratingSceneEngine;

        theMetaNeuralEngine = neuralResponseEngine(@nreNoiseFreeMetaContrast, ...
            @nreNoisyInstancesMetaContrast, ...
            metaNeuralResponseEngineNoiseFreeParams, ...
            metaNeuralResponseEngineNoisyInstanceParams);

        % Update visualizeEachCompute
        theMetaNeuralEngine.visualizeEachCompute = theNeuralEngine.visualizeEachCompute;

        % Compute the threshold for our grating scene with meta scene and
        % and neural response engines. This function does a lot of work,
        % see the function itself, as well as function computePerformance.
        [logThreshold(idx), questObj, ~, para(idx,:)] = ...
            computeThreshold(theMetaSceneEngine, theMetaNeuralEngine, whichClassifierEngine, ...
            classifierPara, thresholdPara, questEnginePara, ...
            'TAFC', true, 'useMetaContrast', useMetaContrast, ...
            'trainFixationalEM', trainFixationalEMObj, ...
            'testFixationalEM', testFixationalEMObj, ...
            'verbose', verbose, ...
            'maxVisualizedNoisyResponseInstanceStimuli', maxVisualizedNoisyResponseInstanceStimuli);
    else
        % Compute the threshold for our grating scene with the previously
        % defined neural and classifier engine.  This function does a lot
        % of work, see the function itself, as well as function
        % computePerformance.
        [logThreshold(idx), questObj, ~, para(idx,:)] = ...
            computeThreshold(gratingSceneEngine, theNeuralEngine, whichClassifierEngine, ...
            classifierPara, thresholdPara, questEnginePara, ...
            'TAFC', true, 'useMetaContrast', useMetaContrast, ...
            'trainFixationalEM', trainFixationalEMObj, ...
            'testFixationalEM', testFixationalEMObj, ...
            'verbose', verbose, ...
            'maxVisualizedNoisyResponseInstanceStimuli', maxVisualizedNoisyResponseInstanceStimuli);
    end

    % Plot stimulus & psychometric curve
    figure(dataFig);

    % This shows one frame of the scene.
    visualizeEachComputeSave = gratingSceneEngine.visualizeEachCompute;
    gratingSceneEngine.visualizeEachCompute = false;
    visualizationContrast = 1.0;
    [theSceneSequence] = gratingSceneEngine.compute(visualizationContrast);
    gratingSceneEngine.visualizeStaticFrame(...
        theSceneSequence, ...
        'frameToVisualize', 1, ...
        'axesHandle', axTop{idx});
    gratingSceneEngine.visualizeEachCompute = visualizeEachComputeSave;

    % Plot data and psychometric curve
    % with a marker size of 2.5
    questObj.plotMLE(2.5,'para',para(idx,:), 'axesHandle', axBottom{idx});
    drawnow;
end
set(dataFig, 'Position',  [0, 0, 800, 800]);
saveas(dataFig,[figureFileBase '_Psychometric.tiff'],'tif');

% Convert returned log threshold to linear threshold
thresholdContrasts = 10 .^ logThreshold;

% Save thresholds
save(fullfile(resultsFileBase,thresholdsDataFileName), 'options', 'spatialFreqs', 'thresholdContrasts');

%% Plot contrast sensitivity function
theCsfFig = figure();
loglog(spatialFreqs, 1 ./ thresholdContrasts, '-ok', 'LineWidth', 2);
xticks(spatialFreqs); xlim([spatialFreqs(1), spatialFreqs(end)]);
if (framesNum == 1)
    yticks([1,2,5,10,20,50]); ylim([1, 50]);
elseif (framesNum == 4)
    yticks([1,2,5,10,20,50,100,200]); ylim([1, 200]);
end
xlabel('Spatial Frequency (cyc/deg)');
ylabel('Sensitivity');
set(theCsfFig, 'Position',  [800, 0, 600, 800]);
saveas(theCsfFig,[figureFileBase,'_CSF.tiff'],'tif');


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
    if (any(abs(thresholdContrasts-validationThresholds)./validationThresholds > validationTolerance))
        thresholdContrasts(:)
        validationThresholds(:)
        error(sprintf('Do not replicate validation thresholds to %d%%. Check that parameters match, or for a bug.',round(100*validationTolerance)));
    else
        fprintf('Validation regression check passes\n');
    end
else
    thresholdContrasts(:)
    fprintf('No validation thresholds, validation regression check not run\n');
end

%% Return threshold values if requested
if (nargout > 0)
    thresholdRet = threshold;
end

%% Save output if desired
if (~isempty(options.resultsFileBase))
    close all;
    drawnow;
    save(resultsFileBase);
end

end


