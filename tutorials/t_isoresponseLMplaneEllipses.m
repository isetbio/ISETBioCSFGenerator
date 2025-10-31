function t_isoresponseLMplaneEllipses(options)
% Compute iso-response contours on the LM-plane
%
% Syntax:
%   t_isoresponseLMplaneEllipses
%
% Description:
%   Illustrates how to compute threshold iso-response ellipses on the LM
%   contrast plane using either CMosaic- or  mRGCmosaic- based response engines. 
%   In the case that an mRGCMosaic response engine is used, this tutorial also 
%   illustrates two additional options;
%   (a) how to simulate an mRGCmosaic with light increment (ON) and 
%       light-decrement (OFF) neurons, although we currently only have ON-center mRGCmosaics
%   (b) how to simulate a non-linear activation function 
%
%   This tutorial was used to generate the examples ellipses that were
%   included in the computational aim of the 2025 R-01 grant.
%
% History:
%   01/27/2025  NPC   Wrote it.

% Examples:
%{

% Run with CMosaic
     t_isoresponseLMplaneEllipses(...
        'useMetaContrast', true, ...
        'whichNoiseFreeNre', 'excitationsCmosaic', ...
        'whichNoisyInstanceNre', 'Poisson', ...
        'whichClassifierEngine', 'rcePoisson', ...
        'useConeContrast', false);

% Run with mRGCMosaic (linear response, only ON-center mosaic)
     t_isoresponseLMplaneEllipses(...
        'useMetaContrast', true, ...
        'whichNoiseFreeNre', 'mRGCMosaic', ...
        'mRGCOutputSignalType', 'mRGCs', ...      % Select between {'cones', 'mRGCs'}
        'whichNoisyInstanceNre', 'Gaussian', ...
        'whichClassifierEngine', 'rceTemplateDistance', ...
        'useConeContrast', true);


%Pts_coneContrast =
%    0.2384    0.3549   -0.0000
%   -0.3058   -0.2796    0.0000
%    0.0000    0.0000         0
%    0.3058    0.2796   -0.0000
%   -0.2384   -0.3549    0.0000


% % % Run with mRGCMosaic (ON mRGC mosaic, non-linear activation function, simulating case in a high saturation regime)
% %      t_isoresponseLMplaneEllipses(...
% %         'useMetaContrast', true, ...
% %         'whichNoiseFreeNre', 'mRGCMosaic', ...
% %         'mosaicEccDegs', [0.0 0.5], ...
% %         'mosaicSizeDegs', [1 0.5], ...
% %         'synthesizedRGCmosaicName', 'PLOSpaperFovealMosaic', ...
% %         'visualizeCroppedRGCmosaicRelationshipToSynthesizedMosaic', true, ...
% %         'referenceLMSconeContrast', [0.2384  0.3549 0.0], ...
% %         'spatialFrequency', 0.0, ...
% %         'spatialPhaseDegs', 0.0, ...
% %         'presentationMode', 'counterphasemodulated', ...
% %         'stimDurationTemporalCycles', 1.0, ...
% %         'frameDurationSeconds', (1/2.5)/8, ...
% %         'temporalFrequencyHz', 2.5, ...
% %         'nTest', 1024, ...
% %         'psychometricCurveSamplesNum', 5, ...
% %         'examinedDirectionsOnLMplane', 0:45:(360-45), ...
% %         'visualizeEachScene', true, ...
% %         'simulateONOFFmosaic', false, ...
% %         'simulateHighSaturationRegime', ~true, ...
% %         'mRGCOutputSignalType', 'mRGCs', ...  
% %         'whichNoisyInstanceNre', 'Gaussian', ...
% %         'gaussianSigma', 0.1, ...
% %         'whichClassifierEngine', 'rceTemplateDistance', ...
% %         'useConeContrast', true);


%Pts_coneContrast =
%    0.2384    0.3549   -0.0000
%   -0.3058   -0.2796    0.0000
%    0.0000    0.0000         0
%    0.3058    0.2796   -0.0000
%   -0.2384   -0.3549    0.0000

% MRGCs (non-linear)
 t_isoresponseLMplaneEllipses(...
        'useMetaContrast', true, ...
        'whichNoiseFreeNre', 'mRGCMosaic', ...
        'mosaicEccDegs', [0.0 0.0], ...
        'mosaicSizeDegs', [1 1], ...
        'synthesizedRGCmosaicName', 'PLOSpaperFovealMosaic', ...
        'visualizeCroppedRGCmosaicRelationshipToSynthesizedMosaic', true, ...
        'backgroundLMSconeExcitations', [0.1529 0.1324 0.0828], ...
        'referenceLMSconeContrast', [0.2384    0.3549   -0.0000], ...
        'spatialFrequency', 0.0, ...
        'spatialPhaseDegs', 0.0, ...
        'presentationMode', 'static', ...
        'nTest', 1024, ...
        'psychometricCurveSamplesNum', 5, ...
        'examinedDirectionsOnLMplane', 0:45:(360-45), ...
        'visualizeEachScene', true, ...
        'simulateONOFFmosaic', false, ...
        'simulateHighSaturationRegime', true, ...
        'mRGCOutputSignalType', 'cones', ...  
        'whichNoisyInstanceNre', 'Gaussian', ...
        'gaussianSigma', 0.5, ...
        'whichClassifierEngine', 'rceTemplateDistance', ...
        'useConeContrast', true);

% MRGCs (linear)
 t_isoresponseLMplaneEllipses(...
        'useMetaContrast', true, ...
        'whichNoiseFreeNre', 'mRGCMosaic', ...
        'mosaicEccDegs', [0.0 0.0], ...
        'mosaicSizeDegs', [1 1], ...
        'synthesizedRGCmosaicName', 'PLOSpaperFovealMosaic', ...
        'visualizeCroppedRGCmosaicRelationshipToSynthesizedMosaic', true, ...
        'backgroundLMSconeExcitations', [0.1529 0.1324 0.0828], ...
        'referenceLMSconeContrast', [0.2384    0.3549   -0.0000], ...
        'spatialFrequency', 0.0, ...
        'spatialPhaseDegs', 0.0, ...
        'presentationMode', 'static', ...
        'nTest', 1024, ...
        'psychometricCurveSamplesNum', 5, ...
        'examinedDirectionsOnLMplane', 0:45:(360-45), ...
        'visualizeEachScene', true, ...
        'simulateONOFFmosaic', false, ...
        'simulateHighSaturationRegime', ~true, ...
        'mRGCOutputSignalType', 'cones', ...  
        'whichNoisyInstanceNre', 'Gaussian', ...
        'gaussianSigma', 0.5, ...
        'whichClassifierEngine', 'rceTemplateDistance', ...
        'useConeContrast', true);

% CONES (input cone mosaic of MRC mosaic, using cone contrast and Gaussian
noise with same sigma as mRGCs)


t_isoresponseLMplaneEllipses(...
        'useMetaContrast', true, ...
        'whichNoiseFreeNre', 'mRGCMosaic', ...
        'mosaicEccDegs', [0.0 0.0], ...
        'mosaicSizeDegs', [1 1], ...
        'synthesizedRGCmosaicName', 'PLOSpaperFovealMosaic', ...
        'visualizeCroppedRGCmosaicRelationshipToSynthesizedMosaic', true, ...
        'backgroundLMSconeExcitations', [0.1529 0.1324 0.0828], ...
        'referenceLMSconeContrast', [-0.2384   -0.3549    0.0000], ...
        'spatialFrequency', 0.0, ...
        'spatialPhaseDegs', 0.0, ...
        'presentationMode', 'static', ...
        'nTest', 1024, ...
        'psychometricCurveSamplesNum', 5, ...
        'examinedDirectionsOnLMplane', 0:45:(360-45), ...
        'mRGCOutputSignalType', 'cones', ...
        'whichNoisyInstanceNre', 'Gaussian', ...
        'gaussianSigma', 1.0, ...
        'whichClassifierEngine', 'rceTemplateDistance', ...
        'useConeContrast', true);


% Run with mRGCMosaic (ON mRGC mosaic, half-wave recrifier activation function, simulating case in subthreshold regime)
t_isoresponseLMplaneEllipses(...
        'useMetaContrast', true, ...
        'whichNoiseFreeNre', 'mRGCMosaic', ...
        'simulateHalfWaveRectification', true, ...
        'simulateONOFFmosaic', false, ...
        'mRGCOutputSignalType', 'mRGCs', ...  
        'whichNoisyInstanceNre', 'Gaussian', ...
        'whichClassifierEngine', 'rceTemplateDistance', ...
        'useConeContrast', true);

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


    % If the neural model is mRGCMosaic, we need to specify the following

    % See RGCMosaicConstructor.helper.utils.initializeRGCMosaicGenerationParameters
    % for what is available and to add new mosaics
    options.synthesizedRGCmosaicName (1,:) char = 'PLOSpaperFovealMosaic';

    % ---- Which species to employ ----
    % Choose between {'macaque', 'human'}. If 'macaque' is chosen, the input
    % cone mosaic has a 1:1 L/M cone ratio.
    options.synthesizedRGCmosaicConeMosaicSpecies  (1,:) char {mustBeMember(options.synthesizedRGCmosaicConeMosaicSpecies,{'human','macaque'})} = 'human';

    % ----- Which subject optics to employ -----
    options.synthesizedRGCmosaicOpticsSubjectName (1,:) ...
        char ...
        {...
        mustBeMember(options.synthesizedRGCmosaicOpticsSubjectName, ...
            { ...
            'PLOSpaperDefaultSubject' ...
            'PLOSpaperSecondSubject' ...
            'VSS2024TalkFirstSubject' ...
            'VSS2024TalkSecondSubject' ...
            'PLOSpaperStrehlRatio_0.87' ...
            'PLOSpaperStrehlRatio_0.72' ...
            'PLOSpaperStrehlRatio_0.59' ...
            'PLOSpaperStrehlRatio_0.60' ...
            'PLOSpaperStrehlRatio_0.27' ...
            'PLOSpaperStrehlRatio_0.23' ...
            'PLOSpaperStrehlRatio_0.21' ...
            'PLOSpaperStrehlRatio_0.19' ...
            'PLOSpaperStrehlRatio_0.09' ...
            } ...
            ) ...
        } ...
        = 'PLOSpaperSecondSubject';


    % ------ targetVisualSTF options ----
    % Options are : {'default', 'x1.3 RsRcRatio'}
    % These are with respect to the macaque data of the Croner & Kaplan '95 study
    % 'default': target the mean Rs/Rc, and the mean Ks/Kc (Rs/Rc)^2
    % See RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.generateTargetVisualSTFmodifiersStruct
    % for all existing options
    options.synthesizedRGCmosaicTargetVisualSTFdescriptor (1,:) char = 'default';

    % Whether to visualize the relationship of the cropped RGC mosaic to
    % the synthesized RGC mosaic
    options.visualizeCroppedRGCmosaicRelationshipToSynthesizedMosaic (1,1) logical = false

    % If the neural model is mRGCMosaic, we can specify the Optics type to use
    options.opticsType (1,:) char  = 'loadComputeReadyRGCMosaic';

    % Choose noise model
    %   Choices: 'Poisson'
    %            'Gaussian'
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
    % For mRGCmosaic, this will be the crop size and eccentricity
    % so it has to be appropriate to the synthesizedRGCmosaicName
    options.mosaicEccDegs (1,2) double = [0 0];
    options.mosaicSizeDegs (1,2) double = [0.5 0.5];

    % Varied stimulus parameter
    options.examinedDirectionsOnLMplane (1,:) double = 0:45:315;

    % The reference LMS cone contrat with respect to which the 
    % the examinedDirectionsOnLMplane are used for compute the
    % tested LMS cone contrsts.
    % For contrasts symmetric around the background,
    % referenceLMSconeContrast should be set to [0 0 0]
    options.referenceLMSconeContrast (1,3) double = [0 0 0];
    
    options.backgroundLMSconeExcitations (1,:) double = [];

    % Fixed stimulus parameters
    options.employMosaicSpecificConeFundamentals (1,1) logical = true;
    options.meanLuminanceCdPerM2 (1,1) double = 100;
    options.meanChromaticityXY (1,2) double = [0.30 0.32];
    options.spatialFrequency(1,1) double = 0.0;
    options.orientationDegs (1,1) double = 90;
    options.spatialPhaseDegs (1,1) double = 0
    options.numberOfFrames (1,:) double = []
    options.stimOnFrameIndices (1,:) double = [];
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
        'slopeRangeHigh', 200/20, ...
        'slopeDelta', 2.5/50, ...
        'thresholdCriterion', 0.81606, ...
        'guessRate', 1/2, ...
        'lapseRate', 0);

    % Number of psychometric curve samples
    options.psychometricCurveSamplesNum (1,:) double = 3
    
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

    % Simulate a high saturation regime
    options.simulateHighSaturationRegime (1,1) logical = false

    % Simulate a high saturation regime
    options.simulateHalfWaveRectification (1,1) logical = false

    % Simulate an ON/OFF mRGC mosaic?
    options.simulateONOFFmosaic (1,1) logical = false

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

% Properties of snthesized mRGCmosaic to use
synthesizedRGCmosaicName = options.synthesizedRGCmosaicName;
synthesizedRGCmosaicOpticsSubjectName = options.synthesizedRGCmosaicOpticsSubjectName;
synthesizedRGCmosaicTargetVisualSTFdescriptor = options.synthesizedRGCmosaicTargetVisualSTFdescriptor;
synthesizedRGCmosaicConeMosaicSpecies = options.synthesizedRGCmosaicConeMosaicSpecies;
visualizeCroppedRGCmosaicRelationshipToSynthesizedMosaic = options.visualizeCroppedRGCmosaicRelationshipToSynthesizedMosaic;

% Mosaic sizes
mosaicEccDegs = options.mosaicEccDegs;
mosaicSizeDegs = options.mosaicSizeDegs;

% Stimulus params
% Custom cone fundamentals
employMosaicSpecificConeFundamentals = options.employMosaicSpecificConeFundamentals;

% Background
meanLuminanceCdPerM2 = options.meanLuminanceCdPerM2;
meanChromaticityXY = options.meanChromaticityXY;

% Background LMS
backgroundLMSconeExcitations = options.backgroundLMSconeExcitations;


% Spatial params
spatialPhaseDegs = options.spatialPhaseDegs;
stimSizeDegs = options.stimSizeDegs;
pixelsNum = options.pixelsNum;

% Temporal params
temporalFrequencyHz = options.temporalFrequencyHz;
stimOnFrameIndices = options.stimOnFrameIndices;
temporalFilterValues = options.temporalFilterValues;
numberOfFrames = options.numberOfFrames;
frameDurationSeconds = options.frameDurationSeconds;
stimDurationTemporalCycles = options.stimDurationTemporalCycles;
presentationMode = options.presentationMode;

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
visualizeNonLinearActivationFunction = ~true;

thresholdsDataFileName = options.thresholdsDataFileName;
parPoolSize = options.parPoolSize;

simulateONOFFmosaic = options.simulateONOFFmosaic;
simulateHighSaturationRegime = options.simulateHighSaturationRegime;
simulateHalfWaveRectification = options.simulateHalfWaveRectification;


%% Freeze rng for replicatbility and validation
rng(1);

%% Close any stray figs
hAllFigs = findall(groot,'Type','figure');

% Close all figures
for i = 1:numel(hAllFigs)
    set(hAllFigs(i), 'HandleVisibility', 'on')
end
close all;

%% Make sure local/figures directory exists so we can write out our figures in peace
[figureFileBaseDir, resultsFileBaseDir] = setupFigureDirectory(mfilename, ...
    useMetaContrast, useConeContrast, useFixationalEMs, ...
    whichNoiseFreeNre, whichNoisyInstanceNre,...
    whichClassifierEngine, mRGCOutputSignalType);

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

%% List of LM angles to be tested
examinedDirectionsOnLMplane = options.examinedDirectionsOnLMplane;
referenceLMSconeContrast = options.referenceLMSconeContrast;


%% Spatial frequency & orientation
examinedSpatialFrequencyCPD = options.spatialFrequency;
orientationDegs = options.orientationDegs;

debugStimulusConfig = ~true;



%% Thresholds filename
if (isempty(thresholdsDataFileName))
    thresholdsDataFileName = ...
        sprintf('%sCSF_SF_%2.2fCPD_Optics_%s_EccDegs_x%2.1f_%2.1f_SizeDegs_%2.1fx%2.1f_OriDegs_%2.0f_RefLMScontrast_%2.1f_%2.1f_%2.1f%s.mat', ...
        mRGCOutputSignalType, ...
        examinedSpatialFrequencyCPD, ...
        opticsType, ...
        mosaicEccDegs(1), mosaicEccDegs(2), ...
        mosaicSizeDegs(1),mosaicSizeDegs(2), ...
        orientationDegs, ...
        referenceLMSconeContrast(1)*100, ...
        referenceLMSconeContrast(2)*100, ...
        referenceLMSconeContrast(3)*100, ...
        presentationMode);
end

fprintf('Thresholds will be saved to %s', fullfile(resultsFileBaseDir,thresholdsDataFileName));

% Max RMS contrast (so as to keep stimuli within the display gamut)
rmsLMconeContrast = 0.05;

if (debugStimulusConfig)
    rmsLMconeContrast = 1.0;
    thresholdPara.logThreshLimitLow = 0.1;
end

% Compute directions based on rmsLMconeContrast
[theLMSconeContrastDirections, examinedDirectionsOnLMplane] = ...
    computeLMSconeContrastDirections(rmsLMconeContrast, examinedDirectionsOnLMplane);

% Add the reference LMScone contrast. 
% For contrasts symmetric around the background,
% referenceLMSconeContrast should be set to [0 0 0]
theLMSconeContrastDirections = bsxfun(@plus, theLMSconeContrastDirections, referenceLMSconeContrast(:));


%% Set grating engine parameters
if (strcmp(presentationMode, 'static'))
    gratingSceneParams = struct( ...
            'meanLuminanceCdPerM2', meanLuminanceCdPerM2, ...
            'meanChromaticityXY', meanChromaticityXY, ...
            'backgroundLMSconeExcitations', backgroundLMSconeExcitations, ...
            'spectralSupport', 400:20:750, ...
            'fovDegs', stimSizeDegs, ...
            'pixelsNum', pixelsNum, ...
            'spatialEnvelope', 'rect', ...
            'spatialEnvelopeRadiusDegs', stimSizeDegs, ...
            'orientation', orientationDegs, ...
            'presentationMode', 'flashed', ... 
            'duration', frameDurationSeconds, ...
            'frameDurationSeconds', frameDurationSeconds, ...
            'spatialPhase', spatialPhaseDegs);
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
        'backgroundLMSconeExcitations', backgroundLMSconeExcitations, ...
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

%% Which stimuli to skip
skippedDirections = 0;
if (length(examinedDirectionsOnLMplane) > 50)
    % Too many, plot every other stimulus
    skippedDirections = 2;
end


if (employMosaicSpecificConeFundamentals)
    setupStepsRequired = 2;
else
    setupStepsRequired = 1;
end

for iStep = 1:setupStepsRequired

    %% Configure our stimulus scene engines 
    [theNullSceneEngine, theTestSceneEngines] = configureStimulusSceneEngines(...
        theLMSconeContrastDirections, examinedSpatialFrequencyCPD, gratingSceneParams);

    if ((employMosaicSpecificConeFundamentals) && (iStep == setupStepsRequired))
        % Do not redo the neural engine installation. Just the stimulus
        % scene engines above (with the customConeFundamentals, the second
        % time around)
        continue;
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
            % Pre-baked mRGC mosaic name
            nreNoiseFreeParams.mosaicParams.rgcMosaicName = synthesizedRGCmosaicName;
    
            % Optics used to synthesize the mRGCMosaic
            nreNoiseFreeParams.mosaicParams.opticsSubjectName = synthesizedRGCmosaicOpticsSubjectName;
    
            % Target visualSTF of synthesized mRGCmosaic
            nreNoiseFreeParams.mosaicParams.targetVisualSTFdescriptor = synthesizedRGCmosaicTargetVisualSTFdescriptor;
    
            % Input cone mosaic with a cone density for human retinas
            nreNoiseFreeParams.mosaicParams.coneMosaicSpecies = synthesizedRGCmosaicConeMosaicSpecies;


            % 2. We can crop the mRGCmosaic to some desired size.
            %    Passing [] for sizeDegs will not crop.
            %    Passing [] for eccentricityDegs will crop the mosaic at its center.
            nreNoiseFreeParams.mRGCMosaicParams.cropParams = struct(...
                'sizeDegs', mosaicSizeDegs, ...
                'eccentricityDegs', mosaicEccDegs, ...
                'visualizeSpatialRelationshipToSourceMosaic', visualizeCroppedRGCmosaicRelationshipToSynthesizedMosaic ...
                );
    
            % 3. Set the input cone mosaic integration time to match the stimulus frame duration
            nreNoiseFreeParams.mRGCMosaicParams.coneIntegrationTimeSeconds = gratingSceneParams.frameDurationSeconds;
    
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
    
    
            % NON-LINEARITY HACKS
            % Simulate incomplete background adaptation - puts mosaic in
            % saturation regime of the activation function
    
            nreNoiseFreeParams.mRGCMosaicParams.responseBias = 0.0;
    
            if (simulateHighSaturationRegime)
                % Push responses into high saturation regime
                nreNoiseFreeParams.mRGCMosaicParams.responseBias = 0.03;
            end
    
            if (simulateHalfWaveRectification)
                % Push responses to the subthreshold regime
                nreNoiseFreeParams.mRGCMosaicParams.responseBias = -0.05;
            end
    
            % Simulate ON-OFF mosaic
            % All odd-indexes mRGCs will be treated as OFF-center, by inverting their 
            % noise-free responses polarity in nreNoiseFreeMidgetRGCMosaic()
            nreNoiseFreeParams.mRGCMosaicParams.simulateONOFFmosaic = simulateONOFFmosaic;
    
    
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
            nreNoiseFreeParams.coneMosaicParams.timeIntegrationSeconds = gratingSceneParams.frameDurationSeconds;
    
            % Handle cone contrast setting
            if (useConeContrast)
                nreNoiseFreeParams.coneMosaicParams.outputSignalType = 'coneContrast';
            else
                nreNoiseFreeParams.coneMosaicParams.outputSignalType = 'coneExcitations';
            end
    
            % No temporal filter
            nreNoiseFreeParams.temporalFilter = temporalFilter;
    
            % Set Gaussian sigma in case we're using Gaussian noise below.
            if (isempty(gaussianSigma))
                gaussianSigma = 50;
            end
   
    
        otherwise
            error('Unsupported noise free neural response engine: ''%s''.', whichNoiseFreeNre);
    end  % switch (whichNoiseFreeNre)

   
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

        [theOI,theMosaic] = generateOpticsAndMosaicFromParams(...
            nreNoiseFreeParams.opticsParams, ...
            [], ...
            mosaicParams);

        switch (whichNoiseFreeNre)
           case 'mRGCMosaic'
               theConeMosaic = theMosaic.inputConeMosaic;
           case 'excitationsCmosaic'
               theConeMosaic = theMosaic;
        end

        % Compute the custom cone fundamentals and store them in gratingSceneParams
        % so in the next pass, we can generate the desired scene
        maxConesNumForAveraging = 3;
        gratingSceneParams.customConeFundamentals = visualStimulusGenerator.coneFundamentalsForPositionWithinConeMosaic(...
            theConeMosaic, theOI, ...
            mosaicEccDegs, gratingSceneParams.fovDegs, maxConesNumForAveraging);

        % Dont need these anymore
        clear 'theConeMosaic'
        clear 'theMosaic'
    end

end %for iStep = 1:setupStepsRequired


% Setup the noisy neural response engine
switch (whichNoisyInstanceNre)
    case 'Poisson'
        nreNoisyInstancesComputeFunction = @nreNoisyInstancesPoisson;
        nreNoisyInstancesParams = nreNoisyInstancesPoisson;
       
        % No activation function for cone mosaic excitations
        %nreNoisyInstancesParams.activationFunctionParams = [];

    case 'Gaussian'
        nreNoisyInstancesComputeFunction = @nreNoisyInstancesGaussian;
        nreNoisyInstancesParams = nreNoisyInstancesGaussian;
        nreNoisyInstancesParams.sigma = gaussianSigma;
        
        if (simulateHighSaturationRegime) && (simulateHalfWaveRectification)
            error('Either select ''simulateHighSaturationRegime'' or ''simulateHalfWaveRectification'', NOT both.')
        end

        if (simulateHighSaturationRegime)
            % Saturating activation function. If the type is specified as
            % 'halfWaveSigmoidalRectifier', the response bias specified above
            % like so: nreNoiseFreeParams.mRGCMosaicParams.responseBias
            % specifies how the noise-free response is pushed into different
            % regimes of the sigmoidal function
            nreNoisyInstancesParams.activationFunctionParams = struct(...
                'type', 'halfwaveSigmoidalRectifier', ...                 % choose between {'linear', 'halfwaveRectifier', 'halfwaveSigmoidalRectifier'}
                'exponent', 2.0, ...                                      % only relevant for 'halfwaveSigmoidalRectifier'
                'semiSaturationReponseAmplitude', 0.15*gaussianSigma, ... % only relevant for 'halfwaveSigmoidalRectifier'
                'gain', 3*gaussianSigma, ...                              % only relevant for 'halfwaveSigmoidalRectifier'
                'visualize', visualizeNonLinearActivationFunction ...
                );
        elseif (simulateHalfWaveRectification)
            nreNoisyInstancesParams.activationFunctionParams = struct(...
                'type', 'halfwaveRectifier', ...                 % choose between {'linear', 'halfwaveRectifier', 'halfwaveSigmoidalRectifier'}
                'visualize', visualizeNonLinearActivationFunction ...
                );
        end

    otherwise
        error('Unsupported noisy instances neural response engine: ''%s''.', whichNoisyInstanceNre);
end % switch (whichNoisyInstanceNre)


%% If we use cone contrast, we will neeed a null scene for normalization.
if (useConeContrast) || (ischar(temporalFilterValues))
   nreNoiseFreeParams.nullStimulusSceneSequence = theNullSceneEngine.compute(0.0);
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
% test. This is the ideal observer for the Poisson noise cone excitations
%
% Set up parameters associated with use of this classifier.
switch (whichClassifierEngine)
    case {'rcePoisson'}
        theClassifierEngine = responseClassifierEngine(@rcePoisson);
        classifierEngineParams = struct(...
            'trainFlag', 'none', ...
            'testFlag', 'random', ...
            'nTrain', 1, 'nTest', nTest);

    case {'rceTemplateDistance'}
        theClassifierEngine = responseClassifierEngine(@rceTemplateDistance);
        classifierEngineParams = struct(...
            'trainFlag', 'none', ...
            'testFlag', 'random', ...
            'nTrain', 1, 'nTest', nTest);

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
        error('Unsupported response classifier engine: ''%s''.', whichClassifierEngine)
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

psychometricCurveSamplesNum = options.psychometricCurveSamplesNum;

questEnginePara = struct( ...
    'qpPF',@qpPFWeibullLog, ...
    'minTrial', nTest*psychometricCurveSamplesNum, ...
    'maxTrial', nTest*psychometricCurveSamplesNum, ...
    'numEstimator', 1, ...
    'stopCriterion', 0.05);


%% Setup figures
psychometricCurvesFig = figure(floor(sum(datevec(datetime('now'))*100))); clf;
set(psychometricCurvesFig, 'Position', [10 10 2000 800], 'Color', [1 1 1]);

subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', length(examinedDirectionsOnLMplane), ...
       'heightMargin',  0.06, ...
       'widthMargin',    0.02, ...
       'leftMargin',     0.02, ...
       'rightMargin',    0.00, ...
       'bottomMargin',   0.04, ...
       'topMargin',      0.01);

for idx = 1:length(examinedDirectionsOnLMplane)
    axTop{idx}  = subplot('Position', subplotPosVectors(1,idx).v);
    axBottom{idx} = subplot('Position', subplotPosVectors(2,idx).v);
end



%% Compute threshold for each chromatic angle
% See toolbox/helpers for functions createGratingSceneEngine, computeThreshold,
% computePeformance

logThreshold = zeros(1, length(examinedDirectionsOnLMplane));
for iChromaDirection = 1:length(examinedDirectionsOnLMplane)
    if (useFixationalEMs)
        % Need to do stuff here
    else
        trainFixationalEMObj = [];
        testFixationalEMObj = [];
    end

    % Get the engine for the tested chromatic direction
    theSceneEngine = theTestSceneEngines{iChromaDirection};
    
    % Set the sceneEngine's visualizeEachCompute property
    theSceneEngine.visualizeEachCompute = visualizeEachScene;
    
     % Create the sceMetaContrast scene engine
     if (useMetaContrast)
            % Instantiate meta contrast neural engine for this spatial
            % frequency and use it to compute threshold
            metaNeuralResponseEngineNoiseFreeParams.sceneEngine = theSceneEngine;
    
            theMetaNeuralEngine = neuralResponseEngine(...
                @nreNoiseFreeMetaContrast, ...
                @nreNoisyInstancesMetaContrast, ...
                metaNeuralResponseEngineNoiseFreeParams, ...
                metaNeuralResponseEngineNoisyInstanceParams);
            
            % Update visualizeEachCompute
            theMetaNeuralEngine.visualizeEachCompute = theNeuralEngine.visualizeEachCompute;
    
            % Compute the threshold for our grating scene with meta scene and
            % and neural response engines. This function does a lot of work,
            % see the function itself, as well as function computePerformance.
            [logThreshold(iChromaDirection), questObj, ~, psychometricCurveParams(iChromaDirection,:)] = ...
                computeThreshold(theMetaSceneEngine, theMetaNeuralEngine, theClassifierEngine, ...
                classifierEngineParams, thresholdPara, questEnginePara, ...
                'TAFC', true, 'useMetaContrast', useMetaContrast, ...
                'trainFixationalEM', trainFixationalEMObj, ...
                'testFixationalEM', testFixationalEMObj, ...
                'verbose', verbose, ...
                'maxVisualizedNoisyResponseInstances', maxVisualizedNoisyResponseInstances);
     else
        % Compute the threshold for our grating scene with the previously
        % defined neural and classifier engine.  This function does a lot
        % of work, see the function itself, as well as function
        % computePerformance.
        [logThreshold(iChromaDirection), questObj, ~, psychometricCurveParams(iChromaDirection,:)] = ...
                computeThreshold(theMetaSceneEngine, theSceneEngine, theClassifierEngine, ...
                classifierEngineParams, thresholdPara, questEnginePara, ...
                'TAFC', true, 'useMetaContrast', useMetaContrast, ...
                'trainFixationalEM', trainFixationalEMObj, ...
                'testFixationalEM', testFixationalEMObj, ...
                'verbose', verbose, ...
                'maxVisualizedNoisyResponseInstances', maxVisualizedNoisyResponseInstances);
     end % useMetaContrast
    
     % Plot stimulus & psychometric curve
     figure(psychometricCurvesFig);
     %set(psychometricCurvesFig, 'HandleVisibility', 'on');
        
     % This shows one frame of the scene.
     visualizeEachComputeSave = theSceneEngine.visualizeEachCompute;
     theSceneEngine.visualizeEachCompute = false;

     theSceneEngine.visualizeStaticFrame(...
        theSceneEngine.compute(1.0), ...
        'frameToVisualize', 1, ...
        'axesHandle', axTop{iChromaDirection});
     theSceneEngine.visualizeEachCompute = visualizeEachComputeSave;


     % Plot data and fitted psychometric curve
     questObj.plotMLE(2.5,'para', psychometricCurveParams(iChromaDirection,:), ...
            'axesHandle', axBottom{iChromaDirection});
     drawnow;
     %set(psychometricCurvesFig, 'HandleVisibility', 'off');
    
end % iChromaDirection

% Convert returned log threshold to linear threshold
thresholdContrasts = 10 .^ logThreshold;

% Save thresholds
stimulusRMSLMconeContrast = sqrt(sum(theLMSconeContrastDirections.^2,1));

save(fullfile(resultsFileBaseDir,thresholdsDataFileName), ...
    'options', 'stimulusRMSLMconeContrast', 'examinedDirectionsOnLMplane', 'thresholdContrasts');

% Figure for plotting the thresholds on the LM plane
% along with the stimuli
hFigStimuliAndThresholds = figure(2346); clf;
set(hFigStimuliAndThresholds, 'Color', [1 1 1], 'Position', [10 10 1200 1200]);
set(hFigStimuliAndThresholds, 'HandleVisibility', 'off');

exportFig = true;
initializeFig = true;
theThresholdAxes = []; 
maxVisualizedThreshold = 0.5;
figName = sprintf('refC_%2.1f_%2.1f_%2.1f', 100*referenceLMSconeContrast(1), 100*referenceLMSconeContrast(2), 100*referenceLMSconeContrast(3));
   
visualizeIsoThresholdEllipsesOnLMplane(...
            stimulusRMSLMconeContrast, ...
            examinedSpatialFrequencyCPD, gratingSceneParams, ...
            examinedDirectionsOnLMplane, skippedDirections, ...
            thresholdContrasts.*stimulusRMSLMconeContrast, maxVisualizedThreshold, referenceLMSconeContrast(1:2), ...
            figureFileBaseDir, hFigStimuliAndThresholds, ...
            exportFig, initializeFig, theThresholdAxes, figName);

%% Restore handle visibility so that a future close will close the figures
hAllFigs = findall(groot,'Type','figure');
for i = 1:numel(hAllFigs)
    set(hAllFigs(i), 'HandleVisibility', 'on')
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

end


%
% SUPPORTING FUNCTIONS
%



function [theNullSceneEngine, theTestSceneEngines] = configureStimulusSceneEngines(...
    theLMSconeContrastDirections, examinedSpatialFrequencyCPD, gratingSceneParams)
 
    % Create the background scene engine
    chromaticDir = [0 0 0.4];
    theNullSceneEngine = createGratingSceneEngine(...
            chromaticDir, examinedSpatialFrequencyCPD, ...
            gratingSceneParams);

    % Store test scene engines
    theTestSceneEngines = cell(1,size(theLMSconeContrastDirections,2));

    % Create the test scene engines
    for iDir = 1:size(theLMSconeContrastDirections,2)
        % Create a grating scene engine with the examined chromatic direction
        theTestSceneEngines{iDir} = createGratingSceneEngine(...
            theLMSconeContrastDirections(:,iDir), ...
            examinedSpatialFrequencyCPD, ...
            gratingSceneParams);
    end % iDir

end




function [figureFileBaseDir, resultsFileBaseDir] = setupFigureDirectory(theScriptName, ...
    useMetaContrast,useConeContrast,useFixationalEMs,...
    whichNoiseFreeNre,whichNoisyInstanceNre,...
    whichClassifierEngine,mRGCOutputSignalType)

    % Make sure local/figures directory exists so we can write out our figures in peace
    projectBaseDir = ISETBioCSFGeneratorRootPath;
    if (~exist(fullfile(projectBaseDir,'local', theScriptName, 'figures'),'dir'))
        mkdir(fullfile(projectBaseDir,'local', theScriptName, 'figures'));
        fprintf('Generated figure directory at %s\n', fullfile(projectBaseDir,'local', theScriptName, 'figures'))
    end

    if (~exist(fullfile(projectBaseDir,'local', theScriptName, 'results'),'dir'))
        mkdir(fullfile(projectBaseDir,'local', theScriptName, 'results'));
        fprintf('Generated results directory at %s\n', fullfile(projectBaseDir,'local', theScriptName, 'figures'))
    end


    figureFileBaseDir = fullfile(projectBaseDir,'local',mfilename,'figures', ...
        sprintf('%s_Meta_%d_ConeContrast_%d_FEMs_%d_%s_%s_%s_%s', mfilename, ...
        useMetaContrast,useConeContrast,useFixationalEMs,whichNoiseFreeNre,whichNoisyInstanceNre,...
        whichClassifierEngine,mRGCOutputSignalType));

    if (~exist(figureFileBaseDir, 'dir'))
        mkdir(figureFileBaseDir);
        fprintf('Generated figure sub-directory at %s\n', figureFileBaseDir);
    end

    resultsFileBaseDir = fullfile(projectBaseDir,'local',mfilename,'results', ...
        sprintf('%s_Meta_%d_ConeContrast_%d_FEMs_%d_%s_%s_%s_%s', mfilename, ...
        useMetaContrast,useConeContrast,useFixationalEMs,whichNoiseFreeNre,whichNoisyInstanceNre,...
        whichClassifierEngine,mRGCOutputSignalType));
    
    if (~exist(resultsFileBaseDir, 'dir'))
        mkdir(resultsFileBaseDir);
        fprintf('Generated figure sub-directory at %s\n', resultsFileBaseDir);
    end



end

