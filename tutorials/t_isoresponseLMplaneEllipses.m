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

%{
% Examples of contours in cone contrast space that
% we can compare with measurements.

% Step 2. Integrate data from all Pts_coneContrasts
% plotMultipleIsoresponseLMplaneEllipseThresholds(Pts_coneContrast)

% Step3. Plot against human data
% compareHumanToMRGCmosaicEllipses()

% Define the reference contrasts
Pts_coneContrast = [
   -0.3000   -0.3000         0; ...
   -0.1500   -0.1500    0.0000; ...
    0.0000    0.0000         0; ...
    0.1500    0.1500         0; ...
    0.3000    0.3000         0; ...
   -0.0500    0.0500         0; ...
    0.0500   -0.0500    0.0000];


% MRGCs (non-linear) and simulating ON/OFF mosaic (odd-numbered mRGCs)

mRGCOutputNonLinearityParams = struct(...
    'visualizeNonLinearActivationFunction', true, ...
    'sourceSignal', 'compositeResponse', ...    % where to apply the non-linerity. choose from {'centerComponentResponse', 'surroundComponentResponse', 'compositeResponse'}
    'type', 'Naka Rushton', ...
    'params', struct(...
        'rectification', 'half', ...            % apply non-linearity to both response polarities % choose between {'half', 'full', 'none'}
        'n', 2.0, ...                           % exponent
        's', 1.0, ...                           % super-saturation exponent (super-saturating if s > 1)
        'c50', 0.05, ...                         % semi-saturation response (in normalized range, based on the passed maxLinearResponse)
        'bias', 0.0, ...                        % bias (in normalized range, based on the passed maxLinearResponse)
        'gain', 1.0, ...                        % post-non-linearity gain
        'maxLinearResponse', 1.0) ...           % max absolute value that the linear excitations-based response can have (clipping after that)
    );

t_isoresponseLMplaneEllipses(...
        'useMetaContrast', true, ...
        'whichNoiseFreeNre', 'mRGCMosaic', ...
        'mosaicEccDegs', [0.0 0.0], ...
        'mosaicSizeDegs', [1 1], ...
        'synthesizedRGCmosaicName', 'PLOSpaperFovealMosaic', ...
        'visualizeCroppedRGCmosaicRelationshipToSynthesizedMosaic', true, ...
        'meanLuminanceCdPerM2', 2000, ...
        'backgroundLMSconeExcitations', [0.1529 0.1324 0.0828], ...
        'referenceLMSconeContrast', Pts_coneContrast, ...
        'contrastCorrection', 'accountForTotalContrastBeing_refC_PLUS_deltaC', ...
        'spatialFrequency', 0.0, ...
        'spatialPhaseDegs', 0.0, ...
        'stimSizeDegs', 1.0, ...
        'spatialEnvelope', 'rect', ...
        'spatialEnvelopeRadiusFraction', 0.95, ...  % for debugging
        'presentationMode', 'static', ...
        'nTest', 1024, ...
        'psychometricCurveSamplesNum', 5, ...
        'examinedDirectionsOnLMplane', 0:22.5:(360-22.5), ...
        'visualizeEachScene', ~true, ...
        'simulateONOFFmosaic', true, ...
        'mRGCOutputNonLinearityParams', mRGCOutputNonLinearityParams, ...
        'mRGCOutputSignalType', 'mRGCs', ...  
        'whichNoisyInstanceNre', 'Gaussian', ...
        'gaussianSigma', 0.1, ...
        'whichClassifierEngine', 'rceTemplateDistance', ...
        'useConeContrast', true, ...
        'thresholdPara', struct( ...
            'logThreshLimitLow', 4.0, ...
            'logThreshLimitHigh', 0.0, ...
            'logThreshLimitDelta', 0.02, ...
            'slopeRangeLow', 1/20, ...
            'slopeRangeHigh', 200/20, ...
            'slopeDelta', 2.5/50, ...
            'thresholdCriterion', 0.81606, ...
            'guessRate', 1/2, ...
            'lapseRate', 0) ...
        );



% MRGCs (linear, based on cone contrast)
t_isoresponseLMplaneEllipses(...
        'useMetaContrast', true, ...
        'whichNoiseFreeNre', 'mRGCMosaic', ...
        'mosaicEccDegs', [0.0 0.0], ...
        'mosaicSizeDegs', [1 1], ...
        'synthesizedRGCmosaicName', 'PLOSpaperFovealMosaic', ...
        'visualizeCroppedRGCmosaicRelationshipToSynthesizedMosaic', true, ...
        'meanLuminanceCdPerM2', 2000, ...
        'backgroundLMSconeExcitations', [0.1529 0.1324 0.0828], ...
        'referenceLMSconeContrast', Pts_coneContrast, ...
        'contrastCorrection', 'accountForTotalContrastBeing_refC_PLUS_deltaC', ...
        'spatialFrequency', 0.0, ...
        'spatialPhaseDegs', 0.0, ...
        'stimSizeDegs', 1.0, ...
        'spatialEnvelope', 'rect', ...
        'spatialEnvelopeRadiusFraction', 0.95, ...  % for debugging
        'presentationMode', 'static', ...
        'nTest', 1024, ...
        'psychometricCurveSamplesNum', 5, ...
        'examinedDirectionsOnLMplane', 0:22.5:(360-22.5), ...
        'visualizeEachScene', ~true, ...
        'simulateONOFFmosaic', false, ...
        'mRGCOutputSignalType', 'mRGCs', ...
        'whichNoisyInstanceNre', 'Gaussian', ...
        'gaussianSigma', 0.1, ...
        'whichClassifierEngine', 'rceTemplateDistance', ...
        'useConeContrast', true);

% CONES (input cone mosaic of MRC mosaic, WITH cone contrast and Gaussian
noise with same sigma as mRGCs)

t_isoresponseLMplaneEllipses(...
        'useMetaContrast', true, ...
        'whichNoiseFreeNre', 'mRGCMosaic', ...
        'mosaicEccDegs', [0.0 0.0], ...
        'mosaicSizeDegs', [1 1], ...
        'synthesizedRGCmosaicName', 'PLOSpaperFovealMosaic', ...
        'visualizeCroppedRGCmosaicRelationshipToSynthesizedMosaic', true, ...
        'meanLuminanceCdPerM2', 2000, ...
        'backgroundLMSconeExcitations', [0.1529 0.1324 0.0828], ...
        'referenceLMSconeContrast', Pts_coneContrast, ...
        'contrastCorrection', 'accountForTotalContrastBeing_refC_PLUS_deltaC', ...
        'spatialFrequency', 0.0, ...
        'spatialPhaseDegs', 0.0, ...
        'stimSizeDegs', 1.0, ...
        'spatialEnvelope', 'rect', ...
        'spatialEnvelopeRadiusFraction', 0.95, ...  % for debugging
        'presentationMode', 'static', ...
        'nTest', 1024, ...
        'psychometricCurveSamplesNum', 5, ...
        'employMosaicSpecificConeFundamentals', true, ...
        'examinedDirectionsOnLMplane', 0:22.5:(360-22.5), ...
        'visualizeEachScene', ~true, ...
        'simulateONOFFmosaic', false, ...
        'mRGCOutputSignalType', 'cones', ...
        'whichNoisyInstanceNre', 'Gaussian', ...
        'gaussianSigma', 0.1, ...
        'whichClassifierEngine', 'rceTemplateDistance', ...
        'useConeContrast', true, ...
        'thresholdPara', struct( ...
            'logThreshLimitLow', 3.0, ...
            'logThreshLimitHigh', 1.0, ...
            'logThreshLimitDelta', 0.02, ...
            'slopeRangeLow', 1/20, ...
            'slopeRangeHigh', 200/20, ...
            'slopeDelta', 2.5/50, ...
            'thresholdCriterion', 0.81606, ...
            'guessRate', 1/2, ...
            'lapseRate', 0) ...
        );

% CONES (input cone mosaic of MRC mosaic, WITHOUT cone contrast and Gaussian
noise with same sigma as mRGCs)

t_isoresponseLMplaneEllipses(...
        'useMetaContrast', true, ...
        'whichNoiseFreeNre', 'mRGCMosaic', ...
        'mosaicEccDegs', [0.0 0.0], ...
        'mosaicSizeDegs', [1 1], ...
        'synthesizedRGCmosaicName', 'PLOSpaperFovealMosaic', ...
        'visualizeCroppedRGCmosaicRelationshipToSynthesizedMosaic', true, ...
        'meanLuminanceCdPerM2', 2000, ...
        'backgroundLMSconeExcitations', [0.1529 0.1324 0.0828], ...
        'referenceLMSconeContrast', Pts_coneContrast, ...
        'contrastCorrection', 'accountForTotalContrastBeing_refC_PLUS_deltaC', ...
        'spatialFrequency', 0.0, ...
        'spatialPhaseDegs', 0.0, ...
        'stimSizeDegs', 1.0, ...
        'spatialEnvelope', 'rect', ...
        'spatialEnvelopeRadiusFraction', 0.95, ...  % for debugging
        'presentationMode', 'static', ...
        'nTest', 1024, ...
        'psychometricCurveSamplesNum', 5, ...
        'employMosaicSpecificConeFundamentals', true, ...
        'examinedDirectionsOnLMplane', 0:22.5:(360-22.5), ...
        'visualizeEachScene', ~true, ...
        'simulateONOFFmosaic', false, ...
        'mRGCOutputSignalType', 'cones', ...
        'whichNoisyInstanceNre', 'Gaussian', ...
        'gaussianSigma', 500, ...
        'whichClassifierEngine', 'rceTemplateDistance', ...
        'useConeContrast', false, ...
        'thresholdPara', struct( ...
            'logThreshLimitLow', 3.0, ...
            'logThreshLimitHigh', 1.0, ...
            'logThreshLimitDelta', 0.02, ...
            'slopeRangeLow', 1/20, ...
            'slopeRangeHigh', 200/20, ...
            'slopeDelta', 2.5/50, ...
            'thresholdCriterion', 0.81606, ...
            'guessRate', 1/2, ...
            'lapseRate', 0) ...
        );

% CONES (input cone mosaic of MRC mosaic, but now Poisson noise)
t_isoresponseLMplaneEllipses(...
        'useMetaContrast', true, ...
        'whichNoiseFreeNre', 'mRGCMosaic', ...
        'mosaicEccDegs', [0.0 0.0], ...
        'mosaicSizeDegs', [1 1], ...
        'synthesizedRGCmosaicName', 'PLOSpaperFovealMosaic', ...
        'visualizeCroppedRGCmosaicRelationshipToSynthesizedMosaic', true, ...
        'meanLuminanceCdPerM2', 2000, ...
        'backgroundLMSconeExcitations', [0.1529 0.1324 0.0828], ...
        'referenceLMSconeContrast', Pts_coneContrast, ...
        'contrastCorrection', 'accountForTotalContrastBeing_refC_PLUS_deltaC', ...
        'spatialFrequency', 0.0, ...
        'spatialPhaseDegs', 0.0, ...
        'stimSizeDegs', 1.0, ...
        'presentationMode', 'static', ...
        'nTest', 1024, ...
        'psychometricCurveSamplesNum', 5, ...
        'examinedDirectionsOnLMplane', 0:22.5:(360-22.5), ...
        'visualizeEachScene', true, ...
        'simulateONOFFmosaic', false, ...
        'mRGCOutputSignalType', 'cones', ...
        'whichNoisyInstanceNre', 'Poisson', ...
        'whichClassifierEngine', 'rcePoisson', ...
        'useConeContrast', false, ...
        'thresholdPara', struct( ...
            'logThreshLimitLow', 4.0, ...
            'logThreshLimitHigh', 1.0, ...
            'logThreshLimitDelta', 0.02, ...
            'slopeRangeLow', 1/20, ...
            'slopeRangeHigh', 200/20, ...
            'slopeDelta', 2.5/50, ...
            'thresholdCriterion', 0.81606, ...
            'guessRate', 1/2, ...
            'lapseRate', 0) ...
       );
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
    %
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
    options.referenceLMSconeContrast (:,3) double = [0 0 0];

    % If we have a non-zero referenceLMSconeContrast, a contrast correction is needed
    % to have the sceneEngine produce a total contrast that is equal to reference contrast + delta contrast
    options.contrastCorrection (1,:) char {mustBeMember(options.contrastCorrection,{'accountForTotalContrastBeing_refC_PLUS_deltaC',''})} = '';


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
    options.spatialEnvelope (1,:) char {mustBeMember(options.spatialEnvelope,{'disk','rect'})} = 'rect';
    options.spatialEnvelopeRadiusFraction (1,1) double = 1.0

    % Presentation mode
    options.presentationMode (1,:) char = 'static';

    % Apply temporal filter?
    %
    % The timebase of the filter is assumed to match the frame rate, so we
    % only need to specify a list of filter values.  Since these canON
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


    options.parPoolSize (1,:) char = 'default'

    % mRGC nonLinearities
    options.mRGCOutputNonLinearityParams = [];

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
spatialEnvelope = options.spatialEnvelope;
spatialEnvelopeRadiusFraction = options.spatialEnvelopeRadiusFraction;

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

parPoolSize = options.parPoolSize;

simulateONOFFmosaic = options.simulateONOFFmosaic;
mRGCOutputNonLinearityParams = options.mRGCOutputNonLinearityParams;

if (~isempty(mRGCOutputNonLinearityParams))
    if (useMetaContrast)
        fprintf('\n***********\nSince an mRGCOutputNonLinearityParams is passed,\nwe override useMetaContrast, which was set to true (and which can be used for linear responses only),\nand set it to false !!! \n*******\n');
        pause(1.0);
    end
    useMetaContrast = false;
end


%% Freeze rng for replicatbility and validation
rng(1);

%% Close any stray figs
hAllFigs = findall(groot,'Type','figure');

% Close all figures
for i = 1:numel(hAllFigs)
    set(hAllFigs(i), 'HandleVisibility', 'on')
end
close all;


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

%% Spatial frequency & orientation
examinedSpatialFrequencyCPD = options.spatialFrequency;
orientationDegs = options.orientationDegs;

debugStimulusConfig = ~true;

%% List of LM angles to be tested
examinedDirectionsOnLMplane = options.examinedDirectionsOnLMplane;

%% Max RMS contrast (so as to keep stimuli within the display gamut)
rmsLMconeContrast = 0.07;

%% List of reference contrasts
referenceLMSconeContrasts = options.referenceLMSconeContrast;

%% Contrast correction
contrastCorrection = options.contrastCorrection;



for iRefContrast = 1:size(referenceLMSconeContrasts,1)

    % The examined reference contrast
    referenceLMSconeContrast = referenceLMSconeContrasts(iRefContrast,:);

    fprintf('Running simulation for reference contrast %d of %d\n', iRefContrast, size(referenceLMSconeContrasts,1));
    
    % Thresholds filename
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

    if (debugStimulusConfig)
        rmsLMconeContrast = 1.0;
        thresholdPara.logThreshLimitLow = 0.1;
    end

    % Compute directions based on rmsLMconeContrast.
    [theDeltaLMSconeContrastDirections, examinedDirectionsOnLMplane] = ...
        computeLMSconeContrastDirections(rmsLMconeContrast, examinedDirectionsOnLMplane);

    % Add the reference LMScone contrast. 
    % For contrasts symmetric around the background,
    % referenceLMSconeContrast should be set to [0 0 0]
    theLMSconeContrastDirections = bsxfun(@plus, theDeltaLMSconeContrastDirections, referenceLMSconeContrast(:));

    % Set grating engine parameters
    gratingSceneParams = struct( ...
        'meanLuminanceCdPerM2', meanLuminanceCdPerM2, ...
        'spectralSupport', 400:20:750, ...
        'fovDegs', stimSizeDegs, ...
        'pixelsNum', pixelsNum, ...
        'spatialPhase', spatialPhaseDegs, ...
        'spatialEnvelope', spatialEnvelope, ...
        'spatialEnvelopeRadiusDegs', spatialEnvelopeRadiusFraction * 0.5* stimSizeDegs, ...
        'orientation', orientationDegs, ...
        'presentationMode', 'flashed', ...
        'duration', frameDurationSeconds, ...
        'frameDurationSeconds', frameDurationSeconds);

    if (isempty(backgroundLMSconeExcitations))
        gratingSceneParams.meanChromaticityXY =  meanChromaticityXY;
    else
        gratingSceneParams.backgroundLMSconeExcitations = backgroundLMSconeExcitations;
    end

    if (~strcmp(presentationMode, 'static'))

        % Non-static. Add temporal parameters
        stimulusDuration = stimDurationTemporalCycles*frameDurationSeconds;
        
        temporalModulationParams.stimOnFrameIndices = [];
        temporalModulationParams.stimDurationFramesNum = [];
        temporalModulationParams.phaseDirection = 1;
        temporalModulationParams.stimDurationTemporalCycles = stimDurationTemporalCycles;
        temporalModulationParams.temporalFrequencyHz = temporalFrequencyHz;
    
        gratingSceneParams.temporalModulationParams = temporalModulationParams;
    end


    % Which stimuli to skip
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
        [theBackgroundSceneEngine, theNullSceneEngine, theTestSceneEngines] = configureStimulusSceneEngines(...
            theDeltaLMSconeContrastDirections, referenceLMSconeContrast, ...
            examinedSpatialFrequencyCPD, ...
            gratingSceneParams, contrastCorrection);
    
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
        
        
                % NON-LINEARITY
                if (~isempty(mRGCOutputNonLinearityParams))
                    nreNoiseFreeParams.mRGCMosaicParams.nonLinearitiesList{1} = mRGCOutputNonLinearityParams;
                else
                    nreNoiseFreeParams.mRGCMosaicParams.nonLinearitiesList = [];
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



    if (iRefContrast == 1)
        mRGCLinearityString = '';
        if strcmp(whichNoiseFreeNre, 'mRGCMosaic')
            if (~isempty(mRGCOutputNonLinearityParams))

                if (nreNoiseFreeParams.mRGCMosaicParams.simulateONOFFmosaic)
                    mRGCLinearityString = sprintf('ON_OFF_Emulation_nLinearRect_''%s''_c50_%1.2f_n%1.1f_s%1.1f', ...
                    mRGCOutputNonLinearityParams.params.rectification, ...
                    mRGCOutputNonLinearityParams.params.c50, ...
                    mRGCOutputNonLinearityParams.params.n, ...
                    mRGCOutputNonLinearityParams.params.s);
                else
                    mRGCLinearityString = sprintf('nonLinearRect_''%s''_c50_%1.2f_n%1.1f_s%1.1f', ...
                        mRGCOutputNonLinearityParams.params.rectification, ...
                        mRGCOutputNonLinearityParams.params.c50, ...
                        mRGCOutputNonLinearityParams.params.n, ...
                        mRGCOutputNonLinearityParams.params.s);
                end

            else
                mRGCLinearityString = 'Linear';
            end
        end % if strcmp(whichNoiseFreeNre, 'mRGCMosaic')

        %% Make sure local/figures directory exists so we can write out our figures in peace
        [figureFileBaseDir, resultsFileBaseDir] = setupFigureDirectory(mfilename, ...
            useMetaContrast, useConeContrast, useFixationalEMs, ...
            whichNoiseFreeNre, whichNoisyInstanceNre,...
            whichClassifierEngine, mRGCOutputSignalType, mRGCLinearityString);
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
            error('Unsupported noisy instances neural response engine: ''%s''.', whichNoisyInstanceNre);
    end % switch (whichNoisyInstanceNre)


    % If we use cone contrast, we will neeed a null scene for normalization.
    if (useConeContrast) || (ischar(temporalFilterValues))
        nreNoiseFreeParams.nullStimulusSceneSequence = theBackgroundSceneEngine.compute(0.0);
    end

    % Create the neural engine
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
        if (isempty(referenceLMSconeContrast))
            metaSceneEngineParams.sceneEngineName = ...
                sprintf('TEST stimulus meta-engine (reference contrast = 0 0 0)');
        else
            metaSceneEngineParams.sceneEngineName = ...
            sprintf('TEST stimulus meta-engine (reference contrast: %2.3f %2.3f %2.3f)', ...
                referenceLMSconeContrast(1), referenceLMSconeContrast(2), referenceLMSconeContrast(3));
        end
        theMetaSceneEngine = sceneEngine(@sceMetaContrast,metaSceneEngineParams);
    
        metaNullSceneEngineParams = sceMetaContrast;
        if (isempty(referenceLMSconeContrast))
            metaNullSceneEngineParams.sceneEngineName = ...
                sprintf('NULL stimulus meta-engine (reference contrast = 0 0 0)');
        else
            metaNullSceneEngineParams.sceneEngineName = ...
                sprintf('NULL stimulus meta engine (reference contrast: %2.3f %2.3f %2.3f)', ...
                    referenceLMSconeContrast(1), referenceLMSconeContrast(2), referenceLMSconeContrast(3));
        end
        theNullSceneMetaEngine = sceneEngine(@sceMetaContrast,metaNullSceneEngineParams);

        % Create nreMetaContrast using the actual scene and neural engines
        metaNeuralResponseEngineNoiseFreeParams = nreNoiseFreeMetaContrast;
        metaNeuralResponseEngineNoiseFreeParams.contrast0 = 0;
        metaNeuralResponseEngineNoiseFreeParams.contrast1 = 1;
        metaNeuralResponseEngineNoiseFreeParams.neuralEngine = theNeuralEngine;
    
        metaNeuralResponseEngineNoisyInstanceParams = nreNoisyInstancesMetaContrast;
        metaNeuralResponseEngineNoisyInstanceParams.neuralEngine = theNeuralEngine;
    end


    % Instantiate the responseClassifierEngine
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
    end % switch (whichClassifierEngine)

    % Parameters for threshold estimation/quest engine
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


    % Setup figures
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


    % Compute threshold for each chromatic angle
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
        theNullSceneEngine.visualizeEachCompute = visualizeEachScene;

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
                    'TAFC', true, ...
                    'nullSceneEngineToEmploy', theNullSceneMetaEngine, ...
                    'useMetaContrast', useMetaContrast, ...
                    'trainFixationalEM', trainFixationalEMObj, ...
                    'testFixationalEM', testFixationalEMObj, ...
                    'visualizeStimulus', visualizeEachScene, ...
                    'verbose', verbose, ...
                    'maxVisualizedNoisyResponseInstances', maxVisualizedNoisyResponseInstances);
         else
            % Compute the threshold for our grating scene with the previously
            % defined neural and classifier engine.  This function does a lot
            % of work, see the function itself, as well as function
            % computePerformance.
            [logThreshold(iChromaDirection), questObj, ~, psychometricCurveParams(iChromaDirection,:)] = ...
                    computeThreshold(theSceneEngine, theNeuralEngine, theClassifierEngine, ...
                    classifierEngineParams, thresholdPara, questEnginePara, ...
                    'TAFC', true, ...
                    'nullSceneEngineToEmploy', theNullSceneEngine, ...
                    'useMetaContrast', useMetaContrast, ...
                    'trainFixationalEM', trainFixationalEMObj, ...
                    'testFixationalEM', testFixationalEMObj, ...
                    'visualizeStimulus', visualizeEachScene, ...
                    'verbose', verbose, ...
                    'maxVisualizedNoisyResponseInstances', maxVisualizedNoisyResponseInstances);
         end % useMetaContrast
        

        % Call computeThreshold one last time with threshold contrast to get the non-linearity metaData at
        % the thresholdContrat
        thresholdParaFinal = thresholdPara;
        thresholdParaFinal.logThreshLimitLow = -logThreshold(iChromaDirection);
        thresholdParaFinal.logThreshLimitHigh = -logThreshold(iChromaDirection);
        questEngineParaFinal = questEnginePara;
        questEngineParaFinal.minTrial = nTest;
        questEngineParaFinal.maxTrial = nTest;

        [~, ~, ~, ~, ~,~,~, trialByTrialWhichMetaData] = ...
            computeThreshold(theSceneEngine, theNeuralEngine, theClassifierEngine, ...
                    classifierEngineParams, thresholdParaFinal, questEngineParaFinal, ...
                    'TAFC', true, ...
                    'nullSceneEngineToEmploy', theNullSceneEngine, ...
                    'useMetaContrast', false, ...
                    'trainFixationalEM', trainFixationalEMObj, ...
                    'testFixationalEM', testFixationalEMObj, ...
                    'visualizeStimulus', visualizeEachScene, ...
                    'verbose', true, ...
                    'visualizeAllComponents', true, ...
                    'maxVisualizedNoisyResponseInstances', maxVisualizedNoisyResponseInstances);

        % Print the non-linearity figure at the threshold contrast
        theKeys = keys(trialByTrialWhichMetaData);
        for i = 1:numel(theKeys)
            ww = trialByTrialWhichMetaData(theKeys{i});
            for jj = 1:numel(ww)
                theMetaData = ww{jj};
                if (isstruct(theMetaData)) && (isfield(theMetaData, 'mRGCnonLinearityFigureHandle'))
                    hFig = figure(theMetaData.mRGCnonLinearityFigureHandle);
                    thePDFfileName = sprintf('NonLinearity_ReferencePos_%d_ChromaDirection_%d_''%s''.pdf', iRefContrast, iChromaDirection, theMetaData.conditionLabel);
                    NicePlot.exportFigToPDF(fullfile(figureFileBaseDir,thePDFfileName), hFig, 300, 'beVerbose');
                    delete(hFig);
                end
            end

        end

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
    
    % The rms cone contrast with respect to the reference
    stimulusRMSLMconeContrast = sqrt(sum(theDeltaLMSconeContrastDirections.^2,1));
    
    % Save thresholds
    fprintf('Saving to %s\n', fullfile(resultsFileBaseDir,thresholdsDataFileName));

    save(fullfile(resultsFileBaseDir,thresholdsDataFileName), ...
        'options', 'theLMSconeContrastDirections', 'theDeltaLMSconeContrastDirections', ...
        'stimulusRMSLMconeContrast', 'examinedDirectionsOnLMplane', 'thresholdContrasts');
    
    % Figure for plotting the thresholds on the LM plane
    % along with the stimuli
    hFigStimuliAndThresholds = figure(2346); clf;
    set(hFigStimuliAndThresholds, 'Color', [1 1 1], 'Position', [10 10 1200 1200]);
    set(hFigStimuliAndThresholds, 'HandleVisibility', 'off');
    
    exportFig = true;
    initializeFig = true;
    theThresholdAxes = []; 
    maxVisualizedThreshold = 0.5;
    figName = sprintf('PsychometricFunctionsReferenceLMScontrast_%2.2f_%2.2f_%2.2f', 100*referenceLMSconeContrast(1), 100*referenceLMSconeContrast(2), 100*referenceLMSconeContrast(3));
   
   % Plot the data
   visualizeIsoThresholdEllipsesOnLMplane(...
            theLMSconeContrastDirections, ...
            theDeltaLMSconeContrastDirections, ...
            stimulusRMSLMconeContrast, ...
            examinedSpatialFrequencyCPD, ...
            gratingSceneParams, ...
            examinedDirectionsOnLMplane, ...
            skippedDirections, ...
            thresholdContrasts, ...
            maxVisualizedThreshold, ...
            referenceLMSconeContrast, ...
            figureFileBaseDir, hFigStimuliAndThresholds, ...
            exportFig, initializeFig, theThresholdAxes, figName);

    % Restore handle visibility so that a future close will close the figures
    hAllFigs = findall(groot,'Type','figure');
    for i = 1:numel(hAllFigs)
        set(hAllFigs(i), 'HandleVisibility', 'on')
    end

    % Do a check on the answer
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
end % for iRefContras

end


%
% SUPPORTING FUNCTIONS
%


function [theBackgroundSceneEngine, theNullSceneEngine, theTestSceneEngines] = configureStimulusSceneEngines(...
    theDeltaLMSconeContrastDirections, referenceLMSconeContrast, examinedSpatialFrequencyCPD, ...
    gratingSceneParams, contrastCorrection)
 
    % Create the NULL scene engine
    gratingSceneParamsForBackgroundScene = gratingSceneParams;
    gratingSceneParamsForBackgroundScene.sceneEngineName = sprintf('background stimulus engine');

    % Create the BKGND stimulus engine
    LMSconeContrastDirection = [0 0 0];

    theBackgroundSceneEngine = createGratingSceneEngine(...
            LMSconeContrastDirection, ...
            examinedSpatialFrequencyCPD, ...
            gratingSceneParamsForBackgroundScene);

    % Create the NULL scene engine
    gratingSceneParamsForNullScene = gratingSceneParams;

    % Add a name for the null stimulus engine
    if (isempty(referenceLMSconeContrast))
        gratingSceneParamsForNullScene.sceneEngineName = ...
            sprintf('null stimulus engine (reference contrast = 0 0 0)');
    else
        gratingSceneParamsForNullScene.sceneEngineName = ...
            sprintf('null stimulus engine (reference contrast: %2.3f %2.3f %2.3f)', ...
                referenceLMSconeContrast(1), referenceLMSconeContrast(2), referenceLMSconeContrast(3));

        % Adjust the background LMS excitations if we
        if (~isempty(gratingSceneParamsForNullScene.backgroundLMSconeExcitations))
            fprintf('\n\n>>>> Null scene background LMS before adjusting for reference contrast (%2.3f %2.3f %2.3f) were: (%2.3f %2.3f %2.3f)',...
                referenceLMSconeContrast(1), referenceLMSconeContrast(2), referenceLMSconeContrast(3), ...
                gratingSceneParamsForNullScene.backgroundLMSconeExcitations(1), ...
                gratingSceneParamsForNullScene.backgroundLMSconeExcitations(2), ...
                gratingSceneParamsForNullScene.backgroundLMSconeExcitations(3));

            gratingSceneParamsForNullScene.backgroundLMSconeExcitations = ...
                gratingSceneParamsForNullScene.backgroundLMSconeExcitations .* (1+referenceLMSconeContrast);

            fprintf('\n>>>> Null scene background LMS after adjusting for reference contrast (%2.3f %2.3f %2.3f) are: (%2.3f %2.3f %2.3f)\n\n\n',...
                referenceLMSconeContrast(1), referenceLMSconeContrast(2), referenceLMSconeContrast(3), ...
                gratingSceneParamsForNullScene.backgroundLMSconeExcitations(1), ...
                gratingSceneParamsForNullScene.backgroundLMSconeExcitations(2), ...
                gratingSceneParamsForNullScene.backgroundLMSconeExcitations(3));
        end
    end

    % Create the null stimulus engine
    LMSconeContrastDirection = [0 0 0];
    theNullSceneEngine = createGratingSceneEngine(...
            LMSconeContrastDirection, ...
            examinedSpatialFrequencyCPD, ...
            gratingSceneParamsForNullScene);

    % Store test scene engines
    theTestSceneEngines = cell(1,size(theDeltaLMSconeContrastDirections,2));

    gratingSceneParams.backgroundLMSconeExcitations = ...
                gratingSceneParams.backgroundLMSconeExcitations .* (1+referenceLMSconeContrast);


    % Create the test scene engines
    for iDir = 1:size(theDeltaLMSconeContrastDirections,2)
        % Add a name to each test stimulus engine
        gratingSceneParams.sceneEngineName = sprintf('test stimulus engine for delta contrast direction: %2.3f %2.3f %2.3f)', ...
            theDeltaLMSconeContrastDirections(1,iDir), theDeltaLMSconeContrastDirections(2,iDir), theDeltaLMSconeContrastDirections(3,iDir));

        % Correct the contrast so that we get (1+refC+deltaC)
        % totalLMSdesired = backLMS * (1+refC+deltaC)
        % totalLMSprovidedBySceneEngine = backLMS * (1+refC) .* (1+deltaC)
        % totalLMSprovidedBySceneEngineWithCorrection = backLMS * (1+refC) .* (1+ deltaC-correctionC ) =
        %                                             = backLMS * (1 + refC + deltaC-correctionC + (deltaC-correctionC)*refC)
        %
        % Equating totalLMSdesired == totalLMSprovidedBySceneEngineWithCorrection (to compute correctionC)               \leadsto
        % backLMS * (1+refC+deltaC) = backLMS * (1 + refC + deltaC-correctionC + (deltaC-correctionC)*refC)              \leadsto
        % backLMS * (1+refC+deltaC) = backLMS * (1 + refC + deltaC) + backLMS*(-correctionC + (deltaC-correctionC)*refC) \leadsto
        % 1+refC+deltaC = 1 + refC + deltaC - correctionC + (deltaC-correctionC)*refC \leadsto
        % -correctionC + (deltaC-correctionC)*refC = 0      \leadsto
        % -correctionC -correctionC*refC = -deltaC*refC     \leadsto
        % correctionC * (1+refC) = deltaC*refC              \leadsto
        %
        %
        %  correctionC = deltaC*refC/(1+refC);
        %
        %
        % Verifying:
        % totalLMSprovidedBySceneEngineWithCorrection = backLMS * (1+refC) .* (1+ deltaC-correctionC )
        % totalLMSprovidedBySceneEngineWithCorrection = backLMS * (1+refC) .* (1+ deltaC-deltaC*refC/(1+refC))
        % totalLMSprovidedBySceneEngineWithCorrection = backLMS * (1+refC) .* (1+ deltaC)- backLMS * (1+refC) * deltaC*refC/(1+refC))
        % totalLMSprovidedBySceneEngineWithCorrection = backLMS * (1+refC) .* (1+ deltaC)- backLMS * deltaC*refC
        % totalLMSprovidedBySceneEngineWithCorrection = backLMS * (1+refC) + backLMS * (1+refC)*deltaC - backLMS * deltaC*refC
        % totalLMSprovidedBySceneEngineWithCorrection = backLMS * (1+refC) + backLMS * deltaC
        % totalLMSprovidedBySceneEngineWithCorrection = backLMS * (1+refC + deltaC)

        switch (contrastCorrection)
            case 'accountForTotalContrastBeing_refC_PLUS_deltaC'
                correctionC = theDeltaLMSconeContrastDirections(:,iDir) .* referenceLMSconeContrast(:) ./ (1+referenceLMSconeContrast(:));
                assert(prod(size(correctionC)) == numel(referenceLMSconeContrast), 'Need to reshape the referenceLMSconeContrast');
            otherwise
                correctionC = [0 0 0]';
        end


        % Create a grating scene engine with the examined chromatic direction
        theTestSceneEngines{iDir} = createGratingSceneEngine(...
            theDeltaLMSconeContrastDirections(:,iDir) - correctionC, ...
            examinedSpatialFrequencyCPD, ...
            gratingSceneParams);
    end % iDir

end




function [figureFileBaseDir, resultsFileBaseDir] = setupFigureDirectory(theScriptName, ...
    useMetaContrast,useConeContrast,useFixationalEMs,...
    whichNoiseFreeNre,whichNoisyInstanceNre,...
    whichClassifierEngine,mRGCOutputSignalType, mRGCLinearityString)

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


    if (~isempty(mRGCLinearityString))
        figureFileBaseDir = fullfile(projectBaseDir,'local',mfilename,'figures', ...
            sprintf('%s_Meta_%d_ConeContrast_%d_FEMs_%d_%s_%s_%s_%s_%s', mfilename, ...
            useMetaContrast,useConeContrast,useFixationalEMs,whichNoiseFreeNre,whichNoisyInstanceNre,...
            whichClassifierEngine,mRGCOutputSignalType, mRGCLinearityString));
    else
        figureFileBaseDir = fullfile(projectBaseDir,'local',mfilename,'figures', ...
            sprintf('%s_Meta_%d_ConeContrast_%d_FEMs_%d_%s_%s_%s_%s', mfilename, ...
            useMetaContrast,useConeContrast,useFixationalEMs,whichNoiseFreeNre,whichNoisyInstanceNre,...
            whichClassifierEngine,mRGCOutputSignalType));
    end

    if (~exist(figureFileBaseDir, 'dir'))
        mkdir(figureFileBaseDir);
        fprintf('Generated figure sub-directory at %s\n', figureFileBaseDir);
    end

    if (~isempty(mRGCLinearityString))
        resultsFileBaseDir = fullfile(projectBaseDir,'local',mfilename,'results', ...
            sprintf('%s_Meta_%d_ConeContrast_%d_FEMs_%d_%s_%s_%s_%s_%s', mfilename, ...
            useMetaContrast,useConeContrast,useFixationalEMs,whichNoiseFreeNre,whichNoisyInstanceNre,...
            whichClassifierEngine,mRGCOutputSignalType, mRGCLinearityString));
    else
        resultsFileBaseDir = fullfile(projectBaseDir,'local',mfilename,'results', ...
            sprintf('%s_Meta_%d_ConeContrast_%d_FEMs_%d_%s_%s_%s_%s', mfilename, ...
            useMetaContrast,useConeContrast,useFixationalEMs,whichNoiseFreeNre,whichNoisyInstanceNre,...
            whichClassifierEngine,mRGCOutputSignalType));
    end

    
    if (~exist(resultsFileBaseDir, 'dir'))
        mkdir(resultsFileBaseDir);
        fprintf('Generated figure sub-directory at %s\n', resultsFileBaseDir);
    end



end

