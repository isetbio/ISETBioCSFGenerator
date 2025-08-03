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

% Run with mRGCMosaic (linear response, only ON-center mosaic), dynamic stimulus
    t_isoresponseLMplaneEllipses(...
        'useMetaContrast', true, ...
        'whichNoiseFreeNre', 'mRGCMosaic', ...
        'mRGCOutputSignalType', 'mRGCs', ...      % Select between {'cones', 'mRGCs'}
        'whichNoisyInstanceNre', 'Gaussian', ...
        'gaussianSigma', 0.1, ...
        'whichClassifierEngine', 'rceTemplateDistance', ...
        'useConeContrast', true, ...
        'useFixationalEMs', false, ...
        'temporalFilterValues', [], ...
        'oiPadMethod', 'zero', ...
        'visualizeEachScene', false, ...
        'visualizeEachResponse', true, ...
        'responseVisualizationFunction', @nreVisualizeMRGCmosaic, ...
        'maxVisualizedNoisyResponseInstances', 2, ...
        'maxVisualizedNoisyResponseInstanceStimuli', 2, ...
        'opticsType', 'loadComputeReadyRGCMosaic', ...
        'mRGCMosaicRawEccDegs', [0 0], ...
        'mRGCMosaicRawSizeDegs', [2 2], ...
        'mosaicEccDegs', [0.0 0], ...
        'mosaicSizeDegs', [1 1], ...
        'stimSizeDegs', 1.2, ...
        'pixelsNum', 512, ...
        'spatialFrequency', 0.0, ...
        'presentationMode', 'counterphasemodulated', ...
        'numberOfFrames', 8, ...
        'stimOnFrameIndices', 1:8, ...
        'frameDurationSeconds', (1/2.5)/8, ...
        'temporalFrequencyHz', 2.5, ...
        'nTest', 1024, ...
        'psychometricCurveSamplesNum', 5, ...
        'examinedDirectionsOnLMplane', 0:45:(360-45), ...
        'validationThresholds', []);

% Run with mRGCMosaic (linear response)
     t_isoresponseLMplaneEllipses(...
        'useMetaContrast', true, ...
        'whichNoiseFreeNre', 'mRGCMosaic', ...
        'mRGCOutputSignalType', 'mRGCs', ... 
        'whichNoisyInstanceNre', 'Gaussian', ...
        'whichClassifierEngine', 'rceTemplateDistance', ...
        'useConeContrast', true);

% Run with mRGCMosaic (ON mRGC mosaic, non-linear activation function, simulating case in a high saturation regime)
     t_isoresponseLMplaneEllipses(...
        'useMetaContrast', true, ...
        'whichNoiseFreeNre', 'mRGCMosaic', ...
        'simulateHighSaturationRegime', true, ...
        'simulateONOFFmosaic', false, ...
        'mRGCOutputSignalType', 'mRGCs', ...  
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
    options.mosaicEccDegs (1,2) double = [0 0];
    options.mosaicSizeDegs (1,2) double = [0.5 0.5];

    % Variable stimulus parameter
    options.examinedDirectionsOnLMplane (1,:) double = 0:45:315;

    % Fixed stimulus parameters
    options.meanLuminanceCdPerM2 (1,1) double = 100;
    options.meanChromaticityXY (1,2) double = [0.30 0.32];
    options.spatialFrequency(1,1) double = 0.0;
    options.orientationDegs (1,1) double = 90;
    options.spatialPhaseDegs (1,1) double = 90
    options.numberOfFrames double = []
    options.frameDurationSeconds (1,1) double = 0.1;
    options.temporalFrequencyHz (1,1) double = 5;
    options.stimOnFrameIndices (1,:) double = [];
    
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

% Mosaic sizes
mosaicEccDegs = options.mosaicEccDegs;
mosaicSizeDegs = options.mosaicSizeDegs;
mRGCRawEccDegs = options.mRGCMosaicRawEccDegs;
mRGCRawSizeDegs = options.mRGCMosaicRawSizeDegs;
mRGCCropSize = mosaicSizeDegs;

meanLuminanceCdPerM2 = options.meanLuminanceCdPerM2;
meanChromaticityXY = options.meanChromaticityXY;
orientationDegs = options.orientationDegs;
spatialPhaseDegs = options.spatialPhaseDegs;
temporalFrequencyHz = options.temporalFrequencyHz;
stimOnFrameIndices = options.stimOnFrameIndices;


temporalFilterValues = options.temporalFilterValues;
numberOfFrames = options.numberOfFrames;
frameDurationSeconds = options.frameDurationSeconds;
presentationMode = options.presentationMode;
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
figureFileBaseDir = setupFigureDirectory(mfilename, ...
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

%% Spatial frequency & orientation
examinedSpatialFrequencyCPD = options.spatialFrequency;
orientationDegs = options.orientationDegs;

debugStimulusConfig = ~true;

% Max RMS contrast (so as to keep stimuli within the display gamut)
rmsLMconeContrast = 0.12;

if (debugStimulusConfig)
    rmsLMconeContrast = 1.0;
    thresholdPara.logThreshLimitLow = 0.1;
end

% Compute directions based on rmsLMconeContrast
[theLMSconeContrastDirections, examinedDirectionsOnLMplane] = ...
    computeLMSconeContrastDirections(rmsLMconeContrast, examinedDirectionsOnLMplane);







%% Set grating engine parameters

if (strcmp(presentationMode, 'static'))
    gratingSceneParams = struct( ...
            'meanLuminanceCdPerM2', meanLuminanceCdPerM2, ...
            'meanChromaticityXY', meanChromaticityXY, ...
            'spectralSupport', 400:20:750, ...
            'fovDegs', stimSizeDegs, ...
            'pixelsNum', pixelsNum, ...
            'spatialEnvelope', 'rect', ...
            'spatialEnvelopeRadiusDegs', stimSizeDegs, ...
            'orientation', orientationDegs, ...
            'presentationMode', 'flashed', ... 
            'duration', 50/1000, ...
            'spatialPhase', spatialPhaseDegs);
else
    % Non-static
    stimulusDuration = framesNum*frameDurationSeconds;
    spatialPhaseAdvanceDegs = 360*temporalFrequencyHz/(framesNum+1);

    gratingSceneParams = struct( ...
        'meanLuminanceCdPerM2', meanLuminanceCdPerM2, ...
        'meanChromaticityXY', meanChromaticityXY, ...
        'spectralSupport', 400:20:750, ...
        'fovDegs', stimSizeDegs, ...
        'pixelsNum', pixelsNum, ...
        'spatialEnvelope', 'rect', ...
        'spatialEnvelopeRadiusDegs', stimSizeDegs, ...
        'orientation', orientationDegs, ...
        'presentationMode', presentationMode, ...
        'duration', stimulusDuration, ...
        'frameDurationSeconds', frameDurationSeconds, ...
        'spatialPhase', spatialPhaseDegs, ...
        'spatialPhaseAdvanceDegs', spatialPhaseAdvanceDegs);
end


% Plot all stimuli
skippedDirections = 0;
if (length(examinedDirectionsOnLMplane) > 50)
    % Too many, plot every other stimulus
    skippedDirections = 2;
end

























%% Configure our stimulus scene engines 
[theNullSceneEngine, theTestSceneEngines] = configureStimulusSceneEngines(...
    theLMSconeContrastDirections, examinedSpatialFrequencyCPD, gratingSceneParams);

%% Thresholds filename
if (isempty(thresholdsDataFileName))
    thresholdsDataFileName = ...
        sprintf('%sCSF_SF_%2.2fCPD_Optics_%s_EccDegs_x%2.1f_%2.1f_SizeDegs_%2.1fx%2.1f_OriDegs_%2.0f_%s.mat', ...
        mRGCOutputSignalType, ...
        examinedSpatialFrequencyCPD, ...
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

        % 2. We can crop the mRGCmosaic to some desired size.
        %    Passing [] for sizeDegs will not crop.
        %    Passing [] for eccentricityDegs will crop the mosaic at its center.
        nreNoiseFreeParams.mRGCMosaicParams.cropParams = struct(...
            'sizeDegs', mRGCCropSize, ...
            'eccentricityDegs', mosaicEccDegs ...
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

        if (useConeContrast) || (ischar(temporalFilterValues))
            % Since we specified that mRGCs will operate on cone-contrast responses,
            % will need a null scene for normalization
            nreNoiseFreeParams.nullStimulusSceneSequence = theNullSceneEngine.compute(0.0);
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
        theVisualizedConeContrastOffset = [0.0 0.0];

        if (simulateHighSaturationRegime)
            % Push responses into high saturation regime
            nreNoiseFreeParams.mRGCMosaicParams.responseBias = 0.03;
        
            % This simply translates the visualized ellipse on the LM plane 
            % (simulating a non-adapted simulation)
            theVisualizedConeContrastOffset = [0.1 0.1];
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

        % NON-LINEARITY HACKS
        % This simply translates the visualized ellipse on the LM plane 
        % (simulating a non-adapted simulation)
        theVisualizedConeContrastOffset = [0.0 0.0];

    otherwise
        error('Unsupported noise free neural response engine: ''%s''.', whichNoiseFreeNre);
end  % switch (whichNoiseFreeNre)



















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
if (useConeContrast) 
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
save(thresholdsDataFileName, 'options', 'examinedDirectionsOnLMplane', 'thresholdContrasts');

% Figure for plotting the thresholds on the LM plane
% along with the stimuli
hFigStimuliAndThresholds = figure(2346); clf;
set(hFigStimuliAndThresholds, 'Color', [1 1 1], 'Position', [10 10 1200 1200]);
set(hFigStimuliAndThresholds, 'HandleVisibility', 'off');

exportFig = true;
initializeFig = true;
theThresholdAxes = []; 
maxVisualizedThreshold = 0.2; % 0.05;
visualizeStimuliAndThresholdsOnLMPlane(...
            rmsLMconeContrast, ...
            examinedSpatialFrequencyCPD, gratingSceneParams, ...
            examinedDirectionsOnLMplane, skippedDirections, ...
            thresholdContrasts, maxVisualizedThreshold, theVisualizedConeContrastOffset, ...
            figureFileBaseDir, hFigStimuliAndThresholds, ...
            exportFig, initializeFig, theThresholdAxes);



%% Restore handle visibility so that a future close will close the figures
hAllFigs = findall(groot,'Type','figure');
for i = 1:numel(hAllFigs)
    set(hAllFigs(i), 'HandleVisibility', 'on')
end

end


%
% SUPPORTING FUNCTIONS
%


function [LMSconeContrasts, examinedDirectionsOnLMplane] = computeLMSconeContrastDirections(...
    rmsLMconeContrast, examinedDirectionsOnLMplane)

    if (isempty(examinedDirectionsOnLMplane))
        examinedDirectionsOnLMplane = [0:10:20 25:5:65 70:10:170];
        examinedDirectionsOnLMplane = cat(2, examinedDirectionsOnLMplane, 180+examinedDirectionsOnLMplane);
    end

    LMSconeContrasts = zeros(3, numel(examinedDirectionsOnLMplane));
    if (numel(rmsLMconeContrast) == numel(examinedDirectionsOnLMplane))
        LMSconeContrasts(1,:) = rmsLMconeContrast .* cosd(examinedDirectionsOnLMplane);
        LMSconeContrasts(2,:) = rmsLMconeContrast .* sind(examinedDirectionsOnLMplane);
    else
        LMSconeContrasts(1,:) = rmsLMconeContrast(1) * cosd(examinedDirectionsOnLMplane);
        LMSconeContrasts(2,:) = rmsLMconeContrast(1) * sind(examinedDirectionsOnLMplane);
    end
end

function theThresholdAxes = visualizeStimuliAndThresholdsOnLMPlane(...
    rmsLMconeContrast, examinedSpatialFrequencyCPD, gratingSceneParams, ...
    theChromaticDirections, skippedDirections, theThresholds, ...
    maxThreshold, theVisualizedConeContrastOffset, ...
    figureFileBaseDir, hFig, ...
    exportFig, initializeFig, theThresholdAxes)

    
    %gratingSceneParams.warningInsteadOfErrorOnOutOfGamut = true;

    visualizedLMangles = 0:15:345;
    theLMSconeContrastDirections = computeLMSconeContrastDirections(rmsLMconeContrast, visualizedLMangles);

    % Create the background scene engine
    theNullSceneEngine = createGratingSceneEngine(...
            [0 0 0.4], examinedSpatialFrequencyCPD, ...
            gratingSceneParams);

    % Compute the scene sequence
    theNullSceneSequence = theNullSceneEngine.compute(0.0);

    % Draw
    set(hFig, 'HandleVisibility', 'on');
    figure(hFig);
    
    plotColorCircle = false;

    if (initializeFig)
        clf;
        
        if (plotColorCircle)
            % Fill the figure with the background scene
            theBackgroundAxes = axes('Position', [0 0 1 1]);
            % Fill the figure with the background scene
            theNullSceneEngine.visualizeStaticFrame(...
                    theNullSceneSequence, ...
                    'frameToVisualize', 1, ...
                    'axesHandle', theBackgroundAxes, ...
                    'sRGBforSceneVisualization', true);
            axis(theBackgroundAxes, 'equal');
    
        
            % Superimpose the examined stimuli
            panelWidth = 0.085;
            panelR = 0.5-0.5*panelWidth;
        
            for iDir = 1:(skippedDirections+1):size(theLMSconeContrastDirections,2)
        
                % Create a grating scene engine with the examined chromatic direction
                theSceneEngine = createGratingSceneEngine(...
                    theLMSconeContrastDirections(:,iDir), examinedSpatialFrequencyCPD, ...
                    gratingSceneParams);
        
                % Compute the stimulus scene at max contrast
                testContrast = 1;
                theSceneSequence = theSceneEngine.compute(testContrast);
        
                % Visualize the stimulus scene
                figure(hFig);
                rmsC = norm(theLMSconeContrastDirections(:,iDir));
                theAxes = axes('Position', [panelR+panelR*0.95*theLMSconeContrastDirections(1,iDir)/rmsC panelR + panelR*0.95*theLMSconeContrastDirections(2,iDir)/rmsC panelWidth panelWidth]);
                theSceneEngine.visualizeStaticFrame(...
                    theSceneSequence, ...
                    'frameToVisualize', 1, ...
                    'axesHandle', theAxes, ...
                    'sRGBforSceneVisualization', true);
                axis(theAxes, 'equal');
                set(theAxes, 'Xcolor', 'none', 'YColor', 'none');
                
            end % iDir
        end

        if (isempty(theThresholdAxes))
            thresholdFigureHalfSize = 0.42;
            theThresholdAxes = axes(...
                'Position', [0.5-thresholdFigureHalfSize 0.5-thresholdFigureHalfSize thresholdFigureHalfSize*2 thresholdFigureHalfSize*2]);
            set(theThresholdAxes, 'XLim', maxThreshold*[-1 1], 'YLim', maxThreshold*[-1 1]);
            axis(theThresholdAxes, 'square');
            box(theThresholdAxes, 'off');
            set(theThresholdAxes, 'Color', 'none', 'Box', 'off', 'XColor', 'none', 'YColor', 'none');
        end
        
    end % if (initializeFig)


    thresholdConeContrasts(1,:) = cosd(theChromaticDirections) .* theThresholds;
    thresholdConeContrasts(2,:) = sind(theChromaticDirections) .* theThresholds;

    
    % Plot the axes
    plot(theThresholdAxes,  maxThreshold*[-1 1], [0 0], 'k-', 'LineWidth', 1.0);
    hold(theThresholdAxes, 'on');
    plot(theThresholdAxes,  [0 0], maxThreshold*[-1 1], 'k-', 'LineWidth', 1.0);

    % Plot the data points
    scatter(theThresholdAxes, theVisualizedConeContrastOffset(1) + thresholdConeContrasts(1,:), ...
        theVisualizedConeContrastOffset(2) + thresholdConeContrasts(2,:), 111, ...
        'MarkerEdgeColor', [0.99 0.0 0.0], 'MarkerFaceColor', [0.95 0.5 0.5], ...
        'LineWidth', 2.0, 'MarkerFaceAlpha', 0.7);

    if (numel(theChromaticDirections)>6)
        [z, a, b, rotationRadians] = fitEllipseToXYpoints(...
            thresholdConeContrasts, ...
            'nonLinear', false);
        %[z, a, b, rotationRadians] = fitellipse(thresholdConeContrasts, 'linear');
    
        % Plot the fitted ellipse
        % form the parameter vector
        npts = 100;
        t = linspace(0, 2*pi, npts);
    
        % Rotation matrix
        Q = [cos(rotationRadians), -sin(rotationRadians); sin(rotationRadians) cos(rotationRadians)];
        % Ellipse points
        X = Q * [a * cos(t); b * sin(t)] + repmat(z, 1, npts);

        % PLot the ellipse points
        h = plot(theThresholdAxes, theVisualizedConeContrastOffset(1) + X(1,:), theVisualizedConeContrastOffset(2) + X(2,:), 'k-', 'LineWidth', 6.0);
        h = plot(theThresholdAxes, theVisualizedConeContrastOffset(1) + X(1,:), theVisualizedConeContrastOffset(2) + X(2,:), 'c-', 'LineWidth',2.0);
    end

    
    
    hold(theThresholdAxes, 'off');
    set(theThresholdAxes, 'XLim', maxThreshold*[-1 1], 'YLim', maxThreshold*[-1 1]);
    axis(theThresholdAxes, 'square');
    box(theThresholdAxes, 'off');
    set(theThresholdAxes, 'Color', 'none', 'Box', 'off', 'XColor', [0.1 0.1 0.1], 'YColor', [0.1 0.1 0.1]);
    xlabel(theThresholdAxes, 'threshold contrast (L-cone)');
    ylabel(theThresholdAxes, 'threshold contrast (M-cone)');
    set(theThresholdAxes, 'FontSize', 24);
    drawnow;

    if ((~isempty(figureFileBaseDir)) && (exportFig))
        NicePlot.exportFigToPNG(fullfile(figureFileBaseDir,'/stimuliOnLMplane.png'), hFig, 300);
    end

    %set(hFig, 'HandleVisibility', 'off');
end


function [theNullSceneEngine, theTestSceneEngines] = configureStimulusSceneEngines(...
    theLMSconeContrastDirections, examinedSpatialFrequencyCPD, gratingSceneParams)

    % Create the background scene engine
    theNullSceneEngine = createGratingSceneEngine(...
            [0 0 0.4], examinedSpatialFrequencyCPD, ...
            gratingSceneParams);

    % Store test scene engines
    theTestSceneEngines = cell(1,size(theLMSconeContrastDirections,2));

    for iDir = 1:size(theLMSconeContrastDirections,2)
        % Create a grating scene engine with the examined chromatic direction
        theTestSceneEngines{iDir} = createGratingSceneEngine(...
            theLMSconeContrastDirections(:,iDir), examinedSpatialFrequencyCPD, ...
            gratingSceneParams);
    end % iDir

end




function figureFileBaseDir = setupFigureDirectory(theScriptName, ...
    useMetaContrast,useConeContrast,useFixationalEMs,...
    whichNoiseFreeNre,whichNoisyInstanceNre,...
    whichClassifierEngine,mRGCOutputSignalType)

    % Make sure local/figures directory exists so we can write out our figures in peace
    projectBaseDir = ISETBioCSFGeneratorRootPath;
    if (~exist(fullfile(projectBaseDir,'local', theScriptName, 'figures'),'dir'))
        mkdir(fullfile(projectBaseDir,'local', theScriptName, 'figures'));
        fprintf('Generated figure directory at %s\n', fullfile(projectBaseDir,'local', theScriptName, 'figures'))
    end

    figureFileBaseDir = fullfile(projectBaseDir,'local',mfilename,'figures', ...
        sprintf('%s_Meta_%d_ConeContrast_%d_FEMs_%d_%s_%s_%s_%s', mfilename, ...
        useMetaContrast,useConeContrast,useFixationalEMs,whichNoiseFreeNre,whichNoisyInstanceNre,...
        whichClassifierEngine,mRGCOutputSignalType));
    if (~exist(figureFileBaseDir, 'dir'))
        mkdir(figureFileBaseDir);
        fprintf('Generated figure sub-directory at %s\n', figureFileBaseDir);
    end


end


function [theThresholds, theVisualizedConeContrastOffset] = previouslyComputedThresholds()

theConeMosaicThresholds = [ ...
        0.0092 ...
    0.0089 ...
    0.0088 ...
    0.0101 ...
    0.0085 ...
    0.0096 ...
    0.0099 ...
    0.0107 ...
    0.0109 ...
    0.0128 ...
    0.0116 ...
    0.0126 ...
    0.0133 ...
    0.0132 ...
    0.0160 ...
    0.0108 ...
    0.0137 ...
    0.0102 ...
    0.0115 ...
    0.0083 ...
    0.0089 ...
    0.0085 ...
    0.0083 ...
    0.0085 ...
    0.0090 ...
    0.0080 ...
    0.0080 ...
    0.0092 ...
    0.0100 ...
    0.0098 ...
    0.0105 ...
    0.0123 ...
    0.0116 ...
    0.0126 ...
    0.0122 ...
    0.0125 ...
    0.0148 ...
    0.0135 ...
    0.0137 ...
    0.0128 ...
    0.0121 ...
    0.0109 ...
    0.0094 ...
    0.0101 ...
    0.0099 ...
    0.0091 ...
        ];

    theThresholdsLowSaturation = [ ...
         0.0077  ...
    0.0086 ...
    0.0109 ...
    0.0142 ...
    0.0152 ...
    0.0185 ...
    0.0175 ...
    0.0206 ...
    0.0201 ...
    0.0185 ...
    0.0177 ...
    0.0130 ...
    0.0122 ...
    0.0093 ...
    0.0078 ...
    0.0064 ...
    0.0059 ...
    0.0057 ...
    0.0048 ...
    0.0056 ...
    0.0060 ...
    0.0065 ...
    0.0068 ...
    0.0071 ...
    0.0093 ...
    0.0115 ...
    0.0155 ...
    0.0149 ...
    0.0174 ...
    0.0187 ...
    0.0234 ...
    0.0212 ...
    0.0198 ...
    0.0161 ...
    0.0144 ...
    0.0114 ...
    0.0083 ...
    0.0078 ...
    0.0061 ...
    0.0063 ...
    0.0061 ...
    0.0052 ...
    0.0048 ...
    0.0057 ...
    0.0059 ...
    0.0074 ...
        ];


    theThresholdsHeavySaturation = [...
        0.0208 ...
    0.0250 ...
    0.0335 ...
    0.0374 ...
    0.0430 ...
    0.0487 ...
    0.0540 ...
    0.0593 ...
    0.0543 ...
    0.0515 ...
    0.0426 ...
    0.0345 ...
    0.0341 ...
    0.0283 ...
    0.0215 ...
    0.0177 ...
    0.0171 ...
    0.0164 ...
    0.0166 ...
    0.0175 ...
    0.0138 ...
    0.0172 ...
    0.0195 ...
    0.0202 ...
    0.0257 ...
    0.0291 ...
    0.0348 ...
    0.0436 ...
    0.0475 ...
    0.0532 ...
    0.0597 ...
    0.0500 ...
    0.0430 ...
    0.0403 ...
    0.0364 ...
    0.0311 ...
    0.0256 ...
    0.0185 ...
    0.0174 ...
    0.0175 ...
    0.0138 ...
    0.0151 ...
    0.0143 ...
    0.0169 ...
    0.0173 ...
    0.0186 ...
        ];

    theThresholds = theThresholdsLowSaturation;
    theVisualizedConeContrastOffset = [0.0 0.0];

    %theThresholds = theThresholdsHeavySaturation;
    %theVisualizedConeContrastOffset = [0.15 0.15];

    %theThresholds = theConeMosaicThresholds
    %theVisualizedConeContrastOffset = [0.0 0.0];

end


function [z, a, b, alpha] = fitEllipseToXYpoints(xyPoints, varargin)

	% Parse input
    p = inputParser;
    p.addParameter('maxIterations', 200, @isnumeric);
    p.addParameter('tolerance', 1e-5, @isnumeric);
    p.addParameter('constraint', 'bookstein', @(x)(ismember(x, {'bookstein', 'trace'})));
    p.addParameter('nonLinear', true, @islogical);
    p.parse(varargin{:});

	rows = size(xyPoints,1);
	cols = size(xyPoints,2);
	assert(rows == 2, 'xyPoints must be a 2 x N matrix');
	assert(cols >=6, 'xyPoints must have at least 6 points');

	fitParams = struct(...
		'nonLinear', p.Results.nonLinear, ...
		'constraint', p.Results.constraint, ...
		'maxIterations', p.Results.maxIterations, ...
		'tolerance', p.Results.tolerance);

	% Remove centroid
	centroid = mean(xyPoints, 2);
	xyPoints = bsxfun(@minus, xyPoints, centroid);

	% Obtain a linear estimate
	switch fitParams.constraint
    	case 'bookstein'
        	[z, a, b, alpha] = fitbookstein(xyPoints);
    	case 'trace'
       		[z, a, b, alpha] = fitggk(xyPoints);
	end % switch


	if (fitParams.nonLinear)
		% Initial conditions
	    z0     = z;
	    a0     = a;
	    b0     = b;
	    alpha0 = alpha;

	    % Fit
    	[z, a, b, alpha, converged, isCircle] = fitNonLinear(xyPoints, z0, a0, b0, alpha0, params);

    	% Return linear estimate if GN doesn't converge or if the data points fall on a circle
	    if (~converged) || (isCircle)
	        fprintf('*** FailureToConverge: Gauss-Newton did not converge, returning linear estimate.\n');
	        z = z0;
	        a = a0;
	        b = b0;
	        alpha = alpha0;
	    end
	end % if (fitParams.nonLinear)

	% Add the centroid back on
	z = z + centroid;
end


function [z, a, b, alpha] = fitbookstein(x)
	%FITBOOKSTEIN   Linear ellipse fit using bookstein constraint
	%   lambda_1^2 + lambda_2^2 = 1, where lambda_i are the eigenvalues of A

	% Convenience variables
	m  = size(x, 2);
	x1 = x(1, :)';
	x2 = x(2, :)';

	% Define the coefficient matrix B, such that we solve the system
	% B *[v; w] = 0, with the constraint norm(w) == 1
	B = [x1, x2, ones(m, 1), x1.^2, sqrt(2) * x1 .* x2, x2.^2];

	% To enforce the constraint, we need to take the QR decomposition
	[Q, R] = qr(B);

	% Decompose R into blocks
	R11 = R(1:3, 1:3);
	R12 = R(1:3, 4:6);
	R22 = R(4:6, 4:6);

	% Solve R22 * w = 0 subject to norm(w) == 1
	[U, S, V] = svd(R22);
	w = V(:, 3);

	% Solve for the remaining variables
	v = -R11 \ R12 * w;

	% Fill in the quadratic form
	A        = zeros(2);
	A(1)     = w(1);
	A([2 3]) = 1 / sqrt(2) * w(2);
	A(4)     = w(3);
	bv       = v(1:2);
	c        = v(3);

	% Find the parameters
	[z, a, b, alpha] = conic2parametric(A, bv, c);
end

function [z, a, b, alpha] = fitggk(x)
	% Linear least squares with the Euclidean-invariant constraint Trace(A) = 1
	% Convenience variables
	m  = size(x, 2);
	x1 = x(1, :)';
	x2 = x(2, :)';

	% Coefficient matrix
	B = [2 * x1 .* x2, x2.^2 - x1.^2, x1, x2, ones(m, 1)];
	v = B \ -x1.^2;

	% For clarity, fill in the quadratic form variables
	A        = zeros(2);
	A(1,1)   = 1 - v(2);
	A([2 3]) = v(1);
	A(2,2)   = v(2);
	bv       = v(3:4);
	c        = v(5);

	% find parameters
	[z, a, b, alpha] = conic2parametric(A, bv, c);
end


function [z, a, b, alpha, converged, isCircle] = fitNonLinear(x, z0, a0, b0, alpha0, params)
	% Gauss-Newton least squares ellipse fit minimising geometric distance 

	% Get initial rotation matrix
	Q0 = [cos(alpha0), -sin(alpha0); sin(alpha0) cos(alpha0)];
	m = size(x, 2);

	% Get initial phase estimates
	phi0 = angle( [1 i] * Q0' * (x - repmat(z0, 1, m)) )';
	u = [phi0; alpha0; a0; b0; z0];

	% Iterate using Gauss Newton
	converged = false;

	for nIts = 1:params.maxIterations
	    % Find the function and Jacobian
	    [f, J, isCircle] = computeJacobian(u);
    
    	if (isCircle)
    		fprintf('Ellipse is near-circular - nonlinear fit may not succeed\n.')
    	end

	    % Solve for the step and update u
	    h = -J \ f;
	    u = u + h;
    
	    % Check for convergence
	    delta = norm(h, inf) / norm(u, inf);
	    if delta < params.tolerance
	        converged = true;
	        break
	    end
	end

	alpha = u(end-4);
	a = u(end-3);
	b = u(end-2);
	z = u(end-1:end);

	% ---- Nested function ---
	function [f, J, isCircle] = computeJacobian(u)
        % Define the system of nonlinear equations and Jacobian. 

        % Tolerance for whether it is a circle
        circTol = 1e-5;
        
        % Unpack parameters from u
        phi   = u(1:end-5);
        alpha = u(end-4);
        a     = u(end-3);
        b     = u(end-2);
        z     = u(end-1:end);
        
        % If it is a circle, the Jacobian will be singular, and the
        % Gauss-Newton step won't work. 
        %TODO: This can be fixed by switching to a Levenberg-Marquardt
        %solver
        if (abs(a - b) / (a + b) < circTol)
            isCircle = true;
        else
        	isCircle = false;
        end

        % Convenience trig variables
        c = cos(phi);
        s = sin(phi);
        ca = cos(alpha);
        sa = sin(alpha);
        
        % Rotation matrices
        Q    = [ca, -sa; sa, ca];
        Qdot = [-sa, -ca; ca, -sa];

        % Preallocate function and Jacobian variables
        f = zeros(2 * m, 1);
        J = zeros(2 * m, m + 5);
        for i = 1:m
            rows = (2*i-1):(2*i);
            % Equation system - vector difference between point on ellipse
            % and data point
            f((2*i-1):(2*i)) = x(:, i) - z - Q * [a * cos(phi(i)); b * sin(phi(i))];
            
            % Jacobian
            J(rows, i) = -Q * [-a * s(i); b * c(i)];
            J(rows, (end-4:end)) = ...
                [-Qdot*[a*c(i); b*s(i)], -Q*[c(i); 0], -Q*[0; s(i)], [-1 0; 0 -1]];
        end
    end % ---- Nested function ---
end

function [z, a, b, alpha] = conic2parametric(A, bv, c)
	% Diagonalise A - find Q, D such at A = Q' * D * Q
	[Q, D] = eig(A);
	Q = Q';

	% If the determinant < 0, it's not an ellipse
	if prod(diag(D)) <= 0 
	    error('NotEllipse', 'Linear fit did not produce an ellipse');
	end

	% We have b_h' = 2 * t' * A + b'
	t = -0.5 * (A \ bv);

	c_h = t' * A * t + bv' * t + c;

	z = t;
	a = sqrt(-c_h / D(1,1));
	b = sqrt(-c_h / D(2,2));
	alpha = atan2(Q(1,2), Q(1,1));
end % conic2parametric
