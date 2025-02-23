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
     
    % Validate with cone contrast, temporal filter, meta contrast.  This
    % is a delta function and should not change threshold relative to the
    % same parameters without a filter.
    t_spatialCSF('useMetaContrast', true, ...
        'whichNoiseFreeNre', 'excitationsCmosaic', ...
        'whichNoisyInstanceNre', 'Poisson', ...
        'whichClassifierEngine', 'rcePoisson', ...
        'useConeContrast', false, ...
        'useFixationalEMs', false, ...
        'numberOfFrames', 4, ...
        'temporalFilterValues', [1 0 0 0], ...
        'oiPadMethod', 'mean', ...
        'validationThresholds',  [0.0208    0.0370    0.0770    0.3300]);

    % Validate with cone contrast, temporal filter based on photocurrent, meta contrast.
    t_spatialCSF('useMetaContrast', true, ...
        'whichNoiseFreeNre', 'excitationsCmosaic', ...
        'whichNoisyInstanceNre', 'Poisson', ...
        'whichClassifierEngine', 'rcePoisson', ...
        'useConeContrast', false, ...
        'useFixationalEMs', false, ...
        'numberOfFrames', 4, ...
        'temporalFilterValues', 'photocurrentImpulseResponseBased', ...
        'oiPadMethod', 'mean');

    % Verify that cMosaic FEMs work without crashing. Same EM paths on
    % train and test.
    t_spatialCSF('useMetaContrast', false, ...
        'whichNoiseFreeNre', 'excitationsCmosaic', ...
        'whichNoisyInstanceNre', 'Poisson', ...
        'whichClassifierEngine', 'rcePoisson', ...
        'useConeContrast', false, ...
        'useFixationalEMs', true, ...
        'testEMsMatchTrain', true, ...
        'nTrainEMs', 4, ...
        'nTestEMs', 4, ...
        'nTrain', 1, ...
        'nTest', 64', ...
        'temporalFilterValues', [], ...
        'oiPadMethod', 'zero', ...
        'validationThresholds', [0.0258    0.0526    0.1576    0.4792], ...
        'visualizeEachScene', false, ...
        'visualizeEachResponse', false, ...
        'responseVisualizationFunction', @nreVisualizeCMosaic, ...
        'maxVisualizedNoisyResponseInstances', 2, ...
        'maxVisualizedNoisyResponseInstanceStimuli',2);

    % Verify that cMosaic FEMs work without crashing. Different EM paths on
    % train and test. 
    %
    % Can't use meta contrast for this case.
    %
    % DHB: THRESHOLD VALUES DO NOT MAKE SENSE TO ME FOR THIS CASE.  EDGE
    % EFFECTS?
    t_spatialCSF('useMetaContrast', false, ...
        'whichNoiseFreeNre', 'excitationsCmosaic', ...
        'whichNoisyInstanceNre', 'Poisson', ...
        'whichClassifierEngine', 'rcePoisson', ...
        'useConeContrast', false, ...
        'useFixationalEMs', true, ...
        'testEMsMatchTrain', false, ...
        'nTrainEMs', 4, ...
        'nTestEMs', 4, ...
        'nTrain', 1, ...
        'nTest', 64', ...
        'temporalFilterValues', [], ...
        'oiPadMethod', 'zero', ...
        'validationThresholds', [1.0000    1.0000    0.1939    0.6992], ...
        'visualizeEachScene', false, ...
        'visualizeEachResponse', false, ...
        'responseVisualizationFunction', @nreVisualizeCMosaic, ...
        'maxVisualizedNoisyResponseInstances', 2, ...
        'maxVisualizedNoisyResponseInstanceStimuli',2);

    % Verify that Gaussian noise works, as well as template classifier
    t_spatialCSF('useMetaContrast', true, ...
        'whichNoiseFreeNre', 'excitationsCmosaic', ...
        'whichNoisyInstanceNre', 'Gaussian', ...
        'whichClassifierEngine', 'rceTemplateDistance', ...
        'useConeContrast', false, ...
        'useFixationalEMs', false, ...
        'temporalFilterValues', [], ...
        'oiPadMethod', 'mean', ...
        'validationThresholds', [0.0993    0.1673    0.3377    1.0000]);

    % Verify that rcePcaSVM works
    t_spatialCSF('useMetaContrast', true, ...
        'whichNoiseFreeNre', 'excitationsCmosaic', ...
        'whichNoisyInstanceNre', 'Poisson', ...
        'whichClassifierEngine', 'rcePcaSVM', ...
        'useConeContrast', false, ...
        'useFixationalEMs', false, ...
        'temporalFilterValues', [], ...
        'oiPadMethod', 'mean', ...
        'validationThresholds', [0.0377    0.0796    0.1427    0.6956]);

    % mRGCMosaic nre basic test
    t_spatialCSF('useMetaContrast', true, ...
        'whichNoiseFreeNre', 'mRGCMosaic', ...
        'mRGCOutputSignalType', 'mRGCs', ...
        'whichNoisyInstanceNre', 'Gaussian', ...
        'whichClassifierEngine', 'rceTemplateDistance', ...
        'useConeContrast', true, ...
        'useFixationalEMs', false, ...
        'temporalFilterValues', [], ...
        'oiPadMethod', 'zero', ...
        'validationThresholds', [0.0507    0.0347    0.0407    0.0512], ...
        'visualizeEachScene', false, ...
        'visualizeEachResponse', false, ...
        'responseVisualizationFunction', @nreVisualizeMRGCmosaic, ...
        'maxVisualizedNoisyResponseInstances', 2, ...
        'maxVisualizedNoisyResponseInstanceStimuli', 2);

    % mRGCMosaic nre basic test with fixational eye movements
    t_spatialCSF('useMetaContrast', true, ...
        'whichNoiseFreeNre', 'mRGCMosaic', ...
        'mRGCOutputSignalType', 'mRGCs', ...
        'whichNoisyInstanceNre', 'Gaussian', ...
        'whichClassifierEngine', 'rceTemplateDistance', ...
        'useConeContrast', true, ...
        'useFixationalEMs', true, ...
        'temporalFilterValues', [], ...
        'oiPadMethod', 'zero', ...
        'validationThresholds', [], ...
        'visualizeEachScene', false, ...
        'visualizeEachResponse', false, ...
        'responseVisualizationFunction', @nreVisualizeMRGCmosaic, ...
        'maxVisualizedNoisyResponseInstances', 2); 

    % mRGCMosaic nre, but pick off cone excitations as output
    t_spatialCSF('useMetaContrast', true, ...
        'whichNoiseFreeNre', 'mRGCMosaic', ...
        'mRGCOutputSignalType', 'cones', ...
        'whichNoisyInstanceNre', 'Gaussian', ...
        'whichClassifierEngine', 'rceTemplateDistance', ...
        'useConeContrast', false, ...
        'useFixationalEMs', false, ...
        'temporalFilterValues', [], ...
        'oiPadMethod', 'zero', ...
        'visualizeEachResponse', false, ...
        'maxVisualizedNoisyResponseInstances', 9, ... 
        'validationThresholds', [0.1522    0.1651    0.2929   0.6418]);

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

    % Choose noise model
    %   Choices: 'Poisson'
    %            'Gaussian'
    options.whichNoisyInstanceNre (1,:) char = 'Poisson'

    % Choose classifier engine
    %    rcePoisson - signal known exactly Poission max likelihood
    %    rceTemplateDistance - signal known exactly nearest L2 template
    %                 distance.
    %    rcePcaSVM  - support vector machine linear classifier after PCA.
    options.whichClassifierEngine (1,:) char = 'rcePoisson'

    % Use cone contrast rather than cone excitations
    options.useConeContrast (1,1) logical = false

    % Use fixational eye movements?
    %
    % The first variable here controls whether we simulate fixational EMs.
    % The second determines whether the training and test EM paths match
    % each other. It can only be true if the number of training and test
    % EMs match each other.
    options.useFixationalEMs (1,1) logical = false
    options.testEMsMatchTrain (1,1) logical = true

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
    % this number must be a multiple of the number of fixationsal EMs
    % specified.
    options.nTrain (1,1) double = 1
    options.nTest (1,1) double = 128

    % Number of temporal frames.  If not empty, overrides default values
    options.numberOfFrames double = []

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
    options.temporalFilterValues (1,:) = []

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
        'slopeDelta', 2.5/20, ...
        'thresholdCriterion', 0.81606)

    % Verbose?
    options.verbose (1,1) logical = true

    % Visualize the output of each compute (useful for debuging code)
    % Scene and response visualization are controlled separately.
    options.visualizeEachScene (1,1) logical = false
    options.visualizeEachResponse (1,1) logical = false
    options.responseVisualizationFunction = []
    options.maxVisualizedNoisyResponseInstances = 1
    options.maxVisualizedNoisyResponseInstanceStimuli = 1
end

%% Close any stray figs
close all;

%% Make sure local/figures directory exists so we can write out our figures in peace
projectBaseDir = ISETBioCSFGeneratorRootPath;
if (~exist(fullfile(projectBaseDir,'local',mfilename,'figures'),'dir'))
    mkdir(fullfile(projectBaseDir,'local',mfilename,'figures'));
end

%% Set flags from key/value pairs
filter = options.filter;
useMetaContrast = options.useMetaContrast;
useConeContrast = options.useConeContrast;
useFixationalEMs = options.useFixationalEMs;
testEMsMatchTrain = options.testEMsMatchTrain;
nTrainEMs = options.nTrainEMs;
nTestEMs = options.nTestEMs;
nTrain = options.nTrain;
nTest = options.nTest;
whichNoiseFreeNre = options.whichNoiseFreeNre;
whichNoisyInstanceNre = options.whichNoisyInstanceNre;
whichClassifierEngine = options.whichClassifierEngine;
validationThresholds = options.validationThresholds;
mRGCOutputSignalType = options.mRGCOutputSignalType;
temporalFilterValues = options.temporalFilterValues;
numberOfFrames = options.numberOfFrames;
fastParameters = options.fastParameters;
oiPadMethod = options.oiPadMethod;
thresholdPara = options.thresholdPara;
verbose = options.verbose;
visualizeEachScene = options.visualizeEachScene;
visualizeEachResponse = options.visualizeEachResponse;
responseVisualizationFunction = options.responseVisualizationFunction;
maxVisualizedNoisyResponseInstances = options.maxVisualizedNoisyResponseInstances;
maxVisualizedNoisyResponseInstanceStimuli = options.maxVisualizedNoisyResponseInstanceStimuli;

%% Freeze rng for replicatbility and validation
rng(1);

%% Figure output base name
figureFileBase = fullfile(projectBaseDir,'local',mfilename,'figures', ...
    sprintf('%s_Meta_%d_ConeContrast_%d_FEMs_%d_%s_%s_%s_%s', mfilename, ...
    useMetaContrast,useConeContrast,useFixationalEMs,whichNoiseFreeNre,whichNoisyInstanceNre,...
    whichClassifierEngine,mRGCOutputSignalType));

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
pixelsNum = 128;
gratingOrientationDegs = 90;
gratingSpatialPhase = 90;
frameDurationSeconds = 0.1;
theTemporalFrequencyHz = 5.0;

% Set up some sizes.  Note that these are small so that the examples run
% fast. 
mosaicEccDegs = [0 0];
mosaicSizeDegs = [0.5 0.5];
mRGCRawSizeDegs = [2 2];
mRGCCropSize = [0.5 0.5];
stimulusSizeDegs = 0.5;

%% Set up temporal filter if we have one. 
%
% Note that the nre returns the same number of frames it was passed.  So if
% you want post-stimulus responses produced by the temporal filter, you
% need to pad your input appropriately. That padding process is not
% illustrated in this tutorial script.
if (~isempty(temporalFilterValues))
    temporalFilter.filterValues = temporalFilterValues;
    if (ischar(temporalFilterValues))
        temporalFilter.temporalSupport = '';
    else
        temporalFilter.temporalSupport = frameDurationSeconds*(0:(length(temporalFilterValues)-1));
    end
else
    temporalFilter = [];
end

%% Null stimulus contrast
%
% Used when we operate with cone contrasts rather than excitations
nullContrast = 0.0;

%% List of spatial frequencies to be tested.
spatialFreqs = [4, 8, 16, 32];
if (length(spatialFreqs) ~= 4 | ~all(spatialFreqs == [4, 8, 16, 32]))
    validationThresholds = [];
end

%% Chromatic direction and contrast
%
% Choose stimulus chromatic direction specified as a 1-by-3 vector
% of L, M, S cone contrast.  These vectors get normalized below, so only
% their direction matters in the specification.
stimType = 'luminance';
switch (stimType)
    case 'luminance'
        chromaDir = [1.0, 1.0, 1.0]';
    case 'red-green'
        chromaDir = [1.0, -1.0, 0.0]';
        validationThresholds = [];
    case 'L-isolating'
        chromaDir = [1.0, 0.0, 0.0]';
        validationThresholds = [];
end

% Set the RMS cone contrast of the stimulus. Things may go badly if you
% exceed the gamut of the monitor, so we are conservative and set this at a
% value that is within gamut of typical monitors and don't worry about it
% further for this tutorial.  A vector length contrast of 0.1 should be OK.
rmsContrast = 0.1;
chromaDir = chromaDir / norm(chromaDir) * rmsContrast;
assert(abs(norm(chromaDir) - rmsContrast) <= 1e-10);

%% Create neural response engine
%
% Determine neural responses to analyze, and set it up.
switch (whichNoiseFreeNre)
    case 'mRGCMosaic'
        % mRGCMosaic responses
        nreNoiseFreeResponse = @nreNoiseFreeMidgetRGCMosaic;
        noiseFreeResponseParams = nreNoiseFreeMidgetRGCMosaic([],[],[],[], ...
            'oiPadMethod',oiPadMethod);

        % Modify certain parameters of interest
        %
        % 1. Select one of the pre-computed mRGC mosaics by specifying its
        % eccentricityDegs & sizeDegs asnd its center type
        noiseFreeResponseParams.mRGCMosaicParams.eccDegs = mosaicEccDegs;
        noiseFreeResponseParams.mRGCMosaicParams.sizeDegs = mRGCRawSizeDegs;
        noiseFreeResponseParams.mRGCMosaicParams.inputSignalType = 'coneContrast';
        noiseFreeResponseParams.mRGCMosaicParams.rgcType = 'ONcenterMidgetRGC';

        % 2. We can crop the mRGCmosaic to some desired size.
        %    Passing [] for sizeDegs will not crop.
        %    Passing [] for eccentricityDegs will crop the mosaic at its center.
        noiseFreeResponseParams.mRGCMosaicParams.cropParams = struct(...
            'sizeDegs', mRGCropSize, ...
            'eccentricityDegs', [] ...
            );

        % 3. Set the input cone mosaic integration time
        noiseFreeResponseParams.mRGCMosaicParams.coneIntegrationTimeSeconds = frameDurationSeconds;

        % 3. Set the noise level.
        %
        % Post-cone summation noise is additive Gaussian noise with a desired
        % sigma. When the input is raw cone excitations, the sigma should be
        % expressed in terms of cone excitations/integration time. When the input
        % to mRGCs is cone modulations with respect to the background, which have a
        % max amplitude of 1.0, the sigma should be scaled appropriately.
        % So if your mRGC mosiac were operating directly on cone excitations, a
        % different value more like 100 or 1000 would be reasonable, depending on
        % light level.
        noiseFreeResponseParams.mRGCMosaicParams.outputSignalType = mRGCOutputSignalType;
        switch (mRGCOutputSignalType)
            case 'mRGCs'
                gaussianSigma = 0.003;
            case 'cones'
                gaussianSigma = 50;
            otherwise
                error('Unknown mRGC output signal type specified');
        end

        % Handle cone contrast setting
        if (useConeContrast)
            noiseFreeResponseParams.mRGCMosaicParams.inputSignalType = 'coneContrast';
        else
            noiseFreeResponseParams.mRGCMosaicParams.inputSignalType = 'coneExcitations';
        end

        % Handle temporal filter
        noiseFreeResponseParams.temporalFilter = temporalFilter;

        % Sanity check on the amount of mRGCMosaicVMembraneGaussianNoiseSigma for
        % the specified noiseFreeParams.mRGCMosaicParams.inputSignalType
        switch (noiseFreeResponseParams.mRGCMosaicParams.inputSignalType)
            case 'coneContrast'
                % Cone contrast is in the range of [-1 1]
                if (gaussianSigma > 1)
                    error('Gaussian noise sigma (%f) seems too large when operating on ''%s''.', ...
                        gaussianSigma, ...
                        noiseFreeResponseParams.mRGCMosaicParams.inputSignalType);
                end

            case 'coneExcitations'
                % Excitations are unlikely to be of order 1.
                if (gaussianSigma < 1)
                    error('Gaussian noise sigma (%f) seems too small when operating on ''%s''.', ...
                        gaussianSigma, ...I
                        noiseFreeResponseParams.mRGCMosaicParams.inputSignalType);
                end

            otherwise
                error('Unknown mRGC signal type specified: ''%s''.', noiseFreeResponseParams.mRGCMosaicParams.inputSignalType)
        end

    case 'excitationsCmosaic'
        nreNoiseFreeResponse = @nreNoiseFreeCMosaic;
        noiseFreeResponseParams = nreNoiseFreeCMosaic([],[],[],[], ...
            'oiPadMethod',oiPadMethod);
        noiseFreeResponseParams.coneMosaicParams.sizeDegs = sizeDegs;
        noiseFreeResponseParams.coneMosaicParams.timeIntegrationSeconds = frameDurationSeconds;
        if (~all(noiseFreeResponseParams.coneMosaicParams.sizeDegs == sizeDegs))
            validationThresholds = [];
        elseif (noiseFreeResponseParams.coneMosaicParams.timeIntegrationSeconds ~= 0.1)
            validationThresholds = [];
        end

        % Handle cone contrast setting
        if (useConeContrast)
            noiseFreeResponseParams.coneMosaicParams.outputSignalType = 'coneContrast';
        else
            noiseFreeResponseParams.coneMosaicParams.outputSignalType = 'coneExcitations';
        end

        % Handle temporal filter
        noiseFreeResponseParams.temporalFilter = temporalFilter;

        % Set Gaussian sigma in case we're using Gaussian noise below.
        gaussianSigma = 50;

    case 'sceneAsResponses'
        % This is image photon counts, and thus provides an estimate of the
        % upper bound on performance for a photon limited system
        %
        % Override size and frame duration to put performance in reasonable range
        sizeDegs = [0.05 0.05];
        frameDurationSeconds = 1e-16;

        % Set up nre
        nreNoiseFreeResponse = @nreNoiseFreeSceneAsResponses;
        noiseFreeResponseParams = nreNoiseFreeSceneAsResponses;
        noiseFreeResponseParams.frameDurationSeconds = frameDurationSeconds;

        % This scene engine should not use cone contrast, no matter what
        useConeContrast = false;

    otherwise
        error('Unsupported noise free nre specified');
end

switch (whichNoisyInstanceNre)
    case 'Poisson'
        nreNoisyInstances = @nreNoisyInstancesPoisson;
        noisyInstancesParams = nreNoisyInstancesPoisson;

    case 'Gaussian'
        nreNoisyInstances = @nreNoisyInstancesGaussian;
        noisyInstancesParams = nreNoisyInstancesGaussian;
        noisyInstancesParams.sigma = gaussianSigma;

    otherwise
        error('Unsupported noisy instances nre specified');
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
        whichClassifierEngine = responseClassifierEngine(@rcePoisson);
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
questEnginePara = struct( ...
    'qpPF',@qpPFWeibullLog, ...
    'minTrial', 1280, ...
    'maxTrial', 1280, ...
    'numEstimator', 1, ...
    'stopCriterion', 0.05);

%% Set grating engine parameters
%
% DHB: I DON'T QUITE UNDERSTAND THE TWO CASES HERE.  WHY DON'T WE JUST
% SPECIFY STATIC OR DYNAMIC AS AN OPTION?
stimulusDuration = framesNum*frameDurationSeconds;
if (~useFixationalEMs & isempty(temporalFilter) & framesNum == 1)
    gratingSceneParams = struct( ...
        'fovDegs', stimSizeDegs, ...
        'presentationMode', 'flashedmultiframe', ...
        'duration', stimulusDuration, ...
        'frameDurationSeconds', stimulusDuration/framesNum, ...
        'orientation', gratingOrientationDegs, ...
        'spatialPhase', gratingSpatialPhase, ...
        'pixelsNum', pixelsNum, ...
        'filter', filter...
        );
else
    % Dynamic stimulus parameters
    presentationMode = 'counterphasemodulated';
    gratingSceneParams = struct( ...
        'fovDegs', stimSizeDegs, ...
        'presentationMode', presentationMode, ... 
        'duration', stimulusDuration, ...
        'frameDurationSeconds', stimulusDuration/framesNum, ...
        'temporalFrequencyHz', theTemporalFrequencyHz, ...
        'orientation', gratingOrientationDegs, ...
        'spatialPhase', gratingSpatialPhase, ...
        'pixelsNum', pixelsNum, ...
        'filter', filter ...
        );
    
    % Here are some other parameters that might get set for a dynamic
    % stimulus.  And we could consider 'drifted' as mode in addition, but
    % may need to special case 0 cpd as 'counterphasemodulated' if we do.
    %
    % 'spatialEnvelope', 'rect', ...
    % 'spatialEnvelopeRadiusDegs', theStimulusSpatialEnvelopeRadiusDegs, ...
    % 'minPixelsNumPerCycle', minPixelsNumPerCycle, ...
    % 'spatialPhaseAdvanceDegs', 360*(frameDurationSeconds*theTemporalFrequencyHz), ...
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
    nullStimulusSceneSequence = nullGratingSceneEngine.compute(nullContrast);
    noiseFreeResponseParams.nullStimulusSceneSequence = nullStimulusSceneSequence;
end

%% Create the neural engine
theNeuralEngine = neuralResponseEngine( ...
    nreNoiseFreeResponse, ...
    nreNoisyInstances, ...
    noiseFreeResponseParams, ...
    noisyInstancesParams ...
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

%% Set up EM object to define fixational eye movement paths
%
% There are lots of EM things we can control. 
%
% Setting the seed to -1 means don't touch the seed or call rng().
%
% DHB: CONSIDER WHETEHR EM PATHS SHOULD BE DIFFERENT FOR EACH SF.
if (useFixationalEMs)
    emMicrosaccadeType = 'none';
    emRandomSeed = -1;
    emComputeVelocitySignal = false;
    emCenterPaths = false;
    emCenterPathsAtSpecificTimeMsec = [];

    % Set up the train EM object
    trainFixationalEMObj = fixationalEM;
    trainFixationalEMObj.microSaccadeType = emMicrosaccadeType;
    trainFixationalEMObj.randomSeed = emRandomSeed ;
    trainFixationalEMObj.compute(stimulusDuration, frameDurationSeconds, nTrainEMs, emComputeVelocitySignal, ...
        'centerPaths', emCenterPaths, 'centerPathsAtSpecificTimeMsec', emCenterPathsAtSpecificTimeMsec);

    % And this one for testing. Here we are doing EM know exactly, so the
    % test object is the same as the training object.
    testFixationalEMObj = fixationalEM;
    testFixationalEMObj.microSaccadeType = emMicrosaccadeType;
    testFixationalEMObj.randomSeed = emRandomSeed;
    testFixationalEMObj.compute(stimulusDuration, frameDurationSeconds, nTestEMs, emComputeVelocitySignal, ...
        'centerPaths', emCenterPaths, 'centerPathsAtSpecificTimeMsec', emCenterPathsAtSpecificTimeMsec);
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

%% Compute threshold for each spatial frequency
% See toolbox/helpers for functions createGratingSceneEngine, computeThreshold,
% computePeformance
dataFig = figure();
for idx = 1:length(spatialFreqs)
    axLeft{idx}  = subplot(length(spatialFreqs), 2, idx * 2 - 1);
    axRight{idx} = subplot(length(spatialFreqs), 2, idx * 2);
end

logThreshold = zeros(1, length(spatialFreqs));
for idx = 1:length(spatialFreqs)
    % Create a static grating scene with a particular chromatic direction,
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
            'maxVisualizedNoisyResponseInstances', maxVisualizedNoisyResponseInstances);
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
            'maxVisualizedNoisyResponseInstances', maxVisualizedNoisyResponseInstances);
    end

    % Plot stimulus
    figure(dataFig);
    
    % This shows one frame of the scene.
    visualizeEachComputeSave = gratingSceneEngine.visualizeEachCompute;
    gratingSceneEngine.visualizeEachCompute = false;
    visualizationContrast = 1.0;
    [theSceneSequence] = gratingSceneEngine.compute(visualizationContrast);
    gratingSceneEngine.visualizeStaticFrame(... 
        theSceneSequence, ...
        'frameToVisualize', 1, ...
        'axesHandle', axLeft{idx});
    gratingSceneEngine.visualizeEachCompute = visualizeEachComputeSave;

    % Plot data and psychometric curve
    % with a marker size of 2.5
    questObj.plotMLE(2.5,'para',para(idx,:), 'axesHandle', axRight{idx});
    drawnow;
end
set(dataFig, 'Position',  [0, 0, 800, 800]);
saveas(dataFig,[figureFileBase '_Psychometric.tiff'],'tif');

% Convert returned log threshold to linear threshold
threshold = 10 .^ logThreshold;

%% Plot contrast sensitivity function
theCsfFig = figure();
loglog(spatialFreqs, 1 ./ threshold, '-ok', 'LineWidth', 2);
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

%% Return threshold values if requested
if (nargout > 0)
    thresholdRet = threshold;
end

end


