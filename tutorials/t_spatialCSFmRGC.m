function t_spatialCSFmRGC(options)
% Compute spatial CSF in different color directions, using the ON-center
% mRGCMosaics. This script illustrates use with and without meta contrast
% method for mRGC calculations.
%
% Description:
%    Use ISETBioCSFGenerator to run out CSFs in different color directions 
%    using mRGCMosaic neural respone engine and (optionally) the meta contrast method.
%
% See also: t_spatialCSFcMosaic, t_metaContrastCSFcMosaic, t_spatialCSFmRGCDynamicStimulus

% History:
%   11/02/2024  FH   Adopted based on t_spatialCSFmRGCMosaic.m
%   12/24/2024  dhb  Update for new architecture

arguments
    % Run the validation check?  This gets overridden to empty if other
    % options change the conditions so that the validation data don't
    % apply.
    options.validationThresholds (1,:) double = [];

    % Apply a filter to the spectra before computing responses?  See
    % t_spatialCSFcMosaicFilter
    options.filter (1,1) = struct('spectralSupport',[],'transmission',[]);

    % Use meta contrast method to speed things up?
    options.useMetaContrast (1,1) logical = true;

    % Use cone contrast rather than cone excitations
    options.useConeContrast (1,1) logical = false;

    % Use fixational eye movements?
    options.useFixationalEMs (1,1) logical = false;

    % Choose noise free neural model
    %   Choices: 'excitationsCmosaic'
    %            'sceneAsResponses'
    options.whichNoiseFreeNre (1,:) char  = 'excitationsCmosaic'

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
end

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
whichNoiseFreeNre = options.whichNoiseFreeNre;
whichNoisyInstanceNre = options.whichNoisyInstanceNre;
whichClassifierEngine = options.whichClassifierEngine;
validationThresholds = options.validationThresholds;

% Clear out stay figures
close all;

%% Figure output base name
figureFileBase = fullfile(projectBaseDir,'local',mfilename,'figures', ...
    sprintf('%s_Meta_%d_ConeContrast_%d_FEMs_%d_%s_%s_%s', mfilename, ...
    useMetaContrast,useConeContrast,useFixationalEMs,whichNoiseFreeNre,whichNoisyInstanceNre,whichClassifierEngine));
figureTypeStr = '.tif';

% Set fastParameters that make this take less time
%
% Setting to false provides more realistic values for real work, but we
% try to keep the demo version run time relatively short.
fastParameters = true;

% Choose stimulus chromatic direction specified as a 1-by-3 vector
% of L, M, S cone contrast.  These vectors get normalized below, so only
% their direction matters in the specification.
stimType = 'achromatic';
switch (stimType)
    case 'achromatic'
        chromaDir = [1.0, 1.0, 1.0]';
    case 'luminance'
        chromaDir = [1.0, 1.0, 0.0]';
    case 'red-green'
        chromaDir = [1.0, -1.0, 0.0]';
    case 'L-isolating'
        chromaDir = [1.0, 0.0, 0.0]';
end

% Set the RMS cone contrast of the stimulus. Things may go badly if you
% exceed the gamut of the monitor, so we are conservative and set this at a
% value that is within gamut of typical monitors and don't worry about it
% further for this tutorial.  A vector length contrast of 0.08 should be OK.
rmsContrast = 0.1;
chromaDir = chromaDir / norm(chromaDir) * rmsContrast;
assert(abs(norm(chromaDir) - rmsContrast) <= 1e-10);

% Set the cone mosaic intergration time. This will also be the stimulus
% frame duration. Here, set it to 50 mseconds.
%
% Match frame duration to integration time
coneIntegrationTimeSeconds = 200/1000;
theFrameDurationSeconds = coneIntegrationTimeSeconds;

%% Create an mRGCMosaic-based neural response engine
%
% nreNoiseFreeMidgetRGCMosaicOISequence calculates the activation of an ON-center mRGCMosaic
% for dynamic stimuli, consisting of multiple frames
noiseFreeComputeFunction = @nreNoiseFreeMidgetRGCMosaic;
noiseFreeParams = noiseFreeComputeFunction();
noisyInstancesComputeFunction = @nreNoisyInstancesGaussian;
noisyInstancesParams = noisyInstancesComputeFunction();

% Modify certain parameters of interest
%
% 1. Select one of the pre-computed mRGC mosaics by specifying its
% eccentricityDegs & sizeDegs asnd its center type
noiseFreeParams.mRGCMosaicParams.eccDegs = [0 0];
noiseFreeParams.mRGCMosaicParams.sizeDegs = [2.0 2.0];
noiseFreeParams.mRGCMosaicParams.inputSignalType = 'coneContrast';
noiseFreeParams.mRGCMosaicParams.rgcType = 'ONcenterMidgetRGC';

% 2. We can crop the mRGCmosaic to some desired size. 
%    Passing [] for sizeDegs will not crop.
%    Passing [] for eccentricityDegs will crop the mosaic at its center.
if (fastParameters)
    cropSize = [0.5 0.5];
else
    cropSize = [1.5 1.5];
end
noiseFreeParams.mRGCMosaicParams.cropParams = struct(...
    'sizeDegs', cropSize, ...
    'eccentricityDegs', [] ...
);

% 3. If we want to use custom optics (not the optics that were used to optimize
% the mRGCMosaic), pass the optics here. This is commented out but
% illustrates where and how a custom oi would be passed.
%noiseFreeParams.customOpticsToEmploy = oiCreate();

% 4. Set the input cone mosaic integration time
noiseFreeParams.mRGCMosaicParams.coneIntegrationTimeSeconds = coneIntegrationTimeSeconds;

% Post-cone summation noise is additive Gaussian noise with a desired
% sigma. When the input is raw cone excitations, the sigma should be
% expressed in terms of cone excitations/integration time. When the input
% to mRGCs is cone modulations with respect to the background, which have a
% max amplitude of 1.0, the sigma should be scaled appropriately.
% So if your mRGC mosiac were operating directly on cone excitations, a
% different value more like 100 or 1000 would be reasonable, depending on
% light level.
noisyInstancesParams.sigma = 0.015;

% Sanity check on the amount of mRGCMosaicVMembraneGaussianNoiseSigma for
% the specified noiseFreeParams.mRGCMosaicParams.inputSignalType 
switch (noiseFreeParams.mRGCMosaicParams.inputSignalType)
    case 'coneContrast'
        % Cone modulations which are in the range of [-1 1]
        if (noisyInstancesParams.sigma > 1)
            error('Gaussian noise sigma (%f) is too large when operating on ''%s''.', ...
                noisyInstancesParams.sigma,...
                noiseFreeParams.mRGCMosaicParams.inputSignalType);
        end

    case 'coneExcitations'
        % Excitations are unlikely to be of order 1.
        if (noisyInstancesParams.sigma < 1)
            error('Gaussian noise sigma (%f) is too small when operating on ''%s''.', ...
                noisyInstancesParams.sigma,...
                noiseFreeParams.mRGCMosaicParams.inputSignalType);
        end

    otherwise
        error('Unknown mRGC signal type specified: ''%s''.', noiseFreeParams.mRGCMosaicParams.inputSignalType)
end

%% Instantiate a dummy neural engine.
%
% We use this to figure out some sizes, as described and done just below.
dummyNoiseFreeParams = noiseFreeParams;
dummyNoiseFreeParams.mRGCMosaicParams.inputSignalType = 'coneExcitations';
dummyNeuralEngine = neuralResponseEngine(noiseFreeComputeFunction, ...
    noisyInstancesComputeFunction, ...
    dummyNoiseFreeParams, ...
    noisyInstancesParams);

% We need access to the generated neuralResponseEngine to determine a stimulus FOV that is matched
% to the size of the inputConeMosaic. This is because the mRGC object creates a cone mosaic
% that is big enough to allow computation of the surrounds, and we don't know how big this
% will be until we try it.
%
% We make the stimulusFOV 25% larger than the input cone mosaic size
% to avoid OI artifacts that arise when the test image has a mean radiance that is
% different than the mean radiance of the null stimulus, and which can
% result in the test stimulus being discriminable from the null stimulus
% just because of differences in the value with which the OI is padded at
% the edges.
%
% Just some dummy params for the grating so we can get the mosaic size.
dummySpatialFrequencyCPD = 1.0;
dummyFOVdegs = 1.0;
theDummyGratingSceneEngine = createGratingScene(chromaDir, dummySpatialFrequencyCPD , ...
        'spatialEnvelope', 'rect', ...
        'fovDegs', dummyFOVdegs, ...
        'minPixelsNumPerCycle', 5, ...
        'pixelsNum', 256);
[dummySceneSequence,dummyTemporalSupport] = theDummyGratingSceneEngine.compute(0);
dummyNeuralEngine.computeNoiseFree(dummySceneSequence,dummyTemporalSupport);
theStimulusFOVdegs = max(dummyNeuralEngine.neuralPipeline.noiseFreeResponse.mRGCMosaic.inputConeMosaic.sizeDegs)*1.25;
theStimulusSpatialEnvelopeRadiusDegs = 0.5*theStimulusFOVdegs;
clear dummyNeuralEngine

%% Create real neural engine
theNeuralEngine = neuralResponseEngine(noiseFreeComputeFunction, ...
    noisyInstancesComputeFunction, ...
    noiseFreeParams, ...
    noisyInstancesParams);

%% Instantiate a PoissonTAFC or a PcaSVMTAFC responseClassifierEngine
%
% Options:
%   'idealObserver' - ideal TAFC observer for Poisson limited signals
%   'computationalObsever' - SVM based learned classifier

% Since we are using a contrast-modulation based input to the mRGCmosaic,
% the Poisson noise is not accurate. So we use a computationalObserver (SVM
% based)
classifierChoice = 'computationalObserver';

switch (classifierChoice) 
    case 'idealObserver'
        % PoissonTAFC makes decision by performing the Poisson likelihood ratio test
        % Also set up parameters associated with use of this classifier.
        theClassifierEngine = responseClassifierEngine(@rcePoisson);
        
        % Train classifier using 1 noise-free instance, 
        % Test performance using a set of 512 noisy instances
        if (fastParameters)
            nTest = 32;
        else
            nTest = 1024;
        end
        nTest = 512;
        classifierParams = struct('trainFlag', 'none', ...
                                'testFlag', 'random', ...
                                'nTrain', 1, 'nTest', nTest);
    
    case 'computationalObserver'
        % Handle fastParameters
        if (fastParameters)
            crossValidationFolds = 2;
            nTrain = 64;
            nTest = 32;
        else
            crossValidationFolds = 10;
            nTrain = 1024;
            nTest = 1024;
        end

        % Instantiate a computational observer consisting of a linear SVM 
        % coupled with a PCA operating on 2 components
        theClassifierEngine = responseClassifierEngine(@rcePcaSVM, ...
            struct(...
                'PCAComponentsNum', 2, ...          % number of PCs used for feature set dimensionality reduction
                'crossValidationFoldsNum', crossValidationFolds, ...  % employ a 10-fold cross-validated linear 
                'kernelFunction', 'linear', ...     % linear
                'classifierType', 'svm' ...         % binary SVM classifier
            ));

        % Train SVM classifier using a set of 4K noisy instances, 
        % Test performance using a set of 4K noisy instances
        classifierParams = struct('trainFlag', 'random', ...
                                'testFlag', 'random', ...
                                'nTrain', nTest, 'nTest', nTrain);
                        
    otherwise
        error('Unknown classifier: ''%s''.', classifierChoice);
end

%% Parameters for threshold estimation/quest engine
% The actual threshold varies enough with the different engines that we
% need to adjust the contrast range that Quest+ searches over, as well as
% the range of psychometric function slopes. Threshold limits are computed
% as 10^-logThreshLimitVal.
thresholdParams = struct('logThreshLimitLow', 2.5, ...
                       'logThreshLimitHigh', 0, ...
                       'logThreshLimitDelta', 0.05, ...
                       'slopeRangeLow', 1, ...
                       'slopeRangeHigh', 200, ...
                       'slopeDelta', 1.0);

% Parameters for running the QUEST+
%
% See t_spatialCSF.m for more on options of the two different mode of
% operation (fixed numer of trials vs. adaptive)
%
% Sample the contrast-response psychometric curve at this number of
% contrast levels.
if (fastParameters)
    contrastLevelsSampled = 8;
else
    contrastLevelsSampled = 10;
end
questEngineParams = struct(...
    'minTrial', contrastLevelsSampled*nTest, ...
    'maxTrial', contrastLevelsSampled*nTest, ...
    'numEstimator', 1, ...
    'stopCriterion', 0.05);

% This will run faster if you reduce theStimulusPixelsNum to something
% smaller, but then you might get an annoying limit about the precision
% of the cone aperture blurring.
if (fastParameters)
    theStimulusPixelsNum = 256;
    minPixelsNumPerCycle = 4;
else
    theStimulusPixelsNum = 512;
    minPixelsNumPerCycle = 12;
end



%% With access to theGratingSceneEngine, we can compute theNullStimulusScene
nullContrast = 0.0;
gratingEngineParams = struct( ...
    'spatialEnvelope', 'soft', ...
    'orientation', theStimulusOrientationDegs, ...
    'fovDegs', theStimulusFOVdegs, ...
    'spatialEnvelopeRadiusDegs', theStimulusSpatialEnvelopeRadiusDegs, ...
    'minPixelsNumPerCycle', minPixelsNumPerCycle, ...
    'pixelsNum', theStimulusPixelsNum ...
);
nullGratingSceneEngine = createGratingScene(chromaDir, dummySpatialFrequencyCPD, ...
        gratingEngineParams);
nullStimulusSceneSequence = nullGratingSceneEngine.compute(nullContrast);
noiseFreeParams.nullStimulusSceneSequence = nullStimulusSceneSequence;

% And instruct the mRGCMosaic neural response engine to operate on cone
% modulations
noiseFreeParams.mRGCMosaicParams.inputSignalType = 'coneContrast';

% Update theNeuralEngine with the new noiseFreeParams
theNeuralEngine.updateParamsStruct(noiseFreeParams,noisyInstancesParams);

% Generate Matlab filename for saving computed data.
% 
% Code that actually does the save commented out at the end.
projectBaseDir = ISETBioCSFGeneratorRootPath;
matFileName = sprintf('mRGCMosaicSpatialCSF_eccDegs_%2.1f_%2.1f_coneContrasts_%2.2f_%2.2f_%2.2f_OrientationDegs_%d_inputSignal_%s_noiseSd_%0.3f.mat', ...
    noiseFreeParams.mRGCMosaicParams.eccDegs(1), ...
    noiseFreeParams.mRGCMosaicParams.eccDegs(2), ...
    chromaDir(1), chromaDir(2), chromaDir(3), ...
    theStimulusOrientationDegs, ...
    regexprep(noiseFreeParams.mRGCMosaicParams.inputSignalType, '_+(\w)', '${upper($1)}'), ...
    noisyInstancesParams.sigma);

%% Ready to compute thresholds at a set of examined spatial frequencies
% Choose the minSF so that it contains 1 full cycle within the smallest
% dimension of the input cone  mosaic
minSF = 0.5/theStimulusFOVdegs;
maxSF = 60;
if (fastParameters)
    spatialFrequenciesSampled = 8;
else
    spatialFrequenciesSampled = 16;
end

% List of spatial frequencies to be tested.
spatialFreqs = logspace(log10(minSF), log10(maxSF), spatialFrequenciesSampled);

% Create the sceMetaContrast scene engine
if (useMetaContrast)
    metaSceneEngineParams = sceMetaContrast;
    metaSceneEngine = sceneEngine(@sceMetaContrast,metaSceneEngineParams);
end

%% Compute threshold for each spatial frequency
% 
logThreshold = zeros(1, length(spatialFreqs));
theComputedQuestObjects = cell(1, length(spatialFreqs));
thePsychometricFunctions = cell(1, length(spatialFreqs));
theFittedPsychometricParams = cell(1, length(spatialFreqs));
theStimulusScenes = cell(1, length(spatialFreqs));

dataFig = figure();
plotRows = 4;
plotCols = 8;
for iSF = 1:length(spatialFreqs)
    % Create a static grating scene with a particular chromatic direction,
    % spatial frequency, orientation, FOV, and size
    gratingSceneEngine = createGratingScene(chromaDir, spatialFreqs(iSF), gratingEngineParams);

    if (useMetaContrast)
        % Create nreMetaContrast using the actual scene and neural engines
        metaNeuralResponseEngineNoiseFreeParams = nreNoiseFreeMetaContrast;
        metaNeuralResponseEngineNoiseFreeParams.contrast0 = 0;
        metaNeuralResponseEngineNoiseFreeParams.contrast1 = 1;
        metaNeuralResponseEngineNoiseFreeParams.sceneEngine = gratingSceneEngine;
        metaNeuralResponseEngineNoiseFreeParams.neuralEngine = theNeuralEngine;

        metaNeuralResponseEngineNoisyInstanceParams = nreNoisyInstancesMetaContrast;
        metaNeuralResponseEngineNoisyInstanceParams.neuralEngine = theNeuralEngine;
        theMetaNeuralEngine = neuralResponseEngine(@nreNoiseFreeMetaContrast, ...
            @nreNoisyInstancesMetaContrast, ...
            metaNeuralResponseEngineNoiseFreeParams, ...
            metaNeuralResponseEngineNoisyInstanceParams);

        % Compute the threshold for our grating scene with the previously
        % defined neural and classifier engine.
        [logThreshold(iSF), questObj, psychometricFunction, fittedPsychometricParams] = ...
            computeThreshold(metaSceneEngine, theMetaNeuralEngine, theClassifierEngine, ...
            classifierParams, thresholdParams, questEngineParams,'TAFC', true);
    else
        [logThreshold(iSF), questObj, psychometricFunction, fittedPsychometricParams] = ...
            computeThreshold(gratingSceneEngine, theNeuralEngine, theClassifierEngine, ...
            classifierParams, thresholdParams, questEngineParams,'TAFC', true);
    end
    
    % Plot stimulus
    figure(dataFig);
    subplot(plotRows, plotCols, iSF * 2 - 1);
    
    visualizationContrast = 1.0;
    [theSceneSequence] = gratingSceneEngine.compute(visualizationContrast);
    gratingSceneEngine.visualizeStaticFrame(theSceneSequence);

    % Plot data and psychometric curve 
    % with a marker size of 2.5
    figure(dataFig);
    subplot(plotRows, plotCols, iSF * 2);
    questObj.plotMLE(2.5);
    drawnow;

    % Save data for off-line visualizations
    theComputedQuestObjects{iSF} = questObj;
    thePsychometricFunctions{iSF} = psychometricFunction;
    theFittedPsychometricParams{iSF}  = fittedPsychometricParams;
    theStimulusScenes{iSF} = theSceneSequence(1);
end

% Pretty up data figure
set(dataFig, 'Position',  [0, 0, 800, 800]);

% Save the figure as a PDF
% set(dataFig, 'PaperSize', [30, 20]);
saveas(dataFig, fullfile(projectBaseDir,'local',myName, ...
    [matFileName(1:10),'PMF', matFileName(21:end-4), figureTypeStr]));

% Convert returned log threshold to linear threshold
threshold = 10 .^ logThreshold;
elapsedTime = toc;

%% Plot Contrast Sensitivity Function
theCsfFig = figure();
loglog(spatialFreqs, 1 ./ threshold, '-ok', 'LineWidth', 2);
xlabel('Spatial Frequency (cyc/deg)');
ylabel('Sensitivity');
set(theCsfFig, 'Position',  [800, 0, 600, 800]);

% Save the figure as a PDF
saveas(theCsfFig, fullfile(projectBaseDir,'local',myName,[matFileName(1:end-4), figureTypeStr]));

%% Export computed data. 
%
% Commented out. If you want to save data like this from a tutorial, put it
% into a local directory and make sure that is gitignored so it doesn't
% clog up the repository.
%
fprintf('Results will be saved in %s.\n', fullfile(projectBaseDir,'local',myName,matFileName));
save(fullfile(projectBaseDir,'local',myName,matFileName), 'spatialFreqs', 'threshold', 'chromaDir', ...
    'theStimulusFOVdegs', 'theStimulusSpatialEnvelopeRadiusDegs', ...
    'noiseFreeComputeFunction', 'noiseFreeParams', ...
    'noisyInstancesComputeFunction', 'noisyInstancesParams', ...    
    'classifierChoice', 'classifierParams', 'thresholdParams', ...
    'theComputedQuestObjects', 'thePsychometricFunctions', 'theFittedPsychometricParams','elapsedTime','-v7.3');
end
