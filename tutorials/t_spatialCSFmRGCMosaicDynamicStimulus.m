function t_spatialCSFmRGCMosaicDynamicStimulus
% Compute spatial CSF in different color directions, using the ON-center
% mRGCMosaics and dynamic stimuli (drifted/counter-phase modulated)
%
% Description:
%    Use ISETBioCSFGenerator to run out CSFs in different color directions 
%    using an mRGCMosaic neural respone engine with dynamic stimuli.
%
% See also: t_spatialCSF, t_modulatedGratingsSceneGeneration,
%           t_chromaticThresholdContour, computeThreshold, computePerformance

% History:
%   05/20/23  NPC   Wrote it
%   12/23/24  dhb   Update for new architecture

% Clear and close
close all;

% Grab root folder for figure output
csfGeneratorRootPath = ISETBioCSFGeneratorRootPath;
figurePath = fullfile(csfGeneratorRootPath,'local','figures');
resultsPath = fullfile(csfGeneratorRootPath,'local','results');
if (~exist(figurePath,'dir'))
    mkdir(figurePath);
end
if (~exist(resultsPath,'dir'))
    mkdir(resultsPath);
end

% Set fastParameters that make this take less time
%
% Setting to false provides more realistic values for real work, but we
% try to keep the demo version relatively run time relatively short.
fastParameters = false;

% Grating orientation
theStimulusOrientationDegs = 0;

% Choose stimulus chromatic direction specified as a 1-by-3 vector
% of L, M, S cone contrast.  These vectors get normalized below, so only
% their direction matters in the specification.
stimType = 'luminance';
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

% Set the RMS cone contrast of the stimulus at max contrast. Things may go
% badly if you exceed the gamut of the monitor, so we are conservative and
% set this at a value that is within gamut of typical monitors and don't
% worry about it further for this tutorial.  A vector length contrast of
% 0.08 should be OK.
rmsContrast = 0.08;
chromaDir = chromaDir / norm(chromaDir) * rmsContrast;
assert(abs(norm(chromaDir) - rmsContrast) <= 1e-10);

% Set the cone mosaic intergration time. This will also be the stimulus
% frame duration. Here, set it to 50 mseconds.
%
% Match frame duration to integration time
coneIntegrationTimeSeconds = 50/1000;
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
noiseFreeParams.mRGCMosaicParams.rgcType = 'ONcenterMidgetRGC';

% 2. We can crop the mRGCmosaic to some desired size. 
%     Passing [] for sizeDegs will not crop.
%     Passing [] for eccentricityDegs will crop the mosaic at its center.
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
%neuralResponsePipelineParams.customOpticsToEmploy = oiCreate();

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
% the specified neuralResponsePipelineParams.mRGCMosaicParams.inputSignalType 
switch (noiseFreeParams.mRGCMosaicParams.inputSignalType)
    case 'cone_modulations'
        % Cone modulations which are in the range of [-1 1]
        if (noisyInstancesParams.sigma > 1)
            error('Gaussian noise sigma (%f) is too large when operating on ''%s''.', ...
                noisyInstancesParams.sigma,...
                noiseFreeParams.mRGCMosaicParams.inputSignalType);
        end

    case 'cone_excitations'
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
dummyNeuralEngine = neuralResponseEngine(noiseFreeComputeFunction, ...
    noisyInstancesComputeFunction, ...
    noiseFreeParams, ...
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

%% Instantiate the responseClassifierEngine
%
% rcePoisson makes decision by performing the Poisson likelihood ratio
% test. This is the ideal observer for the Poisson noice cone excitations
% illustrated in this script.  But you can run other rce's as well.
%    rcePoisson - signal known exactly Poission max likelihood
%    rceTemplateDistance - signal known exactly nearest L2 template
%                 distance.
%    rcePcaSVM  - support vector machine linear classifier after PCA.
%
% Also set up parameters associated with use of this classifier.
classifierEngine = 'rcePoisson';
switch (classifierEngine)
    case {'rcePoisson'}
        nTest = 128;
        classifierEngine = responseClassifierEngine(@rcePoisson);
        classifierUseParams = struct('trainFlag', 'none', ...
            'testFlag', 'random', ...
            'nTrain', 1, 'nTest', nTest);
        doValidationCheck = false;
    case {'rceTemplateDistance'}
        nTest = 128;
        classifierEngine = responseClassifierEngine(@rcePoisson);
        classifierUseParams = struct('trainFlag', 'none', ...
            'testFlag', 'random', ...
            'nTrain', 1, 'nTest', nTest);
        doValidationCheck = false;
    case 'rcePcaSVM'
        nTest = 512;
        pcaSVMParams = struct(...
            'PCAComponentsNum', 2, ...          % number of PCs used for feature set dimensionality reduction
            'crossValidationFoldsNum', 10, ...  % employ a 10-fold cross-validated linear
            'kernelFunction', 'linear', ...     % linear
            'classifierType', 'svm' ...         % binary SVM classifier
            );
        classifierEngine = responseClassifierEngine(@rcePcaSVM,pcaSVMParams);
        classifierUseParams = struct('trainFlag', 'none', ...
            'testFlag', 'random', ...
            'nTrain', 512, 'nTest', nTest);  
        doValidationCheck = false;
    otherwise
        error('Unsupported rce specified')
end

%% Parameters for threshold estimation/quest engine
% The actual threshold varies enough with the different engines that we
% need to adjust the contrast range that Quest+ searches over, as well as
% the range of psychometric function slopes. Threshold limits are computed
% as 10^-logThreshLimitVal.
thresholdParams = struct('logThreshLimitLow', 4.5, ...
                       'logThreshLimitHigh', 1, ...
                       'logThreshLimitDelta', 0.05, ...
                       'slopeRangeLow', 1, ...
                       'slopeRangeHigh', 100, ...
                       'slopeDelta', 1);

% Parameter for running the QUEST+
% See t_spatialCSF.m for more on options of the two different mode of
% operation (fixed numer of trials vs. adaptive)
%
% Sample the contrast-response psychometric curve at this number of
% contrast levels.
%
% Might want to up the number for the non-fastParameters case.
if (fastParameters)
    contrastLevelsSampled = 5;
else
    contrastLevelsSampled = 6;
end
questEngineParams = struct(...
    'minTrial', contrastLevelsSampled*nTest, ...
    'maxTrial', contrastLevelsSampled*nTest, ...
    'numEstimator', 1, ...
    'stopCriterion', 0.05);

% Contrast for the nullStimulusSceneSequence, typically zero
nullContrast = 0.0;

% Dynamic stimulus parameters
thePresentationMode = 'drifted';
theTemporalFrequencyHz = 5.0;

% This will run faster if you reduce theStimulusPixelsNum to something
% smaller, but then you will get an annoying limit about the precision
% of the cone aperture blurring.
if (fastParameters)
    theStimulusPixelsNum = 128;
    minPixelsNumPerCycle = 4;
    spatialFrequenciesSampled = 2;

    % Just two frames
    theStimulusDurationSeconds = 2*theFrameDurationSeconds;
else
    theStimulusPixelsNum = 512;
    minPixelsNumPerCycle = 8;
    spatialFrequenciesSampled = 6;

    % One full temporal cycle
    theStimulusDurationSeconds = 1.0/theTemporalFrequencyHz;
end

% Generate Matlab filename for saving computed data
% 
% Code that actually does the save commented out at the end.
matFileName = sprintf('mRGCMosaicSpatialCSF_eccDegs_%2.1f_%2.1f_coneContrasts_%2.2f_%2.2f_%2.2f_OrientationDegs_%d_temporalFrequencyHz_%2.1f_inputSignal_%s_noiseSd_%0.3f.mat', ...
    noiseFreeParams.mRGCMosaicParams.eccDegs(1), ...
    noiseFreeParams.mRGCMosaicParams.eccDegs(2), ...
    chromaDir(1), chromaDir(2), chromaDir(3), ...
    theStimulusOrientationDegs, ...
    theTemporalFrequencyHz, ...
    regexprep(noiseFreeParams.mRGCMosaicParams.inputSignalType, '_+(\w)', '${upper($1)}'), ...
    noisyInstancesParams.sigma);

%% Ready to compute thresholds at a set of spatial frequencies
minSF = 2;
maxSF = 10;
spatialFreqs = [0 logspace(log10(minSF), log10(maxSF), spatialFrequenciesSampled)];

%% Compute threshold for each spatial frequency
% 
logThreshold = zeros(1, length(spatialFreqs));
theComputedQuestObjects = cell(1, length(spatialFreqs));
thePsychometricFunctions = cell(1, length(spatialFreqs));
theFittedPsychometricParams = cell(1, length(spatialFreqs));
theStimulusScenes = cell(1, length(spatialFreqs));

dataFig = figure;
set(dataFig, 'Position',  [0, 0, 800, 800]);
plotCols = 8;
plotRows = ceil((length(spatialFreqs)*2) / plotCols);

for iSF = 1:length(spatialFreqs)
    if (spatialFreqs(iSF) == 0)
        % For 0 c/deg, we do a counter-phase modulation
        thePresentationModeForThisSF = 'counter phase modulated';
    else
        thePresentationModeForThisSF = thePresentationMode;
    end

    fprintf('Computing contrast sensitivity at %2.1f c/deg using a %s stimulus\n', ...
        spatialFreqs(iSF), thePresentationModeForThisSF);

    % Generate the dynamic grating scene engine for this SF
    theDynamicGratingSceneEngine = createGratingScene(chromaDir, spatialFreqs(iSF), ...
        'spatialEnvelope', 'rect', ...
        'orientation', theStimulusOrientationDegs, ...
        'fovDegs', theStimulusFOVdegs, ...
        'spatialEnvelopeRadiusDegs', theStimulusSpatialEnvelopeRadiusDegs, ...
        'minPixelsNumPerCycle', minPixelsNumPerCycle, ...
        'pixelsNum', theStimulusPixelsNum, ...
        'presentationMode', thePresentationModeForThisSF, ...
        'temporalFrequencyHz', theTemporalFrequencyHz, ...
        'spatialPhaseAdvanceDegs', 360*(theFrameDurationSeconds*theTemporalFrequencyHz), ...
        'duration', theStimulusDurationSeconds);

    % First time through, need to compute the null scene and set up the nre
    % with that null scene.  But we don't need to do this separately for
    % each spatial frequency, just once in the whole script. 
    if (iSF == 1)
        % Compute the nullStimulusSceneSequece and put it in the
        % neuralResponsePipelineParams, so that the neural engine can use
        % it to compute mRGCmosaic responses operating on cone contrast
        % (modulation) responses, instead of operating on raw cone
        % excitation responses
        nullStimulusSceneSequence = theDynamicGratingSceneEngine.compute(nullContrast);
        noiseFreeParams.mRGCMosaicParams.inputSignalType = 'cone_modulations';
        noiseFreeParams.nullStimulusSceneSequence = nullStimulusSceneSequence;

        % Re-instantiate theNeuralEngine with the null scene now computed
        theNeuralEngine = neuralResponseEngine(noiseFreeComputeFunction, ...
            noisyInstancesComputeFunction, ...
            noiseFreeParams, ...
            noisyInstancesParams);
    end

    % Compute the threshold for our grating scene with the previously
    % defined neural and classifier engine.
    [logThreshold(iSF), questObj, psychometricFunction, fittedPsychometricParams] = ...
        computeThreshold(theDynamicGratingSceneEngine, theNeuralEngine, classifierEngine, ...
        classifierUseParams, thresholdParams, questEngineParams, 'TAFC', true);
    
    % Plot stimulus
    figure(dataFig);
    subplot(plotRows, plotCols, iSF * 2 - 1);
    
    visualizationContrast = 1.0;
    [theSceneSequence,theSceneSequenceTemporalSupportSeconds] = theDynamicGratingSceneEngine.compute(visualizationContrast);
    theDynamicGratingSceneEngine.visualizeStaticFrame(theSceneSequence);

    % Plot data and psychometric curve 
    % with a marker size of 2.5
    figure(dataFig);
    subplot(plotRows, plotCols, iSF * 2);
    questObj.plotMLE(2.5);
    drawnow;

    % Visualize the entire scene sequence
    theDynamicGratingSceneEngine.visualizeSceneSequence( ...
        theSceneSequence, theSceneSequenceTemporalSupportSeconds, ...
        'videoFilename', fullfile(figurePath,sprintf('stimulus_%2.2fcpd', spatialFreqs(iSF))));

    % Save data for off-line visualizations
    theComputedQuestObjects{iSF} = questObj;
    thePsychometricFunctions{iSF} = psychometricFunction;
    theFittedPsychometricParams{iSF}  = fittedPsychometricParams;
    theStimulusScenes{iSF} = theSceneSequence(1);
end % iSF

% Convert returned log threshold to linear threshold
threshold = 10 .^ logThreshold;

%% Plot Contrast Sensitivity Function
theCsfFig = figure();
theSpatialFreqs = spatialFreqs;
if (theSpatialFreqs(1) == 0)
    theSpatialFreqs(1) = 0.01;
end

loglog(theSpatialFreqs, 1 ./ threshold, '-ok', 'LineWidth', 2);
xlabel('Spatial Frequency (cyc/deg)');
ylabel('Sensitivity');
set(theCsfFig, 'Position',  [800, 0, 600, 800]);

%% Export computed data. 
fprintf('Results will be saved in %s.\n', matFileName);
save(fullfile(resultsPath,matFileName), 'spatialFreqs', 'threshold', 'chromaDir', ...
    'theStimulusFOVdegs', 'theStimulusSpatialEnvelopeRadiusDegs', 'theStimulusScenes',...
    'theTemporalFrequencyHz', 'theFrameDurationSeconds',  'theStimulusDurationSeconds', ... 
    'noiseFreeComputeFunction', 'noiseFreeParams', ...
    'noisyInstancesComputeFunction', 'noisyInstancesParams', ...
    'classifierEngine', 'classifierUseParams', 'thresholdParams', ...
    'theComputedQuestObjects', 'thePsychometricFunctions', 'theFittedPsychometricParams', '-v7.3');
end
