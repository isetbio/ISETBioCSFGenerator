function t_spatialCSFmRGCMosaicOISequence
% Compute spatial CSF in different color directions, using the ON-center
% mRGCMosaics and dynamic stimuli (drifted/counter-phase modulated)
%
% Description:
%    Use ISETBioCSFGenerator to run out CSFs in different color directions 
%    using an mRGCMosaic neural respone engine with dynamic stimuli.
%
% See also: t_spatialCSFCMosaic, t_thresholdEngine, t_modulatedGratingsSceneGeneration,
%           t_chromaticThresholdContour, computeThresholdTAFC, computePerformanceTAFC
%

% History:
%   05/20/23  NPC   Wrote it

% Clear and close
clear; close all;

% Grating orientation
theStimulusOrientationDegs = 90;

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

% Set the RMS cone contrast of the stimulus. Things may go badly if you
% exceed the gamut of the monitor, so we are conservative and set this at a
% value that is within gamut of typical monitors and don't worry about it
% further for this tutorial.  A vector length contrast of 0.08 should be OK.
rmsContrast = 0.08;
chromaDir = chromaDir / norm(chromaDir) * rmsContrast;
assert(abs(norm(chromaDir) - rmsContrast) <= 1e-10);


% Set the cone mosaic intergration time. This will also be the stimulus
% frame duration. Here, set it to 50 mseconds
coneIntegrationTimeSeconds = 50/1000;

%% Create an mRGCMosaic-based neural response engine
%
% nreMidgetRGCMosaicOISequence calculates the activation of an ON-center mRGCMosaic
% for dynamic stimuli, consisting of multiple frames, without eye movements
theNeuralComputePipelineFunction = @nreMidgetRGCMosaicOISequence;

% Retrieve the default params for this engine
neuralResponsePipelineParams = theNeuralComputePipelineFunction();

% Modify certain params of interest
% 1. Select one of the pre-computed mRGC mosaics by specifying its
% eccentricityDegs & sizeDegs asnd its center type
neuralResponsePipelineParams.eccDegs = [0 0];
neuralResponsePipelineParams.sizeDegs =  [2 2];
neuralResponsePipelineParams.rgcType = 'ONcenterMidgetRGC';

% 2. We can crop the mRGCmosaic to some desired size. 
%     Passing [] for sizeDegs will not crop.
%     Passing [] for eccentricityDegs will crop the mosaic at its center.
neuralResponsePipelineParams.mRGCMosaicParams.cropParams = struct(...
    'sizeDegs', [1.5 1.5], ...
    'eccentricityDegs', [] ...
);


% 3. If we want to use custom optics (not the optics that were used to optimize
% the mRGCMosaic), pass the optics here.
%neuralResponsePipelineParams.customOpticsToEmploy = oiCreate();

% 4. Set the input cone mosaic integration time
neuralResponsePipelineParams.mRGCMosaicParams.coneIntegrationTimeSeconds = coneIntegrationTimeSeconds;

% 5. PRE and POST-CONE SUMMATION NOISE
% Pre- cone summation noise (ie Poisson noise)
neuralResponsePipelineParams.noiseParams.inputConeMosaicNoiseFlag = 'none';

% Post-cone summation noise (Gaussian noise)
neuralResponsePipelineParams.noiseParams.mRGCMosaicNoiseFlag = 'random';

% Post-cone summation noise is additive Gaussian noise with a desired
% sigma. When the input is raw cone excitations, the sigma should be expressed in
% terms of cone excitations/integration time. 
% When the input to mRGCs is cone modulations with respect to the background,
% which have a max amplitude of 1.0, the sigma should be scaled appropriately. 

% Post-cone summation noise when mRGCs are integrating raw cone excitation signals
% neuralResponsePipelineParams.noiseParams.mRGCMosaicVMembraneGaussianNoiseSigma = 1e3 * 0.1;

% Post-cone summation noise when mRGCs are integrating  cone excitation modulations
neuralResponsePipelineParams.noiseParams.mRGCMosaicVMembraneGaussianNoiseSigma = 0.015;

% Sanity check on the amount of mRGCMosaicVMembraneGaussianNoiseSigma for
% the specified neuralResponsePipelineParams.mRGCMosaicParams.inputSignalType 
switch (neuralResponsePipelineParams.mRGCMosaicParams.inputSignalType)
    case 'cone_modulations'
        % Ensure specificed mRGCVmembraneGaussianNoiseSigma is
        % appropriately scaled for cone modulations which are in the range
        % of [-1 1]
        if (neuralResponsePipelineParams.noiseParams.mRGCMosaicVMembraneGaussianNoiseSigma > 1)
            error('mRGC vMembrane Gaussian noise sigma (%f) is too large when operating on ''%s''.', ...
                neuralResponsePipelineParams.noiseParams.mRGCMosaicVMembraneGaussianNoiseSigma,...
                neuralResponsePipelineParams.mRGCMosaicParams.inputSignalType);
        end

    case 'cone_excitations'
        % Ensure specificed mRGCVmembraneGaussianNoiseSigma is
        % appropriately scaled for cone excitations
        if (neuralResponsePipelineParams.noiseParams.mRGCMosaicVMembraneGaussianNoiseSigma < 1)
            error('mRGC vMembrane Gaussian noise sigma (%f) is too small when operating on ''%s''.', ...
                neuralResponsePipelineParams.noiseParams.mRGCMosaicVMembraneGaussianNoiseSigma,...
                neuralResponsePipelineParams.mRGCMosaicParams.inputSignalType);
        end

    otherwise
        error('Unknown input signal type specified: ''%s''.', neuralResponsePipelineParams.mRGCMosaicParams.inputSignalType)
end

% Instantiate theNeuralEngine!
theNeuralEngine = neuralResponseEngine(theNeuralComputePipelineFunction, neuralResponsePipelineParams);



%% Instantiate a PoissonTAFC or a PcaSVMTAFC responseClassifierEngine
%
% Options:
%   'idealObserver' - ideal TAFC observer for Poisson limited signals
%   'computationalObsever' - SVM based learned classifier

% Since we are using a contrast-modulation based input to the mRGCmosaic,
% the Poisson noise is not valid. So we use a computationalObserver (SVM
% based)
classifierChoice = 'computationalObserver';

switch (classifierChoice) 
    case 'idealObserver'
        % PoissonTAFC makes decision by performing the Poisson likelihood ratio test
        % Also set up parameters associated with use of this classifier.
        theClassifierEngine = responseClassifierEngine(@rcePoissonTAFC);
        
        
        % Train classifier using 1 noise-free instance, 
        % Test performance using a set of 512 noisy instances
        nTest = 512;
        classifierParams = struct('trainFlag', 'none', ...
                                'testFlag', 'random', ...
                                'nTrain', 1, 'nTest', nTest);
    case 'computationalObserver'
        % Instantiate a computational observer consisting of a linear SVM 
        % coupled with a PCA operating on 16 components
        theClassifierEngine = responseClassifierEngine(@rcePcaSVMTAFC, ...
            struct(...
                'PCAComponentsNum', 4, ...         % number of PCs used for feature set dimensionality reduction
                'crossValidationFoldsNum', 10, ...  % employ a 10-fold cross-validated linear 
                'kernelFunction', 'linear', ...     % linear
                'classifierType', 'svm' ...         % binary SVM classifier
            ));

        % Train SVM classifier using a set of 4K noisy instances, 
        % Test performance using a set of 4K noisy instances
        nTrain = 1024*4;
        nTest = 1024*4;
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
                       'logThreshLimitHigh', 0.0, ...
                       'logThreshLimitDelta', 0.01, ...
                       'slopeRangeLow', 1, ...
                       'slopeRangeHigh', 200, ...
                       'slopeDelta', 0.25);

% Parameter for running the QUEST+
% See t_thresholdEngine.m for more on options of the two different mode of
% operation (fixed numer of trials vs. adaptive)

% Sample the contrast-response psychometric curve at 5 contrast levels
contrastLevelsSampled = 5;

questEngineParams = struct(...
    'minTrial', contrastLevelsSampled*nTest, ...
    'maxTrial', contrastLevelsSampled*nTest, ...
    'numEstimator', 1, ...
    'stopCriterion', 0.05);


% We need access to the generated neuralResponseEngine to determine a stimulus FOV that is matched
% to the size of the inputConeMosaic. We also need theGratingSceneEngine to
% compute the theNullStimulusScene (which is used by the midgetRGCMosaic 
% neural response engine to compute mRGC responses based on cone mosaic contrast responses)
% To obtain these 2 engines, we call computeThresholdTAFC as with some dummy params as follows:

% Just some dummy params for the grating.
dummySpatialFrequencyCPD = 1.0;
dummyFOVdegs = 1.0;
theDummyGratingSceneEngine = createGratingScene(chromaDir, dummySpatialFrequencyCPD , ...
        'spatialEnvelope', 'rect', ...
        'fovDegs', dummyFOVdegs, ...
        'minPixelsNumPerCycle', 5, ...
        'pixelsNum', 64);

% Some low res quest params to run fast
questEngineParamsDummy = struct(...
    'minTrial', 64, ...
    'maxTrial', 64, ...
    'numEstimator', 1, ...
    'stopCriterion', 0.5);

fprintf('Computing a dummy threshold to get access to theNeuralEngine\n');
% Run the dummy TAFC, just to generate theNeuralEngine and theGratingSceneEngine
computeThresholdTAFC(theDummyGratingSceneEngine, theNeuralEngine, theClassifierEngine, ...
        classifierParams, thresholdParams, questEngineParamsDummy);

% Having ran the computeThresholdTAFC() function, theNeuralEngine has been generated,
% so we can retrieve from it the size of the inputConeMosaic, and therefore match
% the stimulus spatial params to it as follows.
% We make the stimulusFOV 25% larger than the input cone mosaic size
% to avoid OI artifacts that arise when the test image has a mean radiance that is
% different than the mean radiance of the null stimulus, and which can
% result in the test stimulus being discriminable from the null stimulus
% just because of differences in the value with which the OI is padded at
% the edges
theStimulusFOVdegs = max(theNeuralEngine.neuralPipeline.mRGCMosaic.inputConeMosaic.sizeDegs)*1.25;
theStimulusSpatialEnvelopeRadiusDegs = 0.5*theStimulusFOVdegs;


% Enough pixels so that the cone mosaic object does not complain that the
% OI resolution is too low compared to the cone aperture.
theStimulusPixelsNum = 512;
minPixelsNumPerCycle = 8;



% Compute theNullStimulusScene (so we can express cone responses in terms
% of cone modulations)
nullContrast = 0.0;
nullSpatialFrequency = 1.0;
theNullGratingSceneEngine = createGratingScene(chromaDir, nullSpatialFrequency, ...
        'spatialEnvelope', 'rect', ...
        'orientation', theStimulusOrientationDegs, ...
        'fovDegs', theStimulusFOVdegs, ...
        'spatialEnvelopeRadiusDegs', theStimulusSpatialEnvelopeRadiusDegs, ...
        'minPixelsNumPerCycle', minPixelsNumPerCycle, ...
        'pixelsNum', theStimulusPixelsNum);
theNullStimulusSceneSequence = theNullGratingSceneEngine.compute(nullContrast);

% Save theNullStimulusScene in the neuralResponsePipelineParams, so
% that the neural engine can use it to compute mRGCmosaic responses
% operating on cone contrast (modulation) responses, instead of operating
% on raw cone excitation responses
neuralResponsePipelineParams.theNullStimulusScene = theNullStimulusSceneSequence{1};

% Dynamic stimulus parameters
thePresentationMode = 'drifted';
theTemporalFrequencyHz = 5.0;

% Match the frame duration to the cone integration time
theFrameDurationSeconds = coneIntegrationTimeSeconds;

% Make the stimulus last for 1 full temporal cycle
theStimulusDurationSeconds = 1.0/theTemporalFrequencyHz;


% And instruct the mRGCMosaic neural response engine to operate on cone
% modulations
neuralResponsePipelineParams.mRGCMosaicParams.inputSignalType = 'cone_modulations';

% Update theNeuralEngine with the new neuralResponsePipelineParams
theNeuralEngine.updateParamsStruct(neuralResponsePipelineParams);

% Generate Matlab filename for saving computed data
matFileName = sprintf('mRGCMosaicSpatialCSF_eccDegs_%2.1f_%2.1f_coneContrasts_%2.2f_%2.2f_%2.2f_OrientationDegs_%d_temporalFrequencyHz_%2.1f_inputSignal_%s_coneMosaicNoise_%s_mRGCMosaicNoise_%s.mat', ...
    neuralResponsePipelineParams.mRGCMosaicParams.eccDegs(1), ...
    neuralResponsePipelineParams.mRGCMosaicParams.eccDegs(2), ...
    chromaDir(1), chromaDir(2), chromaDir(3), ...
    theStimulusOrientationDegs, ...
    theTemporalFrequencyHz, ...
    regexprep(neuralResponsePipelineParams.mRGCMosaicParams.inputSignalType, '_+(\w)', '${upper($1)}'), ...
    neuralResponsePipelineParams.noiseParams.inputConeMosaicNoiseFlag, ...
    neuralResponsePipelineParams.noiseParams.mRGCMosaicNoiseFlag);

fprintf('Results will be saved in %s.\n', matFileName);

%% Ready to compute thresholds at a set of examined spatial frequencies
% Choose the minSF so that it contains 1 full cycle within the smallest
% dimension of the input cone mosaic
minSF = 0.5/(min(theNeuralEngine.neuralPipeline.mRGCMosaic.inputConeMosaic.sizeDegs));
maxSF = 60;
spatialFrequenciesSampled = 16;

% List of spatial frequencies to be tested.
spatialFreqs = [0 logspace(log10(minSF), log10(maxSF), spatialFrequenciesSampled)];

%% Compute threshold for each spatial frequency
% 
logThreshold = zeros(1, length(spatialFreqs));
theComputedQuestObjects = cell(1, length(spatialFreqs));
thePsychometricFunctions = cell(1, length(spatialFreqs));
theFittedPsychometricParams = cell(1, length(spatialFreqs));
theStimulusScenes = cell(1, length(spatialFreqs));

dataFig = figure();
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

    % Compute the threshold for our grating scene with the previously
    % defined neural and classifier engine.
    [logThreshold(iSF), questObj, psychometricFunction, fittedPsychometricParams] = ...
        computeThresholdTAFC(theDynamicGratingSceneEngine, theNeuralEngine, theClassifierEngine, ...
        classifierParams, thresholdParams, questEngineParams);
    
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


    % Also visualize the entire scene sequence
    theDynamicGratingSceneEngine.visualizeSceneSequence(...
        theSceneSequence, theSceneSequenceTemporalSupportSeconds, ...
        'videoFilename', sprintf('stimulus_%2.2fcpd', spatialFreqs(iSF)));

    % Save data for off-line visualizations
    theComputedQuestObjects{iSF} = questObj;
    thePsychometricFunctions{iSF} = psychometricFunction;
    theFittedPsychometricParams{iSF}  = fittedPsychometricParams;
    theStimulusScenes{iSF} = theSceneSequence(1);
end % iSF

set(dataFig, 'Position',  [0, 0, 800, 800]);

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

% Export computed data
save(matFileName, 'spatialFreqs', 'threshold', 'chromaDir', ...
    'theStimulusFOVdegs', 'theStimulusSpatialEnvelopeRadiusDegs', 'theStimulusScenes',...
    'theTemporalFrequencyHz', 'theFrameDurationSeconds',  'theStimulusDurationSeconds', ... 
    'theNeuralComputePipelineFunction', 'neuralResponsePipelineParams', ...
    'classifierChoice', 'classifierParams', 'thresholdParams', ...
    'theComputedQuestObjects', 'thePsychometricFunctions', 'theFittedPsychometricParams', '-v7.3');
end
