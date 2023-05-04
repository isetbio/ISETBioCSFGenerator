function t_spatialCSFmRGCMosaic
% Compute spatial CSF in different color directions, using the ON-center mRGCMosaics
%
% Description:
%    Use ISETBioCSFGenerator to run out CSFs in different color directions.
%
% See also: t_spatialCSFCMosaic, t_thresholdEngine, t_modulatedGratingsSceneGeneration,
%           t_chromaticThresholdContour, computeThresholdTAFC, computePerformanceTAFC
%

% History:
%   05/03/23  NPC   Wrote it

% Clear and close
clear; close all;

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
% further for this tutorial.  A vector length contrast of 0.08 should be
% OK.
rmsContrast = 0.08;
chromaDir = chromaDir / norm(chromaDir) * rmsContrast;
assert(abs(norm(chromaDir) - rmsContrast) <= 1e-10);


%% Create an mRGCMosaic-based neural response engine
%
% nreMidgetRGCMosaicSingleShot calculates the activation of an ON-center mRGCMosaic
% for stimuli consisting of a single frame and without eye movements
theNeuralComputePipelineFunction = @nreMidgetRGCMosaicSingleShot;

% Retrieve the default params for this engine
neuralResponsePipelineParams = theNeuralComputePipelineFunction();

% Modify certain params of interest
% 1. Crop the mRGCmosaic to size. Passing [] for size will not crop/.
% Passing an empty value for eccentricityDegs will crop
% the mosaic at its center.
neuralResponsePipelineParams.mRGCMosaicParams.cropParams = struct(...
    'sizeDegs', [], ...
    'eccentricityDegs', [] ...
);


% 2. If we want to use custom optics (not the optics that were used to optimize
% the mRGCMosaic), pass the optics here.
%neuralResponsePipelineParams.customOpticsToEmploy = oiCreate();

% 3. Set the input cone mosaic integration time
neuralResponsePipelineParams.mRGCMosaicParams.coneIntegrationTimeSeconds = 200/1000;

% 4. PRE and POST-CONE SUMMATION NOISE
% Pre-summation noise (here just the cone mosaic noise)i
neuralResponsePipelineParams.noiseParams.inputConeMosaicNoiseFlag = 'none';
neuralResponsePipelineParams.noiseParams.mRGCMosaicNoiseFlag = 'random';

% Post-cone summation noise, additive Gaussian noise with desired sigma
% If the input is cone modulations with respect to the background, the max
% response is 1.0, so the sigma should be scaled appropriate
% If the input is raw cone excitations, the sigma should be expressed in
% terms of cone excitations/integration time.

% Post-cone summation noise when mRGCs are integrating raw cone excitation signals
% neuralResponsePipelineParams.noiseParams.mRGCMosaicVMembraneGaussianNoiseSigma = 1e3 * 0.4;

% Post-cone summation noise when mRGCs are integrating  cone excitation modulations
neuralResponsePipelineParams.noiseParams.mRGCMosaicVMembraneGaussianNoiseSigma = 0.007;

% Instantiate theNeuralEngine!
theNeuralEngine = neuralResponseEngine(theNeuralComputePipelineFunction, neuralResponsePipelineParams);

% Matlab filename for saving computed data
matFileName = sprintf('mRGCMosaicSpatialCSF_eccDegs_%2.1f_%2.1f_coneContrasts_%2.2f_%2.2f_%2.2f_coneMosaicNoise_%s_mRGCMosaicNoise_%s.mat', ...
    neuralResponsePipelineParams.mRGCMosaicParams.eccDegs(1), ...
    neuralResponsePipelineParams.mRGCMosaicParams.eccDegs(2), ...
    chromaDir(1), chromaDir(2), chromaDir(3), ...
    neuralResponsePipelineParams.noiseParams.inputConeMosaicNoiseFlag, ...
    neuralResponsePipelineParams.noiseParams.mRGCMosaicNoiseFlag);

fprintf('Results will be saved in %s.\n', matFileName);

%% Instantiate a PoissonTAFC or a PcaSVMTAFC responseClassifierEngine
%
% Options:
%   'idealObserver' - ideal TAFC observer for Poisson limited signals
%   'computationalObsever' - SVM based learned classifier
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
        % coupled with a PCA operating on 2 components
        theClassifierEngine = responseClassifierEngine(@rcePcaSVMTAFC, ...
            struct(...
                'PCAComponentsNum', 2, ...          % number of PCs used for feature set dimensionality reduction
                'crossValidationFoldsNum', 10, ...  % employ a 10-fold cross-validated linear 
                'kernelFunction', 'linear', ...     % linear
                'classifierType', 'svm' ...         % binary SVM classifier
            ));

        % Train SVM classifier using a set of 512 noisy instances, 
        % Test performance using a set of 512 noisy instances
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
                       'slopeRangeHigh', 50, ...
                       'slopeDelta', 2.0);

% Parameter for running the QUEST+
% See t_thresholdEngine.m for more on options of the two different mode of
% operation (fixed numer of trials vs. adaptive)
contrastLevelsSampled = 5;
questEngineParams = struct(...
    'minTrial', contrastLevelsSampled*nTest, ...
    'maxTrial', contrastLevelsSampled*nTest, ...
    'numEstimator', 1, ...
    'stopCriterion', 0.05);


% We need the neuralResponseEngine to determine a stimulus FOV that is matched
% to the size of the inputConeMosaic. We also need theGratingSceneEngine to
% compute the theNullStimulusScene (which is used by the midgetRGCMosaic 
% neural response engine to compute mRGC responses based on cone mosaic contrast responses)
% To obtain these 2 engines, we call computeThresholdTAFC as follows:

% Just some dummy params for the grating.
dummySpatialFrequencyCPD = 4.0;
dummyFOVdegs = 2.0;
theGratingSceneEngine = createGratingScene(chromaDir, dummySpatialFrequencyCPD , ...
        'spatialEnvelope', 'rect', ...
        'fovDegs', dummyFOVdegs, ...
        'minPixelsNumPerCycle', 5, ...
        'pixelsNum', 512);

% Some low res quest params to run fast
questEngineParamsDummy = struct(...
    'minTrial', nTest, ...
    'maxTrial', nTest, ...
    'numEstimator', 1, ...
    'stopCriterion', 0.5);

% Run the dummy TAFC, just to generate theNeuralEngine and theGratingSceneEngine
computeThresholdTAFC(theGratingSceneEngine, theNeuralEngine, theClassifierEngine, ...
        classifierParams, thresholdParams, questEngineParamsDummy);

% Make the stimulus envelope equal to the the input cone mosaic size
theStimulusSpatialEnvelopeRadiusDegs = 0.5*max(theNeuralEngine.neuralPipeline.mRGCMosaic.inputConeMosaic.sizeDegs);

% Make the stimulusFOV 25% larger than the input cone mosaic size
% to avoid OI artifacts which arise when the test image has a mean radiance that is
% different than the mean radiance of the null stimulus, and which can
% result in the test stimulus being discriminable from the null stimulus
% just because of differences in the value with which the OI is padded at
% the edges
theStimulusFOVdegs = max(theNeuralEngine.neuralPipeline.mRGCMosaic.inputConeMosaic.sizeDegs)*1.25;

% Enough pixels so that the cone mosaic object does not complain that the
% OI resolution is too low compared to the cone aperture.
theStimulusPixelsNum = 512;

% Now that we have theGratingSceneEngine, use it to compute theNullStimulusScene
theGratingSceneEngine = createGratingScene(chromaDir, 5.0, ...
        'spatialEnvelope', 'rect', ...
        'fovDegs', theStimulusFOVdegs, ...
        'spatialEnvelopeRadiusDegs', theStimulusSpatialEnvelopeRadiusDegs, ...
        'minPixelsNumPerCycle', 5, ...
        'pixelsNum', theStimulusPixelsNum);

nullContrast = 0.0;
theNullStimulusSceneSequence = theGratingSceneEngine.compute(nullContrast);

% Save theNullStimulusScene in the passed neuralResponsePipelineParams, so
% that the neural engine can use it to compute mRGCmosaic responses
% operating on cone contrast (modulation) responses, instead of operating
% on raw cone excitation responses
neuralResponsePipelineParams.theNullStimulusScene = theNullStimulusSceneSequence{1};

% Re-instantiate theNeuralEngine so we can pass the updated
% neuralResponsePipelineParams (which now contains the null stimulus
% scene). This in turn, enables computation of mRGC mosaic responses based
% on the modulation of cone mosaic responses with respec to the cone mosaic
% response to the null stimulus.
theNeuralEngine.updateParamsStruct(neuralResponsePipelineParams);


%% Ready to compute thresholds at different spatial frequencies
% Choose the minSF so that it contains 1 full cycle within the smallest
% dimension of the input cone mosaic
minSF = 1/(min(theNeuralEngine.neuralPipeline.mRGCMosaic.inputConeMosaic.sizeDegs));
maxSF = 20;

% List of spatial frequencies to be tested.
spatialFreqs = logspace(log10(minSF), log10(maxSF), 10);

%% Compute threshold for each spatial frequency
% 
logThreshold = zeros(1, length(spatialFreqs));
theComputedQuestObjects = cell(1, length(spatialFreqs));
thePsychometricFunctions = cell(1, length(spatialFreqs));
theFittedPsychometricParams = cell(1, length(spatialFreqs));
theStimulusScenes = cell(1, length(spatialFreqs));

dataFig = figure();
for iSF = 1:length(spatialFreqs)
    % Create a static grating scene with a particular chromatic direction,
    % spatial frequency, and temporal duration
    theGratingSceneEngine = createGratingScene(chromaDir, spatialFreqs(iSF), ...
        'spatialEnvelope', 'rect', ...
        'fovDegs', theStimulusFOVdegs, ...
        'spatialEnvelopeRadiusDegs', theStimulusSpatialEnvelopeRadiusDegs, ...
        'minPixelsNumPerCycle', 5, ...
        'pixelsNum', theStimulusPixelsNum);
    
    [logThreshold(iSF), questObj, psychometricFunction, fittedPsychometricParams] = ...
        computeThresholdTAFC(theGratingSceneEngine, theNeuralEngine, theClassifierEngine, ...
        classifierParams, thresholdParams, questEngineParams);
    
    % Plot stimulus
    figure(dataFig);
    subplot(4, 6, iSF * 2 - 1);
    
    visualizationContrast = 1.0;
    [theSceneSequence] = theGratingSceneEngine.compute(visualizationContrast);
    theGratingSceneEngine.visualizeStaticFrame(theSceneSequence);

    % Plot data and psychometric curve 
    % with a marker size of 2.5
    subplot(4, 6, iSF * 2);
    questObj.plotMLE(2.5);
    drawnow;

    % Save for visualizing
    theComputedQuestObjects{iSF} = questObj;
    thePsychometricFunctions{iSF} = psychometricFunction;
    theFittedPsychometricParams{iSF}  = fittedPsychometricParams;
    theStimulusScenes{iSF} = theSceneSequence(1);
    

end
set(dataFig, 'Position',  [0, 0, 800, 800]);

% Convert returned log threshold to linear threshold
threshold = 10 .^ logThreshold;

%% Plot Contrast Sensitivity Function
theCsfFig = figure();
loglog(spatialFreqs, 1 ./ threshold, '-ok', 'LineWidth', 2);
xlabel('Spatial Frequency (cyc/deg)');
ylabel('Sensitivity');
set(theCsfFig, 'Position',  [800, 0, 600, 800]);

% Export computed data
save(matFileName, 'spatialFreqs', 'threshold', 'chromaDir', ...
    'theStimulusFOVdegs', 'theStimulusSpatialEnvelopeRadiusDegs', 'theStimulusScenes',...
    'theNeuralComputePipelineFunction', 'neuralResponsePipelineParams', ...
    'classifierChoice', 'classifierParams', 'thresholdParams', ...
    'theComputedQuestObjects', 'thePsychometricFunctions', 'theFittedPsychometricParams');
end
