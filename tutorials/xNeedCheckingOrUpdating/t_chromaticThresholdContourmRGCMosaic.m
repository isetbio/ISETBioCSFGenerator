function t_chromaticThresholdContourmRGCMosaic
% Compute isothreshold contour in different color directions, using the ON-center mRGCMosaics
%
% Description:
%    Use ISETBioCSFGenerator to run out an isothreshold contour in the LM
%    contrast plane using mRGCMosaic neural respone engines.
%
% See also: t_spatialCSFCMosaic, t_thresholdEngine, t_modulatedGratingsSceneGeneration,
%           t_chromaticThresholdContour, computeThresholdTAFC, computePerformanceTAFC
%

% History:
%   05/05/23  NPC   Wrote it

% Clear and close
clear; close all;

% Set fastParameters that make this take less time
%
% Setting to false provides more realistic values for real work, but we
% try to keep the demo version run time relatively short.
%
% Note that with the fast parameters, this doesn't really show the
% way that the ellipse rotationes with noise at the mRGCs.  
fastParameters = true;

% Choose stimulus spatial frequency, orientation, and spatial phase
theStimulusSpatialFrequencyCPD = 0;
theStimulusSpatialPhaseDegs = 0;
theStimulusOrientationDegs = 0;

% Set the RMS cone contrast of the stimulus. Things may go badly if you
% exceed the gamut of the monitor, so we are conservative and set this at a
% value that is within gamut of typical monitors and don't worry about it
% further for this tutorial.  A vector length contrast of 0.08 should be OK.
rmsContrast = 0.08;

% List of LM chromatic directions to be tested 
if (fastParameters)
    nChromaticDirections = 8;
else
    nChromaticDirections = 16;
end
for ii = 1:nChromaticDirections
    theta = (ii-1)/nChromaticDirections*2*pi;
    theChromaticDirections(:,ii) = [cos(theta) sin(theta) 0]';
    theChromaticDirections(:,ii) = theChromaticDirections(:,ii) / norm(theChromaticDirections(:,ii)) * rmsContrast;
    assert(abs(norm(theChromaticDirections(:,ii)) - rmsContrast) <= 1e-10);
end

%% Create an mRGCMosaic-based neural response engine
%
% nreMidgetRGCMosaicSingleShot calculates the activation of an ON-center mRGCMosaic
% for stimuli consisting of a single frame and without eye movements
theNeuralComputePipelineFunction = @nreMidgetRGCMosaicSingleShot;

% Retrieve the default params for this engine
neuralResponsePipelineParams = theNeuralComputePipelineFunction();

% Modify certain params of interest
% 1. We can crop the mRGCmosaic to some desired size. Passing [] for size will not crop.
% Passing an empty value for eccentricityDegs will crop the mosaic at its center.
neuralResponsePipelineParams.mRGCMosaicParams.cropParams = struct(...
    'sizeDegs', [], ...
    'eccentricityDegs', [] ...
);

% 2. If we want to use custom optics (not the optics that were used to optimize
% the mRGCMosaic), pass the optics here.
%neuralResponsePipelineParams.customOpticsToEmploy = oiCreate();

% 3. Set the input cone mosaic integration time
neuralResponsePipelineParams.mRGCMosaicParams.coneIntegrationTimeSeconds = 200/1000;

% 4. mRGCs operating either on 'cone_excitations' or on 'cone_modulations'
%    If 'cone modulations' is selected, you must provide the null scen.
%    This is done later on.
neuralResponsePipelineParams.mRGCMosaicParams.inputSignalType = 'cone_modulations';
neuralResponsePipelineParams.mRGCMosaicParams.inputSignalType = 'cone_excitations';

% 5. PRE and POST-CONE SUMMATION NOISE
% Pre- cone summation noise (ie Poisson noise)
    neuralResponsePipelineParams.noiseParams.inputConeMosaicNoiseFlag = 'random';

% Post-cone summation noise (Gaussian noise)
neuralResponsePipelineParams.noiseParams.mRGCMosaicNoiseFlag = 'none';

% Post-cone summation noise is additive Gaussian noise with a desired
% sigma. When the input is raw cone excitations, the sigma should be expressed in
% terms of cone excitations/integration time. 
% When the input to mRGCs is cone modulations with respect to the background,
% which have a max amplitude of 1.0, the sigma should be scaled appropriately. 

% Post-cone summation noise when mRGCs are integrating raw cone excitation signals
neuralResponsePipelineParams.noiseParams.mRGCMosaicVMembraneGaussianNoiseSigma = 1e3 * 0.1;

% Post-cone summation noise when mRGCs are integrating  cone excitation modulations
%neuralResponsePipelineParams.noiseParams.mRGCMosaicVMembraneGaussianNoiseSigma = 0.015;

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
%
% Since we are using a contrast-modulation based input to the mRGCmosaic,
% the Poisson noise is a totally correct assumption. But it illustrates
% the inefficiency of the luminance direction.
classifierChoice = 'idealObserver';

switch (classifierChoice) 
    case 'idealObserver'
        % PoissonTAFC makes decision by performing the Poisson likelihood ratio test
        % Also set up parameters associated with use of this classifier.
        theClassifierEngine = responseClassifierEngine(@rcePoisson);
               
        % Train classifier using 1 noise-free instance, 
        % Test performance using a set of 512 noisy instances
        if (fastParameters)
            nTest = 128;
        else
            nTest = 512;
        end
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
            nTrain = 1024*4;
            nTest = 1024*4;
        end

        % Instantiate a computational observer consisting of a linear SVM 
        % coupled with a PCA operating on 2 components
        theClassifierEngine = responseClassifierEngine(@rcePcaSVMTAFC, ...
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
thresholdParams = struct('logThreshLimitLow', 1.5, ...
                       'logThreshLimitHigh', 0, ...
                       'logThreshLimitDelta', 0.01, ...
                       'slopeRangeLow', 1, ...
                       'slopeRangeHigh', 50, ...
                       'slopeDelta', 1.0);

% Parameter for running the QUEST+
% See t_thresholdEngine.m for more on options of the two different mode of
% operation (fixed numer of trials vs. adaptive)

% Sample the contrast-response psychometric curve at this number of
% contrast levels.
%
% Might want to up the number for the non-fastParameters case.
if (fastParameters)
    contrastLevelsSampled = 5;
else
    constrastLevelsSampled = 5;
end

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
dummyFOVdegs = 2.0;
dummyChromaDir = theChromaticDirections(:,1);
theGratingSceneEngine = createGratingScene(dummyChromaDir, theStimulusSpatialFrequencyCPD , ...
        'spatialEnvelope', 'rect', ...
        'fovDegs', dummyFOVdegs, ...
        'minPixelsNumPerCycle', 5, ...
        'pixelsNum', 512);

% Some low res quest params to run fast
questEngineParamsDummy = struct(...
    'minTrial', 16, ...
    'maxTrial', 16, ...
    'numEstimator', 1, ...
    'stopCriterion', 0.5);

% Run the dummy TAFC, just to generate theNeuralEngine and theGratingSceneEngine
computeThreshold(theGratingSceneEngine, theNeuralEngine, theClassifierEngine, ...
        classifierParams, thresholdParams, questEngineParamsDummy,'TAFC',true);

% Having ran the computeThresholdTAFC() function, theNeuralEngine has been generated,
% so we can retrieve from it the size of the inputConeMosaic, and therefore match
% the stimulus spatial params to it as follows.
% The stimulus spatial envelope radius = 1/2 the inputConeMosaic size
theStimulusSpatialEnvelopeRadiusDegs = 0.5*max(theNeuralEngine.neuralPipeline.mRGCMosaic.inputConeMosaic.sizeDegs);

% We make the stimulusFOV 25% larger than the input cone mosaic size
% to avoid OI artifacts that arise when the test image has a mean radiance that is
% different than the mean radiance of the null stimulus, and which can
% result in the test stimulus being discriminable from the null stimulus
% just because of differences in the value with which the OI is padded at
% the edges
theStimulusFOVdegs = max(theNeuralEngine.neuralPipeline.mRGCMosaic.inputConeMosaic.sizeDegs)*1.25;

% This will run faster if you reduce theStimulusPixelsNum to something
% smaller, but then you will get an annoying limit about the precision
% of the cone aperture blurring.
if (fastParameters)
    theStimulusPixelsNum = 512;
    minPixelsNumPerCycle = 4;
else
    theStimulusPixelsNum = 512;
    minPixelsNumPerCycle = 12;
end

% Neural response engine updating
switch (neuralResponsePipelineParams.mRGCMosaicParams.inputSignalType)
    case 'cone_modulations'
        % To compute cone modulations, theNeuralEngine must be provided with theNullStimulusScene
        nullContrast = 0.0;
        theGratingSceneEngine = createGratingScene(dummyChromaDir, theStimulusSpatialFrequencyCPD, ...
            'spatialEnvelope', 'rect', ...
            'orientation', theStimulusOrientationDegs, ...
            'fovDegs', theStimulusFOVdegs, ...
            'spatialEnvelopeRadiusDegs', theStimulusSpatialEnvelopeRadiusDegs, ...
            'spatialPhase', theStimulusSpatialPhaseDegs, ...
            'minPixelsNumPerCycle', minPixelsNumPerCycle, ...
            'pixelsNum', theStimulusPixelsNum);
        theNullStimulusSceneSequence = theGratingSceneEngine.compute(nullContrast);
    
        % Save theNullStimulusScene in the neuralResponsePipelineParams, so
        % that the neural engine can use it to compute mRGCmosaic responses
        % operating on cone contrast (modulation) responses, instead of operating
        % on raw cone excitation responses
        neuralResponsePipelineParams.theNullStimulusScene = theNullStimulusSceneSequence{1};
    
        % Update theNeuralEngine with the new neuralResponsePipelineParams
        theNeuralEngine.updateParamsStruct(neuralResponsePipelineParams);

    case 'cone_excitations'
        % No updating of theNeuralEngine's neuralResponsePipelineParams is needed

    otherwise
        error('Unknown input signal type specified: ''%s''.', neuralResponsePipelineParams.mRGCMosaicParams.inputSignalType)

end

% Generate Matlab filename for saving computed data.
% 
% Code that actually does the save commented out at the end.matFileName = sprintf('mRGCMosaicIsothresholdContour_eccDegs_%2.1f_%2.1f_SpatialFrequencyCPD_%2.1f_OrientationDegs_%d_inputSignal_%s_coneMosaicNoise_%s_mRGCMosaicNoise_%s_vMembraneSigma_%2.3f.mat', ...
matFileName = sprintf('mRGCMosaicIsothresholdContour_eccDegs_%2.1f_%2.1f_SpatialFrequencyCPD_%2.1f_OrientationDegs_%d_inputSignal_%s_coneMosaicNoise_%s_mRGCMosaicNoise_%s_vMembraneSigma_%2.3f.mat', ...
    neuralResponsePipelineParams.mRGCMosaicParams.eccDegs(1), ...
    neuralResponsePipelineParams.mRGCMosaicParams.eccDegs(2), ...
    theStimulusSpatialFrequencyCPD, ...
    theStimulusOrientationDegs, ...
    regexprep(neuralResponsePipelineParams.mRGCMosaicParams.inputSignalType, '_+(\w)', '${upper($1)}'), ...
    neuralResponsePipelineParams.noiseParams.inputConeMosaicNoiseFlag, ...
    neuralResponsePipelineParams.noiseParams.mRGCMosaicNoiseFlag, ...
    neuralResponsePipelineParams.noiseParams.mRGCMosaicVMembraneGaussianNoiseSigma);

%% Ready to compute threshold for each chromatic direction
logThreshold = zeros(1, nChromaticDirections );
theComputedQuestObjects = cell(1,nChromaticDirections);
thePsychometricFunctions = cell(1,nChromaticDirections);
theFittedPsychometricParams = cell(1,nChromaticDirections);
theStimulusScenes = cell(1, nChromaticDirections);

dataFig = figure();
plotRows = 4;
plotCols = 8;
for iChromaDir = 1:nChromaticDirections
    % Create a static grating scene with a particular chromatic direction,
    % spatial frequency, orientation, FOV, and size
    theGratingSceneEngine = createGratingScene(theChromaticDirections(:,iChromaDir), theStimulusSpatialFrequencyCPD, ...
        'spatialEnvelope', 'rect', ...
        'orientation', theStimulusOrientationDegs, ...
        'fovDegs', theStimulusFOVdegs, ...
        'spatialEnvelopeRadiusDegs', theStimulusSpatialEnvelopeRadiusDegs, ...
        'spatialPhase', theStimulusSpatialPhaseDegs, ...
        'minPixelsNumPerCycle', minPixelsNumPerCycle, ...
        'pixelsNum', theStimulusPixelsNum);
   
    % Compute the threshold for our grating scene with the previously
    % defined neural and classifier engine.
    [logThreshold(iChromaDir), questObj, psychometricFunction, fittedPsychometricParams] = ...
        computeThreshold(theGratingSceneEngine, theNeuralEngine, theClassifierEngine, ...
        classifierParams, thresholdParams, questEngineParams, 'TAFC', true);

    % Plot stimulus
    figure(dataFig);
    subplot(plotRows, plotCols, iChromaDir * 2 - 1);

    visualizationContrast = 1.0;
    [theSceneSequence] = theGratingSceneEngine.compute(visualizationContrast);
    theGratingSceneEngine.visualizeStaticFrame(theSceneSequence);

    % Plot data and psychometric curve 
    % with a marker size of 2.5
    figure(dataFig);
    subplot(plotRows, plotCols, iChromaDir * 2);
    questObj.plotMLE(2.5);
    drawnow;

    % Save data for off-line visualizations
    theComputedQuestObjects{iChromaDir} = questObj;
    thePsychometricFunctions{iChromaDir} = psychometricFunction;
    theFittedPsychometricParams{iChromaDir}  = fittedPsychometricParams;
    theStimulusScenes{iChromaDir} = theSceneSequence(1);
end

set(dataFig, 'Position',  [0, 0, 800, 800]);

% Convert returned log threshold to linear threshold
threshold = 10 .^ logThreshold;

% Threshold cone contrasts
thresholdConeContrasts = [threshold.*theChromaticDirections(1,:) ; threshold.*theChromaticDirections(2,:) ; threshold.*theChromaticDirections(3,:)];

%% Fit an ellipse to the data.  See EllipseTest and FitEllipseQ.
%
% The use of scaleFactor to scale up the data and scale down the fit by the
% same amount is fmincon black magic.  Doing this puts the objective
% function into a better range for the default size of search steps.
%
% We constrain the ellipse to line up with the x and y axes.  Change flag
% below to relax this.  Doesn't make very much difference inthis case.
scaleFactor = 1/(max(abs(thresholdConeContrasts(:))));
[fitEllParams,fitA,fitAinv,fitQ] = FitEllipseQ(scaleFactor*thresholdConeContrasts(1:2,:),'lockAngleAt0',true);
nThetaEllipse = 200;
circleIn2D = UnitCircleGenerate(nThetaEllipse);
fitEllipse = PointsOnEllipseQ(fitQ,circleIn2D)/scaleFactor;

%% Plot the fitted ellipse
contrastLim = 0.06;
theContourFig = figure; clf; hold on
plot(thresholdConeContrasts(1,:), thresholdConeContrasts(2,:), 'ok', 'MarkerFaceColor','k', 'MarkerSize',12);
plot(fitEllipse(1,:),fitEllipse(2,:),'r','LineWidth',3);
plot([-contrastLim contrastLim],[0 0],'k:','LineWidth',1);
plot([0 0],[-contrastLim contrastLim],'k:','LineWidth',1);
xlabel('L Cone Contrast');
ylabel('M Cone Contrsast');
set(theContourFig, 'Position',  [800, 0, 600, 800]);
xlim([-contrastLim contrastLim]); ylim([-contrastLim contrastLim]);
axis('square');

%% Export computed data. 
%
% Commented out. If you want to save data like this from a tutorial, put it
% into a local directory and make sure that is gitignored so it doesn't
% clog up the repository.
%
% fprintf('Results will be saved in %s.\n', matFileName);
% save(matFileName, 'theChromaticDirections', 'threshold', ...
%     'theStimulusSpatialFrequencyCPD', 'theStimulusSpatialPhaseDegs', 'theStimulusOrientationDegs', ...
%     'theStimulusFOVdegs', 'theStimulusSpatialEnvelopeRadiusDegs', 'theStimulusScenes',...
%     'theNeuralComputePipelineFunction', 'neuralResponsePipelineParams', ...
%     'classifierChoice', 'classifierParams', 'thresholdParams', ...
%     'theComputedQuestObjects', 'thePsychometricFunctions', 'theFittedPsychometricParams');

end









