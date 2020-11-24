% Compute isothreshold contour in different color directions using a midget RGC mosaic
%
% Description:
%    Use ISETBioCSFGenerator to run out an isothreshold contour in the LM
%    contrast plane. This example uses an ideal Poisson TAFC observer and circularly
%    windowed gratings of constant size and one spatial frequency.
%
% See also: t_spatialCSF, t_thresholdEngine, t_modulatedGratingsSceneGeneration,
%           computeThresholdTAFC, computePerformanceTAFC,
%           createGratingScene.
%

% History:
%   11/04/20  NPC   Wrote it from t_chromaticThresholdContour

% Clear and close
clear; close all;

% Spatial frequencies to be tested, and other parameters.
% Beyond spatial frequency, these override default parameters
% in createGratingScene via key/value pair.
%

% 100 msec stimulus duration
stimulusDurationSeconds = 100/1000;

% Options for presentationMode are {'sampled motion', 'flashed'}
presentationMode = 'flashed';
    
% For 'flashed' presentation mode, present the grating at 90 spatial phase 
% (odd symmetry).  Using 90 degree (sine phase) makes the stimulus symmetric 
% in terms of balanced incremental and decremental components, so that the 
% isothreshold contour is also symmetric.
gratingPhaseDeg = 90;

% Grating spatial frequency in c/deg
spatialFreq = 1.0;

% Set up a set of chromatic directions. These are at constant rms (vector length)
% contrast.
%
% Things may go badly if you exceed the gamut of the monitor, so we are
% conservative and set this at a value that is within gamut of typical
% monitors and don't worry about it further for this tutorial.  A vector
% length contrast of 0.08 should be OK.
rmsContrast = 0.08;
nDirs = 8;
for ii = 1:nDirs
    theta = (ii-1)/nDirs*2*pi;
    theDirs(:,ii) = [cos(theta) sin(theta) 0]';
    theDirs(:,ii) = theDirs(:,ii) / norm(theDirs(:,ii)) * rmsContrast;
    assert(abs(norm(theDirs(:,ii)) - rmsContrast) <= 1e-10);
end

%% Create neural response engine
%
% Obtain default neural response engine params
neuralParams = nreMidgetRGC;

% Modify mRGC mosaic eccentricity and size
neuralParams.mRGCmosaicParams.eccDegs = [1 0];
neuralParams.mRGCmosaicParams.sizeDegs = 0.5*[1 1];

% *** POST-CONE SUMMATION NOISE ***
% Set the mRGC mosaic (post-cone summation) noise flag. If set to 'none',
% the only noise in the computation is that of the coneMosaic. 
%
% If set to 'random', Gaussian noise is added at the final mRGC response.
% The noise sd is noiseFactor*maxResponse
neuralParams.coneMosaicParams.noiseFlag = 'none';
neuralParams.mRGCmosaicParams.noiseFlag = 'random';
neuralParams.mRGCmosaicParams.noiseFactor = 0.25;

% Modify some cone mosaic params
neuralParams.coneMosaicParams.coneMosaicResamplingFactor = 3;
% 25 msec integration time to speed up computation
neuralParams.coneMosaicParams.integrationTime = 25/1000;

% Instantiate the neural response engine
theNeuralEngine = neuralResponseEngine(@nreMidgetRGC, neuralParams);

%% Instantiate the PoissonTAFC or the PcaSVMTAFC responseClassifierEngine
%
% Options:
%   'idealObserver' - ideal TAFC observer for Poisson limited signals
%   'computationalObsever' - SVM based learned classifier
classifierChoice = 'idealObserver';
switch (classifierChoice) 
    case 'idealObserver'
        % PoissonTAFC makes decision by performing the Poisson likelihood ratio test
        % Also set up parameters associated with use of this classifier.
        classifierEngine = responseClassifierEngine(@rcePoissonTAFC);
        % Train classifier using 1 noise-free instance, 
        % Test performance using a set of 128 noisy instances
        classifierPara = struct('trainFlag', 'none', ...
                                'testFlag', 'random', ...
                                'nTrain', 1, 'nTest', 128);
    case 'computationalObserver'
        % Instantiate a computational observer consisting of a linear SVM 
        % coupled with a PCA operating on 2 components
        classifierEngine = responseClassifierEngine(@rcePcaSVMTAFC, ...
            struct(...
                'PCAComponentsNum', 2, ...          % number of PCs used for feature set dimensionality reduction
                'crossValidationFoldsNum', 10, ...  % employ a 10-fold cross-validated linear 
                'kernelFunction', 'linear', ...     % linear
                'classifierType', 'svm' ...         % binary SVM classifier
            ));
        % Train classifier using a set of 256 noisy instances, 
        % Test performance using a set of 128 noisy instances
        classifierPara = struct('trainFlag', 'random', ...
                                'testFlag', 'random', ...
                                'nTrain', 256, 'nTest', 128);
                        
    otherwise
        error('Unknown classifier: ''%s''.', classifierChoice);
end


%% Parameters for threshold estimation/quest engine
% The actual threshold varies enough with the different engines that we
% need to adjust the contrast range that Quest+ searches over, as well as
% the range of psychometric function slopes.
thresholdPara = struct('logThreshLimitLow', 2.4, ...
                       'logThreshLimitHigh', 0.0, ...
                       'logThreshLimitDelta', 0.02, ...
                       'slopeRangeLow', 1, ...
                       'slopeRangeHigh', 50, ...
                       'slopeDelta', 2.5);

% Parameter for running the QUEST+
% See t_thresholdEngine.m for more on options of the two different mode of
% operation (fixed numer of trials vs. adaptive)
questEnginePara = struct('minTrial', 1280, 'maxTrial', 1280, ...
                         'numEstimator', 1, 'stopCriterion', 0.05);

% Visualization params
visualizationPara.visualizeStimulus = true;

% Data saving params
datasavePara.destDir = '~/Desktop/tmpDir';
datasavePara.saveMRGCResponses = true;

%% Compute threshold for each spatial direction
% 
% See toolbox/helpers for functions createGratingScene computeThresholdTAFC
dataFig = figure();
logThreshold = zeros(1, nDirs);

for ii = 1:nDirs
    % Create a static grating scene with a particular chromatic direction,
    % spatial frequency, and temporal duration. Make it twice as large as
    % the mRGC mosaic so that it extends over cone inputs to  the surround
    % subregions of the RGC cells, which are quite large (~7 times the RF
    % center).
    % Options for presentationMode are {'sampled motion', 'flashed'}
    % For 'flashed' we make the duration equal to the
    gratingScene = createGratingScene(theDirs(:,ii), spatialFreq, ...
        'spatialPhase', gratingPhaseDeg, ...
        'duration', stimulusDurationSeconds, ...
        'fovDegs', max(neuralParams.mRGCmosaicParams.sizeDegs)*2, ...
        'spatialEnvelope', 'square', ...
        'presentationMode', presentationMode ...
        );
    
    % Compute the threshold for our grating scene with the previously
    % defined neural and classifier engine.  This function does a lot of
    % work, see t_tresholdEngine and the function itself, as well as
    % function computePerformanceTAFC.
    [logThreshold(ii), questObj] = ...
        computeThresholdTAFC(gratingScene, theNeuralEngine, classifierEngine, ...
        classifierPara, thresholdPara, questEnginePara, visualizationPara, ...
        datasavePara);
    
    % Plot stimulus
    figure(dataFig);
    subplot(ceil(nDirs/2), 4, ii * 2 - 1);
    
    visualizationContrast = 1.0;
    [theSceneSequence] = gratingScene.compute(visualizationContrast);
    gratingScene.visualizeStaticFrame(theSceneSequence);
    
    % Plot data and psychometric curve 
    % with a marker size of 2.5
    subplot(ceil(nDirs/2), 4, ii * 2);
    questObj.plotMLE(2.5);
end

set(dataFig, 'Position',  [0, 0, 800, 800]);

% Convert returned log threshold to linear threshold
threshold = 10 .^ logThreshold;

% Threshold cone contrasts
thresholdConeContrasts = [threshold.*theDirs(1,:) ; threshold.*theDirs(2,:) ; threshold.*theDirs(3,:)];

% Fit an ellipse to the data.  See EllipseTest and EllipsoidFit.
%
% The use of scaleFactor to scale up the data and scale down the fit by the
% same amount is fmincon black magic.  Doing this puts the objective
% function into a better range for the default size of search steps.
scaleFactor = 10;
fitCenter = zeros(3,1);
[fitA,fitAinv,fitQ,fitEllParams] = EllipsoidFit(scaleFactor*thresholdConeContrasts,[],false,true);
nThetaEllipse = 200;
circleIn2D = UnitCircleGenerate(nThetaEllipse);
circleInLMPlane = [circleIn2D(1,:) ; circleIn2D(2,:) ; zeros(size(circleIn2D(1,:)))];
fitEllipse = PointsOnEllipsoidFind(fitQ,circleInLMPlane,fitCenter)/scaleFactor;

% Plot
contrastLim = 0.08;
figure; clf; hold on
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

