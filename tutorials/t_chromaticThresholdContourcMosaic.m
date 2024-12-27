% Compute isothreshold contour in different color directions
%
% Description:
%    Use ISETBioCSFGenerator to run out an isothreshold contour in the LM
%    contrast plane. This example uses an ideal Poisson TAFC observer and circularly
%    windowed gratings of constant size and one spatial frequency.
%
% See also: t_spatialCSFcMosaic, t_modulatedGratingsSceneGeneration,
%           computeThreshold, computePerformance,
%           createGratingScene.
%

% History:
%   10/23/20  dhb   Wrote it from what is now called t_spatialCSFcMosaic

% Clear and close
clear; close all;

% Spatial frequencies to be tested, and other parameters.
% Beyond spatial frequency, these override default parameters
% in createGratingScene via key/value pair.
%
% Using 90 degree (sine phase) makes the stimulus symmetric in terms of
% balanced incremental and decremental components, so that the isothreshold
% contour is also symmetric.
spatialFreq = 4;
gratingPhaseDeg = 90;
gratingFovDegs = 0.3;

% Set up a set of chromatic directions. Passing elevation = 90 puts these
% in the LM contrast plan.  These are at constant rms (vector length)
% contrast.
%
% Things may go badly if you exceed the gamut of the monitor, so we are
% conservative and set this at a value that is within gamut of typical
% monitors and don't worry about it further for this tutorial.  A vector
% length contrast of 0.08 should be OK.
rmsContrast = 0.1;
nDirs = 8;
for ii = 1:nDirs
    theta = (ii-1)/nDirs*2*pi;
    theDirs(:,ii) = [cos(theta) sin(theta) 0]';
    theDirs(:,ii) = theDirs(:,ii) / norm(theDirs(:,ii)) * rmsContrast;
    assert(abs(norm(theDirs(:,ii)) - rmsContrast) <= 1e-10);
end

%% Instantiate the Poisson responseClassifierEngine
%
% PoissonTAFC makes decision by performing the Poisson likelihood ratio test
% Also set up parameters associated with use of this classifier.
classifierEngine = responseClassifierEngine(@rcePoisson);
classifierPara = struct('trainFlag', 'none', ...
                        'testFlag', 'random', ...
                        'nTrain', 1, 'nTest', 128);

%% Parameters for threshold estimation/quest engine
% The actual threshold varies enough with the different engines that we
% need to adjust the contrast range that Quest+ searches over, as well as
% the range of psychometric function slopes. Threshold limits are computed
% as 10^-logThreshLimitVal.
thresholdPara = struct('logThreshLimitLow', 2.4, ...
                       'logThreshLimitHigh', 0.0, ...
                       'logThreshLimitDelta', 0.02, ...
                       'slopeRangeLow', 1, ...
                       'slopeRangeHigh', 50, ...
                       'slopeDelta', 2.5);

% Parameter for running the QUEST+
% See t_spatialCSF.m for more on options of the two different mode of
% operation (fixed numer of trials vs. adaptive)
questEnginePara = struct('minTrial', 1280*2, 'maxTrial', 1280*2, ...
                         'numEstimator', 1, 'stopCriterion', 0.05);

%% Compute threshold for each chromatic direction
% 
% See toolbox/helpers for functions createGratingScene, computeThreshold
dataFig = figure();
logThreshold = zeros(1, nDirs);
for ii = 1:nDirs
    % Create a grating scene engine with a particular chromatic direction,
    % spatial frequency, and temporal duration
    theSceneEngine = createGratingScene(theDirs(:,ii), spatialFreq, 'spatialPhase', gratingPhaseDeg, 'fovDegs', gratingFovDegs);

    %% Create neural response engine on first pass.
    % 
    % Do this after creating the first scene so we can match up integration
    % time with frame duration.
    if (ii == 1)
        % This calculations isomerizations in a patch of cone mosaic with Poisson
        % noise, and includes optical blur.
        noiseFreeResponseParams = nreNoiseFreePhotopigmentExcitationsCMosaic;
        noiseFreeResponseParams.coneMosaicParams.sizeDegs = [0.25 0.25];
        noiseFreeResponseParams.coneMosaicParams.timeIntegrationSeconds = 0.1;
        noisyInstancesParams = nreNoisyInstancesPoisson;
        theNeuralEngine = neuralResponseEngine( ...
            @nreNoiseFreePhotopigmentExcitationsCMosaic, ...
            @nreNoisyInstancesPoisson, ...
            noiseFreeResponseParams, ...
            noisyInstancesParams);
    end
    
    % Compute the threshold for our grating scene with the previously
    % defined neural and classifier engine.  This function does a lot of
    % work, see t_tresholdEngine and the function itself, as well as
    % function computePerformanceTAFC.
    [logThreshold(ii), questObj] = ...
        computeThreshold(theSceneEngine, theNeuralEngine, classifierEngine, classifierPara, ...
        thresholdPara, questEnginePara,'TAFC',true);
    
    % Plot stimulus
    figure(dataFig);
    subplot(ceil(nDirs/2), 4, ii * 2 - 1);
    
    visualizationContrast = 1.0;
    [theSceneSequence] = theSceneEngine.compute(visualizationContrast);
    theSceneEngine.visualizeStaticFrame(theSceneSequence);
    
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

% Fit an ellipse to the data.  See EllipseTest and FitEllipseQ.
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

% Plot
contrastLim = 0.02;
theContourFig = figure; clf; hold on
plot(thresholdConeContrasts(1,:), thresholdConeContrasts(2,:), 'ok', 'MarkerFaceColor','k', 'MarkerSize',12);
plot(fitEllipse(1,:),fitEllipse(2,:),'r','LineWidth',3);
plot([-contrastLim contrastLim],[0 0],'k:','LineWidth',1);
plot([0 0],[-contrastLim contrastLim],'k:','LineWidth',1);
xlabel('L Cone Contrast');
ylabel('M Cone Contrsast');
set(theContourFig, 'Position',  [800, 0, 600, 800]);
% xlim([-contrastLim contrastLim]); ylim([-contrastLim contrastLim]);
%axis('square'); 
axis equal

