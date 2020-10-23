% Compute isothreshold contour in different color directions
%
% Description:
%    Use ISETBioCSFGenerator to run out an isothreshold contour in the LM
%    contrast plane. This example uses an ideal Poisson TAFC observer and circularly
%    windowed gratings of constant size and one spatial frequency.
%
% See also: t_spatialCSF, t_thresholdEngine, t_modulatedGratingsSceneGeneration,
%           computeThresholdTAFC, computePerformanceTAFC
%

% History:
%   10/23/20  dhb   Wrote it from t_spatialCSF

% Clear and close
clear; close all;

% Spatial frequencies to be tested.
spatialFreq = 2;

% Set up a set of chromatic directions. Passing elevation = 90 puts these
% in the LM contrast plan.  These are at constant rms (vector length)
% contrast.
%
% Things may go badly if you exceed the gamut of the monitor, so we are
% conservative and set this at a value that is within gamut of typical
% monitors and don't worry about it further for this tutorial.  A vector
% length contrast of 0.08 should be OK.
rmsContrast = 0.08;
nDirs = 10;
for ii = 1:nDirs
    theta = (ii-1)/nDirs*2*pi;
    theDirs(:,ii) = [cos(theta) sin(theta) 0]';
    theDirs(:,ii) = theDirs(:,ii) / norm(theDirs(:,ii)) * rmsContrast;
    assert(abs(norm(theDirs(:,ii)) - rmsContrast) <= 1e-10);
end

%% Create neural response engine
%
% This calculations isomerizations in a patch of cone mosaic with Poisson
% noise, and includes optical blur.
neuralParams = nrePhotopigmentExcitationsWithNoEyeMovements;
neuralParams.coneMosaicParams.fovDegs = 0.25;
theNeuralEngine = neuralResponseEngine(@nrePhotopigmentExcitationsWithNoEyeMovements, neuralParams);

%% Instantiate the PoissonTAFC responseClassifierEngine
%
% PoissonTAFC makes decision by performing the Poisson likelihood ratio test
% Also set up parameters associated with use of this classifier.
classifierEngine = responseClassifierEngine(@rcePoissonTAFC);
classifierPara = struct('trainFlag', 'none', ...
                        'testFlag', 'random', ...
                        'nTrain', 1, 'nTest', 128);

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

%% Compute threshold for each spatial direction
% 
% See toolbox/helpers for functions createGratingScene computeThresholdTAFC
dataFig = figure();
logThreshold = zeros(1, nDirs);
for ii = 1:nDirs
    % Create a static grating scene with a particular chromatic direction,
    % spatial frequency, and temporal duration
    gratingScene = createGratingScene(theDirs(:,ii), spatialFreq);
    
    % Compute the threshold for our grating scene with the previously
    % defined neural and classifier engine.  This function does a lot of
    % work, see t_tresholdEngine and the function itself, as well as
    % function computePerformanceTAFC.
    [logThreshold(ii), questObj] = ...
        computeThresholdTAFC(gratingScene, theNeuralEngine, classifierEngine, classifierPara, thresholdPara, questEnginePara);
    
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

%% Plot Contrast Sensitivity Function
theCsfFig = figure();
plot(threshold.*theDirs(1,:), threshold.*theDirs(2,:), '-ok', 'MarkerSize',12, 'LineWidth', 2);
xlabel('L Cone Contrast');
ylabel('M Cone Contrsast');
set(theCsfFig, 'Position',  [800, 0, 600, 800]);
xlim([-0.04 0.04]); ylim([-0.04 0.04]);
axis('square');

