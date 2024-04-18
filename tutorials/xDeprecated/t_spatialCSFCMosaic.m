% Compute spatial CSF in different color directions
%
% Description:
%    Use ISETBioCSFGenerator to run out CSFs in different color directions.
%    This example uses the new @cMosaic object, an ideal Poisson TAFC observer,
%    and square windowed gratings of constant size.
%
% See also: t_thresholdEngine, t_modulatedGratingsSceneGeneration,
%           t_chromaticThresholdContour, computeThresholdTAFC, computePerformanceTAFC
%

% History:
%   3/29/21  npc   Wrote it for @cMosaic by adapting t_spatialCSF


% Clear and close
clear; close all;

% List of spatial frequencies to be tested.
spatialFreqs = [2 4 8 12 16 24];

% Choose stimulus chromatic direction specified as a 1-by-3 vector
% of L, M, S cone contrast.  These vectors get normalized below, so only
% their direction matters in the specification.
stimType = 'luminance';
switch (stimType)
    case 'luminance'
        chromaDir = [1.0, 1.0, 1.0]';
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

%% Create neural response engine
%
% This calculates isomerizations in a patch of cone mosaic with Poisson
% noise, and includes optical blur.
theNeuralComputePipelineFunction = @nrePhotopigmentExcitationsCmosaic;
neuralParams = theNeuralComputePipelineFunction();

% Make the mosaic large enough so that it covers 1 cycle of the lowest
% spatial frequecy examined, plus a little extra
headroomFactor = 1.25;
neuralParams.coneMosaicParams.sizeDegs = headroomFactor*1.0/(min(spatialFreqs))*[1 1];
theNeuralEngine = neuralResponseEngine(theNeuralComputePipelineFunction, neuralParams);

%% Instantiate the PoissonTAFC responseClassifierEngine
%
% PoissonTAFC makes decision by performing the Poisson likelihood ratio test
% Also set up parameters associated with use of this classifier.
classifierEngine = responseClassifierEngine(@rcePoissonTAFC);
nTrain = 1;
nTest = 512;
classifierPara = struct('trainFlag', 'none', ...
                        'testFlag', 'random', ...
                        'nTrain', nTrain, ...
                        'nTest', nTest);

%% Parameters for threshold estimation/quest engine
% The actual threshold varies enough with the different engines that we
% need to adjust the contrast range that Quest+ searches over, as well as
% the range of psychometric function slopes. Threshold limits are computed
% as 10^-logThreshLimitVal.
thresholdPara = struct('logThreshLimitLow', 2.4, ...
                       'logThreshLimitHigh', 0.0, ...
                       'logThreshLimitDelta', 0.01, ...
                       'slopeRangeLow', 1, ...
                       'slopeRangeHigh', 50, ...
                       'slopeDelta', 2.0);

% Parameter for running the QUEST+
% See t_thresholdEngine.m for more on options of the two different mode of
% operation (fixed numer of trials vs. adaptive)
contrastLevelsSampled = 10;
questEnginePara = struct(...
    'minTrial', contrastLevelsSampled*nTest, ...
    'maxTrial', contrastLevelsSampled*nTest, ...
    'numEstimator', 1, ...
    'stopCriterion', 0.05);

%% Compute threshold for each spatial frequency
% 
% See toolbox/helpers for functions createGratingScene computeThresholdTAFC
dataFig = figure();
logThreshold = zeros(1, length(spatialFreqs));
for idx = 1:length(spatialFreqs)
    % Create a static grating scene with a particular chromatic direction,
    % spatial frequency, and temporal duration.  Using a spatial phase
    % of 90 degrees makes the stimulus symmetric in the scene, and prevents
    % there from being a DC component that varies as a function of spatial
    % frequency.  We almost always want to use sine phase for this sort of
    % calculation, but somehow the defaults are almost always cosine phase.
    gratingScene = createGratingScene(chromaDir, spatialFreqs(idx), ...
        'spatialEnvelope', 'rect', ...
        'fovDegs', max(neuralParams.coneMosaicParams.fovDegs)*1.02, ...
        'minPixelsNumPerCycle', 5, ...
        'spatialPhase',90);
    
    % Compute the threshold for our grating scene with the previously
    % defined neural and classifier engine.  This function does a lot of
    % work, see t_tresholdEngine and the function itself, as well as
    % function computePerformanceTAFC.
    [logThreshold(idx), questObj] = ...
        computeThresholdTAFC(gratingScene, theNeuralEngine, classifierEngine, ...
        classifierPara, thresholdPara, questEnginePara);
    % [logThreshold(idx), questObj] = ...
    %     computeThreshold(gratingScene, theNeuralEngine, classifierEngine, ...
    %     classifierPara, thresholdPara, questEnginePara);
    
    % Plot stimulus
    figure(dataFig);
    subplot(4, 4, idx * 2 - 1);
    
    visualizationContrast = 1.0;
    [theSceneSequence] = gratingScene.compute(visualizationContrast);
    gratingScene.visualizeStaticFrame(theSceneSequence);
    
    % Plot data and psychometric curve 
    % with a marker size of 2.5
    subplot(4, 4, idx * 2);
    questObj.plotMLE(2.5);
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
