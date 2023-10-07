% Compute spatial CSF in different color directions
%
% Description:
%    Use ISETBioCSFGenerator to run out CSFs in different color directions.
%    This example uses an ideal Poisson observer and circularly
%    windowed gratings of constant size.
%
% See also: t_thresholdEngine, t_modulatedGratingsSceneGeneration,
%           t_chromaticThresholdContour, computeThreshold, computePerformance
%

% History:
%   10/20/20  lqz   Wrote it.
%   10/21/20  dhb   More commments.
%   10/22/20  lqz   Restructure the code
%   10/23/20  dhb   More commments.
%   10/25/20  dhb   Change contrast vectors to column vectors.  This is PTB
%                   convention, and also the convention of my brain.
%   05/10/3   fh    Edited it to call the new functions computeThreshold.m
%                       & computePerformance.m & rcePossion.m

% Clear and close
clear; close all;

% List of spatial frequencies to be tested.
spatialFreqs = 2.^(2:5); %[4, 8, 16, 32]

% Choose stimulus chromatic direction specified as a 1-by-3 vector
% of L, M, S cone contrast.  These vectors get normalized below, so only
% their direction matters in the specification.
stimType = 'red-green';
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
% further for this tutorial.  A vector length contrast of 0.1 should be OK.
rmsContrast = 0.1;
chromaDir = chromaDir / norm(chromaDir) * rmsContrast;
assert(abs(norm(chromaDir) - rmsContrast) <= 1e-10);

%% Create neural response engine
%
% This calculations isomerizations in a patch of cone mosaic with Poisson
% noise, and includes optical blur.
% neuralParams = nrePhotopigmentExcitationsConeMosaicHexWithNoEyeMovements;
neuralParams = nrePhotopigmentExcitationsCmosaic;
neuralParams.coneMosaicParams.fovDegs = 0.25;
neuralParams.coneMosaicParams.timeIntegrationSeconds  = 0.1;
% theNeuralEngine = neuralResponseEngine(@nrePhotopigmentExcitationsConeMosaicHexWithNoEyeMovements, neuralParams);
theNeuralEngine = neuralResponseEngine(@nrePhotopigmentExcitationsCmosaic, neuralParams);

%% Instantiate the Poisson responseClassifierEngine
%
% rcePoisson makes decision by performing the Poisson likelihood ratio test
% Also set up parameters associated with use of this classifier.
useOldWay = false;
if useOldWay
    classifierEngine = responseClassifierEngine(@rcePoissonTAFC);
else
    classifierEngine = responseClassifierEngine(@rcePoisson);
end
classifierPara = struct('trainFlag', 'none', ...
                        'testFlag', 'random', ...
                        'nTrain', 1, 'nTest', 128);

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
% See t_thresholdEngine.m for more on options of the two different mode of
% operation (fixed numer of trials vs. adaptive)
thresholdPara = struct('logThreshLimitLow', 2.4, ...
                       'logThreshLimitHigh', 0.0, ...
                       'logThreshLimitDelta', 0.02, ...
                       'slopeRangeLow', 1/20, ...
                       'slopeRangeHigh', 50/20, ...
                       'slopeDelta', 2.5/20, ...
                       'thresholdCriterion', 0.81606);

questEnginePara = struct( ...
    'qpPF',@qpPFWeibullLog, ...
    'minTrial', 1280, ...
    'maxTrial', 1280, ...
    'numEstimator', 1, ...
    'stopCriterion', 0.05);

%% Compute threshold for each spatial frequency
% 
% See toolbox/helpers for functions createGratingScene computeThreshold
dataFig = figure();
logThreshold = zeros(1, length(spatialFreqs));
for idx = 1:length(spatialFreqs)
    % Create a static grating scene with a particular chromatic direction,
    % spatial frequency, and temporal duration
    gratingScene = createGratingScene(chromaDir, spatialFreqs(idx),...
        'fovDegs', neuralParams.coneMosaicParams.fovDegs,'duration',...
            neuralParams.coneMosaicParams.timeIntegrationSeconds);
    
    % Compute the threshold for our grating scene with the previously
    % defined neural and classifier engine.  This function does a lot of
    % work, see t_tresholdEngine and the function itself, as well as
    % function computePerformance.
    if useOldWay
        [logThreshold(idx), questObj, ~, para(idx,:)] = ...
            computeThresholdTAFC(gratingScene, theNeuralEngine, classifierEngine, ...
            classifierPara, thresholdPara, questEnginePara);
    else
        [logThreshold(idx), questObj, ~, para(idx,:)] = ...
            computeThreshold(gratingScene, theNeuralEngine, classifierEngine, ...
            classifierPara, thresholdPara, questEnginePara, 'TAFC', true);
    end
    
    % Plot stimulus
    figure(dataFig);
    subplot(4, 2, idx * 2 - 1);
    
    visualizationContrast = 1.0;
    [theSceneSequence] = gratingScene.compute(visualizationContrast);
    gratingScene.visualizeStaticFrame(theSceneSequence);
    
    % Plot data and psychometric curve 
    % with a marker size of 2.5
    subplot(4, 2, idx * 2);
    questObj.plotMLE(2.5,'para',para(idx,:));
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
