% Compute contrast sensitivity for different numbers of alternatives
%
% Description:
%    Use ISETBioCSFGenerator to run out sensitivity for different N alternatives.
%    We expect sensitivity to go down, both because stimuli become more
%    similar and because uncertainty increases.
%
%    This example uses an ideal Poisson observer and gratings at different
%    orientations.
%
%    This version illustrates how to model an NWay forced choice task, in
%    which one alternative of N alternatives is presented on each trial,
%    and the subject indicates which.
%
% See also: t_thresholdEngine, t_modulatedGratingsSceneGeneration,
%           t_spatialCSF, t_chromaticThresholdContour, computeThresholdTAFC, computePerformanceTAFC
%

% History:
%    12/07/21  dhb  Wrote it from t_spatialCSF.

%% Clear and close
clear; close all;

%% Set number of alternatives
nAlternativesList = [2 4 8 16];

% Run just one spatial frequency.
spatialFreq = 2;

% Set the RMS cone contrast of the stimulus. Things may go badly if you
% exceed the gamut of the monitor, so we are conservative and set this at a
% value that is within gamut of typical monitors and don't worry about it
% further for this tutorial. You would need to use smaller numbers for
% chromatic directions.
chromaDir = [1.0, 1.0, 1.0]';
rmsContrast = 0.9;
chromaDir = chromaDir / norm(chromaDir) * rmsContrast;
assert(abs(norm(chromaDir) - rmsContrast) <= 1e-10);

%% Create neural response engine
%
% This calculations isomerizations in a patch of cone mosaic with Poisson
% noise, and includes optical blur.
neuralParams = nrePhotopigmentExcitationsWithNoEyeMovements;
neuralParams.coneMosaicParams.fovDegs = 0.25;
neuralParams.coneMosaicParams.timeIntegrationSeconds = 0.005;
theNeuralEngine = neuralResponseEngine(@nrePhotopigmentExcitationsWithNoEyeMovements, neuralParams);

%% Instantiate the PoissonTAFC responseClassifierEngine
%
% PoissonTAFC makes decision by performing the Poisson likelihood ratio test
% Also set up parameters associated with use of this classifier.
classifierEngine = responseClassifierEngineNWay(@rcePoissonNWay_OneStimPerTrial);
classifierPara = struct('trainFlag', 'none', ...
                        'testFlag', 'random', ...
                        'nTrain', 1, 'nTest', 128);

%% Parameters for threshold estimation/quest engine
% The actual threshold varies enough with the different engines that we
% need to adjust the contrast range that Quest+ searches over, as well as
% the range of psychometric function slopes. Threshold limits are computed
% as 10^-logThreshLimitVal. The reason it is log units is that below we
% define the PF for the questEngine as @qpPFWeibullLog. Had we used the
% default (@qpPFWeibull), the units would have been dB.
%
% Also note explicit passing of proportion correct criterion for threshold.
% The default value of 0.81606 matches the parameterization of mQUESTPlus'
% Weibull PFs, when lapse rate is 0 and guess rate is 0.5.  But it seems
% better to pass it explicitly so we know what it is.  Using a lower
% criterion for the NAFC example run here seems fun.
%
% There are two separate structures below. The conceptual distinction
% between them is not entirely clear.  These are interpretted by
% computeThresholdNWay_OneStimulusPerTrial.
thresholdPara = struct('logThreshLimitLow', 2.4, ...
    'logThreshLimitHigh', 0.0, ...
    'logThreshLimitDelta', 0.02, ...
    'slopeRangeLow', 1/20, ...
    'slopeRangeHigh', 100/20, ...
    'slopeDelta', 5/20, ...
    'thresholdCriterion', 0.60, ...
    'guessRate', [], ...
    'lapseRate', [0 0.02]);

questEnginePara = struct( ...
    'qpPF',@qpPFWeibullLog, ...
    'minTrial', 1280, ...
    'maxTrial', 1280, ...
    'numEstimator', 1, ...
    'stopCriterion', 0.05);

%% Compute threshold for each number of alternatives specified
% 
% See toolbox/helpers for functions createGratingScene computeThresholdTAFC
logThreshold = zeros(1, length(nAlternativesList));
dataFig = figure();
for idx = 1:length(nAlternativesList)
    % Open figure

    % Set nAlternatives
    nAlternatives = nAlternativesList(idx);
    
    % Set guess rate
    thresholdPara.guessRate = 1/nAlternatives;

    % Create grating scenes with a particular chromatic direction for each
    % alternative. 
    % spatial frequency, and temporal duration
    orientations = linspace(0,(nAlternatives-1)*180/nAlternatives,nAlternatives);
    for oo = 1:length(orientations)
        gratingScenes{oo} = createGratingScene(chromaDir, spatialFreq,'orientation',orientations(oo));
    end
    
    % Compute the threshold for our grating scene with the previously
    % defined neural and classifier engine.  This function does a lot of
    % work, see t_tresholdEngine and the function itself, as well as
    % function computePerformanceNWay_OneStimulusPerTrial.
    [logThreshold(idx), questObj] = ...
        computeThresholdNWay_OneStimulusPerTrial(gratingScenes, theNeuralEngine, classifierEngine, ...
        classifierPara, thresholdPara, questEnginePara);
    
    % Plot stimulus
    figure(dataFig);
    subplot(4, 4, idx * 2 - 1);
    
%     visualizationContrast = 1.0;
%     [theSceneSequence] = gratingScenes.compute(visualizationContrast);
%     gratingScenes.visualizeStaticFrame(theSceneSequence);
    
    % Plot data and psychometric curve 
    % with a marker size of 2.5
    subplot(4, 4, idx * 2);
    questObj.plotMLE(2.5);
    drawnow;
end
set(dataFig, 'Position',  [0, 0, 800, 800]);

% Convert returned log threshold to linear threshold
threshold = 10 .^ logThreshold;

%% Plot sensitivity against number of alternatives
theNwayFig = figure;
loglog(nAlternativesList, 1 ./ threshold, '-ok', 'LineWidth', 2);
xlabel('Number of Alternatives');
ylabel('Sensitivity');
set(theNwayFig, 'Position',  [800, 0, 600, 800]);
