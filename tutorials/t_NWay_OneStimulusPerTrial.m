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
% See also: t_spatialCSF, t_modulatedGratingsSceneGeneration,
%           t_spatialCSF, t_chromaticThresholdContour, computeThreshold, computePerformance
%

% History:
%    12/07/21  dhb  Wrote it from t_spatialCSF.
%    05/10/23  fh   Edited it to call the new functions computeThreshold.m
%                       & computePerformance.m & rcePossion.m
%    04/17/24  dhb  Deep six oldWay option. Ever forward.  Gratings in sine
%                   phase.
%    12/20/24  dhb  Update for new architecture

%% Clear and close
clear; close all;

% Freeze rng
rng(0);

%% Set number of alternatives
nAlternativesList = [2 4 8];
nAList = length(nAlternativesList);

% Use meta contrast or not
useMetaContrast = true;

% Use fixational eye movements or not
useFixationalEMs = false;

% Run just one spatial frequency.
spatialFreq = 15;

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
 noiseFreeResponseParams = nreNoiseFreeCMosaic;
    noiseFreeResponseParams.coneMosaicParams.sizeDegs = [0.25 0.25];
    noiseFreeResponseParams.coneMosaicParams.timeIntegrationSeconds = 0.1;
    noisyInstancesParams = nreNoisyInstancesPoisson;
    theNeuralEngine = neuralResponseEngine( ...
        @nreNoiseFreeCMosaic, ...
        @nreNoisyInstancesPoisson, ...
        noiseFreeResponseParams, ...
        noisyInstancesParams);

%% Instantiate the PoissonTAFC responseClassifierEngine
%
% PoissonTAFC makes decision by performing the Poisson likelihood ratio test
% Also set up parameters associated with use of this classifier.
classifierPara = struct('trainFlag', 'none', ...
                        'testFlag', 'random', ...
                        'nTrain', 1, 'nTest', 128);
classifierEngine = responseClassifierEngine(@rcePoisson, classifierPara);

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
% computeThreshold.
thresholdPara = struct('logThreshLimitLow', 2.4, ...
    'logThreshLimitHigh', 0.0, ...
    'logThreshLimitDelta', 0.01, ...
    'slopeRangeLow', 1/20, ...
    'slopeRangeHigh', 60/20, ...
    'slopeDelta', 3/20, ...
    'thresholdCriterion', 0.75, ... %will be respecified in the loop below
    'guessRate', 0.5, ... %will be respecified in the loop below
    'lapseRate', 1e-4);

% Define grids for the threshold and slope
estDomain = -thresholdPara.logThreshLimitLow:...
    thresholdPara.logThreshLimitDelta: -thresholdPara.logThreshLimitHigh;
slopeRange = thresholdPara.slopeRangeLow: ...
    thresholdPara.slopeDelta: thresholdPara.slopeRangeHigh; 

questEnginePara = struct( ...
    'employMethodOfConstantStimuli', true,...
    'nTest', 10,...
    'validation', true,...
    'blocked',true,...
    'estDomain',estDomain,...
    'slopeRange',slopeRange,...
    'qpPF',@qpPFWeibullLog, ...
    'numEstimator', 1, ...
    'stopCriterion', 0.05);

%% Set up EM object to define fixational eye movement paths
%
% There are lots of things we can control. In any case, the path
% defined here is used for training.
%
% Setting the seed to -1 means don't touch the seed or call rng().
if (useFixationalEMs)
    trainFixationalEMObj = fixationalEM;
    trainFixationalEMObj.microSaccadeType = 'none';
    trainFixationalEMObj.randomSeed = -1;
    femDuration = stimulusDuration;
    femTimeStep = frameDurationSeconds;
    femNumberFEMs = 1;
    femComputeVelocity = false;
    trainFixationalEMObj.compute(femDuration, femTimeStep, femNumberFEMs, femComputeVelocity);

    % And this one for testing.
    testFixationalEMObj = trainFixationalEMObj;
else
    trainFixationalEMObj = [];
    testFixationalEMObj = [];
end

%% Compute threshold for each number of alternatives specified
% 
% See toolbox/helpers for functions createGratingSceneEngine computeThreshold
logThreshold = zeros(1, nAList);
para         = NaN(nAList, 4); %4 paramters (lapse rate and guess rate are fixed)

for idx = 1:nAList
    % Create grating scenes with a particular chromatic direction for each
    % alternative, spatial frequency, and temporal duration.  Use sine
    % phase gratings to keep mean constant even if some of the grating is
    % truncated given scene size.
    orientations = linspace(0,(nAlternativesList(idx)-1)*...
        180/nAlternativesList(idx),nAlternativesList(idx));
    gratingSceneEngines = cell(1,nAlternativesList(idx));
    for t = 1:nAlternativesList(idx)
        gratingSceneEngines{t} = createGratingSceneEngine(chromaDir, spatialFreq,...
            'orientation',orientations(t), ...
            'fovDegs', min(noiseFreeResponseParams.coneMosaicParams.sizeDegs), ...
            'duration', noiseFreeResponseParams.coneMosaicParams.timeIntegrationSeconds, ...
            'spatialPhase', 90);
    end

    % Respecify the threshold and the guess rate
    thresholdPara.thresholdCriterion = (1- 1/nAlternativesList(idx))/2 + 1/nAlternativesList(idx);
    thresholdPara.guessRate = 1/nAlternativesList(idx);

    % If using meta contrast, set this up.
    if (useMetaContrast)
        metaSceneEngineParams = sceMetaContrast;
        theMetaSceneEngine = sceneEngine(@sceMetaContrast,metaSceneEngineParams);

        % Create nreMetaContrast using the actual scene and neural engines
        for t = 1:nAlternativesList(idx)
            metaNeuralResponseEngineNoiseFreeParams = nreNoiseFreeMetaContrast;
            metaNeuralResponseEngineNoiseFreeParams.contrast0 = 0;
            metaNeuralResponseEngineNoiseFreeParams.contrast1 = 1;
            metaNeuralResponseEngineNoiseFreeParams.neuralEngine = theNeuralEngine;

            metaNeuralResponseEngineNoisyInstanceParams =  nreNoisyInstancesMetaContrast;
            metaNeuralResponseEngineNoisyInstanceParams.neuralEngine = theNeuralEngine;

            metaNeuralResponseEngineNoiseFreeParams.sceneEngine = gratingSceneEngines{t};
            theMetaNeuralEngines{t} = neuralResponseEngine(@nreNoiseFreeMetaContrast, ...
                @nreNoisyInstancesMetaContrast, ...
                metaNeuralResponseEngineNoiseFreeParams, ...
                metaNeuralResponseEngineNoisyInstanceParams);
        end

        % Compute the threshold for our grating scene with meta scene and
        % and neural response engines. This function does a lot of
        % work, well as function computePerformance.
        [logThreshold(idx), questObj, ~, para(idx,:)] = ...
            computeThreshold(theMetaSceneEngine, theMetaNeuralEngines, classifierEngine, ...
            classifierPara, thresholdPara, questEnginePara, ...
            'trainFixationalEM', trainFixationalEMObj, ...
            'testFixationalEM', testFixationalEMObj, ...
            'TAFC', false, 'useMetaContrast', true);
    else
        % Compute the threshold for our grating scene with the previously
        % defined neural and classifier engine.  This function does a lot of
        % work, as well as function computePerformance.
        [logThreshold(idx), questObj, ~, para(idx,:)] = computeThreshold(...
            gratingSceneEngines, theNeuralEngine, classifierEngine,...
            classifierPara, thresholdPara, questEnginePara, ...
            'trainFixationalEM', trainFixationalEMObj, ...
            'testFixationalEM', testFixationalEMObj, ...
            'TAFC',false, 'useMetaContrast', false);
    end
    
    % Plot stimulus
    figure(idx)
    visualizationContrast = 1.0;
    for t = 1:nAlternativesList(idx)
        subplot(8, 4, t);
        theSceneSequence = gratingSceneEngines{t}.compute(visualizationContrast);
        gratingSceneEngines{t}.visualizeStaticFrame(theSceneSequence);
    end
    
    % Plot data and psychometric curve 
    % with a marker size of 2.5
    subplot(8,4,17:32);
    questObj.plotMLE(2.5,'para',para(idx,:));
    drawnow;
    set(gcf, 'Units','normalized','Position',  [0, 0, 0.15, 0.5]);
end

% Convert returned log threshold to linear threshold
threshold = 10 .^ logThreshold;

%% Plot sensitivity against number of alternatives
theNwayFig = figure;
loglog(nAlternativesList, 1 ./ threshold, '-ok', 'LineWidth', 2);
xlabel('Number of Alternatives');
ylabel('Sensitivity');
set(theNwayFig, 'Position',  [800, 0, 600, 800]);

%% These are the thresholds obtained in December 2024
%
% Both with and without meta contrast
%   threshold = 0.0365    0.0503    0.0633

