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
nAList = length(nAlternativesList);

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
neuralParams = nrePhotopigmentExcitationsCmosaicWithNoEyeMovements;
neuralParams.coneMosaicParams.fovDegs = 0.25;
neuralParams.coneMosaicParams.timeIntegrationSeconds = 0.1;
theNeuralEngine = neuralResponseEngine(@nrePhotopigmentExcitationsCmosaicWithNoEyeMovements, neuralParams);

%% Instantiate the PoissonTAFC responseClassifierEngine
%
% PoissonTAFC makes decision by performing the Poisson likelihood ratio test
% Also set up parameters associated with use of this classifier.
classifierEngine = responseClassifierEngineNWay(@rcePoissonNWay_OneStimulusPerTrial);
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
[thresholdPara, questEnginePara] = deal(cell(1, nAList));
%loop through each NWay alternative choiced choice because the number of
%alternatives affect the guess rate as well as the threshold criterion.
for t = 1:nAList
    thresholdPara{t} = struct('logThreshLimitLow', 2.4, ...
        'logThreshLimitHigh', 0.0, ...
        'logThreshLimitDelta', 0.05, ...
        'slopeRangeLow', 1/20, ...
        'slopeRangeHigh', 60/20, ...
        'slopeDelta', 5/20, ...
        'thresholdCriterion', (1- 1/nAlternativesList(t))/2 + 1/nAlternativesList(t), ...
        'guessRate', 1/nAlternativesList(t), ...
        'lapseRate', 1e-4);
    
    estDomain = -thresholdPara{t}.logThreshLimitLow:...
        thresholdPara{t}.logThreshLimitDelta: -thresholdPara{t}.logThreshLimitHigh;
    slopeRange = thresholdPara{t}.slopeRangeLow: ...
        thresholdPara{t}.slopeDelta: thresholdPara{t}.slopeRangeHigh; 
    
    questEnginePara{t} = struct( ...
        'employMethodOfConstantStimuli', true,...
        'nTest', 10,...
        'validation', true,...
        'blocked',true,...
        'estDomain',estDomain,...
        'slopeRange',slopeRange,...
        'qpPF',@qpPFWeibullLog, ...
        'numEstimator', 1, ...
        'stopCriterion', 0.05);
end

%% Compute threshold for each number of alternatives specified
% 
% See toolbox/helpers for functions createGratingScene computeThresholdTAFC
logThreshold = zeros(1, nAList);
dataFig = figure();
for idx = 1:nAList
    % Open figure

    % Set nAlternatives
    nAlternatives = nAlternativesList(idx);

    % Create grating scenes with a particular chromatic direction for each
    % alternative. 
    % spatial frequency, and temporal duration
    orientations = linspace(0,(nAlternatives-1)*180/nAlternatives,nAlternatives);
    for oo = 1:length(orientations)
        gratingScenes{oo} = createGratingScene(chromaDir, spatialFreq,...
            'orientation',orientations(oo), 'fovDegs', neuralParams.coneMosaicParams.fovDegs);
    end
    
    % Compute the threshold for our grating scene with the previously
    % defined neural and classifier engine.  This function does a lot of
    % work, see t_tresholdEngine and the function itself, as well as
    % function computePerformanceNWay_OneStimulusPerTrial.
    [logThreshold(idx), questObj, ~, para(idx,:)] = ...
        computeThresholdNWay_OneStimulusPerTrial(gratingScenes,  ...
        theNeuralEngine, classifierEngine, ...
        classifierPara, thresholdPara{idx}, questEnginePara{idx});
    
    % Plot stimulus
    figure(idx)
    visualizationContrast = 1.0;
    for t = 1:nAlternativesList(idx)
        subplot(8, 4, t);
        theSceneSequence = gratingScenes{t}.compute(visualizationContrast);
        gratingScenes{t}.visualizeStaticFrame(theSceneSequence);
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
