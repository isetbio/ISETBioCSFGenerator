function thresholdRet = t_spatialCSFEyeMovements(varargin)
% Compute spatial CSF in different color directions, with fEM
%
% Syntax:
%   thresholdRet = t_spatialCSF;
%
% Description:
%    Use ISETBioCSFGenerator to run out CSFs in different color directions.
%    This example uses an ideal Poisson observer and circularly
%    windowed gratings of constant size, with fixational eye movements.
%
%
% See also: t_spatialCSF, t_modulatedGratingsSceneGeneration,
%           t_chromaticThresholdContour, computeThreshold, computePerformance

% History:
%   10/20/20  lqz   Wrote it.
%   10/21/20  dhb   More commments.
%   10/22/20  lqz   Restructure the code
%   10/23/20  dhb   More commments.
%   10/25/20  dhb   Change contrast vectors to column vectors.  This is PTB
%                   convention, and also the convention of my brain.
%   05/10/23  fh    Edited it to call the new functions computeThreshold.m
%                       & computePerformance.m & rcePossion.m
%   04/17/24  dhb   Remove oldWay option.  Ever forward.  Enforce sine
%                   phase.
%   12/19/24  dhb   Update for new architecture.

% Close figs
close all

% Parse input
p = inputParser;
p.addParameter('filter', struct('spectralSupport',[],'transmission',[]), @isstruct);
p.addParameter('doValidationCheck', true, @islogical);
parse(p, varargin{:});
filter = p.Results.filter;
doValidationCheck = p.Results.doValidationCheck;

% Freeze rng for replicatbility
rng(1);

% List of spatial frequencies to be tested.
spatialFreqs = [4, 8, 16, 32];
if (length(spatialFreqs) ~= 4 | ~all(spatialFreqs == [4, 8, 16, 32]))
    doValidationCheck = false;
end

% Choose stimulus chromatic direction specified as a 1-by-3 vector
% of L, M, S cone contrast.  These vectors get normalized below, so only
% their direction matters in the specification.
stimType = 'luminance';
switch (stimType)
    case 'luminance'
        chromaDir = [1.0, 1.0, 1.0]';
    case 'red-green'
        chromaDir = [1.0, -1.0, 0.0]';
        doValidationCheck = false;
    case 'L-isolating'
        chromaDir = [1.0, 0.0, 0.0]';
        doValidationCheck = false;
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
noiseFreeResponseParams = nreNoiseFreePhotopigmentExcitationsCmosaic;
noiseFreeResponseParams.coneMosaicParams.sizeDegs = [0.5 0.5];
noiseFreeResponseParams.coneMosaicParams.timeIntegrationSeconds = 0.1;

% Stimulus timing parameters
frameDurationSeconds = noiseFreeResponseParams.coneMosaicParams.timeIntegrationSeconds;
framesNum = 4;
stimulusDuration = framesNum*frameDurationSeconds;

noisyInstancesParams = nreNoisyInstancesPoisson;
theNeuralEngine = neuralResponseEngine( ...
    @nreNoiseFreePhotopigmentExcitationsCmosaic, ...
    @nreNoisyInstancesPoisson, ...
    noiseFreeResponseParams, ...
    noisyInstancesParams);

if (~all(noiseFreeResponseParams.coneMosaicParams.sizeDegs == [0.5 0.5]))
    doValidationCheck = false;
end
if (noiseFreeResponseParams.coneMosaicParams.timeIntegrationSeconds ~= 0.1)
    doValidationCheck = false;
end

%% Instantiate the responseClassifierEngine
%
% rcePoisson makes decision by performing the Poisson likelihood ratio
% test. This is the ideal observer for the Poisson noice cone excitations
% illustrated in this script.  But you can run other rce's as well.
%    rcePoisson - signal known exactly Poission max likelihood
%    rceTemplateDistance - signal known exactly nearest L2 template
%                 distance.
%    rcePcaSVM  - support vector machine linear classifier after PCA.
%
% Also set up parameters associated with use of this classifier.
classifierEngine = 'rcePoisson';
switch (classifierEngine)
    case {'rcePoisson'}
        classifierEngine = responseClassifierEngine(@rcePoisson);
        classifierUsePara = struct('trainFlag', 'none', ...
            'testFlag', 'random', ...
            'nTrain', 1, 'nTest', 128);
    case {'rceTemplateDistance'}
        classifierEngine = responseClassifierEngine(@rcePoisson);
        classifierUsePara = struct('trainFlag', 'none', ...
            'testFlag', 'random', ...
            'nTrain', 1, 'nTest', 128);
        doValidationCheck = false;
    case 'rcePcaSVM'
        pcaSVMParams = struct(...
            'PCAComponentsNum', 2, ...          % number of PCs used for feature set dimensionality reduction
            'crossValidationFoldsNum', 10, ...  % employ a 10-fold cross-validated linear
            'kernelFunction', 'linear', ...     % linear
            'classifierType', 'svm' ...         % binary SVM classifier
            );
        classifierEngine = responseClassifierEngine(@rcePcaSVM,pcaSVMParams);
        classifierUsePara = struct('trainFlag', 'none', ...
            'testFlag', 'random', ...
            'nTrain', 128, 'nTest', 128);  
        doValidationCheck = false;
    otherwise
        error('Unsupported rce specified')
end

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
% See t_spatialCSF.m for more on options of the two different mode of
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
% See toolbox/helpers for functions createGratingScene computeThreshold
dataFig = figure();
logThreshold = zeros(1, length(spatialFreqs));
for idx = 1:length(spatialFreqs)
    % Create a static grating scene with a particular chromatic direction,
    % spatial frequency, and temporal duration.  Put grating in sine phase
    % becuase that keeps the spatial mean constant across spatial
    % frequencies.
    %
    % Create scene produces square scenes.  We use the min of the mosaic
    % field size to pick a reasonable size
    gratingScene = createGratingScene(chromaDir, spatialFreqs(idx),...
        'fovDegs', min(noiseFreeResponseParams.coneMosaicParams.sizeDegs), ...
        'presentationMode', 'flashedmultiframe', ...
        'duration', stimulusDuration, ...
        'temporalFrequency', framesNum/stimulusDuration ...
        );

    % Compute the threshold for our grating scene with the previously
    % defined neural and classifier engine.  This function does a lot of
    % work, see t_tresholdEngine and the function itself, as well as
    % function computePerformance.
    [logThreshold(idx), questObj, ~, para(idx,:)] = ...
        computeThreshold(gratingScene, theNeuralEngine, classifierEngine, ...
        classifierUsePara, thresholdPara, questEnginePara, ...
        'TAFC', true, ...
        'trainFixationalEM', [], ...
        'testFixationalEM', []);

    % Plot stimulus
    figure(dataFig);
    subplot(length(spatialFreqs), 2, idx * 2 - 1);

    visualizationContrast = 1.0;
    [theSceneSequence] = gratingScene.compute(visualizationContrast);
    gratingScene.visualizeStaticFrame(theSceneSequence);

    % Plot data and psychometric curve
    % with a marker size of 2.5
    subplot(length(spatialFreqs), 2, idx * 2);
    questObj.plotMLE(2.5,'para',para(idx,:));
    drawnow;
end
set(dataFig, 'Position',  [0, 0, 800, 800]);

% Convert returned log threshold to linear threshold
threshold = 10 .^ logThreshold;

%% Plot Contrast Sensitivity Function
theCsfFig = figure();
loglog(spatialFreqs, 1 ./ threshold, '-ok', 'LineWidth', 2);
xticks(spatialFreqs); xlim([spatialFreqs(1), spatialFreqs(end)]);
yticks([2,5,10,20,50]); ylim([1, 50]);
xlabel('Spatial Frequency (cyc/deg)');
ylabel('Sensitivity');
set(theCsfFig, 'Position',  [800, 0, 600, 800]);

%% Do a check on the answer
%
% So that if we break something in the future we will have
% a chance of knowing it. The current numbers don't quite
% match the old version, but I think that is because of a change
% away from iePoisson which was not freezing the rng, and also
% other changes somewhere in stochasticity that I have not quite
% tracked down. But this validation generally passes.  Might fail
% sometimes.
if (doValidationCheck)
    % validationThresholds = [0.0418    0.0783    0.1540    0.6759];
    % if (any(abs(threshold-validationThresholds)./validationThresholds > 0.25))
    %     error('Do not replicate validation thresholds to 25%. Check that parameters match, or for a bug.');
    % end
end

%% Return a value if it was requested
if (nargout > 0)
    thresholdRet = threshold;
end
end
