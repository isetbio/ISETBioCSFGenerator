function thresholdRet = t_spatialCSFcMosaic(options)
% Compute spatial CSF in different color directions, with fEM
%
% Syntax:
%   thresholdRet = t_spatialCSFcMosaic;
%
% Description:
%    Use ISETBioCSFGenerator to run out CSFs in different color directions.
%    This example uses an ideal Poisson observer and circularly
%    windowed gratings of constant size, with fixational eye movements.
%
%    This is set up with key/value pairs that demonstate how to mix and
%    match noise free neural response models, noisy instance models, and
%    classifier engines. Different choices are illustrated in the examples
%    in the source code.
%
% Optional key/value pairs
%    See source code arguments block for a list of key/value pairs.

% See also: t_spatialCSFcMosaic, t_modulatedGratingsSceneGeneration,
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
%   04/17/24  dhb   Remove oldWay option.  Ever forward.  Enforce sine phase.
%   12/19/24  dhb   Update for new architecture.

% Examples:
%{
     t_spatialCSFcMosaic('useMetaContrast', false, ...
        'useFixationalEMs', false, ...
        'whichNoiseFreeNre', 'excitationsCmosaic', ...
        'whichNoisyInstanceNre', 'Poisson', ...
        'whichClassifierEngine', 'rcePoisson');

    t_spatialCSFcMosaic('useMetaContrast', true, ...
        'useFixationalEMs', false, ...
        'whichNoiseFreeNre', 'sceneAsResponses', ...
        'whichNoisyInstanceNre', 'Poisson', ...
        'whichClassifierEngine', 'rcePoisson');

    t_spatialCSFcMosaic('useMetaContrast', true, ...
        'useFixationalEMs', true, ...
        'whichNoiseFreeNre', 'excitationsCmosaic', ...
        'whichNoisyInstanceNre', 'Poisson', ...
        'whichClassifierEngine', 'rcePoisson');

     t_spatialCSFcMosaic('useMetaContrast', true, ...
        'useFixationalEMs', false, ...
        'whichNoiseFreeNre', 'excitationsCmosaic', ...
        'whichNoisyInstanceNre', 'Gaussian', ...
        'whichClassifierEngine', 'rceTemplateDistance');

    t_spatialCSFcMosaic('useMetaContrast', true, ...
        'useFixationalEMs', true, ...
        'whichNoiseFreeNre', 'excitationsCmosaic', ...
        'whichNoisyInstanceNre', 'Poisson', ...
        'whichClassifierEngine', 'rcePcaSVM');

%}

arguments
    % Run the validation check?  This gets overridden to false if other
    % options change the conditions so that the validation data don't
    % apply.
    options.doValidationCheck (1,1) logical = true;

    % Apply a filter to the spectra before computing responses?  See
    % t_spatialCSFcMosaicFilter
    options.filter (1,1) = struct('spectralSupport',[],'transmission',[]);

    % Use meta contrast method to speed things up?
    options.useMetaContrast (1,1) logical = true;

    % Use fixational eye movements?
    options.useFixationalEMs (1,1) logical = false;

    % Choose noise free neural model
    %   Choices: 'excitationsCmosaic'
    %            'sceneAsResponses'
    options.whichNoiseFreeNre (1,:) char  = 'excitationsCmosaic'

    % Choose noise model
    %   Choices: 'Poisson'
    %            'Gaussian'
    options.whichNoisyInstanceNre (1,:) char = 'Poisson'

    % Choose classifier engine
    %    rcePoisson - signal known exactly Poission max likelihood
    %    rceTemplateDistance - signal known exactly nearest L2 template
    %                 distance.
    %    rcePcaSVM  - support vector machine linear classifier after PCA.
    options.whichClassifierEngine (1,:) char = 'rcePoisson'
end

%% Close any stray figs
close all;

%% Set flags from key/value pairs
filter = options.filter;
doValidationCheck = options.doValidationCheck;
useMetaContrast = options.useMetaContrast;
useFixationalEMs = options.useFixationalEMs;
whichNoiseFreeNre = options.whichNoiseFreeNre;
whichNoisyInstanceNre = options.whichNoisyInstanceNre;
whichClassifierEngine = options.whichClassifierEngine;


%% Set some parameters that control what we do
if (~useFixationalEMs)
    framesNum = 1;
else
    framesNum = 4;
    doValidationCheck = false;
end
sizeDegs = [0.5 0.5];
frameDurationSeconds = 0.1;

%% Freeze rng for replicatbility
rng(1);

%% List of spatial frequencies to be tested.
spatialFreqs = [4, 8, 16, 32];
if (length(spatialFreqs) ~= 4 | ~all(spatialFreqs == [4, 8, 16, 32]))
    doValidationCheck = false;
end

%% Chromatic direction and contrast
%
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
switch (whichNoiseFreeNre)
    case 'excitationsCmosaic'
        nreNoiseFreeResponse = @nreNoiseFreePhotopigmentExcitationsCMosaic;
        noiseFreeResponseParams = nreNoiseFreePhotopigmentExcitationsCMosaic;
        noiseFreeResponseParams.coneMosaicParams.sizeDegs = sizeDegs;
        noiseFreeResponseParams.coneMosaicParams.timeIntegrationSeconds = frameDurationSeconds;
        if (~all(noiseFreeResponseParams.coneMosaicParams.sizeDegs == [0.5 0.5]))
            doValidationCheck = false;
        end
        if (noiseFreeResponseParams.coneMosaicParams.timeIntegrationSeconds ~= 0.1)
            doValidationCheck = false;
        end

        % These are tuned to match range of nre performance. See Quest
        % structure below for comments on what they mean.
        logThresholdLimitLow = 2.4;
        logThresholdLimitHigh = 0;
        logThresholdLimitDelta = 0.02;
        slopeRangeLow = 1/20;
        slopeRangeHigh = 50/20;
        slopeDelta = 2.5/20;

    case 'sceneAsResponses'
        % This is image photon counts, and thus provides an estimate of the
        % upper bound on performance for a photon limited system.
        nreNoiseFreeResponse = @nreNoiseFreeSceneAsResponses;
        noiseFreeResponseParams = nreNoiseFreeSceneAsResponses;
        doValidationCheck = false;

        % These are tuned to match range of nre performance. See Quest
        % structure below for comments on what they mean.
        logThresholdLimitLow = 12;
        logThresholdLimitHigh = 6;
        logThresholdLimitDelta = 0.05;
        slopeRangeLow = 1e-2/20;
        slopeRangeHigh = 1/20;
        slopeDelta = 1e-1/20;

    otherwise
        error('Unsupported noise free nre specified');
end

switch (whichNoisyInstanceNre)
    case 'Poisson'
        nreNoisyInstances = @nreNoisyInstancesPoisson;
        noisyInstancesParams = nreNoisyInstancesPoisson;

    case 'Gaussian'
        nreNoisyInstances = @nreNoisyInstancesGaussian;
        noisyInstancesParams = nreNoisyInstancesGaussian;
        doValidationCheck = false;

    otherwise
        error('Unsupported noisy instances nre specified');
end

% Stimulus timing parametersf
stimulusDuration = framesNum*frameDurationSeconds;

theNeuralEngine = neuralResponseEngine( ...
    nreNoiseFreeResponse, ...
    nreNoisyInstances, ...
    noiseFreeResponseParams, ...
    noisyInstancesParams);

%% Instantiate the responseClassifierEngine
%
% rcePoisson makes decision by performing the Poisson likelihood ratio
% test. This is the ideal observer for the Poisson noice cone excitations
% illustrated in this script.  But you can run other rce's as well.
%
% Set up parameters associated with use of this classifier.
switch (whichClassifierEngine)
    case {'rcePoisson'}
        whichClassifierEngine = responseClassifierEngine(@rcePoisson);
        classifierPara = struct('trainFlag', 'none', ...
            'testFlag', 'random', ...
            'nTrain', 1, 'nTest', 128);

    case {'rceTemplateDistance'}
        whichClassifierEngine = responseClassifierEngine(@rcePoisson);
        classifierPara = struct('trainFlag', 'none', ...
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
        whichClassifierEngine = responseClassifierEngine(@rcePcaSVM,pcaSVMParams);
        classifierPara = struct('trainFlag', 'none', ...
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
thresholdPara = struct('logThreshLimitLow', logThresholdLimitLow, ...
    'logThreshLimitHigh', logThresholdLimitHigh, ...
    'logThreshLimitDelta', logThresholdLimitDelta, ...
    'slopeRangeLow', slopeRangeLow, ...
    'slopeRangeHigh', slopeRangeHigh, ...
    'slopeDelta', slopeDelta, ...
    'thresholdCriterion', 0.81606);

questEnginePara = struct( ...
    'qpPF',@qpPFWeibullLog, ...
    'minTrial', 1280, ...
    'maxTrial', 1280, ...
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

%% If using meta contrast, set this up.
if (useMetaContrast)
    metaSceneEngineParams = sceMetaContrast;
    theMetaSceneEngine = sceneEngine(@sceMetaContrast,metaSceneEngineParams);

    % Create nreMetaContrast using the actual scene and neural engines
    metaNeuralResponseEngineNoiseFreeParams = nreNoiseFreeMetaContrast;
    metaNeuralResponseEngineNoiseFreeParams.contrast0 = 0;
    metaNeuralResponseEngineNoiseFreeParams.contrast1 = 1;
    metaNeuralResponseEngineNoiseFreeParams.neuralEngine = theNeuralEngine;

    metaNeuralResponseEngineNoisyInstanceParams =  nreNoisyInstancesMetaContrast;
    metaNeuralResponseEngineNoisyInstanceParams.neuralEngine = theNeuralEngine;
end

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
    gratingSceneEngine = createGratingScene(chromaDir, spatialFreqs(idx),...
        'fovDegs', min(sizeDegs), ...
        'presentationMode', 'flashedmultiframe', ...
        'duration', stimulusDuration, ...
        'temporalFrequency', framesNum/stimulusDuration, ...
        'spatialPhase', 90, ...
        'filter', filter...
        );

    % Create the sceMetaContrast scene engine
    if (useMetaContrast)
        % Instantiate meta contrast neural engine for this spatial
        % frequency and use it to compute threshold
        metaNeuralResponseEngineNoiseFreeParams.sceneEngine = gratingSceneEngine;
        theMetaNeuralEngine = neuralResponseEngine(@nreNoiseFreeMetaContrast, ...
            @nreNoisyInstancesMetaContrast, ...
            metaNeuralResponseEngineNoiseFreeParams, ...
            metaNeuralResponseEngineNoisyInstanceParams);

        % Compute the threshold for our grating scene with meta scene and
        % and neural response engines. This function does a lot of work,
        % see the function itself, as well as function computePerformance.
        [logThreshold(idx), questObj, ~, para(idx,:)] = ...
            computeThreshold(theMetaSceneEngine, theMetaNeuralEngine, whichClassifierEngine, ...
            classifierPara, thresholdPara, questEnginePara, ...
            'TAFC', true, 'useMetaContrast', true, ...
            'trainFixationalEM', trainFixationalEMObj, ...
            'testFixationalEM', testFixationalEMObj);
    else
        % Compute the threshold for our grating scene with the previously
        % defined neural and classifier engine.  This function does a lot
        % of work, see the function itself, as well as function
        % computePerformance.
        [logThreshold(idx), questObj, ~, para(idx,:)] = ...
            computeThreshold(gratingSceneEngine, theNeuralEngine, whichClassifierEngine, ...
            classifierPara, thresholdPara, questEnginePara, ...
            'TAFC', true, 'useMetaContrast', false, ...
            'trainFixationalEM', trainFixationalEMObj, ...
            'testFixationalEM', testFixationalEMObj);
    end

    % Plot stimulus
    figure(dataFig);
    subplot(length(spatialFreqs), 2, idx * 2 - 1);

    visualizationContrast = 1.0;
    [theSceneSequence] = gratingSceneEngine.compute(visualizationContrast);
    gratingSceneEngine.visualizeStaticFrame(theSceneSequence);

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
    validationThresholds = [0.0418    0.0783    0.1540    0.6759];
    if (any(abs(threshold-validationThresholds)./validationThresholds > 0.25))
        error('Do not replicate validation thresholds to 25%. Check that parameters match, or for a bug.');
    else
        fprintf('Validation regression check passes\n');
    end
else
    fprintf('Validation regression check not run, presumably because of parameter change\n');

end

%% Return a value if it was requested
if (nargout > 0)
    thresholdRet = threshold;
end

end
