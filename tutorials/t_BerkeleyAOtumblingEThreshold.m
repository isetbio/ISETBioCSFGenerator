function t_BerkeleyAOtumblingEThreshold(options)
% Compute tumbling E threshold with AO optics
%
% This takes a number of key/value pairs that control its detailed
% conditions.  See comments under options block below.
%
% Note, in case you are tempted, that this can't be run with the meta
% contrast method because the parameter is stimulus size rather than
% contrast, and responses do not scale as a linear function of size the way
% they do with contrast.
%
% See also t_BerkeleyAOtumblingESceneEngine.

%% Pick up optional arguments
arguments
    % Add this amount of defocus to the diffraction limited optics
    options.defocusDiopters (1,1) double = 0.05;

    % Pupil diameter in mm.
    options.pupilDiameterMm (1,1) double = 6;

    % Print out more diagnostics, or not
    options.verbose (1,1) logical = false;
end

% Initialize
close all;

% Make sure figures and results directories exist so that output writes
% don't fail
rootPath = ISETBioCSFGeneratorRootPath;
if (~exist(fullfile(rootPath,'local',mfilename,'figures'),'dir'))
    mkdir(fullfile(rootPath,'local',mfilename,'figures'));
end
if (~exist(fullfile(rootPath,'local',mfilename,'results'),'dir'))
    mkdir(fullfile(rootPath,'local',mfilename,'results'));
end

% Define the AO scene parameters for the experiment we are modeling
%
% Get the tumbling E scene engines.
%
% At the moment cannot vary the set of orientations, but would be easy
% enough to allow this with a key/value pair passed to the called tutorial
% function.
%
% The scene engine tutorial returns its parameters, which are used below to
% try to match things up as best as possible.
orientations = [0 90 180 270];
[sce0,sce90,sce180,sce270,backgroundSceneEngine,sceneParams] = t_BerkeleyAOtumblingESceneEngine('VisualizeScene',false);
tumblingEsceneEngines = {sce0, sce90, sce180, sce270};
clear sce0 sce90 sce180 sce270

% Parameters. These control many aspects of what gets done, particular the subject.
%
% To run a different subject or pupil size, change 'psdDataSubdir'
% field below to have the desired subject number in the directory string, and the desired
% pupil string in the name if using other than the default 4 mm.
%
% Parameter fields below allow change of integration time, lens age,
% MPD, L:M:S proportion (although you won't get many S because of the
% tritanopic zone).
%
% Note that that the custom pupil diameter field does not affect the
% optical quality because we are using the PSF read from the PSF data
% file.  What it does is allow you to set a pupil diameter different
% from the one for which the PSF was computed.  What this does is allow
% independent control of retinal illuminance and PSF, if you want to
% separate the two effects.  The full data set does contain PSFs
% computed for different pupil sizes for each TCA/LCA combination which
% may be used to explore the effect of pupil size on optical quality.
% Note again that using a PSF computed for various pupil sizes will not
% affect the retinal illuminance as that is controlled by the pupil
% size set explicitly here.  The PSF file naming convention for the
% various pupil sizes is described in the README.
%
% The code in this script is reasonably clever about creating figure
% and results subdirectories to hold its ouput, that keep the separate
% conditions you might run separate.  But it may not be perfect at
% this, particularly if you start digging deeper into the code and
% customizing more things.
params = struct(...
    'letterSizesNumExamined',  9, ...                           % How many sizes to use for sampling the psychometric curve (9 used in the paper)
    'maxLetterSizeDegs', 0.2, ...                               % The maximum letter size in degrees of visual angle
    'sceneUpSampleFactor', 4, ...                               % Upsample scene, so that the pixel for the smallest scene is < cone aperture
    'mosaicIntegrationTimeSeconds', 1/backgroundSceneEngine.sceneParams.temporalModulationParams.frameRateHz, ... % Integration time, here one scene frame
    'nTest', 512, ...                                           % Number of trial to use for computing Pcorrect
    'thresholdP', 0.781, ...                                    % Probability correct level for estimating threshold performance
    'customLensAgeYears', [], ...                               % Lens age in years (valid range: 20-80), or empty to use the default age of 32.
    'customMacularPigmentDensity', [], ...                      % Cu∂stom MPD, or empty to use the default density of 0.35; example, 0.7
    'customConeDensities', [], ...                              % Custom L-M-S ratio or empty to use default; example [0.6 0.3 0.1]
    'customPupilDiameterMM', [], ...                            % Custom pupil diameter in MM or empty to use the value from the psfDataFile
    'visualizedPSFwavelengths', [], ...                         % Vector with wavelengths for visualizing the PSF. If set to empty[] there is no visualization; example 400:20:700
    'visualizeDisplayCharacteristics', ~true, ...               % Flag, indicating whether to visualize the display characteristics
    'visualizeScene', ~true, ...                                % Flag, indicating whether to visualize one of the scenes
    'visualEsOnMosaic', ~true ...                               % Flag, indicating whether to visualize E's against mosaic as function of their size
    );

% % Set up summary filename and output dir
% summaryFileName = sprintf('Summary_%s_%dms.mat', strrep(params.psfDataSubDir, '.mat', ''), round(1000*params.mosaicIntegrationTimeSeconds));
% if (~isempty(params.customMacularPigmentDensity))
%     summaryFileName = strrep(summaryFileName, '.mat', sprintf('_MPD_%2.2f.mat', params.customMacularPigmentDensity));
% end
% if (~isempty(params.customPupilDiameterMM))
%     summaryFileName = strrep(summaryFileName, '.mat', sprintf('_pupilDiamMM_%2.2f.mat', params.customPupilDiameterMM));
% end
% if (~isempty(params.customConeDensities))
%     summaryFileName = strrep(summaryFileName, '.mat', sprintf('_cones_%2.2f_%2.2f_%2.2f.mat', params.customConeDensities(1), params.customConeDensities(2), params.customConeDensities(3)));
% end
% if (~isempty(params.customLensAgeYears))
%     summaryFileName = strrep(summaryFileName, '.mat', sprintf('_lensAge_%d.mat', params.customLensAgeYears));
% end
% params.outputResultsDir = fullfile(ISETBioCSFGeneratorrootPath,'local',mfilename,'results',strrep(summaryFileName, '.mat',''));
% params.outputFiguresDir =  fullfile(ISETBioCSFGeneratorrootPath,'local',mfilename,'figures',strrep(summaryFileName, '.mat',''));
% if (~exist(params.outputResultsDir,'dir'))
%     mkdir(params.outputResultsDir);
% end
% if (~exist(params.outputFiguresDir,'dir'))
%     mkdir(params.outputFiguresDir);
% end

% Unpack simulation params
letterSizesNumExamined = params.letterSizesNumExamined;
maxLetterSizeDegs = params.maxLetterSizeDegs;
mosaicIntegrationTimeSeconds = params.mosaicIntegrationTimeSeconds;
nTest = params.nTest;
thresholdP = params.thresholdP;

%% Create neural response engine
%
% This calculates excitations in a patch of cone mosaic with Poisson
% noise, and includes optical blur.
noiseFreeResponseParams = nreNoiseFreeCMosaic([],[],[],[],'opticsType','BerkeleyAO');

% Set optics params
wls = sceneParams.spectralSupport;
fieldSizeDegs = sceneParams.displayParams.displayFOVDeg;
accommodatedWl = sceneParams.displayParams.AOPrimaryWls(1);
pupilDiameterMm = options.pupilDiameterMm;
defocusDiopters = options.defocusDiopters;

noiseFreeResponseParams.opticsParams.wls = wls;
noiseFreeResponseParams.opticsParams.opticsType = 'BerkeleyAO';
noiseFreeResponseParams.opticsParams.pupilDiameterMM = pupilDiameterMm;
noiseFreeResponseParams.opticsParams.defocusAmount = defocusDiopters;
noiseFreeResponseParams.opticsParams.accommodatedWl = accommodatedWl;
noiseFreeResponseParams.opticsParams.zCoeffs = zeros(66,1);
noiseFreeResponseParams.opticsParams.defeatLCA = true;
noiseFreeResponseParams.verbose = options.verbose;

% Cone params
noiseFreeResponseParams.coneMosaicParams.wave = wls;
noiseFreeResponseParams.coneMosaicParams.fovDegs = fieldSizeDegs;

noiseFreeResponseParams = nreNoiseFreeCMosaic;
noiseFreeResponseParams.coneMosaicParams.sizeDegs = [0.5 0.5];
noiseFreeResponseParams.coneMosaicParams.timeIntegrationSeconds = mosaicIntegrationTimeSeconds;

noisyInstancesParams = nreNoisyInstancesPoisson;
theNeuralEngine = neuralResponseEngine( ...
    @nreNoiseFreeCMosaic, ...
    @nreNoisyInstancesPoisson, ...
    noiseFreeResponseParams, ...
    noisyInstancesParams);

% Poisson N-way classifier
classifierEngine = responseClassifierEngine(@rcePoisson);
classifierPara = struct('trainFlag', 'none', ...
    'testFlag', 'random', ...
    'nTrain', 1, 'nTest', nTest);

%% Parameters for threshold estimation/quest engine
thresholdPara = struct(...
    'maxParamValue', maxLetterSizeDegs, ...    % The maximum value of the examined param (letter size in degs)
    'logThreshLimitLow', 2.0, ...              % minimum log10(normalized param value)
    'logThreshLimitHigh', 0.0, ...             % maximum log10(normalized param value)
    'logThreshLimitDelta', 0.01, ...
    'slopeRangeLow', 1/20, ...
    'slopeRangeHigh', 500/20, ...
    'slopeDelta', 2/20, ...
    'thresholdCriterion', thresholdP, ...
    'guessRate', 1/numel(orientations), ...
    'lapseRate', [0 0.02]);

% Parameters for Quest
questEnginePara = struct( ...
    'qpPF',@qpPFWeibullLog, ...
    'minTrial', nTest*letterSizesNumExamined, ...
    'maxTrial', nTest*letterSizesNumExamined, ...
    'numEstimator', 1, ...
    'stopCriterion', 0.05);

% Compute psychometric function for the 4AFC paradigm with the 4 E scenes
[threshold, questObj, psychometricFunction, fittedPsychometricParams] = computeThreshold(...
    tumblingEsceneEngines, theNeuralEngine, classifierEngine, ...
    classifierPara, thresholdPara, questEnginePara, ...
    'visualizeAllComponents', ~true, ...
    'verbose', true, ...
    'TAFC', false, 'useMetaContrast', false);

% Plot the derived psychometric function and other things.  The lower
% level routines put this in ISETBioJandJRootPath/figures.
% pdfFileName = sprintf('Performance_%s_Reps_%d.pdf', strrep(params.psfDataFile, '.mat', ''), nTest);
% plotDerivedPsychometricFunction(questObj, threshold, fittedPsychometricParams, ...ISETBio
%     thresholdParameters, fullfile(params.outputFiguresDir,pdfFileName), 'xRange', [0.02 0.2]);
% if (params.visualEsOnMosaic)
%     pdfFileName = sprintf('Simulation_%s_Reps_%d.pdf', strrep(params.psfDataFile, '.mat', ''), nTest);
%     visualizeSimulationResults(questObj, threshold, fittedPsychometricParams, ...
%         thresholdParameters, tumblingEsceneEngines, theNeuralEngine, ...
%         fullfile(params.outputFiguresDir,pdfFileName));
% end

% % Export the results
% exportFileName = sprintf('Results_%s_Reps_%d.mat', strrep(params.psfDataFile, '.mat', ''), nTest);
% if (~isempty(params.customMacularPigmentDensity))
%     exportFileName = strrep(exportFileName, '.mat', sprintf('_MPD_%2.2f.mat', params.customMacularPigmentDensity));
% end
% if (~isempty(params.customPupilDiameterMM))
%     exportFileName = strrep(exportFileName, '.mat', sprintf('_PupilDiamMM_%2.2f.mat', params.customPupilDiameterMM));
% end
% if (~isempty(params.customConeDensities))
%     exportFileName = strrep(exportFileName, '.mat', sprintf('_cones_%2.2f_%2.2f_%2.2f.mat', params.customConeDensities(1), params.customConeDensities(2), params.customConeDensities(3)));
% end
% if (~isempty(params.customLensAgeYears))
%     exportFileName = strrep(exportFileName, '.mat', sprintf('_lensAge_%d.mat', params.customLensAgeYears));
% end
% 
% fprintf('Saving data to %s\n', fullfile(params.outputResultsDir,exportFileName));
% exportSimulation(questObj, threshold, fittedPsychometricParams, ...
%     thresholdParameters, classifierPara, questEnginePara, ...
%     tumblingEsceneEngines, theNeuralEngine, classifierEngine, ...
%     fullfile(params.outputResultsDir,exportFileName));

% logMAR(iPSF) = log10(threshold(iPSF)*60/5);

% Save summary,  This allows examination of the numbers and/or
% replotting.
% save(fullfile(params.outputResultsDir,summaryFileName),"examinedPSFDataFiles","threshold","logMAR","LCA","TCA","theConeMosaic");

end
