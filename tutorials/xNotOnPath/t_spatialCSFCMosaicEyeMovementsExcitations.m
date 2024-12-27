% Compute spatial CSF in the achromatic direction in the presence of fEM
%
% Description:
%    Use ISETBioCSFGenerator to run out CSFs in the achromatic color direction
%    in the presence of fixatinal eye movements using the cone excitations
%    signal.
%
%    This example uses the new @cMosaic object, an SVM computational observer,
%    operating on the output of a quadrature energy spatial pooling mechanism
%    and square windowed gratings of constant size.
%
% See also: t_spatialCSF, t_modulatedGratingsSceneGeneration,
%           t_chromaticThresholdContour, computeThreshold, computePerformance
%

% History:
%   8/24/21  NPC   Wrote it

% Clear and close
clear; close all;

% Stimulus mean luminance
meanLuminanceCdM2 = 34;

% List of spatial frequencies to be tested.
spatialFreqs = [30]; % [2, 4, 8, 12, 16, 25];

% Achromatic color direction specified as a 1-by-3 vector
chromaDir = [0.1, -0.1, 0];%0.7*[1.0, 1.0, 1.0]';

%% Create neural response engine
%
theConeMosaic = cMosaic(...
        'whichEye', 'right eye', ...
        'sizeDegs', [1 1]*0.5*30/spatialFreqs(1), ...
        'eccentricityDegs', [0 0], ...
        'integrationTime', 5/1000, ...
        'eccVaryingConeAperture', true, ...
        'eccVaryingOuterSegmentLength', true, ...
        'eccVaryingConeBlur', false, ...
        'eccVaryingMacularPigmentDensity', true, ...
        'eccVaryingMacularPigmentDensityDynamic', false);
  
zernikeDataBase = 'Polans2015';
rankedSubjectIDs = PolansOptics.constants.subjectRanking;
subjectID = rankedSubjectIDs(1);
subtractCentralRefraction = PolansOptics.constants.subjectRequiresCentralRefractionCorrection(subjectID);

oiEnsemble = theConeMosaic.oiEnsembleGenerate(theConeMosaic.eccentricityDegs, ...
            'zernikeDataBase', zernikeDataBase, ...
            'subjectID', subjectID, ...
            'pupilDiameterMM', 3.0, ...
            'subtractCentralRefraction', subtractCentralRefraction, ...
            'zeroCenterPSF', true, ...
            'wavefrontSpatialSamples', 301 ...
            );
theOptics = oiEnsemble{1};      
  
% The neural compute function to employ
neuralComputeFunction = @nrePbotopigmentExcitationsCMosaic;

% Obtain default neural response engine params
neuralParams = neuralComputeFunction();

% Update the eyeMovementsParams to simulate a different model
neuralParams.eyeMovementsParams.durationSeconds = 100/1000;
neuralParams.eyeMovementsParams.driftModel = 'high velocity';    % Choose between {'low velocity', 'high velocity', 'default'}

% Instantiate the neural engine with custom optics, custom cone mosaic and
% custom eye movement parms using cone photocurrents as the output signal
theNeuralEngine = neuralResponseEngine(neuralComputeFunction, neuralParams);

% Update the cone mosaic and optics using custom cone mosaic and optics
theNeuralEngine.customNeuralPipeline(struct(...
          'coneMosaic', theConeMosaic, ...
          'optics', theOptics));
      
%% Instantiate the PoissonTAFC responseClassifierEngine
%
% The computational observer consists of a spatial pooling mechanism combined with 
% an SVM classifier operating on the output of a quadrature-energy spatial pooling mechanism 

% The spatial pooling struct
poolingParamsStruct = struct(...
    'type', 'quadratureEnergy', ...
    'weights', [] ... % pooling weights are computed separately for each stimulus and test contrast in computeThreshold()
    );

classifierEngineParamsStruct = struct(...
        'pooling', poolingParamsStruct, ... % Spatial pooling params struct
        'crossValidationFoldsNum', 10, ...  % employ a 10-fold cross-validated linear 
        'kernelFunction', 'linear', ...     % linear
        'classifierType', 'svm' ...         % binary SVM classifier
    );
classifierEngine = responseClassifierEngine(@rcePoolingSVMTAFC, classifierEngineParamsStruct);

% Train classifier using a set of 256 noisy instances,
nTrain = 256;
nTest = 256;
classifierPara = struct('trainFlag', 'random', ...
                        'testFlag', 'random', ...
                        'nTrain', nTrain, 'nTest', nTest);
                            
%% Parameters for threshold estimation/quest engine
% The actual threshold varies enough with the different engines that we
% need to adjust the contrast range that Quest+ searches over, as well as
% the range of psychometric function slopes. Threshold limits are computed
% as 10^-logThreshLimitVal.
thresholdPara = struct('logThreshLimitLow', 2.4, ...
                       'logThreshLimitHigh', 0.04, ...
                       'logThreshLimitDelta', 0.02, ...
                       'slopeRangeLow', 10, ...
                       'slopeRangeHigh', 500, ...
                       'slopeDelta', 1.0);

% Parameter for running the QUEST+
% See t_spatialCSF.m for more on options of the two different mode of
% operation (fixed numer of trials vs. adaptive)
contrastLevelsExamined = 8;
questEnginePara = struct('minTrial', nTest*contrastLevelsExamined, ...
                         'maxTrial', nTest*contrastLevelsExamined, ...
                         'numEstimator', 1, 'stopCriterion', 0.05);

%% Compute threshold for each spatial frequency
% 
% See toolbox/helpers for functions createGratingScene computeThreshold
dataFig = figure();
logThreshold = zeros(1, length(spatialFreqs));
for idx = 1:length(spatialFreqs)
    % Create a static grating scene engine with a particular chromatic direction,
    % spatial frequency, and temporal duration
    gratingSceneEngine = createGratingScene(chromaDir, spatialFreqs(idx), ...
        'meanLuminanceCdPerM2', meanLuminanceCdM2, ...
        'spatialEnvelope', 'rect', ...
        'fovDegs', max(theConeMosaic.sizeDegs)*1.02, ...
        'spatialEnvelopeRadiusDegs', max(theConeMosaic.sizeDegs)*1.02, ...
        'warningInsteadOfErrorOnOutOfGamut', false, ...
        'minPixelsNumPerCycle', 5);
    
    % Compute the threshold for our grating scene with the previously
    % defined neural and classifier engine.  This function does a lot of
    % work, see t_tresholdEngine and the function itself, as well as
    % function computePerformance.
    [logThreshold(idx), questObj] = ...
        computeThreshold(gratingSceneEngine, theNeuralEngine, classifierEngine, ...
        classifierPara, thresholdPara, questEnginePara, ...
        'visualizeAllComponents', ~true, ...
        'TAFC', true);
    
    % Plot stimulus
    figure(dataFig);
    subplot(4, 4, idx * 2 - 1);
    
    visualizationContrast = 0.7;
    [theSceneSequence] = gratingSceneEngine.compute(visualizationContrast);
    gratingSceneEngine.visualizeStaticFrame(theSceneSequence);
    
    % Plot data and psychometric curve 
    % with a marker size of 2.5
    subplot(4, 4, idx * 2);
    %questObj.plotMLE(2.5);
    [threshold, para, dataOut{idx}] = questObj.thresholdMLE(...
        'returnData', true, ...
        'showPlot', true, 'newFigure', false, 'pointSize', 2.5);

end
set(dataFig, 'Position',  [0, 0, 800, 800]);

% Convert returned log threshold to linear threshold
threshold = 10 .^ logThreshold;

figure();
d = dataOut{1}; plot(d.examinedContrasts, d.pCorrect, 'ks'); hold on; plot(d.examinedContrastsFit, d.pCorrectFit, 'r-')

%% Plot Contrast Sensitivity Function
theCsfFig = figure();
loglog(spatialFreqs, 1 ./ threshold, '-ok', 'LineWidth', 2);
xlabel('Spatial Frequency (cyc/deg)');
ylabel('Sensitivity');
set(theCsfFig, 'Position',  [800, 0, 600, 800]);
