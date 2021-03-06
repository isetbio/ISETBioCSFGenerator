% Compute spatial CSF in different color directions, using ISETBio midget RGC
%
% Description:
%    Use ISETBioCSFGenerator to run out CSFs in different color directions.
%    This example uses circularly windowed gratings of constant size.
%
% See also: t_chromaticThresholdContourMidgetRGC, t_thresholdEngine, t_modulatedGratingsSceneGeneration,
%           t_chromaticThresholdContour, computeThresholdTAFC, computePerformanceTAFC
%

% History:
%   11/07/20  dhb   Wrote by combining t_spatialCSF and t_chromaticThresholdContourMidgetRGC

% Clear and close
clear; close all;

% List of spatial frequencies to be tested.
spatialFreqs = [0.5, 1, 2, 4, 8, 12, 16];

% 100 msec stimulus duration
stimulusDurationSeconds = 100/1000;

% Options for presentationMode are {'drifted', 'flashed'}
presentationMode = 'drifted';

% How motion is sampled, 45 degs = 8 spatial phases/period
spatialPhaseAdvanceDegs = 45;  

% Temporal frequency in Hz
temporalFrequencyHz = 5;    

% For 'flashed' presentation mode, present the grating at 90 spatial phase 
% (odd symmetry).  Using 90 degree (sine phase) makes the stimulus symmetric 
% in terms of balanced incremental and decremental components, so that the 
% isothreshold contour is also symmetric.
gratingPhaseDeg = 90;

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
% Obtain default neural response engine params
neuralParams = nreMidgetRGC;

% Modify mRGC mosaic eccentricity and size
neuralParams.mRGCmosaicParams.eccDegs = [1 0];
neuralParams.mRGCmosaicParams.sizeDegs = 0.5*[1 1];

% *** POST-CONE SUMMATION NOISE ***
% Set the mRGC mosaic (post-cone summation) noise flag. If set to 'none',
% the only noise in the computation is that of the coneMosaic. 
%
% If set to 'random', Gaussian noise is added at the final mRGC response.
% The noise sd is noiseFactor*maxResponse
neuralParams.coneMosaicParams.noiseFlag = 'none';
neuralParams.mRGCmosaicParams.noiseFlag = 'random';
neuralParams.mRGCmosaicParams.noiseFactor = 0.25;

% Modify some cone mosaic params
neuralParams.coneMosaicParams.coneMosaicResamplingFactor = 3;
neuralParams.coneMosaicParams.integrationTime = 100/1000;

% Instantiate the neural response engine
theNeuralEngine = neuralResponseEngine(@nreMidgetRGC, neuralParams);

%% Instantiate the PoissonTAFC or the PcaSVMTAFC responseClassifierEngine
%
% Options:
%   'idealObserver' - ideal TAFC observer for Poisson limited signals
%   'computationalObsever' - SVM based learned classifier
classifierChoice = 'idealObserver';
switch (classifierChoice) 
    case 'idealObserver'
        % PoissonTAFC makes decision by performing the Poisson likelihood ratio test
        % Also set up parameters associated with use of this classifier.
        classifierEngine = responseClassifierEngine(@rcePoissonTAFC);
        
        % Train classifier using 1 noise-free instance, 
        % Test performance using a set of 128 noisy instances
        classifierPara = struct('trainFlag', 'none', ...
                                'testFlag', 'random', ...
                                'nTrain', 1, 'nTest', 128);
    case 'computationalObserver'
        % Instantiate a computational observer consisting of a linear SVM 
        % coupled with a PCA operating on 2 components
        classifierEngine = responseClassifierEngine(@rcePcaSVMTAFC, ...
            struct(...
                'PCAComponentsNum', 2, ...          % number of PCs used for feature set dimensionality reduction
                'crossValidationFoldsNum', 10, ...  % employ a 10-fold cross-validated linear 
                'kernelFunction', 'linear', ...     % linear
                'classifierType', 'svm' ...         % binary SVM classifier
            ));
        % Train classifier using a set of 256 noisy instances, 
        % Test performance using a set of 128 noisy instances
        classifierPara = struct('trainFlag', 'random', ...
                                'testFlag', 'random', ...
                                'nTrain', 256, 'nTest', 128);
                        
    otherwise
        error('Unknown classifier: ''%s''.', classifierChoice);
end

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


%% Compute threshold for each spatial frequency
% 
% See toolbox/helpers for functions createGratingScene computeThresholdTAFC
dataFig = figure();
logThreshold = zeros(1, length(spatialFreqs));
for idx = 1:length(spatialFreqs)
    % Create a static grating scene with a particular chromatic direction,
    % spatial frequency, and temporal duration. Make it larger than
    % the mRGC mosaic so that it extends over cone inputs to  the surround
    % subregions of the RGC cells, which are quite large (~7 times the RF
    % center).
    maxEccDegs = max(neuralParams.mRGCmosaicParams.eccDegs) + max(0.5*neuralParams.mRGCmosaicParams.sizeDegs);
    extraDegsForRGCSurround = 2.0 * ...
        RGCmodels.CronerKaplan.constants.surroundCharacteristicRadiusFromFitToPandMcells(maxEccDegs);
    stimFOVdegs = max(neuralParams.mRGCmosaicParams.sizeDegs) + extraDegsForRGCSurround;
    
    
    % Options for presentationMode are {'sampled motion', 'flashed'}
    % For 'flashed' we make the duration equal to the
    gratingSceneEngine = createGratingScene(chromaDir, spatialFreqs(idx), ...
        'spatialPhase', gratingPhaseDeg, ...
        'duration', stimulusDurationSeconds, ...
        'temporalFrequencyHz', temporalFrequencyHz, ...
        'spatialPhaseAdvanceDegs', spatialPhaseAdvanceDegs, ...
        'fovDegs', stimFOVdegs, ...
        'spatialEnvelope', 'square', ...
        'presentationMode', presentationMode ...
        );
    
    % Compute the threshold for our grating scene with the previously
    % defined neural and classifier engine.  This function does a lot of
    % work, see t_tresholdEngine and the function itself, as well as
    % function computePerformanceTAFC.
    [logThreshold(idx), questObj] = ...
        computeThresholdTAFC(gratingSceneEngine, theNeuralEngine, classifierEngine, classifierPara, ...
        thresholdPara, questEnginePara);
    
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
