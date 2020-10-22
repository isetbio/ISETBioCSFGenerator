% Compute spatial CSF in different color directions
%
% Description:
%    Use ISETBioCSFGenerator to run out CSFs in different color directions.
%    This example uses an ideal Poisson TAFC observer and circularly
%    windowed gratings of constant size.
%
% See also: t_thresholdEngine, t_modulatedGratingsSceneGeneration
%

% History:
%   10/20/20  lqz   Wrote it.
%   10/21/20  dhb  More commments.
%   10/22/20  lqz   Restructure the code

%% CSF Calculation
%
% List of spatial frequencies to be tested.
spatialFreqs = [0.5, 1, 2, 4, 8, 12, 16, 25];

% Allocate space for thresholds
logThreshold = zeros(1, length(spatialFreqs));

% Choose stimulus chromatic direction specified as a 1-by-3 vector
% of L, M, S cone contrast.  These vectors get normalized below, so only
% their direction matters in the specification.
stimType = 'luminance';
switch (stimType)
    case 'luminance'
        chromaDir = [1.0, 1.0, 1.0];
    case 'red-green'
        chromaDir = [1.0, -1.0, 0.0];
    case 'L-isolating'
        chromaDir = [1.0, 0.0, 0.0];
end

% Set the RMS cone contrast of the stimulus. Things may go badly if you
% exceed the gamut of the monitor, so we are conservative and set this at a
% value that is within gamut of typical monitors and don't worry about it
% further for this tutorial.  A vector length contrast of 0.08 should be
% OK.
rmsContrast = 0.08;
chromaDir = chromaDir / norm(chromaDir) * rmsContrast;
assert(abs(norm(chromaDir) - rmsContrast) <= 1e-10);

%% Create neural engine
% Instantiate a neuralResponseEngine
neuralParams = nrePhotopigmentExcitationsWithNoEyeMovements;
neuralParams.coneMosaicParams.fovDegs = 0.25;
theNeuralEngine = neuralResponseEngine(@nrePhotopigmentExcitationsWithNoEyeMovements, neuralParams);

%% Instantiate the PoissonTAFC responseClassifierEngine
% PoissonTAFC makes decision by performing the Poisson likelihood ratio test
classifierEngine = responseClassifierEngine(@rcePoissonTAFC);

% Parameter associated with this classifier
classifierPara = struct('trainFlag', 'none', ...
                                           'testFlag', 'random', ...
                                            'nTrain', 1, 'nTest', 128);

%% Parameter for threshold estimation/quest engine
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
% See t_thresholdEngine.m for options of the two different mode of
% operation (fixed numer of trials vs. adaptive)
questEnginePara = struct('minTrial', 1280, 'maxTrial', 1280, ...
                                                 'numEstimator', 1, 'stopCriterion', 0.05);

%% Compute threshold for each spatial frequency
% See toolbox/helpers for the definitin of
% function 'createGratingScene' and
% function 'computeThreshold'

dataFig = figure();
for idx = 1:length(spatialFreqs)
    % create a static grating scene with a particular chromatic direction,
    % spatial frequency, and temporal duration
    gratingScene = createGratingScene('chromaDir', chromaDir, 'spatialFreq', spatialFreqs(idx));
    
    % compute the threshold for our garting scene with the previously
    % defined neural and classifier engine
    [logThreshold(idx), questObj] = ...
        computeThreshold(gratingScene, theNeuralEngine, classifierEngine, classifierPara, thresholdPara, questEnginePara);
    
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
