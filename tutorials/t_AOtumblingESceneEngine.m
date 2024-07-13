function [sce0,sce90,sce180,sce270,sceBg] = t_AOtumblingESceneEngine(varargin)

% Initialize
close all;

% Parse optional input
p = inputParser;
p.addParameter('visualizeScene', true, @islogical);
p.parse(varargin{:});
visualizeScene = p.Results.visualizeScene;

% Make sure figures directory exists so that output writes
% don't fail
rootPath = ISETBioCSFGeneratorRootPath;
if (~exist(fullfile(rootPath,'local','figures'),'dir'))
    mkdir(fullfile(rootPath,'local','figures'));
end

% Load in a monitor that mimics the primaries in the Berkely AO system.
theDisplay = load(fullfile(rootPath,'sampledata','monoDisplay.mat'));
presentationDisplay = theDisplay.monoDisplay;

% % These are (or might be) some of the key parameters for the Berkeley
% % AO experiments.
% %
% % Other stimulus parameters
% %
% % Define basic parameters of the AO stimulus,
% % We can put in arbitrary spectra but the
% % actual E we want to model was monochromatic.
% wls = 600:10:800;
% eWavelengthNm = 780;
%
% % Spatial parameters
% nPixels = 128;
% fieldSizeMinutes = 60;
% fieldSizeDegs = fieldSizeMinutes/60;
%
% % Background power in uW.  Depending on how Will ships this to us,
% % we may need to convert into these units.
% fieldPowerUWPerDeg2 = 2000;
% fieldPowerUW = (fieldSizeDegs^2)*fieldPowerUWPerDeg2;
%
% % E power in uW.  This is the power for the stimulus that
% % is added to the background.  Let's assume for now that
% % it is twice the background power.
% ePowerUW = 2*fieldPowerUW;

% Set the parameters for the AO mimicing display
sceneParams = sceTumblingEscene;
sceneParams.presentationDisplay = presentationDisplay;
sceneParams.plotDisplayCharacteristics = false;

% Define the foreground (E) and background RGB wrt the monochromatic
% monitor.  I think the Berkeley experiments use only the green
% channel, so we'll do that here.
sceneParams.chromaSpecification.backgroundRGB = [0 0 0];
sceneParams.chromaSpecification.foregroundRGB = [0 0.5 0];

% Instantiate a tumblingEsceneEngine for 0 deg rotation E
sceneParams.letterRotationDegs = 0;
tumblingEsceneEngine0degs = sceneEngine(@sceTumblingEscene,sceneParams);

% Generate sceneEngine for 90 deg rotation E
sceneParams.letterRotationDegs = 90;
tumblingEsceneEngine90degs = sceneEngine(@sceTumblingEscene,sceneParams);

% Generate sceneEngine for 180 deg rotation E
sceneParams.letterRotationDegs = 180;
tumblingEsceneEngine180degs = sceneEngine(@sceTumblingEscene,sceneParams);

% Generate sceneEngine for 270 deg rotation E
sceneParams.letterRotationDegs = 270;
tumblingEsceneEngine270degs = sceneEngine(@sceTumblingEscene,sceneParams);

% Generate params for the background scene
backgroundSceneParams = sceneParams;
backgroundSceneParams.chromaSpecification.foregroundRGB = sceneParams.chromaSpecification.backgroundRGB;
backgroundSceneEngine = sceneEngine(@sceTumblingEscene,backgroundSceneParams);

% Generate scenes with size of 0.1 deg
sizeDegs = 0.1;
theSmallEsceneSequence0degs = tumblingEsceneEngine0degs.compute(sizeDegs);
theSmallEsceneSequence90degs = tumblingEsceneEngine90degs.compute(sizeDegs);
theSmallEsceneSequence180degs = tumblingEsceneEngine180degs.compute(sizeDegs);
theSmallEsceneSequence270degs = tumblingEsceneEngine270degs.compute(sizeDegs);
theSmallBackgroundSceneSequence = backgroundSceneEngine.compute(sizeDegs);

% Get first frame of the scene sequences
theSmallEscene0degs = theSmallEsceneSequence0degs{1};
theSmallEscene90degs = theSmallEsceneSequence90degs{1};
theSmallEscene180degs = theSmallEsceneSequence180degs{1};
theSmallEscene270degs = theSmallEsceneSequence270degs{1};
theSmallBackgroundScene = theSmallBackgroundSceneSequence{1};

if (visualizeScene)
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
        'rowsNum', 1, ...
        'colsNum', 4, ...
        'heightMargin',  0.0, ...
        'widthMargin',    0.05, ...
        'leftMargin',     0.05, ...
        'rightMargin',    0.00, ...
        'bottomMargin',   0.05, ...
        'topMargin',      0.00);

    domainVisualizationLimits = 0.3*0.5*[-1 1 -1 1];

    hFig = figure(1);
    clf;
    set(hFig, 'Position', [10 10 1550 400], 'Color', [1 1 1]);
    ax = subplot('Position', subplotPosVectors(1,1).v);
    visualizeScene(theSmallEscene0degs, ...
        'spatialSupportInDegs', true, ...
        'crossHairsAtOrigin', true, ...
        'displayRadianceMaps', false, ...
        'avoidAutomaticRGBscaling', true, ...
        'noTitle', true, ...
        'axesHandle', ax);
    set(ax, 'XLim', [domainVisualizationLimits(1) domainVisualizationLimits(2)], ...
        'YLim', [domainVisualizationLimits(3) domainVisualizationLimits(4)], ...
        'XTick', [-0.1 0 0.1], 'YTick', [-0.1 0 0.1] ...
        );

    ax = subplot('Position', subplotPosVectors(1,2).v);
    visualizeScene(theSmallEscene90degs, ...
        'spatialSupportInDegs', true, ...
        'crossHairsAtOrigin', true, ...
        'displayRadianceMaps', false, ...
        'avoidAutomaticRGBscaling', true, ...
        'noYLabel', true, ...
        'noYTicks', true, ...
        'noTitle', true, ...
        'axesHandle', ax);
    set(ax, 'XLim', [domainVisualizationLimits(1) domainVisualizationLimits(2)], ...
        'YLim', [domainVisualizationLimits(3) domainVisualizationLimits(4)], ...
        'XTick', [-0.1 0 0.1], 'YTick', [-0.1 0 0.1] ...
        );

    ax = subplot('Position', subplotPosVectors(1,3).v);
    visualizeScene(theSmallEscene180degs, ...
        'spatialSupportInDegs', true, ...
        'crossHairsAtOrigin', true, ...
        'displayRadianceMaps', false, ...
        'avoidAutomaticRGBscaling', true, ...
        'noYLabel', true, ...
        'noYTicks', true, ...
        'noTitle', true, ...
        'axesHandle', ax);
    set(ax, 'XLim', [domainVisualizationLimits(1) domainVisualizationLimits(2)], ...
        'YLim', [domainVisualizationLimits(3) domainVisualizationLimits(4)], ...
        'XTick', [-0.1 0 0.1], 'YTick', [-0.1 0 0.1] ...
        );

    ax = subplot('Position', subplotPosVectors(1,4).v);
    visualizeScene(theSmallEscene270degs, ...
        'spatialSupportInDegs', true, ...
        'crossHairsAtOrigin', true, ...
        'displayRadianceMaps', false, ...
        'avoidAutomaticRGBscaling', true, ...
        'noYLabel', true, ...
        'noYTicks', true, ...
        'noTitle', true, ...
        'axesHandle', ax);
    set(ax, 'XLim', [domainVisualizationLimits(1) domainVisualizationLimits(2)], ...
        'YLim', [domainVisualizationLimits(3) domainVisualizationLimits(4)], ...
        'XTick', [-0.1 0 0.1], 'YTick', [-0.1 0 0.1] ...
        );

    projectBaseDir = ISETBioCSFGeneratorRootPath;
    pdfFile = fullfile(projectBaseDir,'figures','t_AOTumblingSceneEngine_stimuli.pdf');
    NicePlot.exportFigToPDF(pdfFile,hFig, 300);
end

% Set output if any requested
if (nargout > 0)
    sce0 = tumblingEsceneEngine0degs;
    sce90 = tumblingEsceneEngine90degs;
    sce180 = tumblingEsceneEngine180degs;
    sce270 = tumblingEsceneEngine270degs;
    sceBG = backgroundSceneEngine;
end

end

