function t_tumblingEsceneEngine()

    % Initialize
    clear; close all;
    
    % Make sure figures and results directories exist so that output writes
    % don't fail
    rootPath = ISETBioCSFGeneratorRootPath;
    if (~exist(fullfile(rootPath,'local','figures'),'dir'))
        mkdir(fullfile(rootPath,'local','figures'));
    end

    % Obtain the default params for the tumblingEscene engine
    theSceneEngine = createTumblingEsceneEngine(0);
    defaultParams = theSceneEngine.sceneComputeFunction();
    
    % Set the plotDisplayCharacteristics field to true to visualize the
    % display SPDs and gamma functions
    customSceneParams = defaultParams;
    customSceneParams.plotDisplayCharacteristics = true;
    % Specify datafile name for the display SPDs
    customSceneParams.spdDataFile = 'BVAMS_White_Guns_At_Max.mat';
    
    % Generate sceneEngine for 0 deg rotation E
    letterRotationDegs = 0;
    tumblingEsceneEngine0degs = createTumblingEsceneEngine(letterRotationDegs, 'customSceneParams', customSceneParams);

    % Generate sceneEngine for 90 deg rotation E
    letterRotationDegs = 90;
    tumblingEsceneEngine90degs = createTumblingEsceneEngine(letterRotationDegs);


    % Generate sceneEngine for 180 deg rotation E
    letterRotationDegs = 180;
    tumblingEsceneEngine180degs = createTumblingEsceneEngine(letterRotationDegs);

    % Generate sceneEngine for 270 deg rotation E
    letterRotationDegs = 270;
    tumblingEsceneEngine270degs = createTumblingEsceneEngine(letterRotationDegs);


    % Generate params for the background scene
    sceneParams = tumblingEsceneEngine0degs.sceneComputeFunction();
    backgroundSceneParams = sceneParams;
    backgroundSceneParams.chromaSpecification.foregroundRGB = sceneParams.chromaSpecification.backgroundRGB;
    backgroundSceneEngine = createTumblingEsceneEngine(letterRotationDegs, 'customSceneParams', backgroundSceneParams);

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
    pdfFile = fullfile(projectBaseDir,'figures','t_tumblingEsceneEngine_stimuli.pdf');
    NicePlot.exportFigToPDF(pdfFile,hFig, 300);
end

function tumblingEsceneEngine = createTumblingEsceneEngine(orientation, varargin)

    % Parse
    p = inputParser;
    p.addParameter('customSceneParams', [], @(x)(isempty(x)||(isstruct(x))));
    p.parse(varargin{:});
    customSceneParams = p.Results.customSceneParams;

    % Handle to the compute function which will compute a new scene
    % which varies in the size of 'E' letter (the variable for which 
    % we assess performance)
    sceneComputeFunction = @sceTumblingEscene;

    if (isempty(customSceneParams))
        % Retrieve the default params for the tumbingE scene
        sceneParams = sceneComputeFunction();
    else
        sceneParams = customSceneParams;
    end


    % Change the orientation to the passed orientation
    sceneParams.letterRotationDegs = orientation;

    % Instantiate a tumblingEsceneEngine
    tumblingEsceneEngine = sceneEngine(sceneComputeFunction, sceneParams);
end

