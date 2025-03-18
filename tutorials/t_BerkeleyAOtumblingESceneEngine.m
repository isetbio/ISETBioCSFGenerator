function [sce0,sce90,sce180,sce270,sceBg,sceneParams] = t_BerkeleyAOtumblingESceneEngine(options)
% Show how to describe and control a tumbling E stimulus as presented in
% the Berkeley (Tuten lab) AO system.
%
% Tuten's description of the stimulus parameters for this experiment.
%   - Wavelength: the tumbling-E was modulated in the imaging channel. The
%   center wavelength is 840 nm. We don't have a spectroradiometer that
%   operates in this part of the spectrum, so we haven't measured the
%   spectral bandwidth directly, but the filter used to select light in
%   this channel is this one:
%   https://www.idex-hs.com/store/product-detail/ff01_840_12_25/fl-004460 -
%   Display size: 512 by 512 pixels, 1.413 by 1.413 degrees, 362.3
%   pixels/degree.
%
%   - The E was presented as a decrement against the 840 nm imaging
%   channel. The contrast is nominally 100%, but this in practice this will
%   be slightly less due to (a) imperfect extinction by the AOM and (b)
%   background 940 nm light that is used for wavefront sensing and is not
%   modulated.
%
%   - The average 840 nm power measured at the cornea was 141.4 microwatts.
%   This would be distributed more or less uniformly over the 1.413 degree
%   square field described above. Assuming a 16.67 mm eye length, I'm
%   getting 8.37E-04 W/mm2, but you might want to check my math on that.
%
%   - Temporal characteristics of the AOSLO: 30 fps, 50 ns per pixel. This
%   means each pixel in the raster is illuminated very briefly every ~33
%   ms. Depending on the size of the cone, the effective pulse of light it
%   receives every frame will be slightly longer since the AOSLO pixels are
%   smaller than the cones themselves. At this AOSLO field size, a cone
%   with a 1 arcmin diameter will receive stimulus light over the course of
%   ~0.40 ms per frame. The timing of stimulation might be relevant for
%   this modeling, so happy to discuss the details further.
%
% Other potentially relevant details:
%   - Test eccentricity: 1 degree (usually temporal)
%   - Stimulus duration: 3 AOSLO frames
%   - Stimulus size: 20/40 tumbling-E optotype (width of E bar = 2 arcmin; overall E size is 10 arcmin)
%   - Imposed stimulus motion on the retina
%   - Magnitude: we shifted the letter on each frame by increments of the
%   bar width -- either 0, 0.5, 1.0, or 2.0. "0" imposed motion in this
%   case means it was retinally stabilized
%   direction: shifts were imposed in retinal coordinates along the four
%   cardinal meridians Performance was assessed as correct for a fixed
%   letter size
%
% Spatial parameters
%   nPixels = 512;
%   fieldSizeMinutes = 1.413*60;
%   fieldSizeDegs = fieldSizeMinutes/60;
%
% Imaging (840 nm) power in uW.
%   fieldPowerUW = 141.4;
%   fieldPowerUWPerDeg2 = fieldPowerUW/(fieldSizeDegs^2);
%
% Need pupil diameter to convert corneal power to appropriate equivalent
% radiance.  Probably around 6 mm.
%
% This rouine takes a number of key value pairs. See arguments block below
% for a list and their defaults.

arguments
    options.visualizeScene (1,1) logical = true;
    options.displayNPixels (1,1) double = 512;
    options.displayFOVDeg (1,1) double = 1.413;
    options.wave (:,1) double = (500:5:870)';
    options.AOPrimaryWls (1,3) double = [840 683 543]; % [700 683 54];
    options.AOPrimaryFWHM (1,3) double = [22 27 23];
    options.AOCornealPowersUW (1,3) double = [141.4 10 10];
    options.ambientSpd (:,1) double = zeros(size((500:5:870)'));  % Adjust this if wave changes
    options.pupilSizeMM (1,1) double = 6;
    options.plotDisplayCharacteristics (1,1) logical = false;
    options.chromaSpecification_type (1,:) char = 'RGBsettings';
    options.chromaSpecification_backgroundRGB (1,:) double = [1 0 0];
    options.chromaSpecification_foregroundRGB (1,3) double = [0 0 0];
    options.eHeightMin (1,1) double = 30;
    options.temporalModulationParams_frameRateHz (1,1) double = 60;
    options.temporalModulationParams_numFrame (1,1) double = 3;
    options.temporalModulationParams_xShiftPerFrame (1,:) double = [0 10/60 0];
    options.temporalModulationParams_yShiftPerFrame (1,:) double = [0 0 10/60];
    options.temporalModulationParams_stimOnFrames (1,:) double = [1 1 1];
    options.temporalModulationParams_backgroundRGBPerFrame (:,:) double = [0 0 0; 1 0 0; 0 0 0];
    options.responseFlag (1,:) char = 'excitations';
    options.exportCondition (1,:) char = 'no change';
end

% Initialize
close all;

% Make sure figures directory exists so that output writes
% don't fail
rootPath = ISETBioCSFGeneratorRootPath;
if (~exist(fullfile(rootPath,'local',mfilename,'figures'),'dir'))
    mkdir(fullfile(rootPath,'local',mfilename,'figures'));
end

% Get the basic parameters for the AO scene engine and override according
% to passed key/value pairs.
sceneParams = sceBerkeleyAOTumblingEscene;
sceneParams.spectralSupport = options.wave; 

% Spatial
sceneParams.displayPixelSize = options.displayNPixels;
sceneParams.displayFOVDeg = options.displayFOVDeg; 

% Update parameters according to key/value pairs
sceneParams.displayParams.AOPrimaryWls = options.AOPrimaryWls; 
sceneParams.displayParams.AOPrimaryFWHM = options.AOPrimaryFWHM;
sceneParams.displayParams.AOCornealPowersUW= options.AOCornealPowersUW;
sceneParams.displayParams.ambientSpd = options.ambientSpd;
sceneParams.displayParams.plotDisplayCharacteristics = options.plotDisplayCharacteristics;
sceneParams.pupilSizeMM = options.pupilSizeMM;

% Create spectral power distribution for display from updated parameters
for pp = 1:length(sceneParams.displayParams.AOPrimaryWls)
    sceneParams.displayParams.spd(:,pp) = generateAOPrimarySpd(sceneParams.spectralSupport, ...
        sceneParams.displayParams.AOPrimaryWls(pp), ...
        sceneParams.displayParams.AOPrimaryFWHM(pp), ...
        sceneParams.displayParams.AOCornealPowersUW(pp), ...
        sceneParams.displayParams.displayFOVDeg, ...
        sceneParams.pupilSizeMM);
end

% Define the foreground (E) and background RGB wrt the monochromatic
% monitor.  The Berkeley tumbling E experiments use only the 840
% channel, with background nominally full on and E nominally full off.
sceneParams.chromaSpecification.type = options.chromaSpecification_type;
sceneParams.chromaSpecification.backgroundRGB = options.chromaSpecification_backgroundRGB;
sceneParams.chromaSpecification.foregroundRGB = options.chromaSpecification_foregroundRGB;

% Specifiy E spatial parameters.  The scene engine only produces binary
% text characters so the size of the letter has to be an integer, and we
% want it to be an even integer so that we can put it in the center of an 
% image array that itself has an even number of pixels.
sceneParams.eHeightMin = options.eHeightMin;

% Define the E size temporal sequence of scenes that we want that includes
% how the E should move around across frames
testESizeMin = 10;
sceneParams.temporalModulationParams = struct(...   % temporal: modulation params struct
                'frameRateHz', options.temporalModulationParams_frameRateHz, ...              % frame rate in Hz
                'numFrames', options.temporalModulationParams_numFrame, ...                   % number of frames we want the E on for
                'xShiftPerFrame', options.temporalModulationParams_xShiftPerFrame, ...        % how much the E should be shifted in the x dimension in each frame, in degrees
                'yShiftPerFrame', options.temporalModulationParams_yShiftPerFrame,  ...       % how much the E should be shifted in the y dimension in each frame, in degrees
                'backgroundRGBPerFrame', options.temporalModulationParams_backgroundRGBPerFrame, ... % background color in RGB for each frame
                'stimOnFrames', options.temporalModulationParams_stimOnFrames ...
                );

% Instantiate a tumblingEsceneEngine for 0 deg rotation E
sceneParams.letterRotationDegs = 0;
tumblingEsceneEngine0degs = sceneEngine(@sceBerkeleyAOTumblingEscene,sceneParams);

% Generate sceneEngine for 90 deg rotation E
sceneParams.letterRotationDegs = 90;
tumblingEsceneEngine90degs = sceneEngine(@sceBerkeleyAOTumblingEscene,sceneParams);

% Generate sceneEngine for 180 deg rotation E
sceneParams.letterRotationDegs = 180;
tumblingEsceneEngine180degs = sceneEngine(@sceBerkeleyAOTumblingEscene,sceneParams);

% Generate sceneEngine for 270 deg rotation E
sceneParams.letterRotationDegs = 270;
tumblingEsceneEngine270degs = sceneEngine(@sceBerkeleyAOTumblingEscene,sceneParams);

% Generate params for the background scene
backgroundSceneParams = sceneParams;
backgroundSceneParams.chromaSpecification.foregroundRGB = sceneParams.chromaSpecification.backgroundRGB;
backgroundSceneEngine = sceneEngine(@sceBerkeleyAOTumblingEscene,backgroundSceneParams);

% Generate scenes with size of 0.1 deg
testESizeDeg = testESizeMin/60;
theSmallEsceneSequence0degs = tumblingEsceneEngine0degs.compute(testESizeDeg);
theSmallEsceneSequence90degs = tumblingEsceneEngine90degs.compute(testESizeDeg);
theSmallEsceneSequence180degs = tumblingEsceneEngine180degs.compute(testESizeDeg);
theSmallEsceneSequence270degs = tumblingEsceneEngine270degs.compute(testESizeDeg);
theSmallBackgroundSceneSequence = backgroundSceneEngine.compute(testESizeDeg);

% Visualize each frame ofthe scene sequence
displaySizePixels = sceneParams.displayParams.displayPixelSize;
displaySizeDegrees = sceneParams.displayParams.displayFOVDeg;
pixelsPerDegree = displaySizePixels / displaySizeDegrees;

domainVisualizationLimits = [-displaySizeDegrees/2, displaySizeDegrees/2, -displaySizeDegrees/2, displaySizeDegrees/2];

for ff = 1:length(theSmallEsceneSequence0degs)
    theSmallEscene0degs = theSmallEsceneSequence0degs{ff};
    theSmallEscene90degs = theSmallEsceneSequence90degs{ff};
    theSmallEscene180degs = theSmallEsceneSequence180degs{ff};
    theSmallEscene270degs = theSmallEsceneSequence270degs{ff};
    theSmallBackgroundScene = theSmallBackgroundSceneSequence{ff};

    % Take a look at the scene
    if (options.visualizeScene)
        subplotPosVectors = NicePlot.getSubPlotPosVectors(...
            'rowsNum', 1, ...
            'colsNum', 4, ...
            'heightMargin',  0.0, ...
            'widthMargin',    0.05, ...
            'leftMargin',     0.05, ...
            'rightMargin',    0.00, ...
            'bottomMargin',   0.05, ...
            'topMargin',      0.00);

        domainVisualizationLimits = 1.5*[-1 1 -1 1]; % 0.3*0.5*[-1 1 -1 1];

        hFig = figure;
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
        % set(ax, 'XLim', [domainVisualizationLimits(1) domainVisualizationLimits(2)], ...
        %     'YLim', [domainVisualizationLimits(3) domainVisualizationLimits(4)], ...
        %     'XTick', [-0.1 0 0.1], 'YTick', [-0.1 0 0.1] ...
        %     );
        set(ax, 'XLim', [domainVisualizationLimits(1) domainVisualizationLimits(2)], ...
            'YLim', [domainVisualizationLimits(3) domainVisualizationLimits(4)], ...
            'XTick', linspace(-displaySizeDegrees/2, displaySizeDegrees/2, 5), ...
            'YTick', linspace(-displaySizeDegrees/2, displaySizeDegrees/2, 5));
        xTickLabels = arrayfun(@(x) sprintf('%.2f', x), linspace(-displaySizeDegrees/2, displaySizeDegrees/2, 5), 'UniformOutput', false);
        yTickLabels = arrayfun(@(y) sprintf('%.2f', y), linspace(-displaySizeDegrees/2, displaySizeDegrees/2, 5), 'UniformOutput', false);

        set(ax, 'XTickLabel', xTickLabels, 'YTickLabel', yTickLabels);
        set(ax, 'FontSize', 10, 'FontWeight', 'bold');
        
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
            'XTick', linspace(-displaySizeDegrees/2, displaySizeDegrees/2, 5), ...
            'YTick', linspace(-displaySizeDegrees/2, displaySizeDegrees/2, 5));
        xTickLabels = arrayfun(@(x) sprintf('%.2f', x), linspace(-displaySizeDegrees/2, displaySizeDegrees/2, 5), 'UniformOutput', false);
        yTickLabels = arrayfun(@(y) sprintf('%.2f', y), linspace(-displaySizeDegrees/2, displaySizeDegrees/2, 5), 'UniformOutput', false);

        set(ax, 'XTickLabel', xTickLabels, 'YTickLabel', yTickLabels);
        set(ax, 'FontSize', 10, 'FontWeight', 'bold');

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
            'XTick', linspace(-displaySizeDegrees/2, displaySizeDegrees/2, 5), ...
            'YTick', linspace(-displaySizeDegrees/2, displaySizeDegrees/2, 5));
        xTickLabels = arrayfun(@(x) sprintf('%.2f', x), linspace(-displaySizeDegrees/2, displaySizeDegrees/2, 5), 'UniformOutput', false);
        yTickLabels = arrayfun(@(y) sprintf('%.2f', y), linspace(-displaySizeDegrees/2, displaySizeDegrees/2, 5), 'UniformOutput', false);

        set(ax, 'XTickLabel', xTickLabels, 'YTickLabel', yTickLabels);
        set(ax, 'FontSize', 10, 'FontWeight', 'bold');

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
            'XTick', linspace(-displaySizeDegrees/2, displaySizeDegrees/2, 5), ...
            'YTick', linspace(-displaySizeDegrees/2, displaySizeDegrees/2, 5));
        xTickLabels = arrayfun(@(x) sprintf('%.2f', x), linspace(-displaySizeDegrees/2, displaySizeDegrees/2, 5), 'UniformOutput', false);
        yTickLabels = arrayfun(@(y) sprintf('%.2f', y), linspace(-displaySizeDegrees/2, displaySizeDegrees/2, 5), 'UniformOutput', false);

        set(ax, 'XTickLabel', xTickLabels, 'YTickLabel', yTickLabels);
        set(ax, 'FontSize', 10, 'FontWeight', 'bold');

        projectBaseDir = ISETBioCSFGeneratorRootPath;
        pdfFile = fullfile(projectBaseDir,'local',mfilename,'figures',sprintf('t_AOTumblingSceneEngine_stimuli_frame%d.pdf',ff));
        NicePlot.exportFigToPDF(pdfFile,hFig, 300);
    end
end

% Set output if any requested
if (nargout > 0)
    sce0 = tumblingEsceneEngine0degs;
    sce90 = tumblingEsceneEngine90degs;
    sce180 = tumblingEsceneEngine180degs;
    sce270 = tumblingEsceneEngine270degs;
    sceBg = backgroundSceneEngine;
end

end

