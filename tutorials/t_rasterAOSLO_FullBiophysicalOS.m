function t_rasterAOSLO_FullBiophysicalOS

    % Simulation parameters.
    simulateStimulusRaster = true;
     
    % How many frames / trial. 
    nStimulusFramesPerTrial = 2;

    % How many trials, which also means how many fEMs
    nTrials = 3;

    % Compute cone mosaic and retinal images of stimulus and background
    recomputeRetinalImages = ~true;

    % Compute cone excitations response
    recomputeConeExcitations = ~true;

    % Visualize the stimulus and the cone excitations response
    visualizeStimulusAndConeExcitationSequence = true;

    recomputePhotocurrents = ~true;


    % Where to output results, figures and videos
    projectBaseDir = ISETBioCSFGeneratorRootPath;
    if (~exist(fullfile(projectBaseDir,'local',mfilename,'figures'),'dir'))
        mkdir(fullfile(projectBaseDir,'local',mfilename,'figures'));
    end
    figureFileBase = fullfile(projectBaseDir,'local',mfilename,'figures');

    % Simulation results filename
    if (~exist(fullfile(projectBaseDir,'local',mfilename,'results'),'dir'))
        mkdir(fullfile(projectBaseDir,'local',mfilename,'results'));
    end
    simFileBase = fullfile(projectBaseDir,'local',mfilename,'results');
    

    if (simulateStimulusRaster)
        retinalImagesMatFileName = fullfile(simFileBase, 'coneMosaicAndRetinalImages0degsWithRaster.mat');
        coneExcitationsMatFileName = fullfile(simFileBase, 'coneExcitations0degsWithRaster.mat');
        photocurrentsMatFileName = fullfile(simFileBase, 'photoCurrents0degsWithRaster.mat');
        videoFilename = fullfile(figureFileBase, 'TumblingE0degsMosaicActivationWithRaster');
    else
        retinalImagesMatFileName = fullfile(simFileBase, 'coneMosaicAndRetinalImages0degsNoRaster.mat');
        coneExcitationsMatFileName = fullfile(simFileBase, 'coneExcitations0degsNoRaster.mat');
        photocurrentsMatFileName = fullfile(simFileBase, 'photoCurrents0degsNoRaster.mat');
        videoFilename = fullfile(figureFileBase, 'TumblingE0degsMosaicActivationNoRaster');
    end

    if (recomputeRetinalImages)
        % Generate scenes
        testESizeDeg = 10/60;
        [theNullScene, theEscene0degs] = generateScenes(testESizeDeg);

        pixelTimeOnSeconds = 50 * 1e-9;
        rasterLineDutyCycle = 0.4;

     
        % Simulation time step: here we pick a simulation time step during
        % which the AOSLO will scan through 8 of the 512 raster lines. This
        % computes to around 0.51 milliseconds
        fractionLinesScannedPerSimulationTimeStep = 8/512;
        [simulationTimeStepSeconds, stimulusRefreshIntervalSeconds] = ...
            computeSimulationTimeStep(pixelTimeOnSeconds, rasterLineDutyCycle, ...
                theEscene0degs, fractionLinesScannedPerSimulationTimeStep);


        % Set the mosaic integration time to the simulation time step
        mosaicIntegrationTimeSeconds = simulationTimeStepSeconds;
        
        % Mosaic size and eccentricity
        mosaicSizeDegs = 0.5*[1.0 1.0];
        mosaicEccDegs = [-1 0];

        % Generate the cone mosaic
        theConeMosaic = generateConeMosaic(mosaicEccDegs, mosaicSizeDegs, ...
                mosaicIntegrationTimeSeconds);


        % Generate optics
        simulateAdaptiveOpticsviewingConditions = true;
        theOI = generateOptics(simulateAdaptiveOpticsviewingConditions, theConeMosaic);
    
        % Whether to visualize the retinal images of the rasterized stimulus
        visualizeEachRetinalImage = true;

        % Generate the OIs for one period of the rasterized background
        fprintf('\nGenerating retinal images for the background raster. Please wait ...')
        theListOfBackgroundRasterRetinalImages = computeRasterScanRetinalImagesForOneFullRefresh(theOI, theNullScene, ...
                pixelTimeOnSeconds, rasterLineDutyCycle, ...
                stimulusRefreshIntervalSeconds, simulationTimeStepSeconds, ...
                visualizeEachRetinalImage);
        fprintf('Done !\n');

        fprintf('\nGenerating retinal images for the stimulus raster. Please wait ...')
        % Generate the OIs for one period of the rasterized stimulus
        theListOfStimulusRasterRetinalImages = computeRasterScanRetinalImagesForOneFullRefresh(theOI, theEscene0degs, ...
                pixelTimeOnSeconds, rasterLineDutyCycle, ...
                stimulusRefreshIntervalSeconds, simulationTimeStepSeconds, ...
                visualizeEachRetinalImage);
        fprintf('Done !\n');

        fprintf('\nSaving retinal images for the stimulus and the background raster. Please wait ...')
        save(retinalImagesMatFileName, ...
            'theConeMosaic', ...
            'stimulusRefreshIntervalSeconds', ...
            'simulationTimeStepSeconds', ...
            'theListOfBackgroundRasterRetinalImages', ...
            'theListOfStimulusRasterRetinalImages', ...
            '-v7.3', '-nocompression');
        fprintf('Done !\n');
    end


    if (recomputeConeExcitations)

        fprintf('\nLoading retinal images for the background raster. Please wait ...')
        load(retinalImagesMatFileName, ...
            'theConeMosaic', ...
            'simulationTimeStepSeconds', ...
            'theListOfBackgroundRasterRetinalImages');
        fprintf('Done !\n');

        fprintf('\nGenerating the background OI sequence. Please wait ...');
        theSceneTemporalSupportSeconds = (0:(numel(theListOfBackgroundRasterRetinalImages)-1)) * simulationTimeStepSeconds;
        
        theBackgroundOIsequence = oiArbitrarySequence(...
            theListOfBackgroundRasterRetinalImages, ...
            theSceneTemporalSupportSeconds);
        fprintf('Done ! \n');

        % Cone excitation response time series to the background raster
        fprintf('\nComputing cone excitations response to the background raster OI sequence. Please wait ... ');
        [backgroundConeExcitationResponseTimeSeries, ~, ~,~, backgroundConeExcitationsTimeAxis] = ...
            theConeMosaic.compute(theBackgroundOIsequence, ...
            'withFixationalEyeMovements', false);
        fprintf('Done !\n');

        % Free some memory
        clear 'theBackgroundOIsequence'
        clear 'theListOfBackgroundRasterRetinalImages'



        fprintf('\nLoading retinal images for the stimulus raster. Please wait ...')
        load(retinalImagesMatFileName, ...
            'stimulusRefreshIntervalSeconds', ...
            'simulationTimeStepSeconds', ...
            'theListOfStimulusRasterRetinalImages');
        fprintf('Done !\n');


        % Simulation duration
        simulationDurationSeconds = stimulusRefreshIntervalSeconds * nStimulusFramesPerTrial;

        % Generate  fixational eye movements for nTrials, each lasting for simulationDurationSeconds
        generateFixationalEyeMovements(theConeMosaic, ...
                simulationDurationSeconds, nTrials);


        fprintf('\nGenerating the (periodic) stimulus OI sequence. Please wait ...');
        totalSimulationTimeSteps = size(theConeMosaic.fixEMobj.emPosArcMin,2);
        theSceneTemporalSupportSeconds = (0:(totalSimulationTimeSteps-1)) * simulationTimeStepSeconds;

        theStimulusOIsequence = oiArbitrarySequence(...
            theListOfStimulusRasterRetinalImages, ...
            theSceneTemporalSupportSeconds, ...
            'isPeriodic', true);
        fprintf('Done ! \n');


        % Cone excitation response time series to the stimulus raster 
        fprintf('\nComputing cone excitations response to the stimulus raster OI sequence. Please wait ... ');
        [noiseFreeConeExcitationResponseTimeSeries, noisyConeExcitationResponseTimeSeries, ~,~, timeAxis] = ...
            theConeMosaic.compute(theStimulusOIsequence, ...
            'withFixationalEyeMovements', true);
        fprintf('Done !\n');

        fprintf('\nSaving cone excitation responses to %s ... ', coneExcitationsMatFileName);
        save(coneExcitationsMatFileName, ...
            'theConeMosaic', ...
            'backgroundConeExcitationResponseTimeSeries', ...
            'backgroundConeExcitationsTimeAxis', ...
            'noiseFreeConeExcitationResponseTimeSeries', ...
            'noisyConeExcitationResponseTimeSeries', ...
            'timeAxis', '-v7.3');
        fprintf('Done !\n');
    end

    if (visualizeStimulusAndConeExcitationSequence) && (~recomputePhotocurrents)
        visualizeRetinalImageAndConeExcitations(...
            retinalImagesMatFileName, ...
            coneExcitationsMatFileName, ...
            videoFilename);
    end


    if (recomputePhotocurrents)
        load(coneExcitationsMatFileName, ...
            'theConeMosaic', ...
            'noiseFreeConeExcitationResponseTimeSeries', ...
            'timeAxis');

        warmupPeriodSeconds = 0.0;
        backgroundConeExcitations = zeros(1, theConeMosaic.conesNum);

        [photocurrentResponseTimeSeries, photocurrentResponseTimeAxis] = ...
            computeFullPhotocurrentModelResponse(backgroundConeExcitations, noiseFreeConeExcitationResponseTimeSeries, ...
            theConeMosaic.integrationTime, warmupPeriodSeconds);

        save(photocurrentsMatFileName, ...
            'photocurrentResponseTimeSeries', ...
            'photocurrentResponseTimeAxis', ...
            '-v7.3');

        videoFilename = strrep(videoFilename, 'MosaicActivation', 'MosaicPhotocurrents');
        visualizeRetinalImageAndConePhotoCurrents(coneExcitationsMatFileName, photocurrentsMatFileName, videoFilename);
    end

end


function visualizeRetinalImageAndConePhotoCurrents(coneExcitationsMatFileName, photocurrentsMatFileName, videoFilename)

    fprintf('\nLoading data from %s. Please wait ...', coneExcitationsMatFileName);

    load(coneExcitationsMatFileName, ...
            'theConeMosaic', ...
            'theOIsequence', ...
            'timeAxis');
    fprintf('Done \n');


    fprintf('\nLoading data from %s. Please wait ...', photocurrentsMatFileName);

    load(photocurrentsMatFileName, ...
            'photocurrentResponseTimeSeries', ...
            'photocurrentResponseTimeAxis');
    fprintf('Done \n');



    % Determine the cone indices for which we will visualize their time series responses 
    [LconeIndicesVisualized, MconeIndicesVisualized, SconeIndicesVisualized] = ...
        determineVisualizedConeIndices(theConeMosaic);


    generateMosaicActivationVideo(theConeMosaic, theOIsequence, photocurrentResponseTimeSeries, ...
        LconeIndicesVisualized, MconeIndicesVisualized, SconeIndicesVisualized, ...
        timeAxis, videoFilename);
end


function visualizeRetinalImageAndConeExcitations(...
    retinalImagesMatFileName,coneExcitationsMatFileName, videoFilename)

    fprintf('\nLoading cone excitation data from %s. Please wait ...', coneExcitationsMatFileName);
    load(coneExcitationsMatFileName, ...
            'theConeMosaic', ...
            'backgroundConeExcitationResponseTimeSeries', ...
            'backgroundConeExcitationsTimeAxis');

    fprintf('\nLoading retinal image data from %s. Please wait ...', coneExcitationsMatFileName);
    load(retinalImagesMatFileName, ...
            'simulationTimeStepSeconds', ...
            'theListOfBackgroundRasterRetinalImages');


    % Determine the cone indices for which we will visualize their time series responses 
    [LconeIndicesVisualized, MconeIndicesVisualized, SconeIndicesVisualized] = ...
        determineVisualizedConeIndices(theConeMosaic);

    % Transform excitation counts to excitations rate
    dT = backgroundConeExcitationsTimeAxis(2)-backgroundConeExcitationsTimeAxis(1);
    backgroundConeExcitationResponseTimeSeries = backgroundConeExcitationResponseTimeSeries / dT;


    % First, visualize the response to the background raster
    fprintf('\nGenerating the background raster OI sequence. Please wait ...');
    theSceneTemporalSupportSeconds = (0:(numel(theListOfBackgroundRasterRetinalImages)-1)) * simulationTimeStepSeconds;
        
    theBackgroundOIsequence = oiArbitrarySequence(...
            theListOfBackgroundRasterRetinalImages, ...
            theSceneTemporalSupportSeconds);
    fprintf('Done ! \n');

    % Generate the video of the cone mosaic excitations response to the stimulus raster
    generateMosaicActivationVideo(theConeMosaic, theBackgroundOIsequence , backgroundConeExcitationResponseTimeSeries, ...
        LconeIndicesVisualized, MconeIndicesVisualized, SconeIndicesVisualized, ...
        backgroundConeExcitationsTimeAxis, sprintf('%s-BackgroundRaster', videoFilename));
    
    clear 'theBackgroundOIsequence'
    clear 'theListOfBackgroundRasterRetinalImages'


    fprintf('\nLoading cone excitation data from %s. Please wait ...', coneExcitationsMatFileName);
    load(coneExcitationsMatFileName, ...
            'noiseFreeConeExcitationResponseTimeSeries', ...
            'timeAxis');

    fprintf('\nLoading retinal image data from %s. Please wait ...', coneExcitationsMatFileName);
    load(retinalImagesMatFileName, ...
            'theListOfStimulusRasterRetinalImages');


    % Next, visualize the response to the stimulus raster
    fprintf('\nGenerating the (periodic) stimulus raster OI sequence. Please wait ...');
    totalSimulationTimeSteps = size(theConeMosaic.fixEMobj.emPosArcMin,2);
    theSceneTemporalSupportSeconds = (0:(totalSimulationTimeSteps-1)) * simulationTimeStepSeconds;
      
    theStimulusOIsequence = oiArbitrarySequence(...
            theListOfStimulusRasterRetinalImages, ...
            theSceneTemporalSupportSeconds, ...
            'isPeriodic', true);
    fprintf('Done ! \n');

    % Transform excitation counts to excitations rate
    dT = timeAxis(2)-timeAxis(1);
    noiseFreeConeExcitationResponseTimeSeries = noiseFreeConeExcitationResponseTimeSeries / dT;


    % Generate the video of the cone mosaic excitations response to the stimulus raster
    generateMosaicActivationVideo(theConeMosaic, theStimulusOIsequence, noiseFreeConeExcitationResponseTimeSeries, ...
        LconeIndicesVisualized, MconeIndicesVisualized, SconeIndicesVisualized, ...
        timeAxis, sprintf('%s-StimulusRaster', videoFilename));

end


function generateMosaicActivationVideo(theConeMosaic, theOIsequence, mosaicResponseTimeSeries, ...
    LconeIndicesVisualized, MconeIndicesVisualized, SconeIndicesVisualized, ...
    timeAxis, videoFilename)

    labeledConeIndices = [...
        LconeIndicesVisualized(:); ...
        MconeIndicesVisualized(:); ...
        SconeIndicesVisualized(:)];

    LconeIndicesResponses = mosaicResponseTimeSeries(:,:,LconeIndicesVisualized);
    MconeIndicesResponses = mosaicResponseTimeSeries(:,:,MconeIndicesVisualized);
    SconeIndicesResponses = mosaicResponseTimeSeries(:,:,SconeIndicesVisualized);
    

    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1504 768], 'Color', [1 1 1]);
    
    colsNum = 2; 
    rowsNum = 1;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
               'rowsNum', rowsNum, ...
               'colsNum', colsNum, ...
               'heightMargin',  0.01, ...
               'widthMargin',    0.04, ...
               'leftMargin',     0.04, ...
               'rightMargin',    0.00, ...
               'bottomMargin',   0.05, ...
               'topMargin',      0.02);

    ax1 = subplot('Position', subplotPosVectors(1,1).v);
    ax2 = subplot('Position', subplotPosVectors(1,2).v);

    ax3 = axes('Position', [0.575 0.12 0.4 0.22]);

    videoOBJ = VideoWriter(sprintf('%s_RetinalImageAndConeExcitations',videoFilename), 'MPEG-4');  % H264format (has artifacts)
    videoOBJ.FrameRate = 60;
    videoOBJ.Quality = 100;
    videoOBJ.open();

    nTrials = size(mosaicResponseTimeSeries,1);
    nTimePoints = size(mosaicResponseTimeSeries,2);

    visualizedConeAperture = 'lightCollectingArea5Sigma';
    backgroundColor = [0.2 0.2 0.2];

    activationRange = [min(mosaicResponseTimeSeries(:)) max(mosaicResponseTimeSeries(:))]; 

    theRetinalImage = theOIsequence.frameAtIndex(1);

    % Convert spatial support in degrees
    spatialSupportMM = oiGet(theRetinalImage, 'spatial support', 'mm');
    theOptics = oiGet(theRetinalImage, 'optics');
    focalLength = opticsGet(theOptics, 'focal length');
    mmPerDegree = focalLength*tand(1)*1e3;
    spatialSupportDegs = spatialSupportMM/mmPerDegree;
    spatialSupportX = theConeMosaic.eccentricityDegs(1) + spatialSupportDegs(1,:,1);
    spatialSupportY = theConeMosaic.eccentricityDegs(2) + spatialSupportDegs(:,1,2);

    visualizedFOV = max(spatialSupportX) - min(spatialSupportX);
    domainVisualizationLimits(1:2) = theConeMosaic.eccentricityDegs(1) + 0.5*visualizedFOV*[-1 1];
    domainVisualizationLimits(3:4) = theConeMosaic.eccentricityDegs(2) + 0.5*visualizedFOV*[-1 1];
    domainVisualizationTicks.x = -2:0.2:2;
    domainVisualizationTicks.y = -2:0.2:2;

    targetWavelength = 840;

    irradianceMapWatts480nmPerMMsquared = zeros(nTimePoints, numel(spatialSupportY), numel(spatialSupportX));
    for iTimePoint = 1:nTimePoints
        % Get the current retinal image
        theRetinalImage = theOIsequence.frameAtIndex(iTimePoint);
        irradianceMapWatts480nmPerMMsquared(iTimePoint,:,:) = computeIrradianceInWattsPerMM2(theRetinalImage, targetWavelength);
    end
    irradianceRange = [min(irradianceMapWatts480nmPerMMsquared(:)) max(irradianceMapWatts480nmPerMMsquared(:))];

    dT = timeAxis(2)-timeAxis(1);
    for iTrial = 1:nTrials
        for iTimePoint = 1:nTimePoints

            % Plot retinal irradiance at target wavelengh
            imagesc(ax1, spatialSupportX, spatialSupportY, squeeze(irradianceMapWatts480nmPerMMsquared(iTimePoint,:,:)));
            set(ax1, 'CLim', irradianceRange);
            axis(ax1,'image');
            set(ax1, 'XTick', domainVisualizationTicks.x, 'YTick', domainVisualizationTicks.y);
            set(ax1, 'XLim', domainVisualizationLimits(1:2), 'YLim', domainVisualizationLimits(3:4));

            colormap(ax1, 'gray');
            colorbar(ax1,'north', 'Color', [0.8 0.8 0.8], 'FontSize', 12, 'FontName', 'Spot mono');
    
            set(ax1, 'FontSize', 16);
            title(ax1, sprintf('irradiance, mWatts/mm^2 @ %dnm (simulated Tuten AOSLO, time step: %2.2f msec)',targetWavelength, dT*1000), ...
                'FontSize', 14);
            xlabel(ax1, 'eccentricity, x (degs)');
            ylabel(ax1, 'eccentricity, y (degs)');

            theConeMosaic.visualize('figureHandle', hFig, 'axesHandle', ax2, ...
                'visualizedConeAperture', visualizedConeAperture, ...
                'activation', mosaicResponseTimeSeries(iTrial,iTimePoint,:), ...
                'activationRange', activationRange, ...
                'currentEMposition', squeeze(theConeMosaic.fixEMobj.emPosArcMin(iTrial,iTimePoint,:))/60, ...
                'displayedEyeMovementData', struct('trial', iTrial, 'timePoints', 1:iTimePoint), ...
                'domainVisualizationLimits', domainVisualizationLimits, ...
                'domainVisualizationTicks', domainVisualizationTicks, ...
                'horizontalActivationColorBarInside', true, ...
                'colorbarFontSize', 10, ...
                'labelConesWithIndices', labeledConeIndices, ...
                'noYLabel', true, ...
                'backgroundColor', backgroundColor, ...
                'plotTitleFontSize', 14, ...
                'plotTitle',  sprintf('simulated cone mosaic excitation rate, R*/sec (trial no.%d, t = %2.2f ms)', iTrial, timeAxis(iTimePoint)*1000));
        

            scatter(ax3, timeAxis(1:iTimePoint)*1000, squeeze(LconeIndicesResponses(iTrial, 1:iTimePoint,:)), 50, ...
                'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0.5 0.5], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.4, 'LineWidth', 0.5);
            hold(ax3, 'on');
            scatter(ax3, timeAxis(1:iTimePoint)*1000, squeeze(MconeIndicesResponses(iTrial, 1:iTimePoint,:)), 30, ...
                'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', [0.5 0.8 0.5], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.4, 'LineWidth', 0.5);
            scatter(ax3, timeAxis(1:iTimePoint)*1000, squeeze(SconeIndicesResponses(iTrial, 1:iTimePoint,:)), 10, ...
                'MarkerFaceColor', [0 0 1],  'MarkerEdgeColor', [0.5 0.5 1], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.4, 'LineWidth', 0.5);

            hold(ax3, 'off');
            set(ax3, 'XLim', [timeAxis(1) timeAxis(end)]*1000, 'XTick', 0:3:200, 'YLim', activationRange);
            set(ax3, 'XColor', [0.75 0.75 0.75], 'YColor', [0.75 0.75 0.75], 'Color', 'none', 'FontSize', 14);
            grid(ax3, 'on');
            box(ax3, 'off')
            xlabel(ax3, 'time (ms)');
            ylabel(ax3, 'cone excitation rate (R*/sec)');

            drawnow;
            videoOBJ.writeVideo(getframe(hFig));
        end % for iTimePoint
    end % for iTrial

    videoOBJ.close();
end



function [LconeIndicesVisualized, MconeIndicesVisualized, SconeIndicesVisualized] = ...
        determineVisualizedConeIndices(theConeMosaic)
        
    slicePositionDegs = -0.08;
    sliceThicknessDegs = 0.05;
    [LconeIndicesVisualized, theROI] = ...
            theConeMosaic.coneIndicesOfCertainTypeWithinHorizontalSlice(...
            slicePositionDegs, sliceThicknessDegs, cMosaic.LCONE_ID);
    MconeIndicesVisualized = ...
            theConeMosaic.coneIndicesOfCertainTypeWithinHorizontalSlice(...
            slicePositionDegs, sliceThicknessDegs, cMosaic.MCONE_ID);

    SconeIndicesVisualized = ...
            theConeMosaic.coneIndicesOfCertainTypeWithinHorizontalSlice(...
            slicePositionDegs, sliceThicknessDegs, cMosaic.SCONE_ID);
end



function [theNullScene, theEscene0degs] = generateScenes(testESizeDeg)
    defaultParams = sceBerkeleyAOTumblingEscene();

    temporalModulationParams_xShiftPerFrame = [0 5/60 0];
    temporalModulationParams_yShiftPerFrame = [0 5/60 0];

    [sce0,sce90,sce180,sce270,sceBg,sceneParams] = t_BerkeleyAOTumblingESceneEngine( ...
        'temporalModulationParams_xShiftPerFrame', temporalModulationParams_xShiftPerFrame, ...
        'temporalModulationParams_yShiftPerFrame', temporalModulationParams_yShiftPerFrame);

    theEsceneSequence0degs = sce0.compute(testESizeDeg);
    theNullSequence = sceBg.compute(testESizeDeg);

    frameNo = 2;
    theEscene0degs = theEsceneSequence0degs{frameNo};
    theNullScene = theNullSequence{frameNo};
end

function theConeMosaic = generateConeMosaic(...
    mosaicEccDegs, mosaicSizeDegs, mosaicIntegrationTimeSeconds)
    
    % Set cone aperture parameters for the @cMosaic
    % Use a Gaussian cone aperture with
    % sigma equal to 0.204 x inner segment diameter (cone diameter)
    sigmaGaussian = 0.204;  % From McMahon et al, 2000
    coneApertureModifiers = struct(...
            'smoothLocalVariations', true, ...
            'sigma',  sigmaGaussian, ...
            'shape', 'Gaussian');

    % Instantiate a @sMosaic
    theConeMosaic = cMosaic(...
                'sourceLatticeSizeDegs', 64, ...
                'sizeDegs', mosaicSizeDegs, ...      
                'eccentricityDegs', mosaicEccDegs, ...   
                'eccVaryingConeBlur', true, ...
                'rodIntrusionAdjustedConeAperture', true, ...
                'coneApertureModifiers', coneApertureModifiers, ...
                'integrationTime', mosaicIntegrationTimeSeconds ...
                );
end

function generateFixationalEyeMovements(theConeMosaic, simulationDurationSeconds, nTrials)

    % Fixational eye movements
    eyeMovementsPerTrial = round(simulationDurationSeconds/theConeMosaic.integrationTime);
    
    fixEMobj = fixationalEM();
    fixEMobj.controlGamma = 0.18;
    fixEMobj.feedbackGain = 0.16;
    fixEMobj.microSaccadeType = 'none';    % No micro-saccades
    fixEMobj.randomSeed = 1;

    % Compute fEM paths
    fixEMobj.computeForCmosaic(...
        theConeMosaic, eyeMovementsPerTrial+10, ...  % add some extra time point, which we trim to what we need next
        'nTrials', nTrials, 'rSeed', 857);

    % Only keep what we need
    fixEMobj.trimPathToIncludeSelectTimePoints(1:eyeMovementsPerTrial);
    
    % Set the fMEs to the coneMosaic
    theConeMosaic.emSetFixationalEMObj(fixEMobj);
end

function theOI = generateOptics(simulateAdaptiveOpticsviewingConditions, theConeMosaic)

    if (simulateAdaptiveOpticsviewingConditions)
        testSubjectID = 0;
        subtractCentralRefraction = ~true;
        pupilDiamMM = 6.0;
        noLCA = true;
        wavefrontSpatialSamples = 501;
        psfUpsampleFactor = 2;
    else
        testSubjectID = 1;
        subtractCentralRefraction = true;
        pupilDiamMM = 3.0;
        noLCA = false;
        wavefrontSpatialSamples = 301;
        psfUpsampleFactor = 1;
    end

    [oiEnsemble, psfEnsemble] = theConeMosaic.oiEnsembleGenerate(theConeMosaic.eccentricityDegs, ...
			       'zernikeDataBase', 'Polans2015', ...
			       'subjectID', testSubjectID, ...
			       'pupilDiameterMM', pupilDiamMM , ...
			       'refractiveErrorDiopters', 0, ...
			       'noLCA', noLCA, ...
			       'zeroCenterPSF', true, ...
			       'subtractCentralRefraction', subtractCentralRefraction, ...
			       'wavefrontSpatialSamples', wavefrontSpatialSamples, ...
			       'upsampleFactor', psfUpsampleFactor, ...
			       'warningInsteadOfErrorForBadZernikeCoeffs', false);


    % Extract the optics
    theOI = oiEnsemble{1};
end


function [simulationTimeStepSeconds, stimulusRefreshIntervalSeconds] = ...
    computeSimulationTimeStep(pixelTimeOnSeconds, rasterLineDutyCycle, theTestScene, fractionLinesScannedPerSimulationTimeStep)

    theScenePhotons = sceneGet(theTestScene, 'photons');
    rasterLinesNum = size(theScenePhotons,1);
    pixelsPerRasterLine = size(theScenePhotons,2);

    % Active scan line duration
    rasterLineActiveDurationSeconds = pixelsPerRasterLine * pixelTimeOnSeconds;

    % Total duration of a scan (active + blanking time)
    rasterLineTotalDurationSeconds = rasterLineActiveDurationSeconds / rasterLineDutyCycle;

    % Make the simulation time step equal to the time required to raster
    % through the specified raster lines
    rasterLinesNumScannedPerSimulationTimeStep = round(rasterLinesNum * fractionLinesScannedPerSimulationTimeStep);

    if (mod(rasterLinesNum, rasterLinesNumScannedPerSimulationTimeStep) ~= 0)
        fprintf('With the passed fractionLinesScannedPerSimulationTimeStep (%f)\n, the number of raster lines in the stimulus (%d)\ndoes not divide fully with the number of raster lines scanned per simulation time step (%d)\n', ...
            fractionLinesScannedPerSimulationTimeStep, rasterLinesNum , rasterLinesNumScannedPerSimulationTimeStep);
        error('Specify a different fractionLinesScannedPerSimulationTimeStep')
    end
    simulationTimeStepSeconds = rasterLinesNumScannedPerSimulationTimeStep * rasterLineTotalDurationSeconds;

    % Compute stimulus refresh interval
    stimulusRefreshIntervalSeconds = rasterLinesNum * rasterLineTotalDurationSeconds;
end


function  irradianceEnergyMapWattsPerMMsquared = computeIrradianceInWattsPerMM2(theRetinalImage, targetWavelength)
    % Irradiance in photons
    irradiancePhotons = oiGet(theRetinalImage,'photons');

    wave = oiGet(theRetinalImage,'wave');
    idxTargetWavelength = find(wave == targetWavelength);

    % Irradiance in Watts / m^2 / nm
    irradianceEnergyMapWattsPerMeterSquared = Quanta2Energy(wave,irradiancePhotons);

    % Irradiance at target wavelength, in milliWatts per mm^2
    mmPerMeter = 1e3;
    mmSquaredPerMeterSquared = mmPerMeter^2;
    irradianceEnergyMapWattsPerMMsquared = 1/mmSquaredPerMeterSquared * squeeze(irradianceEnergyMapWattsPerMeterSquared(:,:,idxTargetWavelength));
end



function theListOfRetinalImages = computeRasterScanRetinalImagesForOneFullRefresh(theOI, theTestScene, ...
        pixelTimeOnSeconds, rasterLineDutyCycle,  ...
        stimulusRefreshIntervalSeconds, simulationTimeStepSeconds, ...
        visualizeEachRetinalImage)

    % Compute the retinal image of the non-rasterized scene
    theRetinalImage = oiCompute(theOI, theTestScene, 'pad value', 'mean');

    if (visualizeEachRetinalImage)
        % Compute the spatial support
        spatialSupportMM = oiGet(theRetinalImage, 'spatial support', 'mm');
        theOptics = oiGet(theRetinalImage, 'optics');
        focalLength = opticsGet(theOptics, 'focal length');
        mmPerDegree = focalLength*tand(1)*1e3;
        spatialSupportDegs = spatialSupportMM/mmPerDegree;
        spatialSupportXdegs = spatialSupportDegs(1,:,1);
        spatialSupportYdegs = spatialSupportDegs(:,1,2);
    end


    % Compute the power at 840 nm
    irradiancePowerWatts480nmPerMMsquared = computeIrradianceInWattsPerMM2(theRetinalImage, 840);
    nonRasterizedEscenePower = mean(irradiancePowerWatts480nmPerMMsquared(:));
    fprintf(2,'\nNon-rasterized scene power at 840 nm is: %E W/mm2.', nonRasterizedEscenePower);
    fprintf(2,'\nWill says it should be: 8.37E-04 W/mm2.\n');

    % Retrieve the photons map of the non-raster scene
    theScenePhotons = sceneGet(theTestScene, 'photons');
    rasterLinesNum = size(theScenePhotons,1);
    pixelsPerRasterLine = size(theScenePhotons,2);

    % Active scan line duration
    rasterLineActiveDurationSeconds = pixelsPerRasterLine * pixelTimeOnSeconds;

    % Total duration of a scan (active + blanking time)
    rasterLineTotalDurationSeconds = rasterLineActiveDurationSeconds / rasterLineDutyCycle;

    % Blanking interval
    blankingIntervalDurationSeconds = rasterLineTotalDurationSeconds - rasterLineActiveDurationSeconds;

    % Assuming we have equal left and right blanking intervals
    leftBlankingIntervalDurationSeconds = 0.5 * blankingIntervalDurationSeconds;

    % Compute pixel ON times
    pixelONtimes = zeros(1,rasterLinesNum*pixelsPerRasterLine);
    for iRasterLine = 1:rasterLinesNum
        activeScanONtime = leftBlankingIntervalDurationSeconds + (iRasterLine-1) * rasterLineTotalDurationSeconds;
        pixelIndices = (iRasterLine-1) * pixelsPerRasterLine + (1:pixelsPerRasterLine);
        pixelONtimes(pixelIndices) = activeScanONtime  + (1:numel(pixelIndices)) * pixelTimeOnSeconds;
    end

    % Checks
    assert(stimulusRefreshIntervalSeconds == rasterLinesNum * rasterLineTotalDurationSeconds, ...
        'something is wrong with the temporal properties of the stimulus');
    assert(rem(stimulusRefreshIntervalSeconds,simulationTimeStepSeconds) == 0, ...
        'the stimulus refresh interval is not an integer multiple of the simulationTimeStepSeconds');

    simulationTimeStepsNum = stimulusRefreshIntervalSeconds/simulationTimeStepSeconds;

    % The testScene was generated asumming that all pixels are ON during the stimulusRefreshIntervalSeconds
    % In the rasterized scene, a few raster lines are ON and only during the duration of the simulation step.
    % So we need to amplify their irradiance by the ratio of stimulusRefreshIntervalSeconds / simulationTimeStepSeconds
    radianceMultiplierToAccountForRasterSimulation = simulationTimeStepsNum;

    % Generate a scene for each simulation time step and its associated retinal image
    theListOfRetinalImages = cell(1, simulationTimeStepsNum);
    for iTimeStep = 1:simulationTimeStepsNum

        if (~visualizeEachRetinalImage)
            fprintf('Generating scene for simulation time step %d of %d\n', iTimeStep, simulationTimeStepsNum);
        end

        t1 = (iTimeStep-1)*simulationTimeStepSeconds;
        t2 = t1 + simulationTimeStepSeconds;


        % Find which pixels are ON between t1 and t2
        activePixelsIndices = find( (pixelONtimes >= t1) & (pixelONtimes < t2));
        
        % Find (rows,cols) of the active pixels 
        [activePixelCols, activePixelRows] = ind2sub([rasterLinesNum pixelsPerRasterLine], activePixelsIndices);

        % Zero photons everywhere for the scene corresponding to this iTimeStep
        theFrameScenePhotons = theScenePhotons*0;

        % Adjusted radiance for the rasterized stimulus
        theFrameScenePhotons(activePixelRows, activePixelCols,:) = ...
                theScenePhotons(activePixelRows, activePixelCols,:) * radianceMultiplierToAccountForRasterSimulation;
        
        theFrameScene = sceneSet(theTestScene, 'photons', theFrameScenePhotons);
        theRetinalImageAtThisSimulationTimeStep = oiCompute(theOI, theFrameScene, 'pad value', 'zero');

        if (visualizeEachRetinalImage)
            figure(33); clf;
            image(spatialSupportXdegs, spatialSupportYdegs, oiGet(theRetinalImageAtThisSimulationTimeStep, 'rgbimage'));
            axis 'image'
        end

        theListOfRetinalImages{iTimeStep} = theRetinalImageAtThisSimulationTimeStep;
    end

end



function [photocurrentResponseTimeSeries, photocurrentResponseTimeAxis, coneExcitationResponseTimeSeries] = ...
        computeConeMosaicPhotoCurrentsResponse(theConeMosaic, ...
            coneExcitationResponseBackground, ...
            coneExcitationResponseTimeSeries, ...
            warmupPeriodSeconds, ...
            coolDownPeriodSeconds)

    
    warmUpSamplesNum = warmupPeriodSeconds / theConeMosaic.integrationTime;
    coolDownSamplesNum = coolDownPeriodSeconds / theConeMosaic.integrationTime;

    % Concatenate (pre-stimulus) warm up background
    backgroundConeExcitationsRepMat = repmat(coneExcitationResponseBackground, [size(coneExcitationResponseTimeSeries,1) warmUpSamplesNum 1 ]);
    coneExcitationResponseTimeSeries = cat(2, backgroundConeExcitationsRepMat, coneExcitationResponseTimeSeries);
    size(coneExcitationResponseTimeSeries)
    % Concatenate (post-stimulus) cool down background
    backgroundConeExcitationsRepMat = repmat(coneExcitationResponseBackground, [size(coneExcitationResponseTimeSeries,1) coolDownSamplesNum 1 ]);
    coneExcitationResponseTimeSeries = cat(2, coneExcitationResponseTimeSeries, backgroundConeExcitationsRepMat);


    [photocurrentResponseTimeSeries, photocurrentResponseTimeAxis] = ...
         computeFullPhotocurrentModelResponse(coneExcitationResponseBackground, coneExcitationResponseTimeSeries, ...
         theConeMosaic.integrationTime, warmupPeriodSeconds);
end




function [photocurrentResponseTimeSeries, photocurrentResponseTimeAxis] = ...
    computeFullPhotocurrentModelResponse(backgroundConeExcitations, stimulusConeExcitationResponseTimeSeries, ...
    coneMosaicIntegrationTime, warmupPeriodSeconds)

    osTimeStep = 1e-4;
    integerMultiplier = round(coneMosaicIntegrationTime/osTimeStep);
    osTimeStep = coneMosaicIntegrationTime/integerMultiplier;

    assert(rem(coneMosaicIntegrationTime, osTimeStep) == 0, 'coneMosaic.intergrationTime must be an integral multiple of os.time step');

    upSampleFactor = floor(coneMosaicIntegrationTime / osTimeStep);
    upSampleGain = 1/upSampleFactor;

    nTimeBins = size(stimulusConeExcitationResponseTimeSeries,2);
    theOSphotocurrentResponseTimeBinsNum = nTimeBins*upSampleFactor-1;
    

    nTrials = size(stimulusConeExcitationResponseTimeSeries,1);
    nCones = size(stimulusConeExcitationResponseTimeSeries,3);

    tIn = linspace(0,1,nTimeBins);
    tOut = linspace(0,1,theOSphotocurrentResponseTimeBinsNum);
    photocurrentResponseTimeSeries = zeros(nTrials, nTimeBins, nCones);
    photocurrentResponseTimeAxis = (0:(nTimeBins-1)) * coneMosaicIntegrationTime;

    theOSphotoCurrentResponseTimeAxis = (0:(theOSphotocurrentResponseTimeBinsNum-1))*osTimeStep;


    parfor iCone = 1:nCones

        fprintf('Computing pCurrent response for cone %d of %d\n', iCone, nCones);
        os = osBioPhys('eccentricity',0);
        os.timeStep = osTimeStep;
        os.set('noise flag', 'none');

        theBackgroundExcitationRate = backgroundConeExcitations(iCone)/coneMosaicIntegrationTime;
        theBackgroundExcitationsOStimeFrame = theBackgroundExcitationRate * os.timeStep;

        theState = os.osAdaptSteadyState(theBackgroundExcitationsOStimeFrame, [1 1]);
        theState.timeStep = os.timeStep;
        os.setModelState(theState);

        for iTrial = 1:nTrials
            singleConeSingleTrialConeExcitationsCountResponse = ...
                squeeze(stimulusConeExcitationResponseTimeSeries(iTrial, :, iCone));

            % upsample counts to os.timeStep
            singleConeSingleTrialConeExcitationsCountResponseOStimeSupport = upSampleGain * ...
                interp1(tIn, singleConeSingleTrialConeExcitationsCountResponse, tOut, 'nearest');

            % Convert to rate
            singleConeSingleTrialConeExcitationsRate = ...
                singleConeSingleTrialConeExcitationsCountResponseOStimeSupport / os.timeStep;

            % Compute full photocurrent model
            theOSphotoCurrentResponse = reshape(...
                os.osAdaptTemporal(reshape(singleConeSingleTrialConeExcitationsRate, [1 1 numel(singleConeSingleTrialConeExcitationsRate)])), ...
                [1 theOSphotocurrentResponseTimeBinsNum 1]);
       
            photocurrentResponseTimeSeries(iTrial,:,iCone) = interp1(theOSphotoCurrentResponseTimeAxis, theOSphotoCurrentResponse, photocurrentResponseTimeAxis, 'linear');

        end % iTrial
    end % iCone

    % Only keep response at t >= 0;
    if (warmupPeriodSeconds > 0)
        photocurrentResponseTimeAxis = photocurrentResponseTimeAxis - warmupPeriodSeconds;
        idx = find(photocurrentResponseTimeAxis>=0);
    
        % The delta-photocurrent (pCurrent - background)
        baselinePhotoCurrents = photocurrentResponseTimeSeries(1,idx(1)-1,:);
        photocurrentResponseTimeSeries = bsxfun(@minus, photocurrentResponseTimeSeries, baselinePhotoCurrents);

        photocurrentResponseTimeAxis = photocurrentResponseTimeAxis(idx);
        photocurrentResponseTimeSeries = photocurrentResponseTimeSeries(:,idx,:);
    end

end
