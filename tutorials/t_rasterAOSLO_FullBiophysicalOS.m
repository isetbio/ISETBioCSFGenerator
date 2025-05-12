function t_rasterAOSLO_FullBiophysicalOS

    % Simulation parameters.
    simulateStimulusRaster = true;
     
    % Simulation time step: here we pick a simulation time step during
    % which the AOSLO will scan through 8 of the 512 raster lines. This
    % computes to around 0.51 milliseconds
    fractionLinesScannedPerSimulationTimeStep = 8/512;

    % 2 raster lines
    %fractionLinesScannedPerSimulationTimeStep = 2/512;

    % How many frames / trial. 
    nStimulusFramesPerTrial = 3;

    % How many trials, which also means how many fEMs
    nTrials = 5;

    testIncrementDecrementScenes = true;

    observerFixationalEyeMovementCharacteristics = 'fast';
    observerFixationalEyeMovementCharacteristics = 'slow';
    observerFixationalEyeMovementCharacteristics = 'default';

    % Compute cone mosaic and retinal images of stimulus and background
    recomputeRetinalImages = ~true;

    cropRetinalImagesForConeMosaic = true;
    visualizeTheSceneRadiance = ~true;
    
    % Compute cone excitations response
    recomputeConeExcitations = ~true;

    % Visualize the stimulus and the cone excitations response
    visualizeStimulusAndConeExcitationSequence = ~true;

    % Compute photocurrent response
    recomputePhotocurrents = ~true;
    subtractBackgroundPhotoCurrents = false;

    % Visualize the stimulus and the photocurrents response
    visualizeStimulusAndPhotocurrentSequence = true;


    NDfilterDensity = 1.7;
    %NDfilterDensity = 1.0;
    %NDfilterDensity = 0.5;
    %NDfilterDensity = 0.0;
    pCurrentVisualizedRange = 30;

    mosaicHorizontalEccentricity = -1;
    

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
        retinalImagesMatFileName = fullfile(simFileBase, sprintf('eccDegs_%1.1f_ND_%2.2f_coneMosaicAndRetinalImages.mat', mosaicHorizontalEccentricity, NDfilterDensity));
        coneExcitationsMatFileName = fullfile(simFileBase, sprintf('eccDegs_%1.1f_ND_%2.2f_excitations.mat', mosaicHorizontalEccentricity,NDfilterDensity));
        photocurrentsMatFileName = fullfile(simFileBase, sprintf('eccDegs_%1.1f_ND_%2.2f_photocurrents.mat', mosaicHorizontalEccentricity,NDfilterDensity));
        videoFilename = fullfile(figureFileBase, sprintf('eccDegs_%1.1f_ND_%2.2f_activation', mosaicHorizontalEccentricity,NDfilterDensity));
    else
        retinalImagesMatFileName = fullfile(simFileBase, sprintf('noRaster_eccDegs_%1.1f_ND_%2.2f_coneMosaicAndRetinalImages.mat', mosaicHorizontalEccentricity, NDfilterDensity));
        coneExcitationsMatFileName = fullfile(simFileBase, sprintf('noRaster_eccDegs_%1.1f_ND_%2.2f_excitations.mat', mosaicHorizontalEccentricity,NDfilterDensity));
        photocurrentsMatFileName = fullfile(simFileBase, sprintf('noRaster_eccDegs_%1.1f_ND_%2.2f_photocurrents.mat', mosaicHorizontalEccentricity,NDfilterDensity));
        videoFilename = fullfile(figureFileBase, sprintf('noRaster_eccDegs_%1.1f_ND_%2.2f_activation', mosaicHorizontalEccentricity,NDfilterDensity));
    end

    if (~strcmp(observerFixationalEyeMovementCharacteristics, 'default'))
        postFixString = sprintf('_%sFEM.mat',observerFixationalEyeMovementCharacteristics);
        retinalImagesMatFileName = strrep(retinalImagesMatFileName, '.mat', postFixString);
        coneExcitationsMatFileName = strrep(coneExcitationsMatFileName, '.mat', postFixString);
        photocurrentsMatFileName = strrep(photocurrentsMatFileName, '.mat', postFixString);
        videoFilename = strep(videoFilename, 'activation', sprintf('activation_%sFEM.mat',observerFixationalEyeMovementCharacteristics));
    end

    
    if (recomputeRetinalImages)
        % Generate scenes
        testESizeDeg = 10/60;

        if (testIncrementDecrementScenes)
            defaultParams = sceBerkeleyAOTumblingEscene();
            % Change the second primary to 680 from 650
            AOPrimaryWls = defaultParams.AOPrimaryWls;
            AOAOCornealPowersUW = defaultParams.AOAOCornealPowersUW;
    
            inFocusWavelength = 680;
            AOPrimaryWls(1) = 840; AOAOCornealPowersUW(1) = 170;     % Imaging channel
            AOPrimaryWls(2) = 680; AOAOCornealPowersUW(2) = 170 * 10^-NDfilterDensity;  % Stimulus modulation channel 
            AOAOCornealPowersUW(3) = 170;

            [theNullScene, theIncrementsEscene0degs, theDecrementsEscene0degs, thePresentationDisplay, sceneParams] = ...
                generateIncrementDecrementScenes(...
                testESizeDeg, AOPrimaryWls, AOAOCornealPowersUW, visualizeTheSceneRadiance, figureFileBase);

        else
            inFocusWavelength = 840;
            [theNullScene, theDecrementsEscene0degs, thePresentationDisplay, sceneParams] = generateScenes(testESizeDeg, visualizeTheSceneRadiance, figureFileBase);
            theIncrementsEscene0degs = [];
        end


        % These are fixed params
        pixelTimeOnSeconds = 50 * 1e-9;
        rasterLineDutyCycle = 0.4;

        [simulationTimeStepSeconds, stimulusRefreshIntervalSeconds] = ...
            computeSimulationTimeStep(pixelTimeOnSeconds, rasterLineDutyCycle, ...
                theDecrementsEscene0degs, fractionLinesScannedPerSimulationTimeStep);


        % Set the mosaic integration time to the simulation time step
        mosaicIntegrationTimeSeconds = simulationTimeStepSeconds;
        
        % Mosaic size and eccentricity
        mosaicSizeDegs = 0.5*[1.0 1.0];
        mosaicEccDegs = [mosaicHorizontalEccentricity 0];

        if (cropRetinalImagesForConeMosaic)
            % Add some extra real estate for fMEs
            cropRetinalImagesForConeMosaicSize = mosaicSizeDegs + [0.1 0.1];
        else
            cropRetinalImagesForConeMosaicSize = [];
        end

        % Generate the cone mosaic
        theConeMosaic = generateConeMosaic(mosaicEccDegs, mosaicSizeDegs, ...
                mosaicIntegrationTimeSeconds, sceneParams.spectralSupport);


        % Generate optics
        simulateAdaptiveOpticsviewingConditions = true;
        theOI = generateOptics(simulateAdaptiveOpticsviewingConditions, theConeMosaic, inFocusWavelength);
    
        % Whether to visualize the retinal images of the rasterized stimulus
        visualizeEachRetinalImage = true;

        % Generate the OIs for one period of the rasterized background
        fprintf('\nGenerating retinal images for the background raster. Please wait ...')
        theListOfBackgroundRasterRetinalImages = computeRasterScanRetinalImagesForOneFullRefresh(theOI, ...
                theNullScene, pixelTimeOnSeconds, rasterLineDutyCycle, ...
                stimulusRefreshIntervalSeconds, simulationTimeStepSeconds, ...
                cropRetinalImagesForConeMosaicSize, ...
                visualizeEachRetinalImage, thePresentationDisplay, testIncrementDecrementScenes);
        fprintf('Done !\n');

        if (testIncrementDecrementScenes)
            fprintf('\nGenerating retinal images for the stimulus raster (E-decrements). Please wait ...')
            % Generate the OIs for one period of the rasterized stimulus
            theListOfStimulusDecrementsRasterRetinalImages = computeRasterScanRetinalImagesForOneFullRefresh(theOI, ...
                theDecrementsEscene0degs, pixelTimeOnSeconds, rasterLineDutyCycle, ...
                stimulusRefreshIntervalSeconds, simulationTimeStepSeconds, ...
                cropRetinalImagesForConeMosaicSize, ...
                visualizeEachRetinalImage, thePresentationDisplay, testIncrementDecrementScenes);
            fprintf('Done !\n');

            fprintf('\nGenerating retinal images for the stimulus raster (E-increments). Please wait ...')
            % Generate the OIs for one period of the rasterized stimulus
            theListOfStimulusIncrementsRasterRetinalImages = computeRasterScanRetinalImagesForOneFullRefresh(theOI, ...
                theIncrementsEscene0degs, pixelTimeOnSeconds, rasterLineDutyCycle, ...
                stimulusRefreshIntervalSeconds, simulationTimeStepSeconds, ...
                cropRetinalImagesForConeMosaicSize, ...
                visualizeEachRetinalImage, thePresentationDisplay, testIncrementDecrementScenes);
            fprintf('Done !\n');

            fprintf('\nSaving retinal images for the increments and decrements stimuli and the background raster. Please wait ...')
            save(retinalImagesMatFileName, ...
                'theConeMosaic', ...
                'stimulusRefreshIntervalSeconds', ...
                'simulationTimeStepSeconds', ...
                'thePresentationDisplay', ...
                'theListOfBackgroundRasterRetinalImages', ...
                'theListOfStimulusDecrementsRasterRetinalImages', ...
                'theListOfStimulusIncrementsRasterRetinalImages', ...
                '-v7.3');
            fprintf('Done !\n');

        else
            fprintf('\nGenerating retinal images for the stimulus raster (default E scene). Please wait ...')
            % Generate the OIs for one period of the rasterized stimulus
            theListOfStimulusRasterRetinalImages = computeRasterScanRetinalImagesForOneFullRefresh(theOI, ...
                theDecrementsEscene0degs, pixelTimeOnSeconds, rasterLineDutyCycle, ...
                stimulusRefreshIntervalSeconds, simulationTimeStepSeconds, ...
                cropRetinalImagesForConeMosaicSize, ...
                visualizeEachRetinalImage, thePresentationDisplay, testIncrementDecrementScenes);
            fprintf('Done !\n');


            fprintf('\nSaving retinal images for the stimulus and the background raster. Please wait ...')
            save(retinalImagesMatFileName, ...
                'theConeMosaic', ...
                'stimulusRefreshIntervalSeconds', ...
                'simulationTimeStepSeconds', ...
                'thePresentationDisplay', ...
                'theListOfBackgroundRasterRetinalImages', ...
                'theListOfStimulusRasterRetinalImages', ...
                '-v7.3');
            fprintf('Done !\n');
        end

    end


    if (recomputeConeExcitations)

        fprintf('\nLoading retinal images for the background raster. Please wait ...')
        load(retinalImagesMatFileName, ...
            'theConeMosaic', ...
            'simulationTimeStepSeconds', ...
            'stimulusRefreshIntervalSeconds', ...
            'theListOfBackgroundRasterRetinalImages');
        fprintf('Done !\n');

        fprintf('\nGenerating the background OI sequence. Please wait ...');
        theSceneTemporalSupportSeconds = (0:(numel(theListOfBackgroundRasterRetinalImages)-1)) * simulationTimeStepSeconds;
        
        % Generate OI sequence for the background raster sequence
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


        % Simulation duration
        simulationDurationSeconds = stimulusRefreshIntervalSeconds * nStimulusFramesPerTrial;

        % Generate  fixational eye movements for nTrials, each lasting for simulationDurationSeconds
        generateFixationalEyeMovements(theConeMosaic, ...
                simulationDurationSeconds, nTrials, ...
                observerFixationalEyeMovementCharacteristics);


        if (testIncrementDecrementScenes)
            % Load the Decrements retinal images
            load(retinalImagesMatFileName, ...
                'simulationTimeStepSeconds', ...
                'theListOfStimulusDecrementsRasterRetinalImages');

            % Generate OI sequence for the stimulus raster (periodic)
            fprintf('\nGenerating the (periodic) stimulus OI sequence. Please wait ...');
            totalSimulationTimeSteps = size(theConeMosaic.fixEMobj.emPosArcMin,2);
            theSceneTemporalSupportSeconds = (0:(totalSimulationTimeSteps-1)) * simulationTimeStepSeconds;

            fprintf('\nGenerating the (periodic) stimulus OI sequence for the DECREMENTS stimulus. Please wait ...');
            theStimulusOIsequence = oiArbitrarySequence(...
                theListOfStimulusDecrementsRasterRetinalImages, ...
                theSceneTemporalSupportSeconds, ...
                'isPeriodic', true);
            fprintf('Done ! \n');

            clear 'theListOfStimulusDecrementsRasterRetinalImages';

            % Cone excitation response time series to the stimulus raster 
            fprintf('\nComputing cone excitations response to the DECREMENTS stimulus raster OI sequence. Please wait ... ');
            [noiseFreeConeExcitationDecrementsResponseTimeSeries, noisyConeExcitationDecrementsResponseTimeSeries, ~,~, timeAxis] = ...
                theConeMosaic.compute(theStimulusOIsequence, ...
                'withFixationalEyeMovements', true);
            fprintf('Done !\n');

            fprintf('\nSaving cone excitation DECREMENTS responses to %s ... ', coneExcitationsMatFileName);
            save(coneExcitationsMatFileName, ...
                'theConeMosaic', ...
                'backgroundConeExcitationResponseTimeSeries', ...
                'backgroundConeExcitationsTimeAxis', ...
                'noiseFreeConeExcitationDecrementsResponseTimeSeries', ...
                'noisyConeExcitationDecrementsResponseTimeSeries', ...
                'timeAxis', '-v7.3');

            % Repeat for the increments
            % Clear some RAM
            clear 'theStimulusOIsequence'
            clear 'noiseFreeConeExcitationDecrementsResponseTimeSeries'
            clear 'noisyConeExcitationDecrementsResponseTimeSeries'

            % Load the Decrements retinal images
            load(retinalImagesMatFileName, ...
                'theListOfStimulusIncrementsRasterRetinalImages');


            fprintf('\nGenerating the (periodic) stimulus OI sequence for the INCREMENTS stimulus. Please wait ...');
            theStimulusOIsequence = oiArbitrarySequence(...
                theListOfStimulusIncrementsRasterRetinalImages, ...
                theSceneTemporalSupportSeconds, ...
                'isPeriodic', true);
            fprintf('Done ! \n');

            clear 'theListOfStimulusIncrementsRasterRetinalImages';

            % Cone excitation response time series to the stimulus raster 
            fprintf('\nComputing cone excitations response to the INCREMENTS stimulus raster OI sequence. Please wait ... ');
            [noiseFreeConeExcitationIncrementsResponseTimeSeries, noisyConeExcitationIncrementsResponseTimeSeries, ~,~, timeAxis] = ...
                theConeMosaic.compute(theStimulusOIsequence, ...
                'withFixationalEyeMovements', true);
            fprintf('Done !\n');

            fprintf('\nAppending cone excitation INCREMENTS responses to %s ... ', coneExcitationsMatFileName);
            save(coneExcitationsMatFileName, ...
                'noiseFreeConeExcitationIncrementsResponseTimeSeries', ...
                'noisyConeExcitationIncrementsResponseTimeSeries', ...
                '-append');

        else
            load(retinalImagesMatFileName, ...
                'simulationTimeStepSeconds', ...
                'theListOfStimulusRasterRetinalImages');

            % Generate OI sequence for the stimulus raster (periodic)
            fprintf('\nGenerating the (periodic) stimulus OI sequence. Please wait ...');
            totalSimulationTimeSteps = size(theConeMosaic.fixEMobj.emPosArcMin,2);
            theSceneTemporalSupportSeconds = (0:(totalSimulationTimeSteps-1)) * simulationTimeStepSeconds;

            fprintf('\nGenerating the (periodic) stimulus OI sequence. Please wait ...');
            theStimulusOIsequence = oiArbitrarySequence(...
                theListOfStimulusRasterRetinalImages, ...
                theSceneTemporalSupportSeconds, ...
                'isPeriodic', true);
            fprintf('Done ! \n');
        
            clear 'theListOfStimulusRasterRetinalImages';

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
    end

    
    if (visualizeStimulusAndConeExcitationSequence)
        visualizeRetinalImageAndConeExcitations(...
            retinalImagesMatFileName, ...
            coneExcitationsMatFileName, ...
            testIncrementDecrementScenes, ...
            strrep(videoFilename, 'activation', 'excitations'));
    end


    if (recomputePhotocurrents)

        load(retinalImagesMatFileName, 'stimulusRefreshIntervalSeconds');

        warmupPeriodSeconds = 1.0;
        warmupPeriodStimulusIntervalsNum = ceil(warmupPeriodSeconds /stimulusRefreshIntervalSeconds);
        coolDownPeriodStimulusIntervalsNum = 4;
        warmUpPeriodSecondsToIncludeInVisualization = stimulusRefreshIntervalSeconds;
        debugWarmUpTime = false;

        if (testIncrementDecrementScenes)
            load(coneExcitationsMatFileName, ...
                'theConeMosaic', ...
                'backgroundConeExcitationResponseTimeSeries', ...
                'backgroundConeExcitationsTimeAxis', ...
                'noiseFreeConeExcitationDecrementsResponseTimeSeries', ...
                'timeAxis');  

            [photocurrentResponseDecrementsTimeSeries, photocurrentResponseTimeAxis, photocurrentResponseDecrementsTimeSeriesNoisy] = ...
                computeConeMosaicPhotoCurrentsResponse(...
                    theConeMosaic, ...
                    backgroundConeExcitationResponseTimeSeries, ...
                    noiseFreeConeExcitationDecrementsResponseTimeSeries, ...
                    warmupPeriodStimulusIntervalsNum, warmUpPeriodSecondsToIncludeInVisualization, ...
                    coolDownPeriodStimulusIntervalsNum, ...
                    subtractBackgroundPhotoCurrents, debugWarmUpTime);

             fprintf('\nSaving photocurrent DECREMENTS responses to %s ... ', photocurrentsMatFileName);
             save(photocurrentsMatFileName, ...
                'theConeMosaic', ...
                'photocurrentResponseDecrementsTimeSeries', ...
                'photocurrentResponseDecrementsTimeSeriesNoisy', ...
                'photocurrentResponseTimeAxis', ...
                '-v7.3');

             % Clear some RAM
             clear 'noiseFreeConeExcitationDecrementsResponseTimeSeries';
             clear 'noisyConeExcitationDecrementsResponseTimeSeries';
             clear 'photocurrentResponseDecrementsTimeSeries';
             clear 'photocurrentResponseDecrementsTimeSeriesNoisy'

             % Repeat for increments
             load(coneExcitationsMatFileName, ...
                'noiseFreeConeExcitationIncrementsResponseTimeSeries');

            [photocurrentResponseIncrementsTimeSeries, photocurrentResponseTimeAxis, photocurrentResponseIncrementsTimeSeriesNoisy] = ...
                computeConeMosaicPhotoCurrentsResponse(...
                    theConeMosaic, ...
                    backgroundConeExcitationResponseTimeSeries, ...
                    noiseFreeConeExcitationIncrementsResponseTimeSeries, ...
                    warmupPeriodStimulusIntervalsNum, warmUpPeriodSecondsToIncludeInVisualization, ...
                    coolDownPeriodStimulusIntervalsNum, ...
                    subtractBackgroundPhotoCurrents, debugWarmUpTime);

            fprintf('\nAppending photocurrent INCREMENTS responses to %s ... ', photocurrentsMatFileName);
            save(photocurrentsMatFileName, ...
                'photocurrentResponseIncrementsTimeSeries', ...
                'photocurrentResponseIncrementsTimeSeriesNoisy', ...
                '-append');

        else

            load(coneExcitationsMatFileName, ...
                'theConeMosaic', ...
                'backgroundConeExcitationResponseTimeSeries', ...
                'backgroundConeExcitationsTimeAxis', ...
                'noiseFreeConeExcitationResponseTimeSeries', ...
                'noisyConeExcitationResponseTimeSeries', ...
                'timeAxis');

             load(retinalImagesMatFileName, 'stimulusRefreshIntervalSeconds');

            [photocurrentResponseTimeSeries, photocurrentResponseTimeAxis, photocurrentResponseTimeSeriesNoisy] = ...
                computeConeMosaicPhotoCurrentsResponse(...
                    theConeMosaic, ...
                    backgroundConeExcitationResponseTimeSeries, ...
                    noiseFreeConeExcitationResponseTimeSeries, ...
                    warmupPeriodStimulusIntervalsNum, warmUpPeriodSecondsToIncludeInVisualization, ...
                    coolDownPeriodStimulusIntervalsNum, ...
                    subtractBackgroundPhotoCurrents, debugWarmUpTime);
    
            save(photocurrentsMatFileName, ...
                'theConeMosaic', ...
                'photocurrentResponseTimeSeries', ...
                'photocurrentResponseTimeSeriesNoisy', ...
                'photocurrentResponseTimeAxis', ...
                '-v7.3');
        end

    end % recompute photocurrents


    % Visualize photocurrents
    if (visualizeStimulusAndPhotocurrentSequence)
        % Visualize photocurrents response
        visualizeRetinalImageAndConePhotoCurrents(...
                maxVideoTrials, ...
                nStimulusFramesPerTrial, ...
                testIncrementDecrementScenes, ...
                pCurrentVisualizedRange, ...
                subtractBackgroundPhotoCurrents, ...
                retinalImagesMatFileName, ...
                photocurrentsMatFileName, ...
                strrep(videoFilename, 'activation', 'photocurrents'));
    end
end


function visualizeRetinalImageAndConePhotoCurrents(maxVideoTrials, nStimulusFramesPerTrial, ...
    testIncrementDecrementScenes, pCurrentVisualizedRange, subtractBackgroundPhotoCurrents, ...
    retinalImagesMatFileName, photocurrentsMatFileName, videoFilename)

    fprintf('\nLoading photocurrent response data from %s. Please wait ...', photocurrentsMatFileName);
    if (testIncrementDecrementScenes)
        load(photocurrentsMatFileName, ...
            'theConeMosaic', ...
            'photocurrentResponseDecrementsTimeSeries', ...
            'photocurrentResponseDecrementsTimeSeriesNoisy', ...
            'photocurrentResponseIncrementsTimeSeries', ...
            'photocurrentResponseIncrementsTimeSeriesNoisy', ...
            'photocurrentResponseTimeAxis');

        fprintf('\nLoading background raster images. Please wait...')
        load(retinalImagesMatFileName,'theListOfBackgroundRasterRetinalImages');
        fprintf('Done\n');
        
        % Determine the cone indices for which we will visualize their time series responses 
        [LconeIndicesVisualized, MconeIndicesVisualized, SconeIndicesVisualized] = ...
            determineVisualizedConeIndices(theConeMosaic);
    
        irradianceAtTargetWavelengthInsteadOfRGBimage = true;
        displayEyeMovements = true;
        
        simulationTimeStepSeconds = photocurrentResponseTimeAxis(2)-photocurrentResponseTimeAxis(1);
        theSceneTemporalSupportSeconds = (0:(numel(theListOfBackgroundRasterRetinalImages)-1)) * simulationTimeStepSeconds;
            
        theRetinalImage = theListOfBackgroundRasterRetinalImages{1};
        spatialSupportMM = oiGet(theRetinalImage, 'spatial support', 'mm');
        theOptics = oiGet(theRetinalImage, 'optics');
        focalLength = opticsGet(theOptics, 'focal length');
        mmPerDegree = focalLength*tand(1)*1e3;
        spatialSupportDegs = spatialSupportMM/mmPerDegree;
        spatialSupportX = theConeMosaic.eccentricityDegs(1) + spatialSupportDegs(1,:,1);
        spatialSupportY = theConeMosaic.eccentricityDegs(2) + spatialSupportDegs(:,1,2);

        if (~subtractBackgroundPhotoCurrents)
            pCurrentVisualizedRange = [-95 -25];
        end

        targetWavelength = 680;

        nTimePoints = numel(theListOfBackgroundRasterRetinalImages);
       
        theListOfBackgroundRGBImages = zeros(nTimePoints, numel(spatialSupportY), numel(spatialSupportX),3, 'single');
        theListOfStimulusRGBImages = zeros(nTimePoints, numel(spatialSupportY), numel(spatialSupportX), 3, 'single');

        theListOfBackgroundRetinalIrradianceMaps = zeros(nTimePoints, numel(spatialSupportY), numel(spatialSupportX), 'single');
        theListOfStimulusRetinalIrradianceMaps = zeros(nTimePoints, numel(spatialSupportY), numel(spatialSupportX), 'single'); 
        
        
        for iTimePoint = 1:nTimePoints
            % Get the current retinal image
            theRetinalImage = theListOfBackgroundRasterRetinalImages{iTimePoint};
            theListOfBackgroundRGBImages(iTimePoint,:,:,:) = single(oiGet(theRetinalImage, 'rgbimage'));
            theListOfBackgroundRetinalIrradianceMaps(iTimePoint,:,:) = single(computeIrradianceInWattsPerMM2(theRetinalImage, targetWavelength));
        end
        clear 'theListOfBackgroundRasterRetinalImages';

        fprintf('\nLoading stimulus (increments) raster images. Please wait...')
        load(retinalImagesMatFileName,'theListOfStimulusIncrementsRasterRetinalImages');
        for iTimePoint = 1:nTimePoints
            % Get the current retinal image
            theRetinalImage = theListOfStimulusIncrementsRasterRetinalImages{iTimePoint};
            theListOfStimulusRGBImages(iTimePoint,:,:,:) = single(oiGet(theRetinalImage, 'rgbimage'));
            theListOfStimulusRetinalIrradianceMaps(iTimePoint,:,:) = single(computeIrradianceInWattsPerMM2(theRetinalImage, targetWavelength));
        end
        clear 'theListOfStimulusIncrementsRasterRetinalImages';

        % Generate the video of the cone mosaic photocurrent response to the INCREMENTS stimulus raster
        generateMosaicActivationVideo(theConeMosaic, ...
            maxVideoTrials, ...
            nStimulusFramesPerTrial, ...
            theSceneTemporalSupportSeconds, ...
            spatialSupportX, spatialSupportY, ...
            theListOfStimulusRetinalIrradianceMaps, ...
            theListOfBackgroundRetinalIrradianceMaps, ...
            theListOfStimulusRGBImages, ...
            theListOfBackgroundRGBImages, ...
            photocurrentResponseIncrementsTimeSeries, ...
            photocurrentResponseIncrementsTimeSeriesNoisy, ...
            LconeIndicesVisualized, MconeIndicesVisualized, SconeIndicesVisualized, ...
            displayEyeMovements, ...
            photocurrentResponseTimeAxis, ...
            'photocurrent (pAmps)', ...
            'photocurrent', ...
            pCurrentVisualizedRange, ...
            irradianceAtTargetWavelengthInsteadOfRGBimage, ...
            targetWavelength, ...
            sprintf('%s-StimulusRaster-Increments', videoFilename));

      
       
       % Generate the video of the cone mosaic photocurrent response to the DECREMENTS stimulus raster
       fprintf('\nLoading stimulus (decrements) raster images. Please wait...')
       load(retinalImagesMatFileName, 'theListOfStimulusDecrementsRasterRetinalImages');
        for iTimePoint = 1:nTimePoints
            % Get the current retinal image
            theRetinalImage = theListOfStimulusDecrementsRasterRetinalImages{iTimePoint};
            theListOfStimulusRGBImages(iTimePoint,:,:,:) = single(oiGet(theRetinalImage, 'rgbimage'));
            theListOfStimulusRetinalIrradianceMaps(iTimePoint,:,:) = single(computeIrradianceInWattsPerMM2(theRetinalImage, targetWavelength));
        end
        clear 'theListOfStimulusDecrementsRasterRetinalImages';


        % Generate the video of the cone mosaic photocurrent response to the INCREMENTS stimulus raster
        generateMosaicActivationVideo(theConeMosaic, ...
            maxVideoTrials, ...
            nStimulusFramesPerTrial, ...
            theSceneTemporalSupportSeconds, ...
            spatialSupportX, spatialSupportY, ...
            theListOfStimulusRetinalIrradianceMaps, ...
            theListOfBackgroundRetinalIrradianceMaps, ...
            theListOfStimulusRGBImages, ...
            theListOfBackgroundRGBImages, ...
            photocurrentResponseDecrementsTimeSeries, ...
            photocurrentResponseDecrementsTimeSeriesNoisy, ...
            LconeIndicesVisualized, MconeIndicesVisualized, SconeIndicesVisualized, ...
            displayEyeMovements, ...
            photocurrentResponseTimeAxis, ...
            'photocurrent (pAmps)', ...
            'photocurrent', ...
            pCurrentVisualizedRange, ...
            irradianceAtTargetWavelengthInsteadOfRGBimage, ...
            targetWavelength, ...
            sprintf('%s-StimulusRaster-Decrements', videoFilename));

    else
        load(photocurrentsMatFileName, ...
            'theConeMosaic', ...
            'photocurrentResponseTimeSeries', ...
            'photocurrentResponseTimeSeriesNoisy', ...
            'photocurrentResponseTimeAxis');

        fprintf('\nLoading retinal image data from %s. Please wait ...', retinalImagesMatFileName);
        load(retinalImagesMatFileName, ...
            'theListOfBackgroundRasterRetinalImages', ...
            'theListOfStimulusRasterRetinalImages');
    

        % Determine the cone indices for which we will visualize their time series responses 
        [LconeIndicesVisualized, MconeIndicesVisualized, SconeIndicesVisualized] = ...
            determineVisualizedConeIndices(theConeMosaic);
    
        irradianceAtTargetWavelengthInsteadOfRGBimage = true;
        displayEyeMovements = true;
        yAxisLabel = 'photocurrent (pAmps)';


        % Next, visualize the response to the stimulus raster
        fprintf('\nGenerating the (periodic) stimulus raster OI sequence. Please wait ...');
        totalSimulationTimeSteps = size(theConeMosaic.fixEMobj.emPosArcMin,2);
        simulationTimeStepSeconds = photocurrentResponseTimeAxis(2)-photocurrentResponseTimeAxis(1);
        theSceneTemporalSupportSeconds = (0:(totalSimulationTimeSteps-1)) * simulationTimeStepSeconds;
          
        theStimulusOIsequence = oiArbitrarySequence(...
                theListOfStimulusRasterRetinalImages, ...
                theSceneTemporalSupportSeconds, ...
                'isPeriodic', true);
        fprintf('Done ! \n');


        fprintf('\nGenerating the background raster OI sequence. Please wait ...');
        theBackgroundOIsequence = oiArbitrarySequence(...
                theListOfBackgroundRasterRetinalImages, ...
                theSceneTemporalSupportSeconds);
        fprintf('Done ! \n');

        targetWavelength = 840;

        % Generate the video of the cone mosaic NOISY photocurrent response to the stimulus raster
        generateMosaicActivationVideo(theConeMosaic, theStimulusOIsequence, theBackgroundOIsequence, photocurrentResponseTimeSeriesNoisy, ...
            LconeIndicesVisualized, MconeIndicesVisualized, SconeIndicesVisualized, ...
            displayEyeMovements, ...
            photocurrentResponseTimeAxis, yAxisLabel, ...
            'photocurrent', ...
            pCurrentVisualizedRange, ...
            irradianceAtTargetWavelengthInsteadOfRGBimage, ...
            targetWavelength, ...
            sprintf('%sNoisy-StimulusRaster', videoFilename));
    
        % Generate the video of the cone mosaic NOISE-FREE photocurrent response to the stimulus raster
        generateMosaicActivationVideo(theConeMosaic, theStimulusOIsequence, theBackgroundOIsequence, photocurrentResponseTimeSeries, ...
            LconeIndicesVisualized, MconeIndicesVisualized, SconeIndicesVisualized, ...
            displayEyeMovements, ...
            photocurrentResponseTimeAxis, yAxisLabel, ...
            'photocurrent', ...
            pCurrentVisualizedRange, ...
            irradianceAtTargetWavelengthInsteadOfRGBimage, ...
            targetWavelength, ...
            sprintf('%s-StimulusRaster', videoFilename));
    end

end

function generateMosaicActivationVideo(theConeMosaic, ...
    maxVideoTrials, ...
    nStimulusFramesPerTrial, ...
    retinalImageSequenceTimeAxis, ...
    spatialSupportX, spatialSupportY, ...
    theListOfStimulusRetinalIrradianceMaps, ...
    theListOfBackgroundRetinalIrradianceMaps, ...
    theListOfStimulusRGBImages, ...
    theListOfBackgroundRGBImages, ...
    mosaicNoiseFreeResponseTimeSeries, ...
    mosaicNoisyResponseTimeSeries, ...
    LconeIndicesVisualized, MconeIndicesVisualized, SconeIndicesVisualized, ...
    displayEyeMovements, responseTimeAxis, yAxisLabel, signalType, signalRange, ...
    irradianceAtTargetWavelengthInsteadOfRGBimage, targetWavelength, videoFilename)

    labeledConeIndices = [...
        LconeIndicesVisualized(:); ...
        MconeIndicesVisualized(:); ...
        SconeIndicesVisualized(:)];

    labeledConeIndices = [];

    LconeIndicesResponses = mosaicNoiseFreeResponseTimeSeries(:,:,LconeIndicesVisualized);
    MconeIndicesResponses = mosaicNoiseFreeResponseTimeSeries(:,:,MconeIndicesVisualized);
    SconeIndicesResponses = mosaicNoiseFreeResponseTimeSeries(:,:,SconeIndicesVisualized);
    
    LconeIndicesNoisyResponses = mosaicNoisyResponseTimeSeries(:,:,LconeIndicesVisualized);
    MconeIndicesNoisyResponses = mosaicNoisyResponseTimeSeries(:,:,MconeIndicesVisualized);
    SconeIndicesNoisyResponses = mosaicNoisyResponseTimeSeries(:,:,SconeIndicesVisualized);


    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1600 960], 'Color', [0.1 0.1 0.1]);
    colsNum = 3; 
    rowsNum = 2;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
               'rowsNum', rowsNum, ...
               'colsNum', colsNum, ...
               'heightMargin',  0.04, ...
               'widthMargin',    0.05, ...
               'leftMargin',     0.04, ...
               'rightMargin',    0.00, ...
               'bottomMargin',   0.05, ...
               'topMargin',      0.02);


    ax1 = subplot('Position', subplotPosVectors(1,1).v);
    ax2 = subplot('Position', subplotPosVectors(2,1).v);
    ax3 = axes('Position', [0.38 0.06 0.6 0.92]);

    videoOBJ = VideoWriter(videoFilename, 'MPEG-4');  % H264format 
    videoOBJ.FrameRate = 60;
    videoOBJ.Quality = 100;
    videoOBJ.open();

   
    nTrials = size(mosaicNoiseFreeResponseTimeSeries,1);
    nTrials = min([nTrials maxVideoTrials]);

    visualizedConeAperture = 'lightCollectingArea5Sigma';
    backgroundColor = [0.5 0.5 0.5];

    if (isempty(signalRange))
        activationRange = [min(mosaicNoiseFreeResponseTimeSeries(:)) max(mosaicNoiseFreeResponseTimeSeries(:))]; 
    elseif (numel(signalRange) == 1)
        activationRange = signalRange * [-1 1];
    else
        activationRange = signalRange;
    end

    
    visualizedFOV = max(spatialSupportX) - min(spatialSupportX);
    domainVisualizationLimits(1:2) = theConeMosaic.eccentricityDegs(1) + 0.51*visualizedFOV*[-1 1];
    domainVisualizationLimits(3:4) = theConeMosaic.eccentricityDegs(2) + 0.51*visualizedFOV*[-1 1];
    domainVisualizationTicks.x = -2:0.2:2;
    domainVisualizationTicks.y = -2:0.2:2;

    
    % Subtract the background response (estimated from non-causal times)
    preStimulusTimeIndices = find(responseTimeAxis<0);
    if (~isempty(preStimulusTimeIndices)) && (~strcmp(signalType, 'excitations'))
        % Mean over all trials and non-causal time bins
        backgroundResponses = mean(mean(mosaicNoiseFreeResponseTimeSeries(:,1:preStimulusTimeIndices(end),:),2),1);
    
        % Subtract the bakgroundResponses when visualizing the mosaic activation
        mosaicNoiseFreeResponseTimeSeries = bsxfun(@minus, mosaicNoiseFreeResponseTimeSeries, backgroundResponses);
    end

    mosaicActivationRange = 0.5*max(abs(mosaicNoiseFreeResponseTimeSeries(:)))*[-1 1];

    irradianceRange = [0 max(theListOfStimulusRetinalIrradianceMaps(:))];
    dT = responseTimeAxis(2)-responseTimeAxis(1);
    singleFrameDuration = numel(retinalImageSequenceTimeAxis) * (retinalImageSequenceTimeAxis(2)-retinalImageSequenceTimeAxis(1));
    stimulusDuration = singleFrameDuration * nStimulusFramesPerTrial;
    for iTrial = 1:nTrials
        for iTimePoint = 1:size(LconeIndicesResponses,2) 

            iTimePointMod = iTimePoint;
            if (~strcmp(signalType, 'excitations'))
                % Photocurrents, which have pre- and post-stimulus time bins
                currentTime = responseTimeAxis(iTimePoint);
                
                if (currentTime < 0) || (currentTime > stimulusDuration)
                    % show the background raster
                    inStimulusInterval = false;
                    if (currentTime < 0)
                        sampleNum = numel(retinalImageSequenceTimeAxis)-round(-currentTime/dT);
                    else
                        sampleNum = round(currentTime/dT);
                    end
                    iTimePointMod = mod(sampleNum-1,numel(retinalImageSequenceTimeAxis)) + 1;
                else
                    % show the stimulus raster
                    inStimulusInterval = true;
                    sampleNum = round(currentTime/dT);
                    iTimePointMod = mod(sampleNum-1,numel(retinalImageSequenceTimeAxis)) + 1;
                end
            end % if (~strcmp(signalType, 'excitations'))

            if (irradianceAtTargetWavelengthInsteadOfRGBimage)
                if (inStimulusInterval)
                    % Stimulus radiance
                    imagesc(ax1, spatialSupportX, spatialSupportY, squeeze(theListOfStimulusRetinalIrradianceMaps(iTimePointMod,:,:)));
                else
                    % Background radiance
                    imagesc(ax1, spatialSupportX, spatialSupportY, squeeze(theListOfBackgroundRetinalIrradianceMaps(iTimePointMod,:,:)));
                end
                set(ax1, 'CLim', irradianceRange);
                colormap(ax1, 'gray');
                colorbar(ax1,'north', 'Color', [0.8 0.8 0.8], 'FontSize', 12, 'FontName', 'Spot mono');
                set(ax1, 'FontSize', 16, 'XColor', [0.5 0.5 0.5], 'YColor', [0.5 0.5 0.5]);
                title(ax1, ...
                    sprintf('simulated Tuten AOSLO display (irradiance, mWatts/mm^2 @ %dnm )', targetWavelength), ...
                    'FontSize', 14, 'Color', [0.5 0.5 0.5]);
            else
                if (inStimulusInterval)
                    % Stimulus OI
                    theRetinalImage = squeeze(theListOfStimulusRGBImages(iTimePointMod,:,:,:));
                else
                    % Background OI
                    theRetinalImage = squeeze(theListOfBackgroundRGBImages(iTimePointMod,:,:,:));
                end
                image(ax1, spatialSupportX, spatialSupportY, oiGet(theRetinalImage, 'rgbimage'));
                set(ax1, 'FontSize', 16, 'XColor', [0.5 0.5 0.5], 'YColor', [0.5 0.5 0.5]);
                title(ax1, sprintf('simulated Tuten AOSLO display'),'FontSize', 14, 'Color', [0.5 0.5 0.5]);
            end
            axis(ax1,'image');
            set(ax1, 'XTick', domainVisualizationTicks.x, 'YTick', domainVisualizationTicks.y);
            set(ax1, 'XLim', domainVisualizationLimits(1:2), 'YLim', domainVisualizationLimits(3:4), 'XTickLabel', {});
            set(ax1, 'Color', 'none');
            ylabel(ax1, 'eccentricity, y (degs)');


            emTimePointsVisualized = find(retinalImageSequenceTimeAxis <= responseTimeAxis(iTimePoint));
            if (displayEyeMovements) && (~isempty(emTimePointsVisualized))
                theConeMosaic.visualize('figureHandle', hFig, 'axesHandle', ax2, ...
                    'visualizedConeAperture', visualizedConeAperture, ...
                    'activation', mosaicNoiseFreeResponseTimeSeries(iTrial,iTimePoint,:), ...
                    'activationRange', mosaicActivationRange, ...
                    'currentEMposition', squeeze(theConeMosaic.fixEMobj.emPosArcMin(iTrial,emTimePointsVisualized(end),:))/60, ...
                    'displayedEyeMovementData', struct('trial', iTrial, 'timePoints', emTimePointsVisualized), ...
                    'domainVisualizationLimits', domainVisualizationLimits, ...
                    'domainVisualizationTicks', domainVisualizationTicks, ...
                    'horizontalActivationColorBarInside', true, ...
                    'colorbarFontSize', 10, ...
                    'labelConesWithIndices', labeledConeIndices, ...
                    'noYLabel', true, ...
                    'backgroundColor', backgroundColor, ...
                    'plotTitleFontSize', 14, ...
                    'plotTitleColor', [0.8 0.8 0.8], ...
                    'figureBackgroundColor', [], ...
                    'plotTitle',  sprintf('%s (trial no.%d, t = %2.2f ms)', yAxisLabel, iTrial, responseTimeAxis(iTimePoint)*1000));
            else
                theConeMosaic.visualize('figureHandle', hFig, 'axesHandle', ax2, ...
                    'visualizedConeAperture', visualizedConeAperture, ...
                    'activation', mosaicNoiseFreeResponseTimeSeries(iTrial,iTimePoint,:), ...
                    'activationRange', mosaicActivationRange, ...
                    'domainVisualizationLimits', domainVisualizationLimits, ...
                    'domainVisualizationTicks', domainVisualizationTicks, ...
                    'horizontalActivationColorBarInside', true, ...
                    'colorbarFontSize', 10, ...
                    'labelConesWithIndices', labeledConeIndices, ...
                    'backgroundColor', backgroundColor, ...
                    'plotTitleFontSize', 14, ...
                    'plotTitleColor', [0.8 0.8 0.8], ...
                    'figureBackgroundColor', [], ...
                    'plotTitle',  sprintf('%s (trial no.%d, t = %2.2f ms)', yAxisLabel, iTrial, responseTimeAxis(iTimePoint)*1000));
      
            end
            set(ax2, 'XColor', [0.5 0.5 0.5], 'YColor', [0.5 0.5 0.5]);

            idx = 1:iTimePoint;
            
            if (strcmp(signalType, 'excitations'))
                % Excitations
                scatter(ax3, responseTimeAxis(idx)*1000, squeeze(LconeIndicesResponses(iTrial, idx,:)), 50, ...
                    'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0.5 0.5], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.4, 'LineWidth', 0.5);
                hold(ax3, 'on');
                scatter(ax3, responseTimeAxis(idx)*1000, squeeze(MconeIndicesResponses(iTrial, idx,:)), 30, ...
                    'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', [0.5 0.8 0.5], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.4, 'LineWidth', 0.5);
                scatter(ax3, responseTimeAxis(idx)*1000, squeeze(SconeIndicesResponses(iTrial, idx,:)), 10, ...
                    'MarkerFaceColor', [0 0 1],  'MarkerEdgeColor', [0.5 0.5 1], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.4, 'LineWidth', 0.5);
            else
                % Photocurrents
                plot(ax3, responseTimeAxis(idx)*1000, squeeze(LconeIndicesNoisyResponses(iTrial, idx,:)), '-', ...
                    'Color', 0.5*[1 0 0], 'LineWidth', 0.5);
                hold(ax3, 'on');
                plot(ax3, responseTimeAxis(idx)*1000, squeeze(MconeIndicesNoisyResponses(iTrial, idx,:)), '-', ...
                    'Color', 0.5*[0 1 0], 'LineWidth', 0.5);
                plot(ax3, responseTimeAxis(idx)*1000, squeeze(SconeIndicesNoisyResponses(iTrial, idx,:)), '-', ...
                    'Color', 0.5*[0 0 1], 'LineWidth', 0.5);

                plot(ax3, responseTimeAxis(idx)*1000, squeeze(LconeIndicesResponses(iTrial, idx,:)), '-', ...
                    'Color', [1 0 0], 'LineWidth', 2.0);
                plot(ax3, responseTimeAxis(idx)*1000, squeeze(MconeIndicesResponses(iTrial, idx,:)), '-', ...
                    'Color', [0 1 0], 'LineWidth', 2.0);
                plot(ax3, responseTimeAxis(idx)*1000, squeeze(SconeIndicesResponses(iTrial, idx,:)), '-', ...
                    'Color', [0 0 1], 'LineWidth', 2.0);
                set(ax3, 'YTick', -100:10:100)
            end

            plot(ax3, [0 0], activationRange, 'w--', 'LineWidth', 1.5);
            hold(ax3, 'off');
            set(ax3, 'XLim', [responseTimeAxis(1) responseTimeAxis(end)]*1000, 'XTick', -300:10:500, 'YLim', activationRange, 'YTick', -100:5:0);
            set(ax3, 'XColor', [0.75 0.75 0.75], 'YColor', [0.75 0.75 0.75], 'Color', 'none', 'FontSize', 16);
            grid(ax3, 'on');
            box(ax3, 'off')
            xtickangle(ax3, 90);
            xlabel(ax3, 'time (ms)');
            ylabel(ax3, yAxisLabel);
            drawnow;
            videoOBJ.writeVideo(getframe(hFig));
        end % for iTimePoint

    end % for iTrial

    videoOBJ.close();
end




function visualizeRetinalImageAndConeExcitations(...
    retinalImagesMatFileName, coneExcitationsMatFileName, ...
    testIncrementDecrementScenes, videoFilename)

    fprintf('\nLoading cone excitation data from %s. Please wait ...', coneExcitationsMatFileName);
    load(coneExcitationsMatFileName, ...
                'theConeMosaic', ...
                'backgroundConeExcitationResponseTimeSeries', ...
                'backgroundConeExcitationsTimeAxis');

    fprintf('\nLoading retinal image data from %s. Please wait ...', retinalImagesMatFileName);
        load(retinalImagesMatFileName, ...
                'simulationTimeStepSeconds', ...
                'theListOfBackgroundRasterRetinalImages');


    % Determine the cone indices for which we will visualize their time series responses 
    [LconeIndicesVisualized, MconeIndicesVisualized, SconeIndicesVisualized] = ...
        determineVisualizedConeIndices(theConeMosaic);

    % Transform excitation counts to excitations rate
    dT = backgroundConeExcitationsTimeAxis(2)-backgroundConeExcitationsTimeAxis(1);
    backgroundConeExcitationResponseTimeSeries = backgroundConeExcitationResponseTimeSeries / dT;


    if (testIncrementDecrementScenes)
        targetWavelength = 680;
    else
        targetWavelength = 860;
    end

    % First, visualize the response to the background raster
    fprintf('\nGenerating the background raster OI sequence. Please wait ...');
    theSceneTemporalSupportSeconds = (0:(numel(theListOfBackgroundRasterRetinalImages)-1)) * simulationTimeStepSeconds;
        
    theBackgroundOIsequence = oiArbitrarySequence(...
            theListOfBackgroundRasterRetinalImages, ...
            theSceneTemporalSupportSeconds);
    fprintf('Done ! \n');

    irradianceAtTargetWavelengthInsteadOfRGBimage = true;
    displayEyeMovements = false;
    yAxisLabel = 'cone excitation rate (R*/sec)';

    % Generate the video of the cone mosaic excitations response to the background raster
    generateMosaicActivationVideo(theConeMosaic, theBackgroundOIsequence, [], backgroundConeExcitationResponseTimeSeries, ...
        LconeIndicesVisualized, MconeIndicesVisualized, SconeIndicesVisualized, ...
        displayEyeMovements, ...
        backgroundConeExcitationsTimeAxis, yAxisLabel, ...
        'excitations', ...
        [], ...
        irradianceAtTargetWavelengthInsteadOfRGBimage, ...
        targetWavelength, ...
        sprintf('%s-BackgroundRaster', videoFilename));
    
    clear 'theListOfBackgroundRasterRetinalImages'


    if (testIncrementDecrementScenes)

        fprintf('\nLoading retinal image data from %s. Please wait ...', retinalImagesMatFileName);
        load(retinalImagesMatFileName, 'theListOfStimulusDecrementsRasterRetinalImages');

        % Next, visualize the response to the stimulus raster - decrements
        fprintf('\nGenerating the (periodic) stimulus raster OI sequence (DECREMENTS). Please wait ...');
        totalSimulationTimeSteps = size(theConeMosaic.fixEMobj.emPosArcMin,2);
        theSceneTemporalSupportSeconds = (0:(totalSimulationTimeSteps-1)) * simulationTimeStepSeconds;
          
        theStimulusOIsequence = oiArbitrarySequence(...
                theListOfStimulusDecrementsRasterRetinalImages, ...
                theSceneTemporalSupportSeconds, ...
                'isPeriodic', true);
        fprintf('Done ! \n');
    
        clear 'theListOfStimulusDecrementsRasterRetinalImages';

        fprintf('\nLoading cone excitation (Decrements) data from %s. Please wait ...', coneExcitationsMatFileName);
        load(coneExcitationsMatFileName, ...
            'noiseFreeConeExcitationDecrementsResponseTimeSeries', ...
            'timeAxis');

        % Transform excitation counts to excitations rate
        dT = timeAxis(2)-timeAxis(1);
        noiseFreeConeExcitationDecrementsResponseTimeSeries = noiseFreeConeExcitationDecrementsResponseTimeSeries / dT;
    
    
        displayEyeMovements = true;
        yAxisLabel = 'cone excitation rate (R*/sec)';
        % Generate the video of the cone mosaic excitations response to the stimulus raster
        generateMosaicActivationVideo(theConeMosaic, theStimulusOIsequence, [], noiseFreeConeExcitationDecrementsResponseTimeSeries, ...
            LconeIndicesVisualized, MconeIndicesVisualized, SconeIndicesVisualized, ...
            displayEyeMovements, ...
            timeAxis, yAxisLabel, ...
            'excitations', ...
            [], ...
            irradianceAtTargetWavelengthInsteadOfRGBimage, ...
            targetWavelength, ...
            sprintf('%s-StimulusRaster-Decrements', videoFilename));

        clear 'noiseFreeConeExcitationDecrementsResponseTimeSeries';


        % Repeat for increments
        fprintf('\nLoading retinal image data (increments) from %s. Please wait ...', retinalImagesMatFileName);
        load(retinalImagesMatFileName, 'theListOfStimulusIncrementsRasterRetinalImages');

        % Next, visualize the response to the stimulus raster - increments
        fprintf('\nGenerating the (periodic) stimulus raster OI sequence (INCREMENTS). Please wait ...');
        theStimulusOIsequence = oiArbitrarySequence(...
                theListOfStimulusIncrementsRasterRetinalImages, ...
                theSceneTemporalSupportSeconds, ...
                'isPeriodic', true);
        fprintf('Done ! \n');
    
        clear 'theListOfStimulusIncrementsRasterRetinalImages';
        
        fprintf('\nLoading cone excitation (Increments) data from %s. Please wait ...', coneExcitationsMatFileName);
        load(coneExcitationsMatFileName, ...
            'noiseFreeConeExcitationIncrementsResponseTimeSeries', ...
            'timeAxis');

        % Transform excitation counts to excitations rate
        dT = timeAxis(2)-timeAxis(1);
        noiseFreeConeExcitationIncrementsResponseTimeSeries = noiseFreeConeExcitationIncrementsResponseTimeSeries / dT;
    
    
        % Generate the video of the cone mosaic excitations response to the stimulus raster
        generateMosaicActivationVideo(theConeMosaic, theStimulusOIsequence, [], noiseFreeConeExcitationIncrementsResponseTimeSeries, ...
            LconeIndicesVisualized, MconeIndicesVisualized, SconeIndicesVisualized, ...
            displayEyeMovements, ...
            timeAxis, yAxisLabel, ...
            'excitations', ...
            [], ...
            irradianceAtTargetWavelengthInsteadOfRGBimage, ...
            targetWavelength, ...
            sprintf('%s-StimulusRaster-Increments', videoFilename));

    else

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
    
    
        displayEyeMovements = true;
        yAxisLabel = 'cone excitation rate (R*/sec)';
        % Generate the video of the cone mosaic excitations response to the stimulus raster
        generateMosaicActivationVideo(theConeMosaic, theStimulusOIsequence, [], noiseFreeConeExcitationResponseTimeSeries, ...
            LconeIndicesVisualized, MconeIndicesVisualized, SconeIndicesVisualized, ...
            displayEyeMovements, ...
            timeAxis, yAxisLabel, ...
            'excitations', ...
            [], ...
            irradianceAtTargetWavelengthInsteadOfRGBimage, ...
            targetWavelength, ...
            sprintf('%s-StimulusRaster', videoFilename));
    end


end



function [LconeIndicesVisualized, MconeIndicesVisualized, SconeIndicesVisualized] = ...
        determineVisualizedConeIndices(theConeMosaic)
        

    theROI = regionOfInterest(...
        'geometryStruct', struct(...
            'units', 'degs', ...
            'shape', 'rect', ...
            'center', [-0.9 -0.1], ...
            'width', 0.25, ...
            'height', 0.25, ...
            'rotation', 0.0...
        ));
    

    idx = theROI.indicesOfPointsInside(theConeMosaic.coneRFpositionsDegs(theConeMosaic.lConeIndices,:));
    LconeIndicesVisualized = theConeMosaic.lConeIndices(idx);

    idx = theROI.indicesOfPointsInside(theConeMosaic.coneRFpositionsDegs(theConeMosaic.mConeIndices,:));
    MconeIndicesVisualized = theConeMosaic.mConeIndices(idx);
     
    idx = theROI.indicesOfPointsInside(theConeMosaic.coneRFpositionsDegs(theConeMosaic.sConeIndices,:));
    SconeIndicesVisualized = theConeMosaic.sConeIndices(idx);

end




function [theNullScene, theIncrementsEscene0degs, theDecrementsEscene0degs, thePresentationDisplay, sceneParams] = ...
    generateIncrementDecrementScenes(testESizeDeg, AOPrimaryWls, AOAOCornealPowersUW, ...
    visualizeTheSceneRadiance, figureFileBase)

    temporalModulationParams_xShiftPerFrame = [0 5/60 0];
    temporalModulationParams_yShiftPerFrame = [0 5/60 0];
    
    % We keep the imaging channel at max gain
    imagingChannelGain = 1.0;

    % The background: Stimulus channel at 0.5
    stimulusChannelGain = 0.5;
    backgroundRGB = [imagingChannelGain stimulusChannelGain 0];

    % The E increments engines: Stimulus channel at 1.0
    stimulusChannelGain = 1.0;
    foregroundRGB = [imagingChannelGain stimulusChannelGain 0];
    [sce0Increments, sce90Increments, sce180Increments, sce270Increments, sceBg, sceneParams] = t_BerkeleyAOTumblingESceneEngine( ...
        'temporalModulationParams_xShiftPerFrame', temporalModulationParams_xShiftPerFrame, ...
        'temporalModulationParams_yShiftPerFrame', temporalModulationParams_yShiftPerFrame, ...
        'AOPrimaryWls', AOPrimaryWls, ...
        'AOCornealPowersUW', AOAOCornealPowersUW, ...
        'temporalModulationParams_backgroundRGBPerFrame', [backgroundRGB; backgroundRGB; backgroundRGB], ...
        'temporalModulationParams_foregroundRGBPerFrame', [backgroundRGB; foregroundRGB; backgroundRGB]);

    frameNo = 2;

    theNullSequence = sceBg.compute(testESizeDeg);
    theNullScene = theNullSequence{frameNo};
    if (visualizeTheSceneRadiance)
        hFig = visualizeSceneRadiance(200, theNullScene, 'background', sceBg.presentationDisplay);
        NicePlot.exportFigToPDF(fullfile(figureFileBase,'BackgroundScene.pdf'), hFig, 300);
    end

    theEsceneSequence = sce0Increments.compute(testESizeDeg);
    theIncrementsEscene0degs = theEsceneSequence{frameNo};
    if (visualizeTheSceneRadiance)
        hFig = visualizeSceneRadiance(300, theIncrementsEscene0degs, 'increments', sceBg.presentationDisplay);
        NicePlot.exportFigToPDF(fullfile(figureFileBase,'IncrementsScene.pdf'), hFig, 300);
    end
    

    % The E decrements engines: Stimulus channel at 0.0
    stimulusChannelGain = 0.0;
    foregroundRGB = [imagingChannelGain stimulusChannelGain 0];
    [sce0Decrements, sce90Decrements, sce180Decrements, sce270Decrements] = t_BerkeleyAOTumblingESceneEngine( ...
        'temporalModulationParams_xShiftPerFrame', temporalModulationParams_xShiftPerFrame, ...
        'temporalModulationParams_yShiftPerFrame', temporalModulationParams_yShiftPerFrame, ...
        'AOPrimaryWls', AOPrimaryWls, ...
        'AOCornealPowersUW', AOAOCornealPowersUW, ...
        'temporalModulationParams_backgroundRGBPerFrame', [backgroundRGB; backgroundRGB; backgroundRGB], ...
        'temporalModulationParams_foregroundRGBPerFrame', [backgroundRGB; foregroundRGB; backgroundRGB]);

    
    theEsceneSequence = sce0Decrements.compute(testESizeDeg);
    theDecrementsEscene0degs = theEsceneSequence{frameNo};
    if (visualizeTheSceneRadiance)
        hFig = visualizeSceneRadiance(100, theDecrementsEscene0degs, 'decrements', sceBg.presentationDisplay);
        NicePlot.exportFigToPDF(fullfile(figureFileBase,'DecrementsScene.pdf'), hFig, 300);
    end

    thePresentationDisplay = sceBg.presentationDisplay;
end

function RGBsettingsImage = visualizeRGBimageFromXYZimage(XYZimage,presentationDisplay)

    displayLinearRGBToXYZ = displayGet(presentationDisplay, 'rgb2xyz');
    displayXYZToLinearRGB = inv(displayLinearRGBToXYZ);
    RGBsettingsImage = imageLinearTransform(XYZimage, displayXYZToLinearRGB);

    totalPixels = numel(RGBsettingsImage)
    pixelsNumWithRGBsettingsLessThan0 = numel(find(RGBsettingsImage<0))
    [min(RGBsettingsImage(:)) max(RGBsettingsImage(:))]

    pixelsNumWithRGBsettingsGreaterThan1 = numel(find(RGBsettingsImage>1))
    RGBsettingsImage(RGBsettingsImage<0) = 0;
    RGBsettingsImage(RGBsettingsImage>1) = 1;
end


function hFig = visualizeSceneRadiance(figNo, scene, sceneLabel, presentationDisplay)
    
    if (1==2)
        XYZimage = sceneGet(scene, 'xyz');
        sceneRGBsettings = visualizeRGBimageFromXYZimage(XYZimage,presentationDisplay);
    end

    renderFlag = -2;
    sceneRGBsettings = sceneShowImage(scene,renderFlag);
    rChannel = sceneRGBsettings(:,:,1);
    sceneRGBsettings(:,:,2) = 0;
    sceneRGBsettings(:,:,3) = 0;
    %sceneRGBsettings = sceneGet(scene, 'rgbimage');
    prctile(rChannel(:), [0 1 50 99 100])
    pause
    viewingDistance = sceneGet(scene, 'distance');
    spatialSupportMilliMeters = sceneGet(scene, 'spatial support', 'mm');
    spatialSupportDegs = 2 * atand(spatialSupportMilliMeters/1e3/2/viewingDistance);
    
    spatialSupport = spatialSupportDegs;
    spatialSupportX = squeeze(spatialSupport(1,:,1));
    spatialSupportY = squeeze(spatialSupport(:,1,2));

    [~,targetRow] = min(abs(spatialSupportY-0.01));
    [~,targetCol] = min(abs(spatialSupportX-0.09));

    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 1100 1080], 'Color', [1 1 1]);
    ax = subplot(2,2,1);
    image(ax,spatialSupportX, spatialSupportY, lrgb2srgb(sceneRGBsettings));
    hold(ax, 'on')
    plot(ax, zeros(size(spatialSupportY)) + spatialSupportX(targetCol), spatialSupportY, 'g:', 'LineWidth', 1.0);
    axis(ax,'xy'); axis(ax, 'square')
    ylabel('space, y (degs)')
    title(ax, sprintf('%s scene',sceneLabel));
    set(ax, 'FontSize', 16);

    ax = subplot(2,2,3);
    luminanceMap = sceneGet(scene, 'luminance');
    imagesc(ax,spatialSupportX, spatialSupportY, luminanceMap);
    maxDisplayedLuminance = 15000;
    set(ax, 'CLim', [0 maxDisplayedLuminance]);
    hold(ax, 'on');
    plot(ax, zeros(size(spatialSupportY)) + spatialSupportY(targetCol), spatialSupportY, 'r--');
    colormap(ax,gray);
    axis(ax,'xy'); axis(ax, 'square')
    xlabel('space, x (degs)')
    ylabel('space, y (degs)')
    title(ax, sprintf('luminance (cd/m2), range: [%2.1f - %2.0f]', min(luminanceMap(:)), max(luminanceMap(:))));
    colorbar(ax, 'north', 'Color', [0.8 0.8 0.8])
    set(ax, 'FontSize', 16);

    photons = sceneGet(scene, 'energy');
    wave = sceneGet(scene, 'wave');

    ax = subplot(2,2,2);
    spaceWave = squeeze(photons(:, targetCol, :));
    meanE = mean(spaceWave, 1);
    max(meanE)
    imagesc(ax,  wave, spatialSupportY, log10(spaceWave));
    axis(ax,'xy'); axis(ax, 'square');
    set(ax, 'CLim', log10([0.01 300]), 'XTick', 500:20:900, 'XLim', [500 880]);
    grid(ax, 'on')
    ylabel(ax,'space, y (degs)')
     
    set(ax, 'FontSize', 16);
    colormap(ax,gray);
    colorbar(ax, 'north', 'Color', [1 0.2 0.2]);
    title('log10 (peak power) [Watts]')
    
    ax = subplot(2,2,4);
    plot(ax, wave, mean(spaceWave, 1), 'r-', 'LineWidth', 1.5);
    axis(ax, 'square');
    set(ax, 'XLim', [500 880], 'YScale', 'log', 'YLim', [0.01 350], ...
        'YTick', [0.01 0.03 0.1 0.3 1 3 10 30 100 300], 'XTick', 500:20:900);
    grid(ax, 'on')
    xlabel(ax, 'wavelength (nm)');
    ylabel(ax, 'peak power (Watts)');
    set(ax, 'FontSize', 16);
end

function [theNullScene, theEscene0degs, thePresentationDisplay, sceneParams] = generateScenes(testESizeDeg, visualizeTheSceneRadiance, figureFileBase)

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

    if (visualizeTheSceneRadiance)
        visualizeSceneRadiance(100, theEscene0degs , 'decrements', sceBg.presentationDisplay);
    end

    thePresentationDisplay = sceBg.presentationDisplay;
end

function theConeMosaic = generateConeMosaic(...
    mosaicEccDegs, mosaicSizeDegs, mosaicIntegrationTimeSeconds, wave)
    
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
                'integrationTime', mosaicIntegrationTimeSeconds, ...
                'wave', wave ...
                );
end

function generateFixationalEyeMovements(theConeMosaic, ...
    simulationDurationSeconds, nTrials, observerCharacteristics)

    % Fixational eye movements
    eyeMovementsPerTrial = round(simulationDurationSeconds/theConeMosaic.integrationTime);
    
    fixEMobj = fixationalEM();
    switch (observerCharacteristics)
        case 'fast'
            fixEMobj.controlGamma = 0.17;
            fixEMobj.feedbackGain = 0.15;
        case 'slow'
            fixEMobj.controlGamma = 0.50;
            fixEMobj.feedbackGain = 0.12;
        otherwise
            % default params of gamma, feeedback
    end % switch observerCharacteristics

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

function theOI = generateOptics(simulateAdaptiveOpticsviewingConditions, theConeMosaic, inFocusWavelength)

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
                   'inFocusWavelength', inFocusWavelength, ...
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

    s = sceneGet(theTestScene, 'size');
    rasterLinesNum = s(1);
    pixelsPerRasterLine = s(2);

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
        cropRetinalImagesForConeMosaicSize, ...
        visualizeEachRetinalImage, thePresentationDisplay, testIncrementDecrementScenes)

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
    if (testIncrementDecrementScenes)
        % Compute the power at 680 nm
        irradiancePowerWattsTargetWavelengthPerMMsquared = computeIrradianceInWattsPerMM2(theRetinalImage, 680);
        nonRasterizedEscenePower = mean(irradiancePowerWattsTargetWavelengthPerMMsquared(:));
    else
        % Compute the power at 840 nm
        irradiancePowerWattsTargetWavelengthPerMMsquared = computeIrradianceInWattsPerMM2(theRetinalImage, 840);
        nonRasterizedEscenePower = mean(irradiancePowerWattsTargetWavelengthPerMMsquared(:));
        fprintf(2,'\nNon-rasterized scene power at 840 nm is: %E W/mm2.', nonRasterizedEscenePower);
        fprintf(2,'\nWill says it should be: 8.37E-04 W/mm2.\n');
    end


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

    % Each pixel is ON for pixelTimeOnSeconds. But in the simulation it is
    % on for simulationTimeStepSeconds. Factor to account for that: 
    factorToAccountForPixelOnTime = pixelTimeOnSeconds/simulationTimeStepSeconds;

    radianceMultiplierToAccountForRasterSimulation = ...
        radianceMultiplierToAccountForRasterSimulation * ...
        factorToAccountForPixelOnTime;


    % Generate a scene for each simulation time step and its associated retinal image
    theListOfRetinalImages = cell(1, simulationTimeStepsNum);

    for iTimeStep = 1:simulationTimeStepsNum

        fprintf('Generating scene for simulation time step %d of %d\n', iTimeStep, simulationTimeStepsNum);   

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

        % Set the new radiance
        theFrameScene = sceneSet(theTestScene, 'photons', theFrameScenePhotons);

        % Compute the OI
        theRetinalImageAtThisSimulationTimeStep = oiCompute(theOI, theFrameScene, 'pad value', 'zero');

       
        % Crop the OI
        if (~isempty(cropRetinalImagesForConeMosaicSize))
            
            if (iTimeStep == 1)
                % Compute the cropping rect
                spatialSupportMM = oiGet(theRetinalImageAtThisSimulationTimeStep, 'spatial support', 'mm');
                theOptics = oiGet(theRetinalImageAtThisSimulationTimeStep, 'optics');
                focalLength = opticsGet(theOptics, 'focal length');
                mmPerDegree = focalLength*tand(1)*1e3;
                spatialSupportDegs = spatialSupportMM/mmPerDegree;
                spatialSupportXdegs = spatialSupportDegs(1,:,1);
                spatialSupportYdegs = spatialSupportDegs(:,1,2);
 
                idx = find(abs(spatialSupportXdegs) <= 0.5*cropRetinalImagesForConeMosaicSize(1));
                leftColumn = idx(1);
                colsNum = numel(idx);
    
                idx = find(abs(spatialSupportYdegs) <= 0.5*cropRetinalImagesForConeMosaicSize(2));
                topRow = idx(1);
                rowsNum = numel(idx);

                croppingRect = [topRow leftColumn rowsNum colsNum];
            end
            
            if (visualizeEachRetinalImage)
                theRetinalIlluminance = oiGet(theRetinalImageAtThisSimulationTimeStep, 'illuminance');
                theRetinalIlluminanceRange = [0 max(theRetinalIlluminance(:))];
            end

            % Crop the image
            bytesNumOfRetinalImageBeforeCropping = numel(getByteStreamFromArray(theRetinalImageAtThisSimulationTimeStep));

            theRetinalImageAtThisSimulationTimeStep = oiCrop(...
                theRetinalImageAtThisSimulationTimeStep, croppingRect);

            bytesNumOfRetinalImageAfterCropping = numel(getByteStreamFromArray(theRetinalImageAtThisSimulationTimeStep));
            bytesRatioAfterCropping = bytesNumOfRetinalImageAfterCropping/bytesNumOfRetinalImageBeforeCropping

            spatialSupportMM = oiGet(theRetinalImageAtThisSimulationTimeStep, 'spatial support', 'mm');
            mmPerDegree = focalLength*tand(1)*1e3;
            spatialSupportDegs = spatialSupportMM/mmPerDegree;
            spatialSupportXdegs = spatialSupportDegs(1,:,1);
            spatialSupportYdegs = spatialSupportDegs(:,1,2);
        end


        if (visualizeEachRetinalImage)
            figure(33); clf;
            %XYZimage = oiGet(theRetinalImageAtThisSimulationTimeStep, 'XYZ');
            %RGBsettingsImage = visualizeRGBimageFromXYZimage(XYZimage,thePresentationDisplay);
            %image(spatialSupportXdegs, spatialSupportYdegs, lrgb2srgb(RGBsettingsImage));
            theRetinalIlluminance = oiGet(theRetinalImageAtThisSimulationTimeStep, 'illuminance');
            imagesc(spatialSupportXdegs, spatialSupportYdegs, theRetinalIlluminance);
            set(gca, 'CLim', theRetinalIlluminanceRange);
            axis 'image'
            colormap(gray(1024));
        end

        theListOfRetinalImages{iTimeStep} = theRetinalImageAtThisSimulationTimeStep;
    end
end


function [photocurrentResponseTimeSeries, photocurrentResponseTimeAxis, photocurrentResponseTimeSeriesNoisy] = ...
        computeConeMosaicPhotoCurrentsResponse(theConeMosaic, ...
            coneExcitationResponseBackground, ...
            coneExcitationResponseTimeSeries, ...
            warmBackgroundRefeshCyclesNum, warmUpPeriodSecondsToIncludeInVisualization, ...
            coolDownBackgroundRefeshCyclesNum, subtractBackgroundPhotoCurrents, debugWarmUpTime)


    % Concatenate (pre-stimulus) warm up background
    backgroundConeExcitationsRepMatPreStimulus = repmat(coneExcitationResponseBackground, [size(coneExcitationResponseTimeSeries,1) warmBackgroundRefeshCyclesNum 1 ]);
    coneExcitationResponseTimeSeries = cat(2, backgroundConeExcitationsRepMatPreStimulus, coneExcitationResponseTimeSeries);

    % Concatenate (post-stimulus) cool down background
    backgroundConeExcitationsRepMatPostStimulus = repmat(coneExcitationResponseBackground, [size(coneExcitationResponseTimeSeries,1) coolDownBackgroundRefeshCyclesNum 1 ]);
    coneExcitationResponseTimeSeries = cat(2, coneExcitationResponseTimeSeries, backgroundConeExcitationsRepMatPostStimulus);

    % mean background excitation count over warm-up period
    backgroundConeExcitationCounts = squeeze(mean(coneExcitationResponseBackground,2));
    warmUpSamplesToIncludeInVisualization = ceil(warmUpPeriodSecondsToIncludeInVisualization / theConeMosaic.integrationTime);
    actualWarmUpPeriodSecondsToIncludeInVisualization = warmUpSamplesToIncludeInVisualization * theConeMosaic.integrationTime;

    warmupPeriodSeconds = size(backgroundConeExcitationsRepMatPreStimulus,2)*theConeMosaic.integrationTime;
    
    [photocurrentResponseTimeSeries, photocurrentResponseTimeAxis, photocurrentResponseTimeSeriesNoisy] = ...
         computeFullPhotocurrentModelResponse(theConeMosaic,  backgroundConeExcitationCounts, coneExcitationResponseTimeSeries, ...
         theConeMosaic.integrationTime, warmupPeriodSeconds, actualWarmUpPeriodSecondsToIncludeInVisualization, ...
         subtractBackgroundPhotoCurrents, debugWarmUpTime);
end




function [photocurrentResponseTimeSeries, photocurrentResponseTimeAxis, photocurrentResponseTimeSeriesNoisy] = ...
    computeFullPhotocurrentModelResponse(theConeMosaic, backgroundConeExcitationCounts, stimulusConeExcitationResponseTimeSeries, ...
    coneMosaicIntegrationTime, warmupPeriodSeconds, warmUpPeriodSecondsToIncludeInVisualization, ...
    subtractBackgroundPhotoCurrents, debugWarmUpTime)

    osTimeStep = 1e-5;
    integerMultiplier = round(coneMosaicIntegrationTime/osTimeStep);
    osTimeStep = coneMosaicIntegrationTime/integerMultiplier;

    assert(rem(coneMosaicIntegrationTime, osTimeStep) == 0, 'coneMosaic.intergrationTime must be an integral multiple of os.time step');

    upSampleFactor = floor(coneMosaicIntegrationTime / osTimeStep);
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

        if (theConeMosaic.coneTypes(iCone) == cMosaic.SCONE_ID)
            %fprintf('S-cone. Skipping photocurrent computation.\n')
            continue
        end
        
        % Setup biophysical model of outer segment for each cone 
        os = osBioPhys('eccentricity',0);
        os.timeStep = osTimeStep;
        os.set('noise flag', 'none');

        % Set the state
        backgroundConeExcitationRate = backgroundConeExcitationCounts(iCone)/coneMosaicIntegrationTime;
        theState = os.osAdaptSteadyState(backgroundConeExcitationRate, [1 1]);

        theState.timeStep = os.timeStep;
        os.setModelState(theState);

        for iTrial = 1:nTrials

            % Get the cone excitations count time series
            singleConeSingleTrialConeExcitationsCountResponse = ...
                squeeze(stimulusConeExcitationResponseTimeSeries(iTrial, :, iCone));

           
            % upsample counts time series to the os.timeStep timebase
            singleConeSingleTrialConeExcitationsCountResponse = ...
                interp1(tIn, singleConeSingleTrialConeExcitationsCountResponse, tOut, 'nearest');
     

            % Convert excitation counts to excitation rates
            singleConeSingleTrialConeExcitationsRate = ...
                singleConeSingleTrialConeExcitationsCountResponse / coneMosaicIntegrationTime;
              

            % Compute full photocurrent model
            theOSphotoCurrentResponse = reshape(...
                os.osAdaptTemporal(reshape(singleConeSingleTrialConeExcitationsRate, [1 1 numel(singleConeSingleTrialConeExcitationsRate)])), ...
                [1 theOSphotocurrentResponseTimeBinsNum 1]);

            fprintf(' pCurrent response for cone %d (type:%d) of %d: range: [%2.2f %2.2f]\n', iCone, theConeMosaic.coneTypes(iCone), nCones, min(theOSphotoCurrentResponse), max(theOSphotoCurrentResponse));

            % Downsample to the original cone mosaic time base
            thePcurrentResponse = interp1(theOSphotoCurrentResponseTimeAxis, theOSphotoCurrentResponse, photocurrentResponseTimeAxis, 'linear');
            photocurrentResponseTimeSeries(iTrial,:,iCone) = thePcurrentResponse;
            
        end % iTrial
        
    end % iCone

    % Add noise

    % osAddNoise() assumes that the noise-free photocurrent signal, here photocurrentResponseTimeSeries,
    % is arranged in a 3D matrix with the following format: [row, col, time]
    % So we have to permute the photocurrentResponseTimeSeries which is in [trial, time, mCones] format

    % Move the time dimension to the 3rd place
    photocurrentResponseTimeSeries = permute(photocurrentResponseTimeSeries, [1 3 2]);

    % Add photocurrent noise
    photocurrentResponseTimeSeriesNoisy = osAddNoise(...
        photocurrentResponseTimeSeries, ...
        'sampTime', coneMosaicIntegrationTime);

    % Move the time dimension back to the 2nd place
    photocurrentResponseTimeSeries = permute(photocurrentResponseTimeSeries, [1 3 2]);
    photocurrentResponseTimeSeriesNoisy = permute(photocurrentResponseTimeSeriesNoisy, [1 3 2]);

    % Only keep response at t >= -warmUpPeriodSecondsToIncludeInVisualization;
    if (warmupPeriodSeconds > 0)

        photocurrentResponseTimeAxis = photocurrentResponseTimeAxis - warmupPeriodSeconds;
        
        if (subtractBackgroundPhotoCurrents)
            % The delta-photocurrent (pCurrent - background)
            idx = find(photocurrentResponseTimeAxis<=0);
            baselinePhotoCurrents = photocurrentResponseTimeSeries(1,idx(end),:);

            photocurrentResponseTimeSeries = bsxfun(@minus, photocurrentResponseTimeSeries, baselinePhotoCurrents);
            photocurrentResponseTimeSeriesNoisy = bsxfun(@minus, photocurrentResponseTimeSeriesNoisy, baselinePhotoCurrents);
        end

        idx = find(photocurrentResponseTimeAxis >= -warmUpPeriodSecondsToIncludeInVisualization);
        photocurrentResponseTimeAxis = photocurrentResponseTimeAxis(idx);
        photocurrentResponseTimeSeries = photocurrentResponseTimeSeries(:,idx,:);
        photocurrentResponseTimeSeriesNoisy = photocurrentResponseTimeSeriesNoisy(:,idx,:);

        if (debugWarmUpTime)
            figure(3333); clf;
            plot(photocurrentResponseTimeAxis, squeeze(photocurrentResponseTimeSeries(1,:,:)), 'k-');
            xlabel('time');
            ylabel('cone no')
            colormap(gray)
            pause
        end
    end

end
