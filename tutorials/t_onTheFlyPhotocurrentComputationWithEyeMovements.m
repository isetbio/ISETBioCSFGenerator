function t_onTheFlyPhotocurrentComputationWithEyeMovements(options)
% Demonstrates how to compute photocurrent responses on the fly
%
% Syntax:
%    t_onTheFlyPhotocurrentComputationWithEyeMovements
%
% Description:
%    Demonstrates how to generate an expanding bulleye stimulus, present it
%    to a cone mosaic under fixational eye movements and compute 
%    the cone mosaic excitations and outer segment photocurrent response.
%    The code shows explicity how to compute noise-free and noisy photocurrent responses 
%    by calls to the low-level os routines.
%    The computed responses are visualized and animation videos of the
%    spatiotemporal activation of the cone mosaic are generated. 
%
%    This code was used for generate material for Randolph Blake's
%    chapter on motion perception, for which he asked us for help with
%    ISETBio-based illustrations
%
% Inputs:
%    None.
%
% Outputs:
%   None.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%   t_neuralResponseCompute, t_thresholdEngine, t_neuralResponseCompute,
%   t_responseClassifier.
%

% History:
%    02/19/2025  NPC  Wrote it.
%    03/06/2025  NPC  Fast parameters

% Examples:
%{

    % Run with defaults, computing new responses 
    t_onTheFlyPhotocurrentComputationWithEyeMovements();

    % Do not re-run the simulation, simply load results from a
    % previous simulation so we can generate figures and videos
    t_onTheFlyPhotocurrentComputationWithEyeMovements(...
        'loadResponsesFromPreviousSimulation', true);
%}
%{
    % ETTBSkip
    % Set to skip on autorun because it takes a long time
    %
    % Compute new responses using high-res parameters which take 
    % longer to compute.     
    t_onTheFlyPhotocurrentComputationWithEyeMovements(...
        'fastParameters', false);
%}

    arguments
        % Fast parameters?
        options.fastParameters (1,1) logical = true;
    
        % Whether to load responses already computed 
        % If not, we compute new responses
        options.loadResponsesFromPreviousSimulation (1,1) logical = false;
    end

    % Parse input
    fastParameters = options.fastParameters;

    % Whether to re-run the simulation or simply load results from a
    % previous simulation so we can generate figures and videos
    loadResponsesFromPreviousSimulation = options.loadResponsesFromPreviousSimulation;

    % Close figs
    close all;

    % Whether to include fixational eye movements or not
    fixationalEyeMovements = true;

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

    if (fixationalEyeMovements)
        simFileName = sprintf('simulationWithFixationalEM.mat');
    else
        simFileName = sprintf('simulationWithoutFixationalEM.mat');
    end
    simFileName = fullfile(projectBaseDir,'local',mfilename,'results',simFileName);
    
    % Whether to reformat the exported AVI videos to MP4 format
    reformatExportedAVIvideoToMP4format = ~true;

    if (loadResponsesFromPreviousSimulation)
        % Load previous results
        load(simFileName, ...
            'theConeMosaic', ...
            'theExemplarConeIndices', ...
            'thePhotocurrentImpulseResponseStruct', ...
            'responseTimeAxis', ...
            'coneMosaicPhotocurrentSpatiotemporalActivation', ...
            'coneMosaicNoisyPhotocurrentSpatiotemporalActivation', ...
            'coneMosaicSpatiotemporalActivation', ...
            'coneMosaicNoisySpatiotemporalActivation');

        % Generate videos of responses
        generateResponseVideos(...
            responseTimeAxis, ...
            coneMosaicSpatiotemporalActivation, ...
            coneMosaicPhotocurrentSpatiotemporalActivation, ...
            coneMosaicNoisySpatiotemporalActivation, ...
            coneMosaicNoisyPhotocurrentSpatiotemporalActivation, ...
            theConeMosaic, theExemplarConeIndices, ...
            thePhotocurrentImpulseResponseStruct, ...
            fixationalEyeMovements, figureFileBase, ...
            reformatExportedAVIvideoToMP4format);

        return;
    end

    % Run the simulation from scratch

    % Mosaic parameters
    mosaicSizeDegs = 0.6*[1 1];
    mosaicEccDegs = [0 0];
    mosaicIntegrationTimeSeconds = 2/1000;
    
    % Stimulus params
    contrast = 0.75;
    spatialFrequencyCPD = 10;
    temporalFrequencyHz = 8;
    stimulationDurationTemporalCycles = 6;
    simulationDurationSeconds = stimulationDurationTemporalCycles/temporalFrequencyHz;  
    if (fixationalEyeMovements)
        nTrials = 10;
    else
        nTrials = 2;
    end

   % Override parameters if we want to go through the paces fast
    if (fastParameters)
        mosaicSizeDegs = 0.4*[1 1];
        mosaicIntegrationTimeSeconds = 15/1000;
        nTrials = 1;
        stimulationDurationTemporalCycles = 4;
        simulationDurationSeconds = stimulationDurationTemporalCycles/temporalFrequencyHz;  
    end

    % Generate scene sequence representing a polar grating
    sceneComputeFunction = @sceGrating;
    customSceneParams = sceneComputeFunction();

    % Spatial params
    customSceneParams.fovDegs = max(mosaicSizeDegs)*1.2;
    customSceneParams.spatialModulationDomain = 'polar';
    customSceneParams.temporalModulation = 'drifted';
    customSceneParams.coneContrastModulation = contrast * [1 1 1];
    customSceneParams.spatialFrequencyCyclesPerDeg = spatialFrequencyCPD;
    customSceneParams.spatialEnvelopeRadiusDegs = max(mosaicSizeDegs)/3;

    % Temporal params
    customSceneParams.frameDurationSeconds = mosaicIntegrationTimeSeconds;
    customSceneParams.temporalModulationParams.phaseDirection = -1;
    customSceneParams.temporalModulationParams.temporalFrequencyHz = temporalFrequencyHz;
    customSceneParams.temporalModulationParams.stimDurationFramesNum = round(simulationDurationSeconds/customSceneParams.frameDurationSeconds);
   
    % Instantiate the scene engine
    thePolarGratingEngine = sceneEngine(sceneComputeFunction, customSceneParams);

    % Generate the null stimulus
    testContrast = 0.0;
    theNullSceneSequence = thePolarGratingEngine.compute(testContrast);

    % Generate the test stimulus
    testContrast = 1.0;
    [theSceneSequence, theSceneTemporalSupportSeconds] = thePolarGratingEngine.compute(testContrast);

    % Visualize the test stimulus
    thePolarGratingEngine.visualizeSceneSequence(...
            theSceneSequence, theSceneTemporalSupportSeconds);

    
    % Generate and export a movie of the stimulus
    theVideoFileName0 = sprintf('%s/Stimulus',figureFileBase);
    generateStimulusMovie(theVideoFileName0, theSceneSequence, theSceneTemporalSupportSeconds, reformatExportedAVIvideoToMP4format);


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

    % Determine indices of exemplar L, M and S-cone near a desired position
    exemplarConePosDegs = [0.13 -0.06];
    theExemplarConeIndices = determineConeIndicesNearPosition(theConeMosaic, exemplarConePosDegs);
   
    % Generate optics params
    rankedSujectIDs = PolansOptics.constants.subjectRanking;
    subjectRankOrder = 3;
    testSubjectID = rankedSujectIDs(subjectRankOrder);
    subtractCentralRefraction = ...
        PolansOptics.constants.subjectRequiresCentralRefractionCorrection(testSubjectID);


    % Generate optics appropriate for the mosaic's eccentricity
    [oiEnsemble, psfEnsemble] = theConeMosaic.oiEnsembleGenerate(mosaicEccDegs, ...
            'zernikeDataBase', 'Polans2015', ...
            'subjectID', testSubjectID, ...
            'pupilDiameterMM', 3.0, ...
            'subtractCentralRefraction', subtractCentralRefraction);

    % Extract the optics
    theOI = oiEnsemble{1};

    % Extract the PSF
    thePSF = psfEnsemble{1};

    % Visualize the cone mosaic and the PSF of the employed optics
    visualizeConeMosaicAndOptics(theConeMosaic, thePSF, figureFileBase);

    % Compute the sequence of optical images corresponding to the polar grating
    fprintf('Computing the optical image sequences (null stimulus + test stimulus)');
    framesNum = numel(theSceneSequence);
    theListOfOpticalImages = cell(1, framesNum);
    listOfNullOpticalImages = cell(1, framesNum);
    for frame = 1:framesNum
        theListOfOpticalImages{frame} = oiCompute(theOI, theSceneSequence{frame}, 'pad value', 'mean');
        listOfNullOpticalImages{frame} = oiCompute(theOI, theNullSceneSequence{frame}, 'pad value', 'mean');
    end

    % Generate an @oiSequence object from the list of computed optical images
    theOIsequence = oiArbitrarySequence(theListOfOpticalImages, theSceneTemporalSupportSeconds);
    nullOIsequence = oiArbitrarySequence(listOfNullOpticalImages, theSceneTemporalSupportSeconds);

    % Compute theConeMosaicNullResponse, i.e., the input cone mosaic response to the NULL scene
    coneMosaicNullResponse = theConeMosaic.compute(...
                            nullOIsequence, ...
                            'nTrials', 1);
    
    % Compute the photocurrent impulse response function
    % The contrast of the impulse stimulus
    theConeFlashImpulseContrast = 0.01*[1 1 1];

    % Whether to visualize the computed impulse responses
    visualizePhotocurrentImpulseResponses = true;

    % Compute the impulse responses (full temporal support, not just limited to the framesNum, hence the [] argument)
    thePhotocurrentImpulseResponseStruct = CMosaicNrePhotocurrentImpulseResponses(...
            theConeMosaic, coneMosaicNullResponse, theConeFlashImpulseContrast, [], ...
            visualizePhotocurrentImpulseResponses);
    thePhotocurrentImpulseResponseStruct.LMSconeImpulseResponses = thePhotocurrentImpulseResponseStruct.LMSconeImpulseResponses';

    % Generate fixational eye movement path lasting the entire stimulus duration
    if (fixationalEyeMovements)
        % Compute stimulus duration
        stimulusDurationSeconds = ...
            customSceneParams.temporalModulationParams.stimDurationFramesNum * customSceneParams.frameDurationSeconds;
        
        % Generate fEMs (fix the seed to have the same trajectories)
        theConeMosaic.emGenSequence(stimulusDurationSeconds, ...
            'microsaccadeType', 'none', ...
            'nTrials', nTrials, ...
            'randomSeed', 10);
    end

    % Compute the spatiotemporal cone-mosaic activation
    fprintf('Computing cone mosaic response\n');
    [coneMosaicSpatiotemporalActivation, coneMosaicNoisySpatiotemporalActivation, ~, ~, responseTimeAxis] = ...
        theConeMosaic.compute(theOIsequence, ...
        'withFixationalEyeMovements', fixationalEyeMovements, ...
        'nTrials', nTrials);

    if (size(coneMosaicSpatiotemporalActivation,1) == 1) && (nTrials>1)
        disp('replicating cone mosaic responses')
        coneMosaicSpatiotemporalActivation = repmat(coneMosaicSpatiotemporalActivation, [nTrials 1 1]);
    end

    % Compute the noise-free photocurrent responses by convolving the 
    % differential excitation response (response - background response) with the photocurrent impulse response
    % and adding the steady-state current
    tBinsPhotocurrentResponse = size(coneMosaicSpatiotemporalActivation,2) + numel(thePhotocurrentImpulseResponseStruct.temporalSupportSeconds) - 1;
    coneMosaicPhotocurrentSpatiotemporalActivation = zeros(nTrials, tBinsPhotocurrentResponse, theConeMosaic.conesNum);

    theConeTypes = theConeMosaic.coneTypes;
    theSteadyCurrents = thePhotocurrentImpulseResponseStruct.steadyStateCurrents;
    thePcurrentImpulseResponses = thePhotocurrentImpulseResponseStruct.LMSconeImpulseResponses;

    % Compute photocurrent responses
    parfor coneIndex = 1:theConeMosaic.conesNum

        % Retrieve the cone class
        theConeClass = theConeTypes(coneIndex);

        % Retrieve the steady-state photocurrent for this cone class
        theConeClassSpecificMeanPhotocurrent = theSteadyCurrents(theConeClass);

        % Retrieve the photocurrent impulse response function for this cone class
        theConeClassSpecificPhotocurrentImpulseResponse = thePcurrentImpulseResponses(theConeClass,:);
       
        for iTrial = 1:nTrials
            % The differential cone excitation time series (stimulus-background)
            theConeExcitationResponseTimeSeries = ...
                squeeze(coneMosaicSpatiotemporalActivation(iTrial,:, coneIndex)) - ...
                coneMosaicNullResponse(1,:, coneIndex);

            % Insert the null cone excitation differential response (which is 0) at the begining 
            theConeExcitationResponseTimeSeries = ...
                cat(2, 0*coneMosaicNullResponse(1,:, coneIndex), theConeExcitationResponseTimeSeries);

            % Convolve with the photocurrent impulse response
            thePhotocurrentResponse = ...
                conv(theConeExcitationResponseTimeSeries, theConeClassSpecificPhotocurrentImpulseResponse);

            % Retrieve the part of the response following the inserted null cone excitation response
            idx = size(coneMosaicNullResponse,2)+1 : size(thePhotocurrentResponse,2);
            coneMosaicPhotocurrentSpatiotemporalActivation(iTrial, :, coneIndex) = thePhotocurrentResponse(1,idx);

            % Add the mean (background) photocurrent
            coneMosaicPhotocurrentSpatiotemporalActivation(iTrial, :, coneIndex) = ...
                theConeClassSpecificMeanPhotocurrent + ...
                coneMosaicPhotocurrentSpatiotemporalActivation(iTrial, :, coneIndex);
        end % iTrial
    end  % parfor mCone
    
    % Keep the part of the photocurrent response that extends over the temporal extent of the cone mosaic excitation 
    tBins = 1:size(coneMosaicSpatiotemporalActivation,2);
    coneMosaicPhotocurrentSpatiotemporalActivation = coneMosaicPhotocurrentSpatiotemporalActivation(:,tBins,:);
    
    % osAddNoise() assumes that the noise-free photocurrent signal, here coneMosaicPhotocurrentSpatiotemporalActivation,
    % is arranged in a 3D matrix with the following format: [row, col, time]
    % So we have to permute the coneMosaicPhotocurrentSpatiotemporalActivation which is in [trial, time, mCones] format

    % Move the time dimension to the 3rd place
    osNoiseFormatNoiseFreeSpatiotemporalPhotocurrentSignal = ...
        permute(coneMosaicPhotocurrentSpatiotemporalActivation, [1 3 2]);

    % Add photocurrent noise
    osNoiseFormatNoisySpatiotemporalPhotocurrentSignal = ...
        osAddNoise(osNoiseFormatNoiseFreeSpatiotemporalPhotocurrentSignal, ...
            'sampTime', theConeMosaic.integrationTime);

    % Move the time dimension back to the 2nd place
    coneMosaicNoisyPhotocurrentSpatiotemporalActivation = ...
        permute(osNoiseFormatNoisySpatiotemporalPhotocurrentSignal, [1 3 2]);

    % Save the results
    save(simFileName, ...
        'theConeMosaic', ...
        'theExemplarConeIndices', ...
        'thePhotocurrentImpulseResponseStruct', ...
        'responseTimeAxis', ...
        'coneMosaicPhotocurrentSpatiotemporalActivation', ...
        'coneMosaicNoisyPhotocurrentSpatiotemporalActivation', ...
        'coneMosaicSpatiotemporalActivation', ...
        'coneMosaicNoisySpatiotemporalActivation', ...
        '-v7.3');

    % Generate videos of responses
    generateResponseVideos(...
            responseTimeAxis, ...
            coneMosaicSpatiotemporalActivation, ...
            coneMosaicPhotocurrentSpatiotemporalActivation, ...
            coneMosaicNoisySpatiotemporalActivation, ...
            coneMosaicNoisyPhotocurrentSpatiotemporalActivation, ...
            theConeMosaic, theExemplarConeIndices, ...
            thePhotocurrentImpulseResponseStruct, ...
            fixationalEyeMovements, figureFileBase, ...
            reformatExportedAVIvideoToMP4format);
end

%%  VISUALIZATION ROUTINES
function generateResponseVideos(responseTimeAxis, ...
            coneMosaicSpatiotemporalActivation, ...
            coneMosaicPhotocurrentSpatiotemporalActivation, ...
            coneMosaicNoisySpatiotemporalActivation, ...
            coneMosaicNoisyPhotocurrentSpatiotemporalActivation, ...
            theConeMosaic, theExemplarConeIndices, ...
            thePhotocurrentImpulseResponseStruct, ...
            fixationalEyeMovements, figureFileBase, ...
            reformatExportedAVIvideoToMP4format)


    % Plots of response time-series for an exemplar L-, M- and S-cone
    theLconeIndex = theExemplarConeIndices(1);
    theMconeIndex = theExemplarConeIndices(2);
    theSconeIndex = theExemplarConeIndices(3);

    theResponses = coneMosaicNoisySpatiotemporalActivation(:,:,[theLconeIndex theMconeIndex theSconeIndex])/theConeMosaic.integrationTime;
    coneExcitationResponseRange = [0 max(theResponses(:))];
   
    % theResponses = coneMosaicNoisyPhotocurrentSpatiotemporalActivation(:,:,[theLconeIndex theMconeIndex theSconeIndex]);
    pCurrentResponseRange = [-90 -40];  % pAmps
    
    fontSize = 20;

    nTrials = size(coneMosaicSpatiotemporalActivation,1);
    for iTrial = 1:nTrials
        hFig = figure(1000+iTrial); clf;
        set(hFig, 'Position', [10 10 750 1150], 'Color', [1 1 1]);

        ax = axes('Position', [0.11 0.58 0.86 0.40]);
        
        % The noisy excitations
        bar(ax,responseTimeAxis*1e3, squeeze(coneMosaicNoisySpatiotemporalActivation(:,:,theLconeIndex))/theConeMosaic.integrationTime,  1, 'EdgeColor', [1 0 0],   'FaceColor', [1 0.5 0.5],   'FaceAlpha', 0.5, 'EdgeAlpha', 0.5, 'LineWidth', 0.5);
        hold(ax, 'on')
        bar(ax,responseTimeAxis*1e3, squeeze(coneMosaicNoisySpatiotemporalActivation(:,:,theMconeIndex))/theConeMosaic.integrationTime, 1, 'EdgeColor', [0 0.8 0], 'FaceColor', [0.5 0.8 0.5], 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5, 'LineWidth', 0.5);
        bar(ax,responseTimeAxis*1e3, squeeze(coneMosaicNoisySpatiotemporalActivation(:,:,theSconeIndex))/theConeMosaic.integrationTime, 1, 'EdgeColor', [0 0 1],   'FaceColor', [0.5 0.5 1],   'FaceAlpha', 0.5, 'EdgeAlpha', 0.5, 'LineWidth', 0.5);
        
        % The noise-free excitations
        plot(ax,responseTimeAxis*1e3, squeeze(coneMosaicSpatiotemporalActivation(iTrial,:,theLconeIndex))/theConeMosaic.integrationTime, '-', 'Color', 'k', 'LineWidth', 4);
        plot(ax,responseTimeAxis*1e3, squeeze(coneMosaicSpatiotemporalActivation(iTrial,:,theMconeIndex))/theConeMosaic.integrationTime, '-', 'Color', 'k', 'LineWidth', 4);
        plot(ax,responseTimeAxis*1e3, squeeze(coneMosaicSpatiotemporalActivation(iTrial,:,theSconeIndex))/theConeMosaic.integrationTime, '-', 'Color', 'k', 'LineWidth', 4);
       
        p1 = plot(ax,responseTimeAxis*1e3, squeeze(coneMosaicSpatiotemporalActivation(iTrial,:,theLconeIndex))/theConeMosaic.integrationTime, '-', 'Color', 'r', 'LineWidth', 1.5);
        p2 = plot(ax,responseTimeAxis*1e3, squeeze(coneMosaicSpatiotemporalActivation(iTrial,:,theMconeIndex))/theConeMosaic.integrationTime, '-', 'Color', [0 0.8 0], 'LineWidth', 1.5);
        p3 = plot(ax,responseTimeAxis*1e3, squeeze(coneMosaicSpatiotemporalActivation(iTrial,:,theSconeIndex))/theConeMosaic.integrationTime, '-', 'Color', 'b', 'LineWidth', 1.5);
       
        legend(ax, [p1 p2 p3], {'L-cone', 'M-cone', 'S-cone'}, 'Location', 'NorthWest');
        set(ax, 'XLim', [responseTimeAxis(1) responseTimeAxis(end)]*1e3, 'FontSize', fontSize);
        set(ax, 'YTick', [0 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000]);
        set(ax, 'YTickLabel', {'0', '1k', '2k', '3k', '4k', '5k', '6k', '7k', '8k', '9k', '10k'});
        set(ax, 'YLim', coneExcitationResponseRange);
        xlabel(ax,'time (msec)');
        box(ax, 'off');
        grid(ax, 'on');
        ylabel(ax,'cone excitations (isomerizations/sec)');

        ax = axes('Position', [0.11 0.05 0.86 0.40]);
        dT = responseTimeAxis(2)-responseTimeAxis(1);
        coneMosaicPhotocurrentResponseTimeAxis = (0:(size(coneMosaicPhotocurrentSpatiotemporalActivation,2)-1))*dT;
        
        % The noisy photocurrents
        plot(ax,coneMosaicPhotocurrentResponseTimeAxis*1e3, squeeze(coneMosaicNoisyPhotocurrentSpatiotemporalActivation(iTrial,:,theLconeIndex)), '-', 'Color', 'r', 'LineWidth', 1);
        hold(ax, 'on')
        plot(ax,coneMosaicPhotocurrentResponseTimeAxis*1e3, squeeze(coneMosaicNoisyPhotocurrentSpatiotemporalActivation(iTrial,:,theMconeIndex)), '-', 'Color', 'g', 'LineWidth', 1);
        plot(ax,coneMosaicPhotocurrentResponseTimeAxis*1e3, squeeze(coneMosaicNoisyPhotocurrentSpatiotemporalActivation(iTrial,:,theSconeIndex)), '-', 'Color', 'b', 'LineWidth', 1);
        
        % The noise-free photocurrents
        plot(ax,coneMosaicPhotocurrentResponseTimeAxis*1e3, squeeze(coneMosaicPhotocurrentSpatiotemporalActivation(iTrial,:,theLconeIndex)), '-', 'Color', 'k', 'LineWidth', 4);
        plot(ax,coneMosaicPhotocurrentResponseTimeAxis*1e3, squeeze(coneMosaicPhotocurrentSpatiotemporalActivation(iTrial,:,theMconeIndex)), '-', 'Color', 'k', 'LineWidth', 4);
        plot(ax,coneMosaicPhotocurrentResponseTimeAxis*1e3, squeeze(coneMosaicPhotocurrentSpatiotemporalActivation(iTrial,:,theSconeIndex)), '-', 'Color', 'k', 'LineWidth', 4);
        
        p1 = plot(ax,coneMosaicPhotocurrentResponseTimeAxis*1e3, squeeze(coneMosaicPhotocurrentSpatiotemporalActivation(iTrial,:,theLconeIndex)), '-', 'Color', 'r', 'LineWidth', 1.5);
        p2 = plot(ax,coneMosaicPhotocurrentResponseTimeAxis*1e3, squeeze(coneMosaicPhotocurrentSpatiotemporalActivation(iTrial,:,theMconeIndex)), '-', 'Color', 'g', 'LineWidth', 1.5);
        p3 = plot(ax,coneMosaicPhotocurrentResponseTimeAxis*1e3, squeeze(coneMosaicPhotocurrentSpatiotemporalActivation(iTrial,:,theSconeIndex)), '-', 'Color', 'b', 'LineWidth', 1.5);
        legend(ax, [p1 p2 p3], {'L-cone', 'M-cone', 'S-cone'}, 'Location', 'NorthWest');
        xlabel(ax,'time (msec)');
        ylabel(ax, '');
        %set(ax, 'YTickLabel', {});
        box(ax, 'off');
        grid(ax, 'on');
        set(ax, 'XLim', [responseTimeAxis(1) responseTimeAxis(end)]*1e3, 'FontSize', fontSize);
        set(ax, 'YLim', pCurrentResponseRange);
        ylabel(ax, 'photocurrent (pAmps)');

        % Plot the pCurrent impulse response
        ax = axes('Position', [0.72 0.43 0.25 0.1]);
        plot(ax, thePhotocurrentImpulseResponseStruct.temporalSupportSeconds*1e3, thePhotocurrentImpulseResponseStruct.LMSconeImpulseResponses(cMosaic.LCONE_ID,:), 'k-', 'LineWidth', 4);
        hold on;
        p1 = plot(ax, thePhotocurrentImpulseResponseStruct.temporalSupportSeconds*1e3, thePhotocurrentImpulseResponseStruct.LMSconeImpulseResponses(cMosaic.LCONE_ID,:), 'r-', 'LineWidth', 1.5);
        plot(ax, thePhotocurrentImpulseResponseStruct.temporalSupportSeconds*1e3, thePhotocurrentImpulseResponseStruct.LMSconeImpulseResponses(cMosaic.MCONE_ID,:), 'k-', 'LineWidth', 4);
        p2 = plot(ax, thePhotocurrentImpulseResponseStruct.temporalSupportSeconds*1e3, thePhotocurrentImpulseResponseStruct.LMSconeImpulseResponses(cMosaic.MCONE_ID,:), 'g-', 'LineWidth', 1.5);
        plot(ax, thePhotocurrentImpulseResponseStruct.temporalSupportSeconds*1e3, thePhotocurrentImpulseResponseStruct.LMSconeImpulseResponses(cMosaic.SCONE_ID,:), 'b-', 'LineWidth', 4);
        p3 = plot(ax, thePhotocurrentImpulseResponseStruct.temporalSupportSeconds*1e3, thePhotocurrentImpulseResponseStruct.LMSconeImpulseResponses(cMosaic.SCONE_ID,:), 'b-', 'LineWidth', 1.5);
        ylabel(ax, 'pCurrent IR');
        legend(ax, [p1 p2 p3], {'L-cone', 'M-cone', 'S-cone'});
        set(ax, 'fontSize', 12);
        box(ax, 'off');
        grid(ax, 'on')
        set(ax, 'YTickLabel', {}, 'XLim', [0 500], 'XTick', 0:100:1000);

        NicePlot.exportFigToPDF(sprintf('%s/responseTimeSeriesTrial%d.pdf',figureFileBase, iTrial), hFig, 300);
    end

    % Videos of cone mosaic spatiotemporal activatio
    colorbarTickLabelColor = [0.1 0.1 0.1];
    backgroundColor = [0.1 0.1 0.1];
    pCurrentResponseRange = [-61 -49];            % outer segment current, pAmps
    coneExcitationsResponseRange = [1000 3500];   % isomerizations/sec (R*/sec)

    domainvisualizationlimits(1) = theConeMosaic.eccentricityDegs(1) - 0.51*theConeMosaic.sizeDegs(1);
    domainvisualizationlimits(2) = theConeMosaic.eccentricityDegs(1) + 0.51*theConeMosaic.sizeDegs(1);
    domainvisualizationlimits(3) = theConeMosaic.eccentricityDegs(2) - 0.51*theConeMosaic.sizeDegs(2);
    domainvisualizationlimits(4) = theConeMosaic.eccentricityDegs(2) + 0.51*theConeMosaic.sizeDegs(2);

    for iTrial = 1:size(coneMosaicNoisySpatiotemporalActivation,1)

        theVideoFileName1 = sprintf('%s/NoiseFreeExcitationsTrial%d',figureFileBase, iTrial);
	    %videoOBJ1 = VideoWriter(theVideoFileName1, 'M PEG-4');  % H264format (has artifacts)
        videoOBJ1 = VideoWriter(theVideoFileName1, 'Uncompressed AVI');
	    videoOBJ1.FrameRate = 30;
	    %videoOBJ1.Quality = 100;
	    videoOBJ1.open();

        theVideoFileName2 = sprintf('%s/NoiseFreePhotocurrentsTrial%d',figureFileBase, iTrial);
	    % videoOBJ2 = VideoWriter(theVideoFileName2, 'MPEG-4'); % H264format (has artifacts)
        videoOBJ2 = VideoWriter(theVideoFileName2, 'Uncompressed AVI');
	    videoOBJ2.FrameRate = 30;
	    %videoOBJ2.Quality = 100;
	    videoOBJ2.open();
    
        theVideoFileName3 = sprintf('%s/NoisyExcitationsTrial%d',figureFileBase, iTrial);
	    % videoOBJ3 = VideoWriter(theVideoFileName3, 'MPEG-4'); % H264format (has artifacts)
        videoOBJ3 = VideoWriter(theVideoFileName3, 'Uncompressed AVI');
	    videoOBJ3.FrameRate = 30;
	    %videoOBJ3.Quality = 100;
	    videoOBJ3.open();

        theVideoFileName4 = sprintf('%s/NoisyPhotocurrentsTrial%d',figureFileBase, iTrial);
	    % videoOBJ4 = VideoWriter(theVideoFileName3, 'MPEG-4'); % H264format (has artifacts)
        videoOBJ4 = VideoWriter(theVideoFileName4, 'Uncompressed AVI');
	    videoOBJ4.FrameRate = 30;
	    %videoOBJ4.Quality = 100;
	    videoOBJ4.open();


        hFig1 = figure(21); clf;
        set(hFig1, 'Position', [10 10 600 650], 'Color', [1 1 1]);
        ax1 = subplot('Position', [0.09 0.09 0.91 0.89]);
    
        hFig2 = figure(22); clf;
        set(hFig2, 'Position', [10 10 600 650], 'Color', [1 1 1]);
        ax2 = subplot('Position', [0.09 0.09 0.91 0.89]);
    
        hFig3 = figure(23); clf;
        set(hFig3, 'Position', [10 10 600 650], 'Color', [1 1 1]);
        ax3 = subplot('Position', [0.09 0.09 0.91 0.89]);
    
        hFig4 = figure(24); clf;
        set(hFig4, 'Position', [10 10 600 650], 'Color', [1 1 1]);
        ax4 = subplot('Position', [0.09 0.09 0.91 0.89]);

    
        for t = 1:size(coneMosaicNoisySpatiotemporalActivation,2)
            if (fixationalEyeMovements)
                % Visualize cone mosaic response with fEM
                theConeMosaic.visualize('figureHandle', hFig1, 'axesHandle', ax1, ...
                    'activation', coneMosaicSpatiotemporalActivation(1,t,:)/theConeMosaic.integrationTime, ...
                    'activationRange', coneExcitationsResponseRange, ...
                    'domainvisualizationlimits', domainvisualizationlimits, ...
                    'labelConesWithIndices', theExemplarConeIndices, ...
                    'displayedEyeMovementData', struct(...
                        'trial', iTrial, ...
                        'timePoints', 1:t), ...
                    'horizontalActivationColorbar', true, ...
                    'colorbarTickLabelColor', colorbarTickLabelColor, ...
                    'backgroundColor', backgroundColor, ...
                    'fontSize', fontSize, ...
                    'plotTitle', sprintf('noise-free excitations (R*/sec)\nt = %2.0f msec', 1000*responseTimeAxis(t)));
            else
                 % Visualize cone mosaic response without fEM
                theConeMosaic.visualize('figureHandle', hFig1, 'axesHandle', ax1, ...
                    'activation', coneMosaicSpatiotemporalActivation(1,t,:)/theConeMosaic.integrationTime, ...
                    'activationRange', coneExcitationsResponseRange, ...
                    'domainvisualizationlimits', domainvisualizationlimits, ...
                    'labelConesWithIndices', theExemplarConeIndices, ...
                    'horizontalActivationColorbar', true, ...
                    'colorbarTickLabelColor', colorbarTickLabelColor, ...
                    'backgroundColor', backgroundColor, ...
                    'fontSize', fontSize, ...
                    'plotTitle', sprintf('noise-free excitations (R*/sec)\nt = %2.0f msec', 1000*responseTimeAxis(t)));
 
            end

            drawnow;
		    videoOBJ1.writeVideo(getframe(hFig1));

            if (fixationalEyeMovements)
                % Visualize cone mosaic response with fEM
                theConeMosaic.visualize('figureHandle', hFig2, 'axesHandle', ax2, ...
                    'activation', squeeze(coneMosaicPhotocurrentSpatiotemporalActivation(iTrial,t,:)), ...
                    'activationRange', pCurrentResponseRange, ...
                    'domainvisualizationlimits', domainvisualizationlimits, ...
                    'labelConesWithIndices', theExemplarConeIndices, ...
                    'displayedEyeMovementData', struct(...
                        'trial', iTrial, ...
                        'timePoints', 1:t), ...
                    'horizontalActivationColorbar', true, ...
                    'colorbarTickLabelColor', colorbarTickLabelColor, ...
                    'backgroundColor', backgroundColor, ...
                    'fontSize', fontSize, ...
                    'plotTitle', sprintf('noise-free photocurrents (pAmps)\nt = %2.0f msec', 1000*responseTimeAxis(t)));
            else
                % Visualize cone mosaic response without fEM
                theConeMosaic.visualize('figureHandle', hFig2, 'axesHandle', ax2, ...
                    'activation', squeeze(coneMosaicPhotocurrentSpatiotemporalActivation(iTrial,t,:)), ...
                    'activationRange', pCurrentResponseRange, ...
                    'domainvisualizationlimits', domainvisualizationlimits, ...
                    'labelConesWithIndices', theExemplarConeIndices, ...
                    'horizontalActivationColorbar', true, ...
                    'colorbarTickLabelColor', colorbarTickLabelColor, ...
                    'backgroundColor', backgroundColor, ...
                    'fontSize', fontSize, ...
                    'plotTitle', sprintf('noise-free photocurrents (pAmps)\nt = %2.0f msec', 1000*responseTimeAxis(t)));
     
            end

            drawnow;
		    videoOBJ2.writeVideo(getframe(hFig2));

            if (fixationalEyeMovements)
                % Visualize cone mosaic response with fEM
                theConeMosaic.visualize('figureHandle', hFig3, 'axesHandle', ax3, ...
                    'activation', squeeze(coneMosaicNoisySpatiotemporalActivation(iTrial,t,:))/theConeMosaic.integrationTime, ...
                    'activationRange', coneExcitationsResponseRange, ...
                    'domainvisualizationlimits', domainvisualizationlimits, ...
                    'labelConesWithIndices', theExemplarConeIndices, ...
                    'displayedEyeMovementData', struct(...
                        'trial', iTrial, ...
                        'timePoints', 1:t), ...
                    'horizontalActivationColorbar', true, ...
                    'colorbarTickLabelColor', colorbarTickLabelColor, ...
                    'backgroundColor', backgroundColor, ...
                    'fontSize', fontSize, ...
                    'plotTitle', sprintf('noisy excitations (R*/sec)\nt = %2.0f msec', 1000*responseTimeAxis(t)));
            else
                % Visualize cone mosaic response without fEM
                theConeMosaic.visualize('figureHandle', hFig3, 'axesHandle', ax3, ...
                    'activation', squeeze(coneMosaicNoisySpatiotemporalActivation(iTrial,t,:))/theConeMosaic.integrationTime, ...
                    'activationRange', coneExcitationsResponseRange, ...
                    'domainvisualizationlimits', domainvisualizationlimits, ...
                    'labelConesWithIndices', theExemplarConeIndices, ...
                    'horizontalActivationColorbar', true, ...
                    'colorbarTickLabelColor', colorbarTickLabelColor, ...
                    'backgroundColor', backgroundColor, ...
                    'fontSize', fontSize, ...
                    'plotTitle', sprintf('noisy excitations (R*/sec)\nt = %2.0f msec', 1000*responseTimeAxis(t)));
            end

            drawnow;
		    videoOBJ3.writeVideo(getframe(hFig3));

            if (fixationalEyeMovements)
                % Visualize cone mosaic response with fEM
                theConeMosaic.visualize('figureHandle', hFig4, 'axesHandle', ax4, ...
                    'activation', squeeze(coneMosaicNoisyPhotocurrentSpatiotemporalActivation(iTrial,t,:)), ...
                    'activationRange', pCurrentResponseRange, ...
                    'domainvisualizationlimits', domainvisualizationlimits, ...
                    'labelConesWithIndices', theExemplarConeIndices, ...
                    'displayedEyeMovementData', struct(...
                        'trial', iTrial, ...
                        'timePoints', 1:t), ...
                    'horizontalActivationColorbar', true, ...
                    'colorbarTickLabelColor', colorbarTickLabelColor, ...
                    'backgroundColor', backgroundColor, ...
                    'fontSize', fontSize, ...
                    'plotTitle', sprintf('noisy photocurrents (pAmps)\nt = %2.0f msec', 1000*responseTimeAxis(t)));
            else
                % Visualize cone mosaic response without fEM
                theConeMosaic.visualize('figureHandle', hFig4, 'axesHandle', ax4, ...
                    'activation', squeeze(coneMosaicNoisyPhotocurrentSpatiotemporalActivation(iTrial,t,:)), ...
                    'activationRange', pCurrentResponseRange, ...
                    'domainvisualizationlimits', domainvisualizationlimits, ...
                    'labelConesWithIndices', theExemplarConeIndices, ...
                    'horizontalActivationColorbar', true, ...
                    'colorbarTickLabelColor', colorbarTickLabelColor, ...
                    'backgroundColor', backgroundColor, ...
                    'fontSize', fontSize, ...
                    'plotTitle', sprintf('noisy photocurrents (pAmps)\nt = %2.0f msec', 1000*responseTimeAxis(t)));
     
            end

            drawnow;
		    videoOBJ4.writeVideo(getframe(hFig4));
        end
 
        videoOBJ1.close;
        videoOBJ2.close;
        videoOBJ3.close;
        videoOBJ4.close;
    
        % Reformat video to AVI 
        if (reformatExportedAVIvideoToMP4format)
            reformatVideo(theVideoFileName1);
            reformatVideo(theVideoFileName2);
            reformatVideo(theVideoFileName3);
            reformatVideo(theVideoFileName4);
        end

    end % iTrial
end

function generateStimulusMovie(theVideoFileName, theSceneSequence, theSceneTemporalSupportSeconds, reformatExportedAVIvideoToMP4format)
    % Generate stimulus video
    videoOBJ0 = VideoWriter(theVideoFileName, 'Uncompressed AVI');
	videoOBJ0.FrameRate = 30;
	videoOBJ0.open();

    hFig = figure(20);
    set(hFig, 'Position', [10 10 600 650]);
    ax = subplot('Position', [0.02 0.02 0.96 0.96]);
    
    for frame = 1:numel(theSceneSequence)
        image(ax, sceneGet(theSceneSequence{frame}, 'rgbimage'));
        axis(ax, 'image');
        title(ax, sprintf('t = %2.0f msec', theSceneTemporalSupportSeconds(frame)*1e3));
        set(ax, 'XTick', [], 'YTick', [], 'FontSize', 16);
        drawnow;
        videoOBJ0.writeVideo(getframe(hFig));
    end
    videoOBJ0.close();

    if (reformatExportedAVIvideoToMP4format)
        reformatVideo(theVideoFileName);
    end
end

function reformatVideo(videoFileName)
    if (exist('/opt/homebrew/bin/ffmpeg'))
        sysCommand = sprintf('/opt/homebrew/bin/ffmpeg -i %s.avi -c:v libx264 -crf 22 -pix_fmt yuv420p %s.mp4', ...
            videoFileName, videoFileName);
        system(sysCommand);
    
        sysCommand = sprintf('rm %s.avi ', videoFileName);
        system(sysCommand);
    end
end

function theExemplarConeIndices = determineConeIndicesNearPosition(theConeMosaic, exemplarConePosDegs)
    
    exemplarConeDistances = sqrt(sum((bsxfun(@minus, theConeMosaic.coneRFpositionsDegs, exemplarConePosDegs)).^2,2));
    [~,idx] = sort(exemplarConeDistances, 'ascend');
    theLconeIndex = [];
    theMconeIndex = [];
    theSconeIndex = [];

    i = 0;
    while ((isempty(theLconeIndex)) || (isempty(theMconeIndex)) || (isempty(theSconeIndex))) && (i <= theConeMosaic.conesNum)
        i = i + 1;
        theConeIndex = idx(i);
        theConeType = theConeMosaic.coneTypes(theConeIndex);
        switch (theConeType)
            case cMosaic.LCONE_ID
                if (isempty(theLconeIndex))
                    theLconeIndex = theConeIndex;
                    fprintf('Found exemplar Lcone\n')
                end
            case cMosaic.MCONE_ID
                if (isempty(theMconeIndex))
                    theMconeIndex = theConeIndex;
                    fprintf('Found exemplar Mcone\n')
                end
            case cMosaic.SCONE_ID
                if (isempty(theSconeIndex))
                    theSconeIndex = theConeIndex;
                    fprintf('Found exemplar Scone\n')
                end
        end % switch
    end % while
    
    theExemplarConeIndices = [theLconeIndex theMconeIndex theSconeIndex];
end

function visualizeConeMosaicAndOptics(theConeMosaic, thePSF, figureFileBase)

    % Extract the PSF at 550 nm
    targetWavelength = 550;
    [~,idx] = min(abs(thePSF.supportWavelength-targetWavelength));
    thePSFstuct = struct(...
         'supportXdegs', thePSF.supportX/60, ...
         'supportYdegs', thePSF.supportY/60, ...
         'data', thePSF.data(:,:,idx));

    % Visualize the full cone mosaic
    hFig0 = figure(20); clf;
    set(hFig0, 'Position', [10 10 600 650], 'Color', [1 1 1]);
    ax = subplot('Position', [0.1 0.06 0.89 0.89]);
    theConeMosaic.visualize(...
        'figureHandle', hFig0, ...
        'axesHandle', ax, ...
        'withSuperimposedPSF', thePSFstuct, ...
        'visualizedConeAperture', 'lightCollectingArea4sigma', ...
        'domainVisualizationTicks', struct('x', -0.2:0.1:0.2, 'y', -0.2:0.1:0.2), ...
        'domainVisualizationLimits', 0.16*[-1 1 -1 1], ...
        'fontSize', 20, ...
        'plotTitle', 'central mosaic + PSF');
    NicePlot.exportFigToPDF(sprintf('%s/centralRegionOfConeMosaicWithOpticalPSF.pdf',figureFileBase), hFig0, 300);

    % Visualize the central part of the cone mosaic with the PSF at 550 nm
    hFig0 = figure(20); clf;
    set(hFig0, 'Position', [10 10 600 650], 'Color', [1 1 1]);
    ax = subplot('Position', [0.1 0.06 0.89 0.89]);
    theConeMosaic.visualize(...
        'figureHandle', hFig0, ...
        'domainVisualizationTicks', struct('x', -1:0.2:1, 'y', -1:0.2:1), ...
        'axesHandle', ax, ...
        'fontSize', 20, ...
        'plotTitle', 'full mosaic');
    NicePlot.exportFigToPDF(sprintf('%s/fullConeMosaic.pdf',figureFileBase), hFig0, 300);
end
