function thePhotocurrentImpulseResponseStruct = CMosaicNrePhotocurrentImpulseResponses(...
    theConeMosaic, coneMosaicNullResponse, theConeFlashImpulseContrast, temporalSamplesNum, ...
    visualizePhotocurrentImpulseResponses)
% Function for computing cone outer segment photocurrent impulse responses
%
% Syntax:
%   thePhotocurrentImpulseResponseStruct = CMosaicNrePhotocurrentImpulseResponses(...
%       theConeMosaic, coneMosaicNullResponse, theConeFlashImpulseContrast, temporalSamplesNum, ...
%       visualizePhotocurrentImpulseResponses)
%
% Description:
%   Compute cone outersegment photocurrent impulse response (pCurrent IR)
%   functions at a temporal resolution equal to theConeMosaic.integrationTime 
%   and a duration equal to temporalSamplesNum x theConeMosaic.integrationTime
%
% Inputs:
%   - theConeMosaic                         - the @cMosaic for which to generate outer segment pCurrent IR functions
%
%   - coneMosaicNullResponse                - a 3D matrix, [kTrials x tFrames x nCones], of cone responses to the background stimulus (zero contrast)
%                                            (kTrials and tFrames can both be equal to 1)
%
%   - theConeFlashImpulseContrast           - a 3-element vector with values in (0..1) specifying the 
%                                             L-, M-, and S-cone contrast of the stimulus impulse for which to compute the 
%                                             L-cone, M-cone and S-cone pCurrent IR functions
%
%   - temporalSamplesNum                    - the number of temporal samples in the returned pCurrent IRs, 
%                                             in which each sample is equal to theConeMosaic.integrationTime
%
%   - visualizePhotocurrentImpulseResponses - boolean, whether to visualize the computed pCurrent IRs
% 
% Outputs:
%   - thePhotocurrentImpulseResponseStruct - struct with the following fields
%       - LMSconeIncrementImpulseResponses                - [temporalSamplesNum, 3] matrix of L-, M- and S-cone impulse increment responses
%       - LMSconeDecrementImpulseResponses                - [temporalSamplesNum, 3] matrix of L-, M- and S-cone impulse decrement responses
%       - LMSconeImpulseResponses                         - [temporalSamplesNum, 3] matrix of L-, M- and S-cone impulse responses = 0.5*(LMSconeIncrementImpulseResponses-LMSconeDecrementImpulseResponses);
%       - steadyStateCurrents                             - [1 3] vector of steady-stage (background) L-, M-, and S-cone photocurrents
%       - coneDensityWeightedPhotocurrentImpulseResponse  - [temporalSamplesNum, 1] L-, M-, S-cone density weighted pCurrent IR function
%       - temporalSupportSeconds                          - [temporalSamplesNum, 1] vector containing the temporal support of the computed pCurrent IRs
%
%
% Example usage:
% 
%  coneMosaicNullResponse = ...
%    theConeMosaic.compute(nullOIsequence, ...
%                          'nTrials', 1);
%    
%   theConeFlashImpulseContrast = 0.01*[1 1 1];
%   visualizePhotocurrentImpulseResponses = true;
%   temporalSamplesNum = [];  % full duration
%
%    thePhotocurrentImpulseResponseStruct = CMosaicNrePhotocurrentImpulseResponses(...
%            theConeMosaic, coneMosaicNullResponse, theConeFlashImpulseContrast, temporalSamplesNum, ...
%            visualizePhotocurrentImpulseResponses);
%
%
% For more usage examples, see t_onTheFlyPhotocurrentComputationWithEyeMovements
%

% History:
%    02/21/2025  NPC  Wrote it


    % Get the dimensions of the null response
    [instancesNum, timeBinsNum, mConesNum] = size(coneMosaicNullResponse);
    assert(mConesNum == theConeMosaic.conesNum, 'cMosaicResponse dimensionality error');

    theLconeExcitations = coneMosaicNullResponse(1:instancesNum, 1:timeBinsNum, theConeMosaic.lConeIndices);
    theMconeExcitations = coneMosaicNullResponse(1:instancesNum, 1:timeBinsNum, theConeMosaic.mConeIndices);
    theSconeExcitations = coneMosaicNullResponse(1:instancesNum, 1:timeBinsNum, theConeMosaic.sConeIndices);

    % Compute the mean cone excitations for each cone class based on the null response
    meanConeExcitationsPerIntegrationTime(cMosaic.LCONE_ID) = mean(theLconeExcitations(:));
    meanConeExcitationsPerIntegrationTime(cMosaic.MCONE_ID) = mean(theMconeExcitations(:));
    meanConeExcitationsPerIntegrationTime(cMosaic.SCONE_ID) = mean(theSconeExcitations(:));

    % Compute the background photon rate (R*/sec) for each cone class
    theBackgroundConeExcitationRates = meanConeExcitationsPerIntegrationTime / theConeMosaic.integrationTime;

    % Compute impulse responses for a flash duration that equals the cone mosaic integration time
    % theFlashDurationSeconds = theConeMosaic.integrationTime;

    % Or compute impulse responses for a very short flash (0.1milliseconds)
    theFlashDurationSeconds = 0.1/1000;

    [LMSincrementImpulseResponses, LMSdecrementImpulseResponses, temporalSupport, steadyStateCurrents] = ...
        computePhotocurrentResponses(theBackgroundConeExcitationRates, theConeFlashImpulseContrast,  ...
        theFlashDurationSeconds, sqrt(sum(theConeMosaic.eccentricityDegs.^2,2)));

    % Decimate at the integrationTime ensuring with sample the peak response
    downsampleMethod = 'linearInterpolationPeakTimeAligned';

    for iCone = cMosaic.LCONE_ID:cMosaic.SCONE_ID
        [decimatedLMSincrementImpulseResponses(:,iCone), decimatedTemporalSupport] = ...
            resampleResponse(temporalSupport, LMSincrementImpulseResponses(:,iCone), theConeMosaic.integrationTime, downsampleMethod);
        decimatedLMSdecrementImpulseResponses(:,iCone) = ...
            resampleResponse(temporalSupport, LMSdecrementImpulseResponses(:,iCone), theConeMosaic.integrationTime, downsampleMethod);
    end

    if (visualizePhotocurrentImpulseResponses)
        hFig = visualizeResponses(100, theFlashDurationSeconds, ...
            temporalSupport, LMSincrementImpulseResponses, LMSdecrementImpulseResponses, ...
            decimatedTemporalSupport, decimatedLMSincrementImpulseResponses, decimatedLMSdecrementImpulseResponses);
         NicePlot.exportFigToPDF('100.pdf', hFig, 300);
    end

    % Compute cone density weighted impulse response
    coneDensityWeightedPhotocurrentImpulseResponse = ...
        theConeMosaic.achievedConeDensities(cMosaic.LCONE_ID) * decimatedLMSincrementImpulseResponses(:, cMosaic.LCONE_ID) + ...
        theConeMosaic.achievedConeDensities(cMosaic.MCONE_ID) * decimatedLMSincrementImpulseResponses(:, cMosaic.MCONE_ID) + ...
        theConeMosaic.achievedConeDensities(cMosaic.SCONE_ID) * decimatedLMSincrementImpulseResponses(:, cMosaic.SCONE_ID);

    % Only return the requested frames num
    if (~isempty(temporalSamplesNum))
        coneDensityWeightedPhotocurrentImpulseResponse = coneDensityWeightedPhotocurrentImpulseResponse(1:temporalSamplesNum);
        decimatedTemporalSupport = decimatedTemporalSupport(1:temporalSamplesNum);
    end


    % Fill the returned struct
    thePhotocurrentImpulseResponseStruct.LMSconeIncrementImpulseResponses = decimatedLMSincrementImpulseResponses;
    thePhotocurrentImpulseResponseStruct.LMSconeDecrementImpulseResponses = decimatedLMSdecrementImpulseResponses;
    thePhotocurrentImpulseResponseStruct.LMSconeImpulseResponses = 0.5*(decimatedLMSincrementImpulseResponses-decimatedLMSdecrementImpulseResponses);
    thePhotocurrentImpulseResponseStruct.steadyStateCurrents = steadyStateCurrents;
    thePhotocurrentImpulseResponseStruct.coneDensityWeightedPhotocurrentImpulseResponse = coneDensityWeightedPhotocurrentImpulseResponse;
    thePhotocurrentImpulseResponseStruct.temporalSupportSeconds = decimatedTemporalSupport;

    if (visualizePhotocurrentImpulseResponses)
        figure(11); clf;
        stem(thePhotocurrentImpulseResponseStruct.temporalSupportSeconds, ...
            thePhotocurrentImpulseResponseStruct.coneDensityWeightedPhotocurrentImpulseResponse, ...
            'k-', 'filled', ...
            'LineWidth', 2.0, 'MarkerSize', 16);
        set(gca, 'XTick', thePhotocurrentImpulseResponseStruct.temporalSupportSeconds);
        set(gca, 'YLim', max(abs(thePhotocurrentImpulseResponseStruct.coneDensityWeightedPhotocurrentImpulseResponse(:))) * [-1 1]);
        grid on
        set(gca, 'FontSize', 16)
        xlabel('time (sec)');
        axis 'square'
        title(sprintf('cone density weighted\nphotocurrent impulse response'));
    end

end


function [theResampledResponse, theResampledResponseTemporalSupport] = resampleResponse(...
    temporalSupport, theResponse, resampledTimeInterval, downsampleMethod)

    theResponse = reshape(theResponse, [1 numel(theResponse)]);

    % The resampled signal time support
    theResampledResponseTemporalSupport = temporalSupport(1):resampledTimeInterval:temporalSupport(end);

    switch (downsampleMethod)
        case 'resample'
            [p,q] = rat(numel(theResampledResponseTemporalSupport) / numel(temporalSupport));
            theResampledResponse = resample(theResponse,p,q);
        
        case 'linearInterpolation'
            theResampledResponse = interp1(temporalSupport, theResponse, theResampledResponseTemporalSupport);
        
        case 'linearInterpolationPeakTimeAligned'
            % Shift the temporal support so that theResampledResponse samples the peak of the original response
            [~,idxOfPeakResponse] = max(abs(theResponse));
            theTemporalSupportShift = min(abs(theResampledResponseTemporalSupport-temporalSupport(idxOfPeakResponse)));
            [~,idx] = min(abs(temporalSupport-theTemporalSupportShift));
            theResampledResponse = interp1(temporalSupport, theResponse, theResampledResponseTemporalSupport - temporalSupport(idx));
            theResampledResponse(isnan(theResampledResponse)) = 0;

        otherwise
            error('Unknown resampling method: ''%s''.', resampleMethod)
    end % switch downSampleMethod
end


function hFig = visualizeResponses(figNo, theFlashDurationSeconds, ...
    temporalSupport, LMSincrementImpulseResponses, LMSdecrementImpulseResponses, ...
    decimatedTemporalSupport, decimatedLMSincrementImpulseResponses, decimatedLMSdecrementImpulseResponses)


    maxAmplitude = max([max(abs(LMSincrementImpulseResponses(:))) max(abs(LMSdecrementImpulseResponses(:)))]);
    amplitudeRangeOriginal = maxAmplitude * 1.05*[-1 1];

    maxAmplitude = max([max(abs(decimatedLMSincrementImpulseResponses(:))) max(abs(decimatedLMSdecrementImpulseResponses(:)))]);
    amplitudeRangeDecimated = maxAmplitude * 1.05*[-1 1];

    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 1210 1285], 'Color', [1 1 1]);
    set(hFig, 'Name', sprintf('Flash Duration (msec): %2.1f', theFlashDurationSeconds*1e3));

    wideLineWidth = 6;
    narrowLineWidth = 2;
    ax = subplot(1,2,1);
    plot(ax,temporalSupport*1e3, LMSincrementImpulseResponses(:,cMosaic.LCONE_ID), 'k-', 'LineWidth', wideLineWidth);
    hold(ax, 'on');
    p1 = plot(ax,temporalSupport*1e3, LMSincrementImpulseResponses(:,cMosaic.LCONE_ID), 'r-', 'LineWidth', narrowLineWidth);

    plot(ax,temporalSupport*1e3, LMSdecrementImpulseResponses(:,cMosaic.LCONE_ID), 'k-', 'LineWidth', wideLineWidth)
    plot(ax,temporalSupport*1e3, LMSdecrementImpulseResponses(:,cMosaic.LCONE_ID), 'r:', 'LineWidth', narrowLineWidth);

    plot(ax,temporalSupport*1e3, LMSincrementImpulseResponses(:,cMosaic.MCONE_ID), 'k-', 'LineWidth', wideLineWidth);
    p2 = plot(ax,temporalSupport*1e3, LMSincrementImpulseResponses(:,cMosaic.MCONE_ID), 'g-', 'LineWidth', narrowLineWidth);

    plot(ax,temporalSupport*1e3, LMSdecrementImpulseResponses(:,cMosaic.MCONE_ID), 'k-', 'LineWidth', wideLineWidth);
    plot(ax,temporalSupport*1e3, LMSdecrementImpulseResponses(:,cMosaic.MCONE_ID), 'g:', 'LineWidth', narrowLineWidth);

    plot(ax,temporalSupport*1e3, LMSincrementImpulseResponses(:,cMosaic.SCONE_ID), 'k-', 'LineWidth', wideLineWidth);
    p3 = plot(ax,temporalSupport*1e3, LMSincrementImpulseResponses(:,cMosaic.SCONE_ID), 'c-', 'LineWidth', narrowLineWidth);

    plot(ax,temporalSupport*1e3, LMSdecrementImpulseResponses(:,cMosaic.SCONE_ID), 'k-', 'LineWidth', wideLineWidth);
    plot(ax,temporalSupport*1e3, LMSdecrementImpulseResponses(:,cMosaic.SCONE_ID), 'c:', 'LineWidth', narrowLineWidth);

    legend(ax, [p1 p2 p3], {'L-cone', 'M-cone', 'S-cone'});
    xlabel(ax, 'time (msec)')
    set(ax, 'YLim', amplitudeRangeOriginal, 'XLim', [temporalSupport(1) temporalSupport(end)]*1e3);
    set(ax, 'XTick', 0:100:2000, 'FontSize', 16);
    grid(ax, 'on');
    title(ax, sprintf('pCurrent impulse responses\n(dT = %1.3f msec)', (temporalSupport(2)-temporalSupport(1))*1e3));


    idx = find(temporalSupport>=0 & temporalSupport <=0.5);
    dtOriginal = (temporalSupport(2)-temporalSupport(1))*1e3;
    volOriginal = dtOriginal * trapz(abs(squeeze(LMSincrementImpulseResponses(idx,cMosaic.LCONE_ID))));

    idx = find(decimatedTemporalSupport>=0 & decimatedTemporalSupport <=0.5);
    dtDecimated = (decimatedTemporalSupport(2)-decimatedTemporalSupport(1))*1e3;
    volDecimated = dtDecimated * trapz(abs(squeeze(decimatedLMSincrementImpulseResponses(idx,cMosaic.LCONE_ID))));
    
    sprintf('\nVolume of original pCurrent IR: %f (dt: %2.1f msec), vol of decimated pCurrent IR: %f (dt: %2.1f msec)\n', ...
        volOriginal, dtOriginal, volDecimated, dtDecimated);


    wideLineWidth = 4;
    narrowLineWidth = 1.5;
    ax = subplot(1,2,2);
    
    plot(ax,decimatedTemporalSupport*1e3, decimatedLMSincrementImpulseResponses(:,cMosaic.LCONE_ID), ...
        'k-', 'LineWidth', wideLineWidth);
    hold(ax, 'on');
    p1 = plot(ax,decimatedTemporalSupport*1e3, decimatedLMSincrementImpulseResponses(:,cMosaic.LCONE_ID), ...
        'ro-', 'LineWidth', narrowLineWidth, 'MarkerSize', 16, 'MarkerFaceColor', [0.85 0 0.0], 'MarkerEdgeColor', 'k');
    
    plot(ax,decimatedTemporalSupport*1e3, decimatedLMSdecrementImpulseResponses(:,cMosaic.LCONE_ID), ...
        'k-', 'LineWidth', wideLineWidth);
    plot(ax,decimatedTemporalSupport*1e3, decimatedLMSdecrementImpulseResponses(:,cMosaic.LCONE_ID), ...
        'ro:', 'LineWidth', narrowLineWidth, 'MarkerSize', 16, 'MarkerFaceColor', [0.85 0 0.0], 'MarkerEdgeColor', 'k');

    plot(ax,decimatedTemporalSupport*1e3, decimatedLMSincrementImpulseResponses(:,cMosaic.MCONE_ID), ...
        'k-', 'LineWidth', wideLineWidth);
    p2 = plot(ax,decimatedTemporalSupport*1e3, decimatedLMSincrementImpulseResponses(:,cMosaic.MCONE_ID), ...
        'go-', 'LineWidth', narrowLineWidth, 'MarkerSize', 16, 'MarkerFaceColor', [0 0.85 0.0], 'MarkerEdgeColor', 'k');

    plot(ax,decimatedTemporalSupport*1e3, decimatedLMSdecrementImpulseResponses(:,cMosaic.MCONE_ID), ...
        'k-', 'LineWidth', wideLineWidth);
    plot(ax,decimatedTemporalSupport*1e3, decimatedLMSdecrementImpulseResponses(:,cMosaic.MCONE_ID), ...
        'go:', 'LineWidth', narrowLineWidth, 'MarkerSize', 16, 'MarkerFaceColor', [0 0.85 0.0], 'MarkerEdgeColor', 'k');

    plot(ax,decimatedTemporalSupport*1e3, decimatedLMSincrementImpulseResponses(:,cMosaic.SCONE_ID), ...
        'k-', 'LineWidth', wideLineWidth);
    p3 = plot(ax,decimatedTemporalSupport*1e3, decimatedLMSincrementImpulseResponses(:,cMosaic.SCONE_ID), ...
        'co-', 'LineWidth', narrowLineWidth, 'MarkerSize', 16, 'MarkerFaceColor', [0 0.85 0.85], 'MarkerEdgeColor', 'k');

    plot(ax,decimatedTemporalSupport*1e3, decimatedLMSdecrementImpulseResponses(:,cMosaic.SCONE_ID), ...
        'k-', 'LineWidth', wideLineWidth);
    plot(ax,decimatedTemporalSupport*1e3, decimatedLMSdecrementImpulseResponses(:,cMosaic.SCONE_ID), ...
        'co:', 'LineWidth', narrowLineWidth, 'MarkerSize', 16, 'MarkerFaceColor', [0 0.85 0.85], 'MarkerEdgeColor', 'k');
    xlabel(ax, 'time (msec)')
    legend(ax, [p1 p2 p3], {'L-cone', 'M-cone', 'S-cone'});
    set(ax, 'YLim', amplitudeRangeDecimated, 'XLim', [temporalSupport(1) temporalSupport(end)]*1e3);
    set(ax, 'XTick', 0:100:2000, 'FontSize', 16);
    grid(ax, 'on');
    title(ax, sprintf('pCurrent impulse responses\n(dT = %1.3f msec)', (decimatedTemporalSupport(2)-decimatedTemporalSupport(1))*1e3));
end



function [LMSincrementImpulseResponses, LMSdecrementImpulseResponses, temporalSupport, steadyStateCurrents] = ...
    computePhotocurrentResponses(backgroundLMSexcitationRates, flashIntensityLMSContrasts, flashDurationSeconds, radialEccDegs)

    os = osBioPhys('eccentricity', radialEccDegs);
    os.timeStep = 1e-4;
    os.set('noise flag', 'none');

    % Duration of photocurrent impulse response
    warmUpDurationSeconds = 0.5;
    impulseResponseDurationSeconds = 1.0 + warmUpDurationSeconds;
    nSamples = round(impulseResponseDurationSeconds / os.timeStep) + 1;

    warmupTimeSamples = round(warmUpDurationSeconds / os.timeStep);
    flashTimeSamples = round(flashDurationSeconds  / os.timeStep);

    coneTypesNum = length(backgroundLMSexcitationRates);
    LMSincrementImpulseResponses = zeros(nSamples, coneTypesNum);
    LMSdecrementImpulseResponses = LMSincrementImpulseResponses;
    steadyStateCurrents = zeros(1,coneTypesNum);

    for iConeClass = 1:length(backgroundLMSexcitationRates)
        % Get the isomerization count per time step for this cone class
        backgroundExcitationPhotons  = backgroundLMSexcitationRates(iConeClass) * os.timeStep;

        % The flash intensity contrast for this cone class
        flashIntensityContrast = flashIntensityLMSContrasts(iConeClass);

        % Create a constant stimulus at this rate
        backgroundStimulusPhotons = backgroundExcitationPhotons*ones(nSamples, 1);
        backgroundStimulusPhotonRate = backgroundStimulusPhotons/os.timeStep;

        theState = os.osAdaptSteadyState(backgroundExcitationPhotons, [1 1]);
        theState.timeStep = os.timeStep;
        os.setModelState(theState);
        backgroundCurrentTimeSeries = os.osAdaptTemporal(reshape(backgroundStimulusPhotonRate, [1 1 numel(backgroundStimulusPhotonRate)]));

        % The increment flash
        flashIntensityPhotons = flashIntensityContrast * backgroundExcitationPhotons;
        flashIncrementStimulusPhotons = backgroundStimulusPhotons;
        flashIncrementStimulusPhotons(warmupTimeSamples+(1:flashTimeSamples)) = ...
            flashIncrementStimulusPhotons(warmupTimeSamples+(1:flashTimeSamples)) + flashIntensityPhotons;
        flashIncrementStimulusPhotonRate = flashIncrementStimulusPhotons/os.timeStep;

        % The increment photocurrent response
        currentResponse = os.osAdaptTemporal(reshape(flashIncrementStimulusPhotonRate, [1 1 numel(flashIncrementStimulusPhotonRate)]));
        
        % The differential photocurrent (increment - steady state) divided
        % by the flash intensity
        deltaCurrent = squeeze(currentResponse - backgroundCurrentTimeSeries)/abs(flashIntensityPhotons);

        % Trim in time to get the impulse response
        LMSincrementImpulseResponses(:, iConeClass) = deltaCurrent;

        % The decrement flash
        flashIntensityPhotons = -flashIntensityPhotons;
        flashDecrementStimulusPhotons = backgroundStimulusPhotons;
        flashDecrementStimulusPhotons(warmupTimeSamples+(1:flashTimeSamples)) = ...
            flashDecrementStimulusPhotons(warmupTimeSamples+(1:flashTimeSamples)) + flashIntensityPhotons;
        flashDecrementStimulusPhotonRate = flashDecrementStimulusPhotons/os.timeStep;

        % The decrement photocurrent response
        currentResponse = os.osAdaptTemporal(reshape(flashDecrementStimulusPhotonRate, [1 1 numel(flashDecrementStimulusPhotonRate)]));
        
        % The differential photocurrent (decrement - steady state) divided
        % by the flash intensity
        deltaCurrent = squeeze(currentResponse - backgroundCurrentTimeSeries)/abs(flashIntensityPhotons);

        % Trim in time to get the impulse response
        LMSdecrementImpulseResponses(:, iConeClass) = deltaCurrent;

        steadyStateCurrents(iConeClass) = backgroundCurrentTimeSeries(1,1,end);
    end % iConeClass

    temporalSupport = (1:size(LMSdecrementImpulseResponses,1)) * os.timeStep - warmUpDurationSeconds;

    % Only keep positive time data
    idx = find(temporalSupport>=0);
    temporalSupport = temporalSupport(idx);
    LMSdecrementImpulseResponses = LMSdecrementImpulseResponses(idx,:);
    LMSincrementImpulseResponses = LMSincrementImpulseResponses(idx,:);
end
