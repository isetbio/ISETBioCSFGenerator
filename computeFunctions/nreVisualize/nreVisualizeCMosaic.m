% Visualization function for CMosaic based neural responses 
%
% Add comments
%

% History:
%    01/11/2025  NPC  Wrote it

function nreVisualizeCMosaic(theConeMosaic, theNeuralResponses, temporalSupportSeconds, responseLabel)

    % Back to cMosaic's native response format because we will be calling
    % cMosaic visualizer functions which expect it
    if (size(theNeuralResponses,2) == theConeMosaic.conesNum)
        theNeuralResponses = permute(theNeuralResponses,[1 3 2]);
    end
    

    hFig = figure(1999);
    clf;
    set(hFig, 'Position', [10 10 1700 500], 'Color', [1 1 1]);
    axCmosaicActivation = subplot(1,2,1);
    axConeResponseTimeResponses = subplot(1,2,2);

    [nTrials, nTimePoints, nCones] = size(theNeuralResponses);
    activationRange = prctile(theNeuralResponses(:),[5 95]);

    theLconeResponses = theNeuralResponses(:,:,theConeMosaic.lConeIndices);
    theMconeResponses = theNeuralResponses(:,:,theConeMosaic.mConeIndices);
    theSconeResponses = theNeuralResponses(:,:,theConeMosaic.sConeIndices);


    if (~isempty(theConeMosaic.fixEMobj)) 
        emPathsDegs = theConeMosaic.fixEMobj.emPosArcMin/60;
    else
        emPathsDegs = [];
    end

    for iTrial = 1:nTrials
        for iPoint = 1:nTimePoints
            if (isempty(emPathsDegs))
                theConeMosaic.visualize(...
                    'figureHandle', hFig, ...
                    'axesHandle', axCmosaicActivation, ...
                    'activation', theNeuralResponses(iTrial, iPoint,:), ...
                    'activationRange', activationRange, ...
                    'plotTitle', sprintf('%s (t: %2.1f msec)', responseLabel, temporalSupportSeconds(iPoint)*1e3));
                
            else
                theEMtrial = min([iTrial size(emPathsDegs,1)]);

                theConeMosaic.visualize(...
                    'figureHandle', hFig, ...
                    'axesHandle', axCmosaicActivation, ...
                    'activation', theNeuralResponses(iTrial, iPoint,:), ...
                    'activationRange', activationRange, ...
                    'currentEMposition', squeeze(emPathsDegs(theEMtrial,iPoint,:)), ...
                    'displayedEyeMovementData', struct('trial', theEMtrial, 'timePoints', 1:iPoint), ...
                    'plotTitle', sprintf('%s (t: %2.1f msec, trial %d of %d)', responseLabel, temporalSupportSeconds(iPoint)*1e3, iTrial, nTrials));
     
            end 
            drawnow;
        end % iPoint

        mosaicSpatioTemporalActivation = [];
        mosaicSpatioTemporalActivation(:, 1:numel(theConeMosaic.lConeIndices)) = squeeze(theLconeResponses(iTrial,:,:));
        mosaicSpatioTemporalActivation(:, size(mosaicSpatioTemporalActivation,2)+(1:numel(theConeMosaic.mConeIndices))) = squeeze(theMconeResponses(iTrial,:,:));
        mosaicSpatioTemporalActivation(:, size(mosaicSpatioTemporalActivation,2)+(1:numel(theConeMosaic.sConeIndices))) = squeeze(theSconeResponses(iTrial,:,:));
        
        dt = 0.5*(temporalSupportSeconds(2)-temporalSupportSeconds(1));
        imagesc(axConeResponseTimeResponses, temporalSupportSeconds*1e3, 1:theConeMosaic.conesNum, mosaicSpatioTemporalActivation');
        hold(axConeResponseTimeResponses, 'on');
        LconeRect.x = [temporalSupportSeconds(1)-dt temporalSupportSeconds(end)+dt temporalSupportSeconds(end)+dt temporalSupportSeconds(1)-dt temporalSupportSeconds(1)-dt]*1e3;
        LconeRect.y = [1  1 numel(theConeMosaic.lConeIndices) numel(theConeMosaic.lConeIndices) 1];
        MconeRect.x = [temporalSupportSeconds(1)-dt temporalSupportSeconds(end)+dt temporalSupportSeconds(end)+dt temporalSupportSeconds(1)-dt temporalSupportSeconds(1)-dt]*1e3;
        MconeRect.y = numel(theConeMosaic.lConeIndices) + [1  1 numel(theConeMosaic.mConeIndices) numel(theConeMosaic.mConeIndices) 1];
        SconeRect.x = [temporalSupportSeconds(1)-dt temporalSupportSeconds(end)+dt temporalSupportSeconds(end)+dt temporalSupportSeconds(1)-dt temporalSupportSeconds(1)-dt]*1e3;
        SconeRect.y = numel(theConeMosaic.lConeIndices) + numel(theConeMosaic.mConeIndices) + [1  1 numel(theConeMosaic.sConeIndices) numel(theConeMosaic.sConeIndices) 1];
        
        plot(axConeResponseTimeResponses, LconeRect.x, LconeRect.y, 'r-', 'LineWidth', 2);
        plot(axConeResponseTimeResponses, MconeRect.x, MconeRect.y, 'g-', 'LineWidth', 2);
        plot(axConeResponseTimeResponses, SconeRect.x, SconeRect.y, 'c-', 'LineWidth', 2);
        for i = 1:numel(temporalSupportSeconds)
            plot(axConeResponseTimeResponses, (temporalSupportSeconds(i)-dt)*[1 1]*1e3, [1 theConeMosaic.conesNum], 'k-', 'LineWidth', 1.0);
        end

        colormap(brewermap(1024, '*greys'));
    
        hold(axConeResponseTimeResponses, 'off');
        set(axConeResponseTimeResponses, ...
            'XLim', [temporalSupportSeconds(1)-dt temporalSupportSeconds(end)+dt]*1e3, ...
            'YLim', [1 theConeMosaic.conesNum], ...
            'CLim', [activationRange(1) activationRange(2)]);
        axis(axConeResponseTimeResponses, 'xy');
        xlabel(axConeResponseTimeResponses, 'time (msec)');
        ylabel(axConeResponseTimeResponses, sprintf('cone index (%d L-cones, %d M-cones, %d S-cones', numel(theConeMosaic.lConeIndices), numel(theConeMosaic.mConeIndices), numel(theConeMosaic.sConeIndices)));
        set(axConeResponseTimeResponses, 'FontSize', 16, 'Color', [1 1 1], 'CLim', activationRange);
        title(axConeResponseTimeResponses, sprintf('spatio-temporal %s (trial %d of %d)', responseLabel, iTrial, nTrials));
        drawnow;

    end % iTrial
end