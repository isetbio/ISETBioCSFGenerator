% Visualization function for CMosaic based neural responses 
%
% Add comments
%

% History:
%    01/11/2025  NPC  Wrote it

function nreVisualizeMRGCmosaic(theMRGCmosaic, theNeuralResponses, theConeMosaicResponses, temporalSupportSeconds, responseLabel)
    
    % Back to mRGCMosaic's native response format because we will be calling
    % cMosaic visualizer functions which expect it
    if (size(theNeuralResponses,2) == theMRGCmosaic.rgcsNum)
        theNeuralResponses = permute(theNeuralResponses,[1 3 2]);
    end

    hFig = figure(1999);
    clf;
    set(hFig, 'Position', [10 10 1700 500], 'Color', [1 1 1]);
    axInputConeMosaicActivation = subplot('Position', [0.05 0.1 0.3 0.85]);
    axMRGCMosaicActivation = subplot('Position', [0.37 0.1 0.3 0.85]);
    axResponseTimeResponses = subplot('Position', [0.7 0.1 0.3 0.85]);

    [nTrials, nTimePoints, nCones] = size(theNeuralResponses);
    activationRangeMRGCMosaic = prctile(theNeuralResponses(:),[1 99]);
    activationRangeConeMosaic = prctile(theConeMosaicResponses(:),[1 99]);

    if (~isempty(theMRGCmosaic.inputConeMosaic.fixEMobj)) 
        emPathsDegs = theMRGCmosaic.inputConeMosaic.fixEMobj.emPosArcMin/60;
    else
        emPathsDegs = [];
    end

    for iTrial = 1:nTrials
        for iPoint = 1:nTimePoints

            theConeMosaicReponseTrial = min([iTrial size(theConeMosaicResponses,1)]);

            if (isempty(emPathsDegs))
                theMRGCmosaic.inputConeMosaic.visualize(...
                    'figureHandle', hFig, ...
                    'axesHandle', axInputConeMosaicActivation, ...
                    'activation', theConeMosaicResponses(theConeMosaicReponseTrial, iPoint,:), ...
                    'activationRange', activationRangeConeMosaic , ...
                    'plotTitle', sprintf('cMosaic %s (t: %2.1f msec)', responseLabel, temporalSupportSeconds(iPoint)*1e3));
            else
                theEMtrial = min([iTrial size(emPathsDegs,1)]);

                theMRGCmosaic.inputConeMosaic.visualize(...
                    'figureHandle', hFig, ...
                    'axesHandle', axInputConeMosaicActivation, ...
                    'activation', theConeMosaicResponses(theConeMosaicReponseTrial, iPoint,:), ...
                    'activationRange', activationRangeConeMosaic , ...
                    'currentEMposition', squeeze(emPathsDegs(theEMtrial,iPoint,:)), ...
                    'displayedEyeMovementData', struct('trial', theEMtrial, 'timePoints', 1:iPoint), ...
                    'plotTitle', sprintf('cMosaic %s (t: %2.1f msec, trial %d of %d)', responseLabel, temporalSupportSeconds(iPoint)*1e3, iTrial, nTrials));
     
            end
            
            theMRGCmosaic.visualize(...
                'figureHandle', hFig, ...
                'axesHandle', axMRGCMosaicActivation, ...
                'activation', theNeuralResponses(iTrial, iPoint,:), ...
                'activationRange', activationRangeMRGCMosaic, ...
                'plotTitle', sprintf('mRGCMosaic %s (t: %2.1f msec)', responseLabel, temporalSupportSeconds(iPoint)*1e3));
            drawnow;
        end % iPoint

        
        mosaicSpatioTemporalActivation(:, 1:theMRGCmosaic.rgcsNum) = squeeze(theNeuralResponses(iTrial,:,:));
        dt = 0.5*(temporalSupportSeconds(2)-temporalSupportSeconds(1));
        imagesc(axResponseTimeResponses, temporalSupportSeconds*1e3, 1:theMRGCmosaic.rgcsNum, mosaicSpatioTemporalActivation');
        hold(axResponseTimeResponses, 'on');
        for i = 1:numel(temporalSupportSeconds)
            plot(axResponseTimeResponses, (temporalSupportSeconds(i)-dt)*[1 1]*1e3, [1 theMRGCmosaic.rgcsNum], 'k-', 'LineWidth', 1.0);
        end
        hold(axResponseTimeResponses, 'off');

        linearRamp = linspace(0,1,1024);
        linearRamp = linearRamp(:);
        cMap = [linearRamp*0 linearRamp linearRamp*0];

        colormap(axResponseTimeResponses, cMap);

        axis(axResponseTimeResponses, 'xy');
        xlabel(axResponseTimeResponses, 'time (msec)');
        ylabel(axResponseTimeResponses, 'mRGC index');
        set(axResponseTimeResponses, 'FontSize', 16, 'Color', [1 1 1], 'CLim', activationRangeMRGCMosaic);
        title(axResponseTimeResponses, sprintf('spatio-temporal %s (trial %d of %d)', responseLabel, iTrial, nTrials));
        drawnow;

    end % iTrial

end