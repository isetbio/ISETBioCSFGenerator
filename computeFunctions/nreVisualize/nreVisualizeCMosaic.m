% Visualization function for CMosaic based neural responses 
%
% Add comments
%

% History:
%    01/11/2025  NPC  Wrote it

function nreVisualizeCMosaic(neuralPipeline, neuralResponses, temporalSupportSeconds, varargin)

    [figureHandle, axesHandle, clearAxesBeforeDrawing, responseLabel, maxVisualizedInstances, visualizationMetaData] = ...
        neuralResponseEngine.parseVisualizationOptions(varargin{:});

    if (isempty(figureHandle))
        figureHandle = figure(); clf;
        set(figureHandle, 'Position', [10 10 1700 500]);
    end
    if (isempty(axesHandle))
        axCmosaicActivation = subplot(1,2,1);
        axConeTimeResponses = subplot(1,2,2);
    else
        if (iscell(axesHandle)&&(numel(axesHandle) == 2))
            axCmosaicActivation = axesHandle{1};
            axConeTimeResponses = axesHandle{2};
        else
            warning('nreVisualizeCMosaic:expected a cell array with 2 axes handles but did not receive it', 'Will generate new axes.');
            axCmosaicActivation = subplot(1,2,1);
            axConeTimeResponses = subplot(1,2,2);
        end
    end


    % Retrieve the visualization data we need from the pipeline
    theConeMosaic = neuralPipeline.noiseFreeResponse.coneMosaic;

    % Back to cMosaic's native response format because we will be calling
    % cMosaic visualizer functions which expect it
    if (size(neuralResponses,2) == theConeMosaic.conesNum)
        neuralResponses = permute(neuralResponses,[1 3 2]);
    end
    

    % Get the response dimensins and range
    [nInstances, nTimePoints, ~] = size(neuralResponses);
    activationRange = [0 max(neuralResponses(:))];


    if (~isempty(theConeMosaic.fixEMobj)) 
        emPathsDegs = theConeMosaic.fixEMobj.emPosArcMin/60;
    else
        emPathsDegs = [];
    end

    if (numel(temporalSupportSeconds)>1)
        dt = temporalSupportSeconds(2)-temporalSupportSeconds(1);
        XLim = [temporalSupportSeconds(1)-dt/2 temporalSupportSeconds(end)+dt/2];
        XTick = temporalSupportSeconds(1):dt:temporalSupportSeconds(end);
    else
        dt = 0;
        XLim = [temporalSupportSeconds(1)-0.01 temporalSupportSeconds(1)+0.01];
        XTick = temporalSupportSeconds(1);
    end

    for iTrial = 1:min([nInstances maxVisualizedInstances])

        % The instantaneous spatial activation
        for iPoint = 1:nTimePoints

            % Retrieve spatiotemporal response up to this time point
            [mosaicSpatioTemporalActivation, LconeRect, MconeRect, SconeRect] = ...
                spatioTemporalResponseComponents(theConeMosaic, neuralResponses, temporalSupportSeconds, iTrial, iPoint);
            
            % The spatiotemporal mosaic activation up to this time point
            imagesc(axConeTimeResponses, temporalSupportSeconds, 1:theConeMosaic.conesNum, ...
                mosaicSpatioTemporalActivation);
            hold(axConeTimeResponses, 'on');

            % The L,M,S cone rectangles
            plot(axConeTimeResponses, LconeRect.x, LconeRect.y, 'r-', 'LineWidth', 2);
            plot(axConeTimeResponses, MconeRect.x, MconeRect.y, 'g-', 'LineWidth', 2);
            plot(axConeTimeResponses, SconeRect.x, SconeRect.y, 'c-', 'LineWidth', 2);

            % The stimulus frames
            for i = 1:numel(temporalSupportSeconds)
                plot(axConeTimeResponses, (temporalSupportSeconds(i)-dt/2)*[1 1], [1 theConeMosaic.conesNum], 'k-', 'LineWidth', 1.0);
            end

            hold(axConeTimeResponses, 'off');
            axis(axConeTimeResponses, 'xy');

            colormap(axConeTimeResponses, brewermap(1024, '*greys'));

            set(axConeTimeResponses, ...
                'XLim', XLim, ...
                'XTick', XTick, ...
                'YLim', [1 theConeMosaic.conesNum], ...
                'CLim', activationRange, ...
                'Color', [0 0 0], ...
                'FontSize', 16);

            colorbar(axConeTimeResponses, 'NorthOutside');

            xlabel(axConeTimeResponses, 'time (sec)');
            ylabel(axConeTimeResponses, sprintf('cone index (%d L-cones, %d M-cones, %d S-cones', numel(theConeMosaic.lConeIndices), numel(theConeMosaic.mConeIndices), numel(theConeMosaic.sConeIndices)));
            title(axConeTimeResponses, sprintf('spatio-temporal %s (trial %d of %d)', responseLabel, iTrial, nInstances));
            


            if (isempty(emPathsDegs))
                theConeMosaic.visualize(...
                    'figureHandle', figureHandle, ...
                    'axesHandle', axCmosaicActivation, ...
                    'activation', neuralResponses(iTrial, iPoint,:), ...
                    'activationRange', activationRange, ...
                    'clearAxesBeforeDrawing', clearAxesBeforeDrawing, ...
                    'plotTitle', sprintf('%s (t: %2.1f msec)', responseLabel, temporalSupportSeconds(iPoint)*1e3));
                
            else
                theEMtrial = min([iTrial size(emPathsDegs,1)]);

                theConeMosaic.visualize(...
                    'figureHandle', figureHandle, ...
                    'axesHandle', axCmosaicActivation, ...
                    'activation', neuralResponses(iTrial, iPoint,:), ...
                    'activationRange', activationRange, ...
                    'currentEMposition', squeeze(emPathsDegs(theEMtrial,iPoint,:)), ...
                    'displayedEyeMovementData', struct('trial', theEMtrial, 'timePoints', 1:iPoint), ...
                    'clearAxesBeforeDrawing', clearAxesBeforeDrawing, ...
                    'plotTitle', sprintf('%s (t: %2.1f msec, trial %d of %d)', responseLabel, temporalSupportSeconds(iPoint)*1e3, iTrial, nInstances));
     
            end 

            drawnow;

        end % iPoint
    end % iTrial
end



function [mosaicSpatioTemporalActivation, LconeRect, MconeRect, SconeRect] = spatioTemporalResponseComponents(theConeMosaic, neuralResponses, temporalSupportSeconds, iTrial, iPoint)

    % Retrieve the spatiotemporal L-, M- and S-cone responses
    % so we can plot them separately along the y-axis
    theLconeResponses = neuralResponses(iTrial,1:iPoint,theConeMosaic.lConeIndices);
    theMconeResponses = neuralResponses(iTrial,1:iPoint,theConeMosaic.mConeIndices);
    theSconeResponses = neuralResponses(iTrial,1:iPoint,theConeMosaic.sConeIndices);

    % The spatiotemporal activation
    lConesNum = numel(theConeMosaic.lConeIndices);
    mConesNum = numel(theConeMosaic.mConeIndices);
    sConesNum = numel(theConeMosaic.sConeIndices);
    mosaicSpatioTemporalActivation = zeros(theConeMosaic.conesNum, numel(temporalSupportSeconds));
    mosaicSpatioTemporalActivation(1:lConesNum, 1:iPoint) = (squeeze(theLconeResponses))';
    mosaicSpatioTemporalActivation((lConesNum)+(1:mConesNum), 1:iPoint) = (squeeze(theMconeResponses))';
    mosaicSpatioTemporalActivation((lConesNum+mConesNum)+(1:sConesNum), 1:iPoint) = (squeeze(theSconeResponses))';

    % The L-, M-, and S-cone rects
    if (numel(temporalSupportSeconds) == 1)
        dt = 1/1000;
    else
        dt = 0.5*(temporalSupportSeconds(2)-temporalSupportSeconds(1));
    end
    LconeRect.x = [temporalSupportSeconds(1)-dt temporalSupportSeconds(end)+dt temporalSupportSeconds(end)+dt temporalSupportSeconds(1)-dt temporalSupportSeconds(1)-dt]*1e3;
    LconeRect.y = [1  1 numel(theConeMosaic.lConeIndices) numel(theConeMosaic.lConeIndices) 1];
    MconeRect.x = LconeRect.x;
    MconeRect.y = numel(theConeMosaic.lConeIndices) + [1  1 numel(theConeMosaic.mConeIndices) numel(theConeMosaic.mConeIndices) 1];
    SconeRect.x = LconeRect.x;
    SconeRect.y = numel(theConeMosaic.lConeIndices) + numel(theConeMosaic.mConeIndices) + [1  1 numel(theConeMosaic.sConeIndices) numel(theConeMosaic.sConeIndices) 1];
end