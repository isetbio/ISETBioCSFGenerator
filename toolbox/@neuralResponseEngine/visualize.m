function visualize(obj, neuralResponses, temporalSupportSeconds, varargin)
% Gateway for visualization of @neuralResponseEngine computed data.

    if (~isempty(obj.customVisualizationFunctionHandle))
        if (~isempty(obj.neuralPipeline))
            % User specified custom visualization function and we have a non-empty neural pipeline
            obj.customVisualizationFunctionHandle(obj.neuralPipeline, neuralResponses, temporalSupportSeconds, varargin{:});
        end
    else
        % Just the default visualization
        defaultVisualization(neuralResponses, temporalSupportSeconds, varargin{:});
    end
end

function defaultVisualization(neuralResponses, temporalSupportSeconds, varargin)
   
    [figureHandle, axesHandle, clearAxesBeforeDrawing, responseLabel, maxVisualizedInstances] = ...
        neuralResponseEngine.parseVisualizationOptions(varargin{:});

    %% Set figure size
    if (isempty(figureHandle))
        figureHandle = figure(); clf;
        set(figureHandle, 'Position', [10 10 700 700], 'Color', [1 1 1]);
    end
    if (isempty(axesHandle))
        axesHandle = subplot('Position', [0.09 0.07 0.85 0.90]);
    end

    % If we were passed more than one axes handles, just use the first one
    if (iscell(axesHandle))&&(numel(axesHandle) > 1)
            axesHandle = axesHandle{1};
    end

    if (clearAxesBeforeDrawing)
        cla(axesHandle);
    end


    [nInstances, mNeurons, tBins] = size(neuralResponses);
    activationRange = [min(neuralResponses(:)) max(neuralResponses(:))];
    assert(tBins == size(temporalSupportSeconds,2), 'Inconsistency in time dimension');

    visualizedInstancesNum = min([nInstances maxVisualizedInstances]);

    if (numel(temporalSupportSeconds)>1)
        dt = temporalSupportSeconds(2)-temporalSupportSeconds(1);
        XLim = [temporalSupportSeconds(1)-dt/2 temporalSupportSeconds(end)+dt/2];
        XTick = temporalSupportSeconds(1):dt:temporalSupportSeconds(end);
    else
        XLim = [temporalSupportSeconds(1)-0.01 temporalSupportSeconds(1)+0.01];
        XTick = temporalSupportSeconds(1);
    end

    if (mNeurons == 1)
        plot(axesHandle, temporalSupportSeconds, neuralResponses(1:visualizedInstancesNum,1,:), 'k-', 'LineWidth', 1.5);
        xlabel(axesHandle, 'time (msec)');
        ylabel(axesHandle, 'response amplitude');
        title(axesHandle, responseLabel);
        set(axesHandle, 'XLim', XLim, 'YLim', activationRange, 'XTick', XTick);
        set(axesHandle, 'fontSize', 16);
        drawnow;
    elseif (mNeurons > 1)
        set(axesHandle, 'fontSize', 16);
        xlabel(axesHandle, 'time (msec)');
        ylabel(axesHandle, 'neuron index');
        
        for iTrial = 1:visualizedInstancesNum
            theNeuralResponseMap = squeeze(neuralResponses(iTrial,:,:));
            if (size(theNeuralResponseMap,2) == mNeurons)
                 theNeuralResponseMap =  theNeuralResponseMap';
            end

            imagesc(axesHandle, temporalSupportSeconds, 1:mNeurons, theNeuralResponseMap);
            set(axesHandle, ...
                'CLim', activationRange, ...
                'XLim', XLim, ...
                'XTick', XTick, ...
                'YLim', [1 mNeurons]);
            set(axesHandle, 'fontSize', 16);
            title(axesHandle, sprintf('%s (trial: %d)',responseLabel, iTrial));
            colormap(gray);
            drawnow;
        end
    end
    
end


% Validation function
function mustBeEmptyOrHandle(argumentName, x)
    if (isempty(x) || ishandle(x))
    else
        error(neuralResponseEngine:visualize,'must be EITHER empty OR a %s handle', argumentName);
    end
end