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
        neuralResponseEngine.parseDefaultVisualizationOptionsStruct(varargin{:});

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
    assert(tBins == size(temporalSupportSeconds,2), 'Inconsistency in time dimension');

    visualizedInstancesNum = min([nInstances maxVisualizedInstances]);

    if (mNeurons == 1)
        plot(axesHandle, temporalSupportSeconds*1e3, neuralResponses(1:visualizedInstancesNum,1,:), 'k-', 'LineWidth', 1.5);
        xlabel(axesHandle, 'time (msec)');
        ylabel(axesHandle, 'response amplitude');
        title(axesHandle, responseLabel);
        set(axesHandle, 'fontSize', 16);
        drawnow;
    elseif (mNeurons > 1)
        set(axesHandle, 'fontSize', 16);
        xlabel(axesHandle, 'time (msec)');
        ylabel(axesHandle, 'neuron index');
        for iTrial = 1:visualizedInstancesNum
            imagesc(axesHandle, temporalSupportSeconds*1e3, 1:mNeurons, squeeze(neuralResponses(iTrial,:,:)));
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