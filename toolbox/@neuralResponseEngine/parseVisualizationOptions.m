function [figureHandle, axesHandle, clearAxesBeforeDrawing, responseLabel, maxVisualizedNoisyResponseInstances, visualizationMetaData] = ...
    parseVisualizationOptions(options)
    
    arguments
        options.figureHandle = [];
        options.axesHandle = [];
        options.responseLabel (1,:) char  = ''
        options.clearAxesBeforeDrawing (1,1) logical = true;
        options.maxVisualizedNoisyResponseInstances (1,1) double = inf;
        options.visualizationMetaData = [];
    end

    figureHandle = options.figureHandle;
    axesHandle = options.axesHandle;
    responseLabel = options.responseLabel;
    clearAxesBeforeDrawing = options.clearAxesBeforeDrawing;
    maxVisualizedNoisyResponseInstances = options.maxVisualizedNoisyResponseInstances;
    visualizationMetaData = options.visualizationMetaData;

    assert(isEmptyOrHandle('figure', figureHandle), 'Figure handle is not valid');
    
    if (numel(axesHandle) <= 1)
        assert(isEmptyOrHandle('axes', axesHandle), 'Axes handle is not valid');
    elseif iscell(axesHandle)
         for iAxisHandle = 1:numel(axesHandle)
             theAxisHandle = axesHandle{iAxisHandle};
             assert(mustBeEmptyOrHandle('axes', theAxisHandle), sprintf('Axes handle %d is not valid', iAxisHandle));
         end
    else
      error('invalid passed axesHandle: nor an axis handle, neither a cell array of axis handles.')
    end

end


% Validation function
function val = isEmptyOrHandle(argumentName, x)
    if (isempty(x) || ishandle(x))
        val = true;
    else
        val = false;
        sprintf('must be EITHER empty OR a %s handle\n', argumentName);
    end
end