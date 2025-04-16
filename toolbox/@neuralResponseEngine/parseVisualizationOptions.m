function [figureHandle, axesHandle, clearAxesBeforeDrawing, responseLabel, ...
    responseVideoFileName, neuralPipelineID, visualizeResponsesAsModulations] = ...
    parseVisualizationOptions(options)
    
    arguments
        options.figureHandle = [];
        options.axesHandle = [];
        options.responseLabel (1,:) char  = ''
        options.clearAxesBeforeDrawing (1,1) logical = true;
        options.responseVideoFileName (1,:) char  = ''
        options.neuralPipelineID (1,:) char  = ''
        options.visualizeResponsesAsModulations (1,1) logical = false;
    end

    figureHandle = options.figureHandle;
    axesHandle = options.axesHandle;
    responseLabel = options.responseLabel;
    clearAxesBeforeDrawing = options.clearAxesBeforeDrawing;
    responseVideoFileName = options.responseVideoFileName;
    neuralPipelineID = options.neuralPipelineID;
    visualizeResponsesAsModulations = options.visualizeResponsesAsModulations;
    
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