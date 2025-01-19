function [figureHandle, axesHandle, clearAxesBeforeDrawing, responseLabel, maxVisualizedNoisyResponseInstances] = ...
    parseVisualizationOptionsStruct(optionsStruct)
    
    arguments
        optionsStruct.figureHandle = [];
        optionsStruct.axesHandle = [];
        optionsStruct.responseLabel (1,:) char  = ''
        optionsStruct.clearAxesBeforeDrawing (1,1) logical = true;
        optionsStruct.maxVisualizedNoisyResponseInstances (1,1) double = inf;
    end

    figureHandle = optionsStruct.figureHandle;
    axesHandle = optionsStruct.axesHandle;
    responseLabel = optionsStruct.responseLabel;
    clearAxesBeforeDrawing = optionsStruct.clearAxesBeforeDrawing;
    maxVisualizedNoisyResponseInstances = optionsStruct.maxVisualizedNoisyResponseInstances;

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



 % if (isstruct(optionsStruct))
 %        % Default values
 %        figureHandle = [];
 %        axesHandle = [];
 %        responseLabel = '';
 %        maxVisualizedInstances = inf;
 %        clearAxesBeforeDrawing = true;
 % 
 %        % Parse options struct
 %        fNames = fieldnames(optionsStruct);
 % 
 %        for iField = 1:numel(fNames)
 %            switch(fNames{iField})
 %                case 'figureHandle'
 %                    assert(isEmptyOrHandle('figure', optionsStruct.figureHandle), 'Figure handle is not valid');
 %                    figureHandle = optionsStruct.figureHandle;
 % 
 %                case 'axesHandle'
 %                    if (numel(optionsStruct.axesHandle) == 1)
 %                        assert(isEmptyOrHandle('axes', optionsStruct.axesHandle), 'Axes handle is not valid');
 %                    elseif iscell(optionsStruct.axesHandle)
 %                        for iAxisHandle = 1:numel(optionsStruct.axesHandle)
 %                            theAxisHandle = optionsStruct.axesHandle{iAxisHandle};
 %                            assert(mustBeEmptyOrHandle('axes', theAxisHandle), sprintf('Axes handle %d is not valid', iAxisHandle));
 %                        end
 %                    else
 %                        error('invalid passes axesHandle: nor an axis handle, neither a cell array of axis handles.')
 %                    end
 %                    axesHandle = optionsStruct.axesHandle;
 % 
 %                case 'responseLabel'
 %                    assert(ischar(optionsStruct.responseLabel), 'Response label is not a char');
 %                    responseLabel = optionsStruct.responseLabel;
 % 
 %                case 'maxVisualizedNoisyResponseInstances'
 %                    assert(isnumeric(optionsStruct.maxVisualizedNoisyResponseInstances), 'maxVisualizedNoisyResponseInstances is not  numeric');
 %                    maxVisualizedInstances = optionsStruct.maxVisualizedNoisyResponseInstances;
 % 
 %                case 'clearAxesBeforeDrawing'
 %                    assert(islogical(optionsStruct.clearAxesBeforeDrawing), 'clearAxesBeforeDrawing is not boolean');
 %                    clearAxesBeforeDrawing = optionsStruct.clearAxesBeforeDrawing;
 % 
 %                otherwise
 %                    warning('Uknown visualization option: ''%s'' is ignored,', fNames{iField});
 % 
 %            end % switch
 %        end % for iField
 %    else
 %        warning('Passed visualization options argument is not a struct. Ignoring it.')
 %    end % if (isstruct(optionsStruct))