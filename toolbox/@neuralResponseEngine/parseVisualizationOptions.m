function [figureHandle, axesHandle, clearAxesBeforeDrawing, responseLabel, maxVisualizedNoisyResponseInstances] = ...
    parseVisualizationOptions(options)
    
    arguments
        options.figureHandle = [];
        options.axesHandle = [];
        options.responseLabel (1,:) char  = ''
        options.clearAxesBeforeDrawing (1,1) logical = true;
        options.maxVisualizedNoisyResponseInstances (1,1) double = inf;
    end

    figureHandle = options.figureHandle;
    axesHandle = options.axesHandle;
    responseLabel = options.responseLabel;
    clearAxesBeforeDrawing = options.clearAxesBeforeDrawing;
    maxVisualizedNoisyResponseInstances = options.maxVisualizedNoisyResponseInstances;

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



 % if (isstruct(options))
 %        % Default values
 %        figureHandle = [];
 %        axesHandle = [];
 %        responseLabel = '';
 %        maxVisualizedInstances = inf;
 %        clearAxesBeforeDrawing = true;
 % 
 %        % Parse options struct
 %        fNames = fieldnames(options);
 % 
 %        for iField = 1:numel(fNames)
 %            switch(fNames{iField})
 %                case 'figureHandle'
 %                    assert(isEmptyOrHandle('figure', options.figureHandle), 'Figure handle is not valid');
 %                    figureHandle = options.figureHandle;
 % 
 %                case 'axesHandle'
 %                    if (numel(options.axesHandle) == 1)
 %                        assert(isEmptyOrHandle('axes', options.axesHandle), 'Axes handle is not valid');
 %                    elseif iscell(options.axesHandle)
 %                        for iAxisHandle = 1:numel(options.axesHandle)
 %                            theAxisHandle = options.axesHandle{iAxisHandle};
 %                            assert(mustBeEmptyOrHandle('axes', theAxisHandle), sprintf('Axes handle %d is not valid', iAxisHandle));
 %                        end
 %                    else
 %                        error('invalid passes axesHandle: nor an axis handle, neither a cell array of axis handles.')
 %                    end
 %                    axesHandle = options.axesHandle;
 % 
 %                case 'responseLabel'
 %                    assert(ischar(options.responseLabel), 'Response label is not a char');
 %                    responseLabel = options.responseLabel;
 % 
 %                case 'maxVisualizedNoisyResponseInstances'
 %                    assert(isnumeric(options.maxVisualizedNoisyResponseInstances), 'maxVisualizedNoisyResponseInstances is not  numeric');
 %                    maxVisualizedInstances = options.maxVisualizedNoisyResponseInstances;
 % 
 %                case 'clearAxesBeforeDrawing'
 %                    assert(islogical(options.clearAxesBeforeDrawing), 'clearAxesBeforeDrawing is not boolean');
 %                    clearAxesBeforeDrawing = options.clearAxesBeforeDrawing;
 % 
 %                otherwise
 %                    warning('Uknown visualization option: ''%s'' is ignored,', fNames{iField});
 % 
 %            end % switch
 %        end % for iField
 %    else
 %        warning('Passed visualization options argument is not a struct. Ignoring it.')
 %    end % if (isstruct(options))