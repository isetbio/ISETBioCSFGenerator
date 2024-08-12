function validateAndSetParamsStruct(obj, paramsStruct)
    
    % Accept an empty params struct
    if (isempty(paramsStruct))
        % Set the neuralParams to the default params struct, which is returned
        % by the neuralComputeFunction when it is called with no input arguments
        obj.neuralParams = obj.neuralComputeFunction();
        return;
    end
    
    % Assert that the input is a struct
    assert(isstruct(paramsStruct), ...
        'Expected a struct during neuralResponseEngine instantiation');
    
    % DHB: This type checking too restrictive as our view of what is
    % possible with this general compute structure expands.  Commented out.
    %
%    % Obligatory sub-truct: 'opticsParams'
%    % Assert that we have a sub-struct named 'opticsParams'
%    assert((isfield(paramsStruct, 'opticsParams'))&&(isstruct(paramsStruct.opticsParams)), ...
%        'Expected an ''opticsParams'' substruct during neuralResponseEngine instantiation.');
%    
%     % Obligatory sub-struct: 'coneMosaicParams'
%     % Assert that we have a sub-struct named 'coneMosaicParams'
%     assert((isfield(paramsStruct, 'coneMosaicParams'))&&(isstruct(paramsStruct.coneMosaicParams)), ...
%         'Expected a ''coneMosaicParams'' substruct during neuralResponseEngine instantiation.');
%     
%     % Optional sub-struct: 'rgcMosaicParams'
%     % Assert that we have a sub-struct named 'coneMosaicParams'
%     if (isfield(paramsStruct, 'rgcMosaicParams'))
%         assert((isstruct(paramsStruct.rgcMosaicParams)), ...
%             'Optional ''rgcMosaicParams''passed during neuralResponseEngine instantiation is not a struct.');
%    end
    
    % Set the neural params
    obj.neuralParams = paramsStruct;
end
