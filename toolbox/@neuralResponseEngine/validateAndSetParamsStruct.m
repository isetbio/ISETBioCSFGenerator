function validateAndSetParamsStruct(obj, paramsStruct)
    
    % Assert that the input is a struct
    assert(isstruct(paramsStruct), ...
        'Expected a struct during neuralResponseEngine instantiation');
    
    % Obligatory sub-truct: 'opticsParams'
    % Assert that we have a sub-struct named 'opticsParams'
    assert((isfield(paramsStruct, 'opticsParams'))&&(isstruct(paramsStruct.opticsParams)), ...
        'Expected an ''opticsParams'' substruct during neuralResponseEngine instantiation.');
    
    % Obligatory sub-struct: 'coneMosaicParams'
    % Assert that we have a sub-struct named 'coneMosaicParams'
    assert((isfield(paramsStruct, 'coneMosaicParams'))&&(isstruct(paramsStruct.coneMosaicParams)), ...
        'Expected a ''coneMosaicParams'' substruct during neuralResponseEngine instantiation.');
    
    % Optional sub-struct: 'rgcMosaicParams'
    % Assert that we have a sub-struct named 'coneMosaicParams'
    if (isfield(paramsStruct, 'rgcMosaicParams'))
        assert((isstruct(paramsStruct.rgcMosaicParams)), ...
            'Optional ''rgcMosaicParams''passed during neuralResponseEngine instantiation is not a struct.');
    end
    
    % Set the neural params
    obj.neuralParams = paramsStruct;
end
