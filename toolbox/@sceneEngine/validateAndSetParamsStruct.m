function validateAndSetParamsStruct(obj, paramsStruct)
    assert(isstruct(paramsStruct), 'Expected a struct during sceneGenerationEngine instantiation.');
    obj.sceneParams = paramsStruct;

    if (isfield(paramsStruct, 'sceneEngineName'))
        obj.name = paramsStruct.sceneEngineName;
    end

end

