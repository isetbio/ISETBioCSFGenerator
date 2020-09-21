function validateAndSetParamsStruct(obj, paramsStruct)
    assert(isstruct(paramsStruct), 'Expected a struct during sceneGenerationEngine instantiation.');
    obj.sceneParams = paramsStruct;
end

