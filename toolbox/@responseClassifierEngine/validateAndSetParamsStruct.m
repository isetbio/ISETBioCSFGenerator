function validateAndSetParamsStruct(obj, paramsStruct)
    assert(isstruct(paramsStruct), 'Expected a struct during responseClassifierEngine instantiation.');
    obj.classifierParams = paramsStruct;
end