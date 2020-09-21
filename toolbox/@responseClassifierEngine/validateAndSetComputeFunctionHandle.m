function validateAndSetComputeFunctionHandle(obj,classifierComputeFunctionHandle)
    assert(isa(classifierComputeFunctionHandle, 'function_handle'), 'Expected a function handle during responseClassifierEngine instantiation');
    obj.classifierComputeFunction = classifierComputeFunctionHandle;
end



