function validateAndSetComputeFunctionHandle(obj,sceneComputeFunctionHandle)
    assert(isa(sceneComputeFunctionHandle, 'function_handle'), 'Expected a function handle during sceneGenerationEngine instantiation');
    obj.sceneComputeFunction = sceneComputeFunctionHandle;
end

