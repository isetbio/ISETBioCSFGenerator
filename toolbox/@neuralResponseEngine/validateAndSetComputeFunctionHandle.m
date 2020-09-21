function validateAndSetComputeFunctionHandle(obj,sceneComputeFunctionHandle)
    assert(isa(sceneComputeFunctionHandle, 'function_handle'), 'Expected a function handle during neuralResponseEngine instantiation');
    obj.neuralComputeFunction = sceneComputeFunctionHandle;
end
