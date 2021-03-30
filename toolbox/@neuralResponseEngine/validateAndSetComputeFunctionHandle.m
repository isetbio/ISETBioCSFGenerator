function validateAndSetComputeFunctionHandle(obj,neuralComputeFunctionHandle)
    assert(isa(neuralComputeFunctionHandle, 'function_handle'), 'Expected a function handle during neuralResponseEngine instantiation');
    obj.neuralComputeFunction = neuralComputeFunctionHandle;
end
