function validateAndSetComputeFunctionHandle(obj, noiseFreeComputeFunctionHandle, noisyInstancesComputeFunctionHandle)
    assert(isa(noiseFreeComputeFunctionHandle, 'function_handle'), 'Expected a function handle during neuralResponseEngine instantiation');
    obj.noiseFreeComputeFunction = noiseFreeComputeFunctionHandle;

    assert(isa(noisyInstancesComputeFunctionHandle, 'function_handle'), 'Expected a function handle during neuralResponseEngine instantiation');
    obj.noisyInstancesComputeFunction = noisyInstancesComputeFunctionHandle;
end
