function validateAndSetParamsStruct(obj, noiseFreeComputeParams,noisyInstancesComputeParams)
    
    % Accept an empty params struct
    if (isempty(noiseFreeComputeParams))
        % Set the noiseFreeNeuralParams to the default params struct, which is returned
        % by the neuralComputeFunction when it is called with no input arguments
        noiseFreeComputeParams = obj.noiseFreeComputeFunction();
    end
    
    % Assert that the input is a struct
    assert(isstruct(noiseFreeComputeParams), ...
        'Expected a struct during neuralResponseEngine instantiation');
    
    % Set the neural params
    obj.noiseFreeComputeParams = noiseFreeComputeParams;

        % Accept an empty params struct
    if (isempty(noisyInstancesComputeParams))
        % Set the noisyInstancesNeuralParams to the default params struct, which is returned
        % by the neuralComputeFunction when it is called with no input arguments
        noisyInstancesComputeParams = obj.noisyInstancesComputeFunction();
    end
    
    % Assert that the input is a struct
    assert(isstruct(noisyInstancesComputeParams), ...
        'Expected a struct during neuralResponseEngine instantiation');
    
    % Set the neural params
    obj.noisyInstancesComputeParams = noisyInstancesComputeParams;
end
