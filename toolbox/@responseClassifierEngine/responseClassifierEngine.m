classdef responseClassifierEngine < handle
    %% Public properties
    properties

    end
    
    %% Private properties
    properties (SetAccess=private)
        % User-passed function handle to the classifier computation routine
        classifierComputeFunction
        % User-passed struct with all classifier computation params
        classifierParams
        
        % The trained classifier
        trainedClassifier = [];
        
        % Preprocessing params from training dataset
        preProcessingConstants = [];
        
        % Valid modes of operation
        validOperationModes = {'train', 'predict'};
    end
    
    % Public methods
    methods
        % Constructor
        function obj = responseClassifierEngine(classifierComputeFunctionHandle, classifierParamsStruct)
            % Validate and set the classifier compute function handle
            obj.validateAndSetComputeFunctionHandle(classifierComputeFunctionHandle);
            
            % If we dont receice a paramsStruct as the second argument use
            % the default params returned by the
            % classifierComputeFunctionHandle
            if (nargin == 1)
                classifierParamsStruct = obj.classifierComputeFunction();
            end
            
            % Validate and set the classifier params
            obj.validateAndSetParamsStruct(classifierParamsStruct);
        end
        
        % Compute method
        function dataOut = compute(obj, operationMode, nullResponses, testResponses)
            % Validate the operationMode
            assert(ismember(operationMode, obj.validOperationModes), sprintf('The passed responseClassifierEngine.compute() ''%s'' is invalid.', operationMode));
            
            % Call the user-supplied compute function
            dataOut = obj.classifierComputeFunction(obj, operationMode, nullResponses, testResponses, obj.classifierParams);
        
            % Set the trainedClassifier property for future predictions
            if (isfield(dataOut, 'trainedClassifier'))
                obj.trainedClassifier = dataOut.trainedClassifier;
                obj.preProcessingConstants = dataOut.preProcessingConstants;
            end
        end
    end
    
    % Private methods
    methods (Access = private)
        % Method to validate and set the classifier compute function handle
        validateAndSetComputeFunctionHandle(obj,classifierComputeFunctionHandle);
        % % Method to validate and set the classifier params struct
        validateAndSetParamsStruct(obj, classifierParamsStruct);
    end
    
end
