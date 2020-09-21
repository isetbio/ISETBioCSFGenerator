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
    end
    
    % Public methods
    methods
        % Constructor
        function obj = responseClassifierEngine(classifierComputeFunctionHandle, classifierParamsStruct)
            % Validate and set the scene compute function handle
            obj.validateAndSetComputeFunctionHandle(classifierComputeFunctionHandle);
            
            % Validate and set the scene params
            obj.validateAndSetParamsStruct(classifierParamsStruct);
        end
        
        % Compute method
        function [pCorrect, classificationData, decisionBoundary] = compute(obj, nullResponses, testResponses)
            % Run the classifier for these null and test responses
            [pCorrect, classificationData, decisionBoundary] = obj.classifierComputeFunction(obj.classifierParams, nullResponses, testResponses);
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
