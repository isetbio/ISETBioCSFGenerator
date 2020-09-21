classdef sceneGenerationEngine < handle
    %% Public properties
    properties

    end
    
    %% Private properties
    properties (SetAccess=private)
        % User-passed function handle to the scene computation routine
        sceneComputeFunction
        % User-passed struct with all scene computation params except contrast
        sceneParams
    end
    
    % Public methods
    methods
        % Constructor
        function obj = sceneGenerationEngine(sceneComputeFunctionHandle, sceneParamsStruct)
            % Validate and set the scene compute function handle
            obj.validateAndSetComputeFunctionHandle(sceneComputeFunctionHandle);
            
            % Validate and set the scene params
            obj.validateAndSetParamsStruct(sceneParamsStruct);
        end
        
        % Compute method
        function [theSceneSequence, temporalSupportSeconds] = compute(obj, sceneContrast)
            % Recompute the scene for the new scene contrast
            [theSceneSequence, temporalSupportSeconds] = obj.sceneComputeFunction(obj.sceneParams, sceneContrast);
        end
    end
    
    % Private methods
    methods (Access = private)
        % Method to validate and set the scene compute function handle
        validateAndSetComputeFunctionHandle(obj,sceneComputeFunctionHandle);
        % % Method to validate and set the scene params struct
        validateAndSetParamsStruct(obj,sceneParamsStruct);
    end
    
end