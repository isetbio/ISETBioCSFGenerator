classdef sceneEngine < handle
% Define a sceneEngine class
%
% Syntax:
%   theSceneEngine = sceneEngine(sceneComputeFunctionHandle, sceneParamsStruct)
%
% Description:
%    The sceneEngine provides a unified way of generating scenes which are 
%    identical in all aspects except for one parameter (contrast for now), 
%    which is the parameter that we vary in order to obtain a detection
%    threshold.
%
%
% Inputs:
%    sceneComputeFunctionHandle     - Function handle to the computeFunction that defines 
%                                     the spatiotemporal aspects of the generated scene sequence
%
%    sceneParamsStruct              - Struct with parameters specific to the computeFunction. 
%                                     Optional. If not defined, the default params
%                                     defined in the computeFunction are used
% Outputs:
%    The created sceneEngine object.
%
% Optional key/value pairs: None
%
%
% See Also:
%    t_sceneGeneration.m, uniformFIeldTemporalModulation.m
%

% History:
%    9/20/2020  NPC Wrote it

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
        function obj = sceneEngine(sceneComputeFunctionHandle, sceneParamsStruct)
            % Validate and set the scene compute function handle
            obj.validateAndSetComputeFunctionHandle(sceneComputeFunctionHandle);
            
            % If we dont receice a paramsStruct as the second argument use
            % the default params returned by the sceneComputeFunctionHandle
            if (nargin == 1)
                sceneParamsStruct = obj.sceneComputeFunction();
            end
            
            % Validate and set the scene params struct
            obj.validateAndSetParamsStruct(sceneParamsStruct);
        end
        
        % Compute method
        function [theSceneSequence, temporalSupportSeconds] = compute(obj, sceneContrast)
            % Call the user-supplied compute function
            dataOut = obj.sceneComputeFunction(sceneContrast, obj.sceneParams);
        
            % Parse dataOut struct
            theSceneSequence = dataOut.sceneSequence;
            temporalSupportSeconds = dataOut.temporalSupport;
            
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