classdef sceneEngine < handle
% Define a sceneEngine class
%
% Syntax:
%   theSceneEngine = sceneEngine(sceneComputeFunctionHandle,
%   sceneParamsStruct)
%
% Description:
%    The sceneEngine provides a unified way of generating scenes which are
%    identical in all aspects except for one parameter (contrast for now),
%    which is the parameter that we vary in order to obtain a detection
%    threshold.
%
%    This object has a useful writable property, visualizeEachCompute.
%    When set to true, it will produce a figure that shows the scene
%    sequence produced by the computeFunction, which is very useful for
%    checking and debugging.  It does slow things down though, so should be
%    turned off once everything is working.
%
% Inputs:
%    sceneComputeFunctionHandle     - Function handle to the
%                                     computeFunction that defines the
%                                     spatiotemporal aspects of the
%                                     generated scene sequence.  See the
%                                     example compute function
%                                     sceUniformTemporalModulation provided
%                                     with this toolbox for a description
%                                     of the full API for the compute
%                                     function.
%    sceneParamsStruct              - Struct with parameters specific to
%                                     the computeFunction. Optional. If not
%                                     defined or passed as empty, the default
%                                     params defined in the computeFunction
%                                     are used.
% Outputs:
%    The created sceneEngine object.
%
% Optional key/value pairs: None
%
%
% See Also:
%    t_sceneGeneration.m, sceUniformFieldTemporalModulation.m
%

% History:
%    09/20/2020  NPC Wrote it
%    10/17/2020  dhb Comments.

    %% Public properties
    properties
        % User-settable flag for visualizing the output of each compute() call
        visualizeEachCompute = false;
    end
    
    %% Private properties
    properties (SetAccess=private)
        % User-passed function handle to the scene computation routine
        sceneComputeFunction
        
        % User-passed struct with all scene computation params except contrast
        sceneParams
        
        % Required dataOut struct fields returned from the compute function
        requiredFieldsForDataOutStruct = {'sceneSequence', 'temporalSupport'};
        
        % The presentation display
        presentationDisplay = [];
    end
    
    % Public methods
    methods
        % Constructor
        function obj = sceneEngine(sceneComputeFunctionHandle, sceneParamsStruct)
            % Validate and set the scene compute function handle
            obj.validateAndSetComputeFunctionHandle(sceneComputeFunctionHandle);
            
            % If we dont receice a paramsStruct as the second argument use
            % the default params returned by the sceneComputeFunctionHandle
            if (nargin == 1 || isempty(sceneParamsStruct))
                sceneParamsStruct = obj.sceneComputeFunction();
            end
            
            % Validate and set the scene params struct
            obj.validateAndSetParamsStruct(sceneParamsStruct);
        end
        
        % Compute method
        [theSceneSequence, temporalSupportSeconds, statusReport] = compute(obj, sceneContrast);
        
        % Visualization method
        visualizeSceneSequence(obj, sceneSequence, temporalSupportSeconds, varargin);
        
        % Visualize one frame of the sequence 
        visualizeStaticFrame(obj, sceneSequence, temporalSupportSeconds, varargin);
        
    end
    
    % Private methods
    methods (Access = private)
        % Method to validate and set the scene compute function handle
        validateAndSetComputeFunctionHandle(obj,sceneComputeFunctionHandle);
        % % Method to validate and set the scene params struct
        validateAndSetParamsStruct(obj,sceneParamsStruct);
    end
    
    % Class methods (static)
    methods (Static)
    end
    
end