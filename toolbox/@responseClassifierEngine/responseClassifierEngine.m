classdef responseClassifierEngine < handle
% Define a responseClassifierEngine class
%
% Syntax:
%   theClassifierEngine =
%      responseClassifierEngine(classifierComputeFunctionHandle, classifierParamsStruct)
%
% Description:
%    The responseClassifierEngine stores the trained classifier and any proprocessing 
%    constants returned by its computeFunction while in 'train' mode, so that when its 
%    computeFunction is called in the 'predict' operation mode, it can have access to them.
%
% Inputs:
%    classifierComputeFunctionHandle  - Function handle to the computeFunction that defines the
%                                       operation of the classifier
%
%    classifierParamsStruct           - Struct with parameters specific to the computeFunction. 
%                                       Optional. If not defined, the default params
%                                       defined in the computeFunction are used
% Outputs:
%    The created responseClassifierEngine object.
%
% Optional key/value pairs: None
%
%
% See Also:
%    t_responseClassifier.m, rcePcaSVMClassifier.m
%

% History:
%    9/20/2020  NPC Wrote it

    %% Public properties
    properties

    end
    
    %% Private properties
    properties (SetAccess=private)
        % User-supplied compute function handle to the classifier computation routine
        classifierComputeFunction
        % User-supplied struct with all classifier computation params
        classifierParams
        
        % The trained classifier
        trainedClassifier = [];
        
        % Preprocessing params from training dataset
        preProcessingConstants = [];
        
        % Valid modes of operation
        validOperationModes = {'train', 'predict'};
        
        % Required dataOut struct fields returned from the compute function during a 'train' operation mode
        requiredFieldsForTrainDataOutStruct = {'trainedClassifier', 'preProcessingConstants'};
        
        % Required dataOut struct fields returned from the compute function during a 'predict' operation mode
        requiredFieldsForPredictDataOutStruct = {'trialPredictions', 'pCorrect'};
    end
    
    % Public methods
    methods
        % Constructor
        function obj = responseClassifierEngine(classifierComputeFunctionHandle, classifierParamsStruct)
            % Validate and set the classifier compute function handle
            obj.validateAndSetComputeFunctionHandle(classifierComputeFunctionHandle);
            
            % If we dont receice a paramsStruct as the second argument use
            % the default params returned by the classifierComputeFunctionHandle
            if (nargin == 1)
                classifierParamsStruct = obj.classifierComputeFunction();
            end
            
            % Validate and set the classifier params
            obj.validateAndSetParamsStruct(classifierParamsStruct);
        end
        
        % Compute method
        dataOut = compute(obj, operationMode, nullResponses, testResponses);

    end
    
    % Private methods
    methods (Access = private)
        % Method to validate and set the classifier compute function handle
        validateAndSetComputeFunctionHandle(obj,classifierComputeFunctionHandle);
        % % Method to validate and set the classifier params struct
        validateAndSetParamsStruct(obj, classifierParamsStruct);
    end
    
end
