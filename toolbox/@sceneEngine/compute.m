function [theSceneSequence, temporalSupportSeconds] = compute(obj, sceneContrast)
% Generic compute method for the @sceneEngine class.
%
% Syntax:
%   [theSceneSequence, temporalSupportSeconds] = compute(obj, sceneContrast)
%
% Description:
%    Generic compute method for the @sceneEngine class. Its purpose is to 
%    provide a unified interface for computing a temporal sequence of scenes
%    independent of scene specifics, at a given contrast. The code for generating
%    the scene and the scene parameters are defined in the computeFunctionHandle 
%    and the sceneParamsStruct, respectively, both of which are set when 
%    instantiating the @sceneEngine object.
%
% Inputs:
%    obj                      - the parent @sceneEngine object                              
%    sceneContrast            - the contrast of the scene to be generated
%
% Optional key/value input arguments: none 
%
% Outputs:
%    theSceneSequence        - a cell array of scenes, representing a spatio-temporal stimulus
%    temporalSupportSeconds  - a vector of time stamps for each frame of the scene sequence 
%
% See Also:
%     t_sceneGeneration, t_modulatedGratingsSceneGeneration,
%     sceUniformFieldModulation, sceGrating

% History:
%    09/20/2020  NPC Wrote it
%    10/18/2020  dhb Pass obj to scene compute function, just in case we
%                    ever need to give the function access to object
%                    properties.

    % Call the user-supplied compute function
    dataOut = obj.sceneComputeFunction(obj, sceneContrast, obj.sceneParams);

    % Retrieve the returned dataOut fields
    theFieldNames = fieldnames(dataOut);
    
    % Validate the dataOut struct
    for k = 1:numel(obj.requiredFieldsForDataOutStruct)
        assert(ismember(obj.requiredFieldsForDataOutStruct{k}, theFieldNames), sprintf('dataOut struct does not contain the ''%s'' field', obj.requiredFieldsForDataOutStruct{k}));
    end
        
    % Parse dataOut struct
    theSceneSequence = dataOut.sceneSequence;
    temporalSupportSeconds = dataOut.temporalSupport;
    
    if (isfield(dataOut, 'presentationDisplay'))
        obj.presentationDisplay = dataOut.presentationDisplay;
    end
end