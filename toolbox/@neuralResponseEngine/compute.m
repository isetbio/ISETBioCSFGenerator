function [neuralResponses, temporalSupportSeconds] = compute(obj, ...
                sceneSequence, sceneTemporalSupportSeconds, instancesNum, varargin)
% Generic compute method for the @neuralResponseEngine class.
%
% Syntax:
%   [neuralResponses, temporalSupportSeconds] = compute(obj, ...
%                sceneSequence, sceneTemporalSupportSeconds, instancesNum, varargin)
%
% Description:
%    Generic compute method for the @neuralResponseEngine class. Its purpose 
%    is to provide a unified interface for computing the response of a neural 
%    processing pipeline independent of the pipeline's specifics. The code 
%    for building up and executing the neural response pipeline, e.g.
%         optics -> photon absorption -> phototransduction -> retinal ganglion cells
%    and the parameters of this pipeline (e.g, type of optics, pupil diameter, 
%    cone mosaic size, etc.) are defined in the computeFunctionHandle and the 
%    neuralResponseParamsStruct, respectively, both of  which are set when  
%    instantiating the @neuralResponse objethet.
%    Apart from providing a unified interface, this method also sets the
%    neuralPipeline struct of the parent @neuralResponseEngine object.
%
% Inputs:
%    obj                            - the parent @neuralResponseEngine object
%                               
%    sceneSequence                  - a cell array of scenes, representing a 
%                                     spatiotemporal stimulus
%
%    sceneTemporalSupportSeconds    - a vector of time stamps for each frame 
%                                     of the scene sequence 
%
%    instancesNum                   - a scalar, specifying the number of response 
%                                     instances to generate
%
%
% Optional key/value input arguments:
%    optional key/value pairs       - these are passed directly to the
%                                     computeFunction of the @neuralResponseEngine object
%
% Outputs:
%    neuralResponses                - computed neural responses
%
%    temporalSupportSeconds         - a vector of time stamps for each time bin of the computed response 
%
%
% See Also:
%     t_sceneGeneration

% History:
%    9/20/2020  NPC Wrote it

    % Call the user-supplied compute function
    dataOut = obj.neuralComputeFunction(obj, obj.neuralParams, sceneSequence, sceneTemporalSupportSeconds, instancesNum, varargin{:});

    % Parse dataOut struct
    neuralResponses = dataOut.neuralResponses;
    temporalSupportSeconds = dataOut.temporalSupport;

    % Set the neural pipeline struct for future computations
    if (isfield(dataOut, 'neuralPipeline'))
        obj.neuralPipeline = dataOut.neuralPipeline;
    end
             
end
        