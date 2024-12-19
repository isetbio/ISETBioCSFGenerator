function [noiseFreeResponses, temporalSupportSeconds] = computeNoiseFree(obj, ...
                sceneSequence, sceneTemporalSupportSeconds, varargin)
% Generic compute method for the @neuralResponseEngine class.
%
% Syntax:
%   [noiseFreeResponses, temporalSupportSeconds] = compute(obj, ...
%                sceneSequence, sceneTemporalSupportSeconds, varargin)
%
% Description:
%    Generic compute noise free method for the @neuralResponseEngine class. Its purpose 
%    is to provide a unified interface for computing the noise free response of a neural 
%    processing pipeline independent of the pipeline's specifics. The code 
%    for building up and executing the neural response pipeline, e.g.
%         optics -> photon absorption -> phototransduction -> retinal ganglion cells
%    and the parameters of this pipeline (e.g, type of optics, pupil diameter, 
%    cone mosaic size, etc.) are defined in the computeFunctionHandle and the 
%    noiseFreeResponseParams, respectively, both of  which are set when  
%    instantiating the @neuralResponseEngine object.
%
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
% Optional key/value input arguments:
%    optional key/value pairs       - these are passed directly to the
%                                     computeFunction of the @neuralResponseEngine object
%
% Outputs:
%    noiseFreeResponses             - computed noise-free neural responses
%
%    temporalSupportSeconds         - a vector of time stamps for each time bin of the computed response 
%
% See Also:
%     t_sceneGeneration

% History:
%    9/20/2020  NPC Wrote it

    % Call the user-supplied compute function
    dataOut = obj.noiseFreeComputeFunction(obj, obj.noiseFreeComputeParams, sceneSequence, sceneTemporalSupportSeconds, varargin{:});

    % Parse dataOut struct
    noiseFreeResponses = dataOut.neuralResponses;
    temporalSupportSeconds = dataOut.temporalSupport;

    % Set the neural pipeline struct for future computations
    if (isfield(dataOut, 'noiseFreeResponsePipeline'))
        obj.neuralPipeline.noiseFreeResponses = dataOut.noiseFreeResponsePipeline;
    end
             
end
        