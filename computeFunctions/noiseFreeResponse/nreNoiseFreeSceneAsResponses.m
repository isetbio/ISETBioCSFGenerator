function dataOut = nreNoiseFreeSceneAsResposes(...
    neuralEngine, noiseFreeComputeParams, sceneSequence, ...
    sceneSequenceTemporalSupport, varargin)
% Compute function for computation of cone excitations witout eye movements 
%
% Syntax:
%   dataOut = nreSceneAsResponses(...
%    neuralEngine, neuralResponseParamsStruct, sceneSequence, ...
%    sceneSequenceTemporalSupport, varargin);
%
% Description:
%    Function serving as the computeFunctionHandle for a @neuralResponseEngine
%    object.  This is an illustrative function that just adds Poisson photon noise
%    to each pixel in the scene.
%
%    There are 2 ways to use this function.
%
%       [1] If called directly and with no arguments, 
%           dataOut = nreScenePhotonNoise()
%       it does not compute anything and simply returns a struct with the 
%       defaultParams (optics and coneMosaic params) that define the neural 
%       compute pipeline for this computation.
%
%       [2] If called from a parent @neuralResponseEngine object, 
%       it computes 'instancesNum' of cone photopigment excitation sequences 
%       in response to the passed 'sceneSequence'.
%
%    It is not a good idea to try to call this function with arguments
%    directly - it should be called by the compute method of its parent
%    @neuralResponseEngine.
%
% Inputs:
%    neuralEngine                   - the parent @neuralResponseEngine object that
%                                     is calling this function as its computeFunctionHandle
%    noiseFreeComputeParams         - a struct containing properties of the employed neural chain.
%                                     This should just be empty for this
%                                     compute function.
%    sceneSequence                  - a cell array of scenes defining the frames of a stimulus
%    sceneSequenceTemporalSupport   - the temporal support for the stimulus frames, in seconds
%
% Optional key/value input arguments:
%   None
%
% Outputs:
%    dataOut  - A struct that depends on the input arguments. 
%
%               If called directly with no input arguments, the returned struct contains
%               the defaultParams, which here is the empty struct.
%
%             - If called from a parent @neuralResponseEngine), the returned
%               struct is organized as follows:
%
%              .neuralResponses : matrix of neural responses.  The first dimension
%                                 is number of instances.  Often there is just one.
%                                 Each "column" of the second dimension 
%                                 is the response vector for one frame of the input.
%                                 The third dimension indexes the frames.
%              .temporalSupport : the temporal support of the neural
%                                   responses, in seconds
%              .noiseFreeResponsePipeline : a struct that the parent
%                                   @neuralResponseEngine
%                                   can tuck away.  Only returned if the
%                                   parent object has an empty parameters
%                                   struct. For this routine, it is just a
%                                   dummy struct as it isn't needed by this
%                                   routine.
% See Also:
%     t_neuralResponseCompute

% History:
%    09/26/2020  dhb  Wrote it.
%    10/19/2020  dhb  Fix comment to reflect fact that we now return
%                     instancesNum instances in noise free case.
%    12/18/2024  dhb  Rewrite for major architecture redo.

% Examples:
%{
    % Usage case #1. Just return the default neural response params
    % Here those are the empty matrix
    clear;
    defaultParams = nreNoiseFreeSceneAsResponses();

    % Usage case #2. Compute noise free, noisy, and repeatable (seed: 346) noisy response instances
    % using a parent @neuralResponseEngine object and the default neural response params

    % Instantiate the parent @neuralResponseEngine object.  Also need a
    % noise generating function
    theNeuralEngine = neuralResponseEngine(@nreNoiseFreeSceneAsResponses,@nreNoisyInstancesPoisson);

    % Compute a second version that doesn't add noise
    theNeuralEngineNoNoise = neuralResponseEngine(@nreNoiseFreeSceneAsResponses,@nreNoisyInstancesNoNoise);

    % Instantiate a @sceneEngine object and generate a test scene
    sceneParams = sceUniformFieldTemporalModulation;
    sceneParams.sizePixels = 5;
    theSceneEngine = sceneEngine(@sceUniformFieldTemporalModulation,sceneParams);
    testContrast = 0.1;
    [sceneSequence, sceneTemporalSupportSeconds] = ...
        theSceneEngine.compute(testContrast);
    
    % Compute noise free response
    [noiseFreeResponse, temporalSupportSeconds] = theNeuralEngine.computeNoiseFree(...
            sceneSequence, ...
            sceneTemporalSupportSeconds ...
            );

    % Add Poisson noise
    instancesNum = 16;
    [noisyInstances, temporalSupportSeconds] = theNeuralEngine.computeNoisyInstances(...
            noiseFreeResponse, ...
            temporalSupportSeconds, ...
            instancesNum, ...
            'random' ...
            );

    % Create response instances with no noise, really just to show that we
    % can do this with the nreNoisyInstancesNoNoise compute function.
    % Note that we can use the noise free responses generated by the version of
    % the nre that can add noise. We could also pass 'none' as the flag to
    % either this call or the one above to get instances without noise.
    [noNoiseInstances, ~] = theNeuralEngineNoNoise.computeNoisyInstances(...
            noiseFreeResponse, ...
            temporalSupportSeconds, ...
            instancesNum, ...
            'random' ...
            );
%}

    % Check input arguments. If called with zero input arguments, just return the default params struct
    if (nargin == 0)
        dataOut = generateDefaultParams();
        return;
    end
    
    % Parse the input arguments
    %   %
    % Allow for possibility that other nre's take key/value pairs that we
    % can ignore, so set KeepUnmatched to true.
    p = inputParser;
    p.KeepUnmatched = true;
    varargin = ieParamFormat(varargin);
    p.parse(varargin{:});
    
    % Return pipeline struct with key parameters  
    if (isempty(neuralEngine.neuralPipeline) | ~isfield(neuralEngine.neuralPipeline,'noiseFreeResponse'))
        noiseFreeResponsePipeline.dummyParameter = 1;
        frameDurationSeconds = noiseFreeComputeParams.frameDurationSeconds;
        returnTheNoiseFreePipeline = true;
    else
        frameDurationSeconds = neuralEngine.noiseFreeComputeParams.frameDurationSeconds;
        returnTheNoiseFreePipeline =  false;
    end

    % Temporal support for the neural response
    temporalSupportSeconds = sceneSequenceTemporalSupport; 
    if (length(temporalSupportSeconds) > 1)
        if (abs( (temporalSupportSeconds(2)-temporalSupportSeconds(1)) - frameDurationSeconds) > 1e-5)
            error('Specified frame duration not consistent with passed temporal support');
        end
    end
    
    % Compute the noise-free response.  The scene photons are photons/sec,
    % so we multiply everthing by the frame duration.
    responseDim = length(sceneSequence{1}.data.photons(:));
    framesNum = numel(sceneSequence);
    theNeuralResponses = zeros(1,responseDim,framesNum);
    for jj = 1:framesNum
        theNeuralResponses(1,:,jj) = sceneSequence{jj}.data.photons(:);
    end
    theNeuralResponses = theNeuralResponses*frameDurationSeconds;
    
    % Assemble the dataOut struct
    dataOut = struct(...
        'neuralResponses', theNeuralResponses, ...
    	'temporalSupport', temporalSupportSeconds ...
    );

    if (returnTheNoiseFreePipeline)
        dataOut.noiseFreeResponsePipeline = noiseFreeResponsePipeline;
    end
end

function p = generateDefaultParams()
    % Default params for this compute function
    p.frameDurationSeconds = 1;
end
