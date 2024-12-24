function dataOut = nreNoiseFreeMetaContrast(...
    neuralEngine, metaNeuralParams, metaSceneSequence, ...
    metaSceneSequenceTemporalSupport, varargin)
% Meta contrast nre noise free compute function
%
% Syntax:
%   dataOut = nreNoiseFreeMetaContrast(...
%    metaNeuralEngine, noiseFreeResponseParams, sceneSequence, ...
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
%    metaNeuralEngine               - the parent @neuralResponseEngine object that
%                                     is calling this function as its computeFunctionHandle
%    noiseFreeComputeParams         - a struct containing properties of the employed neural chain.
%                                     This should just be empty for this
%                                     compute function.
%    sceneSequence                  - a cell array of scenes defining the frames of a stimulus
%    sceneSequenceTemporalSupport   - the temporal support for the stimulus frames, in seconds
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
%              .neuralResponses : matrix of neural responses.  Each column
%                                 is the response vector for one frame of the input.
%              .temporalSupport : the temporal support of the neural
%                                   responses, in seconds
%              .noiseFreeResponsePipeline : a struct that the parent
%                                   @neuralResponseEngine
%                                   can tuck away.  Only returned if the
%                                   parent object has an empty parameters
%                                   struct. For this routine, it has
%                                   information about the oi and cMosaic.
%
%
% Optional key/value input arguments:
%   None.
% 
% See Also:
%     nreNoisyInstancesMetaContrast, t_metaContrastCSF

% History:
%    08/11/24    dhb  Wrote it.
%    12/20/2024  dhb  Rewrite for major architecture redo.

% Check input arguments. If called with zero input arguments, just return the default params struct
if (nargin == 0)
    dataOut = generateDefaultParams();
    return;
end

% Parse the input arguments
p = inputParser;
p.KeepUnmatched = true;
varargin = ieParamFormat(varargin);
p.parse(varargin{:});

% We will use the neural pipeline to store the precomputed response
% vectors we need to then compute for any contrast
if (isempty(neuralEngine.neuralPipeline) | ~isfield(neuralEngine.neuralPipeline,'noiseFreeResponse'))
    if (metaNeuralParams.contrast0 ~= 0)
        error('We assume that contrast0 is 0 in our meta contrast logic');
    end

    % Compute the contrast 0 response
    noiseFreeResponsePipeline.contrast0 = metaNeuralParams.contrast0;
    [actualSceneSequence, actualTemporalSupport] = metaNeuralParams.sceneEngine.compute(metaNeuralParams.contrast0);
    [contrast0Response, ~] = metaNeuralParams.neuralEngine.computeNoiseFree(...
        actualSceneSequence, ...
        actualTemporalSupport, ...
        varargin{:});

    noiseFreeResponsePipeline.response0 = contrast0Response;
    clear contrast0Response

    % Compute the passed contrast1 response
    noiseFreeResponsePipeline.contrast1 = metaNeuralParams.contrast1;
    [actualSceneSequence, actualTemporalSupport] = metaNeuralParams.sceneEngine.compute(metaNeuralParams.contrast1);
    [contrast1Response, ~] = metaNeuralParams.neuralEngine.computeNoiseFree(...
        actualSceneSequence, ...
        actualTemporalSupport, ...
        varargin{:});
    
    % Store in canonical form that we can compute the real response quickly as done
    % below.
    noiseFreeResponsePipeline.response1 = ...
        (contrast1Response-noiseFreeResponsePipeline.response0)/noiseFreeResponsePipeline.contrast1;
    clear contrast1Response
    
    % Save actual temporal support
    noiseFreeResponsePipeline.temporalSupportSecs = actualTemporalSupport;

    % Flag that we need to store what we computed here
    returnTheNeuralPipeline = true;
else
    returnTheNeuralPipeline = false;
end

%% Pipeline is initialized, can compute fast based on what has been precomputed.
%
% Compute noise free responses
theContrast = sceneGet(metaSceneSequence{1},'photons');
if (numel(theContrast) ~= 1)
    error('Something went wrong, scene photons should contain one number');
end
if (returnTheNeuralPipeline)
    theNoiseFreeNeuralResponses = noiseFreeResponsePipeline.response0 + ...
        theContrast*noiseFreeResponsePipeline.response1;
    theNoiseFreeTemporalSupport = noiseFreeResponsePipeline.temporalSupportSecs;
else
    theNoiseFreeNeuralResponses = neuralEngine.neuralPipeline.noiseFreeResponse.response0 + ...
        theContrast*neuralEngine.neuralPipeline.noiseFreeResponse.response1;
    theNoiseFreeTemporalSupport = neuralEngine.neuralPipeline.noiseFreeResponse.temporalSupportSecs;
end

% Assemble the dataOut struct
if (returnTheNeuralPipeline)
    dataOut = struct(...
    'neuralResponses', theNoiseFreeNeuralResponses, ...
    'temporalSupport', noiseFreeResponsePipeline.temporalSupportSecs);
    dataOut.noiseFreeResponsePipeline = noiseFreeResponsePipeline;

else
     dataOut = struct(...
    'neuralResponses', theNoiseFreeNeuralResponses);
     dataOut.temporalSupport = neuralEngine.neuralPipeline.noiseFreeResponse.temporalSupportSecs;
end

end

% Default params for this compute function
function p = generateDefaultParams()
p = struct( ...
    'contrast0', 0, ...
    'contrast1', 0.5, ...
    'sceneEngine',[], ...
    'neuralEngine', [] ...
    );
end

