function dataOut = nreNoiseFreeMetaContrast(...
    metaNeuralEngineOBJ, metaNeuralParams, metaSceneSequence, ...
    metaSceneSequenceTemporalSupport, instancesNum, varargin)
% Meta contrast scene engine, noise free compute
%
% Syntax:
%   dataOut = nreNoiseFreeMetaContrast(...
%    neuralEngineOBJ, noiseFreeResponseParams, sceneSequence, ...
%    sceneSequenceTemporalSupport, instancesNum, varargin);
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
%   'theBackgroundRetinalImage'       - oi containing the computed retinal
%                                     image of the background.  Used by
%                                     some nre's to convert cone
%                                     excitations to modulation/contrast.
%                                     By default this is set to a stub oi
%                                     that doesn't do anything, so this
%                                     only needs to be passed explicitly if
%                                     the underlying nre needs it.
% 
% See Also:
%     t_metaContrastCSF

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
p.addParameter('theBackgroundRetinalImage', struct('type', 'opticalimage'), @isstruct);
varargin = ieParamFormat(varargin);
p.parse(varargin{:});

% Retrieve the response noiseFlag labels and validate them.
rngSeed = p.Results.rngSeed;
metaNeuralEngineOBJ.validateNoiseFlags(noiseFlags);

% For each noise flag we generate a corresponing neural response, and all
% neural responses are stored in a dictionary indexed by the noiseFlag label.
% Setup theNeuralResponses dictionary, loading empty responses for now
theNeuralResponses = containers.Map();
for idx = 1:length(noiseFlags)
    theNeuralResponses(noiseFlags{idx}) = [];
end

% We will use the neural pipeline to store the precomputed response
% vectors we need to then compute for any contrast
if (isempty(metaNeuralEngineOBJ.neuralPipeline))
    if (metaNeuralParams.contrast0 ~= 0)
        error('We assume that contrast0 is 0 in our meta contrast logic');
    end

    % Compute the contrast 0 response
    neuralPipeline.contrast0 = metaNeuralParams.contrast0;
    [actualSceneSequence, actualTemporalSupport] = metaNeuralParams.sceneEngine.compute(metaNeuralParams.contrast0);
    [contrast0Response, ~] = metaNeuralParams.neuralEngine.computeNoiseFree(...
        actualSceneSequence, ...
        actualTemporalSupport, ...
        'theBackgroundRetinalImage',p.Results.theBackgroundRetinalImage);

    neuralPipeline.response0 = contrast0Response;
    clear contrast0Response

    % Compute the passed contrast1 response
    neuralPipeline.contrast1 = metaNeuralParams.contrast1;
    [actualSceneSequence, actualTemporalSupport] = metaNeuralParams.sceneEngine.compute(metaNeuralParams.contrast1);
    [contrast1Response, ~] = metaNeuralParams.neuralEngine.computeNoiseFree(...
        actualSceneSequence, ...
        actualTemporalSupport, ...
        'theBackgroundRetinalImage',p.Results.theBackgroundRetinalImage);
    
    % Store in canonical form that we can compute the real response quickly as done
    % below.
    neuralPipeline.response1 = (contrast1Response-neuralPipeline.response0)/neuralPipeline.contrast1;
    clear contrast1Response
    
    % Save actual temporal support
    neuralPipeline.temporalSupportSecs = actualTemporalSupport;

    % Flag that we need to store what we computed here
    returnTheNeuralPipeline = true;
else
    returnTheNeuralPipeline = false;
end

%% Pipeline is initialized, can compute fast based on what has been precomputed.

% Set rng seed if one was passed. Not clear we need to do this because
% all the randomness is in the @coneMosaic compute object, but it
% doesn't hurt to do so, if we ever choose a random number at this
% level.
if (~isempty(rngSeed))
    oldSeed = rng(rngSeed);
end

% Compute noise free responses
theContrast = sceneGet(metaSceneSequence{1},'photons');
if (numel(theContrast) ~= 1)
    error('Something went wrong, scene photons should contain one number');
end
if (returnTheNeuralPipeline)
    theNoiseFreeNeuralResponses = neuralPipeline.response0 + ...
        theContrast*neuralPipeline.response1;
    theNoiseFreeTemporalSupport = neuralPipeline.temporalSupportSecs;
else
    theNoiseFreeNeuralResponses = metaNeuralEngineOBJ.neuralPipeline.response0 + ...
        theContrast*metaNeuralEngineOBJ.neuralPipeline.response1;
    theNoiseFreeTemporalSupport = metaNeuralEngineOBJ.neuralPipeline.temporalSupportSecs;
end

% Compute responses for each type of noise flag requested
for idx = 1:length(noiseFlags)
    % Return noise free response
    if (contains(ieParamFormat(noiseFlags{idx}), 'none'))

        % Store noise-free response instances.  By convention, we only 
        % return one instance here, independent of how many are requested.
        theNeuralResponses(noiseFlags{idx}) = theNoiseFreeNeuralResponses;

    else
        % We already have the noise free response, so now just need to add
        % noise to it to compute noisy response.
        %
        % If passed, set the rng seed in the hope that doing so will freeze the
        % noise.
        if (~isempty(rngSeed))
            rng(rngSeed);
        end
 
        % Add noise to each response instance
        [~,m,n] = size(theNoiseFreeNeuralResponses);
        theNoisyNeuralResponses = zeros(instancesNum,m,n);
        for iii = 1:instancesNum
            % useSeed = randi(32000,1,1);
            [theNoisyNeuralResponsesOneInstance] = metaNeuralParams.neuralEngine.compute(...
                theNoiseFreeNeuralResponses, ...
                theNoiseFreeTemporalSupport, ...
                1, ...
                'noiseFlags', {'random'}, ...
                'theBackgroundRetinalImage',p.Results.theBackgroundRetinalImage, ...
                'justAddNoise', true, 'rngSeed',[]);
            theNoisyNeuralResponses(iii,:,:) = theNoisyNeuralResponsesOneInstance('random');
        end
        theNeuralResponses(noiseFlags{idx}) = theNoisyNeuralResponses;
    end
end

% Restore rng seed if we set it
if (~isempty(rngSeed))
    rng(oldSeed);
end

% Assemble the dataOut struct
if (returnTheNeuralPipeline)
    dataOut = struct(...
    'neuralResponses', theNeuralResponses, ...
    'temporalSupport', neuralPipeline.temporalSupportSecs);
    dataOut.neuralPipeline = neuralPipeline;

else
     dataOut = struct(...
    'neuralResponses', theNeuralResponses);
     dataOut.temporalSupport = metaNeuralEngineOBJ.neuralPipeline.temporalSupportSecs;
end

end

% Default params for this compute function
function p = generateDefaultParams()
p = struct( ...
    'contrast0', 0, ...
    'contrast1', 0.5, ...
    'sceneEngine',[], ...
    'neuralEngine', [], ...
    'noiseAddMethod', 'nre' ...
    );
end

