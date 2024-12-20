function dataOut = nreNoiseFreeMetaContrast(...
    metaNeuralEngineOBJ, metaNeuralParams, metaSceneSequence, ...
    metaSceneSequenceTemporalSupport, instancesNum, varargin)
% Meta contrast scene engine
%
% Syntax:
%   dataOut = nreNoiseFreeMetaContrast(...
%    neuralEngineOBJ, neuralResponseParamsStruct, sceneSequence, ...
%    sceneSequenceTemporalSupport, instancesNum, varargin);
%
% Description:
%    Function serving as the noise free compute function handle for a @neuralResponseEngine
%    object for the meta contrast scheme. There are 2 ways to use this function.
%
%       [1] If called directly and with no arguments,
%           dataOut = nrePhotopigmentExcitationsConeMosaicHexWithNoEyeMovements()
%       it does not compute anything and simply returns a struct with the
%       defaultParams that define the neural
%       compute pipeline for this computation.
%
%       [2] If called from a parent @neuralResponseEngine object,
%       it computes 'instancesNum' of response sequences
%       in response to the passed 'sceneSequence'.
%
%    It is not a good idea to try to call this function with arguments
%    directly - it should be called by the compute method of its parent
%    @neuralResponseEngine.
%
% Inputs:
%    neuralEngineOBJ                - the parent @neuralResponseEngine object that
%                                     is calling this function as its computeFunctionHandle
%    neuralResponseParamsStruct     - a struct containing properties of the
%                                     employed neural chain.
%    sceneSequence                  - a cell array of scenes defining the frames of a stimulus
%    sceneSequenceTemporalSupport   - the temporal support for the stimulus frames, in seconds
%    instancesNum                   - the number of response instances to compute
%
% Optional key/value input arguments:
%    'noiseFlags'                   - Cell array of strings containing labels
%                                     that encode the type of noise to be included
%                                     Valid values are:
%                                        - 'none' (noise-free responses)
%                                        - 'random' (noisy response instances)
%                                     Default is {'random'}.
%   'rngSeed'                       - Integer.  Set rng seed. Empty (default) means don't touch the
%                                     seed.
%
% Outputs:
%    dataOut  - A struct that depends on the input arguments.
%
%               If called directly with no input arguments, the returned struct contains
%               the defaultParams (optics and coneMosaic) that define the neural
%               compute pipeline for this computation.  This can be useful
%               for a user interested in knowing what needs to be supplied
%               to this.
%
%             - If called from a parent @neuralResponseEngine), the returned
%               struct is organized as follows:
%                .neuralResponses : dictionary of responses indexed with
%                                   labels corresponding to the entries of
%                                   the 'noiseFlags'  optional argument
%                .temporalSupport : the temporal support of the neural
%                                   responses, in seconds
%                .neuralPipeline  : a struct containing the optics and cone mosaic
%                                   employed in the computation (only returned if
%                                   the parent @neuralResponseEngine object has
%                                   an empty neuralPipeline property)
%
%       The computed neural responses can be extracted as:
%           neuralResponses('one of the entries of noiseFlags') 
%       and are arranged in a matrix of:
%           [instancesNum x nTimeBins x mResponses]
%       where mResponses is the dimensionality of the response vector returned
%       by the actual underlying nre (e.g. mCones in the case of a cMosaic).
%
%       BUT: If you ask for multiple instances in the noise free case, you
%            only get one instance.  An exception to may be if the
%            underlying nre involves a set of eye movement paths that are
%            not all zeros, where one noise free instance is returned per
%            eye movement path.
%
%       NOTE: MATLAB always drops the last dimension of an matrix that has
%             more than two dimensions, if that dimension has only 1 entry.
%             So if mCones is 1, the returned array will be [instancesNum x
%             nTimeBins], NOT [instancesNum x nTimeBins x mResponses].
%
% See Also:
%     t_metaContrastCSF

% History:
%    08/11/24    dhb  Wrote it.

% Check input arguments. If called with zero input arguments, just return the default params struct
if (nargin == 0)
    dataOut = generateDefaultParams();
    return;
end

% Parse the input arguments
p = inputParser;
p.addParameter('rngSeed',[],@(x) (isempty(x) | isnumeric(x)));
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

    neuralPipeline.response0 = contrast0Response
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

