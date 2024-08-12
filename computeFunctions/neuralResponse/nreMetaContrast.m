function dataOut = nreMetaContrast(...
    metaNeuralEngineOBJ, metaNeuralParams, metaSceneSequence, ...
    metaSceneSequenceTemporalSupport, instancesNum, varargin)
% Meta contrast scene engine
%
% Syntax:
%   dataOut = nreMetaContrast(...
%    neuralEngineOBJ, neuralResponseParamsStruct, sceneSequence, ...
%    sceneSequenceTemporalSupport, instancesNum, varargin);
%
% Description:
%    Function serving as the computeFunctionHandle for a @neuralResponseEngine
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
%    'amputateScenes'               - Logical. Whether we want to just
%                                       select the first scene and generate
%                                       cone excitations and amputate the
%                                       rest. Default: false.
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
%           [instancesNum x mCones x tTimeBins]
%
%       NOTE: MATLAB always drops the last dimension of an matrix if that  is 1.
%             So if tBins is 1, the returned array will be [instancesNum x  mCones],
%             NOT [instancesNum x mCones x 1].
%
% See Also:
%     t_metaContrast

% History:
%    08/11/24    dhb  Wrote it.

% Examples:
%{
    % Instantiate a @sceneEngine object and generate a test scene
    theSceneEngineOBJ = sceneEngine(@sceUniformFieldTemporalModulation);
    testContrast = 0.1;
    [theTestSceneSequence, theTestSceneTemporalSupportSeconds] = ...
        theSceneEngineOBJ.compute(testContrast);

    % Usage case #1. Just return the default neural response params
    nreParams = nrePhotopigmentExcitationsCmosaic()
    nreParams.coneMosaicParams.timeIntegrationSeconds = ...
        theTestSceneTemporalSupportSeconds(2)-theTestSceneTemporalSupportSeconds(1);
    
    % Usage case #2. Compute noise free, noisy, and repeatable (seed: 346) noisy response instances
    % using a parent @neuralResponseEngine object and the default neural response params

    % Instantiate the parent @neuralResponseEngine object
    theNeuralEngineOBJ = neuralResponseEngine(@nrePhotopigmentExcitationsCmosaic,nreParams);
    
    % Compute 16 response instances for a number of different noise flags
    instancesNum = 16;
    noiseFlags = {'random', 'none'};
    [theResponses, theResponseTemporalSupportSeconds] = theNeuralEngineOBJ.compute(...
            theTestSceneSequence, ...
            theTestSceneTemporalSupportSeconds, ...
            instancesNum, ...
            'noiseFlags', noiseFlags ...
            );

    % Retrieve the different computed responses
    noiseFreeResponses = theResponses('none');
    randomNoiseResponseInstances = theResponses('random');
%}


%% NOISE, BACKGROUNDWTFIMAGE


% Check input arguments. If called with zero input arguments, just return the default params struct
if (nargin == 0)
    dataOut = generateDefaultParams();
    return;
end

% Parse the input arguments
p = inputParser;
p.addParameter('noiseFlags', {'random'});
p.addParameter('rngSeed',[],@(x) (isempty(x) | isnumeric(x)));
p.addParameter('amputateScenes', false, @islogical);
p.addParameter('theBackgroundRetinalImage', struct('type', 'opticalimage'), @isstruct);
varargin = ieParamFormat(varargin);
p.parse(varargin{:});

% We need amputation scenes parameter to match calling routines, but we
% never want it to be true
if (p.Results.amputateScenes)
    error('The amputateScenes option was a temporary hack.  Rewrite your code so you can set it to false.');
end

% Retrieve the response noiseFlag labels and validate them.
noiseFlags = p.Results.noiseFlags;
rngSeed = p.Results.rngSeed;
metaNeuralEngineOBJ.validateNoiseFlags(noiseFlags);

% Get the number of scene sequences
framesNum = numel(metaSceneSequence);

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
    [contrast0Response, ~] = metaNeuralParams.neuralEngine.compute(...
        actualSceneSequence, ...
        actualTemporalSupport, ...
        1, ...
        'noiseFlags', {'none'}, ...
        'amputateScenes', false, ...
        'theBackgroundRetinalImage',p.Results.theBackgroundRetinalImage);

    neuralPipeline.response0 = contrast0Response('none');
    clear contrast0Response

    % Compute the passed contrast response
    neuralPipeline.contrast1 = metaNeuralParams.contrast1;
    [actualSceneSequence, actualTemporalSupport] = metaNeuralParams.sceneEngine.compute(metaNeuralParams.contrast1);
    [contrast1Response, ~] = metaNeuralParams.neuralEngine.compute(...
        actualSceneSequence, ...
        actualTemporalSupport, ...
        1, ...
        'noiseFlags', {'none'}, ...
        'amputateScenes', false);
    %'theBackgroundRetinalImage',theBackgroundRetinalImage);
    
    % Store in canonical form that we can compute the real response quickly as done
    % below.
    neuralPipeline.response1 = contrast1Response('none')/neuralPipeline.contrast1 - neuralPipeline.response0;
    clear contrast1Response1
    
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
else
    theNoiseFreeNeuralResponses = metaNeuralEngineOBJ.neuralPipeline.response0 + ...
        theContrast*metaNeuralEngineOBJ.neuralPipeline.response1;
end

% Compute responses for each type of noise flag requested
for idx = 1:length(noiseFlags)
    % Return noise free response
    if (contains(ieParamFormat(noiseFlags{idx}), 'none'))

        % Store noise-free response instances
        theNeuralResponses(noiseFlags{idx}) = theNoiseFreeNeuralResponses;

    else
        if (~isempty(rngSeed))
            % Compute noisy response instances with a specified random noise seed for repeatability
            theNeuralResponses(noiseFlags{idx}) = theNoiseFreeNeuralResponses; 
        else
            % Because computeForOISequence freezes noise, if we want
            % unfrozen noise (which is the case if we are here),
            % we have to pass it a randomly chosen seed.
            useSeed = randi(32000,1,1);

            % Compute noisy response instances with a random random noise seed 
            theNeuralResponses(noiseFlags{idx}) = theNoiseFreeNeuralResponses;
        end
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

