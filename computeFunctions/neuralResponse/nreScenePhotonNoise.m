function dataOut = nreScenePhotonNoise(...
    neuralEngineOBJ, neuralResponseParamsStruct, sceneSequence, ...
    sceneSequenceTemporalSupport, instancesNum, varargin)
% Compute function for computation of cone excitations witout eye movements 
%
% Syntax:
%   dataOut = nreScenePhotonNoise(...
%    neuralEngineOBJ, neuralResponseParamsStruct, sceneSequence, ...
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
%    neuralEngineOBJ                - the parent @neuralResponseEngine object that
%                                     is calling this function as its computeFunctionHandle
%    neuralResponseParamsStruct     - a struct containing properties of the employed neural chain.
%                                     This should just be empty.
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
%   'amputateScenes'                - Logical. Must now always be false.
%                                     Will go away sooner or later.
%   'theBackgroundRetinalImage'     - OI describing the retinal image to
%                                     the background stimulus.  Default is
%                                     an empty struct. This is not used by
%                                     this nre, but the key/value pair is
%                                     here because some nre's use it.
%   'justAddNoise'                  - Boolean.  Default false. If true,
%                                     treat the sceneSequence as a noise
%                                     free instance and add noise to it.
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
%              .neuralResponses : dictionary of responses indexed with 
%                                   labels corresponding to the entries of
%                                   the 'noiseFlags'  optional argument
%              .temporalSupport : the temporal support of the neural
%                                   responses, in seconds
%              .neuralPipeline  : a dummy struct that the parent
%                                   @neuralResponseEngine
%                                   can tuck away.  Only returned if the
%                                   parent object has an empty parameters
%                                   struct.
%
%       The computed neural responses can be extracted as:
%           neuralResponses('one of the entries of noiseFlags') 
%       and are arranged in a matrix of:
%           [instancesNum x mResponses x tTimeBins] 
%       where the number of responses is the number of pixels times the
%       number of wavelengths in the scene.
%
% See Also:
%     t_neuralResponseCompute

% History:
%    09/26/2020  dhb  Wrote it.
%    10/19/2020  dhb  Fix comment to reflect fact that we now return
%                     instancesNum instances in noise free case.

% Examples:
%{
    % Usage case #1. Just return the default neural response params
    defaultParams = nreScenePhotonNoise()

    % Usage case #2. Compute noise free, noisy, and repeatable (seed: 346) noisy response instances
    % using a parent @neuralResponseEngine object and the default neural response params

    % Instantiate the parent @neuralResponseEngine object
    theNeuralEngineOBJ = neuralResponseEngine(@nreScenePhotonNoise);

    % Instantiate a @sceneEngine object and generate a test scene
    sceneParams = sceUniformFieldTemporalModulation;
    sceneParams.sizePixels = 5;
    theSceneEngineOBJ = sceneEngine(@sceUniformFieldTemporalModulation,sceneParams);
    testContrast = 0.1;
    [theTestSceneSequence, theTestSceneTemporalSupportSeconds] = ...
        theSceneEngineOBJ.compute(testContrast);
    
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
    p.addParameter('justAddNoise',false,@islogical)
    varargin = ieParamFormat(varargin);
    p.parse(varargin{:});

    % Scene amputation really complicates the code, and it is probably not
    % conceptually good to allow it at the nre level. Disallowing and
    % providing message.
    if (p.Results.amputateScenes)
        fprintf('No longer supporting amputation of scene sequences\n');
        fprintf('If scene sequence amputation is needed, do it application specific calling code.\n');
        error('Throwing error to motivate this change.');
    end
    
    % Retrieve the response noiseFlag labels and validate them.
    noiseFlags = p.Results.noiseFlags;
    rngSeed = p.Results.rngSeed;
    neuralEngineOBJ.validateNoiseFlags(noiseFlags);
    
    % For each noise flag we generate a corresponing neural response, and all 
    % neural responses are stored in a dictionary indexed by the noiseFlag label.
    % Setup theNeuralResponses dictionary, loading empty responses for now
    theNeuralResponses = containers.Map();
    for idx = 1:length(noiseFlags)
        theNeuralResponses(noiseFlags{idx}) = [];
    end
    
    % Return a dummy neural pipeline struct to keep parent object happy if
    % it hasn't yet stored one.
    if (isempty(neuralEngineOBJ.neuralPipeline))
        neuralPipeline.dummyfield = 1;
        returnTheNeuralPipeline = true;
    else
        returnTheNeuralPipeline =  false;
    end

    % Set rng seed
    if (~isempty(rngSeed))
        oldSeed = rng(rngSeed);
    end
    
    % Compute responses for each type of noise flag requested
    responseDim = length(sceneSequence{1}.data.photons(:));
    framesNum = numel(sceneSequence);
    for idx = 1:length(noiseFlags)
        switch (ieParamFormat(noiseFlags{idx}))
            case 'none'
                % Compute the noise-free response
                theResponses = zeros(instancesNum,responseDim,framesNum);
                for jj = 1:framesNum
                    for ii = 1:instancesNum
                        theResponses(ii,:,jj) = sceneSequence{jj}.data.photons(:);
                    end
                end
                
            case 'random'
                % Compute noisy response instances
                theResponses = zeros(instancesNum,responseDim,framesNum);
                for jj = 1:framesNum
                    noiseFreeResponses = sceneSequence{jj}.data.photons(:);
                    for ii = 1:instancesNum
                        theResponses(ii,:,jj) = iePoisson(noiseFreeResponses);
                    end
                end
                
            otherwise
                error(sprintf('Unknown noise flag %s passed',ieParamFormat(noiseFlags{idx})));
        end
        
        % Tuck responses into container
        theNeuralResponses(noiseFlags{idx}) = theResponses;
    end
    
    % Restore
    if (~isempty(rngSeed))
        rng(oldSeed);
    end
    
    % Temporal support for the neural response
    temporalSupportSeconds = sceneSequenceTemporalSupport; 
    
    % Assemble the dataOut struct
    dataOut = struct(...
        'neuralResponses', theNeuralResponses, ...
    	'temporalSupport', temporalSupportSeconds);
    if (returnTheNeuralPipeline)
        dataOut.neuralPipeline = neuralPipeline;
    end
end

function p = generateDefaultParams()
    % Default params for this compute function
    p = struct([]);
end
