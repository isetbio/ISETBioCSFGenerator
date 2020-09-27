function dataOut = photopigmentExcitationsWithNoEyeMovements(...
    neuralEngineOBJ, neuralResponseParamsStruct, sceneSequence, ...
    sceneSequenceTemporalSupport, instancesNum, varargin)
% Compute function for computation of cone excitations witout eye movements 
%
% Syntax:
%   dataOut = photopigmentExcitationsWithNoEyeMovements(...
%    neuralEngineOBJ, neuralResponseParamsStruct, sceneSequence, ...
%    sceneSequenceTemporalSupport, instancesNum, varargin);
%
% Description:
%    Compute function to be used as a computeFunctionHandle for a @neuralResponseEngine
%    object. There are 2 ways to use this function.
%
%       [1] If called directly and with no arguments, 
%
%           dataOut = photopigmentExcitationsWithNoEyeMovements()
%
%       it does not compute anything and simply returns a struct with the 
%       defaultParams (optics and coneMosaic params) that define the neural 
%       compute pipeline for this computation.
%
%       [2] If called from a parent @neuralResponseEngine object, 
%       it computes 'instancesNum' of cone photopigment excitation sequences 
%       in response to the passed 'sceneSequence'.
%
% Inputs:
%    neuralEngineOBJ                - the parent @neuralResponseEngine object,
%                                     calling this function as its computeFunctionHandle
%    neuralResponseParamsStruct     - a struct containing properties of the
%                                     employed neural chain.
%    sceneSequence                  - a cell array of scenes defining the frames of a stimulus
%
%    sceneSequenceTemporalSupport   - the temporal support for the stimulus frames, in seconds
%
%    instancesNum                   - the number of response instances to compute
%
%
% Optional key/value input arguments:
%    'noiseFlags'                   - Cell array of strings containing labels
%                                     that encode the type of noise to be included
%                                     Valid values are: 
%                                          - 'none' (noise-free responses)
%                                          - 'random' (random noisy response instances)
%                                          - 'rngSeed_someInt' (repeatable noisy response
%                                             instances controlled by the someInt rng seed) 
%
% Outputs:
%    dataOut  - A struct that depends on the input arguments. 
%
%               If called directly with no input arguments, the returned struct contains
%               the defaultParams (optics and coneMosaic) that define the neural 
%               compute pipeline for this computation.
%
%             - If called from a parent @neuralResponseEngine, the returned
%               struct is organized as follows:
%
%              .neuralResponses : dictionary of responses indexed with 
%                                   labels corresponding to the entries of
%                                   the 'noiseFlags'  optional argument
%              .temporalSupport : the temporal support of the neural
%                                   responses, in seconds
%              .theOptics       : the optics employed in the computation
%                                   (only returned if the parent @neuralResponseEngine
%                                   object does not yet have a value for its 'theOptics' property)
%              .theConeMosaic   : the coneMosaic employed in the computation
%                                   (only returned if the parent @neuralResponseEngine
%                                   object does not yet have a value for its 'theConeMosaic' property)
%
%
%       The computed neural responses can be extracted as:
%
%           neuralResponses('one of the entries of the noiseFlags input argument') 
%
%       and are arranged in a matrix of:
%
%           [instancesNum x mCones x tTimeBins] 
%
%
%
% See Also:
%     t_neuralResponseCompute

% History:
%    09/26/2020  NPC  Wrote it.
%
%   Examples:
%{
    % Usage case #1. Just return the default neural response params
    defaultParams = photopigmentExcitationsWithNoEyeMovements()

    % Usage case #2. Compute noise free, noisy, and repeatable (seed: 346) noisy response instances
    % using a parent @neuralResponseEngine object and the default neural response params

    % Instantiate the parent @neuralResponseEngine object
    theNeuralEngineOBJ = neuralResponseEngine(@photopigmentExcitationsWithNoEyeMovements);

    % Instantiate a @sceneEngine object and generate a test scene
    theSceneEngineOBJ = sceneEngine(@uniformFieldTemporalModulation);
    testContrast = 0.1;
    [theTestSceneSequence, theTestSceneTemporalSupportSeconds] = ...
        theSceneEngineOBJ.compute(testContrast);
    
    % Compute 16 response instances for a number of different noise flags
    instancesNum = 16;
    noiseFlags = {'random', 'none','rNgSeed346', 'rng Seed 100'};
    [theResponses, theResponseTemporalSupportSeconds] = theNeuralEngineOBJ.compute(...
            theTestSceneSequence, ...
            theTestSceneTemporalSupportSeconds, ...
            instancesNum, ...
            'noiseFlags', noiseFlags ...
            );

    % Retrieve the different computed responses
    noiseFreeResponses = theResponses('none');
    repeatableNoisyResponseInstances346 = theResponses('rNgSeed346');
    repeatableNoisyResponseInstances100 = theResponses('rng Seed 100');
    randomNoiseResponseInstances = theResponses('random');
%}

    % Default params for this compute function
    defaultNeuralResponseParams = generateDefaultParams();

    % Check input arguments. If called with zero input arguments, just return the default params struct
    if (nargin == 0)
        dataOut = defaultNeuralResponseParams;
        return;
    end
    
    % Parse the input arguments
    p = inputParser;
    p.addParameter('noiseFlags', {'random'});
    p.parse(varargin{:});
    
    % Retrieve the response noiseFlag labels and validate them.
    noiseFlags = p.Results.noiseFlags;
    neuralEngineOBJ.validateNoiseFlags(noiseFlags);
    
    % For each noise flag we generate a corresponing neural response, and all 
    % neural responses are stored in a dictionary indexed by the noiseFlag label.
    % Setup theNeuralResponses dictionary, loading empty responses for now
    theNeuralResponses = containers.Map();
    for idx = 1:length(noiseFlags)
        theNeuralResponses(noiseFlags{idx}) = [];
    end
    
    if (isempty(neuralEngineOBJ.theOptics))
        % Generate the optics
        theOptics = oiCreate(neuralResponseParamsStruct.opticsParams.type, neuralResponseParamsStruct.opticsParams.pupilDiameterMM);  
        returnTheOptics = true;
    else
        % Load the optics from the previous computations
        theOptics = neuralEngineOBJ.theOptics;
        returnTheOptics = false;
    end
    
    if (isempty(neuralEngineOBJ.theConeMosaic))
        % Generate the cone mosaic
        theConeMosaic = coneMosaicHex(neuralResponseParamsStruct.coneMosaicParams.upsampleFactor, ...
            'fovDegs', neuralResponseParamsStruct.coneMosaicParams.fovDegs, ...
            'integrationTime', neuralResponseParamsStruct.coneMosaicParams.timeIntegrationSeconds ...
        );
        returnTheConeMosaic = true;
    else
        % Load the cone mosaic from the previous computations
        theConeMosaic = neuralEngineOBJ.theConeMosaic;
        returnTheConeMosaic = false;
    end

    % Compute the sequence of optical images corresponding to the sequence of scenes
    framesNum = numel(sceneSequence);
    theListOfOpticalImages = cell(1, framesNum);
    for frame = 1:framesNum
        theListOfOpticalImages{frame} = oiCompute(sceneSequence{frame}, theOptics);
    end
    
    % Generate an @oiSequence object containing the list of computed optical images
    theOIsequence = oiArbitrarySequence(theListOfOpticalImages, sceneSequenceTemporalSupport);
    
    % Zero eye movements
    eyeMovementsNum = theOIsequence.maxEyeMovementsNumGivenIntegrationTime(theConeMosaic.integrationTime);
    emPaths = zeros(instancesNum, eyeMovementsNum, 2);
    
    % Compute responses for each type of noise flag requested
    for idx = 1:length(noiseFlags)
        
        if (contains(ieParamFormat(noiseFlags{idx}), 'none'))
            % Compute the noise-free response
            % To do so, first save the current mosaic noiseFlag
            lastConeMosaicNoiseFlag = theConeMosaic.noiseFlag;
            % Set the coneMosaic.noiseFlag to 'none';
            theConeMosaic.noiseFlag = 'none';
            % Compute noise-free response instances
            theNeuralResponses(noiseFlags{idx}) = theConeMosaic.computeForOISequence(theOIsequence, ...
                'emPaths', emPaths, ...   % the emPaths
                'currentFlag', false ...  % no photocurrent
            );
            % Restore the original noise flag
            theConeMosaic.noiseFlag = lastConeMosaicNoiseFlag;
            
        elseif (contains(ieParamFormat(noiseFlags{idx}), 'rngseed'))
            % Extract the seed from the noise flag
            rngSeed = str2double(strrep(ieParamFormat(noiseFlags{idx}), 'rngseed', ''));
            % Compute noisy response instances with a specified random noise seed for repeatability
            theNeuralResponses(noiseFlags{idx}) = theConeMosaic.computeForOISequence(theOIsequence, ...
                'emPaths', emPaths, ...    % the emPaths
                'currentFlag', false, ...  % no photocurrent
                'seed', rngSeed ...        % random seed
            );
        
        elseif (contains(ieParamFormat(noiseFlags{idx}), 'random'))
            % Compute noisy response instances
            theNeuralResponses(noiseFlags{idx}) = theConeMosaic.computeForOISequence(theOIsequence, ...
                'emPaths', emPaths, ...   % the emPaths
                'currentFlag', false ...  % no photocurrent
            );
        end
    end
    
    % Temporal support for the neural response
    temporalSupportSeconds = theConeMosaic.timeAxis; 
    
    % Assemble dataOut struct
    dataOut = struct(...
        'neuralResponses', theNeuralResponses, ...
    	'temporalSupport', temporalSupportSeconds);
    if (returnTheOptics)
        dataOut.theOptics = theOptics;
    end
    if (returnTheConeMosaic)
        dataOut.theConeMosaic = theConeMosaic;
    end

end

function p = generateDefaultParams()
    % Default params for this compute function
    p = struct(...
        'opticsParams', struct(...
            'type', 'wvf human', ...
            'pupilDiameterMM', 3.0 ...
        ), ...
        'coneMosaicParams', struct(...
            'upsampleFactor', 5, ...
            'fovDegs', 0.3, ...
            'timeIntegrationSeconds', 5/1000 ...
        ) ...
    );
end
