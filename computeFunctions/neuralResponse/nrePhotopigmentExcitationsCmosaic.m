function dataOut = nrePhotopigmentExcitationsCmosaic(...
    neuralEngineOBJ, neuralResponseParamsStruct, sceneSequence, ...
    sceneSequenceTemporalSupport, instancesNum, varargin)
% Compute function for computation of cone excitations witout eye movements 
%
% Syntax:
%   dataOut = nrePhotopigmentExcitationsCmosaic(...
%    neuralEngineOBJ, neuralResponseParamsStruct, sceneSequence, ...
%    sceneSequenceTemporalSupport, instancesNum, varargin);
%
% Description:
%    Function serving as the computeFunctionHandle for a @neuralResponseEngine
%    object using the new @cMosaic object. There are 2 ways to use this function.
%
%       [1] If called directly and with no arguments, 
%           dataOut = nrePhotopigmentExcitationsConeMosaicHexWithNoEyeMovements()
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
%     t_neuralResponseCompute

% History:
%    03/29/2021  npc  Wrote it by adapting nrePhotopigmentExcitationsConeMosaicHexWithNoEyeMovements
%    04/10/2023  fh   Edited it so that with correctly set varargin, this
%                       function behaves the same as
%                       nrePhotopigmentExcitationsCmosaicWithNoEyeMovements.m

% Examples:
%{
    % Usage case #1. Just return the default neural response params
    defaultParams = nrePhotopigmentExcitationsCmosaic()

    % Usage case #2. Compute noise free, noisy, and repeatable (seed: 346) noisy response instances
    % using a parent @neuralResponseEngine object and the default neural response params

    % Instantiate the parent @neuralResponseEngine object
    theNeuralEngineOBJ = neuralResponseEngine(@nrePhotopigmentExcitationsCmosaic);

    % Instantiate a @sceneEngine object and generate a test scene
    theSceneEngineOBJ = sceneEngine(@sceUniformFieldTemporalModulation);
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

    % Get the number of scene sequences
    framesNum = numel(sceneSequence);
    % if there are more than one frames (scenes), then we check if they are
    % identical
    if framesNum ~= 1
        flag_detectDiff = checkSceneSequences_isidentical(sceneSequence);
        % if the scenes are different, throw an error message
        if flag_detectDiff
            warning('Note that the input scenes are temporally modulated.');
        end
    end
    
    % Parse the input arguments
    p = inputParser;
    p.addParameter('noiseFlags', {'random'});
    p.addParameter('rngSeed',[],@(x) (isempty(x) | isnumeric(x)));
    p.addParameter('amputateScenes', false, @islogical);
    varargin = ieParamFormat(varargin);
    p.parse(varargin{:});
    
    % Retrieve the response noiseFlag labels and validate them.
    noiseFlags = p.Results.noiseFlags;
    rngSeed = p.Results.rngSeed;
    amputateScenes = p.Results.amputateScenes;
    neuralEngineOBJ.validateNoiseFlags(noiseFlags);
    
    % For each noise flag we generate a corresponing neural response, and all 
    % neural responses are stored in a dictionary indexed by the noiseFlag label.
    % Setup theNeuralResponses dictionary, loading empty responses for now
    theNeuralResponses = containers.Map();
    for idx = 1:length(noiseFlags)
        theNeuralResponses(noiseFlags{idx}) = [];
    end
    
    if (isempty(neuralEngineOBJ.neuralPipeline))
        % Generate the @cMosaic object
        theConeMosaic = cMosaic(...
            'sizeDegs', neuralResponseParamsStruct.coneMosaicParams.sizeDegs, ...
            'eccentricityDegs', neuralResponseParamsStruct.coneMosaicParams.eccDegs, ...
            'integrationTime', neuralResponseParamsStruct.coneMosaicParams.timeIntegrationSeconds ...
            );
        
        % Generate optics appropriate for the mosaic's eccentricity
        oiEnsemble = theConeMosaic.oiEnsembleGenerate(neuralResponseParamsStruct.coneMosaicParams.eccDegs, ...
            'zernikeDataBase', 'Polans2015', ...
            'subjectID', neuralResponseParamsStruct.opticsParams.PolansSubject, ...
            'pupilDiameterMM', neuralResponseParamsStruct.opticsParams.pupilDiameterMM);
        theOptics = oiEnsemble{1};
        returnTheNeuralPipeline = true;
    else
        % Load the optics from the previously computed neural pipeline
        theOptics = neuralEngineOBJ.neuralPipeline.optics;
        % Load the cone mosaic from the previously computed neural pipeline
        theConeMosaic = neuralEngineOBJ.neuralPipeline.coneMosaic;
        returnTheNeuralPipeline =  false;
    end    

    % Compute the sequence of optical images corresponding to the sequence of scenes
    if framesNum == 1 
        theOIsequence = oiCompute(sceneSequence{1}, theOptics);
    else
        theListOfOpticalImages = cell(1, framesNum);
        for frame = 1:framesNum
            theListOfOpticalImages{frame} = oiCompute(sceneSequence{frame}, theOptics);
        end
        
        % Generate an @oiSequence object containing the list of computed optical images
        theOIsequence = oiArbitrarySequence(theListOfOpticalImages, sceneSequenceTemporalSupport);
    end
    
    % Set rng seed if one was passed. Not clear we need to do this because
    % all the randomness is in the @coneMosaic compute object, but it
    % doesn't hurt to do so, if we ever choose a random number at this
    % level.
    if (~isempty(rngSeed))
        oldSeed = rng(rngSeed);
    end
    
    % Compute responses for each type of noise flag requested
    for idx = 1:length(noiseFlags)
        if (contains(ieParamFormat(noiseFlags{idx}), 'none'))
            % Compute the noise-free response
            % To do so, first save the current mosaic noiseFlag
            lastConeMosaicNoiseFlag = theConeMosaic.noiseFlag;
            
            % Set the coneMosaic.noiseFlag to 'none';
            theConeMosaic.noiseFlag = 'none';
            
            % Compute noise-free response instances
            if ~amputateScenes
                [theNeuralResponses(noiseFlags{idx}), ~, ~, ~, temporalSupportSeconds] = ...
                    theConeMosaic.compute(theOIsequence, ...
                    'nTrials', instancesNum, ...
                    'withFixationalEyeMovements', false);
            else
                [theNeuralResponses(noiseFlags{idx}), ~, ~, ~, temporalSupportSeconds] = ...
                    theConeMosaic.compute(theOIsequence.frameAtIndex(1), ...
                    'nTrials', instancesNum, ...
                    'withFixationalEyeMovements', false, ...
                    'nTimePoints', framesNum);
            end
        
            % Restore the original noise flag
            theConeMosaic.noiseFlag = lastConeMosaicNoiseFlag;
            
        elseif (~isempty(rngSeed))
            % Compute noisy response instances with a specified random noise seed for repeatability
            if ~amputateScenes
                [~,theNeuralResponses(noiseFlags{idx}), ~, ~, temporalSupportSeconds] = ...
                    theConeMosaic.compute(theOIsequence, ...
                    'nTrials', instancesNum, ...
                    'withFixationalEyeMovements', false, ...
                    'seed', rngSeed);        % random seed
            else
                [~,theNeuralResponses(noiseFlags{idx}), ~, ~, temporalSupportSeconds] = ...
                    theConeMosaic.compute(theOIsequence.frameAtIndex(1), ...
                    'nTrials', instancesNum, ...
                    'withFixationalEyeMovements', false, ...
                    'nTimePoints', framesNum,...
                    'seed', rngSeed);        % random seed
            end
        elseif (contains(ieParamFormat(noiseFlags{idx}), 'random'))
            % Because computeForOISequence freezes noise, if we want
            % unfrozen noise (which is the case if we are here), 
            % we have to pass it a randomly chosen seed.
            useSeed = randi(32000,1,1);
            % Compute noisy response instances
            if ~amputateScenes
                [~,theNeuralResponses(noiseFlags{idx}), ~, ~, temporalSupportSeconds] = ...
                    theConeMosaic.compute(theOIsequence, ...
                    'nTrials', instancesNum, ...
                    'withFixationalEyeMovements', false, ...
                    'seed', useSeed);        % random seed
            else
                [~,theNeuralResponses(noiseFlags{idx}), ~, ~, temporalSupportSeconds] = ...
                    theConeMosaic.compute(theOIsequence.frameAtIndex(1), ...
                    'nTrials', instancesNum, ...
                    'withFixationalEyeMovements', false, ...
                    'nTimePoints', framesNum,...
                    'seed', useSeed);        % random seed
            end
        end
    end
    
    % Restore rng seed if we set it
    if (~isempty(rngSeed))
        rng(oldSeed);
    end
    
    % Assemble the dataOut struct
    dataOut = struct(...
        'neuralResponses', theNeuralResponses, ...
    	'temporalSupport', temporalSupportSeconds);
    
    if (returnTheNeuralPipeline)
        dataOut.neuralPipeline.optics = theOptics;
        dataOut.neuralPipeline.coneMosaic = theConeMosaic;
    end
end

function p = generateDefaultParams()
    % Default params for this compute function
    p = struct(...
        'opticsParams', struct(...
            'PolansSubject', 10, ...
            'pupilDiameterMM', 3.0 ...
        ), ...
        'coneMosaicParams', struct(...
            'sizeDegs', 0.3*[1 1], ...
            'eccDegs', [0 0], ...
            'timeIntegrationSeconds', 5/1000 ...
        ) ...
    );
end

function flag_detectDiff = checkSceneSequences_isidentical(sceneSequence)
    %the total number of frames
    total_seq = length(sceneSequence);
    %get the first scene (use it as reference)
    scenePhotons_ref = sceneSequence{1}.data.photons(:);
    %initialize the flag to be 0 for detecting any difference between 
    % scene sequences
    flag_detectDiff = 0; 
    %initialize the counter for the comparison scene
    counter_seq = 2;
    
    %select how many elements we'd like to randomly check
    nSamples    = 100;
    idx_samples = randi(length(scenePhotons_ref), [1,nSamples]);

    while (~flag_detectDiff) && (counter_seq <= total_seq)
        %check if all randomly selected elements are equal
        scenePhotons_comp = sceneSequence{counter_seq}.data.photons(:);
        meanAbsDiff = mean(abs(scenePhotons_ref(:) - scenePhotons_comp(:)));
        if any(meanAbsDiff) && (sum(abs(scenePhotons_ref(idx_samples)-...
                scenePhotons_comp(idx_samples)) < 1e-20*meanAbsDiff) ~= nSamples)
            flag_detectDiff = 1;
        end
        counter_seq = counter_seq+1;
    end
end
