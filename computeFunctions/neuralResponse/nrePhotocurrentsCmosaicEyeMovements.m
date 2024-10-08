function dataOut = nrePhotocurrentsCmosaicEyeMovements(...
    neuralEngineOBJ, neuralResponseParamsStruct, sceneSequence, ...
    sceneSequenceTemporalSupport, instancesNum, varargin)
% Compute function for computation of photocurrent response in the presence of fixational eye movements 
%
% Syntax:
%   dataOut = nrePhotocurrentsCmosaicEyeMovements(...
%    neuralEngineOBJ, neuralResponseParamsStruct, sceneSequence, ...
%    sceneSequenceTemporalSupport, instancesNum, varargin);
%
% Description:
%    Function serving as the computeFunctionHandle for a @neuralResponseEngine
%    object using the new @cMosaic object. There are 2 ways to use this function.
%
%       [1] If called directly and with no arguments, 
%           dataOut = nrePhotocurrentsCmosaicEyeMovements()
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
%
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
% See Also:
%     t_neuralResponseCompute

% History:
%    03/29/2021  npc  Wrote it by adapting nrePhotopigmentExcitationsConeMosaicHexWithNoEyeMovements

% Examples:
%{
    % Usage case #1. Just return the default neural response params
    defaultParams = nrePhotocurrentsCmosaicEyeMovements()

    % Usage case #2. Compute noise free, noisy, and repeatable (seed: 346) noisy response instances
    % using a parent @neuralResponseEngine object and the default neural response params

    % Instantiate the parent @neuralResponseEngine object
    theNeuralEngineOBJ = neuralResponseEngine(@nrePhotocurrentsCmosaicEyeMovements);

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
    
    % Parse the input arguments
    p = inputParser;
    p.addParameter('noiseFlags', {'random'});
    p.addParameter('rngSeed',[],@(x) (isempty(x) | isnumeric(x)));
    varargin = ieParamFormat(varargin);
    p.parse(varargin{:});
    
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

    fprintf('Computing %d instances of eye movements (drift model: ''%s'')\n', ...
        instancesNum, neuralResponseParamsStruct.eyeMovementsParams.driftModel);
    
    switch (neuralResponseParamsStruct.eyeMovementsParams.driftModel)
        case 'default'
            theConeMosaic.emGenSequence(neuralResponseParamsStruct.eyeMovementsParams.durationSeconds, ...
                'centerPaths', true, ...
                'microsaccadeType', 'none', ...
                'nTrials', instancesNum);
    
        case 'low velocity'
            % Generate fixational eye movements with mugh lower drift speed
            driftModelPositionNoiseStd = 0.0;
            theConeMosaic.emGenSequence(neuralResponseParamsStruct.eyeMovementsParams.durationSeconds, ...
                'centerPaths', true, ...
                'microsaccadeType', 'none', ...
                'driftModelPositionNoiseStd', driftModelPositionNoiseStd, ...
                'nTrials', instancesNum);
       
        case 'high velocity'
            % Generate fixational eye movements with mugh higher drift speed
            driftModelPositionNoiseStd = 0.70;
            theConeMosaic.emGenSequence(neuralResponseParamsStruct.eyeMovementsParams.durationSeconds, ...
                'centerPaths', true, ...
                'microsaccadeType', 'none', ...
                'driftModelPositionNoiseStd', driftModelPositionNoiseStd, ...
                'nTrials', instancesNum);
            
        otherwise
            error('Unknown driftModel: ''%s''.', neuralResponseParamsStruct.eyeMovementsParams.driftModel);
    end
    
    % Compute the sequence of optical images corresponding to the sequence of scenes
    framesNum = numel(sceneSequence);
    theListOfOpticalImages = cell(1, framesNum);
    for frame = 1:framesNum
        theListOfOpticalImages{frame} = oiCompute(theOptics, sceneSequence{frame},'padvalue','mean');
    end
    
    % Generate an @oiSequence object containing the list of computed optical images
    theOIsequence = oiArbitrarySequence(theListOfOpticalImages, sceneSequenceTemporalSupport);
    
    % Set rng seed if one was passed. Not clear we need to do this because
    % all the randomness is in the @coneMosaic compute object, but it
    % doesn't hurt to do so, if we ever choose a random number at this
    % level.
    if (~isempty(rngSeed))
        oldSeed = rng(rngSeed);
    end   
    timeSamplesNum = numel(theConeMosaic.fixEMobj.timeAxis);
    
    % Compute responses for each type of noise flag requested
    for idx = 1:length(noiseFlags)
        if (contains(ieParamFormat(noiseFlags{idx}), 'none'))
            % Compute the noise-free response
            % To do so, first save the current mosaic noiseFlag
            lastConeMosaicNoiseFlag = theConeMosaic.noiseFlag;
            
            % Set the coneMosaic.noiseFlag to 'none';
            theConeMosaic.noiseFlag = 'none';
            
            fprintf('Computing 1 noise-free response instance\n');
            % Compute noise-free response instances WITHOUT eye movements
            [noiseFreeConeExcitationResponses, ~, ~, ~, temporalSupportSeconds] = ...
                theConeMosaic.compute(theOIsequence.frameAtIndex(1), ...
                'nTrials', 1, ...
                'nTimePoints', timeSamplesNum ...
            );
        
            % Compute photocurrent responses from the noiseFreeConeExcitationResponses
            theNeuralResponses(noiseFlags{idx}) = CMosaicNreComputePhotocurrent(...
                noiseFreeConeExcitationResponses, temporalSupportSeconds, theConeMosaic.noiseFlag, ...
                neuralResponseParamsStruct.eyeMovementsParams.keptResponsesDurationSeconds);
            
            % Restore the original noise flag
            theConeMosaic.noiseFlag = lastConeMosaicNoiseFlag;
            
        elseif (~isempty(rngSeed))
            
            fprintf('Computing %d response instances\n', instancesNum);
            % Compute noise-free cone excitation response instances under 
            % fixational eye movements, with a specified random noise seed for repeatability
            [noiseFreeConeExcitationResponses, ~, ~, ~, temporalSupportSeconds] = theConeMosaic.compute(theOIsequence.frameAtIndex(1), ...
                'withFixationalEyeMovements', true, ...
                'seed', rngSeed ...        % random seed
            );
            
           % Compute photocurrent responses from the noiseFreeConeExcitationResponses
            theNeuralResponses(noiseFlags{idx}) = cMosaicNreComputePhotocurrent(...
                noiseFreeConeExcitationResponses, temporalSupportSeconds, theConeMosaic.noiseFlag, ...
                neuralResponseParamsStruct.eyeMovementsParams.keptResponsesDurationSeconds);         
            
        elseif (contains(ieParamFormat(noiseFlags{idx}), 'random'))
            % Because computeForOISequence freezes noise, if we want
            % unfrozen noise (which is the case if we are here), 
            % we have to pass it a randomly chosen seed.
            useSeed = randi(32000,1,1);
            
            fprintf('Computing %d response instances\n', instancesNum);
            % Compute noise-free cone excitation response instances under
            % fixational eye movements
            [noiseFreeConeExcitationResponses, ~, ~, ~, temporalSupportSeconds] = theConeMosaic.compute(theOIsequence.frameAtIndex(1), ...
                'withFixationalEyeMovements', true, ...
                'seed', useSeed ...        % random seed
            );
        
            % Compute photocurrent responses from the noiseFreeConeExcitationResponses
            theNeuralResponses(noiseFlags{idx}) = cMosaicNreComputePhotocurrent(...
                noiseFreeConeExcitationResponses, temporalSupportSeconds, theConeMosaic.noiseFlag, ...
                neuralResponseParamsStruct.eyeMovementsParams.keptResponsesDurationSeconds);
           
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
        ), ...
        'eyeMovementsParams', struct(...
            'durationSeconds', 300/1000, ...
            'driftModel', 'default', ... 
            'keptResponsesDurationSeconds', 100/1000)...
    );
end
