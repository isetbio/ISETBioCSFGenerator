function [theNeuralResponses, temporalSupportSeconds, theOptics, theConeMosaic] = photopigmentExcitationsWithNoEyeMovements(...
    neuralEngineOBJ, theNeuralResponseParamsStruct, theSceneSequence, theSceneTemporalSupportSeconds, instancesNum, varargin)

    % Parse the passed key/value pairs
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
        theOptics = oiCreate(theNeuralResponseParamsStruct.opticsParams.type, theNeuralResponseParamsStruct.opticsParams.pupilDiameterMM);  
    else
        % Load the optics from the previous computations
        theOptics = neuralEngineOBJ.theOptics;
    end
    
    if (isempty(neuralEngineOBJ.theConeMosaic))
        % Generate the cone mosaic
        theConeMosaic = coneMosaicHex(theNeuralResponseParamsStruct.coneMosaicParams.upsampleFactor, ...
            'fovDegs', theNeuralResponseParamsStruct.coneMosaicParams.fovDegs, ...
            'integrationTime', theNeuralResponseParamsStruct.coneMosaicParams.timeIntegrationSeconds ...
        );  
    else
        % Load the cone mosaic from the previous computations
        theConeMosaic = neuralEngineOBJ.theConeMosaic;
    end

    % Compute the sequence of optical images corresponding to the sequence of scenes
    framesNum = numel(theSceneSequence);
    theListOfOpticalImages = cell(1, framesNum);
    for frame = 1:framesNum
        theListOfOpticalImages{frame} = oiCompute(theSceneSequence{frame}, theOptics);
    end
    
    % Generate an @oiSequence object containing the list of computed optical images
    theOIsequence = oiArbitrarySequence(theListOfOpticalImages, theSceneTemporalSupportSeconds);
    
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
    
    % Return the temporal support for the neural response
    temporalSupportSeconds = theConeMosaic.timeAxis; 
end