function [theNeuralResponse, temporalSupportSeconds, theOptics, theConeMosaic] = photopigmentExcitationsWithNoEyeMovements(...
    neuralEngineOBJ, theNeuralResponseParamsStruct, theSceneSequence, theSceneTemporalSupportSeconds, instancesNum)

    if (isempty(neuralEngineOBJ.theOptics))
        % Generate the optics
        theOptics = oiCreate(theNeuralResponseParamsStruct.opticsParams.type, theNeuralResponseParamsStruct.opticsParams.pupilDiameterMM);
    else
        theOptics = neuralEngineOBJ.theOptics;
    end
    
    if (isempty(neuralEngineOBJ.theConeMosaic))
        % Generate the cone mosaic
        theConeMosaic = coneMosaicHex(theNeuralResponseParamsStruct.coneMosaicParams.upsampleFactor, ...
            'fovDegs', theNeuralResponseParamsStruct.coneMosaicParams.fovDegs, ...
            'integrationTime', theNeuralResponseParamsStruct.coneMosaicParams.timeIntegrationSeconds ...
        );
    else
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
    
    % Compute the neural response (sequence of isomerizations here)
    theNeuralResponse = theConeMosaic.computeForOISequence(theOIsequence, ...
            'emPaths', emPaths, ...
            'currentFlag', false);
        
    % Return the temporal support for the neural response
    temporalSupportSeconds = theConeMosaic.timeAxis; 
end