function [theSceneSequence, temporalSupportSeconds] = uniformFieldTemporalModulation(theSceneParamsStruct, testContrast)

    % Create an equal photon uniform field
    uniformScene = sceneCreate('uniform equal photon', theSceneParamsStruct.sizePixels);

    % Set the scene width in degrees
    uniformScene = sceneSet(uniformScene, 'wAngular', theSceneParamsStruct.fovDegs);

    % 1 meter away
    uniformScene = sceneSet(uniformScene, 'distance', 1.0);

    % background: adjust radiance according to desired  mean luminance
    meanLuminance = theSceneParamsStruct.meanLuminanceCdPerM2;
    backgroundScene = sceneAdjustLuminance(uniformScene, meanLuminance);

    % pedestal: adjust radiance according to desired  mean luminance and test contrast
    testLuminance = meanLuminance * (1.0 + testContrast);
    pedestalScene = sceneAdjustLuminance(uniformScene, testLuminance);
    
    % Generate temporal support for the scene sequence
    temporalSupportSeconds = (0:(theSceneParamsStruct.stimDurationFramesNum-1))*(theSceneParamsStruct.frameDurationSeconds);
    
    % Generate the scene sequence
    theSceneSequence = cell(1, theSceneParamsStruct.stimDurationFramesNum);
    for frameIndex = 1:theSceneParamsStruct.stimDurationFramesNum
        if (ismember(frameIndex, theSceneParamsStruct.stimOnsetFramesIndices))
            theSceneSequence{frameIndex} = pedestalScene;
        else
            theSceneSequence{frameIndex} = backgroundScene;
        end
    end
    
end