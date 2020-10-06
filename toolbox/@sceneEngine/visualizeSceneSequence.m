function visualizeSceneSequence(obj, sceneSequence, temporalSupportSeconds)
    figure; clf;
    scenesNum = numel(sceneSequence);
    maxColsNum = 8;
    if (scenesNum <= maxColsNum)
        rowsNum = 1;
        colsNum = scenesNum;
    else
        colsNum = maxColsNum;
        rowsNum = ceil(scenesNum/colsNum);
    end
    
    if (isempty(obj.presentationDisplay))
       % Generate generic display
       presentationDisplay = displayCreate('LCD-Apple');
    else
       % Use employed display
       presentationDisplay = obj.presentationDisplay;
    end
    % Compute the RGB settings for the display
    displayLinearRGBToXYZ = displayGet(presentationDisplay, 'rgb2xyz');
    displayXYZToLinearRGB = inv(displayLinearRGBToXYZ);
    
        
    for frameIndex = 1:scenesNum
        % Extract the XYZ image representation
        xyzImage = sceneGet(sceneSequence{frameIndex}, 'xyz');
        % Linear RGB image
        displayLinearRGBimage = imageLinearTransform(xyzImage, displayXYZToLinearRGB);
        % Settings RGB image
        displaySettingsImage = (ieLUTLinear(displayLinearRGBimage, displayGet(presentationDisplay, 'inverse gamma'))) / displayGet(presentationDisplay, 'nLevels');

        % Display
        image(displaySettingsImage);
        axis 'image';
        title(sprintf('frame %d\n(%2.0f msec)', frameIndex,temporalSupportSeconds(frameIndex)*1000));
        set(gca, 'XTick', [], 'YTick', []);
        drawnow;
    end
    
end