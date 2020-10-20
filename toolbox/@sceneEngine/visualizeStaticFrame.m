function visualizeStaticFrame(obj, sceneSequence)

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
    if (frameIndex == 1)
        xPixels = size(xyzImage,2);
        yPixels = size(xyzImage,1);
        x = 1:xPixels;
        y = 1:yPixels;
        x = x-mean(x);
        y = y-mean(y);
    end
    
    % Linear RGB image
    displayLinearRGBimage = imageLinearTransform(xyzImage, displayXYZToLinearRGB);
    % Settings RGB image
    displaySettingsImage = (ieLUTLinear(displayLinearRGBimage, displayGet(presentationDisplay, 'inverse gamma'))) / displayGet(presentationDisplay, 'nLevels');
    
    % Render image
    image(x,y,displaySettingsImage);
    % Cross hairs
    hold on;
    set(gca, 'XTick', [], 'YTick', []);
    drawnow;
end

end