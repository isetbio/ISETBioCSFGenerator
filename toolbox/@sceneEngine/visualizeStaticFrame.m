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
    
    debugMode = ~true;
    if (debugMode)
        displaySettingsImage = sceneGet(sceneSequence{frameIndex}, 'rgb');
    else
        % Linear RGB image
        displayLinearRGBimage = imageLinearTransform(xyzImage, displayXYZToLinearRGB);
        
        % Find out-of-gamut pixels
        pixelsNum = numel(displayLinearRGBimage);
        belowZeroIdx = find(displayLinearRGBimage<0);
        overOneIdx = find(displayLinearRGBimage>1);
        
        % Deal with out-of-gamut pixels
        displayLinearRGBimage(belowZeroIdx) = 0;
        displayLinearRGBimage(overOneIdx) = 1;
        
        % Inform the user about the out-of-gamut pixels
        if ((numel(belowZeroIdx)>0) || (numel(overOneIdx)>0))
             fprintf(2,'Warning: out-of-gamut pixels in frame %d.\n', frameIndex);
        end
        if (numel(belowZeroIdx)>0)
            fprintf(2,'\t%2.1f of the image pixels  were < 0.\n', 100*numel(belowZeroIdx)/pixelsNum);
        end
        if (numel(overOneIdx)>0)
            fprintf(2,'\t%2.1f of the image pixels were > 1)\n', 100*numel(overOneIdx)/pixelsNum);
        end
        
        % Settings RGB image
        displaySettingsImage = (ieLUTLinear(displayLinearRGBimage, displayGet(presentationDisplay, 'inverse gamma'))) / displayGet(presentationDisplay, 'nLevels');
    end
    
    % Render image
    image(x,y,displaySettingsImage);
    % Cross hairs
    hold on;
    set(gca, 'XTick', [], 'YTick', []);
    drawnow;
end

end