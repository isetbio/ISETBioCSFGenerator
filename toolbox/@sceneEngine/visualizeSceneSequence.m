function visualizeSceneSequence(obj, sceneSequence, temporalSupportSeconds)
   
    scenesNum = numel(sceneSequence);
    RGBgunTrace = zeros(scenesNum,3);
    
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
     
    hFig = figure(); clf;
    set(hFig, 'Position', [100 400 1400 640], 'Color', [1 1 1]); 
    for frameIndex = 1:scenesNum
        subplot('Position', [0.01 0.07 0.45 0.9]);
        % Extract the XYZ image representation
        xyzImage = sceneGet(sceneSequence{frameIndex}, 'xyz');
        if (frameIndex == 1)
            xPixels = size(xyzImage,2);
            yPixels = size(xyzImage,1);
            x = 1:xPixels;
            y = 1:yPixels;
            x = x-mean(x);
            y = y-mean(y);
            no = round(xPixels/2);
            mo = round(yPixels/2);
        end
        
        % Linear RGB image
        displayLinearRGBimage = imageLinearTransform(xyzImage, displayXYZToLinearRGB);
        RGBgunTrace(frameIndex,:) = squeeze(displayLinearRGBimage(mo,no,:));
        % Settings RGB image
        displaySettingsImage = (ieLUTLinear(displayLinearRGBimage, displayGet(presentationDisplay, 'inverse gamma'))) / displayGet(presentationDisplay, 'nLevels');

        % Render image
        image(x,y,displaySettingsImage);
        % Cross hairs
        hold on;
        plot([0 0], max(abs(y))*[-1 1], 'k-');
        plot(max(abs(x))*[-1 1], [0 0],'k-');
        axis 'image';
        set(gca, 'FontSize', 16,  'XTick', [], 'YTick', []);
        title(sprintf('frame %d (%2.0f msec)', frameIndex,temporalSupportSeconds(frameIndex)*1000));
        

        % Render RGB gun modulation at center of stimulus
        subplot('Position', [0.50 0.07 0.45 0.9]); 
        plot(temporalSupportSeconds(1:frameIndex), RGBgunTrace(1:frameIndex,1), 'ro-', 'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5], 'LineWidth', 1.5); hold on;
        plot(temporalSupportSeconds(1:frameIndex), RGBgunTrace(1:frameIndex,2), 'go-', 'MarkerSize', 12, 'MarkerFaceColor', [0 1 0.5], 'LineWidth', 1.5);
        plot(temporalSupportSeconds(1:frameIndex), RGBgunTrace(1:frameIndex,3), 'bo-', 'MarkerSize', 12, 'MarkerFaceColor', [0.5 0.5 1], 'LineWidth', 1.5);
        hold off;
        set(gca, 'FontSize', 16, 'YLim', [0 1], 'XLim', [temporalSupportSeconds(1) temporalSupportSeconds(end)]);
        xlabel('time (seconds)');
        title('gun modulation at center of stimulus');
        drawnow;
    end
end