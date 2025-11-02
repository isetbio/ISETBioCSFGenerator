function visualizeSceneSequence(obj, sceneSequence, temporalSupportSeconds, varargin)
    p = inputParser;
    p.addParameter('videoFilename', [], @(x)(isempty(x)||ischar(x)));
    p.addParameter('sRGBforSceneVisualization', false, @islogical);
    parse(p, varargin{:});
    videoFileName = p.Results.videoFilename;
    sRGBforSceneVisualization = p.Results.sRGBforSceneVisualization;

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
    displayLinearRGBToLMS = displayGet(presentationDisplay, 'rgb2lms');

    if (~isempty(videoFileName))
        videoOBJ = VideoWriter(videoFileName, 'MPEG-4');
        videoOBJ.FrameRate = 10;
        videoOBJ.Quality = 100;
        videoOBJ.open();
    end

    hFig = figure(999); clf;
    set(hFig, 'Position', [100 400 1400 640], 'Color', [1 1 1]); 
    ax = subplot('Position', [0.01 0.07 0.45 0.85]);
    axModulation = subplot('Position', [0.50 0.07 0.45 0.85]);
        

    for frameIndex = 1:scenesNum
        

        % Extract the XYZ image representation
        xyzImage = sceneGet(sceneSequence{frameIndex}, 'xyz');
        displaySettingsImage = sceneGet(sceneSequence{frameIndex}, 'rgb');

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

        displayLinearLMSimage = imageLinearTransform(displayLinearRGBimage,displayLinearRGBToLMS);
        Limage = squeeze(displayLinearLMSimage(:,:,1));
        Mimage = squeeze(displayLinearLMSimage(:,:,2));
        Simage = squeeze(displayLinearLMSimage(:,:,3));
        meanL = mean(Limage(:));
        meanM = mean(Mimage(:));
        meanS = mean(Simage(:));

        xChroma = xyzImage(:,:,1)./sum(xyzImage,3);
        yChroma = xyzImage(:,:,2)./sum(xyzImage,3);

        % Render image
        stimProfile = squeeze(sum(displayLinearRGBimage(mo,:,:),3));
        m1 = min(stimProfile);
        m2 = max(stimProfile);
        if (m1 == m2)
            normalizedStimProfile = 2*((stimProfile-m1)/(m2-m1)-0.5);
        else
            normalizedStimProfile = stimProfile*0;
        end

        if (sRGBforSceneVisualization)
            image(ax,x,y,lrgb2srgb(displaySettingsImage));
        else
            image(ax,x,y,displaySettingsImage);
        end

        hold(ax, 'on');
        
        % Profile of the R+B+G guns
        plot(ax, x, yPixels/2 * normalizedStimProfile, 'k.-');
        
        % Cross hairs
        %plot(ax, [0 0], max(abs(y))*[-1 1], 'k-');
        %plot(ax, max(abs(x))*[-1 1], [0 0],'k-');
        axis(ax, 'image');
        hold(ax, 'off');
        set(ax, 'CLim', [0 1], 'XLim', 0.5*xPixels*[-1 1], 'YLim', 0.5*yPixels*[-1 1]);
        set(ax, 'FontSize', 16,  'XTick', [], 'YTick', []);
        title(ax,sprintf('%s\nLMS: %2.3f,%2.3f,%2.3f, xyY = (%2.2f,%2.2f,%2.1f cd/m2) (t:%2.0f msec)', ...
            sceneGet(sceneSequence{frameIndex}, 'name'), ...
            meanL, meanM, meanS, ...
            mean(xChroma(:)), mean(yChroma(:)), sceneGet(sceneSequence{frameIndex}, 'mean luminance'), ...
            temporalSupportSeconds(frameIndex)*1000));
        
        % Render RGB gun modulation at center of stimulus
        plot(axModulation, temporalSupportSeconds(1:frameIndex), RGBgunTrace(1:frameIndex,1), 'ro-', 'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5], 'LineWidth', 1.5); 
        hold(axModulation, 'on');
        plot(axModulation, temporalSupportSeconds(1:frameIndex), RGBgunTrace(1:frameIndex,2), 'go-', 'MarkerSize', 12, 'MarkerFaceColor', [0 1 0.5], 'LineWidth', 1.5);
        plot(axModulation, temporalSupportSeconds(1:frameIndex), RGBgunTrace(1:frameIndex,3), 'bo-', 'MarkerSize', 12, 'MarkerFaceColor', [0.5 0.5 1], 'LineWidth', 1.5);
        hold(axModulation, 'off');
        set(axModulation, 'FontSize', 16, 'YLim', [min(min(RGBgunTrace(1:frameIndex,:))) max(max(RGBgunTrace(1:frameIndex,:)))], ...
                                 'XLim', [temporalSupportSeconds(1) temporalSupportSeconds(end)+eps]);
        xlabel(axModulation,'time (seconds)');
        title(axModulation,'RGB gun modulation at center of stimulus');
        drawnow;

        if (~isempty(videoFileName))
            videoOBJ.writeVideo(getframe(hFig));
        end
    end
    
    if (~isempty(videoFileName))
    	videoOBJ.close()
    end

end