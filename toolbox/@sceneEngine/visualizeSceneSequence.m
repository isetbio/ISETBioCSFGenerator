function visualizeSceneSequence(obj, sceneSequence, temporalSupportSeconds, varargin)
% Render a multi-frame scene sequence as an animated figure.
%
% Left panel  - sRGB (or linear-RGB) image of the scene at each frame,
%               with an optional RGB-sum row-profile overlay (see
%               'showProfile').
% Right panel - RGB-gun modulation at a chosen pixel, grown frame by
%               frame (see 'whichPixel').
%
% Optionally writes the animation to an MPEG-4 file at 10 fps (quality 100).
%
% Syntax:
%   obj.visualizeSceneSequence(sceneSequence, temporalSupportSeconds)
%   obj.visualizeSceneSequence(..., 'Name', Value, ...)
%
% Inputs:
%   sceneSequence          - {1 x nFrames} cell array of iset scene structs
%   temporalSupportSeconds - [1 x nFrames] time stamp of each frame (sec)
%
% Optional key/value pairs:
%   'videoFilename'             - char path for MPEG-4 output (no extension
%                                 needed); default [] writes no file.
%   'sRGBforSceneVisualization' - logical; true (default) gamma-encodes the
%                                 displayed image via lrgb2srgb.
%   'showProfile'               - logical; true (default) overlays a
%                                 normalised RGB-sum profile along the
%                                 sampled row on the scene image.  Set to
%                                 false to keep the scene image uncluttered.
%   'whichPixel'                - [row col] pixel used for the RGB-gun
%                                 modulation trace on the right panel.
%                                 Default [] uses the image centre pixel.
%                                 Pass the stimulus centre pixel when the
%                                 stimulus is offset from image centre so
%                                 that the trace reflects actual modulation.

    p = inputParser;
    p.addParameter('videoFilename',             [],   @(x)(isempty(x)||ischar(x)));
    p.addParameter('sRGBforSceneVisualization', true, @islogical);
    p.addParameter('showProfile',               true, @islogical);
    p.addParameter('whichPixel',                [],   @(x)(isempty(x)||(isnumeric(x)&&numel(x)==2)));
    parse(p, varargin{:});

    videoFileName             = p.Results.videoFilename;
    sRGBforSceneVisualization = p.Results.sRGBforSceneVisualization;
    showProfile               = p.Results.showProfile;
    whichPixel                = p.Results.whichPixel;

    scenesNum   = numel(sceneSequence);
    RGBgunTrace = zeros(scenesNum, 3);

    % Display model: converts scene XYZ to linear RGB gun values
    if (isempty(obj.presentationDisplay))
        presentationDisplay = displayCreate('LCD-Apple');
    else
        presentationDisplay = obj.presentationDisplay;
    end
    displayLinearRGBToXYZ = displayGet(presentationDisplay, 'rgb2xyz');
    displayXYZToLinearRGB = inv(displayLinearRGBToXYZ);
    displayLinearRGBToLMS = displayGet(presentationDisplay, 'rgb2lms');

    if (~isempty(videoFileName))
        videoOBJ = VideoWriter(videoFileName, 'MPEG-4');
        videoOBJ.FrameRate = 10;
        videoOBJ.Quality   = 100;
        videoOBJ.open();
    end

    % Figure layout: scene image on left, modulation trace on right.
    % Narrow-band scenes use a separate figure number to avoid collisions.
    theSceneEngineName = obj.name;
    theSceneName = sceneGet(sceneSequence{1}, 'name');
    if (strcmp(theSceneName, 'narrow band'))
        hFig = figure(9000); clf;
        set(hFig, 'Position', [5000 10 1400 640], 'Color', [1 1 1]);
    else
        hFig = figure(1000); clf;
        set(hFig, 'Position', [100 400 1400 640], 'Color', [1 1 1]);
    end
    ax           = subplot('Position', [0.01 0.07 0.45 0.80]);
    axModulation = subplot('Position', [0.50 0.07 0.45 0.80]);

    for frameIndex = 1:scenesNum

        theSceneName         = sceneGet(sceneSequence{frameIndex}, 'name');
        xyzImage             = sceneGet(sceneSequence{frameIndex}, 'xyz');
        displaySettingsImage = sceneGet(sceneSequence{frameIndex}, 'rgb');

        if (frameIndex == 1)
            xPixels = size(xyzImage, 2);
            yPixels = size(xyzImage, 1);

            % Pixel coordinate axes centred at (0,0) for image/plot overlay alignment
            x = (1:xPixels) - mean(1:xPixels);
            y = (1:yPixels) - mean(1:yPixels);

            % Pixel sampled for both the profile overlay and the right-panel trace
            if isempty(whichPixel)
                sampleRow = round(yPixels / 2);   % image centre (default)
                sampleCol = round(xPixels / 2);
            else
                sampleRow = whichPixel(1);
                sampleCol = whichPixel(2);
            end
        end

        % Convert scene XYZ to display linear RGB
        displayLinearRGBimage = imageLinearTransform(xyzImage, displayXYZToLinearRGB);

        % Record RGB gun values at the chosen pixel for the right-panel trace
        RGBgunTrace(frameIndex, :) = squeeze(displayLinearRGBimage(sampleRow, sampleCol, :));

        % Mean LMS and xy chromaticity: used only in the frame title
        displayLinearLMSimage = imageLinearTransform(displayLinearRGBimage, displayLinearRGBToLMS);
        Limage = squeeze(displayLinearLMSimage(:,:,1));
        Mimage = squeeze(displayLinearLMSimage(:,:,2));
        Simage = squeeze(displayLinearLMSimage(:,:,3));
        meanL  = mean(Limage(:));
        meanM  = mean(Mimage(:));
        meanS  = mean(Simage(:));

        xChroma = xyzImage(:,:,1) ./ sum(xyzImage, 3);
        yChroma = xyzImage(:,:,2) ./ sum(xyzImage, 3);
        luma    = xyzImage(:,:,3);

        % Image-centre stats for the frame title (independent of whichPixel)
        midRow = round(yPixels / 2);
        midCol = round(xPixels / 2);
        centerL         = Limage(midRow, midCol);
        centerM         = Mimage(midRow, midCol);
        centerS         = Simage(midRow, midCol);
        centerChromaX   = xChroma(midRow, midCol, 1);
        centerChromaY   = yChroma(midRow, midCol, 1);
        centerLuminance = luma(midRow, midCol, 1);

        % --- Left panel: scene image ----------------------------------------
        if (sRGBforSceneVisualization)
            image(ax, x, y, lrgb2srgb(displaySettingsImage));
        else
            image(ax, x, y, displaySettingsImage);
        end
        hold(ax, 'on');

        % Optional row-profile overlay along the sampled row.
        % The profile is the column-wise RGB sum normalised to [-1, +1] and
        % scaled to the image half-height, so it spans the full image vertically.
        if showProfile
            stimProfile = squeeze(sum(displayLinearRGBimage(sampleRow, :, :), 3));
            m1 = min(stimProfile);
            m2 = max(stimProfile);
            if (m1 == m2)
                % Uniform row: flat profile at zero (no overlay visible)
                normalizedStimProfile = stimProfile * 0;
            else
                % Map [m1, m2] -> [-1, +1]
                normalizedStimProfile = 2 * ((stimProfile - m1) / (m2 - m1) - 0.5);
            end
            plot(ax, x, yPixels/2 * normalizedStimProfile, 'k.-');
        end

        axis(ax, 'image');
        hold(ax, 'off');
        set(ax, 'CLim', [0 1], 'XLim', 0.5*xPixels*[-1 1], 'YLim', 0.5*yPixels*[-1 1]);
        set(ax, 'FontSize', 16, 'XTick', [], 'YTick', []);
        title(ax, sprintf( ...
            '%s-%s\nmean(LMS): %2.3f,%2.3f,%2.3f, mean(xyY) = (%2.2f,%2.2f,%2.1f cd/m2)\ncenter(LMS): %2.3f,%2.3f,%2.3f, center(xyY) = (%2.2f,%2.2f,%2.1f cd/m2)\n(t:%2.0f msec)', ...
            theSceneEngineName, theSceneName, ...
            meanL, meanM, meanS, ...
            mean(xChroma(:)), mean(yChroma(:)), sceneGet(sceneSequence{frameIndex}, 'mean luminance'), ...
            centerL, centerM, centerS, ...
            centerChromaX, centerChromaY, centerLuminance, ...
            temporalSupportSeconds(frameIndex) * 1000));

        % --- Right panel: RGB-gun modulation trace, grown frame by frame ----
        plot(axModulation, temporalSupportSeconds(1:frameIndex), RGBgunTrace(1:frameIndex,1), ...
            'ro-', 'MarkerSize', 12, 'MarkerFaceColor', [1 0.5 0.5], 'LineWidth', 1.5);
        hold(axModulation, 'on');
        plot(axModulation, temporalSupportSeconds(1:frameIndex), RGBgunTrace(1:frameIndex,2), ...
            'go-', 'MarkerSize', 12, 'MarkerFaceColor', [0 1 0.5], 'LineWidth', 1.5);
        plot(axModulation, temporalSupportSeconds(1:frameIndex), RGBgunTrace(1:frameIndex,3), ...
            'bo-', 'MarkerSize', 12, 'MarkerFaceColor', [0.5 0.5 1], 'LineWidth', 1.5);
        hold(axModulation, 'off');
        set(axModulation, 'FontSize', 16, ...
            'YLim', [min(min(RGBgunTrace(1:frameIndex,:))) max(max(RGBgunTrace(1:frameIndex,:)))], ...
            'XLim', [temporalSupportSeconds(1) temporalSupportSeconds(end) + eps]);
        xlabel(axModulation, 'time (seconds)');
        title(axModulation, sprintf('RGB gun modulation at pixel (%d,%d)', sampleRow, sampleCol));
        drawnow;

        if (~isempty(videoFileName))
            videoOBJ.writeVideo(getframe(hFig));
        end
    end

    if (~isempty(videoFileName))
        videoOBJ.close();
    end

end
