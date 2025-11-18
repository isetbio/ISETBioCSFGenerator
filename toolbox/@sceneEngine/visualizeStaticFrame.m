function visualizeStaticFrame(obj, sceneSequence, varargin)

% Parse input
p = inputParser;
p.addParameter('skipOutOfGamutCheck', false, @islogical);
p.addParameter('frameToVisualize', 1, @isscalar);
p.addParameter('opticalImageInsteadOfScene', []);
p.addParameter('sRGBforSceneVisualization', false, @islogical);
p.addParameter('axesHandle', [], @(x)(isempty(x)||(ishandle(x))));

p.parse(varargin{:});
skipOutOfGamutCheck = p.Results.skipOutOfGamutCheck;
frameIndex = p.Results.frameToVisualize;
theOI = p.Results.opticalImageInsteadOfScene;
axesHandle = p.Results.axesHandle;
sRGBforSceneVisualization = p.Results.sRGBforSceneVisualization;

% Get the scene at the frame index
theScene = sceneSequence{frameIndex};

% Compute the optical image of the scene if 'opticalImageInteadOfScene' is set
if (~isempty(theOI))
    theOI = oiCompute(theOI, theScene);
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
displayLinearRGBToLMS = displayGet(presentationDisplay, 'rgb2lms');

% Extract the XYZ image representation
if (isempty(theOI))
    xyzImage = sceneGet(theScene, 'xyz');
else
    xyzImage = oiGet(theOI, 'xyz');
end


% Linear RGB image
displayLinearRGBimage = imageLinearTransform(xyzImage, displayXYZToLinearRGB);
displayLinearLMSimage = imageLinearTransform(displayLinearRGBimage,displayLinearRGBToLMS);
Limage = squeeze(displayLinearLMSimage(:,:,1));
Mimage = squeeze(displayLinearLMSimage(:,:,2));
Simage = squeeze(displayLinearLMSimage(:,:,3));
meanL = mean(Limage(:));
meanM = mean(Mimage(:));
meanS = mean(Simage(:));

xChroma = xyzImage(:,:,1)./sum(xyzImage,3);
yChroma = xyzImage(:,:,2)./sum(xyzImage,3);
luma = xyzImage(:,:,3);

middleRow = round(size(Limage,1)/2);
middleCol = round(size(Limage,2)/2);
centerL = Limage(middleRow, middleCol);
centerM = Mimage(middleRow, middleCol);
centerS = Simage(middleRow, middleCol);
centerChromaX = xChroma(middleRow, middleCol,1);
centerChromaY = yChroma(middleRow, middleCol,1);
centerLuminance = luma(middleRow, middleCol,1);


if (frameIndex == 1)
    xPixels = size(xyzImage,2);
    yPixels = size(xyzImage,1);
    x = 1:xPixels;
    y = 1:yPixels;
    x = x-mean(x);
    y = y-mean(y);
end

if (~skipOutOfGamutCheck)
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

end

if (isempty(theOI))
    displaySettingsImage = sceneGet(theScene, 'rgb');
else
    displaySettingsImage = oiGet(theOI, 'rgb');
end

% Render image
if (isempty(axesHandle))
    axesHandle = gca;
end

if (sRGBforSceneVisualization)
    image(axesHandle,x,y,lrgb2srgb(displaySettingsImage));
else
    image(axesHandle,x,y,displaySettingsImage);
end

axis(axesHandle, 'image');

% Cross hairs
hold(axesHandle,'on');
set(axesHandle, 'XTick', [], 'YTick', []);

title(axesHandle,sprintf('mean LMS: %2.3f,%2.3f,%2.3f\nmean xyY: (%2.2f,%2.2f,%2.1f(%2.1f) cd/m2)\ncenter LMS: %2.3f,%2.3f,%2.3f\ncenter xyY: (%2.2f,%2.2f,%2.1f cd/m2)', ...
    meanL, meanM, meanS, ...
    mean(xChroma(:)), mean(yChroma(:)), sceneGet(sceneSequence{frameIndex}, 'mean luminance'), mean(luma(:)), ...
    centerL, centerM, centerS, ...
    centerChromaX, centerChromaY, centerLuminance));


drawnow;

end