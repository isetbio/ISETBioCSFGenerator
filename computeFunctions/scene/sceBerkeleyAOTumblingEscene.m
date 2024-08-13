function dataOut = sceBerkeleyAOTumblingEscene(sceneEngineOBJ, testESizeDeg, sceneParams)
    % Check input arguments. If called with zero input arguments, just return the default params struct
    if (nargin == 0)
        dataOut = generateDefaultParams();
        return;
    end

    % Validate params
    assert(rem(sceneParams.displayPixelSize,2) == 0, 'display assumed to have even number of pixels');
    assert(rem(sceneParams.letterHeightPixels,2) == 0, 'letterHeight must be an even number of pixels');
    assert(rem(sceneParams.letterWidthPixels,2) == 0, 'letterWidth must be an even number of pixels');
    assert(ismember(sceneParams.letterRotationDegs, [0 90 180 270]), 'letterRotationDegs must be in 0, 90, 180, or 270');

    % Assign passed sceneParams to struct paramsForTextRendering, to avoid
    % confusion.  We will then add fields to the latter and pass it on.
    paramsForTextRendering = sceneParams;

    % Generate the display on which the tumbing E scenes are presented
    presentationDisplay = generateBerkeleyAOPresentationDisplay(sceneParams);

    % Decide bitmap to use.  Here we define a custom E bitmap in a function
    % at the bottom of this script, that returns a 15 by 15 pixel square E.
    theCustomBitMap = textToCustomBitmap('E');

    % Letter pixel aspect ratio. This is the pixel 
    letterPixelAspectRatio = size(theCustomBitMap,2)/size(theCustomBitMap,1);

    % If we define our letter bitmap to be exactly the letter size, then
    % this can be 1.  If the letter has margins around it in the bitmap,
    % this factor can adjust the size in pixels in the code below.
    sizeScalingFactor = 1.0;

    % We can define an extra shift of the letter in pixels relative to its
    % nominal center, to try to get the letter as centered as possible when
    % it is nominally at 0,0. 
    %
    % Allowing this was an idea we had, but because the bitmapped letter
    % gets scaled below, quite what values we want may depend on the letter
    % size.
    paramsForTextRendering.yOffset = 0; % shift up (negative) or down (positive)
    paramsForTextRendering.xOffset = 0; % shift left (negative) or right (positive)

    % Figure out how many rows and columns we want the E bitmap to be
    testESizeMin = testESizeDeg*60;
    letterHeightUnquantizedPixels = sizeScalingFactor*((testESizeMin/60)/(paramsForTextRendering.displayFOVDeg))*paramsForTextRendering.displayPixelSize;
    paramsForTextRendering.letterHeightPixels = 2*round(letterHeightUnquantizedPixels/2);
    paramsForTextRendering.letterWidthPixels = 2*round(letterHeightUnquantizedPixels*(letterPixelAspectRatio)/2);
    paramsForTextRendering.yPixelsNumMargin = (paramsForTextRendering.displayPixelSize-paramsForTextRendering.letterHeightPixels)/2;
    paramsForTextRendering.xPixelsNumMargin =  (paramsForTextRendering.displayPixelSize-paramsForTextRendering.letterWidthPixels)/2;
    paramsForTextRendering.upSampleFactor = uint8(1);

    % Check pixel consistency
    rowsNum = paramsForTextRendering.letterHeightPixels + paramsForTextRendering.yPixelsNumMargin*2;
    colsNum = paramsForTextRendering.letterWidthPixels + paramsForTextRendering.xPixelsNumMargin*2;
    if (rowsNum ~= paramsForTextRendering.displayPixelSize || colsNum ~= paramsForTextRendering.displayPixelSize)
        error('Have not computed rowsNum and colsNum or some other parameter consistently.');
    end

    [theTestSceneSequence, temporalSupportSeconds] = generateTumblingEsceneSequence(...
        presentationDisplay, theCustomBitMap, paramsForTextRendering, ...
        'visualizeScene', paramsForTextRendering.visualizeScene);

    % Assemble dataOut struct - required fields
    dataOut.sceneSequence = theTestSceneSequence;
    dataOut.temporalSupport = temporalSupportSeconds;
    dataOut.presentationDisplay = presentationDisplay;
    dataOut.statusReport = 'done';
end

function [theSceneSequence, temporalSupportSeconds] = generateTumblingEsceneSequence(...
    presentationDisplay, theChar, paramsForTextRendering ...
    , varargin)

    % Parse optional input
    p = inputParser;
    p.addParameter('visualizeScene', false, @islogical);
    p.parse(varargin{:});
    visualizeScene = p.Results.visualizeScene;

    % Extract sequence parameters.  For backwards compatibility, there is
    % just one frame and no shift if nothing is specfied (and in that case,
    % the frame rate has no effect).
    if (~isfield(paramsForTextRendering,'temporalModulationParams') || isempty(paramsForTextRendering.temporalModulationParams))
        frameRateHz = 60;
        numFrames = 1;
        xShiftPerFrame = 0;
        yShiftPerFrame = 0;
    else
        frameRateHz = paramsForTextRendering.temporalModulationParams.frameRateHz;
        numFrames = paramsForTextRendering.temporalModulationParams.numFrames;
        xShiftDegrees = paramsForTextRendering.temporalModulationParams.xShiftPerFrame;
        yShiftDegrees = paramsForTextRendering.temporalModulationParams.yShiftPerFrame;
    end
    frameDurationSec = 1 / frameRateHz;

    % Calculate degrees per pixel
    displayFOVDeg = paramsForTextRendering.displayFOVDeg;
    displayPixelSize = paramsForTextRendering.displayPixelSize;
    degreesPerPixel = displayFOVDeg / displayPixelSize;

    % Convert shifts from degrees to pixels
    xShiftPixels = round(xShiftDegrees / degreesPerPixel);
    yShiftPixels = round(yShiftDegrees / degreesPerPixel);

    % Initialize output
    theSceneSequence = cell(1, numFrames);
    temporalSupportSeconds = zeros(1, numFrames);

    xPixelsNumMargin0 = paramsForTextRendering.xPixelsNumMargin;
    yPixelsNumMargin0 = paramsForTextRendering.yPixelsNumMargin;

    % Generate each frame
    for frameIndex = 1:numFrames
        % Update scene parameters for current frame.  This handles the
        % specified shift of the letter on each frame.
        paramsForTextRendering.xPixelsNumMargin = xPixelsNumMargin0 + xShiftPixels(frameIndex) + paramsForTextRendering.xOffset;
        paramsForTextRendering.yPixelsNumMargin = yPixelsNumMargin0  + yShiftPixels(frameIndex) + paramsForTextRendering.yOffset;

        % Generate the scene frame
        theSceneFrame = generateTumblingEscene(presentationDisplay, theChar, paramsForTextRendering, 'visualizeScene', visualizeScene);
        
        % Store the scene frame
        theSceneSequence{frameIndex} = theSceneFrame;

        % Store the timestamp of the scene frame
        temporalSupportSeconds(frameIndex) = (frameIndex - 1) * frameDurationSec;
    end
end

function theScene = generateTumblingEscene(...
    presentationDisplay, theChar, paramsForTextRendering, varargin)

    % Parse optional input
    p = inputParser;
    p.addParameter('visualizeScene', false, @islogical);
    p.parse(varargin{:});
    visualizeScene = p.Results.visualizeScene;

    % Generate stimulus scene
    textSceneParams = struct(...
        'textString', theChar, ...                                                         % Text to display, can be letter or bitmap. 
        'textRotation', paramsForTextRendering.letterRotationDegs, ...                     % Rotation (0,90,180,270 only)
        'rowsNum', paramsForTextRendering.displayPixelSize, ...                            % Pixels along the vertical (y) dimension
        'colsNum', paramsForTextRendering.displayPixelSize, ...                            % Pixels along the horizontal (x) dimension
        'textBitMapRescaledRowsCols', [paramsForTextRendering.letterHeightPixels paramsForTextRendering.letterWidthPixels], ...
        'targetRow', paramsForTextRendering.yPixelsNumMargin, ...                          % Y-pixel offset 
        'targetCol', paramsForTextRendering.xPixelsNumMargin, ...                          % X-pixel offset 
        'upSampleFactor', uint8(1), ...                                                    % Upsample the scene to increase the retinal image resolution
        'chromaSpecification', paramsForTextRendering.chromaSpecification, ...             % Background and stimulus rgb values
        'temporalModulationParams', paramsForTextRendering.temporalModulationParams, ...   % Parameters describing temporal sequence
        'centerLetter', false ...                                                          % Use rotatedTextSceneRealizedOnDisplay?  (We center manually so false here).
    );

    theScene = rotatedTextSceneRealizedOnDisplay(presentationDisplay, ...
        textSceneParams, visualizeScene);
end

%% Return a structure with the fields and default parameters for this scene engine
function p = generateDefaultParams()

    p = struct(...
        'displayPixelSize', 512, ...            % Linear display pixel size
        'displayFOVDeg', 1.413, ...             % Linear field size in degrees
        'viewingDistanceMeters', 3, ...         % Far enough away to be in good focus
        'letterRotationDegs', 0, ...            % Letter rotation (0,90,180,270 only)
        'letterHeightPixels', 20, ...           % Letter height in pixels - must be 20
        'letterWidthPixels', 18, ...            % Letter width in pixels - must be 18
        'yPixelsNumMargin', 25, ...             % Y-margin
        'xPixelsNumMargin', 25, ...             % X-margin
        'upSampleFactor', uint8(3), ...         % Upsampling for better centering
        'chromaSpecification', struct(...
                'type', 'RGBsettings', ...
                'backgroundRGB', [0.5 0.5 0.5], ...   
                'foregroundRGB',  [0.4 0.4 0.4]), ...
        'temporalModulationParams', struct(...  % temporal: modulation params struct
                'frameRateHz', 60, ...              % frame rate in Hz
                'numFrames', 1, ...                 % number of frames we want the E on for
                'xShiftPerFrame', 0, ...            % shift E in the x dimension in each frame in degree
                'yShiftPerFrame', 0), ...           % shift E in the y dimension in each frame in degree
        'wave', (400:10:860)', ...              % Wavelength sampling for primaries
        'AOPrimaryWls', [840 650 540], ...      % Display spd center wavelengths
        'AOPrimaryFWHM', [10 10 10], ...        % Display spd FWHM in nm
        'AOAOCornealPowersUW', [141.4 0 0], ... % Display spd power full on
        'ambientSpd', [], ...                   % Display ambientSpd
        'pupilSizeMM', 6, ...                   % Need this here to do appropriate primary power conversion to radiance
        'visualizeScene', false, ...            % Whether to visualize the generated scene
        'plotDisplayCharacteristics', false ... % Whether to visualize the display characteristics
       );
end

%% Generate a display that mimics the Berkeley AO system
function presentationDisplay = generateBerkeleyAOPresentationDisplay(sceneParams)

   % We know the dispaly FOV and number of pixels.  We choose a distance
   % far enough way and compute DPI to make it work out.
   inchesPerMeter = 39.3701;
   displaySizeMeters = 2*sceneParams.viewingDistanceMeters*tand(sceneParams.displayFOVDeg/2);
   displayDotsPerMeter = sceneParams.displayPixelSize/displaySizeMeters;
   displayDPI = displayDotsPerMeter/inchesPerMeter;

   % Now build the display given all the wonderful things we know about it.
   presentationDisplay = generateCustomDisplay(...
       'viewingDistanceMeters', sceneParams.viewingDistanceMeters, ...
       'dotsPerInch', displayDPI, ...
       'wavelengthSupportNanoMeters', sceneParams.wave, ...
       'spectralPowerDistributionWattsPerSteradianM2NanoMeter', sceneParams.spd, ...
       'ambientSPDWattsPerSteradianM2NanoMeter', sceneParams.ambientSpd, ...
       'gammaTable', repmat((linspace(0,1,1024)').^2, [1 3]), ...
       'plotCharacteristics', sceneParams.plotDisplayCharacteristics);
    
end

function TxtIm = textToCustomBitmap(text)

switch (text)
    case 'E'
        % 20 x 18 E, with actual E 15 x 15
        % TxtIm=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
        %     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
        %     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
        %     1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
        %     1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
        %     1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
        %     1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
        %     1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
        %     1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
        %     1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
        %     1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
        %     1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
        %     1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
        %     1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
        %     1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
        %     1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
        %     1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
        %     1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
        %     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
        %     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];

        % 15 by 15 E
        TxtIm=[
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
            0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
            0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
            0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
            0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
            0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
            0,0,0,1,1,1,1,1,1,1,1,1,1,1,1;
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
            ];
    otherwise
        error('Do not know about passed text');
end
end