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

    % Generate the display on which the tumbing E scenes are presented
    presentationDisplay = generateBerkeleyAOPresentationDisplay(sceneParams);

    % The height of the actual E in the bitmap is only 15 rows, not the
    % nominal 20.  We care about the E size, not the bitmap size.  So we
    % bump up the desired number of rows by this factor in the computation
    % just below.
    sizeScalingFactor = 1.0/0.75;

    % Figure out how many rows and columns we want the E bitmap to be
    testESizeMin = testESizeDeg*60;
    letterHeightUnquantizedPixels = sizeScalingFactor*((testESizeMin/60)/(sceneParams.displayFOVDeg))*sceneParams.displayPixelSize;
    sceneParams.letterHeightPixels = 2*round(letterHeightUnquantizedPixels/2);
    sceneParams.letterWidthPixels = 2*round(letterHeightUnquantizedPixels*(18/20)/2);
    sceneParams.yPixelsNumMargin = (sceneParams.displayPixelSize-sceneParams.letterHeightPixels)/2;
    sceneParams.xPixelsNumMargin =  (sceneParams.displayPixelSize-sceneParams.letterWidthPixels)/2;
    sceneParams.upSampleFactor = uint8(1);

    %% CONVERT SHIFT TO FROM DEGS TO PIXELS.  SPECIFY AN OFFSET IN PIXELS IN ADDITION
    % TO SCALE FACTOR SO THAT THE E ENDS UP CLOSER TO CENTERED.

    % Check pixel consistency
    rowsNum = sceneParams.letterHeightPixels + sceneParams.yPixelsNumMargin*2;
    colsNum = sceneParams.letterWidthPixels + sceneParams.xPixelsNumMargin*2;
    if (rowsNum ~= sceneParams.displayPixelSize || colsNum ~= sceneParams.displayPixelSize)
        error('Have not computed rowsNum and colsNum or some other parameter consistently.');
    end

    % Generate the E scene frame sequence 
    [theTestSceneSequence, temporalSupportSeconds] = generateTumblingEsceneSequence(...
        presentationDisplay, 'E', sceneParams, ...
        'visualizeScene', sceneParams.visualizeScene);

    % Assemble dataOut struct - required fields
    dataOut.sceneSequence = theTestSceneSequence;
    dataOut.temporalSupport = temporalSupportSeconds;
    dataOut.presentationDisplay = presentationDisplay;
    dataOut.statusReport = 'done';
end

function [theSceneSequence, temporalSupportSeconds] = generateTumblingEsceneSequence(...
    presentationDisplay, theChar, sceneParams, varargin)

    % Parse optional input
    p = inputParser;
    p.addParameter('visualizeScene', false, @islogical);
    p.parse(varargin{:});
    visualizeScene = p.Results.visualizeScene;

    % Extract sequence parameters.  For backwards compatibility, there is
    % just one frame and no shift if nothing is specfied (and in that case,
    % the frame rate has no effect).
    if (~isfield(sceneParams,'temporalModulationParams') || isempty(sceneParams.temporalModulationParams))
        frameRateHz = 60;
        numFrames = 1;
        xShiftPerFrame = 0;
        yShiftPerFrame = 0;
    else
        frameRateHz = sceneParams.temporalModulationParams.frameRateHz;
        numFrames = sceneParams.temporalModulationParams.numFrames;
        xShiftDegrees = sceneParams.temporalModulationParams.xShiftPerFrame;
        yShiftDegrees = sceneParams.temporalModulationParams.yShiftPerFrame;
    end
    frameDurationSec = 1 / frameRateHz;

    % Calculate degrees per pixel
    displayFOVDeg = sceneParams.displayFOVDeg;
    displayPixelSize = sceneParams.displayPixelSize;
    degreesPerPixel = displayFOVDeg / displayPixelSize;

    % Convert shifts from degrees to pixels
    xShiftPixels = xShiftDegrees / degreesPerPixel;
    yShiftPixels = yShiftDegrees / degreesPerPixel;

    % Round the shift to the nearest integer of pixels
    xShiftPixels = round(xShiftPixels);
    yShiftPixels = round(yShiftPixels);

    % Initialize output
    theSceneSequence = cell(1, numFrames);
    temporalSupportSeconds = zeros(1, numFrames);

    xPixelsNumMargin0 = sceneParams.xPixelsNumMargin;
    yPixelsNumMargin0 = sceneParams.yPixelsNumMargin;

    % shift extra pixels by default to let E be on the center
    yOffset = 0; % shift up -0.5
    xOffset = 1; % shift to right

    % Generate each frame
    for frameIndex = 1:numFrames
        % Update scene parameters for current frame.  This handles the
        % specified shift of the letter on each frame.
        sceneParams.xPixelsNumMargin = xPixelsNumMargin0 + (frameIndex - 1) * xShiftPixels(frameIndex) + xOffset;
        sceneParams.yPixelsNumMargin = yPixelsNumMargin0  + (frameIndex - 1) * yShiftPixels(frameIndex) + yOffset;

        % Generate the scene frame
        theSceneFrame = generateTumblingEscene(presentationDisplay, theChar, sceneParams, 'visualizeScene', visualizeScene);
        
        % Store the scene frame
        theSceneSequence{frameIndex} = theSceneFrame;

        % Store the timestamp of the scene frame
        temporalSupportSeconds(frameIndex) = (frameIndex - 1) * frameDurationSec;
    end
end

function theScene = generateTumblingEscene(...
    presentationDisplay, theChar, sceneParams, varargin)

    % Parse optional input
    p = inputParser;
    p.addParameter('visualizeScene', false, @islogical);
    p.parse(varargin{:});
    visualizeScene = p.Results.visualizeScene;

    % Generate stimulus scene
    textSceneParams = struct(...
        'textString', theChar, ...                                              % Text to display
        'textRotation', sceneParams.letterRotationDegs, ...                     % Rotation (0,90,180,270 only)
        'rowsNum', sceneParams.displayPixelSize, ...                            % Pixels along the vertical (y) dimension
        'colsNum', sceneParams.displayPixelSize, ...                            % Pixels along the horizontal (x) dimension
        'textBitMapRescaledRowsCols', [sceneParams.letterHeightPixels sceneParams.letterWidthPixels], ...
        'targetRow', sceneParams.yPixelsNumMargin, ...                          % Y-pixel offset 
        'targetCol', sceneParams.xPixelsNumMargin, ...                          % X-pixel offset 
        'upSampleFactor', uint8(1), ...                                         % Upsample the scene to increase the retinal image resolution
        'chromaSpecification', sceneParams.chromaSpecification, ...             % Background and stimulus rgb values
        'temporalModulationParams', sceneParams.temporalModulationParams, ...   % Parameters describing temporal sequence
        'centerLetter', false ...                                               % Use rotatedTextSceneRealizedOnDisplay?  (We center manually so false here).
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