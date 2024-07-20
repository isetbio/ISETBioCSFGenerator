% The tumblingEsceneEngine.compute(Esize) compute function
function dataOut = sceBerkeleyAOTumblingEscene(sceneEngineOBJ, testEsizeDegs, sceneParams)
    % Check input arguments. If called with zero input arguments, just return the default params struct
    if (nargin == 0)
        dataOut = generateDefaultParams();
        return;
    end

    % Validate params
    assert(sceneParams.letterHeightPixels == 20, 'letterHeight must be 20 pixels');
    assert(sceneParams.letterWidthPixels == 18, 'letterWidth must be 18 pixels');
    assert(ismember(sceneParams.letterRotationDegs, [0 90 180 270]), 'letterRotationDegs must be in 0, 90, 180, or 270');

    % Account for the fact that the E does not occupy all 20 pixels along the
    % y-dimension
    sizeScalingFactor = 1.0/0.75;

    % Generate the display on which the tumbing E scenes are presented
    presentationDisplay = generateBerkeleyAOPresentationDisplay(sceneParams);

    % Generate the E scene with 0 deg rotation
    theTestScene = generateTumblingEscene(...
        presentationDisplay, 'E', sceneParams, ...
        'visualizeScene', sceneParams.visualizeScene);

    % Assemble dataOut struct - required fields
    dataOut.sceneSequence{1} = theTestScene;
    dataOut.temporalSupport(1) = 0;
    dataOut.presentationDisplay = presentationDisplay;
    dataOut.statusReport = 'done';
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
        'rowsNum', sceneParams.letterHeightPixels + sceneParams.yPixelsNumMargin*2, ... % Pixels along the vertical (y) dimension
        'colsNum', sceneParams.letterWidthPixels + sceneParams.xPixelsNumMargin*2, ...  % Pixels along the horizontal (x) dimension
        'targetRow', sceneParams.yPixelsNumMargin, ...                          % Y-pixel offset 
        'targetCol', sceneParams.xPixelsNumMargin, ...                          % X-pixel offset 
        'upSampleFactor', uint8(1), ...                       % Upsample the scene to increase the retinal image resolution
        'horizontalFOVDegs', sceneParams.eWidthMin/60, ...
        'chromaSpecification', sceneParams.chromaSpecification ...              % Background and stimulus rgb values
    );
    theScene = rotatedTextSceneRealizedOnDisplay(presentationDisplay, ...
        textSceneParams, visualizeScene);
end

%% Return a structure with the fields and default parameters for this scene engine
function p = generateDefaultParams()

    p = struct(...
        'displayPixelSize', 512, ...            % Linear display pixel size
        'displayFOVDeg', 1.413, ...             % Linear field size in degrees
        'viewingDistanceMeters', 3, ...            % Far enough away to be in good focus
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