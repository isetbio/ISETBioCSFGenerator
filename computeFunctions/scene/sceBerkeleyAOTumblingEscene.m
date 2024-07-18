% The tublingEsceneEngine.compute(Esize) compute function
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
    %
    % WE ARE NO LONGER GOING TO READ A DISPLAY BUT INSTEAD WILL GENERATE IT
    % HERE BASED ON THE BERKELEY PARAMETERS.
    %
    sceneParams.wave
    sceneParams.spd
    sceneParams.ambiendSpd
    if (~isfield(sceneParams,'presentationDisplay') | isempty(sceneParams.presentationDisplay))
        presentationDisplay = generateBerkeleyAOPresentationDisplay(...
            testEsizeDegs*sizeScalingFactor, sceneParams.letterHeightPixels, ...
            sceneParams.wave, ...
            sceneParams.spd, ...
            sceneParams.ambientSpd, ...
            sceneParams.plotDisplayCharacteristics);
    else
        presentationDisplay  = sceneParams.presentationDisplay;
        presentationDisplay = displaySet(presentationDisplay,'dpi',220);
        pixelSizeMeters = displayGet(presentationDisplay, 'meters per dot');
        letterSizeMeters = sceneParams.letterHeightPixels*pixelSizeMeters;
        desiredViewingDistance = 0.5*letterSizeMeters/(tand(testEsizeDegs*sizeScalingFactor/2));
        presentationDisplay = displaySet(presentationDisplay, 'viewing distance', desiredViewingDistance);    
    end

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

    textSceneParams = struct(...
        'textString', theChar, ...                              % Text to display
        'textRotation', sceneParams.letterRotationDegs, ...                     % Rotation (0,90,180,270 only)
        'rowsNum', sceneParams.letterHeightPixels + sceneParams.yPixelsNumMargin*2, ... % Pixels along the vertical (y) dimension
        'colsNum', sceneParams.letterWidthPixels + sceneParams.xPixelsNumMargin*2, ...  % Pixels along the horizontal (x) dimension
        'targetRow', sceneParams.yPixelsNumMargin, ...                      % Y-pixel offset 
        'targetCol', sceneParams.xPixelsNumMargin, ...                      % X-pixel offset 
        'upSampleFactor', sceneParams.upSampleFactor, ...                         % Upsample the scene to increase the retinal image resolution
        'chromaSpecification', sceneParams.chromaSpecification ...          % Background and stimulus chromaticity
    );
    
    % Generate stimulus scene
    theScene = rotatedTextSceneRealizedOnDisplay(presentationDisplay, ...
        textSceneParams, visualizeScene);
end


function p = generateDefaultParams()

    p = struct(...
        'displayPixelSize', 512, ...            % Linear display pixel size
        'displayFOVDeg', 1.413, ...             % Linear field size in degrees
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
        'wave', (400:10:860)',                  % Wavelength sampling for primaries
        'AOPrimaryWls', [840 650 540], ...      % Display spd center wavelengths
        'AOPrimaryFWHM', [10 10 10], ...        % Display spd FWHM in nm
        'AOAOCornealPowersUW', [141.4 0 0], ... % Display spd power full on
        'ambientSpd', [], ...                   % Display ambientSpd
        'visualizeScene', false, ...            % Whether to visualize the generated scene
        'plotDisplayCharacteristics', false ... % Whether to visualize the display characteristics
       );
end

%% WE NEED TO UPDATE THIS TO PRODUCE A DISPLAY THAT MATCHES THE BERKELEY SYSTEM.
%
% spdDataFile -> spd: Instead of reading a date file which has an nWls by 3
% matrix of spectral power distributions, we will pass such a matrix along
% with the wavelengths.  We'll make the spectra in the first column of the matrix match up with
% the imaging channel light specified above.  The other two columns can be
% anything since for this experiment there is just one wavelength involved.
% So, for example, all zeros should work OK.
%
% See lines 180 to 194 of sceAOTumblingESceneFromTwoSpot for how to make
% the spd for the first column.
%
% ambientSPDDataFile -> ambientSpd instead of a file, we'll pass a one column vector
% with nWls samples.  This describes the background light in the AO system.
% So I guess create a spectrum at 940 nm since that is what Will says is
% the other light source.  Some logic as making the 840 spectrum.

function presentationDisplay = generateBerkeleyAOPresentationDisplay(sceneParams)

   % First thing we need to do is generate the primary spectra.  We know
   % the center wavelength, FWHM, and power for each.
s

   % We need to figure out a distance and pixel meters per do that will
   % make the display have the desired visual angle on the retina
   pixelSizeMeters = displayGet(presentationDisplay, 'meters per dot');
   letterSizeMeters = letterSizePixels*pixelSizeMeters;
   desiredViewingDistance = 0.5*letterSizeMeters/(tand(letterSizeDegs/2));
   presentationDisplay = displaySet(presentationDisplay, 'viewing distance', desiredViewingDistance);

   % Now build the display given all the wonderful things we know about it.
   presentationDisplay = generateCustomDisplay(...
       'dotsPerInch', 220, ...
       'wavelengthSupportNanoMeters', wave, ...
       'spectralPowerDistributionWattsPerSteradianM2NanoMeter', spd, ...
       'ambientSPDWattsPerSteradianM2NanoMeter', ambientspd, ...
       'gammaTable', repmat((linspace(0,1,1024)').^2, [1 3]), ...
       'plotCharacteristics', plotCharacteristics);


    
end