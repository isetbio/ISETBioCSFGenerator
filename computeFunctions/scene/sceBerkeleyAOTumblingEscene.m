function dataOut = sceBerkeleyAOTumblingEscene(sceneEngineOBJ, testESizeDeg, sceneParams)
% Compute function for generating a sequence of scenes depicting a
% scene sequence of tumbling E's.
%
% Syntax:
%   dataOut = sceBerkeleyAOTumblingEscene(sceneEngineOBJ, testESizeDeg, sceneParams)
%
% Description:
%    Compute function to be used as a computeFunctionHandle for a @sceneEngine
%    object. There are 2 ways to use this function.
%
%       [1] If called wihtout arguments, itt does not compute anything and
%           it does not compute anything and simply returns a struct with the 
%           defaultParams that define the scene.
%
%       [2] If called with arguments, as it is from a parent @sceneEngine object,
%           it computes a cell array of scenes defining the frames of a
%           stimulus and the temporal support of the frames. These are
%           returned as named fields of the returned dataOut struct.
%
%    All scene functions used with the sceneEngine class must conform to
%    this API.
%
%    In addition to computing, this function checks the
%    `visualizeEachCompute` flag of the sceneEngineOBJ and, if it is set,
%    calls its visualizeSceneSequence() method of the sceneEngineOBJ. This
%    causes figures to appear that visualize the generated scene sequence,
%    which is helpful for debugging. Note that everything runs much more
%    slowly in this case.  You can also call the visualizeSceneSequence
%    method directly if you want to see the scene sequence at some
%    particular point in your code.  See computeThreshold, for example.
%
% Inputs:
%    sceneEngineOBJ              - Calling @sceneEngine object.  This is
%                                  currently unused, but passing it allows us
%                                  flexibility in the future and matches
%                                  conventions for the other classes in
%                                  this toolbox.
%    testESizeDegs               - Scalar providing the size of the test E in degrees.                          
%    sceneParamsStruct           - Struct containing properties of the
%                                  scene understood by this function.
%                                  As noted above, execute sceGrating at
%                                  the command line to see the structure's
%                                  fields and default values.
% 
%                                  The fields and their meanings are
%                                  provided in the function that generates
%                                  the default parameters at the end of
%                                  this source file.
%
% Outputs:
%    dataOut  - A struct that depends on the input arguments. 
%
%               If called directly with no input arguments, the returned struct contains
%               the defaultParams that define the grating scene
%
%             - If called from a parent @sceneEngine, the returned
%               struct is organized as follows:
%                 .sceneSequence : a cell array of scenes defining the frames of the generated grating scene sequence                            
%                 .temporalSupport : the temporal support of the frames of the generated grating scene sequence, in seconds            
%
% Optional key/value input arguments:
%     None.
%
% See Also:
%     t_sceneGeneration, t_spatialCSF, createGratingSceneEngine.

% History:
%  01/12/24  dhb  Added comments and started work on generalizing

    % Check input arguments. If called with zero input arguments, just return the default params struct
    if (nargin == 0)
        dataOut = generateDefaultParams();
        return;
    end

    % Validate params
    assert(rem(sceneParams.displayParams.displayPixelSize,2) == 0, 'display assumed to have even number of pixels');
    assert(rem(sceneParams.letterHeightPixels,2) == 0, 'letterHeight must be an even number of pixels');
    assert(rem(sceneParams.letterWidthPixels,2) == 0, 'letterWidth must be an even number of pixels');
    assert(ismember(sceneParams.letterRotationDegs, [0 90 180 270]), 'letterRotationDegs must be in 0, 90, 180, or 270');

    % Generate the display on which the tumbing E scenes are presented
    sceneParams.displayParams.spectralSupport = sceneParams.spectralSupport;
    presentationDisplay = generateBerkeleyAODisplay(sceneParams.displayParams);

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
    sceneParams.yOffset = 0; % shift up (negative) or down (positive)
    sceneParams.xOffset = 0; % shift left (negative) or right (positive)

    % Figure out how many rows and columns we want the E bitmap to be
    testESizeMin = testESizeDeg*60;
    letterHeightUnquantizedPixels = sizeScalingFactor*((testESizeMin/60)/(sceneParams.displayParams.displayFOVDeg))*sceneParams.displayParams.displayPixelSize;
    sceneParams.letterHeightPixels = 2*round(letterHeightUnquantizedPixels/2);
    sceneParams.letterWidthPixels = 2*round(letterHeightUnquantizedPixels*(letterPixelAspectRatio)/2);
    sceneParams.yPixelsNumMargin = (sceneParams.displayParams.displayPixelSize-sceneParams.letterHeightPixels)/2;
    sceneParams.xPixelsNumMargin =  (sceneParams.displayParams.displayPixelSize-sceneParams.letterWidthPixels)/2;
    sceneParams.upSampleFactor = uint8(1);

    % Check pixel consistency
    rowsNum = sceneParams.letterHeightPixels + sceneParams.yPixelsNumMargin*2;
    colsNum = sceneParams.letterWidthPixels + sceneParams.xPixelsNumMargin*2;
    if (rowsNum ~= sceneParams.displayParams.displayPixelSize || colsNum ~= sceneParams.displayParams.displayPixelSize)
        error('Have not computed rowsNum and colsNum or some other parameter consistently.');
    end

    [theTestSceneSequence, temporalSupportSeconds] = generateTumblingEsceneSequence(...
        presentationDisplay, theCustomBitMap, sceneParams, ...
        'visualizeScene', sceneParams.visualizeScene);

    % Check the visualizeEachCompute flag of the sceneEngineOBJ, and if set to true,
    % call its visualizeSceneSequence() method to visualize the generate scene sequence.
    if (sceneEngineOBJ.visualizeEachCompute)
        sceneEngineOBJ.visualizeSceneSequence(theSceneSequence, temporalSupportSeconds);
    end

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
        stimOnFrames = 1;
    else
        frameRateHz = paramsForTextRendering.temporalModulationParams.frameRateHz;
        numFrames = paramsForTextRendering.temporalModulationParams.numFrames;
        xShiftDegrees = paramsForTextRendering.temporalModulationParams.xShiftPerFrame;
        yShiftDegrees = paramsForTextRendering.temporalModulationParams.yShiftPerFrame;
        backgroundRGBPerFrame = paramsForTextRendering.temporalModulationParams.backgroundRGBPerFrame;
        stimOnFrames = paramsForTextRendering.temporalModulationParams.stimOnFrames;
    end
    frameDurationSec = 1 / frameRateHz;

    % Calculate degrees per pixel
    displayFOVDeg = paramsForTextRendering.displayParams.displayFOVDeg;
    displayPixelSize = paramsForTextRendering.displayParams.displayPixelSize;
    degreesPerPixel = displayFOVDeg / displayPixelSize;

    % Convert shifts from degrees to pixels
    xShiftPixels = round(xShiftDegrees / degreesPerPixel);
    yShiftPixels = round(yShiftDegrees / degreesPerPixel);

    % Initialize output
    theSceneSequence = cell(1, numFrames);
    temporalSupportSeconds = zeros(1, numFrames);

    xPixelsNumMargin0 = paramsForTextRendering.xPixelsNumMargin;
    yPixelsNumMargin0 = paramsForTextRendering.yPixelsNumMargin;

    % Save actual foregroundRGB
    foregroundRGB = paramsForTextRendering.chromaSpecification.foregroundRGB;

    % Generate each frame
    for frameIndex = 1:numFrames
        % Update scene parameters for current frame.  This handles the
        % specified shift of the letter on each frame.
        paramsForTextRendering.xPixelsNumMargin = xPixelsNumMargin0 + xShiftPixels(frameIndex) + paramsForTextRendering.xOffset;
        paramsForTextRendering.yPixelsNumMargin = yPixelsNumMargin0  + yShiftPixels(frameIndex) + paramsForTextRendering.yOffset;

        % change background for each frame
        paramsForTextRendering.chromaSpecification.backgroundRGB = backgroundRGBPerFrame(frameIndex, :);
        if (stimOnFrames(frameIndex))
            paramsForTextRendering.chromaSpecification.foregroundRGB = foregroundRGB;
        else
            paramsForTextRendering.chromaSpecification.foregroundRGB = backgroundRGBPerFrame(frameIndex, :);
        end

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
        'rowsNum', paramsForTextRendering.displayParams.displayPixelSize, ...                            % Pixels along the vertical (y) dimension
        'colsNum', paramsForTextRendering.displayParams.displayPixelSize, ...                            % Pixels along the horizontal (x) dimension
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

    displayParams = generateBerkeleyAODisplayDefaultParams;

    p = struct(...
        'displayParams', displayParams, ...     % Display parameters
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
                'yShiftPerFrame', 0, ...
                'backgroundRGBPerFrame', [1 0 0]), ...           % shift E in the y dimension in each frame in degree
        'spectralSupport', (400:10:860)', ...   % Wavelength sampling for primaries
        'AOPrimaryWls', [840 650 540], ...      % Display spd center wavelengths
        'AOPrimaryFWHM', [10 10 10], ...        % Display spd FWHM in nm
        'AOAOCornealPowersUW', [141.4 0 0], ... % Display spd power full on
        'ambientSpd', [], ...                   % Display ambientSpd
        'pupilSizeMM', 6, ...                   % Need this here to do appropriate primary power conversion to radiance
        'visualizeScene', false, ...            % Whether to visualize the generated scene
        'plotDisplayCharacteristics', false ... % Whether to visualize the display characteristics
       );
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