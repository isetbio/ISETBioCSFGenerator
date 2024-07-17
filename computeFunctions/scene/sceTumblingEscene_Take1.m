% The tublingEsceneEngine.compute(Esize) compute function
function dataOut = sceTumblingEscene(sceneEngineOBJ, testEsizeDegs, sceneParams)
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
    if (~isfield(sceneParams,'presentationDisplay') | isempty(sceneParams.presentationDisplay))
        presentationDisplay = generatePresentationDisplay(...
            testEsizeDegs*sizeScalingFactor, sceneParams.letterHeightPixels, ...
            sceneParams.spdDataFile, ...
            sceneParams.ambientSPDDataFile, ...
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
        'letterRotationDegs', 0, ...            % Letter rotation (0,90,180,270 only)
        'letterHeightPixels', 20, ...           % Letter height in pixels - must be 20
        'letterWidthPixels', 18, ...            % Letter width in pixels - must be 18
        'yPixelsNumMargin', 25, ...             % Y-margin
        'xPixelsNumMargin', 25, ...             % X-margin
        'presentationDisplay',[], ...           % If this not empty, it is the display file to use
        'upSampleFactor', uint8(3), ...         % Upsampling for better centering
        'chromaSpecification', struct(...
                'type', 'RGBsettings', ...
                'backgroundRGB', [0.5 0.5 0.5], ...   
                'foregroundRGB',  [0.4 0.4 0.4]), ...
        'spdDataFile', 'BVAMS_White_Guns_At_Max.mat', ...  % Name of file containing the display SPDs
        'ambientSPDDataFile', 'BVAMS_White_Background.mat', ... % Name of the file containing the ambient SPD
        'visualizeScene', false, ...                % Whether to visualize the generated scene
        'plotDisplayCharacteristics', false ...     % Whether to visualize the display characteristics
       );
end

% Create an ISETBio display given parameters we are interested in.  This
% was originally for the JandJISETBio chromatic aberration project
function presentationDisplay = generatePresentationDisplay(...
    letterSizeDegs, letterSizePixels, spdDataFile, ambientSPDDataFile, plotCharacteristics)

    % Load the ambient SPD
    projectBaseDir = ISETBioCSFGeneratorRootPath();

    fprintf('Loading ambient SPD from %s\n', fullfile(projectBaseDir,'sampledata',ambientSPDDataFile));
    load(fullfile(projectBaseDir,'sampledata',ambientSPDDataFile),'spd');
    ambientSPD = spd;
    clear 'spd'
    
    % Check data consistency
    assert(size(ambientSPD, 2) == 2, 'The ambient SPD matrix must be an N x 2 matrix, with the first column being the spectral support');
    if (size(ambientSPD,2) > 2)
        fprintf(2,'\nThe ambient SPD matrix must be an N x 2 matrix, with the first column being the spectral support and the second column being the ambient energy.\n');
        fprintf(2,'The data retrieved from ''%s'', contain %d columns. Ignoring all but the first 2 columns.\n\n', ambientSPDDataFile, size(ambientSPD,2));
    end
    
    ambientWave = ambientSPD(:,1);
    ambientSPD = ambientSPD(:,2);
    ambientSPD = ambientSPD/ (ambientWave(2)-ambientWave(1));
    
    % Load the RGB SPDs
    fprintf('Loading SPDs from %s\n', fullfile(getpref('ISETBioJandJ','dataDir'),spdDataFile));
    load(fullfile(getpref('ISETBioJandJ','dataDir'),spdDataFile), 'spd');
    
    % Check data consistency
    assert(size(spd, 2) == 4, 'The SPD matrix must be an N x 4 matrix, with the first column being the spectral support');
    wave = spd(:,1);
    spd = spd(:,2:4);
    spd = spd / (wave(2)-wave(1));
    assert(size(spd,1) == size(ambientSPD,1), 'The ambient SPD must have the same wavelength entries as the display SPD');
    assert(all(ambientWave == wave), 'The ambient wavelength support must match the wavelength support of the display SPD');
    
    presentationDisplay = generateCustomDisplay(...
           'dotsPerInch', 220, ...
           'wavelengthSupportNanoMeters', wave, ...
           'spectralPowerDistributionWattsPerSteradianM2NanoMeter', spd, ...
           'ambientSPDWattsPerSteradianM2NanoMeter', ambientSPD, ...
           'gammaTable', repmat((linspace(0,1,1024)').^2, [1 3]), ...
           'plotCharacteristics', plotCharacteristics);
    
    pixelSizeMeters = displayGet(presentationDisplay, 'meters per dot');
    letterSizeMeters = letterSizePixels*pixelSizeMeters;
    desiredViewingDistance = 0.5*letterSizeMeters/(tand(letterSizeDegs/2));
    presentationDisplay = displaySet(presentationDisplay, 'viewing distance', desiredViewingDistance);
    
end