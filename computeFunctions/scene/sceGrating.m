function dataOut = sceGrating(testContrast, gratingParams)
% Compute function for generating a sequence of scenes depicting a
% a grating that is either drifting, counter-phase modulated, or flashed.
%
% Syntax:
%   dataOut = sceGrating(testContrast, gratingParams);
%
% Description:
%    Compute function to be used as a computeFunctionHandle for a @sceneEngine
%    object. There are 2 ways to use this function.
%
%       [1] If called directly and with no arguments, 
%
%               dataOut = sceGrating()
%
%           it does not compute anything and simply returns a struct with the 
%           defaultParams that define the scene.
%
%       [2] If called from a parent @sceneEngine object, it computes a cell array 
%           of scenes defining the frames of a grating stimulus and the temporal
%           support of the frames
%
% Inputs:
%    testContrast                - the contrast for the scene to be generated
%                               
%    sceneParamsStruct           - a struct containing spatio-temporo-chromatic 
%                                  properties of the grating
%
%
% Optional key/value input arguments: none 
%
% Outputs:
%    dataOut  - A struct that depends on the input arguments. 
%
%               If called directly with no input arguments, the returned struct contains
%               the defaultParams that define the grating scene
%
%             - If called from a parent @sceneEngine, the returned
%               struct is organized as follows:
%
%              .sceneSequence : a cell array of scenes defining the frames of the generated grating scene sequence
%                               
%              .temporalSupport : the temporal support of the frames of the generated grating scene sequence, in seconds            
%
%
%
% See Also:
%     t_sceneGeneration

% History:
%    10/05/2020  NPC  Wrote it.
%
%   Examples:
%{
    % Usage case #1. Just return the default scene params
    defaultParams = sceGrating()

    % Usage case #2. Compute a grating sequence using a parent @sceneEngine 
    % object and the default scene params

    % Instantiate the parent @sceneEngine object
    theSceneEngineOBJ = sceneEngine(@sceGrating);

    % Generate a test scene sequence
    testContrast = 1.0;
    [theTestSceneSequence, temporalSupportSeconds] = ...
        theSceneEngineOBJ.compute(testContrast);

    % Visualize the generated scene frames
    theSceneEngineOBJ.visualizeSceneSequence(theTestSceneSequence, temporalSupportSeconds);
    
%}

    % Check input arguments. If called with zero input arguments, just return the default params struct
    if (nargin == 0)
        dataOut = generateDefaultParams();
        return;
    end
    
    % Validate the passed grating params
    validateParams(gratingParams);
    
    % Generate the display on which the gratings are rendered
    presentationDisplay = generatePresentationDisplay(gratingParams);
    
    % Generate ths scene sequence depicting frames of the modulated grating
    [theSceneSequence, temporalSupportSeconds] = generateGratingSequence(presentationDisplay, gratingParams, testContrast);
    
    % Assemble dataOut struct - required fields
    dataOut.sceneSequence = theSceneSequence;
    dataOut.temporalSupport = temporalSupportSeconds;
    
    % Add optional fields
    dataOut.presentationDisplay = presentationDisplay;
end


function [theSceneSequence, temporalSupportSeconds] = generateGratingSequence(presentationDisplay, gratingParams, testContrast)

    switch (gratingParams.temporalModulation)
        case 'flashed'
            for frameIndex = 1:gratingParams.temporalModulationParams.stimDurationFramesNum
                % Spatial phase is kept constant throughout all frames
                frameSpatialPhaseSequence(frameIndex) = gratingParams.spatialPhaseDegs;
                
                % Contrast is modulated during the flash frames
                frameContrastSequence(frameIndex) = testContrast;
                if (~ismember(frameIndex, gratingParams.temporalModulationParams.stimOnFrameIndices))
                    frameContrastSequence(frameIndex) = 0;
                end
            end
            
        case 'drifted'
            stimDurationOneCycleSeconds = 1.0/gratingParams.temporalModulationParams.temporalFrequencyHz;
            stimDurationOneCycleFrames = stimDurationOneCycleSeconds / gratingParams.frameDurationSeconds;
            deltaSpatialPhaseDegs = 360/stimDurationOneCycleFrames;
            stimDurationFramesNum = ceil(gratingParams.temporalModulationParams.stimDurationTemporalCycles * stimDurationOneCycleFrames);
            
            for frameIndex = 1:stimDurationFramesNum
                % Contrast is kept constant throughout all frames
                frameContrastSequence(frameIndex) = testContrast;
                % Spatial phase is advanced
                frameSpatialPhaseSequence(frameIndex) = (frameIndex-1)*deltaSpatialPhaseDegs;
            end
    
    
        case 'counter phase modulated'
            stimDurationOneCycleSeconds = 1.0/gratingParams.temporalModulationParams.temporalFrequencyHz;
            stimDurationOneCycleFrames = stimDurationOneCycleSeconds / gratingParams.frameDurationSeconds;
            deltaTemporalPhaseDegs = 360/stimDurationOneCycleFrames;
            stimDurationFramesNum = ceil(gratingParams.temporalModulationParams.stimDurationTemporalCycles * stimDurationOneCycleFrames);
            
            for frameIndex = 1:stimDurationFramesNum
                % Contrast is modulated sinusoidally
                frameContrastSequence(frameIndex) = testContrast * sin(360*(frameIndex-1)*deltaTemporalPhaseDegs);
                % Spatial pahse is kept constant throughout the frames
                frameSpatialPhaseSequence(frameIndex) = gratingParams.spatialPhaseDegs;
            end
            
        otherwise
            error('Unknown temporal modulation: ''%s''.\n', gratingParams.temporalModulation);
    end

    % Open progress bar
    hProgressBar = waitbar(0,'Generating scene sequence...');

    % Generate each frame
    for frameIndex = 1:numel(frameContrastSequence)
        % Update progress bar
        waitbar(frameIndex/numel(frameContrastSequence),hProgressBar,...
            sprintf('Calculating scene frame %d of %d', frameIndex, numel(frameContrastSequence)));
        
        % Generate the scene frame
        theSceneFrame = generateGratingSequenceFrame(presentationDisplay, gratingParams, ...
            frameContrastSequence(frameIndex), frameSpatialPhaseSequence(frameIndex));
        
        % Store the scene frame
        theSceneSequence{frameIndex} = theSceneFrame;
        % Store the timestamp of the scene frame
        temporalSupportSeconds(frameIndex) = (frameIndex-1)*gratingParams.frameDurationSeconds;
    end
    
    % Close progress bar
    close(hProgressBar);
end


function theSceneFrame = generateGratingSequenceFrame(presentationDisplay, gratingParams, frameContrast, frameSpatialPhaseDegs)
    % Compute the color transformation matrices for this display
    displayLinearRGBToLMS = displayGet(presentationDisplay, 'rgb2lms');
    displayLMSToLinearRGB = inv(displayLinearRGBToLMS);
    displayLinearRGBToXYZ = displayGet(presentationDisplay, 'rgb2xyz');
    displayXYZToLinearRGB = inv(displayLinearRGBToXYZ);
    
    % Background chromaticity and mean luminance vector
    xyY = [gratingParams.meanChromaticityXY(1) gratingParams.meanChromaticityXY(2)  gratingParams.meanLuminanceCdPerM2];
    
    % Background XYZ tri-stimulus values
    backgroundXYZ = (xyYToXYZ(xyY(:)))';
    
    % Background linear RGB primary values for the presentation display
    backgroundRGB = imageLinearTransform(backgroundXYZ, displayXYZToLinearRGB);
    
    % Background LMS excitations
    backgroundLMS = imageLinearTransform(backgroundRGB, displayLinearRGBToLMS);
    
    % Compute the spatial contrast modulation pattern
    contrastPattern = frameContrast * generateSpatialModulationPattern(gratingParams, frameSpatialPhaseDegs);
    
    % Compute the LMS-cone contrast spatial pattern
    LMScontrastImage = zeros(size(contrastPattern,1), size(contrastPattern,2), 3);
    for coneIndex = 1:3
        LMScontrastImage(:,:,coneIndex) = gratingParams.coneContrastModulation(coneIndex) * contrastPattern;
    end
    
    % Compute the LMS excitations image
    LMSexcitationImage = bsxfun(@times, (1.0+LMScontrastImage), reshape(backgroundLMS, [1 1 3]));
        
    % Compute the linear RGB primaries image
    RGBimage = imageLinearTransform(LMSexcitationImage, displayLMSToLinearRGB);

    % Make sure we are in gamut (no subpixels with primary values outside of [0 1]
    outOfGamutPixels = numel(find((RGBimage(:)<0)|(RGBimage(:)>1)));
    assert(outOfGamutPixels==0, ...
            sprintf('%d subpixels with primary values > 1; %d subpixels with primary values < 0', ...
            numel(find(RGBimage>1)), numel(find(RGBimage<0))));
        
    % Generate a gamma corrected RGB image (RGBsettings) that we can pop in the isetbio scene straightforward
    RGBsettingsImage = (ieLUTLinear(RGBimage, displayGet(presentationDisplay, 'inverse gamma'))) / displayGet(presentationDisplay, 'nLevels');
    
    % Generate scene corresponding to the test stimulus on the presentation display
    format = 'rgb';
    meanLuminance = []; % EMPTY, so that mean luminance is determined from the rgb settings values we pass
    theScene = sceneFromFile(RGBsettingsImage, format, meanLuminance, presentationDisplay);
        
    % Set the desired FOV
    theSceneFrame = sceneSet(theScene, 'h fov', gratingParams.fovDegs);
end

function contrastPattern = generateSpatialModulationPattern(gratingParams, frameSpatialPhaseDegs)

    % Compute pixelsNum
    cyclesWithinFOV = gratingParams.spatialFrequencyCyclesPerDeg * gratingParams.fovDegs;
    pixelsNumPerCycle = round(gratingParams.pixelsNum /cyclesWithinFOV); 
    if (pixelsNumPerCycle < gratingParams.minPixelsNumPerCycle)
        requiredPixelsNum = round(cyclesWithinFOV * gratingParams.minPixelsNumPerCycle);
        %fprintf('Requested %d pixel stim width, but to satisfy the minPixelsNumPerCycle we need %d pixels\n', ...
        %    gratingParams.pixelsNum, requiredPixelsNum);
        gratingParams.pixelsNum = requiredPixelsNum;
    end
    
    % Spatial support
    x = linspace(-gratingParams.fovDegs/2, gratingParams.fovDegs/2, gratingParams.pixelsNum);
    [X,Y] = meshgrid(x);
    
    % Spatial pattern
    X = X - gratingParams.spatialPositionDegs(1);
    Y = Y - gratingParams.spatialPositionDegs(2);
    Xp =  X * cosd(gratingParams.orientationDegs) + Y * sind(gratingParams.orientationDegs);
    Yp = -X * sind(gratingParams.orientationDegs) + Y * cosd(gratingParams.orientationDegs);
    
    switch (gratingParams.spatialModulationDomain)
        case 'cartesian'
            contrastPattern = cosd(360*gratingParams.spatialFrequencyCyclesPerDeg*Yp + frameSpatialPhaseDegs);
        case 'polar'
            R = sqrt(Xp.^2+Yp.^2);
            contrastPattern = cosd(360*gratingParams.spatialFrequencyCyclesPerDeg*R + frameSpatialPhaseDegs);
        otherwise
            error('Unknown spatial modulationDomain: ''%s''.\n', gratingParams.spatialModulationDomain)
    end
    
    
    
    % Harmonic or Square qave
    switch (gratingParams.spatialModulation)
        case 'square'
            contrastPattern = sign(contrastPattern);
        case 'harmonic'
            ; % do nothing
        otherwise
            error('Unknown spatial modulation: ''%s''.\n', gratingParams.spatialModulation)
    end


    % Envelope
    switch (gratingParams.spatialEnvelope)
        case 'disk'
            idx = find(sqrt(Xp.^2 + Yp.^2) < gratingParams.spatialEnvelopeRadiusDegs);
            envelope = X * 0;
            envelope(idx) = 1;
        case 'Gaussian'
            envelope = exp(-(0.5*(Xp/gratingParams.spatialEnvelopeRadiusDegs).^2)) .* ...
                       exp(-(0.5*(Yp/gratingParams.spatialEnvelopeRadiusDegs).^2));
        case 'square'
            idx = find(...
                (abs(Xp) < gratingParams.spatialEnvelopeRadiusDegs) & ...
                (abs(Yp) < gratingParams.spatialEnvelopeRadiusDegs));
            envelope = Xp * 0;
            envelope(idx) = 1;
        otherwise
            error('Unknown spatial envelope: ''%s''.\n', gratingParams.spatialEnvelope)
    end

    % Mask with envelope
    contrastPattern = contrastPattern .* envelope;
end


function presentationDisplay = generatePresentationDisplay(gratingParams)
    % Generic LCD display
    presentationDisplay = displayCreate('LCD-Apple', ...
        'viewing distance', gratingParams.viewingDistanceMeters, ...
        'wave', gratingParams.spectralSupport ...        % Custom spectral support
        );
    % Custom LUT and bit depth
    if (gratingParams.gammaTableExponent == 1)
        % Linear LUT
        N = 2^gratingParams.bitDepth;
        gTable = repmat(linspace(0, 1, N), 3, 1)';
    else
        x = linspace(0,1,2^gratingParams.bitDepth);
        gTable = x(:).^gratingParams.gammaTableExponent;
        gTable = repmat(gTable, [1,3]);
    end
    % Set the gamma table
    presentationDisplay = displaySet(presentationDisplay, 'gTable', gTable);
end

function validateParams(gratingParamsStruct)
    defaultParamsStruct = generateDefaultParams();
    defaultFieldNames = fieldnames(defaultParamsStruct);
    % Assert that each of the field names of the gratingParamsStruct match
    % one field of the defaultParamsStruct
    
    inputFieldNames = fieldnames(gratingParamsStruct);
    for k = 1:numel(inputFieldNames)
        assert(ismember(inputFieldNames{k}, defaultFieldNames), ...
            sprintf('>> Can not deal with grating parameter name ''%s''.\n>> Inspect the generateDefaultParams() function of %s.m for all available parameter names.', ...
            inputFieldNames{k}, mfilename()));
    end
    
    % Further validations
    % spatialEnvelope
    if (isfield(gratingParamsStruct, 'spatialEnvelope'))
        assert(ismember(gratingParamsStruct.spatialEnvelope, {'disk', 'square', 'Gaussian'}), ...
            sprintf('>> Invalid value for spatialModulation: ''%s''.\n>> Inspect the generateDefaultParams() function of %s.m for valid parameter values.', ...
            gratingParamsStruct.spatialModulation, mfilename()));
    end
    
    % spatialModulation
    if (isfield(gratingParamsStruct, 'spatialModulation'))
        assert(ismember(gratingParamsStruct.spatialModulation, {'harmonic', 'square'}), ...
            sprintf('>> Invalid value for spatialModulation: ''%s''.\n>> Inspect the generateDefaultParams() function of %s.m for valid parameter values.', ...
            gratingParamsStruct.spatialModulation, mfilename()));
    end
    
    % temporalModulationParams
    if (isfield(gratingParamsStruct, 'temporalModulation'))
        assert(ismember(gratingParamsStruct.temporalModulation, {'flashed', 'drifted', 'counter phase modulated'}), ...
            sprintf('>> Invalid value for temporalModulation: ''%s''.\n>> Inspect the generateDefaultParams() function of %s.m for valid parameter values.', ...
            gratingParamsStruct.temporalModulation, mfilename()));
    end
end

function p = generateDefaultParams()

    p = struct(...
        'viewingDistanceMeters', 0.57, ...              % display: viewing distance
        'bitDepth', 8, ...                              % display: length of LUT
        'gammaTableExponent', 2.0, ...                  % display: shape of LUT, 1 = Linear
        'spectralSupport', 400:10:750, ...              % display: spectral support of the primary SPDs, in nanometers
        'meanLuminanceCdPerM2', 40, ...                 % background: mean luminance, in candellas per meter squared
        'meanChromaticityXY', [0.3 0.32], ...           % background: neutral mean chromaticity
        'coneContrastModulation', [0.09 -0.09 0.0], ...   % chromatic direction: LMS cone contrasts
        'fovDegs', 1.0, ...                             % spatial: field of view, in degrees
        'minPixelsNumPerCycle', 10, ...                 % spatial: min number of pixels per grating spatial period (to minimize aliasing effects)
        'pixelsNum', 256, ...                           % spatial: desired size of stimulus in pixels - this could be larger depending on the minPixelsNumPerCycle, the SF, and the FOV
        'spatialFrequencyCyclesPerDeg', 5, ...          % spatial: grating spatial freequency, in cycles/deg
        'orientationDegs', 0, ...                       % spatial: grating orientation, in degrees
        'spatialPhaseDegs', 90, ...                     % spatial: grating spatial phase, in degrees
        'spatialPositionDegs', [0.0 0.0], ...           % spatial: center of grating, in degrees
        'spatialModulation', 'harmonic', ...            % spatial: contrast modulation - choose between {'harmonic', 'square'} for sinusoidal/square spatial modulation 
        'spatialModulationDomain', 'cartesian', ...     % spatial: domain of spatial modulation - choose between {'cartesian', 'polar'}
        'spatialEnvelope', 'Gaussian', ...              % spatial: envelope - choose between {'disk', 'square', 'Gaussian'}
        'spatialEnvelopeRadiusDegs', 0.15, ...          % spatial: radius of the spatial envelope, in degs    
        'temporalModulation', 'flashed', ...            %   temporal modulation mode: choose between {'flashed', 'drifted', 'counter phase modulated'}
        'temporalModulationParams', struct(...          % temporal: modulation params struct
            'stimOnFrameIndices', [2 3], ...            %   params relevant to the temporalModulationMode
            'stimDurationFramesNum', 4), ...            %   params relevant to the temporalModulationMode
        'frameDurationSeconds', 20/1000 ...            % temporal: frame duration, in seconds
    );

end