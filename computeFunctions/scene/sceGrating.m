function dataOut = sceGrating(sceneEngineOBJ, testContrast, gratingParams)
% Compute function for generating a sequence of scenes depicting a
% a grating that is either drifting, counter-phase modulated, or flashed.
%
% Syntax:
%   dataOut = sceGrating(obj, testContrast, gratingParams);
%
% Description:
%    Compute function to be used as a computeFunctionHandle for a @sceneEngine
%    object. There are 2 ways to use this function.
%
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
% Inputs:
%    sceneEngineOBJ              - Calling @sceneEngine object.  This is
%                                  currently used, but passing it allows us
%                                  flexibility in the future and matches
%                                  conventions for the other classes in
%                                  this toolbox.
%    testContrast                - Scalar providing the contrast for the
%                                  scene to be generated.                           
%    sceneParamsStruct           - Struct containing properties of the
%                                  scene understood by this function.
%                                  As noted above, execute sceGrating at
%                                  the command line to see the structure's
%                                  fields and default values.
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
% Examples:
%    The source code contains examples.
%
% See Also:
%     t_sceneGeneration, t_thresholdEngine

% History:
%    10/05/2020  NPC  Wrote it.

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
    [theSceneSequence, temporalSupportSeconds, statusReport] = generateGratingSequence(presentationDisplay, gratingParams, testContrast);
    
    % Assemble dataOut struct - required fields
    dataOut.sceneSequence = theSceneSequence;
    dataOut.temporalSupport = temporalSupportSeconds;
    dataOut.statusReport = statusReport;
    
    % Add optional fields
    dataOut.presentationDisplay = presentationDisplay;
end

function [theSceneSequence, temporalSupportSeconds, statusReport] = generateGratingSequence(presentationDisplay, gratingParams, testContrast)

    switch (gratingParams.temporalModulation)
        case 'flashed'
            for frameIndex = 1:gratingParams.temporalModulationParams.stimDurationFramesNum
                % Spatial phase is kept constant throughout all frames
                frameSpatialPhaseSequence(frameIndex) = gratingParams.spatialPhaseDegs;
                
                % Contrast is modulated during the flash frames
                % Initially we set the stimulus contrast of all the frames
                % to testContrast
                frameContrastSequence(frameIndex) = testContrast;
                % if the stimulus should disappear in any given frame, then
                % the contrast is reset to 0
                if (~ismember(frameIndex, gratingParams.temporalModulationParams.stimOnFrameIndices))
                    frameContrastSequence(frameIndex) = 0;
                end
            end
            
        case 'drifted'
            % Compute the time (in sec) needed to complete 1 cycle
            stimDurationOneCycleSeconds = 1.0/gratingParams.temporalModulationParams.temporalFrequencyHz;
            % Compute the time (in frame) needed to compute 1 cycle
            stimDurationOneCycleFrames = stimDurationOneCycleSeconds / gratingParams.frameDurationSeconds;
            % Compute the phase needed to shift for each frame
            deltaSpatialPhaseDegs = 360/stimDurationOneCycleFrames;
            % Compute the total number of frames 
            % = total drifting cycles x number of frames needed for 1 cycle
            stimDurationFramesNum = ceil(gratingParams.temporalModulationParams.stimDurationTemporalCycles * stimDurationOneCycleFrames);
            
            for frameIndex = 1:stimDurationFramesNum
                % Contrast is kept constant throughout all frames
                frameContrastSequence(frameIndex) = testContrast;
                % Spatial phase is advanced
                frameSpatialPhaseSequence(frameIndex) = (frameIndex-1)*deltaSpatialPhaseDegs;
            end
    
    
        case 'counter phase modulated'
            % See the comments for case 'drifted' for the four lines below
            stimDurationOneCycleSeconds = 1.0/gratingParams.temporalModulationParams.temporalFrequencyHz;
            stimDurationOneCycleFrames = stimDurationOneCycleSeconds / gratingParams.frameDurationSeconds;
            deltaTemporalPhaseDegs = 360/stimDurationOneCycleFrames;
            stimDurationFramesNum = ceil(gratingParams.temporalModulationParams.stimDurationTemporalCycles * stimDurationOneCycleFrames);

            for frameIndex = 1:stimDurationFramesNum
                % Contrast is modulated sinusoidally 
                frameContrastSequence(frameIndex) = testContrast * sind((frameIndex-1)*deltaTemporalPhaseDegs);
                % Spatial pahse is kept constant throughout the frames
                frameSpatialPhaseSequence(frameIndex) = gratingParams.spatialPhaseDegs;
            end

        otherwise
            error('Unknown temporal modulation: ''%s''.\n', gratingParams.temporalModulation);
    end

    displayProgressBar = false; % getpref('ISET','waitbar');
    if (displayProgressBar)
        % Open progress bar
        hProgressBar = waitbar(0,'Generating scene sequence...');
    end
    
    theSceneSequence = cell(1, numel(frameContrastSequence));
    temporalSupportSeconds = zeros(1, numel(frameContrastSequence));
    outOfGamutFlag = zeros(1, numel(frameContrastSequence));
    
    % Generate each frame
    for frameIndex = 1:numel(frameContrastSequence)
        % Update progress bar
        if (displayProgressBar)
            waitbar(frameIndex/numel(frameContrastSequence),hProgressBar,...
                sprintf('Calculating scene frame %d of %d', frameIndex, numel(frameContrastSequence)));
        end
        
        % Generate the scene frame
        [theSceneFrame, outOfGamutFlag(frameIndex)] = generateGratingSequenceFrame(presentationDisplay, gratingParams, ...
            frameContrastSequence(frameIndex), frameSpatialPhaseSequence(frameIndex));
        
        % Store the scene frame
        theSceneSequence{frameIndex} = theSceneFrame;
        % Store the timestamp of the scene frame
        temporalSupportSeconds(frameIndex) = (frameIndex-1)*gratingParams.frameDurationSeconds;
    end
    
    % Set the outOfGamutFlag if we went out of gamut in any of the frames
    if (any(outOfGamutFlag))
        statusReport.OutOfGamut = true;
    else
        statusReport = struct();
    end

    if (displayProgressBar)
        % Close progress bar
        close(hProgressBar);
    end
end

function [theSceneFrame, outOfGamutFlag] = generateGratingSequenceFrame(presentationDisplay, gratingParams, frameContrast, frameSpatialPhaseDegs)
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

    % Flag on whether to debug the LUT inversion visually
    debugLUT = ~true;
    if (debugLUT)
        figure(256)
        subplot(1,2,1)
        imagesc(squeeze(sum(RGBimage,3)))
        colormap(gray(1024))
    end
    
    % Make sure we are in gamut (no subpixels with primary values outside of [0 1]
    outOfGamutFlag = false;
    outOfGamutPixels = numel(find((RGBimage(:)<0)|(RGBimage(:)>1)));

    if (isfield(gratingParams, 'warningInsteadOfErrorOnOutOfGamut'))
        if (outOfGamutPixels>0)
            if (gratingParams.warningInsteadOfErrorOnOutOfGamut)
               % Tone map and do not throw an error
               outOfGamutFlag = true;
               RGBimage = tonemap(RGBimage);
            else
               % Throw an error
               error('%d subpixels with primary values > 1; %d subpixels with primary values < 0', ...
                    numel(find(RGBimage>1)), numel(find(RGBimage<0)));
            end
        end
    else
        % Assert we are in-gamut
        assert(outOfGamutPixels==0, ...
            sprintf('%d subpixels with primary values > 1; %d subpixels with primary values < 0', ...
            numel(find(RGBimage>1)), numel(find(RGBimage<0))));
    end
        
    % Generate a gamma corrected RGB image (RGBsettings) that we can pop in the isetbio scene straightforward
    lutLevels = displayGet(presentationDisplay, 'nLevels');
    RGBsettingsImage = (ieLUTLinear(RGBimage, ieLUTInvert(displayGet(presentationDisplay, 'gamma'), lutLevels))) / lutLevels;
    
    if (debugLUT)
        subplot(1,2,2)
        imagesc(squeeze(sum(RGBsettingsImage,3)))
        colormap(gray)
        pause
    end
    
    % Generate scene corresponding to the test stimulus on the presentation display
    format = 'rgb';
    meanLuminance = []; % EMPTY, so that mean luminance is determined from the rgb settings values we pass
    theScene = sceneFromFile(RGBsettingsImage, format, meanLuminance, presentationDisplay);
        
    % Set the desired FOV
    theSceneFrame = sceneSet(theScene, 'h fov', gratingParams.fovDegs);

    % Check if the spectral support and transmission data are available for the filter
    if ~isempty(gratingParams.filter.spectralSupport) && ~isempty(gratingParams.filter.transmission)
        % Replicate the spectral support of the gratingParams across rows equal to the 
        % length of the filter's spectral support
        wvl = repmat(gratingParams.spectralSupport(:)', ...     %row vector
            [length(gratingParams.filter.spectralSupport), 1]);

        % Replicate the spectral support of the filter across columns equal to the 
        % length of the grating's spectral support
        wvl_filter = repmat(gratingParams.filter.spectralSupport(:),...  %column vector
            [1, length(gratingParams.spectralSupport)]);

        % Compute the absolute difference between the replicated wavelength grids
        Diff = abs(wvl_filter - wvl);

        % Find the indices of the minimum values in each column of the difference matrix
        [~, minDiff_idx] = min(Diff, [], 1);

        % Use the indices to select corresponding transmission values and convert 
        % from percentage to decimal
        transmission_wvl = gratingParams.filter.transmission(minDiff_idx)./100; 

        % Retrieve the photon data from the current scene frame
        photons = sceneGet(theSceneFrame,'photons');

        % Replicate the adjusted transmission values across a grid defined by 
        % pixelsNum x pixelsNum, and then permute to adjust dimensions
        transmission_wvl_rep = permute(repmat(transmission_wvl(:),...
            [1, size(photons,1), size(photons,2)]),[2,3,1]);

        % Apply the filter's transmission to the photon data
        photons_filtered = photons.*transmission_wvl_rep;

        % Update the scene frame with the filtered photons
        theSceneFrame = sceneSet(theSceneFrame,'photons',photons_filtered);
    end
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
        % Make the sinusoidal grating in the cartesian coordinate system
        case 'cartesian'
            contrastPattern = cosd(360*gratingParams.spatialFrequencyCyclesPerDeg*Yp + frameSpatialPhaseDegs);
        % Make the sinusoidal grating in the polar coordinate system
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
            % do nothing
        otherwise
            error('Unknown spatial modulation: ''%s''.\n', gratingParams.spatialModulation)
    end

    % Envelope
    switch (gratingParams.spatialEnvelope)
        case 'disk'
            % select all the indices that are within the circle
            idx = find(sqrt(Xp.^2 + Yp.^2) < gratingParams.spatialEnvelopeRadiusDegs);
            % make what's within the circle 1 and everywhere else 0
            % this leads to sharp edges
            envelope = X * 0;
            envelope(idx) = 1;
        case 'soft'
            % soft Gaussian envelope
            envelope = exp(-(0.5*(Xp/gratingParams.spatialEnvelopeRadiusDegs).^2)) .* ...
                       exp(-(0.5*(Yp/gratingParams.spatialEnvelopeRadiusDegs).^2));
        case 'rect'
            % rectangular envelope with sharp edges
            idx = find(...
                (abs(X) < gratingParams.spatialEnvelopeRadiusDegs) & ...
                (abs(Y) < gratingParams.spatialEnvelopeRadiusDegs));
            envelope = Xp * 0;
            envelope(idx) = 1;
        case 'halfcos'
            % the positive half of cosine wave is selected (-pi/2 to pi/2)
            Xp = (X + gratingParams.spatialPositionDegs(1)).* (pi/(2*gratingParams.spatialEnvelopeRadiusDegs));
            Yp = (Y + gratingParams.spatialPositionDegs(2)).* (pi/(2*gratingParams.spatialEnvelopeRadiusDegs));
            envelope = cos(Xp).* cos(Yp);
            idx = find(Xp < -pi / 2 | Xp > pi / 2 ...
                 | Yp < -pi/2 | Yp > pi / 2);
            envelope(idx) = 0;
            envelope = envelope / max(envelope(:));
        case 'none'
            envelope = Xp * 0 + 1;
        otherwise
            error('Unknown spatial envelope: ''%s''.\n', gratingParams.spatialEnvelope)
    end

    % Mask with envelope
    contrastPattern = contrastPattern .* envelope;
end

function presentationDisplay = generatePresentationDisplay(gratingParams)
    % Generic LCD display
    presentationDisplay = displayCreate(gratingParams.screenDisplay, ...
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
    
    % Adjust display SPDs so we can generate desired mean luminance.
    %
    % Why the factor of 2.1 here?  What we (DHB/FH) think is that this
    % means is that when the mean luminance is set, we can present a 100%
    % contrast luminance grating in gamut, with a little headroom.  To put
    % it another way, you need the factor of 2, and the factor of 2.1
    % provides a little bit of safety margin.
    desiredPeakLuminance = gratingParams.meanLuminanceCdPerM2 * 2.1;
    peakLuminance = displayGet(presentationDisplay, 'peak luminance');
    scalingFactor = desiredPeakLuminance / peakLuminance;
    spds = displayGet(presentationDisplay, 'spd primaries');
    presentationDisplay = displaySet(presentationDisplay, 'spd primaries', spds*scalingFactor);
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
        assert(ismember(gratingParamsStruct.spatialEnvelope, {'none', 'disk', 'rect', 'soft','halfcos'}), ...
            sprintf('>> Invalid value for spatialEnvelope: ''%s''.\n>> Inspect the generateDefaultParams() function of %s.m for valid parameter values.', ...
            gratingParamsStruct.spatialEnvelope, mfilename()));
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
        'screenDisplay', 'LCD-Apple',...                % display: the screen display type
        'filter',struct(...                             % spatial: apply a filter in front of the display 
            'spectralSupport',[],...                    %   A vector containing the wavelengths (nm) at which the filter's transmission properties are defined 
            'transmission',[]),...                      %   A corresponding vector of transmission values (%; ranging from 0 - 100) for each wavelength specified in 'spectralSupport' 
        'viewingDistanceMeters', 0.57, ...              % display: viewing distance
        'bitDepth', 20, ...                             % display: length of LUT
        'gammaTableExponent', 2.0, ...                  % display: shape of LUT, 1 = Linear
        'spectralSupport', 400:10:750, ...              % display: spectral support of the primary SPDs, in nanometers
        'meanLuminanceCdPerM2', 40, ...                 % background: mean luminance, in candellas per meter squared
        'meanChromaticityXY', [0.3 0.32], ...           % background: neutral mean chromaticity
        'coneContrastModulation', [0.09 -0.09 0.0], ... % chromatic direction: LMS cone contrasts
        'warningInsteadOfErrorOnOutOfGamut', false, ... % chromatic: whether to throw an error or a warning if out-of-gamut pixels are generated
        'fovDegs', 1.0, ...                             % spatial: field of view, in degrees
        'minPixelsNumPerCycle', 10, ...                 % spatial: min number of pixels per grating spatial period (to minimize aliasing effects)
        'pixelsNum', 256, ...                           % spatial: desired size of stimulus in pixels - this could be larger depending on the minPixelsNumPerCycle, the SF, and the FOV
        'spatialFrequencyCyclesPerDeg', 5, ...          % spatial: grating spatial freequency, in cycles/deg
        'orientationDegs', 0, ...                       % spatial: grating orientation, in degrees
        'spatialPhaseDegs', 90, ...                     % spatial: grating spatial phase, in degrees
        'spatialPositionDegs', [0.0 0.0], ...           % spatial: center of grating, in degrees
        'spatialModulation', 'harmonic', ...            % spatial: contrast modulation - choose between {'harmonic', 'square'} for sinusoidal/square spatial modulation 
        'spatialModulationDomain', 'cartesian', ...     % spatial: domain of spatial modulation - choose between {'cartesian', 'polar'}
        'spatialEnvelope', 'soft', ...                  % spatial: envelope - choose between {'disk', 'rect', 'soft'}
        'spatialEnvelopeRadiusDegs', 0.15, ...          % spatial: radius of the spatial envelope, in degs    
        'temporalModulation', 'flashed', ...            % temporal modulation mode: choose between {'flashed', 'drifted', 'counter phase modulated'}
        'temporalModulationParams', struct(...          % temporal: modulation params struct
            'stimOnFrameIndices', [2 3], ...            %   params relevant to the temporalModulationMode
            'stimDurationFramesNum', 4), ...            %   params relevant to the temporalModulationMode
        'frameDurationSeconds', 20/1000 ...            % temporal: frame duration, in seconds
    );

end

function RGBimage = tonemap(RGBimage)
    % Transform the RGB image [Rows x Cols x 3] into a [3 x N] matrix for faster computations
    [rgbPrimariesImageCalFormat,nCols,mRows] = ImageToCalFormat(RGBimage);
    
    % Tonemap
    maxes = max(rgbPrimariesImageCalFormat,[],1);
    maxes = max(maxes, maxes*0+1);
    mins = min(rgbPrimariesImageCalFormat,[],1);
    mins = min(mins, maxes*0);
    rgbPrimariesImageCalFormat = bsxfun(@times,bsxfun(@minus,rgbPrimariesImageCalFormat,mins), 1./(maxes - mins));
    
    % Back to image format
    RGBimage = CalFormatToImage(rgbPrimariesImageCalFormat, nCols, mRows);
end

