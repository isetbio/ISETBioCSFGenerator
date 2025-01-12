function dataOut = sceAOTwoSpot(sceneEngineOBJ, testContrast, twoSpotParams)
% Compute function for generating a sequence of scenes for two AO spots
%
% Syntax:
%   dataOut = sceAOTwoSpot(obj, testContrast, twoSpotParams)
%
% Description:
%    Compute function to be used as a computeFunctionHandle for a @sceneEngine
%    object. There are 2 ways to use this function.
%
%       [1] If called wihtout arguments, itt does not compute anything and
%           simply returns a struct with the defaultParams that define the
%           scene.
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
%                                  currently not used, but passing it allows us
%                                  flexibility in the future and matches
%                                  conventions for the other classes in
%                                  this toolbox.
%    testContrast                - Scalar providing the contrast for the
%                                  scene to be generated.
%    twoSpotParams               - Struct containing properties of the
%                                  scene understood by this function.
%                                  As noted above, execute sceGrating at
%                                  the command line to see the structure's
%                                  fields and default values.
%
% Outputs:
%    dataOut  - A struct that depends on the input arguments.
%
%             - If called directly with no input arguments, the returned struct contains
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
%     t_sceneGeneration, t_spatialCSF

% History:
%    03/16/2021  dhb  Wrote it.

%   Examples:
%{
    % Usage case #1. Just return the default scene params
    defaultParams = sceAOTwoSpot()

    % Usage case #2. Compute a grating sequence using a parent @sceneEngine
    % object and the default scene params

    % Instantiate the parent @sceneEngine object
    theSceneEngineOBJ = sceneEngine(@sceAOTwoSpot);

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
validateParams(twoSpotParams);

% Generate ths scene sequence depicting frames of the modulated grating
[theSceneSequence, temporalSupportSeconds] = generateTwoSpotSequence(twoSpotParams, testContrast);

% Check the visualizeEachCompute flag of the sceneEngineOBJ, and if set to true,
% call its visualizeSceneSequence() method to visualize the generate scene sequence.
if (sceneEngineOBJ.visualizeEachCompute)
    sceneEngineOBJ.visualizeSceneSequence(theSceneSequence, temporalSupportSeconds);
end

% Assemble dataOut struct - required fields
dataOut.sceneSequence = theSceneSequence;
dataOut.temporalSupport = temporalSupportSeconds;
dataOut.statusReport = struct();

% Add optional fields to dataOut structure. None presently.

end

function [theSceneSequence, temporalSupportSeconds] = generateTwoSpotSequence(twoSpotParams, testContrast)

switch (twoSpotParams.temporalModulation)
    case 'flashed'
        for frameIndex = 1:twoSpotParams.temporalModulationParams.stimDurationFramesNum
            % Can set frame specific values here
            %
            % Contrast is modulated during the flash frames
            frameContrastSequence(frameIndex) = testContrast;
            if (~ismember(frameIndex, twoSpotParams.temporalModulationParams.stimOnFrameIndices))
                frameContrastSequence(frameIndex) = 0;
            end
        end
        
    otherwise
        error('Unknown temporal modulation: ''%s''.\n', twoSpotParams.temporalModulation);
end

% Visualize scene
% vcAddAndSelectObject(theScene);
% sceneWindow;

% Set up
theSceneSequence = cell(1, numel(frameContrastSequence));
temporalSupportSeconds = zeros(1, numel(frameContrastSequence));

% Generate each frame
for frameIndex = 1:numel(frameContrastSequence)
    % Generate the scene frame
    theSceneFrame = generateTwoSpotSequenceFrame(twoSpotParams,frameContrastSequence(frameIndex));
    
    % Store the scene frame
    theSceneSequence{frameIndex} = theSceneFrame;
    
    % Store the timestamp of the scene frame
    temporalSupportSeconds(frameIndex) = (frameIndex-1)*twoSpotParams.frameDurationSeconds;
end

end

function theSceneFrame = generateTwoSpotSequenceFrame(twoSpotParams,frameContrast)

% Make relative spectral power distributions. Each is approximated by a
% Gaussian with specified center wavlength and FWHM, and with total power
% given by the corneal power specified above.  The call to trapz takes
% wavelength spacing into account when normalizing the power.
imagingRelSpd = normpdf(twoSpotParams.wls,twoSpotParams.imagingWl,FWHMToStd(twoSpotParams.imagingFWHM));
imagingUnitSpd = imagingRelSpd/trapz(twoSpotParams.wls,imagingRelSpd);
spotRelSpd = normpdf(twoSpotParams.wls,twoSpotParams.spotWl,FWHMToStd(twoSpotParams.spotFWHM));
spotUnitSpd = spotRelSpd/trapz(twoSpotParams.wls,spotRelSpd);

% Stimulus power
%
% Specified here as the power entering the eye (going through
% the pupil).  If the beam just fills the pupil, then this is just the power
% measured at the cornea.
%
% If you have instead corneal irradiance in UW/mm^2 and this fills or overfills
% the pupil, then multiply by pupil area in mm^2 to get the power entering the eye.
%
% This code assumes that the power for both background and test is spread
% out over the background area. This is because typically for an AOSLO we
% measure power for the whole raster, and then in an experiment switch it
% off for part of the time to make the spatial pattern we want.
%
% The test should be specified as the amount of power added to
% the background (its incremental power).  This is what you'll get if
% background and test are in separate channels of the system and you
% measure them separately.
imagingBgCornealPowerUW = twoSpotParams.imagingBgPowerUW;
spotBgCornealPowerUW = twoSpotParams.spotBgPowerUW;

% Make spds that give desired full power.
backgroundSpdCornealPowerUW = imagingBgCornealPowerUW*imagingUnitSpd;
spotSpdCornealPowerUW = spotBgCornealPowerUW*spotUnitSpd;

% Get equivalent spectral radiance of background and test increment as a
% function of the wavelength support.
%
% The routine here finds the radiance on an external conventional display that
% produces the same retinal illuminance as the corneal power specified
% above.  This is purely geometric calculation; attenuation of light by
% occular media is not taken into account at this stage.  Note that this
% conversion routine expects power per wavelength band, not power per nm,
% as its input, but returns power per nm as its output.  A little
% confusing, but that is why the spd being passed in is multiplied by the
% wavelength spacing.
deltaWl = twoSpotParams.wls(2)-twoSpotParams.wls(1);
backgroundSizeDegs = twoSpotParams.spotBgDegs;
pupilDiameterMm = twoSpotParams.pupilDiameterMm;
pupilAreaMm2 = pi*((pupilDiameterMm/2)^2);
imagingBgSpdRadiance = AOMonochromaticCornealPowerToRadiance(twoSpotParams.wls,twoSpotParams.wls,backgroundSpdCornealPowerUW*deltaWl,pupilDiameterMm,backgroundSizeDegs^2);
spotBgSpdRadiance = AOMonochromaticCornealPowerToRadiance(twoSpotParams.wls,twoSpotParams.wls,spotSpdCornealPowerUW*deltaWl,pupilDiameterMm,backgroundSizeDegs^2);

% Make sure our computed radiance yields the desired corneal
% irradiance when we go in the other direction.
imagingBgSpdCornealIrradianceUWMm2Check = RadianceAndDegrees2ToCornIrradiance(imagingBgSpdRadiance,backgroundSizeDegs^2)*(1e6)*((1e-3)^2);
imagingCornealPowerUWCheck = trapz(twoSpotParams.wls,imagingBgSpdCornealIrradianceUWMm2Check)*pupilAreaMm2;
if (abs(imagingCornealPowerUWCheck-imagingBgCornealPowerUW)/imagingBgCornealPowerUW > 1e-4)
    error('Do not get right cornal power back from computed radiance');
end
spotSpdCornealIrradianceUWMm2Check = RadianceAndDegrees2ToCornIrradiance(spotBgSpdRadiance,backgroundSizeDegs^2)*(1e6)*((1e-3)^2);
spotCornealPowerUWCheck = trapz(twoSpotParams.wls,spotSpdCornealIrradianceUWMm2Check)*pupilAreaMm2;
if (abs(spotCornealPowerUWCheck-spotBgCornealPowerUW)/spotBgCornealPowerUW > 1e-4)
    error('Do not get right cornal power back from computed radiance');
end

% Create scene
%
% Create an empty scene to use for the spot.  Put it far enough away so it
% is basically in focus for an emmetropic eye accommodated to infinity.
theScene = sceneCreate('empty');
theScene = sceneSet(theScene,'wavelength',twoSpotParams.wls);
theScene = sceneSet(theScene,'distance',twoSpotParams.viewingDistanceMeters);

% Make an image with the background + two spots
spotPattern = drawTwoSpot(twoSpotParams);

% Get stimulus 1 and 2 contrast
spot1Contrast = frameContrast*cosd(twoSpotParams.stimAngle);
spot2Contrast = frameContrast*sind(twoSpotParams.stimAngle);
spot1SpdRadiance = spotBgSpdRadiance + spot1Contrast*spotBgSpdRadiance;
spot2SpdRadiance = spotBgSpdRadiance + spot2Contrast*spotBgSpdRadiance;

% Then fill in appropriate radiance at each pixel
nWls = length(twoSpotParams.wls);
radianceEnergySpot = zeros(twoSpotParams.pixelsNum,twoSpotParams.pixelsNum,nWls);
for i = 1:twoSpotParams.pixelsNum
    for j = 1:twoSpotParams.pixelsNum
        % Background pixels are 1, spot1 pixels are 2, spot2pixels are 3
        if (spotPattern(i,j) == 0)
            radianceEnergySpot(i,j,:) = imagingBgSpdRadiance;
        elseif (spotPattern(i,j) == 1)
            radianceEnergySpot(i,j,:) = imagingBgSpdRadiance + spotBgSpdRadiance;
        elseif (spotPattern(i,j) == 2)
            radianceEnergySpot(i,j,:) = imagingBgSpdRadiance + spot1SpdRadiance;
        elseif (spotPattern(i,j) == 3)
            radianceEnergySpot(i,j,:) = imagingBgSpdRadiance + spot2SpdRadiance;
        end
    end
end

% Convert radiance to quantal units
radiancePhotonsSpot = Energy2Quanta(twoSpotParams.wls,radianceEnergySpot);

% Put in the image in quantal units (photons)
theScene = sceneSet(theScene,'photons',radiancePhotonsSpot);

% Set the desired FOV
theSceneFrame = sceneSet(theScene, 'h fov', twoSpotParams.fovDegs);

end

function validateParams(twoSpotParamsStruct)
defaultParamsStruct = generateDefaultParams();
defaultFieldNames = fieldnames(defaultParamsStruct);
% Assert that each of the field names of the twoSpotParamsStruct match
% one field of the defaultParamsStruct

inputFieldNames = fieldnames(twoSpotParamsStruct);
for k = 1:numel(inputFieldNames)
    assert(ismember(inputFieldNames{k}, defaultFieldNames), ...
        sprintf('>> Can not deal with grating parameter name ''%s''.\n>> Inspect the generateDefaultParams() function of %s.m for all available parameter names.', ...
        inputFieldNames{k}, mfilename()));
end

% Further validations
%
% type
if (isfield(twoSpotParamsStruct, 'type'))
    assert(ismember(twoSpotParamsStruct.type, {'basic'}), ...
        sprintf('>> Invalid value for temporalModulation: ''%s''.\n>> Inspect the generateDefaultParams() function of %s.m for valid parameter values.', ...
        twoSpotParamsStruct.type, mfilename()));
end

% temporalModulationParams
if (isfield(twoSpotParamsStruct, 'temporalModulation'))
    assert(ismember(twoSpotParamsStruct.temporalModulation, {'flashed'}), ...
        sprintf('>> Invalid value for temporalModulation: ''%s''.\n>> Inspect the generateDefaultParams() function of %s.m for valid parameter values.', ...
        twoSpotParamsStruct.temporalModulation, mfilename()));
end

end

function p = generateDefaultParams()

p = struct(...
    'type', 'basic', ...                            % type
    'viewingDistanceMeters', 2, ...                 % viewing distance: in meters
    'wls', 400:10:750, ...                          % wavelength support of the primary SPDs: in nanometers
    'stimAngle', 0, ...                             % stimulus angle in incr/decr plane
    'spotWl', 550, ...                              % spot wavelength: in nm
    'spotFWHM', 20, ...                             % spot full width at half max: in nm
    'spotWidthDegs', 1.3/60, ...                    % spot width: in degrees
    'spotHeightDegs', 1/60, ...                     % spot height: in degrees
    'spotVerticalSepDegs', 1/60, ...                % spot center to center vertical separation: in degrees
    'spotHorizontalSepDegs', 0, ...                 % spot center to center horizontal separation: in degrees
    'spotBgDegs', 0.25, ...                         % spot background: in degrees
    'spotBgPowerUW', (1/16)*76/1000*(2/3), ...      % spot background power: in uW
    'imagingWl', 790, ...                           % imaging beam wavelength: in nm
    'imagingFWHM', 20, ...                          % imaging beam full width at half max: in nm
    'imagingBgPowerUW', 0.5, ...                      % imaging beam power entering eye: in uW
    'fovDegs', 0.25, ...                            % spatial: scene field of view, in degrees
    'pixelsNum', 128, ...                           % spatial: desired size of stimulus in pixels
    'temporalModulation', 'flashed', ...            % temporal modulation mode: choose between {'flashed'}
    'temporalModulationParams', struct(...          % temporal: modulation params struct
    'stimOnFrameIndices', [2 3 4], ...              %   params relevant to the temporalModulationMode
    'stimDurationFramesNum', 5), ...                %   params relevant to the temporalModulationMode
    'frameDurationSeconds', 1/16, ...               % temporal: frame duration, in seconds
    'pupilDiameterMm', 6 ...                        % pupil diameter mm
    );

end
