% t_AOTumblingEScene
%
% Description:
%   Produce an AO apparatus scene with a tumbling E stimulus whose
%   orientation we can specify and where the variable that is changed is
%   the size.
%
% See also:
%    t_AORectangleScene.

% History:
%    06/25/24  dhb, qf  Wrote it.

%% Clear and close
clear; close all;

%% Set parameters 
defocusDiopters = 0.05;
pupilDiameterMm = 7;

% Orientation of the E.
%
% Natural choices here are 0, 90, 180 and 270.
letterRotationDegrees = 0;

% Other stimulus parameters
%
% Define basic parameters of the AO stimulus,
% We can put in arbitrary spectra but the
% actual E we want to model was monochromatic.
wls = 600:10:800;
eWavelengthNm = 780;

% Spatial parameters
nPixels = 128;
fieldSizeMinutes = 60;
fieldSizeDegs = fieldSizeMinutes/60;

% Stimulus angle
stimAngle = 0;

% Background power in uW.  Depending on how Will ships this to us,
% we may need to convert into these units.
fieldPowerUWPerDeg2 = 2000;
fieldPowerUW = (fieldSizeDegs^2)*fieldPowerUWPerDeg2;

% E power in uW.  This is the power for the stimulus that
% is added to the background.  Let's assume for now that
% it is twice the background power.
ePowerUW = 2*fieldPowerUW;

% Set up two spot params
rectParams = struct(...
    'type', 'basic', ...                            % type
    'viewingDistanceMeters', 2, ...                 % viewing distance: in meters
    'wls', 400:10:750, ...                          % wavelength support of the primary SPDs: in nanometers
    'stimAngle', stimAngle, ...                     % stimulus angle in incr/decr plane
    'spotWl', eWavelengthNm, ...                    % spot wavelength: in nm
    'spotFWHM', 20, ...                             % spot full width at half max: in nm
    'spotWidthDegs', spotWidthMinutes/60, ...       % spot width: in degrees
    'spotHeightDegs', spotHeightMinutes/60, ...     % spot height: in degrees
    'spotVerticalSepDegs', spotVerticalSepMinutes/60 , ... % spot center to center vertical separation: in degrees
    'spotHorizontalSepDegs', spotHorizontalSepMinutes/60, ... % spot center to center horizontal separation: in degrees
    'spotBgDegs', fieldSizeDegs, ...                % spot background: in degrees
    'spotBgPowerUW', spotBgPowerUW, ...             % spot background power: in uW
    'imagingWl', 750, ...                           % imaging beam wavelength: in nm
    'imagingFWHM', 20, ...                          % imaging beam full width at half max: in nm
    'imagingBgPowerUW', 0, ...                      % imaging beam power entering eye: in uW
    'fovDegs', fieldSizeDegs, ...                   % spatial: full scene field of view, in degrees
    'pixelsNum', nPixels, ...                       % spatial: desired size of stimulus in pixels
    'temporalModulation', 'flashed', ...            % temporal modulation mode: choose between {'flashed'}
    'temporalModulationParams', struct(...          % temporal: modulation params struct
      'stimOnFrameIndices', [1], ...                %   params relevant to the temporalModulationMode
      'stimDurationFramesNum', 1), ...              %   params relevant to the temporalModulationMode
    'frameDurationSeconds', 3/16, ...               % temporal: frame duration, in seconds
    'pupilDiameterMm', pupilDiameterMm ...          % pupil diameter mm
    );
                     
% Create a static two-spot AO scene with a particular incr-decr direction,
% and other relevant parameters. This uses compute function handle for
% two-spot AO stimuli, as commented above the parameters have been cooked
% to make one rectangular spot on a dark background.
rectComputeFunction = @sceAOTwoSpot;

% Instantiate a sceneEngine with the above sceneComputeFunctionHandle
% and the custom params.
rectScene = sceneEngine(rectComputeFunction, rectParams);

% Create a scene sequence at a particular contrast.  Here there is just one
% frame in the returned cell array of scenes.
visualizationContrast = 1.0;
rectSceneSequence = rectScene.compute(bgFactor*visualizationContrast);

% Visualize
ieAddObject(rectSceneSequence{1});
sceneWindow;
