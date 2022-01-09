% t_AORectangleScene
%
% Description:
%   Use the machinery of sceTwoSpot to produce an AO apparatus scene with a
%   rectanglur stimuls.

%% Set parameters 
defocusDiopters = 0.05;
pupilDiameterMm = 7;
 
% Rectangle horizontal and vertical size in arcminutes.
spotHeightMinutes = 2.5;
spotWidthMinutes= 2.5;

% Position of rectangle. These numbers give the center position in arcmin.
%    Vertical - positive is up 
%    Horizontal - positive is right
spotVerticalPositionMinutes = 2.5;
spotHorizontalPositionMinutes = 5;

% Other stimulus parameters
%
% Define basic parameters of the AO stimulus
wls = 400:10:750;
spotWavelengthNm = 580;

nPixels = 128;
fieldSizeMinutes = 60;
fieldSizeDegs = fieldSizeMinutes/60;

% Measured power 550 nm visible for an AO experiment at Penn was 76 nW for
% 1.5 by 1 degree field. This gives us a ballpark reasonable number.
%
% We don't bother to model IR imaging light effect, it's swamped by the
% visible.
spotPowerUW = (fieldSizeDegs^2)*76/1000*(2/3);

% The underlying routine starts with a background and works in contrast.
% To deal with this, we specify a background down from desired power by
% some large factor and then take contrast as this factor.
bgFactor = 5000;
spotBgPowerUW = spotPowerUW/bgFactor;

% Convert our parameters to those that matter to the two spot routine.
% Using stimAngle of 90 means we only see spot 2, which is drawn second and
% thus always over spot 1.
stimAngle = 90;

% This is center-to-center separation of the two spots. So with values of
% zero they overlap. Because the separation is accomplished by applying
% half of it to each underlying spot, we multiply by two to convert
% position offsets to separations.  The minus sign makes the sign
% convention as described above.
spotVerticalSepMinutes = -2*spotVerticalPositionMinutes;
spotHorizontalSepMinutes = -2*spotHorizontalPositionMinutes;

% Set up two spot params
rectParams = struct(...
    'type', 'basic', ...                            % type
    'viewingDistanceMeters', 2, ...                 % viewing distance: in meters
    'wls', 400:10:750, ...                          % wavelength support of the primary SPDs: in nanometers
    'stimAngle', stimAngle, ...                     % stimulus angle in incr/decr plane
    'spotWl', spotWavelengthNm, ...                 % spot wavelength: in nm
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
% and other relevant parameters
% Compute function handle for two-spot AO stimuli
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
