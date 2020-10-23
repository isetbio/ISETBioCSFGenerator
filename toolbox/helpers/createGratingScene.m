function [gratingScene] = createGratingScene(chromaticDir, spatialFrequency, varargin)
% Create a grating scene using sceGrating
%
% Syntax:
%   [gratingScene] = createGratingScene(chromaticDir, spatialFreqency)
%
% Description:
%    Create a @sceneEngine object representing a chromatic grating scene.  Key/value pairs
%    control parameters.
%
% Inputs:
%   chromaticDir      - Three dimensional vector giving L, M, and S
%                       contrasts for the grating at 100% contrast.
%   spatialFrequency  - Grating spatial frequency.
%
% Outputs:
%   gratingScene      - Created sceneEngine
%
% Optional key/value pairs:
%   'spatialPhase'    - The spatial phase of the stimulus, in degrees
%                       Default 0.
%   'orientation'     - The orientation of the stimulus, in degrees.
%                       Default 90 (vertical grating).
%   'duration'        - The duration of the stimulus, in seconds.  Default
%                       0.1.
%
% See also: t_modulatedGratingsSceneGenerateion, t_spatialCSF
%

% History:
%   10/23/20  dhb  Added comments.

% Set up parameters with defaults
p = inputParser;
p.addParameter('spatialPhase', 0, @(x)(isnumeric(x) && numel(x) == 1));
p.addParameter('orientation', 90, @(x)(isnumeric(x) && numel(x) == 1));
p.addParameter('duration', 0.1, @(x)(isnumeric(x) && numel(x) == 1));
parse(p, varargin{:});

% Compute function handle for grating stimuli
sceneComputeFunction = @sceGrating;

% Retrieve the default params for the grating stimulus
gratingParams = sceGrating();

% Configure those parameters that we adjust through key/value pairs.
% chromatic direction and and spatial frequency of the grating
% with a 90 deg orientation, and a cosine spatial phase
gratingParams.coneContrastModulation = chromaticDir;
gratingParams.spatialFrequencyCyclesPerDeg = spatialFrequency;
gratingParams.spatialPhaseDegs = p.Results.spatialPhase;
gratingParams.orientationDegs = p.Results.orientation;

% Configure a disk spatial envelope
gratingParams.spatialEnvelope = 'disk';
gratingParams.minPixelsNumPerCycle = 30;
gratingParams.spatialEnvelopeRadiusDegs = 0.4;

% Configure temporal modulation: 100 ms duration for only 1 frame
gratingParams.frameDurationSeconds = p.Results.duration;
gratingParams.temporalModulation = 'flashed';
gratingParams.temporalModulationParams =  struct(...
    'stimOnFrameIndices', 1, 'stimDurationFramesNum', 1);

% Instantiate a sceneEngine with the above sceneComputeFunctionHandle
% and the custom grating params.
gratingScene = sceneEngine(sceneComputeFunction, gratingParams);

end