
% Create a static grating scene with specific chromatic direction,
% spatial frequency, and duration

function [gratingScene] = createGratingScene(varargin)
%createGratingScene  Create a static grating scene
%
% Usage:
%   createGratingScene('chromaDir', chromaDir, 'spatialFreq', spatialFreq);
%   See t_spatialCSF.m
%
% Inputs:
%   None.
%
% Outputs:
%   sceneEngine Object.
%
% Optional key/value pairs:
%   'chromaDir'       - 1-by-3 vector specifying the chromatic direction
%                     of the stimulus as L, M, S cone contrast
%
%   'spatialFreq'     - The spatial frequency of the stimulus, in cyc/deg
%
%   'spatialPhase'    - The spatial phase of the stimulus, in degree
%
%   'orientation'     - The orientation of the stimulus, in degree
%
%   'duration'        - The duration of the stimulus, in second

p = inputParser;

p.addParameter('chromaDir', [0.1, 0.1, 0.1], @(x)(isnumeric(x) && numel(x) == 3));
p.addParameter('spatialFreq', 1, @(x)(isnumeric(x) && numel(x) == 1));
p.addParameter('spatialPhase', 0, @(x)(isnumeric(x) && numel(x) == 1));
p.addParameter('orientation', 90, @(x)(isnumeric(x) && numel(x) == 1));
p.addParameter('duration', 0.1, @(x)(isnumeric(x) && numel(x) == 1));

parse(p, varargin{:});

% Compute function handle for grating stimuli
sceneComputeFunction = @sceGrating;

% Retrieve the default params for the grating stimulus
gratingParams = sceGrating();

% Configure chromatic direction and and spatial frequency of the grating
% with a 90 deg orientation, and a cosine spatial phase
gratingParams.coneContrastModulation = p.Results.chromaDir;
gratingParams.spatialFrequencyCyclesPerDeg = p.Results.spatialFreq;
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