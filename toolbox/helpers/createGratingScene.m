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
%   'spatialPhase'      - The spatial phase of the stimulus, in degrees
%                         Default 0.
%   'orientation'       - The orientation of the stimulus, in degrees.
%                         Default 90 (vertical grating).
%   'duration'          - The duration of the stimulus, in seconds.
%                         Default: 0.1.
%   'presentationMode'  - Presentation mode, for now either 'flashed' (1
%                         frame), or 'sampled motion', (4 frames, with
%                         spatial phase advancing by 90 degs in each frame)
%
% See also: t_modulatedGratingsSceneGenerateion, t_spatialCSF
%

% History:
%   10/23/20  dhb  Added comments.
%   11/07/20  npc  Added 'presentationMode' key/value pair   

% Set up parameters with defaults
p = inputParser;
p.addParameter('spatialPhase', 0, @(x)(isnumeric(x) && numel(x) == 1));
p.addParameter('spatialEnvelope', 'disk', @(x)(ischar(x) && ismember(x, {'disk', 'square', 'Gaussian'})));
p.addParameter('orientation', 90, @(x)(isnumeric(x) && numel(x) == 1));
p.addParameter('duration', 0.1, @(x)(isnumeric(x) && numel(x) == 1));
p.addParameter('presentationMode', 'flashed', @(x)(ischar(x) && ismember(x,{'flashed', 'sampled motion'})));
p.addParameter('pixelsNum', 0, @(x)(isnumeric(x) && numel(x) == 1));
p.addParameter('fovDegs', 1.0, @(x)(isnumeric(x) && numel(x) == 1));
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
gratingParams.fovDegs = p.Results.fovDegs;

% Set pixel size
pixelsNum = p.Results.pixelsNum;
if(pixelsNum >= 1)
    gratingParams.pixelsNum = pixelsNum;
end

% Configure a disk spatial envelope
gratingParams.spatialEnvelope = p.Results.spatialEnvelope;
gratingParams.minPixelsNumPerCycle = 30;
gratingParams.spatialEnvelopeRadiusDegs = 0.4;

% Configure temporal modulation:
switch (p.Results.presentationMode)
    case 'flashed'
        % Single frame presentation
        gratingParams.frameDurationSeconds = p.Results.duration;
        gratingParams.temporalModulation = 'flashed';
        gratingParams.temporalModulationParams =  struct(...
            'stimOnFrameIndices', 1, 'stimDurationFramesNum', 1);
    case 'sampled motion'
        % N-frame presentation, each frame advancing spatial phase by 360/N degs
        % with the N-frames lasting for p.Results.duration. Here N = 4, so
        % each frame is advancing spatial phase by 90 degs.
        spatialPhaseSamplesNum = 4;
        temporalFrequencyHz = 1.0/p.Results.duration;
        gratingParams.frameDurationSeconds = p.Results.duration/spatialPhaseSamplesNum;
        gratingParams.temporalModulation = 'drifted';
        gratingParams.temporalModulationParams =  struct(...
            'temporalFrequencyHz', temporalFrequencyHz, ...
            'stimDurationTemporalCycles', 1);
    otherwise
        error('Unknown presentationMode: ''%s''.', p.Results.presentationMode);
        
end

% Instantiate a sceneEngine with the above sceneComputeFunctionHandle
% and the custom grating params.
gratingScene = sceneEngine(sceneComputeFunction, gratingParams);

end