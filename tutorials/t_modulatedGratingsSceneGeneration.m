function t_modulatedGratingsSceneGeneration
% Demonstrates how to use the @sceneEngine object with the @sceGrating compute function to generate a variety of grating stimuli
%
% Syntax:
%    t_sceneGeneration
%
% Description:
%   Demonstrates how to use the @sceneEngine object with the @sceGrating scene compute function 
%   to generate a variety of grating stimuli that are commonly employed in neurophysiological and 
%   psychophysical experiments.
%
% Inputs:
%    None.
%
% Outputs:
%   None.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%   t_sceneGeneration

% History:
%    10/06/2020  NPC  Wrote it.

    % Close figs
    close all;

    % Compute function handle for grating stimuli
    sceneComputeFunction = @sceGrating;
    
    % Retrieve the default params for the grating stimulus
    defaultGratingParams = sceGrating();

    % STIMULUS #1
    % A flashed grating sequence in which the grating is flashed ON
    % 3 times during the stimulus duration, each time for a 10 frames, with
    % an inter-flash duration of 20 frames, or 200 msec, since each frame is
    % lasting for 10 msec. Make it an L-only, 4 c/deg grating, with a 45
    % deg orientation, and a circular aperture with a radius of 0.3 degs
    
    % Start with the default params
    customGratingParams = defaultGratingParams;
    
    % Configure an L-only grating
    customGratingParams.coneContrastModulation = [0.1 0.0 0.0];
    
    % Configure a 4 c/deg grating, with a 45 deg orientation
    customGratingParams.spatialFrequencyCyclesPerDeg = 4.0;
    customGratingParams.orientationDegs = 45;
    
    % Configure a disk spatial envelope
    customGratingParams.spatialEnvelope = 'disk';
    customGratingParams.spatialEnvelopeRadiusDegs = 0.3;
    
    % Configure a 100 msec frame duration
    customGratingParams.frameDurationSeconds = 10/1000;
    
    % Configure temporal modulation
    customGratingParams.temporalModulationParams =  struct(... 
            'mode', 'flashed', ...                      
            'stimOnFrameIndices', 1+[0:9 30:39 60:69], ...          
            'stimDurationFramesNum', 80);
            
    % Instantiate a sceneEngine with the above sceneComputeFunctionHandle 
    % and the custom grating params.
    theSceneEngine = sceneEngine(sceneComputeFunction, customGratingParams);
    
    % Specify a test with 100% contrast
    testContrast = 1.0;

    % Compute the scene sequence
    [theSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(testContrast);
    
    % Visualize the generated scene sequence
    theSceneEngine.visualizeSceneSequence(theSceneSequence, theSceneTemporalSupportSeconds);
    disp('Hit enter to continue');
    pause

    
    % STIMULUS #2 
    % An  L+M grating drifting at 8 Hz
    customGratingParams = defaultGratingParams;
    customGratingParams.coneContrastModulation = [0.6 0.6 0.0];
    
    customGratingParams.temporalModulationParams =  struct(... 
            'mode', 'drifted', ...                       
            'temporalFrequencyHz', 4, ...           
            'stimDurationTemporalCycles', 5);
        
    % Re-instantiate the sceneEngine with the customGratingParams
    theSceneEngine = sceneEngine(sceneComputeFunction, customGratingParams);
    
    % Compute the scene sequence
    [theSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(testContrast);
    
    % Visualize the generated scene sequence
    theSceneEngine.visualizeSceneSequence(theSceneSequence, theSceneTemporalSupportSeconds);
    disp('Hit enter to continue');
    pause
    
    
    % STIMULUS #3 
    % An achromatic square-wave grating counterphased at 8 Hz
    % positioned at a spatially-offset location
    customGratingParams = defaultGratingParams;
    customGratingParams.coneContrastModulation = [0.6 0.6 0.6];
    customGratingParams.spatialPositionDegs =  [0.1 -0.1];
    customGratingParams.spatialEnvelope = 'square';
    customGratingParams.spatialModulation =  'square';
    customGratingParams.spatialPhaseDegs = 0;
    customGratingParams.spatialEnvelopeRadiusDegs = 0.3;
    customGratingParams.temporalModulationParams =  struct(... 
            'mode', 'counter phase modulated', ...                       
            'temporalFrequencyHz', 8, ...           
            'stimDurationTemporalCycles', 10);
        
    % Re-instantiate the sceneEngine with the customGratingParams
    theSceneEngine = sceneEngine(sceneComputeFunction, customGratingParams);
    
    % Compute the scene sequence
    [theSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(testContrast);
    
    % Visualize the generated scene sequence
    theSceneEngine.visualizeSceneSequence(theSceneSequence, theSceneTemporalSupportSeconds);
end


