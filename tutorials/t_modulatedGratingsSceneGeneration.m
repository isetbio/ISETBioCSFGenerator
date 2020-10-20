function t_modulatedGratingsSceneGeneration
% Demonstrates how to use the @sceneEngine class with the sceGrating scene compute function to generate a variety of grating stimuli
%
% Syntax:
%    t_modulatedGratingsSceneGeneration
%
% Description:
%   Demonstrates how to use the @sceneEngine class with the sceGrating scene 
%   compute function to generate a variety of grating stimuli that are commonly 
%   employed in neurophysiological and  psychophysical experiments. The 
%   sceGrating compute function enables extensive manipulation of stimulus 
%   chromatic, spatial and temporal parameters. 
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
%    10/17/2020  dhb  Remove pause commands so we can autotest with this.

    % Close figs
    close all;

    % Compute function handle for grating stimuli
    sceneComputeFunction = @sceGrating;
    
    % Retrieve the default params for the grating stimulus
    defaultGratingParams = sceGrating();

    % Specify a test with 100% of the achievable contrast
    testContrast = 1.0;
    
    % STIMULUS #1
    % A flashed grating sequence in which the grating is flashed ON
    % 3 times during the stimulus duration, each time for a 10 frames, with
    % an inter-flash duration of 20 frames, or 200 msec, since each frame is
    % lasting for 10 msec. Make it an L-only, 4 c/deg grating, with a 45
    % deg orientation, and a circular aperture with a radius of 0.3 degs
    
    % Start with the default params
    customGratingParams = defaultGratingParams;
    
    % Configure an L-M grating
    customGratingParams.coneContrastModulation = [0.09 -0.09 0.0];
    
    % Configure a 3 c/deg grating, with a 90 deg orientation, and a cosine spatial phase
    customGratingParams.spatialPhaseDegs = 0;
    customGratingParams.spatialFrequencyCyclesPerDeg = 3.0;
    customGratingParams.orientationDegs = 90;
    
    % Configure a disk spatial envelope
    customGratingParams.spatialEnvelope = 'disk';
    customGratingParams.minPixelsNumPerCycle = 30;
    customGratingParams.spatialEnvelopeRadiusDegs = 0.4;
    
    % Configure a 100 msec frame duration
    customGratingParams.frameDurationSeconds = 100/1000;
    
    % Configure temporal modulation
    customGratingParams.temporalModulation = 'flashed';
    customGratingParams.temporalModulationParams =  struct(...                      
            'stimOnFrameIndices', 1+[0:9 30:39 60:69], ...          
            'stimDurationFramesNum', 70);
            
    % Instantiate a sceneEngine with the above sceneComputeFunctionHandle 
    % and the custom grating params.
    theSceneEngine = sceneEngine(sceneComputeFunction, customGratingParams);
    
    % Compute the scene sequence
    [theSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(testContrast);
    
    % Visualize the generated scene sequence
    theSceneEngine.visualizeSceneSequence(theSceneSequence, theSceneTemporalSupportSeconds);
 
    % STIMULUS #2 
    % An L+M Gabor grating with a 60 deg orientation drifting at 8 Hz
    customGratingParams = defaultGratingParams;
    customGratingParams.coneContrastModulation = [0.6 0.6 0.0];
    customGratingParams.orientationDegs = 60;
    customGratingParams.temporalModulation = 'drifted';
    customGratingParams.temporalModulationParams =  struct(...                        
            'temporalFrequencyHz', 4, ...           
            'stimDurationTemporalCycles', 5);
        
    % Re-instantiate the sceneEngine with the customGratingParams
    theSceneEngine = sceneEngine(sceneComputeFunction, customGratingParams);
    
    % Compute the scene sequence
    [theSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(testContrast);
    
    % Visualize the generated scene sequence
    theSceneEngine.visualizeSceneSequence(theSceneSequence, theSceneTemporalSupportSeconds);
       
    % STIMULUS #3 
    % An achromatic grating with an orientation of 30 deg,
    % square aperture and a spatial phase of 0 degs, positioned at an 
    % off-center location, counterphased at 2 Hz
    %
    customGratingParams = defaultGratingParams;
    customGratingParams.coneContrastModulation = [0.6 0.6 0.6];
    customGratingParams.spatialPositionDegs =  [0.08 -0.07];
    customGratingParams.orientationDegs = 30;
    customGratingParams.spatialEnvelope = 'square';
    customGratingParams.spatialModulation =  'harmonic';
    customGratingParams.minPixelsNumPerCycle = 100;
    customGratingParams.spatialPhaseDegs = 0;
    customGratingParams.spatialEnvelopeRadiusDegs = 0.25;
    customGratingParams.temporalModulation = 'counter phase modulated';
    customGratingParams.temporalModulationParams =  struct(...                 
            'temporalFrequencyHz', 2, ...           
            'stimDurationTemporalCycles', 4);
        
    % Re-instantiate the sceneEngine with the customGratingParams
    theSceneEngine = sceneEngine(sceneComputeFunction, customGratingParams);
    
    % Compute the scene sequence
    [theSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(testContrast);
    
    % Visualize the generated scene sequence
    theSceneEngine.visualizeSceneSequence(theSceneSequence, theSceneTemporalSupportSeconds);
    
    
    % STIMULUS #4 
    % A M-cone isolating disk, pulsing at 4 Hz
    %
    customGratingParams = defaultGratingParams;
    customGratingParams.coneContrastModulation = [0.0 0.1 0.0];
    customGratingParams.spatialFrequencyCyclesPerDeg = 0;
    customGratingParams.spatialPositionDegs =  [-0.1 0.05];
    customGratingParams.spatialEnvelope = 'disk';
    customGratingParams.spatialModulation =  'square';
    customGratingParams.minPixelsNumPerCycle = 30;
    customGratingParams.spatialPhaseDegs = 0;
    customGratingParams.spatialEnvelopeRadiusDegs = 0.2;
    customGratingParams.temporalModulation = 'counter phase modulated';
    customGratingParams.temporalModulationParams =  struct(...                        
            'temporalFrequencyHz', 4, ...           
            'stimDurationTemporalCycles', 5);
        
    % Re-instantiate the sceneEngine with the customGratingParams
    theSceneEngine = sceneEngine(sceneComputeFunction, customGratingParams);
    
    % Compute the scene sequence
    [theSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(testContrast);
    
    % Visualize the generated scene sequence
    theSceneEngine.visualizeSceneSequence(theSceneSequence, theSceneTemporalSupportSeconds);
    
    % STIMULUS #5 
    % A L+M-S isolating polar sinusoid grating, drifting at 2 Hz
    %
    customGratingParams = defaultGratingParams;
    customGratingParams.coneContrastModulation = [0.7 0.7 -0.7];
    customGratingParams.spatialFrequencyCyclesPerDeg = 4;
    customGratingParams.spatialPositionDegs =  [00.1 -0.05];
    customGratingParams.spatialEnvelope = 'Gaussian';
    customGratingParams.spatialModulationDomain = 'polar';
    customGratingParams.spatialEnvelopeRadiusDegs = 0.2;
    customGratingParams.temporalModulation = 'drifted';
    customGratingParams.temporalModulationParams =  struct(...                        
            'temporalFrequencyHz', 2, ...           
            'stimDurationTemporalCycles', 4);
        
    % Re-instantiate the sceneEngine with the customGratingParams
    theSceneEngine = sceneEngine(sceneComputeFunction, customGratingParams);
    
    % Compute the scene sequence
    [theSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(testContrast);
    
    % Visualize the generated scene sequence
    theSceneEngine.visualizeSceneSequence(theSceneSequence, theSceneTemporalSupportSeconds);
    
end


