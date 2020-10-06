function t_modulatedGratingsSceneGeneration
% Demonstrates how to use the @sceneEngine object to generate a variety of grating stimuli
%
% Syntax:
%    t_sceneGeneration
%
% Description:
%   Demonstrates how to use the @sceneEngine object to generate a variety of grating stimuli
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

    % Configure the function handle and the params for the @sceneEngine
    % user supplied compute function
    sceneComputeFunction = @sceGrating;
    
    % Retrieve the default params for the grating stimulus
    defaultGratingParams = sceGrating();

    customGratingParams = defaultGratingParams;
    customGratingParams.temporalModulationParams =  struct(... 
            'mode', 'flashed', ...                      
            'stimOnFrameIndices', [1:10 20:30 40:50 60:80], ...          
            'stimDurationFramesNum', 100);
            
    % Instantiate a sceneEngine with the above sceneComputeFunctionHandle 
    % and the default grating scene params.
    theSceneEngine = sceneEngine(sceneComputeFunction, customGratingParams);
    
    % Specify a test with 100% contrast
    testContrast = 1.0;

    % Compute the scene sequence
    [theSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(testContrast);
    
    % Visualize the generated scene sequence
    theSceneEngine.visualizeSceneSequence(theSceneSequence, theSceneTemporalSupportSeconds);
    pause

    % Params for a drifting grating at 8 Hz
    customGratingParams = defaultGratingParams;
    customGratingParams.temporalModulationParams =  struct(... 
            'mode', 'drifted', ...                       
            'temporalFrequencyHz', 4, ...           
            'stimDurationTemporalCycles', 4);
        
    % Re-instantiate the sceneEngine with the customGratingParams
    theSceneEngine = sceneEngine(sceneComputeFunction, customGratingParams);
    
    % Compute the scene sequence
    [theSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(testContrast);
    
    % Visualize the generated scene sequence
    theSceneEngine.visualizeSceneSequence(theSceneSequence, theSceneTemporalSupportSeconds);
    pause
    
    
    % Params for a counterphased grating at 4 Hz
    customGratingParams.temporalModulationParams =  struct(... 
            'mode', 'counter phase modulated', ...                       
            'temporalFrequencyHz', 2, ...           
            'stimDurationTemporalCycles', 4);
        
    % Re-instantiate the sceneEngine with the customGratingParams
    theSceneEngine = sceneEngine(sceneComputeFunction, customGratingParams);
    
    % Compute the scene sequence
    [theSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(testContrast);
    
    % Visualize the generated scene sequence
    theSceneEngine.visualizeSceneSequence(theSceneSequence, theSceneTemporalSupportSeconds);
end


