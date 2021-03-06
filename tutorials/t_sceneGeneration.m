function t_sceneGeneration
% Demonstrates how to use the @sceneEngine object to generate a simple scene sequence
%
% Syntax:
%    t_sceneGeneration
%
% Description:
%    Demonstrates how to generate a stimulus sequence using a
%    @sceneEngine object and a scene compute function. The scene computed by
%    the 'sceUniformFieldTemporalModulation' compute function is a uniform field whose
%    luminance is stepped up during some interval. Here we are passing a
%    custom scene params struct during instantiation of the @sceneEngine object.
%    If no scene params struct were passed, the default scene params 
%    defined in the 'sceUniformFieldTemporalModulation' compute function would be used.
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
%   t_neuralResponseCompute, t_thresholdEngine, t_neuralResponseCompute,
%   t_responseClassifier.
%

% History:
%    09/21/2020  NPC  Wrote it.

    % Close figs
    close all;

    % Configure the function handle and the params for the @sceneEngine
    % user supplied compute function.  This example function is provided with 
    % the toolbox.  See its help text for a description of the API for the
    % compute function.
    sceneComputeFunction = @sceUniformFieldTemporalModulation;
    
    % User supplied struct with params appropriate for the
    % sceneComputeFunction.  If you didn't know what fields were there,
    % you could execute the compute function with no input arguments. This will
    % return the default parameter struct for this function. That is,
    % enter:
    %   sceUniformFieldTemporalModulation
    customSceneParams = struct(...
        'fovDegs', 0.25, ...                        % 0.25 degs across
        'meanLuminanceCdPerM2', 200, ...            % 200 cd/m2 mean luminance
        'frameDurationSeconds', 50/1000, ...        % 50 msec frame duration
        'stimDurationFramesNum', 4, ...             % total time: 200 msec
        'stimOnsetFramesIndices', [2 3], ...        % modulate luminance at frames 1 and 2, so between 50 and 150 msec
        'sizePixels', 64 ...                        % 64 x 64 pixels
    );

    % Instantiate a sceneEngine with the above sceneComputeFunctionHandle 
    % and the custom scene params.  In the general case, the
    % customSceneParams structure should bewhat the sceneComputeFunction
    % expects to be passed.  If you instantiate without a parameters
    % structure, the default parameter structure is used.
    theSceneEngine = sceneEngine(sceneComputeFunction, customSceneParams);
    
    % Specify a test with 70% contrast
    testContrast = 0.7;

    % Compute the scene sequence
    [theSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(testContrast);
    
    % Visualize the generated scene sequence
    debugSceneGeneration = true;
    if (debugSceneGeneration)
        theSceneEngine.visualizeSceneSequence(theSceneSequence, theSceneTemporalSupportSeconds);
    end
    
    % Set a different parameters structure and generate a new sequence.
    % Here we lower the background luminance so everything looks darker,
    % and change the temporal sequence of the modulation. 
    customSceneParams.meanLuminanceCdPerM2 = 90;
    customSceneParams.stimOnsetFramesIndices = [1 4];
    
    % Note that we need to re-instantiate the @sceneEngine object to do
    % this, since the convention is that a given @sceneEngine object only
    % generates one scene, perhaps with a different contrast.
    theSceneEngine = sceneEngine(sceneComputeFunction, customSceneParams);
    [theSceneSequence1, theSceneTemporalSupportSeconds1] = theSceneEngine.compute(testContrast);
    if (debugSceneGeneration)
        theSceneEngine.visualizeSceneSequence(theSceneSequence1, theSceneTemporalSupportSeconds1);
    end
    
    % As noted above, if you didn't know what fields the parameters struct
    % the sceUniformFieldTemporalModulation function expects, you can get a
    % default parameters structure by calling
    % sceUniformFieldTemporalModulation without any arguments
    defaultSceneParams = sceUniformFieldTemporalModulation
    
    % Also as noted above, you can instantiate a scene generation object with the default
    % parameters by not providing them:
    theSceneEngine = sceneEngine(sceneComputeFunction);
    theSceneEngine.sceneParams
end


