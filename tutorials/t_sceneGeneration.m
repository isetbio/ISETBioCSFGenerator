function t_sceneGeneration
% How the @sceneEngine object to generate a simple scene sequence
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
%   t_neuralResponseCompute

% History:
%    09/21/2020  NPC  Wrote it.

    % Close figs
    close all;

    % Configure the function handle and the params for the @sceneEngine
    % user supplied compute function
    sceneComputeFunction = @sceUniformFieldTemporalModulation;
    
    % User supplied struct with params appropriate for the sceneComputeFunction
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
    % customSceneParams structure should what the sceneComputeFunction
    % expects to be passed.
    theSceneEngine = sceneEngine(sceneComputeFunction, customSceneParams);
    
    % Specify a test with 70% contrast
    testContrast = 0.7;

    % Compute the scene sequence
    [theSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(testContrast);
    
    % Visualize the generated scene sequence
    debugSceneGeneration = true;
    if (debugSceneGeneration)
        renderSceneSequence(theSceneSequence, theSceneTemporalSupportSeconds);
    end
    
    % Set a different parameters structure and generate a new sequence.
    % Here we lower the background luminance so everything looks darker,
    % and change the temporal sequence of the modulation. 
    customSceneParams.meanLuminanceCdPerM2 = 90;
    customSceneParams.stimOnsetFramesIndices = [1 4];
    
    % Note that we need to re-instantiate the @sceneEngine object to do this, since the 
    % convention is that a given @sceneEngine object only generates one scene, perhaps
    % with a different contrast.
    theSceneEngine = sceneEngine(sceneComputeFunction, customSceneParams);
    
    [theSceneSequence1, theSceneTemporalSupportSeconds1] = theSceneEngine.compute(testContrast);
    if (debugSceneGeneration)
        renderSceneSequence(theSceneSequence1, theSceneTemporalSupportSeconds1);
    end
    
    % If you didn't know what fields the parameters struct the
    % sceUniformFieldTemporalModulation function expects, you can get a
    % default parameters structure by calling
    % sceUniformFieldTemporalModulation without any arguments
    defaultSceneParams = sceUniformFieldTemporalModulation
    
    % You can also instantiate a scene generation object with the default
    % parameters by not providing them:
    theSceneEngine = sceneEngine(sceneComputeFunction);
    theSceneEngine.sceneParams

end

function renderSceneSequence(sceneSequence, temporalSupportSeconds)
    figure; clf;
    scenesNum = numel(sceneSequence);
    for frameIndex = 1:numel(sceneSequence)
        subplot(1, scenesNum, frameIndex);
        xyzImage = sceneGet(sceneSequence{frameIndex}, 'xyz');
        luminanceMap = squeeze(xyzImage(:,:,2));
        image(luminanceMap);
        axis 'image';
        set(gca, 'CLim', [0 500]);
        title(sprintf('frame %d\n(%2.0f msec)', frameIndex,temporalSupportSeconds(frameIndex)*1000));
        colormap(gray);
    end
    
end

