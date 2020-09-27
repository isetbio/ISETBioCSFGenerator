function t_sceneGeneration
% How the @sceneEngine object to generate a simple scene sequence
%
% Syntax:
%    t_sceneGeneration
%
% Description:
%    Demonstrates how to generate a stimulus sequence using a
%    @sceneEngine object and a scene compute function. The scene computed by
%    the 'uniformFieldTemporalModulation' compute function is a uniform field whose
%    luminance is stepped up during some interval. Here we are passing a
%    custom scene params struct during instantiation of the @sceneEngine object.
%    If no scene params struct were passed, the default scene params 
%    defined in the 'uniformFieldTemporalModulation' compute function would be used.
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
%
% See Also:
%   t_neuralResponseCompute

% History:
%    09/21/2020  NPC  Wrote it.

    % Configure the function handle and the params for the @sceneEngine
    % User supplied compute function
    sceneComputeFunction = @uniformFieldTemporalModulation;
    % User supplied struct with params appropriate for the @sceneEngine sceneComputeFunction
    customSceneParams = struct(...
        'fovDegs', 0.25, ...                        % 0.25 degs across
        'meanLuminanceCdPerM2', 100, ...            % 100 cd/m2 mean luminance
        'frameDurationSeconds', 50/1000, ...        % 50 msec frame duration
        'stimDurationFramesNum', 4, ...             % total time: 200 msec
        'stimOnsetFramesIndices', [2 3], ...        % modulate luminance at frames 1 and 2, so between 50 and 150 msec
        'sizePixels', 64 ...                        % 64 x 64 pixels
    );

    % Instantiate a sceneEngine with the above sceneComputeFunctionHandle 
    % and the custom scene params
    theSceneEngine = sceneEngine(sceneComputeFunction, customSceneParams);
    
    % Specify a pedestal luminance with 70% contrast
    testContrast = 0.7;

    % Compute the scene sequence
    [theSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(testContrast);
    
    % Visualize the generated scene sequence
    debugSceneGeneration = true;
    if (debugSceneGeneration)
        renderSceneSequence(theSceneSequence, theSceneTemporalSupportSeconds);
    end

end

function renderSceneSequence(sceneSequence, temporalSupportSeconds)
    figure(1); clf;
    scenesNum = numel(sceneSequence);
    for frameIndex = 1:numel(sceneSequence)
        subplot(1, scenesNum, frameIndex);
        xyzImage = sceneGet(sceneSequence{frameIndex}, 'xyz');
        luminanceMap = squeeze(xyzImage(:,:,2));
        imagesc(luminanceMap);
        axis 'image';
        set(gca, 'CLim', [0 200]);
        title(sprintf('frame %d\n(%2.0f msec)', frameIndex,temporalSupportSeconds(frameIndex)*1000));
        colormap(gray);
    end
    
end

