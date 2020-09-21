function t_testSceneGeneration
% How to generate a simple scene sequence using @sceneGenerationEngine objects
%
% Syntax:
%    t_testSceneGeneration
%
% Description:
%    Demonstrates how to generate a stimulus sequence using a
%    @sceneGeneration object. The scene represents a uniform field whose
%    luminance is stepped up during some interval.
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
%   t_testNeuralResponseCompute

% History:
%    09/21/2020  NPC  Wrote it.

    % Configure the function handle and the params for the @sceneGenerationEngine
    % This is a function that the USER has to supply
    sceneComputeFunction = @uniformFieldTemporalModulation;
    % This is a struct that the USER has to supply and which is to work
    % with the user-supplied function handle
    sceneParams = struct(...
        'fovDegs', 0.25, ...                        % 0.25 degs across
        'meanLuminanceCdPerM2', 100, ...            % 100 cd/m2 mean luminance
        'frameDurationSeconds', 50/1000, ...        % 50 msec frame duration
        'stimDurationFramesNum', 4, ...             % total time: 200 msec
        'stimOnsetFramesIndices', [2 3], ...        % modulate luminance at frames 1 and 2, so between 50 and 150 msec
        'sizePixels', 64 ...                        % 64 x 64 pixels
    );

    % Instantiate a sceneGenerationEngine with the above sceneComputeFunctionHandle and sceneParams
    theSceneEngine = sceneGenerationEngine(sceneComputeFunction, sceneParams);
    
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

