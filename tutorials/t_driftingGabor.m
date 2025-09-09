function t_driftingGabor()
% Demonstrates how to generate a stimulus sequence representing a drifting Gabor
%
% Syntax:
%    t_sceneGeneration
%
% Description:
%    Demonstrates how to generate a stimulus sequence representing a
%    drifting Gabor
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
%    05/12/2023  NPC  Wrote it.

    % Close figs
    close all;

    sceneComputeFunction = @sceGrating;

    customSceneParams = sceneComputeFunction();
    %customSceneParams.temporalModulationParams.stimDurationFramesNum = 6;

    customSceneParams.coneContrastModulation = [0.7 0.7 0.7];
    customSceneParams.spatialFrequencyCyclesPerDeg = 6.0;
    customSceneParams.orientationDegs = 90;
    customSceneParams.spatialEnvelopeRadiusDegs = 0.3;
    customSceneParams.temporalModulation = 'counterphasemodulated';
    customSceneParams.temporalModulation = 'drifted';

    customSceneParams.frameDurationSeconds = 0.1;
    customSceneParams.temporalModulationParams.stimOnFrameIndices = [];
    customSceneParams.temporalModulationParams.stimDurationFramesNum = [];
    customSceneParams.temporalModulationParams.stimDurationTemporalCycles = 2.0;
    customSceneParams.temporalModulationParams.temporalFrequencyHz = 1.0;

    customSceneParams.temporalModulationParams

    theSceneEngine = sceneEngine(sceneComputeFunction, customSceneParams);
    testContrast = 1.0;
    [theSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(testContrast);
    figure(1); clf;

    % Determine max luminance
    m2 = [];
    for iFrame = 1:(numel(theSceneSequence))
        lumMap = sceneGet(theSceneSequence{iFrame}, 'luminance');
        m2(iFrame) = max(max(lumMap(:)));
    end
    maxLuminance = max(m2);

    for iFrame = 1:(numel(theSceneSequence))
        imagesc(sceneGet(theSceneSequence{iFrame}, 'luminance'));
        set(gca, 'XTick', [], 'YTick', [], 'CLim', [0 maxLuminance]);
        colormap(gray)
        axis 'image'
        title(sprintf('t:%2.1f milliseconds', theSceneTemporalSupportSeconds(iFrame)*1000))
        drawnow
        pause(0.1)
    end

end
