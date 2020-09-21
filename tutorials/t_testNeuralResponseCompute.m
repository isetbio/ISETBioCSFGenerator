function t_testNeuralResponseCompute
% Use the @neuralResponseEngine object to compute cone mosaic isomerization responses
%
% Syntax:
%    t_testNeuralResponseCompute
%
% Description:
%    Demonstrates how to generate a stimulus sequence using a
%    @sceneGeneration object and how to compute cone mosaic excitation 
%    resoponses used a @neuralResponseEngine object. 
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

    % Configure the function handle and the params for the @neuralResponseEngine
    % This is a function that the USER has to supply
    neuralComputeFunction = @photopigmentExcitationsWithNoEyeMovements;

    % Struct that the USER has to supply which defines the compute pipeline
    opticsParams = struct(...
        'type', 'wvf human', ...
        'pupilDiameterMM', 3.0 ...
    );

    coneMosaicParams = struct(...
        'upsampleFactor', 5, ...
        'fovDegs', 0.3, ...
        'timeIntegrationSeconds', 5/1000 ...
    );

    neuralResponseParams = struct(...
        'opticsParams', opticsParams, ...
        'coneMosaicParams', coneMosaicParams ...
    );

    % Instantiate a sceneGenerationEngine with the above sceneComputeFunctionHandle and sceneParams
    theSceneEngine = sceneGenerationEngine(sceneComputeFunction, sceneParams);
    
    % Instantiate a neuralResponseEngine with the desired neural response params
    theNeuralEngine = neuralResponseEngine(neuralComputeFunction, neuralResponseParams);
    
    % Specify a pedestal luminance with 70% contrast
    testContrast = 0.7;

    % Compute the scene sequence
    [theSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(testContrast);

    % Compute 8 instances of neural responses to the input scene sequence
    instancesNum = 8;
    [theIsomerizationsSequence, theIsomerizationsTemporalSupportSeconds] = ...
        theNeuralEngine.compute(theSceneSequence, theSceneTemporalSupportSeconds, instancesNum);

    % Visualize result
    debugNeuralResponseGeneration = true;
    if (debugNeuralResponseGeneration)
        renderNeuralResponseSequence(theIsomerizationsSequence, theIsomerizationsTemporalSupportSeconds);
    end
    
end

function renderNeuralResponseSequence(theResponseSequence, theResponseTemporalSupportSeconds)
    figure(2);
    meanResponse = squeeze(mean(theResponseSequence,1));
    cellsNum = size(meanResponse ,1);
    imagesc(theResponseTemporalSupportSeconds, 1:cellsNum, meanResponse);
    xlabel('time (sec)');
    ylabel('cells');
    colormap(gray);
end
