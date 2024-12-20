function t_neuralResponseCompute
% Use the @neuralResponseEngine object to compute cone mosaic isomerization responses
%
% Syntax:
%    t_neuralResponseCompute
%
% Description:
%    Demonstrates how to generate a stimulus sequence using a
%    @sceneEngine object and a neural compute function.
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
% See also: t_thresholdEngine, t_sceneGeneration, t_responseClassifier
%

% History:
%    09/21/2020  NPC  Wrote it.

    % Close figures
    close all;

    % Instantiate the scene engine with a compute function.
    % This is a function that the USER has to supply.
    sceneComputeFunction = @sceUniformFieldTemporalModulation;
    sceneParams = sceUniformFieldTemporalModulation;
    theSceneEngine = sceneEngine(sceneComputeFunction,sceneParams);
    
    % Configure the function handles and the params for the @neuralResponseEngine
    % These are functions that the user has to specify, and write if they are not using
    % one of our provided ones.  One funciton computes the noise free
    % responses, while the other adds noise.
    %
    % We illustrate overriding some but not all of the default response
    % parameters in the code below.
    noiseFreeResponseParams = nreNoiseFreePhotopigmentExcitationsCmosaic;
    noiseFreeResponseParams.opticsParams.pupilDiameterMM = 2.0;
    noiseFreeResponseParams.coneMosaicParams.sizeDegs = [0.25 0.25];
    noiseFreeResponseParams.coneMosaicParams.timeIntegrationSeconds = sceneParams.frameDurationSeconds;;
    noisyInstancesParams = nreNoisyInstancesPoisson;
    theNeuralEngine = neuralResponseEngine( ...
        @nreNoiseFreePhotopigmentExcitationsCmosaic, ...
        @nreNoisyInstancesPoisson, ...
        noiseFreeResponseParams, ...
        noisyInstancesParams);
    
    % Specify a pedestal luminance with 70% contrast
    testContrast = 0.7;

    % Compute the scene sequence
    [theSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(testContrast);

    % Compute instances of neural responses to the input scene sequence.
    %
    % We specify the types of noise to be applied to the computed responses.
    % This has to be a cell array.  All nre compute functions need to understand 
    % 'none' and 'random'. Specific nre compute functions are allowed understand additional options
    % as appropriate to the specific model.
    %
    % It is possible to freeze the noise by specifying a seed for the
    % randome number generator through the 'rngSeed' key/value pair.  The
    % compute function should restore the rng to its current state if this
    % is passed, but not otherwise.
    instancesNum = 8;
    [noiseFreeResponse, theResponseTemporalSupportSeconds] = theNeuralEngine.computeNoiseFree(...
            theSceneSequence, ...
            theSceneTemporalSupportSeconds ...
            );
   [noisyInstances, ~] = theNeuralEngine.computeNoisyInstances(...
            noiseFreeResponse, ...
            theResponseTemporalSupportSeconds, ...
            instancesNum, ...
            'random' ...
            );

    assert(size(noisyInstances,1) == instancesNum);
    assert(size(noisyInstances,3) == length(theResponseTemporalSupportSeconds));
        
    % Visualize responses 
    debugNeuralResponseGeneration = true;
    if (debugNeuralResponseGeneration)
        renderNeuralResponseSequence(1, noiseFreeResponse, theResponseTemporalSupportSeconds,'noise free');
        renderNeuralResponseSequence(2, noisyInstances, theResponseTemporalSupportSeconds,'noisy instances');
    end
    
    % Generate 1 noisy instance with a specified rng seed.  We do this
    % a few times to verify that things work as they should
    instancesNum = 1;
    [noisyInstancesFixedSeed1, ~] = theNeuralEngine.computeNoisyInstances(...
        noiseFreeResponse, ...
        theResponseTemporalSupportSeconds, ...
        instancesNum, ...
        'random', ...
        'rngSeed', 10 ...
        );
    [noisyInstancesFixedSeed2, ~] = theNeuralEngine.computeNoisyInstances(...
        noiseFreeResponse, ...
        theResponseTemporalSupportSeconds, ...
        instancesNum, ...
        'random', ...
        'rngSeed', 10 ...
        );
    [noisyInstancesFixedSeed3, ~] = theNeuralEngine.computeNoisyInstances(...
        noiseFreeResponse, ...
        theResponseTemporalSupportSeconds, ...
        instancesNum, ...
        'random', ...
        'rngSeed', 11 ...
        );
    if (any(noisyInstancesFixedSeed2 ~= noisyInstancesFixedSeed1))
        error('Fixing seed does not result in identical noisy response');
    end
      if (all(noisyInstancesFixedSeed3 == noisyInstancesFixedSeed1))
        error('Changing seed does not change noisy response');
    end
end

function renderNeuralResponseSequence(figNo, theResponseSequence, theResponseTemporalSupportSeconds, titleLabel)
    figure(figNo); clf;
    meanResponse = squeeze(mean(theResponseSequence,1));
    cellsNum = size(meanResponse ,1);
    subplot(1,2,1);
    imagesc(theResponseTemporalSupportSeconds, 1:cellsNum, squeeze(theResponseSequence(1,:,:)));
    xlabel('time (sec)');
    ylabel('cells');
    title(sprintf('1st response instance (%s)', titleLabel));
    
    subplot(1,2,2);
    imagesc(theResponseTemporalSupportSeconds, 1:cellsNum, meanResponse);
    xlabel('time (sec)');
    ylabel('cells');
    title(sprintf('mean response (%s)', titleLabel));
    colormap(gray);
end
