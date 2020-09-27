function t_neuralResponseCompute
% Use the @neuralResponseEngine object to compute cone mosaic isomerization responses
%
% Syntax:
%    t_neuralResponseCompute
%
% Description:
%    Demonstrates how to generate a stimulus sequence using a
%    @sceneEngine object and a neural compute function. The neural pipeline defined
%    by the 'photopigmentExcitationsWithNoEyeMovements' compute function represents
%    the cone excitations in the absence of fixational eye movements. Here we are
%    passing a custom neural response params struct during instantiation of the
%    @neuralResponseEngine object. If no neural response params struct were passed, 
%    the default neural response params defined in the 'photopigmentExcitationsWithNoEyeMovements' 
%    compute function would be used.
%
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
% See also: t_sceneGeneration
%

% History:
%    09/21/2020  NPC  Wrote it.

    % Configure the function handle and the params for the @sceneGenerationEngine
    % This is a function that the USER has to supply
    sceneComputeFunction = @uniformFieldTemporalModulation;
    % Instantiate a sceneGenerationEngine with the above sceneComputeFunctionHandle
    theSceneEngine = sceneEngine(sceneComputeFunction);
    
    % Configure the function handle and the params for the @neuralResponseEngine
    % This is a function that the USER has to supply
    neuralComputeFunction = @photopigmentExcitationsWithNoEyeMovements;
    
    % Custom neural response params struct 
    customNeuralResponseParams = struct(...
        'opticsParams',  struct(...
            'type', 'wvf human', ...
            'pupilDiameterMM', 2.0 ...
        ), ...
        'coneMosaicParams',  struct(...
            'upsampleFactor', 5, ...
            'fovDegs', 0.1, ...
            'timeIntegrationSeconds', 5/1000) ...
    );

    % Instantiate a neuralResponseEngine with custom neural response params
    theNeuralEngine = neuralResponseEngine(neuralComputeFunction, customNeuralResponseParams);
    
    % Specify a pedestal luminance with 70% contrast
    testContrast = 0.7;

    % Compute the scene sequence
    [theSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(testContrast);

    % Compute 8 instances of neural responses to the input scene sequence
    instancesNum = 8;
    % Types of noise for thr computed responses
    % This has to be a cell array list with valid entries: {'none', 'random', 'rngSeed_someInt'}
    noiseFlags = {'rNgSeed345', 'rNgSeed345', 'rNgSeed346', 'random', 'none'};
    [theResponses, theResponseTemporalSupportSeconds] = theNeuralEngine.compute(...
            theSceneSequence, ...
            theSceneTemporalSupportSeconds, ...
            instancesNum, ...
            'noiseFlags', noiseFlags ...
            );
    
    [theResponses, theResponseTemporalSupportSeconds] = theNeuralEngine.compute(...
            theSceneSequence, ...
            theSceneTemporalSupportSeconds, ...
            instancesNum, ...
            'noiseFlags', noiseFlags ...
            );
        
    % Visualize all responses computes (different figures for different
    % noise flags)
    debugNeuralResponseGeneration = true;
    if (debugNeuralResponseGeneration)
        for idx = 1:length(noiseFlags)
            renderNeuralResponseSequence(idx, theResponses(noiseFlags{idx}), theResponseTemporalSupportSeconds, noiseFlags{idx});
        end
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
