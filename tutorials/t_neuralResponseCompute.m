function t_neuralResponseCompute
% Use the @neuralResponseEngine object to compute cone mosaic isomerization responses
%
% Syntax:
%    t_neuralResponseCompute
%
% Description:
%    Demonstrates how to generate a stimulus sequence using a
%    @sceneEngine object and a neural compute function. The neural pipeline defined
%    by the 'nrePhotopigmentExcitationsConeMosaicHexWithNoEyeMovements' compute function represents
%    the cone excitations in the absence of fixational eye movements. Here we are
%    passing a custom neural response params struct during instantiation of the
%    @neuralResponseEngine object. If no neural response params struct were passed, 
%    the default neural response params defined in the 'nrePhotopigmentExcitationsConeMosaicHexWithNoEyeMovements' 
%    compute function would be used.
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
    theSceneEngine = sceneEngine(sceneComputeFunction);
    
    % Configure the function handle and the params for the @neuralResponseEngine
    % This is a function that the USER has to supply
    neuralComputeFunction = @nrePhotopigmentExcitationsConeMosaicHexWithNoEyeMovements;
    
    % Custom neural response params struct. The form of this structure is
    % defined by and is specific to the
    % nrePhotopigmentExcitationsConeMosaicHexWithNoEyeMovements compute function.  You
    % can write your own compute function for the model visual system you
    % are interested in, and it can take parameters defined by a structure
    % that you set up as part of writing that function.
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

    % Instantiate a neuralResponseEngine with the custom neural response params
    theNeuralEngine = neuralResponseEngine(neuralComputeFunction, customNeuralResponseParams);
    
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
    noiseFlags = {'random', 'none'};
    [theResponses, theResponseTemporalSupportSeconds] = theNeuralEngine.compute(...
            theSceneSequence, ...
            theSceneTemporalSupportSeconds, ...
            instancesNum, ...
            'noiseFlags', noiseFlags, ...
            'rngSeed', [] ...
            );
        
    % The responses come back as a Matlab container, with one container
    % entry per noise flag passed.  Extract the response matrix from the
    % container by indexing it with the noise flag.   These are arranged as
    % an instancesNum x nDim x nTimepoints matrix, where nDim is the
    % dimension of the neural response at a single time point and
    % nTimespoints is the number of time points in the sequence.
    noiseFreeResponses = theResponses('none');
    
    % The noisy ('random') responses instances also comeback in the
    % container.  These are arranged as an instancesNum x nDim x nTimepoints
    % matrix, where nDim is the dimension of the neural response at a
    % single time point and nTimespoints is the number of time points in
    % the sequence.
    noisyInstances = theResponses('random');
    assert(size(noisyInstances,1) == instancesNum);
    assert(size(noisyInstances,3) == length(theResponseTemporalSupportSeconds));
        
    % Visualize all responses computes (different figures for different
    % noise flags)
    debugNeuralResponseGeneration = true;
    if (debugNeuralResponseGeneration)
        for idx = 1:length(noiseFlags)
            renderNeuralResponseSequence(idx, theResponses(noiseFlags{idx}), theResponseTemporalSupportSeconds, noiseFlags{idx});
        end
    end
    
    % Generate 1 noisy instance with a specified rng seed
    instancesNum = 1;
    % the noise flag must contain the 'rngSeed' substring
    noiseFlags = {'rngSeed whatever'};
    [theResponses, theResponseTemporalSupportSeconds] = theNeuralEngine.compute(...
            theSceneSequence, ...
            theSceneTemporalSupportSeconds, ...
            instancesNum, ...
            'noiseFlags', noiseFlags, ...
            'rngSeed', 123456 ...
            );
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
