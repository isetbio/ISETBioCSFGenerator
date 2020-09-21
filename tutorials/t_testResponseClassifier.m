function t_testResponseClassifier
% Use the @responseClassifierEngine object to compute NULL/TEST discriminability 
%
% Syntax:
%    t_testResponseClassifier
%
% Description:
%    Demonstrates how to generate a stimulus sequence using a
%    @sceneGeneration object, how to compute cone mosaic excitation 
%    responses to a test stimulus and the null stimulus using a @neuralResponseEngine 
%    object and how to compute the probability with which the rest stimulus
%    can be discriminated from the null stimulus based on the respective
%    responses using a @responseClassifierEngine
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

    % =======  STIMULUS SCENE PARAMS ===========
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

    % ======= NEURAL RESPONSE PARAMS ===========
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
        'fovDegs', 0.1, ...
        'timeIntegrationSeconds', 5/1000 ...
    );

    neuralResponseParams = struct(...
        'opticsParams', opticsParams, ...
        'coneMosaicParams', coneMosaicParams ...
    );

    % =======  RESPONSE CLASSIFIER PARAMS  ======= 
    % Configure the function handle and the params for the @responseClassifierEngine
    % This is a function that the USER has to supply
    classifierComputeFunction = @pcaSVMClassifier;
    % This is a struct that the USER has to supply and which is to work
    % with the user-supplied function handle
    classifierParams = struct(...
        'PCAComponentsNum', 2, ...
        'taskIntervals', 1, ...
        'classifierType', 'svm', ...
        'kernelFunction', 'linear', ...
        'crossValidationFoldsNum', 10);
    
    
    % STEP 0. Instantiate the engines

    % Instantiate a sceneGenerationEngine with the above sceneComputeFunctionHandle and sceneParams
    theSceneEngine = sceneGenerationEngine(sceneComputeFunction, sceneParams);
    
    % Instantiate a neuralResponseEngine with the desired neural response params
    theNeuralEngine = neuralResponseEngine(neuralComputeFunction, neuralResponseParams);
    
    % Instantiate a responseClassifierEngine with the above classifierComputeFunctionHandle and classifierParams
    theClassifierEngine = responseClassifierEngine(classifierComputeFunction, classifierParams);
   
    
    % STEP 1. Generate the stimulus scenes
    
    % The NULL stimulus scene
    nullContrast = 0.0;
    % Compute the TEST stimulus
    [theNullSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(nullContrast);

    % The TEST stimulus scene
    testContrast = 1.5/100;
    % Compute the TEST stimulus
    [theTestSceneSequence, ~] = theSceneEngine.compute(testContrast);
    
    
    % Step 2. Compute responses to the null and test stimuli
    % Compute 100 instances of neural responses to the TEST stimulus
    instancesNum = 100;
    
    [theNULLResponses, theIsomerizationsTemporalSupportSeconds] = ...
        theNeuralEngine.compute(theNullSceneSequence, theSceneTemporalSupportSeconds, instancesNum);
    
    [theTESTResponses, ~] = ...
        theNeuralEngine.compute(theTestSceneSequence, theSceneTemporalSupportSeconds, instancesNum);

    
    % Step 3. Run the classifier to obtain pCorrect
    [pCorrect, features, decisionBoundary] = theClassifierEngine.compute(theNULLResponses, theTESTResponses);
    
    % Visualize
    debugClassifier = true;
    if (debugClassifier)
        renderClassifierResults(features, decisionBoundary, pCorrect);
    end
    
end

function renderClassifierResults(features, decisionBoundary, pCorrect)
    figure(4); clf;
    hold on;
    N = length(decisionBoundary.x);
    % Decision boundary as a density plot
    imagesc(decisionBoundary.x,decisionBoundary.y,reshape(decisionBoundary.z,[N N]));
    % The decision boundary as a line
    [C,h] = contour(decisionBoundary.x,decisionBoundary.y,reshape(decisionBoundary.z,[N N]), [0 0]);
    h.LineColor = [0 0 0];
    h.LineWidth = 2.0;
    
    scatter(features(1:100,1), features(1:100,2), 144, 'MarkerFaceColor', [0.8 0.8 0.8]); 
    hold('on')
    scatter(features(101:200,1), features(101:200,2), 144, 'MarkerFaceColor', [0.4 0.4 0.4]);
    xlabel('PCA #1 score');
    ylabel('PCA #2 score');
    title(sprintf('percent correct: %2.3f', pCorrect));
    axis 'square';
    set(gca, 'XLim', [min(features(:)) max(features(:))], ...
             'YLim', [min(features(:)) max(features(:))]);
    colormap(brewermap(1024, 'RdYlGn'));
end
