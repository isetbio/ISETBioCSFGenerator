function t_testResponseClassifier
% Use the @responseClassifierEngine object to compute NULL/TEST discriminability 
%
% Syntax:
%    t_testResponseClassifier
%
% Description:
%    Demonstrates how to generate a stimulus sequence using a
%    @sceneEngine object, how to compute cone mosaic excitation 
%    responses to a test stimulus and the null stimulus using a @neuralResponseEngine 
%    object and how to compute the probability with which the test stimulus
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
% See Also:
%   t_neuralResponseCompute
%   t_sceneGeneration

% History:
%    09/21/2020  NPC  Wrote it.

    % Configure the function handle and the params for the @sceneGenerationEngine
    sceneComputeFunction = @uniformFieldTemporalModulation;

    % Configure the function handle and the params for the @neuralResponseEnginey
    neuralComputeFunction = @photopigmentExcitationsWithNoEyeMovements;


    % =======  RESPONSE CLASSIFIER PARAMS  ======= 
    % Configure the function handle and the params for the @responseClassifierEngine
    % This is a function that the USER has to supply
    classifierComputeFunction = @pcaSVMClassifier;
    % This is a struct that the USER has to supply and which is to work
    % with the user-supplied function handle
    customClassifierParams = struct(...
        'PCAComponentsNum', 2, ...
        'taskIntervals', 2, ...
        'classifierType', 'svm', ...
        'kernelFunction', 'linear', ...
        'crossValidationFoldsNum', 10);
    
    
    % STEP 0. Instantiate the engines

    % Instantiate a sceneGenerationEngine with the above sceneComputeFunctionHandle and sceneParams
    theSceneEngine = sceneEngine(sceneComputeFunction);
    
    % Instantiate a neuralResponseEngine with the desired neural response params
    theNeuralEngine = neuralResponseEngine(neuralComputeFunction);
    
    % Instantiate a responseClassifierEngine with the above classifierComputeFunctionHandle and classifierParams
    theClassifierEngine = responseClassifierEngine(classifierComputeFunction, customClassifierParams);
   
    
    % STEP 1. Generate the NULL and TEST stimulus scene sequences
    
    % The NULL stimulus scene
    nullContrast = 0.0;
    % Compute the TEST stimulus
    [theNullSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(nullContrast);

    % The TEST stimulus scene
    testContrast = 0.7/100;
    % Compute the TEST stimulus
    [theTestSceneSequence, ~] = theSceneEngine.compute(testContrast);
    
    
    % Step 2. Compute NULL and TEST response data for training the classifier
    trainingInstancesNum = 256;
    
    [inSampleNullStimResponses, theIsomerizationsTemporalSupportSeconds] = ...
        theNeuralEngine.compute(theNullSceneSequence, ...
        theSceneTemporalSupportSeconds, ...
        trainingInstancesNum, ...
        'noiseFlags', {'random'});
    
    [inSampleTestStimResponses, ~] = theNeuralEngine.compute(...
        theTestSceneSequence, ...
        theSceneTemporalSupportSeconds, ...
        trainingInstancesNum, ...
        'noiseFlags', {'random'});

    
    % Step 3. Train the classifier on the inSample responses
    trainingData = theClassifierEngine.compute(...
        inSampleNullStimResponses('random'), ...
        inSampleTestStimResponses('random'), ...
        'train');
    
    
    outOfSampleInstancesNum = 2;
    
    for trial = 1:100
        % Step 4. Compute NULL and TEST response data for testing the classifier
        outOfSampleNullStimResponses = theNeuralEngine.compute(...
            theNullSceneSequence, ...
            theSceneTemporalSupportSeconds, ...
            outOfSampleInstancesNum, ...
            'noiseFlags', {'random'});

        outOfSampleTestStimResponses = theNeuralEngine.compute(...
            theTestSceneSequence, ...
            theSceneTemporalSupportSeconds, ...
            outOfSampleInstancesNum, ...
            'noiseFlags', {'random'});

        predictedData = theClassifierEngine.compute(...
            outOfSampleNullStimResponses('random'), ...
            outOfSampleTestStimResponses('random'), ...
            'predict');

        % Visualize
         debugClassifier = true;
         if (debugClassifier)
             plotClassifierResults(theClassifierEngine.classifierParams, trainingData, predictedData);
             drawnow;
         end
   end
        
end

function plotClassifierResults(classifierParams, trainingData, predictedData)
    figure(1); clf;
    
    minFeature = min([min(trainingData.features(:)) min(predictedData.features(:))]);
    maxFeature = max([max(trainingData.features(:)) max(predictedData.features(:))]);
    
    % The training data
    ax = subplot(1,2,1);
    hold(ax, 'on');
    renderDecisionBoundary(ax,trainingData.decisionBoundary); 
    renderFeatures(ax, trainingData.features, classifierParams.taskIntervals);
    xlabel(ax,'PCA #1 score');
    ylabel(ax,'PCA #2 score');
    title(ax,sprintf('In-sample percent correct: %2.3f', trainingData.pCorrectInSample));
    axis(ax,'square');
    set(ax, 'XLim', [minFeature maxFeature], 'YLim', [minFeature maxFeature]);
    colormap(ax,brewermap(1024, 'RdYlGn'));

    
    % The predicted data
    ax = subplot(1,2,2);
    hold(ax, 'on');
    renderDecisionBoundary(ax,trainingData.decisionBoundary); 
    renderFeatures(ax, predictedData.features, classifierParams.taskIntervals);
    xlabel(ax,'PCA #1 score');
    ylabel(ax,'PCA #2 score');
    title(ax,sprintf('Out-of-sample percent correct: %2.3f', predictedData.pCorrectOutOfSample));
    axis(ax,'square');
    set(ax, 'XLim', [minFeature maxFeature], 'YLim', [minFeature maxFeature]);
    colormap(ax,brewermap(1024, 'RdYlGn'));
    
end

function renderFeatures(ax, features, taskIntervals)

    N = max([1 floor(size(features,1)/2)]);

    if (taskIntervals == 2)
            scatter(ax,features(1:N,1), features(1:N,2), 64, ...
                'MarkerFaceColor', [0.8 0.8 0.8], ...
                'MarkerEdgeColor', [0.2 0.2 0.2]); 
            hold('on')
            if (size(features,1)<=2*N)
                scatter(ax,features(N+(1:N),1), features(N+(1:N),2), 64, ...
                'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerEdgeColor', [0.9 0.9 0.9]);
            end
    else
        
        scatter(ax,features(1:N,1), features(1:N,2), 64, ...
            'MarkerFaceColor', [0.8 0.8 0.8], ...
            'MarkerEdgeColor', [0.2 0.2 0.2]); 
        hold('on')
        scatter(ax,features(N+(1:N),1), features(N+(1:N),2), 64, ...
            'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerEdgeColor', [0.9 0.9 0.9]);
    end

end


function renderDecisionBoundary(ax, decisionBoundary)
    if (~isempty(decisionBoundary))
        N = length(decisionBoundary.x);
        % Decision boundary as a density plot
        imagesc(ax,decisionBoundary.x,decisionBoundary.y,reshape(decisionBoundary.z,[N N]));
        % The decision boundary as a line
        [C,h] = contour(ax,decisionBoundary.x,decisionBoundary.y,reshape(decisionBoundary.z,[N N]), [0 0]);
        h.LineColor = [0 0 0];
        h.LineWidth = 2.0;
    end
end