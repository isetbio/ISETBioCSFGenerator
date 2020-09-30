function t_responseClassifier
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

    % Instantiate a sceneGenerationEngine with the above sceneComputeFunctionHandle
    % No sceneParams passed, so we are using the default params specified
    % in the sceneComputeFunction
    theSceneEngine = sceneEngine(sceneComputeFunction);
    
    % Configure the function handle and the params for the @neuralResponseEnginey
    neuralComputeFunction = @photopigmentExcitationsWithNoEyeMovements;

    % Instantiate a neuralResponseEngine. No responseParams passed, so we
    % are using the default params specified in the neuralComputeFunction
    theNeuralEngine = neuralResponseEngine(neuralComputeFunction);
    
    % User-supplied computeFunction for the @responseClassifierEngine
    classifierComputeFunction = @pcaSVMClassifier;
    % User-supplied struct with params appropriate for the @responseClassifierEngine computeFunction
    customClassifierParams = struct(...
        'PCAComponentsNum', 2, ...          % number of PCs used for feature set dimensionality reduction
        'taskIntervals', 2, ...             % simulate a 2-interval task
        'crossValidationFoldsNum', 10, ...  % employ a 10-fold cross-validated linear 
        'kernelFunction', 'linear', ...     % linear
        'classifierType', 'svm' ...         % binary SVM classifier
        );
   
    % Instantiate a responseClassifierEngine with the above classifierComputeFunctionHandle and custom classifierParams
    theClassifierEngine = responseClassifierEngine(classifierComputeFunction);
   
    % Generate the NULL stimulus sequence
    nullContrast = 0.0;
    % Compute the NULL stimulus
    [theNullSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(nullContrast);

    % Generate the TEST stimulus sequence
    testContrast = 0.008;
    % Compute the TEST stimulus
    [theTestSceneSequence, ~] = theSceneEngine.compute(testContrast);
    
    % Compute respose instances to the NULL and TEST stimuli for training the classifier
    trainingInstancesNum = 128;
    
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

    % Train the binary classifier on the above NULL/TEST response set
    trainingData = theClassifierEngine.compute('train',...
        inSampleNullStimResponses('random'), ...
        inSampleTestStimResponses('random'));
    
    % Test predictions of the trained classifier on 10 out of sample
    % response instance at a time
    outOfSampleInstancesNum = 50;
    outOfSampleInstancesNum = max([outOfSampleInstancesNum theClassifierEngine.classifierParams.taskIntervals]);
    
    % Repeat out-of-sample predictions a total of N times
    N = 1;
    for trial = 1:N
        % Compute new respose instances to the NULL and TEST stimuli
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

        % Run the classifier on the new response instances
        predictedData = theClassifierEngine.compute('predict',...
            outOfSampleNullStimResponses('random'), ...
            outOfSampleTestStimResponses('random'));

        % Visualize the classifier performance on the in-sample (training) and 
        % the out of sample responses. Notice that the second component
        % captures almost variance in the test data. Thats because the principal
        % components are computed on the training data. In that data set the second
        % component is capturing the noise. In the second data set, the
        % second component captures almost zero variance because the noise component
        % is different.
        debugClassifier = true;
        if (debugClassifier)
             plotClassifierResults(theClassifierEngine.classifierParams, trainingData, predictedData);
             drawnow;
        end
   end
        
end

function plotClassifierResults(classifierParams, trainingData, predictedData)
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 950 500], 'Color', [1 1 1]);
    minFeature = min([min(trainingData.features(:))]);
    maxFeature = max([max(trainingData.features(:))]);
    
    % The training data
    ax = subplot(1,2,1);
    hold(ax, 'on');
    renderDecisionBoundary(ax,trainingData.decisionBoundary, true); 
    renderFeatures(ax, trainingData.features, classifierParams.taskIntervals);
    xlabel(ax,'PCA #1 score');
    ylabel(ax,'PCA #2 score');
    title(ax,sprintf('In-sample percent correct: %2.3f', trainingData.pCorrectInSample));
    axis(ax,'square');
    set(ax, 'XLim', [minFeature maxFeature], 'YLim', [minFeature maxFeature], 'FontSize', 12);
    colormap(ax,brewermap(1024, 'RdYlGn'));

    
    minFeature = min([min(predictedData.features(:))]);
    maxFeature = max([max(predictedData.features(:))]);
    
    % The predicted data
    ax = subplot(1,2,2);
    hold(ax, 'on');
    renderDecisionBoundary(ax,trainingData.decisionBoundary, false); 
    renderFeatures(ax, predictedData.features, classifierParams.taskIntervals);
    xlabel(ax,'PCA #1 score');
    ylabel(ax,'PCA #2 score');
    title(ax,sprintf('Out-of-sample percent correct: %2.3f', predictedData.pCorrectOutOfSample));
    axis(ax,'square');
    set(ax, 'XLim', [minFeature maxFeature], 'YLim', [minFeature maxFeature],  'FontSize', 12);
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


function renderDecisionBoundary(ax, decisionBoundary, depictStrength)
    if (~isempty(decisionBoundary))
        N = length(decisionBoundary.x);
        if (depictStrength)
        % Decision boundary as a density plot
            imagesc(ax,decisionBoundary.x,decisionBoundary.y,reshape(decisionBoundary.z,[N N]));
        end
        % The decision boundary as a line
        [C,h] = contour(ax,decisionBoundary.x,decisionBoundary.y,reshape(decisionBoundary.z,[N N]), [0 0]);
        h.LineColor = [0 0 0];
        h.LineWidth = 2.0;
    end
end