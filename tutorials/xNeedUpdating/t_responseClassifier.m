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
%   t_thresholdEngine, t_sceneGeneration, t_neuralResponseCompute
%

% History:
%    09/21/2020  NPC  Wrote it.

    % Close all figs
    close all;

    % Configure the function handle and the params for the @sceneGenerationEngine
    sceneParams = sceUniformFieldTemporalModulation;
    sceneComputeFunction = @sceUniformFieldTemporalModulation;

    % Instantiate a sceneGenerationEngine with the above sceneComputeFunctionHandle
    % No sceneParams passed, so we are using the default params specified
    % in the sceneComputeFunction
    theSceneEngine = sceneEngine(sceneComputeFunction);
    
    % Configure the function handle and the params for the
    % @neuralResponseEngine.
    % 
    % Set integration time to match scene sequence frame duration.  These
    % need to match.
    nreParams = nrePhotopigmentExcitationsCmosaic;
    nreParams.coneMosaicParams.timeIntegrationSeconds = sceneParams.frameDurationSeconds;
    neuralComputeFunction = @nrePhotopigmentExcitationsCmosaic;

    % Instantiate a neuralResponseEngine. No responseParams passed, so we
    % are using the default params specified in the neuralComputeFunction
    theNeuralEngine = neuralResponseEngine(neuralComputeFunction,nreParams);
    
    % User-supplied computeFunction for the @responseClassifierEngine
    classifierComputeFunction = @rcePcaSVMTAFC;
    
    % User-supplied struct with params appropriate for the @responseClassifierEngine computeFunction
    customClassifierParams = struct(...
        'PCAComponentsNum', 2, ...          % number of PCs used for feature set dimensionality reduction
        'crossValidationFoldsNum', 10, ...  % employ a 10-fold cross-validated linear 
        'kernelFunction', 'linear', ...     % linear
        'classifierType', 'svm' ...         % binary SVM classifier
        );
   
    % Instantiate a responseClassifierEngine with the above classifierComputeFunctionHandle and custom classifierParams
    theClassifierEngine = responseClassifierEngine(classifierComputeFunction, customClassifierParams);
   
    % Generate the NULL stimulus sequence
    nullContrast = 0.0;
    [theNullSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(nullContrast);

    % Generate the TEST stimulus sequence
    testContrast = 0.008;
    [theTestSceneSequence, ~] = theSceneEngine.compute(testContrast);
    
    % Compute respose instances to the NULL and TEST stimuli for training the classifier
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

    % Train the binary classifier on the above NULL/TEST response set
    trainingData = theClassifierEngine.compute('train',...
        inSampleNullStimResponses('random'), ...
        inSampleTestStimResponses('random'));
    
    % Test predictions of the trained classifier on out of sample
    % response instances.  This parameter tells us how many out of sample 
    % stimulus instances to use.
    outOfSampleInstancesNum = 100;
  
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
        % captures almost no variance in the test data. Thats because the principal
        % components are computed on the training data. In that data set the second
        % component is capturing the noise. In the second data set, the
        % second component captures almost zero variance because the noise component
        % is different.
        debugClassifier = true;
        if (debugClassifier)
             plotClassifierResults(trainingData, predictedData);
             drawnow;
        end
   end
        
end

function plotClassifierResults(trainingData, predictedData)
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 950 500], 'Color', [1 1 1]);
    minFeature = min([ min(trainingData.features(:)) min(predictedData.features(:)) ]);
    maxFeature = max([ max(trainingData.features(:)) max(predictedData.features(:)) ]);
    
    % The training data
    ax = subplot(1,2,1);
    hold(ax, 'on');
    renderDecisionBoundary(ax,trainingData.decisionBoundary, true); 
    renderFeatures(ax, trainingData.features, trainingData.nominalClassLabels);
    xlabel(ax,'PCA #1 score');
    ylabel(ax,'PCA #2 score');
    title(ax,sprintf('In-sample percent correct: %2.3f', trainingData.pCorrect));
    axis(ax,'square');
    set(ax, 'XLim', [minFeature maxFeature], 'YLim', [minFeature maxFeature], 'FontSize', 12);
    colormap(ax,brewermap(1024, 'RdYlGn'));


    % The predicted data
    ax = subplot(1,2,2);
    hold(ax, 'on');
    renderDecisionBoundary(ax,trainingData.decisionBoundary, false); 
    renderFeatures(ax, predictedData.features, predictedData.nominalClassLabels);
    xlabel(ax,'PCA #1 score');
    ylabel(ax,'PCA #2 score');
    title(ax,sprintf('Out-of-sample percent correct: %2.3f', predictedData.pCorrect));
    axis(ax,'square');
    set(ax, 'XLim', [minFeature maxFeature], 'YLim', [minFeature maxFeature],  'FontSize', 12);
    colormap(ax,brewermap(1024, 'RdYlGn'));
    
end

function renderFeatures(ax, features, nominalLabels)

    idx = find(nominalLabels == 0);
    scatter(ax,features(idx,1), features(idx,2), 64, ...
        'MarkerFaceColor', [0.8 0.8 0.8], ...
        'MarkerEdgeColor', [0.2 0.2 0.2]); 
    hold('on')
    idx = find(nominalLabels == 1);
    scatter(ax,features(idx,1), features(idx,2), 64, ...
        'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerEdgeColor', [0.9 0.9 0.9]);

    legend({'decision boundary', 'nominal class 0', 'nomimal class 1'})
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