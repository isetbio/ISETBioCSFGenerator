function dataOut = pcaSVMClassifier(classifierOBJ, theNULLResponses, theTESTResponses, classifierParams, operationMode)

    % Default params for this compute function
    defaultClassifierParams = generateDefaultParams();

    % Check input arguments. If called with zero input arguments, just return the default params struct
    if (nargin == 0)
        dataOut = defaultClassifierParams;
        return;
    end
    
    % Feature assembly phase
    [features, classLabels] = computeFeatures(classifierParams, theNULLResponses, theTESTResponses);
    
    % Feature dimensionality reduction phase
    if (strcmp(operationMode, 'train'))
        % Pre-processing: standardizing
        % Determine centering and scaling constants
        m = mean(features,1);
        s = std(features,1,1);
        % Standardize
        features = bsxfun(@times, bsxfun(@minus, features, m), 1./s);

        % Pre-processing: projection to PC space
        % Determine PCs
        principalComponents = pca(features,'NumComponents', classifierParams.PCAComponentsNum);
        % Project all responses to the space defined by the first classifierParams.PCAComponentsNum
        features = features * principalComponents;
    else
        % Pre-processing: standardizing using the centering/scaling
        % constants derived during the training phase
        m = classifierOBJ.preProcessingConstants.centering;
        s = classifierOBJ.preProcessingConstants.scaling;
        features = bsxfun(@times, bsxfun(@minus, features, m), 1./s);
        
        % Pre-processing: projection to PC space using the PC's derived during the training phase
        features = features * classifierOBJ.preProcessingConstants.principalComponents;    
    end

    % Classification phase
    if (strcmp(operationMode, 'train'))
        % Train an SVM classifier on the data
        [trainedSVM, pCorrectInSample, decisionBoundary] = classifierTrain(classifierParams, features, classLabels);
    else
        % Employ the trained classifier to classify the new data
        [predictedClassLabels, pCorrectOutOfSample] = classifierPredict(classifierOBJ.trainedClassifier, features, classLabels);
    end
    
     % Assemble dataOut struct
     dataOut.features = features;
     if (strcmp(operationMode, 'train'))
        dataOut.trainedClassifier = trainedSVM;
        dataOut.preProcessingConstants = struct(...
            'centering', m, ...
        	'scaling', s, ...
        	'principalComponents', principalComponents);
        dataOut.pCorrectInSample = pCorrectInSample;
        dataOut.decisionBoundary = decisionBoundary;
     else
        dataOut.pCorrectOutOfSample = pCorrectOutOfSample;
        dataOut.predictedClassLabels = predictedClassLabels;
     end
     
end


function [predictedClassLabels, pCorrectOutOfSample] = classifierPredict(trainedSVM, features, actualClassLabels)
    N = size(features,1);
    predictedClassLabels = predict(trainedSVM, features);
    pCorrectOutOfSample = sum(predictedClassLabels(:)==actualClassLabels(:))/N;
end

function [compactSVMModel, pCorrectInSample, decisionBoundary] = classifierTrain(classifierParams, features, classLabels)
 
     % Train the SVM classifier
     svmModel = fitcsvm(features, classLabels, ...
                'KernelFunction', classifierParams.kernelFunction);
            
     % Cross-validate the SVM
     crossValidatedSVM = crossval(svmModel,'KFold', classifierParams.crossValidationFoldsNum);
     
     % Compute percent correct
     percentCorrect = 1 - kfoldLoss(crossValidatedSVM,'lossfun','classiferror','mode','individual');
     stdErr = std(percentCorrect)/sqrt(classifierParams.crossValidationFoldsNum);
     pCorrectInSample = mean(percentCorrect);
            
     % Extract the trained, compact classifier. The compact SVM saves space
     % by not including the train data.
     compactSVMModel = compact(svmModel);
     
     % Compute the 2D decision boundary
     if (size(features,2) == 2)
         N = 200;
         x = linspace(min(features(:)), max(features(:)), N);
         y = linspace(min(features(:)), max(features(:)), N);
         [X,Y] = meshgrid(x, y);

         % Run the SVM to get assigned classes for all the xy grid points
         highResFeatureGrid = [X(:),Y(:)];
         [~, score] = predict(compactSVMModel,highResFeatureGrid);

         % Get our decision boundary (score)
         decisionBoundary.z = score(:,1);
         decisionBoundary.x = x;
         decisionBoundary.y = y;
     else
         decisionBoundary = [];
     end
     
end

function [features, classLabels] = computeFeatures(classifierParams, theNULLResponses, theTESTResponses)
    
    % Reshape data to [nTrials x mResponse dimensions]
    nTrials = size(theNULLResponses,1);
    cellsNum = size(theNULLResponses,2);
    timeBins = size(theNULLResponses,3);
    responseSize = cellsNum*timeBins;
    
    theNULLResponses = reshape(theNULLResponses, [nTrials responseSize]);
    theTESTResponses = reshape(theTESTResponses, [nTrials responseSize]);
    
    % Arange responses according to whether we are simulating a 1- or a 2-interval task
    if (classifierParams.taskIntervals == 1)
        features = zeros(2*nTrials, responseSize);
        classLabels = zeros(2*nTrials, 1);
        
        % Class 1 data and labels
        features(1:nTrials,:) = theNULLResponses;
        classLabels(1:nTrials,1) = 0;
        
        % Class 2 data and labels
        features(nTrials+(1:nTrials),:) = theTESTResponses;
        classLabels(nTrials+(1:nTrials),1) = 1;
    else
        features = zeros(nTrials, 2*responseSize);
        classLabels = zeros(nTrials, 1);
        halfTrials = floor(nTrials/2);
        
        % Class 1 data and labels
        features(1:halfTrials,:) = cat(2, ...
            theNULLResponses(1:halfTrials,:), ...
            theTESTResponses(1:halfTrials,:));
        classLabels(1:halfTrials,1) = 0;
        
        % Class 2 data and labels
        features(halfTrials+(1:halfTrials),:) = cat(2, ...
            theTESTResponses(halfTrials+(1:halfTrials),:), ...
            theNULLResponses(halfTrials+(1:halfTrials),:));
        classLabels(halfTrials+(1:halfTrials),1) = 1;
    end
end


function p = generateDefaultParams()
    p = struct(...
        'PCAComponentsNum', 2, ...
        'taskIntervals', 1, ...
        'classifierType', 'svm', ...
        'kernelFunction', 'linear', ...
        'crossValidationFoldsNum', 10);
end
