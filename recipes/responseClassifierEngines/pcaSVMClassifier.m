function [pCorrect, features, decisionBoundary] = pcaSVMClassifier(classifierParams, theNULLResponses, theTESTResponses)

    % Reshape data to [nTrials x mResponse dimensions]
    nTrials = size(theNULLResponses,1);
    cellsNum = size(theNULLResponses,2);
    timeBins = size(theNULLResponses,3);
    responseSize = cellsNum*timeBins;
    theNULLResponses = reshape(theNULLResponses, [nTrials responseSize]);
    theTESTResponses = reshape(theTESTResponses, [nTrials responseSize]);
    
    % Run PCA on all responses
    allResponses = cat(1, theNULLResponses, theTESTResponses);
        
    % Standardize responses
    m = mean(allResponses,1);
    s = std(allResponses,1,1);
    allResponses = (allResponses - repmat(m,size(allResponses,1),1)) ./ repmat(s,size(allResponses,1),1);
   
    % Find principal components
    principalComponents = pca(allResponses,'NumComponents', classifierParams.PCAComponentsNum);
    % Project all responses to the space defined by the first classifierParams.PCAComponentsNum
    allResponses = allResponses * principalComponents;
    % Extract the projections of the null/test stimulus response instances
    theNULLResponses = allResponses(1:nTrials,:);
    theTESTResponses = allResponses(nTrials+(1:nTrials),:); 
    
    % Assemble classification data according to whether we are simulating a
    % single interval or a 2-interval task
    % Assemble the classification data
    if (classifierParams.taskIntervals == 1)
        features = zeros(2*nTrials, size(theNULLResponses,2));
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
        features(1:halfTrials,:) = cat(2, theNULLResponses(1:halfTrials,:), theTESTResponses(1:halfTrials,:));
        classLabels(1:halfTrials,1) = 0;
        
        % Class 2 data and labels
        features(halfTrials+(1:halfTrials),:) = cat(2, theTESTResponses(halfTrials+(1:halfTrials),:), theNULLResponses(halfTrials+(1:halfTrials),:));
        classLabels(halfTrials+(1:halfTrials),1) = 1;
    end
    
     % Train the SVM classifier
     svm = fitcsvm(features, classLabels, ...
                'KernelFunction', classifierParams.kernelFunction);
     % Cross-validate the SVM
     CVSVM = crossval(svm,'KFold', classifierParams.crossValidationFoldsNum);
     % Compute percent correct
     percentCorrect = 1 - kfoldLoss(CVSVM,'lossfun','classiferror','mode','individual');
     stdErr = std(percentCorrect)/sqrt(classifierParams.crossValidationFoldsNum);
     pCorrect = mean(percentCorrect);
            
     N = 200;
     x = linspace(min(features(:)), max(features(:)), N);
     y = linspace(min(features(:)), max(features(:)), N);
     [X,Y] = meshgrid(x, y);

     % Run the SVM to get assigned classes for all the xy grid points
     pred = [X(:),Y(:)];
     [~, score] = predict(svm,pred);
    
     % Get our decision boundary (score)
     decisionBoundary.z = score(:,1);
     decisionBoundary.x = x;
     decisionBoundary.y = y;
end
