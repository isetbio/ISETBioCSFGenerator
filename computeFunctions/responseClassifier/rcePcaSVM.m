function dataOut = rcePcaSVMTAFC(obj, operationMode, classifierParamsStruct, theResponses, whichAlternatives)
% Compute function for training a binary SVM classifier to predict response classes from response instances.
%
% Syntax:
%   dataOut = rcePcaSVMClassifier(obj, operationMode,classifierParamsStruct, theReponses, whichAlternatives)

% Description:
%    Compute function to be used as a computeFunctionHandle for a @responseClassifierEngine
%    object for training a binary, two-alternative forced choise classifier. The classifier
%    consists of a PCA-data dimensionality stage coupled with a binary SVM classifier.
%
%    When called with no arguments, this returns a default parameters
%    struct for itself. Currently, there are no parameters used for this
%    routine, so a struct with an unintersting field is returned just to
%    keep the calling machinery happy.
%
%    When called from a parent @responseClassifierEngine object with
%    operationMode set to 'train', it sets up the SVM, using the passed
%    training instances.
%
%    When this is called from a parent @responseClassifierEngine object
%    with operationMode set to 'predict', it makes correct/incorrect
%    predictions for the passed instances. The predictions are made using
%    the SVM.
%
%    This function models an N-way forced choice task with one stimulus
%    presented per trial, so that each trial is has responses to be one of the
%    N possible alternatives.  However, N is currently restricted to 2,
%    because we haven't implemented an N-way SVM linear classifer yet.  But
%    the inferace is set up for N to ease that transition in the future.
%
%    It can be used to handle TAFC if you stack the two alternatives into
%    one long response, which is done (e.g.) in computePerformance when the
%    TAFC flag is set.
%
%    Typically, training will be with one noisy instances of each of
%    the N alternatives.
%
% Inputs:
%    obj                      - The calling @responseClassifierEngine.
%    operationMode            - String, may be either 'train' or
%                               'predict'. When set to 'train', the
%                               function stores into the passed object the
%                               information it needs to predict, using the
%                               passed instances.  When set to 'predict'
%                               the function predicts correct/incorrect on
%                               the passed instances.
%    classifierParamsStruct   - Paramters struct passed by the object.  Can
%                               be used to control details of how the classifier
%                               is set up. Currently this function does not
%                               use any parameters.
%               In 'train' mode:
%                    theResponses   - an N-dimensional cell array, with each entry
%                                     a [mInstances x nDims x tTimePoints] matrix
%                                     of responses for one of the N alternatives.
%                    whichAlternatives - Optional. Not used and can be
%                                     omitted or passed as empty.
%
%               In 'predict' mode:
%                    theResponses   - an [mTrials x nDims x nTimePoints]
%                                     matrix of responses for each trial.
%                    whichAlternatives - Vector of dimension mTrials whose
%                                     entries are integers in the range [1,N] and which
%                                     specify the alternative presented on each trial.
%                                     If this is omitted or passed as
%                                     empty, the routine still predicts the
%                                     alternative presented on each trial,
%                                     but returns responses and pCorrect
%                                     fields below as NaN.
%
% Outputs:
%    dataOut                  - If function called with no input arguments,
%                               this is the default parameter structure for
%                               this funciton.
%
%                               If function is called from a
%                               @responseClassifierEngine object, what this
%                               contains depends on whether operationMode
%                               is 'train' or 'predict'.
%
%                               For 'train', the struct has two fields as
%                               required by responseClassifierEngine
%                               objects.
%                                 .trainedClassifier: Set to empty matrix
%                                 .preProcessingConstants: Structure with
%                                   templates to be used by classifier to
%                                   predict.
%
%                               For 'predict', the struct has two fields as
%                               required by responseClassifierEngine
%                               objects.
%                                 .trialPredictions : Vector with trial-by-trial
%                                                     correct/incorrect
%                                                     predictions. 1 means
%                                                     correct, 0 means incorrect.
%                                 .pCorrect         : Probability correct, which is
%                                                     just the mean of the trialPredictions
%                                                     field.
%                                 .whichAlternativesPredicted : Vector with
%                                                   trial-by-trial alternative predicted.
%                                The trialPredictions and pCorrect fields
%                                come back as NaN if whichAlternatives is
%                                not passed, or passed as empty.
%
% Optional key/value pairs:
%   None.
%
% See Also:
%     t_spatialCSF, computePerformance

% History:
%    09/26/2020  NPC  Wrote it.
%    12/21/2024  dhb  Update not to specialize for TAFC here.

% For consistency with the interface
if (nargin == 0)
    dataOut = struct('Classifier', 'Poisson Forced Choice Ideal Observer');
    return;
elseif (nargin < 5 | isempty(whichAlternatives))
    whichAlternatives = [];
end

% Check operation mode
if (~strcmp(operationMode,'train') && ~strcmp(operationMode,'predict'))
    error('Unknown operation mode passed.  Must be ''train'' or ''predict''');
end

% Feature preprocessing analysis
if (strcmp(operationMode, 'train'))
    % Make sure we are dealing with just 2 alternatives training time
    if (length(theResponses) ~= 2)
        error('Currently we are set up to deal with just 2 alternatives');
    end

    % Feature assembly for classifier training.  This concatenates all the
    % instances and returns a vector of labels according to instance
    [features, whichAlternativesTraining] = assembleTrainingFeatures(theResponses);

    % Compute pre-processing constants: centering factor & principal components
    m = mean(features,1);

    % Center the in-sample eatures
    features = bsxfun(@minus, features, m);

    % Determine PCs without centering (we're doing the centering so we can
    % do it the same way for the out-of-sample data during the prediction
    % phase)
    principalComponents = pca(features,'Centered', false,'NumComponents', classifierParamsStruct.PCAComponentsNum);
else
    % If which alternatives passed, make sure max is 2
    if (~isempty(whichAlternatives))
        if (max(whichAlternatives) > 2)
            error('Currently we are set up to deal with just 2 alternatives');
        end
    end

    % Convert passed responses to features format
    features = theResponses(:,:);

    % Retrieve the in-sample pre-processing constants
    m = obj.preProcessingConstants.centering;
    principalComponents = obj.preProcessingConstants.principalComponents;

    % Center the out-of-sample features
    features = bsxfun(@minus, features, m);
end

% Feature dimensionality reduction via projection to the space defined by
% the first classifierParamsStruct.PCAComponentsNum
features = features * principalComponents;

% Classify
if (strcmp(operationMode, 'train'))
    % Train an SVM classifier on the data
    [trainedSVM, whichAlternativesPredicted, pCorrectTraining, decisionBoundary] = ...
        classifierTrain(classifierParamsStruct, features, whichAlternativesTraining);
    responseTraining = whichAlternativesPredicted == whichAlternativesTraining;
else
    % Employ the trained classifier to predict classes for the responses in this data set
    [response, pCorrect, whichAlternativesPredicted] = classifierPredict(obj.trainedClassifier, features, whichAlternatives);
end

if (strcmp(operationMode, 'train'))
    % Return the trained SVM and the preprocessing constants.
    % These fields are required by the calling object.
    dataOut.trainedClassifier = trainedSVM;
    dataOut.preProcessingConstants = struct(...
        'centering', m, ...
        'principalComponents', principalComponents);
    dataOut.decisionBoundary = decisionBoundary;

    % These are optional and specific to this compute function.
    dataOut.trialPredictionsTraining = responseTraining;
    dataOut.pCorrectTraining = pCorrectTraining;
    dataOut.whichAlternativesdTraining = whichAlternativesTraining;
    dataOut.whichAlternativesPredictedTraining = whichAlternativesPredicted;
else
    % Return the out-of-sample pCorrect, and the trial-by-trial vector
    % of predictions (0 == correct, 1 == correct)
    dataOut.trialPredictions = response;
    dataOut.pCorrect = pCorrect;
    dataOut.whichAlternativePredicted = whichAlternativesPredicted;
end
end


function [response, pCorrect, whichAlternativesPredicted] = classifierPredict(trainedSVM, features, actualClassLabels)

N = size(features,1);
whichAlternativesPredicted  = predict(trainedSVM, features);
if (~isempty(actualClassLabels))
    response = whichAlternativesPredicted(:)==actualClassLabels(:);
else
    response = NaN*ones(size(whichAlternativesPredicted));
end
pCorrect = mean(response);

end

function [compactSVMModel, predictedClassLabels, pCorrectInSample, decisionBoundary] = classifierTrain(classifierParamsStruct, features, classLabels)

% Train the SVM classifier
svmModel = fitcsvm(features, classLabels, ...
    'KernelFunction', classifierParamsStruct.kernelFunction);

% Cross-validate the SVM yp obtain pCorrect for the in-Sample data
crossValidatedSVM = crossval(svmModel,'KFold', classifierParamsStruct.crossValidationFoldsNum);

% Compute pCorrect via cross-validation
percentCorrect = 1 - kfoldLoss(crossValidatedSVM,'lossfun','classiferror','mode','individual');
stdErr = std(percentCorrect)/sqrt(classifierParamsStruct.crossValidationFoldsNum);
pCorrectInSample = mean(percentCorrect);

% Extract the trained, compact classifier. The compact SVM saves space
% by not including the train data.
compactSVMModel = compact(svmModel);

% Generated predicted class labels for the training data set
predictedClassLabels = predict(compactSVMModel, features);

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

function [features, classLabels] = assembleTrainingFeatures(theResponses)

% Reshape data to [nTrials x mResponse dimensions]
% nTrials = size(nullResponses,1);
% cellsNum = size(nullResponses,2);
% timeBins = size(nullResponses,3);
% responseSize = cellsNum*timeBins;
% 
% nullResponses = reshape(nullResponses, [nTrials responseSize]);
% testResponses = reshape(testResponses, [nTrials responseSize]);

features = [];
classLabels = [];
for nn = 1:length(theResponses)
    nTrialsHere = size(theResponses{nn},1);
    features = cat(1, features, theResponses{nn}(:,:));
    classLabels = cat(1, classLabels, nn*ones(nTrialsHere,1));
end

end


function p = generateDefaultParams()
p = struct(...
    'PCAComponentsNum', 2, ...
    'classifierType', 'svm', ...
    'kernelFunction', 'linear', ...
    'crossValidationFoldsNum', 10);
end
