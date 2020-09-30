function dataOut = compute(obj, operationMode, nullResponses, testResponses)
% Generic compute method for the @responseClassifierEngine class.
%
% Syntax:
%   dataOut = compute(obj, operationMode, nullResponses, testResponses)
%
% Description:
%    Generic compute method for the @responseClassifierEngine class. Its purpose 
%    is to provide a unified interface for either training or employing a 
%    trained response classifier independent of the response classifier 
%    specifics. The code for building up and executing the response classifier 
%    and the classifier parameters  (e.g, type of classifier, type of data preprocessing) 
%    are defined in the computeFunctionHandle and the responseClassifierParamsStruct, 
%    respectively, both of  which are set when  instantiating the
%    @responseClassifier object.
%    
%    Apart from providing a unified interface, this method also stores the
%    trained classifier model and any preprocessing constants computed during the
%    training phase using the training (in-sample) data set. These are employed
%    in subsequent calls for predicting class labels to novel (out-of-sample) data sets.
%
% Inputs:
%    obj                            - the parent @responseClassifierEngine object
%                               
%    operationMode                  - a string in {'train', 'predict'}
%                                     defining the operation to be performed: either
%                                     training the classifier on a (nullResponses, testResponses) 
%                                     training data set, or using the trained classifier to
%                                     derive predictions on a previously unseen (out-of-sample)
%                                     (nullResponses, testResponses) data set
%
%    nullResponses                  - an [mTrials x nDims] matrix of responses to the null stimulus
%
%    testResponses                  - an [mTrials x nDims] matrix of responses to the test stimulus
%
%
% Optional key/value input arguments: none
%
% Outputs:
%    dataOut                        - a struct that depends on the operationMode selected
%                                     If the operatonMode is set to 'train', dataOut contains the
%                                     following fields:
%                                     .features                : the features used for classification
%                                     .trainedClassifier       : the trained binary SCV classifer
%                                     .preProcessingConstants  : constants computed during the dimensionality reduction preprocessing phase
%                                     .pCorrectInSample        : probability of correct classification for the in-sample (training data)
%                                     .decisionBounday         : if the feature set is 2D, the 2D decision boundary, otherwise []
%
%                                     If the operatonMode is set to 'predict', dataOut contains the
%                                     following fields:
%                                     .features                : the features used for classification
%                                     .pCorrectOutOfSample     : probability of correct classification for the out-of-sample (testing data)
%                                     .predictedClassLabels    : the predicted labels for the out-of-sample responses
%
% See Also:
%     t_responseClassifier

% History:
%    9/20/2020  NPC Wrote it
    % Validate the operationMode
    assert(ismember(operationMode, obj.validOperationModes), sprintf('The passed responseClassifierEngine.compute() ''%s'' is invalid.', operationMode));

    % Call the user-supplied compute function
    dataOut = obj.classifierComputeFunction(obj, operationMode, obj.classifierParams, nullResponses, testResponses);

    % Set the trainedClassifier property for future predictions
    if (isfield(dataOut, 'trainedClassifier'))
        obj.trainedClassifier = dataOut.trainedClassifier;
        obj.preProcessingConstants = dataOut.preProcessingConstants;
    end
end