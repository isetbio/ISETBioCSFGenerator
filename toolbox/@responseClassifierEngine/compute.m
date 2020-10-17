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
%    operationMode                  - a string, either 'train' or 'predict'
%                                     defining the operation to be performed: either
%                                     training the classifier on a (nullResponses, testResponses) 
%                                     training data set, or using the trained classifier to
%                                     derive predictions on a previously unseen (out-of-sample)
%                                     (nullResponses, testResponses) data set
%    nullResponses                  - an [mTrials x nDims x nTimePoints] matrix of responses to the null stimulus
%    testResponses                  - an [mTrials x nDims x nTimePoints] matrix of responses to the test stimulus
%
% Optional key/value input arguments: none
%
% Outputs:
%    dataOut  - If called from a parent @responseClassifierEngine object, the returned
%               struct is organized as follows: 
%
%               In 'train' mode, the struct must have the following two
%               fields.  It may can additional fields that are classifier
%               dependent.
%                   .trainedClassifier       : the trained binary SCV classifer
%                   .preProcessingConstants  : constants computed during the dimensionality reduction preprocessing phase
%
%               In 'predict' mode, the struct must have the following two
%               fields.  It may can additional fields that are classifier
%               dependent.
%                   .pCorrect                : probability of correct classification for the out-of-sample trials (testing data)
%                   .trialPredictions        : vector of the trial-by-trial predictions 
%                                              (0 == incorrectly predicting nominal class, 1 == correctly predicting nominal class)
% See Also:
%     t_responseClassifier

% History:
%    9/20/2020  NPC Wrote it

    % Validate the operationMode
    assert(ismember(operationMode, obj.validOperationModes), sprintf('The passed responseClassifierEngine.compute() ''%s'' is invalid.', operationMode));

    % Call the user-supplied compute function
    dataOut = obj.classifierComputeFunction(obj, operationMode, obj.classifierParams, nullResponses, testResponses);

    % Retrieve the returned dataOut fields
    theFieldNames = fieldnames(dataOut);

    % Set the trainedClassifier property for future predictions
    if (strcmp(operationMode, 'train'))
        % Validate the dataOut struct
        for k = 1:numel(obj.requiredFieldsForTrainDataOutStruct)
            assert(ismember(obj.requiredFieldsForTrainDataOutStruct{k}, theFieldNames), sprintf('dataOut struct does not contain the ''%s'' field', obj.requiredFieldsForTrainDataOutStruct{k}));
        end
        
        % Store the trained classfier
        obj.trainedClassifier = dataOut.trainedClassifier;
        
        % Store the pre-processing constants
        obj.preProcessingConstants = dataOut.preProcessingConstants;
    else
        % Validate the dataOut struct
        for k = 1:numel(obj.requiredFieldsForPredictDataOutStruct)
            assert(ismember(obj.requiredFieldsForPredictDataOutStruct{k}, theFieldNames), sprintf('dataOut struct does not contain the ''%s'' field', obj.requiredFieldsForPredictDataOutStruct{k}));
        end
    end
end