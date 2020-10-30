function dataOut = rceTemplateTAFC(obj, operationMode, ~, nullResponses, testResponses)
% Compute function for ideal signal-known Poisson noise classifier for TAFC
%
% Syntax:
%     dataOut = rceTemplateTAFC(obj, operationMode, classifierParamsStruct, nullResponses, testResponses)
%
% Description:
%    Compute function to be used as a computeFunctionHandle for a
%    @responseClassifierEngine object.
%
%    When called with no arguments, this returns a default parameters
%    struct for itself. Currently, there are no parameters used for this
%    routine, so a struct with an unintersting field is returned just to
%    keep the calling machinery happy.
%
%    When called from a parent @responseClassifierEngine object with
%    operationMOde set to 'train', it sets up template (nearest neighbor)
%    classifier , using the mean of the passed null and test responses as a
%    template.  Pass noise free responses to get the signal known exactly
%    version.
%
%    When this is called from a parent @responseClassifierEngine object
%    with operationMode set to 'predict', it makes correct/incorrect
%    predictions for the passed instances. The predictions are made using a
%    the template classifier.
%
%    This function models a TAFC task, so that each trial is considered to
%    be either null/test across the two alternatives, or test/null.
%
%    Typically, training will be with one noise free instances of the null
%    and test responses, while testing will be on noisy instances.  For
%    training, however, the means of the passed null and test instances are
%    taken as the noise free template, so that you can in fact pass
%    multiple instances and they can be noisy.
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
%    nullResponses            - an [mInstances x nDims x tTimePoints] matrix of responses to the null stimulus
%    testResponses            - an [mInstances x nDims x tTimePoints] matrix of responses to the test stimulus
%
% Outputs:
%    dataOut                  - If function called with no input arguments,
%                               this is the default parameter structure for
%                               this funciton.
%                              
%                               If function is called from a @responseClassifierEngine
%                               object, what this contains depends on
%                               whether operationMode is 'train' or
%                               'predict'.
% 
%                               For 'train', the struct has two fields as
%                               required by responseClassifierEngine
%                               objects.
%                                 .trainedClassifier: Set to empty matrix
%                                 .preProcessingConstants: Structure with
%                                   templates to be used by classifier to
%                                   predict.
%             
%                              For 'predict', the struct has two fields as
%                               required by responseClassifierEngine
%                               objects.
%                                 .trialPredictions : Vector with trial-by-trial
%                                                     correct/incorrect
%                                                     predictions. 1 means
%                                                     correct, 0 means incorrect.
%                                 .pCorrect         : Probability correct, which is
%                                                     just the mean of the trialPredictions
%                                                     field.
%
% Optional key/value pairs:
%   None.
%
% See also: t_thresholdEngine, t_responseClassifer.

% History:
%   10/17/20  dhb  Wrote it from rcePoissonTAFC.


% For consistency with the interface
if (nargin == 0)
    dataOut = struct('Classifier', 'Template 2AFC Classifier');
    return;
end

if (strcmp(operationMode, 'train'))
    
    % Store template for test/null stimulus
    nullTemplate = mean(nullResponses(:, :), 1);
    testTemplate = mean(testResponses(:, :), 1);
    
    dataOut.trainedClassifier = [];
    dataOut.preProcessingConstants = struct('nullTemplate', nullTemplate, 'testTemplate', testTemplate);
    
    return;
end

if (strcmp(operationMode, 'predict'))
    % Get template we tucked away at training time.
    nullTemplate = obj.preProcessingConstants.nullTemplate;
    testTemplate = obj.preProcessingConstants.testTemplate;
    nullTestTemplate = [nullTemplate ; testTemplate];
    testNullTemplate = [testTemplate ; nullTemplate];
    
    % Make sure number of null and test instances matches.
    nTrials = size(nullResponses, 1);
    assert(nTrials == size(testResponses, 1));
    
    % Compute response {0, 1} with nearest neighbor classification. 
    % Evaluate on assumption that stimulus order is null-test.  
    response = zeros(1, nTrials);
    for idx = 1:nTrials
        nullTestResponse = [nullResponses(idx, :) ; testResponses(idx, :)];
        distance_cr = norm(nullTestResponse - nullTestTemplate);
        distance_ic = norm(nullTestResponse - testNullTemplate);
        
        % Correct if distance of null-test response to null-test template
        % is less than distance to test-null template.
        response(idx) = ((distance_cr - distance_ic) < 0);     
    end
    
    % Set up return
    dataOut.trialPredictions = response;
    dataOut.pCorrect = mean(response);
    
    return;
end

end
