function dataOut = rcePoissonTAFC(obj, operationMode, classifierParamsStruct, nullResponses, testResponses)
% Compute function for ideal signal-known Poisson noise classifier for TAFC
%
% Syntax:
%     dataOut = rcePoissonTAFC(obj, operationMode, classifierParamsStruct, nullResponses, testResponses)
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
%    When called from a parent @responseClassifierEngine object with operationMode set to 'predict', it
%    sets up the Poission log-likehood ratio (LL) of test stimulus, using
%    the mean of the passed null and test responses as a template.  Pass
%    noise free responses to get the signal known exactly version.
%
%    When this is called from a parent @responseClassifierEngine object
%    with operationMode set to 'predict', it makes correct/incorrect
%    predictions for the passed instances. The predictions are made using a
%    maximum likelihood classifier for TAFC and Poisson noise.
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
%   10/17/20  dhb  Lots of work on comments.

% For consistency with the interface
if (nargin == 0)
    dataOut = struct('Classifier', 'Poisson 2AFC Ideal Observer');
    return;
end

% Check operation mode
if (~strcmp(operationMode,'train') && ~strcmp(operationMode,'predict'))
    error('Unknown operation mode passed.  Must be ''train'' or ''predict''');
end
    
if (strcmp(operationMode, 'train'))  
    % We simulate the observer with a "detection" protocol
    % No noise response template for test/null stimulus 
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
    
    % Make sure number of null and test instances matches.
    nTrials = size(nullResponses, 1);
    assert(nTrials == size(testResponses, 1));
    
    % Compute response {0, 1} with log likelihood ratio.  Assume alternative order on
    % each trial is null/test, compute likelihood of each order, and call
    % it correct if likelihood for null/test order is higher.  Because each
    % entry of the response vectors is an indendent Poisson observation, we
    % can just sum the log likelihood of the null and test components of
    % the TAFC response.
    response = zeros(1, nTrials);
    for idx = 1:nTrials
        llhdCr = llhd(nullResponses(idx, :), nullTemplate) + llhd(testResponses(idx, :), testTemplate);
        llhdIc = llhd(nullResponses(idx, :), testTemplate) + llhd(testResponses(idx, :), nullTemplate);
        
        % For likelihood ratio extremely close to 1, do a coin flip
        threshold = 1e-10;
        if (abs(llhdCr - llhdIc) <= threshold)
            response(idx) = (rand() > 0.5);
        else
            response(idx) = ((llhdCr - llhdIc) > 0);
        end
    end
    
    % Set up return
    dataOut.trialPredictions = response;
    dataOut.pCorrect = mean(response);
    
    return;
end

end

% Log-likelihood for Poisson response vectors 
%
% Syntax:
%    ll = llhd(response, template)
%
% Description:
%    Compute log likelihood of a vector of independent
%    Poisson observations, given Poisson mean for each.
%    See for example:
%      https://online.stat.psu.edu/stat504/node/28/
%
% Inputs:
%    response  -   Vector of observed responses.
%    template  -   Poisson mean for each entry of response.
%
% Outputs:
%    ll        -   Log likelihood.

% History:
%    10/17/20  dhb  Added comments.

function ll = llhd(response, template)
    response = double(response);
    template = double(template);
    ll = sum(response .* log(template) - template);
end
