function dataOut = rcePoissonNWayFC(obj, operationMode, classifierParamsStruct, theResponses, whichAlternatives)
% Compute function for ideal signal-known Poisson noise classifier for N-way forced choice.
%
% Syntax:
%     dataOut = rcePoissonNWayFC(obj, operationMode,classifierParamsStruct, theReponses, whichAlternatives)
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
%    sets up the Poission log-likehood ratio (LL) for the decision, using
%    the mean of the passed N responses as a template.  Pass
%    noise free responses to get the signal known exactly version.
%
%    When this is called from a parent @responseClassifierEngine object
%    with operationMode set to 'predict', it makes correct/incorrect
%    predictions for the passed instances. The predictions are made using a
%    maximum likelihood classifier for N-way forced choice and Poisson noise.
%
%    This function models an N-way forced choice task with one stimulus
%    presented per trial, so that each trial is considered to be one of the
%    four passed alternatives.
%
%    Typically, training will be with one noise free instances four alternatives,
%    while testing will be on noisy instances.  For training, however, the means
%    of the passed instances are taken as the noise free template, so that you
%    can in fact pass multiple instances per alternative and they can be noisy.
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
%    theResponses             - an N-dimensional cell array, with each entry
%                               a [mInstances x nDims x tTimePoints] matrix
%                               of responses for one of the N alternatives.
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
% See also: t_thresholdEngine, t_responseClassifer, rcePoissonTAFC,
%           PoissonDecisionLogLikelihood,
%           PoissonIdealObserverNAlternativeFC.
%            

% History:
%   12/03/21  dhb  Started on this.

% For consistency with the interface
if (nargin == 0)
    dataOut = struct('Classifier', 'Poisson N-way Forced Choice Ideal Observer');
    return;
end

% Check operation mode
if (~strcmp(operationMode,'train') && ~strcmp(operationMode,'predict'))
    error('Unknown operation mode passed.  Must be ''train'' or ''predict''');
end
    
if (strcmp(operationMode, 'train'))  
    % No noise response template for each alternative stimulus 
    for ii = 1:length(theResponses)
        theTemplates{ii} = mean(theResponses{ii}(:, :), 1);
    end
    dataOut.trainedClassifier = [];
    dataOut.preProcessingConstants = struct('theTemplates', theTemplate);
    
    return;
end

if (strcmp(operationMode, 'predict'))
    % Get templates we tucked away at training time.
    theTemplates = obj.preProcessingConstants.theTemplates;
    
    % Make sure number of instances matches across passed responses.
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
        llhdCr = llhd(nullResponses(idx, :), nullTemplate) + PoissonDecisionLogLikelihoood(testResponses(idx, :), testTemplate);
        llhdIc = llhd(nullResponses(idx, :), testTemplate) + PoissonDecisionLogLikelihoood(testResponses(idx, :), nullTemplate);
        
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


