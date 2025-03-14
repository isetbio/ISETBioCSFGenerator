function dataOut = rceTemplateDistanceTAFC(obj, operationMode, ~, theResponses, whichAlternatives)
% Compute function for ideal signal-known Poisson noise classifier for TAFC
%
% Syntax:
%     dataOut = rceTemplateDistanceTAFC(obj, operationMode, classifierParamsStruct, theResponses, whichAlternatives)
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
%    operationMode set to 'train', it sets up template (nearest neighbor)
%    classifier, using the mean of the passed null and test responses as a
%    template.  Pass noise free responses to get the signal known exactly
%    version.
%
%    When this is called from a parent @responseClassifierEngine object
%    with operationMode set to 'predict', it makes correct/incorrect
%    predictions for the passed instances. The predictions are made using a
%    minimum L2-norm distance to the templates for each alternative as its
%    decision rule.
%
%    This function models an N-way forced choice task with one stimulus
%    presented per trial, so that each trial is has responses to be one of the
%    N possible alternatives.
%
%    It can be used to handle TAFC if you stack the two alternatives into
%    one long response, which is done (e.g.) in computePerformance when the
%    TAFC flag is set.
%
%    Typically, training will be with one noise free instances of the null
%    and test responses, while testing will be on noisy instances.  For
%    training, however, the means of the passed null and test instances are
%    taken as the noise free template, so that you can in fact pass
%    multiple instances and they can be noisy.
%
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
% See also: computePerformance, t_spatialCSF

% History:
%   10/17/20  dhb  Wrote it from rcePoissonTAFC.
%   12/21/d4  dhb  Generalize away from TAFC.  We handle TAFC in the call now.

% For consistency with the interface
if (nargin == 0)
    dataOut = struct('Classifier', 'Template Distance Classifier Ideal');
    return;
elseif (nargin < 5 | isempty(whichAlternatives))
    whichAlternatives = [];
end

% Check operation mode
if (~strcmp(operationMode,'train') && ~strcmp(operationMode,'predict'))
    error('Unknown operation mode passed.  Must be ''train'' or ''predict''');
end

if (strcmp(operationMode, 'train'))
   % No noise response template for each alternative stimulus 
    %
    % This code stretches the responses x time out along the columns,
    % with each instances a row, and then takes the mean over rows.
    %
    % Usually we just pass one instance, but if you have lots of noisy
    % instances you can get the template as the mean of all of them.
    for ii = 1:length(theResponses)
        theTemplates{ii} = mean(theResponses{ii}(:, :), 1);
    end
    dataOut.trainedClassifier = [];
    dataOut.preProcessingConstants.theTemplates = theTemplates;
    
    return;
end

if (strcmp(operationMode, 'predict'))
    % Get templates we tucked away at training time.
    theTemplates = obj.preProcessingConstants.theTemplates;
    
    % Make sure number of instances matches across passed responses.
    nTrials = size(theResponses, 1);
    if (~isempty(whichAlternatives))
        assert(nTrials == length(whichAlternatives));
    end
    
    % Compute response finding which of the possible alternatives has the
    % highest likelihood.
    response = zeros(1, nTrials);
    whichAlternativesPredicted = zeros(1, nTrials);
    nAlternatives = length(theTemplates);
    trialDecisionLikelihoods = zeros(nTrials,nAlternatives);
    for tt = 1:nTrials
        % Find the decision likelihood
        for aa = 1:nAlternatives
            trialDecisionDistances(tt,aa) = norm(theResponses(tt, :) - theTemplates{aa});
        end

        % Pick the winner
        [~,whichAlternativeList] = min(trialDecisionDistances(tt,:));
        whichAlternativesPredicted(tt) = whichAlternativeList(1);

        % See if it was correct
        if (~isempty(whichAlternatives))
            if (whichAlternativesPredicted(tt) == whichAlternatives(tt))
                response(tt) = 1;
            else
                response(tt) = 0;
            end
        else
            response(tt) = NaN;
        end
    end
    
    % Set up return
    dataOut.trialPredictions = response;
    dataOut.pCorrect = mean(response);
    dataOut.whichAlternativePredicted = whichAlternativesPredicted;
    
    return;
end



end
