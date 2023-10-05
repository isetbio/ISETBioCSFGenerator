function dataOut = rcePoissonNWay_OneStimulusPerTrial(obj, operationMode, classifierParamsStruct, theResponses, whichAlternatives)
% Compute function for ideal signal-known Poisson noise classifier for 
% N-way forced choice. This function has been deprecated because there is a
% more general function rcePoisson.m that is applicable for both N-way
% forced choice and TAFC tasks.
%
% History:
%   12/03/21  dhb  Started on this.
%   05/10/23  fh   Edited the code to direct people to the more general
%                   function rcePoisson.m and moved this to the deprecated
%                   folder

% For consistency with the interface
if (nargin == 0)
    warning(['This function has been deprecated. Consider using a more ',...
        'general function rcePoisson.m. ']);
    dataOut = struct('Classifier', 'Poisson Ideal Observer');
    return;
end

dataOut = rcePoisson(obj, operationMode, classifierParamsStruct,...
    theResponses, whichAlternatives);

% % Check operation mode
% if (~strcmp(operationMode,'train') && ~strcmp(operationMode,'predict'))
%     error('Unknown operation mode passed.  Must be ''train'' or ''predict''');
% end
% 
% if (strcmp(operationMode, 'train'))  
%     % No noise response template for each alternative stimulus 
%     for ii = 1:length(theResponses)
%         theTemplates{ii} = mean(theResponses{ii}(:, :), 1);
%     end
%     dataOut.trainedClassifier = [];
%     dataOut.preProcessingConstants.theTemplates = theTemplates;
% 
%     return;
% end
% 
% if (strcmp(operationMode, 'predict'))
%     % Get templates we tucked away at training time.
%     theTemplates = obj.preProcessingConstants.theTemplates;
% 
%     % Make sure number of instances matches across passed responses.
%     nTrials = size(theResponses, 1);
%     assert(nTrials == length(whichAlternatives));
% 
%     % Compute response finding which of the possible alternatives has the
%     % highest likelihood.
%     response = zeros(1, nTrials);
%     nAlternatives = length(theTemplates);
%     trialDecisionLikelihoods = zeros(nTrials,nAlternatives);
%     for tt = 1:nTrials
%         % Find the decision likelihood
%         for aa = 1:nAlternatives
%             trialDecisionLikelihoods(tt,aa) = PoissonDecisionLogLikelihoood(theResponses(tt, :), theTemplates{aa});
%         end
% 
%         % Pick the winner
%         [~,whichAlternativeList] = max(trialDecisionLikelihoods(tt,:));
%         whichAlternative = whichAlternativeList(1);
% 
%         % See if it was correct
%         if (whichAlternative == whichAlternatives(tt))
%             response(tt) = 1;
%         else
%             response(tt) = 0;
%         end
%     end
% 
%     % Set up return
%     dataOut.trialPredictions = response;
%     dataOut.pCorrect = mean(response);
% 
%     return;
% end
% 
% end
% 
% 
