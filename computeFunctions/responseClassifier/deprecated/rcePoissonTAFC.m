function dataOut = rcePoissonTAFC(obj, operationMode, classifierParamsStruct, nullResponses, testResponses)
% Compute function for ideal signal-known Poisson noise classifier for
% TAFC. This function has been deprecated because there is a more general 
% function rcePoisson.m that is applicable for both N-way forced choice 
% and TAFC tasks.
%
% History:
%   10/17/20  dhb  Lots of work on comments.
%   12/03/21  dhb  Pull out Poisson log likelihood as function into
%                  isetbio.
%   05/10/23  fh   Edited the code to direct people to the more general
%                   function rcePoisson.m and moved this to the deprecated
%                   folder

if (nargin == 0)
    warning(['This function has been deprecated. Consider using a more ',...
        'general function rcePoisson.m. ']);
    dataOut = struct('Classifier', 'Poisson Forced Choice Ideal Observer');
    return;
end

%reorganize the null and the test responses
if strcmp(operationMode, 'train')
    %This method / function will only be called in computePerformance.m,
    %which already reorganizes the data by concatenating the cone responses
    %given the test and the null stimuli. So when this function is called
    %inside computePerformance.m, passed 'nullResponses' is actually the
    %concatenated test + null responses.
    theResponses = nullResponses;
    dataOut = rcePoisson(obj, operationMode, classifierParamsStruct,...
        theResponses,[]);
else
    %Analogously, when this function is called inside computePerformance.m
    %and it's 'predict' mode, then passed 'nullResponses' is actually the
    %concatenated test + null responses, and 'testResponses' is actually
    %the label (correct answer)
    theResponses = nullResponses;
    whichAlternatives = testResponses;
    dataOut = rcePoisson(obj, operationMode, classifierParamsStruct,...
        theResponses,whichAlternatives);
end

% if (nargin == 0)
%     dataOut = struct('Classifier', 'Poisson Forced Choice Ideal Observer');
%     return;
% end
% 
% % Check operation mode
% if (~strcmp(operationMode,'train') && ~strcmp(operationMode,'predict'))
%     error('Unknown operation mode passed.  Must be ''train'' or ''predict''');
% end
% 
% if (strcmp(operationMode, 'train'))  
%     % We simulate the observer with a "detection" protocol
%     % No noise response template for test/null stimulus 
%     nullTemplate = mean(nullResponses(:, :), 1);
%     testTemplate = mean(testResponses(:, :), 1);
%     dataOut.trainedClassifier = [];
%     dataOut.preProcessingConstants = struct('nullTemplate', nullTemplate, 'testTemplate', testTemplate);
% 
%     return;
% end
% 
% if (strcmp(operationMode, 'predict'))
%     % Get template we tucked away at training time.
%     nullTemplate = obj.preProcessingConstants.nullTemplate;
%     testTemplate = obj.preProcessingConstants.testTemplate;
% 
%     % Make sure number of null and test instances matches.
%     nTrials = size(nullResponses, 1);
%     assert(nTrials == size(testResponses, 1));
% 
%     % Compute response {0, 1} with log likelihood ratio.  Assume without
%     % loss of generality that alternative order on each trial is null/test,
%     % compute likelihood of each order, and call it correct if likelihood
%     % for null/test order is higher.  Because each entry of the response
%     % vectors is an indendent Poisson observation, we can just sum the log
%     % likelihood of the null and test components of the TAFC response.
%     response = zeros(1, nTrials);
%     for idx = 1:nTrials
%         llhdCr = PoissonDecisionLogLikelihoood(nullResponses(idx, :), nullTemplate) + PoissonDecisionLogLikelihoood(testResponses(idx, :), testTemplate);
%         llhdIc = PoissonDecisionLogLikelihoood(nullResponses(idx, :), testTemplate) + PoissonDecisionLogLikelihoood(testResponses(idx, :), nullTemplate);
% 
%         % For likelihood ratio extremely close to 1, do a coin flip
%         threshold = 1e-10;
%         if (abs(llhdCr - llhdIc) <= threshold)
%             response(idx) = (rand() > 0.5);
%         else
%             response(idx) = ((llhdCr - llhdIc) > 0);
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
