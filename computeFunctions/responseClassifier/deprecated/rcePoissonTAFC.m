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

% For consistency with the interface
if (nargin == 0)
    warning(['This function has been deprecated. Consider using a more ',...
        'general function rcePoisson.m. ']);
    dataOut = struct('Classifier', 'Poisson Ideal Observer');
    return;
end

% Check operation mode
if (~strcmp(operationMode,'train') && ~strcmp(operationMode,'predict'))
    error('Unknown operation mode passed.  Must be ''train'' or ''predict''');
end

%reorganize the null and the test responses
if strcmp(operationMode, 'train')
    %concatenate them along the 3rd dimension (cones)
    cat1 = cat(3, testResponses, nullResponses);
    cat2 = cat(3, nullResponses, testResponses);
    %nicely put the concatenated cells back to the container
    theResponses = containers.Map(trainNoiseFlag, {cat1, cat2});
    dataOut = rcePoisson(obj, operationMode, classifierParamsStruct,...
        theResponses,[]);
else
    %combine the 2nd (time) and the 3rd dimensions (cones)
    outSampleTestStimResponses = testResponses(:,:);
    outSampleNullStimResponses = nullResponses(:,:);
    %concatenate them together along the 2nd dimension (cones)
    cat_resp = cat(2, outSampleTestStimResponses, outSampleNullStimResponses);
    %nicely put them back to the container
    theResponses = containers.Map(testNoiseFlag, cat_resp);
    %define the label for correct responses
    whichAlternatives = ones(size(theResponses,1),1);
    dataOut = rcePoisson(obj, operationMode, classifierParamsStruct,...
        theResponses,whichAlternatives);
end

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
