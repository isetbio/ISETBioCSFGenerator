function dataOut = poissonTemplateClassifier(obj, operationMode, ~, nullResponses, testResponses)
% Compute function for training a binary SVM classifier to predict classes from
% responses.

% Description:
%    Compute function to be used as a computeFunctionHandle for a @responseClassifierEngine
%    object.
%
%    When called from a parent @responseClassifierEngine object, it
%    computes the Poission log-likehood ratio (LL) of test stimulus, using the
%    noise-free response as a template. A response (class label) is then
%    generated as (LL > 0).
%
% Inputs:
%    This function has five input argument s.t. it is consistent with the
%    responseClassifierEngine interface. It actually only uses
%    'nullResponses' and 'testResponses', each contain the noise-free
%     template and noisy trial-by-trial instances.
%
%    nullResponses                  - an [mTrials x nDims] matrix of responses to the null stimulus
%
%    testResponses                  - an [mTrials x nDims] matrix of responses to the test stimulus


% For consistency with the interface
if (nargin == 0)
    dataOut = struct('Classifier', 'Poisson Template Matching');
    return;
end


if (strcmp(operationMode, 'train'))
    
    % We simulate the observer with a "detection" protocol
    % No noise response template for test/null stimulus
    
    nullTemplate = nullResponses(1, :);
    testTemplate = testResponses(1, :);
    
    dataOut.trainedClassifier = 'N/A';
    dataOut.preProcessingConstants = struct('nullTemplate', nullTemplate, 'testTemplate', testTemplate);
    
    return;
end

if (strcmp(operationMode, 'predict'))
        
    testTemplate = obj.preProcessingConstants.testTemplate;
    
    nTrial = size(nullResponses, 1);
    assert(nTrial == size(testResponses, 1));
    
    % Compute response {0, 1} with log likelihood ratio
    response = zeros(1, nTrial);
    for idx = 1:nTrial
        ll = logLikelihood(testTemplate, testResponses(idx, :), nullResponses(idx, :));
        response(idx) = (ll > 0);
    end
    
    dataOut.response = response;
    dataOut.pCorrect = mean(response);
    
    return;
end

end

% Log-likelihood ratio for Poisson R.V.
function ll = logLikelihood(rate, test, null)

zeroThreshold = 1e-20;
test = test + zeroThreshold;
null = null + zeroThreshold;

ll = rate .* log(test ./ null) + null - test;
ll = sum(ll);

end
