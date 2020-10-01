function dataOut = poisson2AFC (~, ~, ~, nullResponses, testResponses)
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
    dataOut = struct('Classifier', 'Poisson 2AFC Ideal Observer');
    return;
end

% We simulate the observer with a "detection" protocol
% No noise response template for test/null stimulus
nullTemplate = nullResponses('none');
nullTemplate = nullTemplate(1, :);

testTemplate = testResponses('none');
testTemplate = testTemplate(1, :);

dataOut.nullTemplate = nullTemplate;
dataOut.testTemplate = testTemplate;

% Obtain trial-by-trial noisy response
nullResponses = nullResponses('random');
testResponses = testResponses('random');

nTrial = size(nullResponses, 1);
assert(nTrial == size(testResponses, 1));

% Compute response {0, 1} with log likelihood ratio
response = zeros(1, nTrial);
for idx = 1:nTrial    
    llhd_cr = llhd(nullResponses(idx, :), nullTemplate) + llhd(testResponses(idx, :), testTemplate);
    llhd_ic = llhd(nullResponses(idx, :), testTemplate) + llhd(testResponses(idx, :), nullTemplate);
    
    response(idx) = ((llhd_cr - llhd_ic) > 0);
end

dataOut.response = response;

end

% % Log-likelihood for Poisson R.V.
function ll = llhd(sample, rate)

ll = sum(sample .* log(rate) - rate);

end