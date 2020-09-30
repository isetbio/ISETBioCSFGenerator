function dataOut = poissonTemplateClassifier(~, ~, ~, nullResponses, testResponses)

if (nargin == 0)
    dataOut = struct('Classifier', 'Poisson Template Ideal Observer');
    return;
end

% we simulate the observer with a "detection" protocol
% no noise response template for test/null stimulus
nullTemplate = nullResponses('none');
nullTemplate = nullTemplate(1, :);

testTemplate = testResponses('none');
testTemplate = testTemplate(1, :);

dataOut.nullTemplate = nullTemplate;
dataOut.testTemplate = testTemplate;

% trial-by-trial noisy response
nullResponses = nullResponses('random');
testResponses = testResponses('random');

nTrial = size(nullResponses, 1);
assert(nTrial == size(testResponses, 1));

% compute response {0, 1} with log likelihood ratio
response = zeros(1, nTrial);
for idx = 1:nTrial
    ll = logLikelihood(testTemplate, testResponses(idx, :), nullResponses(idx, :));
    response(idx) = (ll > 0);
end

dataOut.response = response;

end

function ll = logLikelihood(rate, test, null)
ll = rate .* log(test ./ null) + null - test;
ll = sum(ll);
end