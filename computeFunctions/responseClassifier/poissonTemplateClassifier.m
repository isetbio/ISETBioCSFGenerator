function dataOut = poissonTemplateClassifier(~, ~, ~, nullResponses, testResponses)

if (nargin == 0)
    dataOut = struct('Classifier', 'Poisson Template Ideal Observer');
    return;
end

% no noise response template for null stimulus
nullTemplate = nullResponses('none');
nullTemplate = nullTemplate(1, :);

% no noise response template for test stimulus
testTemplate = testResponses('none');
testTemplate = testTemplate(1, :);

dataOut.null = nullTemplate;
dataOut.test = testTemplate;

end