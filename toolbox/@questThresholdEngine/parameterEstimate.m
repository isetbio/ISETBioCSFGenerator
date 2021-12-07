function estimates = parameterEstimate(this)

% For now, hard code the number of parameters to be 4
estimates = zeros(this.numEstimator, 4);

for idx = 1 : this.numEstimator
    estimator = this.estimators{idx};
    psiParamsIndex = qpListMaxArg(estimator.posterior);
    estimates(idx, :) = estimator.psiParamsDomain(psiParamsIndex, :);
end

end