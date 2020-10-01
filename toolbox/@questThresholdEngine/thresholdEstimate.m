function [threshold, stderr] = thresholdEstimate(this)

estimates = this.parameterEstimate();

% Assume threshold parameter is the 1st one
threshold = mean(estimates(:, 1));
stderr = std(estimates(:, 1)) / sqrt(this.numEstimator);

end