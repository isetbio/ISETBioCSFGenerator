%% setup variables and constants
projectName = 'ISETImagePipeline';
dataBaseDir = getpref(projectName, 'dataDir');

displayFile = 'CRT12BitDisplay.mat';
display = load(fullfile(dataBaseDir, displayFile));

retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 1.0, 'display', display.CRT12BitDisplay, 'integrationTime', 0.05);
retina.visualizeMosaic();

% estimate on the log contrast domain
estDomain  = -6 : 0.025 : 0;
slopeRange = 0.1 : 1 : 100;

observer = PoissonTemplateObserver(retina, display.CRT12BitDisplay, 'L+M+S', 1);

%% Single QUEST+ object with fixed number of trials
estimator = QuestThresholdEstimator('minTrial', 128, 'maxTrial', 256, ...
    'estDomain', estDomain, 'numEstimator', 1, 'slopeRange', slopeRange);

[crst, flag] = estimator.nextStimulus();
while (flag)
    
    % log contrast -> contrast
    stimCrst = 10 ^ crst;
    
    % code for scene engine, neural engine, and classifier engine
    % here we combined them in the Observer class object
    response = observer.singleTrial(stimCrst) + 1;
    
    [crst, flag] = estimator.singleTrial(crst, response);
    [threshold, ~] = estimator.thresholdEstimate();
    
    fprintf('Trial count: %d, threshold estimate: %.3f \n', estimator.nTrial, threshold);
end

% Show results
[threshold, stderr] = estimator.thresholdEstimate();
fprintf('%d trials recorded, (log) threshold estimate: %.2f +/- %.2f \n', estimator.nTrial, threshold, stderr);

[threshold, para] = estimator.thresholdMLE('showPlot', true);
fprintf('Maximum likelihood fit parameters: %0.1f, %0.1f, %0.1f, %0.2f\n', ...
    para(1), para(2), para(3), para(4));

%% Single QUEST+ object but when it's easy to run multiple trials per contrast
estimator = QuestThresholdEstimator('minTrial', 1024, 'maxTrial', 1024, ...
    'estDomain', estDomain, 'numEstimator', 1, 'slopeRange', slopeRange);

[crst, flag] = estimator.nextStimulus();

% 32 trial for each contrast level
nRepeat = 32;
while (flag)
    
    % log contrast -> contrast
    stimCrst = 10 ^ crst;
    
    % code for scene engine, neural engine, and classifier engine
    % here we combined them in the Observer class object
    [~, response] = observer.multiTrial(stimCrst, nRepeat);
    response = response + 1;
    
    [crst, flag] = estimator.multiTrial(crst * ones(1, nRepeat), response);
    [threshold, ~] = estimator.thresholdEstimate();
    
    fprintf('Trial count: %d, threshold estimate: %.3f \n', estimator.nTrial, threshold);
end

% Show results
[threshold, stderr] = estimator.thresholdEstimate();
fprintf('%d trials recorded, (log) threshold estimate: %.2f +/- %.2f \n', estimator.nTrial, threshold, stderr);

[threshold, para] = estimator.thresholdMLE('showPlot', true);
fprintf('Maximum likelihood fit parameters: %0.1f, %0.1f, %0.1f, %0.2f\n', ...
    para(1), para(2), para(3), para(4));

%% Multiple QUEST+ object for adaptive procedure
estimator = QuestThresholdEstimator('minTrial', 32, 'maxTrial', 4096, 'stopCriterion', 0.05, ...
    'estDomain', estDomain, 'numEstimator', 3, 'slopeRange', slopeRange);

[crst, flag] = estimator.nextStimulus();
while (flag)
    
    % log contrast -> contrast
    stimCrst = 10 ^ crst;
    
    % code for scene engine, neural engine, and classifier engine
    % here we combined them in the Observer class object
    response = observer.singleTrial(stimCrst) + 1;
    
    [crst, flag] = estimator.singleTrial(crst, response);
    [threshold, stderr] = estimator.thresholdEstimate();
    
    fprintf('Trial count: %d, threshold estimate: %.3f, stderr: %.3f \n', estimator.nTrial, threshold, stderr);
end

% Show results
[threshold, stderr] = estimator.thresholdEstimate();
fprintf('%d trials recorded, (log) threshold estimate: %.2f +/- %.2f \n', estimator.nTrial, threshold, stderr);

[threshold, para] = estimator.thresholdMLE('showPlot', true);
fprintf('Maximum likelihood fit parameters: %0.1f, %0.1f, %0.1f, %0.2f\n', ...
    para(1), para(2), para(3), para(4));

%% Multiple QUEST+ object for adaptive procedure
estimator = QuestThresholdEstimator('minTrial', 256, 'maxTrial', 4096, 'stopCriterion', 0.01, ...
    'estDomain', estDomain, 'numEstimator', 3, 'slopeRange', slopeRange);

[crst, flag] = estimator.nextStimulus();

% 32 trial for each contrast level
nRepeat = 32;
while (flag)
    
    % log contrast -> contrast
    stimCrst = 10 ^ crst;
    
    % code for scene engine, neural engine, and classifier engine
    % here we combined them in the Observer class object
    [~, response] = observer.multiTrial(stimCrst, nRepeat);
    response = response + 1;
    
    [crst, flag] = estimator.multiTrial(crst * ones(1, nRepeat), response);
    [threshold, stderr] = estimator.thresholdEstimate();
    
    fprintf('Trial count: %d, threshold estimate: %.3f, stderr: %.3f \n', estimator.nTrial, threshold, stderr);
end

% Show results
[threshold, stderr] = estimator.thresholdEstimate();
fprintf('%d trials recorded, (log) threshold estimate: %.2f +/- %.2f \n', estimator.nTrial, threshold, stderr);

[threshold, para] = estimator.thresholdMLE('showPlot', true);
fprintf('Maximum likelihood fit parameters: %0.1f, %0.1f, %0.1f, %0.2f\n', ...
    para(1), para(2), para(3), para(4));


%% CSF
plotlabOBJ = plotlab();
%    
% Apply the default plotlab recipe 
% overriding just the figure size
plotlabOBJ.applyRecipe(...
  'figureWidthInches', 10, ...
  'figureHeightInches', 10);

spatialFreq = [1, 2, 4, 8, 12, 16, 20, 30, 40];
threshold = zeros(1, length(spatialFreq));

figure();
for idx = 1:length(spatialFreq)
    threshold(idx) = adaptiveQUEST(retina, display.CRT12BitDisplay, spatialFreq(idx), estDomain, slopeRange, idx);
end

%% Plot CSF
figure();
crstThreshold = 10 .^ threshold;
plot(spatialFreq * 2, log(1 ./ crstThreshold), '-ok', 'LineWidth', 2);
box off;

ytickVal = yticks();
yticklabels(floor(exp(ytickVal)));

xlabel('Spatial Frequency');
ylabel('Sensitivity');

%% helper function
function threshold = adaptiveQUEST(retina, display, spatialFreq, estDomain, slopeRange, figIdx)

observer = PoissonTemplateObserver(retina, display, 'L+M+S', spatialFreq);

% Multiple QUEST+ object for adaptive procedure
estimator = QuestThresholdEstimator('minTrial', 80, 'maxTrial', 256, 'stopCriterion', 0.05, ...
    'estDomain', estDomain, 'numEstimator', 3, 'slopeRange', slopeRange);

[nextCrst, nextFlag] = estimator.nextStimulus();
while (nextFlag)
    
    % log contrast -> contrast
    stimCrst = 10 ^ nextCrst;
    
    % code for scene engine, neural engine, and classifier engine
    % here we combined them in the Observer class object
    response = observer.singleTrial(stimCrst) + 1;
    
    [nextCrst, nextFlag] = estimator.singleTrial(nextCrst, response);
end

% Show results
fprintf('%d trials recorded \n', estimator.nTrial);

subplot(3, 3, figIdx);
% title(sprintf('spatial frequency: %.2f cyc/deg \n', spatialFreq * 2));
[threshold, para] = estimator.thresholdMLE('showPlot', true);

fprintf('Maximum likelihood fit parameters: %0.1f, %0.1f, %0.1f, %0.2f\n', ...
    para(1), para(2), para(3), para(4));

end

function threshold = adaptiveQUEST_mul(retina, display, spatialFreq, estDomain, slopeRange, figIdx)

observer = PoissonTemplateObserver(retina, display, 'L+M+S', spatialFreq);

% Multiple QUEST+ object for adaptive procedure
estimator = QuestThresholdEstimator('minTrial', 256, 'maxTrial', 4096, 'stopCriterion', 0.025, ...
    'estDomain', estDomain, 'numEstimator', 3, 'slopeRange', slopeRange);

[crst, flag] = estimator.nextStimulus();
% 32 trial for each contrast level
nRepeat = 32;
while (flag)
    
    % log contrast -> contrast
    stimCrst = 10 ^ crst;
    
    % code for scene engine, neural engine, and classifier engine
    % here we combined them in the Observer class object
    [~, response] = observer.multiTrial(stimCrst, nRepeat);
    response = response + 1;
    
    [crst, flag] = estimator.multiTrial(crst * ones(1, nRepeat), response);
    
    % [threshold, stderr] = estimator.thresholdEstimate();
    % fprintf('Trial count: %d, threshold estimate: %.3f, stderr: %.3f \n', estimator.nTrial, threshold, stderr);
end

% Show results
fprintf('%d trials recorded \n', estimator.nTrial);

subplot(3, 3, figIdx);
% title(sprintf('spatial frequency: %.2f cyc/deg \n', spatialFreq * 2));
[threshold, para] = estimator.thresholdMLE('showPlot', true);

fprintf('Maximum likelihood fit parameters: %0.1f, %0.1f, %0.1f, %0.2f\n', ...
    para(1), para(2), para(3), para(4));

end
