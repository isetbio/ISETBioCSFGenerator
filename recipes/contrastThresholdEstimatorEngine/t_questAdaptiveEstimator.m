%% setup variables and constants
projectName = 'ISETImagePipeline';
dataBaseDir = getpref(projectName, 'dataDir');

displayFile = 'CRT12BitDisplay.mat';
display = load(fullfile(dataBaseDir, displayFile));

retina = ConeResponse('eccBasedConeDensity', true, 'eccBasedConeQuantal', true, ...
    'fovealDegree', 1.0, 'display', display.CRT12BitDisplay, 'integrationTime', 0.05);
retina.visualizeMosaic();

% estimate on the log contrast domain
estDomain  = -5 : 0.025 : 0;
slopeRange = 0.1 : 0.5 : 50;

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

% Apply the default plotlab recipe 
% overriding just the figure size
plotlabOBJ.applyRecipe(...
  'figureWidthInches', 10, ...
  'figureHeightInches', 10);

% estimate on the log contrast domain
estDomain  = -5 : 0.05 : 0;
slopeRange = 1.0 : 0.5 : 50;

% spatialFreq = [1, 2, 4, 6, 8, 12, 16, 20, 30];
spatialFreq = [0.25, 0.5, 1, 2, 4, 6, 8, 12, 14];
threshold = zeros(1, length(spatialFreq));

figure();
for idx = 1:length(spatialFreq)
    threshold(idx) = QUEST(retina, display.CRT12BitDisplay, spatialFreq(idx), estDomain, slopeRange, idx);
    % threshold(idx) = adaptiveQUEST_mul(retina, display.CRT12BitDisplay, spatialFreq(idx), estDomain, slopeRange, idx);
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
function threshold = QUEST(retina, display, spatialFreq, estDomain, slopeRange, figIdx)

observer = PoissonTemplateObserver(retina, display, 'L+M+S', spatialFreq);

estimator = QuestThresholdEstimator('minTrial', 100, 'maxTrial', 1024, 'stopCriterion', 0.05, ...
    'estDomain', estDomain, 'numEstimator', 1, 'slopeRange', slopeRange);

[crst, flag] = estimator.nextStimulus();

% 10 trial for each contrast level
nRepeat = 10;
while (flag)
    
    % log contrast -> contrast
    stimCrst = 10 ^ crst;
    
    % code for scene engine, neural engine, and classifier engine
    % here we combined them in the Observer class object
    [~, response] = observer.multiTrial(stimCrst, nRepeat);
    response = response + 1;
    
    [crst, flag] = estimator.multiTrial(crst * ones(1, nRepeat), response);
end

% Show results
fprintf('%d trials recorded \n', estimator.nTrial);

subplot(3, 3, figIdx);
% title(sprintf('spatial frequency: %.2f cyc/deg \n', spatialFreq * 2));
[threshold, para] = estimator.thresholdMLE('showPlot', true, 'pointSize', 8);

fprintf('Maximum likelihood fit parameters: %0.2f, %0.2f, %0.2f, %0.2f\n', ...
    para(1), para(2), para(3), para(4));

end

% Trial-by-trial protocol
function threshold = adaptiveQUEST(retina, display, spatialFreq, estDomain, slopeRange, figIdx)

observer = PoissonTemplateObserver(retina, display, 'L+M+S', spatialFreq);

% Multiple QUEST+ object for adaptive procedure
estimator = QuestThresholdEstimator('minTrial', 64, 'maxTrial', 512, 'stopCriterion', 0.04, ...
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

fprintf('Maximum likelihood fit parameters: %0.2f, %0.2f, %0.2f, %0.2f\n', ...
    para(1), para(2), para(3), para(4));

end

% Run many trials as possible for each contrast level
function threshold = adaptiveQUEST_mul(retina, display, spatialFreq, estDomain, slopeRange, figIdx)

observer = PoissonTemplateObserver(retina, display, 'L+M+S', spatialFreq);

% Multiple QUEST+ object for adaptive procedure
estimator = QuestThresholdEstimator('minTrial', 256, 'maxTrial', 4096, 'stopCriterion', 0.025, ...
    'estDomain', estDomain, 'numEstimator', 3, 'slopeRange', slopeRange);

[crst, flag] = estimator.nextStimulus();
% 64 trial for each contrast level
nRepeat = 64;
while (flag)
    
    % log contrast -> contrast
    stimCrst = 10 ^ crst;
    
    % code for scene engine, neural engine, and classifier engine
    % here we combined them in the Observer class object
    [~, response] = observer.multiTrial(stimCrst, nRepeat);
    response = response + 1;
    
    [crst, flag] = estimator.multiTrial(crst * ones(1, nRepeat), response);
end

% Show results
fprintf('%d trials recorded \n', estimator.nTrial);

subplot(3, 3, figIdx);
% title(sprintf('spatial frequency: %.2f cyc/deg \n', spatialFreq * 2));
[threshold, para] = estimator.thresholdMLE('showPlot', true, 'pointSize', 20);

fprintf('Maximum likelihood fit parameters: %0.2f, %0.2f, %0.2f, %0.2f\n', ...
    para(1), para(2), para(3), para(4));

end
