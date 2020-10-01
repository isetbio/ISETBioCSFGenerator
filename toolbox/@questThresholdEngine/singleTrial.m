function [nextCrst, nextFlag] = singleTrial(this, stim, response)

% Convert from {0, 1} to {1, 2} response encoding for QUEST procedure
response = response + 1;

% Update the current QUEST object in use
this.estimators{this.estIdx} = qpUpdate(this.estimators{this.estIdx}, stim, response);

% Select next QUEST object
this.nTrial = this.nTrial + 1;
this.estIdx = mod(this.estIdx, this.numEstimator) + 1;

if this.nTrial >= this.maxTrial
    this.nextFlag = false;
end

% Running estimate of threshold and its standard error
[threshold, stderr] = thresholdEstimate(this);

% Stop only if  stderr drop below criterion and we have at least minTrial # of trials
if ((stderr / abs(threshold)) < this.stopCriterion) && this.nTrial >= this.minTrial
    this.nextFlag = false;
end

% Set next stimulus to ask for
this.testCrst = qpQuery(this.estimators{this.estIdx});

[nextCrst, nextFlag] = this.nextStimulus();

end