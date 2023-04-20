function [nextCrst, nextFlag] = multiTrialFast(this, stimVec, responseVec)

% Convert from {0, 1} to {1, 2} response encoding for QUEST procedure
responseVec = responseVec + 1;

% Update for all the trials without updating entropy.
saveNoentropy = this.estimators{this.estIdx}.noentropy;
this.estimators{this.estIdx}.noentropy = 1;
for idx = 1:size(stimVec,1)
    this.estimators{this.estIdx} = qpUpdate(this.estimators{this.estIdx}, stimVec(idx), responseVec(idx));
end
this.estimators{this.estIdx}.noentropy = saveNoentropy;

% Update the entropy just once
this.estimators{this.estIdx} = qpUpdate(this.estimators{this.estIdx}, [], []);

% Select next QUEST object
this.nTrial = this.nTrial + size(stimVec,1);
this.estIdx = mod(this.estIdx, this.numEstimator) + 1;

% Check that we're done  
if this.nTrial >= this.maxTrial
    this.nextFlag = false;
end

% Running estimate of threshold and its standard error
[threshold, stderr] = thresholdEstimate(this);

% Stop only if stderr drop below criterion and we have at least minTrial # of trials
if ((this.stopCriterion(threshold, stderr)) && this.nTrial >= this.minTrial)
    this.nextFlag = false;
else
    this.nextFlag = true;
end

% Set next stimulus to ask for
this.testCrst = qpQuery(this.estimators{this.estIdx});

[nextCrst, nextFlag] = this.nextStimulus();
    
end
