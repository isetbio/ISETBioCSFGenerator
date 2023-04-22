function [nextCrst, nextFlag] = singleTrial(this, stim, response)

% Check that we're not doing method of constant stimuli with
% multiple trials, in which case we should not be here.
if (this.validation & this.nRepeat > 1)
    error('This routine is not for validation with nTest > 1');
end

% Convert from {0, 1} to {1, 2} response encoding for QUEST procedure
response = response + 1;

% Update the current QUEST object in use
this.estimators{this.estIdx} = qpUpdate(this.estimators{this.estIdx}, stim, response);

% Select next QUEST object
this.nTrial = this.nTrial + 1;
this.estIdx = mod(this.estIdx, this.numEstimator) + 1;

% Run the entire psychometric curve for validation mode.
%
% With OLDWAY set to false, this goes through blocks of contrast
% on trial at a time, in random order.  Good for real psychophysical
% experiments, and some simulations.
OLDWAY = false;
if this.validation
    % The old way went through all of the contrasts in fixed order,
    % running all trials for one contrast before moving on to the next.
    % It is still here in case the new way breaks something.
    if (OLDWAY)
        error('Should not be here');

        crstIdx = floor(this.nTrial / this.nRepeat) + 1;
        
        if crstIdx > length(this.estDomain)
            this.testCrst = NaN;
            this.nextFlag = false;
        else
            this.testCrst = this.estDomain(crstIdx);
            this.nextFlag = true;
        end
    
    % The new way draws from a pre-randomized order set up when
    % we created the questThresholdObject.  The two ways should
    % give the same result for a computational observer, but the
    % new way is better for a human psychophysical experiment.
    else
        if (this.nTrial > length(this.validationTrialContrasts))
            this.testCrst = NaN;
            this.nextFlag = false;
        else
            this.testCrst = this.validationTrialContrasts(this.nTrial);
            this.nextFlag = true;
        end
    end
    
else
    
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


end