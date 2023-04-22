function [nextCrst, nextFlag] = singleTrialValidationBlocked(this, stim, response)

% Check that we're not doing method of constant stimuli, in which case we
% should not be here.
if (~this.validation)
    error('singleTrialValidationBlocked is only for validation method');
end
if (this.nRepeat <= 1)
    error('singleTrialValidationBlocked is only for nTest > 1');
end

% Convert from {0, 1} to {1, 2} response encoding for QUEST procedure
response = response + 1;

% Update the current QUEST object in use
this.estimators{this.estIdx} = qpUpdate(this.estimators{this.estIdx}, stim, response);

% Select next QUEST object
this.nTrial = this.nTrial + 1;
this.estIdx = mod(this.estIdx, this.numEstimator) + 1;

% Run the entire psychometric curve for validation mode
%
% With OLDWAY set to true, as it is here, each desired contrast
% is run once as a block. The code is pretty obscure, but we
% believe it is right.
OLDWAY = true;
if this.validation
    % The old way went through all of the contrasts in fixed order,
    % running all trials for one contrast before moving on to the next.
    % It is still here in case the new way breaks something.
    if (OLDWAY)
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
        error('This routine should not get here.');

        if (this.nTrial > length(this.validationTrialContrasts))
            this.testCrst = NaN;
            this.nextFlag = false;
        else
            this.testCrst = this.validationTrialContrasts(this.nTrial);
            this.nextFlag = true;
        end
    end
    
else
    error('This routine is only for the validation method')
    
end


end