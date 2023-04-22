function [nextCrst, nextFlag] = multiTrial(this, stimVec, responseVec)

% Check whether we should be calling multiTrialFast instead.
if (~this.validation & this.nRepeat > 1)
    error('You want to be calling the multiTrialQuestBlocked method instead of multiTrial');
end

for idx = 1:length(stimVec)
    if (this.validation & this.nRepeat > 1)
        this.singleTrialValidationBlocked(stimVec(idx), responseVec(idx));
    else
        this.singleTrial(stimVec(idx), responseVec(idx));
    end
end

[nextCrst, nextFlag] = this.nextStimulus();

end