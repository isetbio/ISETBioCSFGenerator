function [nextCrst, nextFlag] = multiTrial(this, stimVec, responseVec)

for idx = 1:length(stimVec)
    this.singleTrial(stimVec(idx), responseVec(idx));
end

[nextCrst, nextFlag] = this.nextStimulus();

end