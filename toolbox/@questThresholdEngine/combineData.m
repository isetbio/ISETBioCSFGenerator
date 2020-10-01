function [stimVec, responseVec, structVec] = combineData(this)

structVec = [];
for idx = 1 : this.numEstimator
    questData = this.estimators{idx};
    structVec = [structVec; questData.trialData];
end

nTotal = length(structVec);

stimVec = zeros(1, nTotal);
responseVec = zeros(1, nTotal);

for idx = 1:nTotal
    stimVec(idx) = structVec(idx).stim;
    responseVec(idx) = structVec(idx).outcome - 1;
end

end