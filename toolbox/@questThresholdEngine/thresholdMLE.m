% Combine all data to make a final MLE
function [threshold, para, dataOut] = thresholdMLE(this, varargin)
% Find maximum likelihood fit to psychometric data and return threshold estimate.
%
% Synopsis:
%    [threshold, para, dataOut] = thresholdMLE(this, varargin) 
%
% Description:
%
% Inputs:
%
% Outputs:
%
% Optional key/value pairs.
%    'thresholdCriterion'   - Threshold fraction correct to which threshold
%                             should correspond.

% History:
%  12/8/21: dhb  Add thresholdCriterion key/value pair, trying to keep
%                default behavior fixed.  Added skeleton header comment.
%           dhb  Move input parsing to the top, where I was looking for it.

% Parse inputs.
p = inputParser;
p.addParameter('showPlot',  false);
p.addParameter('newFigure', false);
p.addParameter('pointSize', 25);
p.addParameter('returnData',  false);\
parse(p, varargin{:});

[stimVec, responseVec, structVec] = this.combineData();

questData = this.estimators{1};
psiParamsIndex = qpListMaxArg(questData.posterior);
psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);

para = qpFit(structVec, questData.qpPF, psiParamsQuest, questData.nOutcomes, ...
    'lowerBounds', [min(this.estDomain) min(this.slopeRange) min(this.guessRate) min(this.lapseRate)], ...
    'upperBounds', [max(this.estDomain) max(this.slopeRange) max(this.guessRate) max(this.lapseRate)]);

% Return the thresholds.  The default, perhaps not well chosen, is to
% return the first parameter of the psychometric function.  If a criterion
% percent correct is set via the key/value pair 'thresholdCriterion', then
% threshold is returned as the stimulus that leads to this percent correct.
% stimContrast  = qpPFWeibullInv(proportionCorrect,psiParams)
% predictedProportions = qpPFWeibull(stimParams,psiParams)

threshold = para(1);



if ((p.Results.showPlot) || (p.Results.returnData))
    if p.Results.newFigure
        figure();
    end
    stimVal = unique(stimVec);
    pCorrect = zeros(1,length(stimVal));
    
    for idx = 1:length(stimVal)
        prop = responseVec(stimVec == stimVal(idx));
        pCorrect(idx) = sum(prop) / length(prop);
        
        if (p.Results.showPlot)
            scatter(stimVal(idx), pCorrect(idx), p.Results.pointSize * 100 / length(stimVec) * length(prop), ...
            'MarkerEdgeColor', zeros(1, 3), 'MarkerFaceColor', ones(1, 3) * 0.5, 'MarkerFaceAlpha', 0.5);
            hold on;
        end
    end
    
    stimSpace = this.estDomain;
    fitCurve = qpPFWeibull(stimSpace', para);
    
    if (p.Results.showPlot)
        plot(stimSpace, fitCurve(:, 2), '-', 'Color', [1.0 0.2 0.0], 'LineWidth', 3);
    
        xlim([min(this.estDomain), max(this.estDomain)]);
        xlabel('Contrast');

        ylim([0, 1]);
        ylabel('P Correct');
    end
    
    if (p.Results.returnData)
        dataOut.examinedContrasts = stimVal;
        dataOut.pCorrect = pCorrect;
        dataOut.examinedContrastsFit = stimSpace;
        dataOut.pCorrectFit = fitCurve(:,2);
    else 
        dataOut = [];
    end
    
end