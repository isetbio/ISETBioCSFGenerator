% Combine all data to make a final MLE
function [threshold, para, dataOut] = thresholdMLE(this, varargin)

[stimVec, responseVec, structVec] = this.combineData();

questData = this.estimators{1};
psiParamsIndex = qpListMaxArg(questData.posterior);
psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);

para = qpFit(structVec, questData.qpPF, psiParamsQuest, questData.nOutcomes, ...
    'lowerBounds', [min(this.estDomain) min(this.slopeRange) min(this.guessRate) min(this.lapseRate)], ...
    'upperBounds', [max(this.estDomain) max(this.slopeRange) max(this.guessRate) max(this.lapseRate)]);

threshold = para(1);

p = inputParser;
p.addParameter('showPlot',  false);
p.addParameter('newFigure', false);
p.addParameter('pointSize', 25);
p.addParameter('returnData',  false);
parse(p, varargin{:});

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