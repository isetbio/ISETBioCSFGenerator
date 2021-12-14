function [threshold, para, dataOut] = thresholdMLE(this, varargin)
% Find maximum likelihood fit to psychometric data and return threshold estimate.
%
% Synopsis:
%    [threshold, para, dataOut] = thresholdMLE(this, varargin) 
%
% Description:
%    This method combines the data, fits the psychometric function, and
%    returns a threshold estimate.  Units of threshold are determined by
%    the PF specified in the object.
% 
%    A key/value pair determines the criterion proportion correct
%    threshold.  The default of 0.81606 corresponds to the performance
%    level in the mQUESTPlus parameterization of the Weibull, when the
%    lapse rate is 0 and the guess rate is 0.5.  Better to specify it
%    explicitly, however.
%
% Inputs:
%
% Outputs:
%    threshold              - Threshold estimate in the units of the
%                             psychometric function specified in the object.
%    para                   - Psychometric function parameters from fit.
%                             If these are passed via key/value pair then
%                             the passed values are returned.
%    dataOut                - Optional structure with more products of this
%                             routine. Use dataOut key/value pair for this
%                             to be returned.
%
% Optional key/value pairs.
%    'thresholdCriterion'   - Threshold fraction correct to which threshold
%                             should correspond.
%    'para'                 - Psychometric function parameters.  If passed,
%                             these are used and the fit is skipped. Must
%                             be matched to what the object's PF expects.
%                             Default empty, so that fitting is done.
%    'showPlot'             - Show a plot of data and fit. Default false.
%    'newFigure'            - Create a new figure window for the plot. 
%                             Otherwise the plot goes into the current
%                             figure. Default false.
%    'pointSize'            - Base point size for plot.  Default 25.
%    'returnedData'         - dataOut structure is only filled if this is
%                             true. Default false.

% History:
%  12/8/21: dhb  Add thresholdCriterion key/value pair, trying to keep
%                default behavior fixed.  Added skeleton header comment.
%           dhb  Move input parsing to the top, where I was looking for it.

% Parse inputs.
p = inputParser;
p.addParameter('showPlot',  false);
p.addParameter('newFigure', false);
p.addParameter('pointSize', 25);
p.addParameter('returnData',  false);
p.addParameter('para', []);
p.addParameter('thresholdCriterion',0.81606);
parse(p, varargin{:});

[stimVec, responseVec, structVec] = this.combineData();

questData = this.estimators{1};
psiParamsIndex = qpListMaxArg(questData.posterior);
psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);

% Make sure object and the quest data structure have the same PF.
if (~isequal(this.qpPF,questData.qpPF))
    error('Psychometric function not matched between object and QuestPlus data structure');
end

% Fit the PF to the data
para = qpFit(structVec, questData.qpPF, psiParamsQuest, questData.nOutcomes, ...
    'lowerBounds', [min(this.estDomain) min(this.slopeRange) min(this.guessRate) min(this.lapseRate)], ...
    'upperBounds', [max(this.estDomain) max(this.slopeRange) max(this.guessRate) max(this.lapseRate)]);

% Return the threshold.  The default, perhaps not well chosen, is to
% return the first parameter of the psychometric function.  If a criterion
% percent correct is set via the key/value pair 'thresholdCriterion', then
% threshold is returned as the stimulus that leads to this percent correct.
%
% The old way was to use the psychometric function mean as the threshold.
% This corresponds to proportion correct 0.81606 when guess rate is 0.5 and
% lapse rate is 0.  The new way explicitly finds the threshold
% corresponding to a passed proportion correct, with default 0.81606.
threshold = this.qpPFInv(p.Results.thresholdCriterion,para);

% Check that we get the same answer as the old way if enough else is
% matched to the way we were computing then.
if (p.Results.thresholdCriterion == 0.81606 && ...
        (isequal(this.qpPF,@qpPFWeibull) || isequal(this.qpPF,@qpPFWeibullLog) || isequal(this.qpPF,@qpPFStandardWeibull)))
    if (para(3) == 0.5 && para(4) == 0 )
        thresholdOld = para(1);
        if (abs(threshold-thresholdOld)/threshold > 1e-4)
            error('Not successfully producing old threshold with new method');
        end
    end
end

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
    fitCurve = this.qpPF(stimSpace', para);
    
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