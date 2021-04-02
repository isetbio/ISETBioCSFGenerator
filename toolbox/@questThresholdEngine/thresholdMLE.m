% Combine all data to make a final MLE
function [threshold, para, dataOut] = thresholdMLE(this, varargin)

[stimVec, responseVec, structVec] = this.combineData();

% Simple MLE for initial estimate
questData = this.estimators{1};
psiParamsIndex = qpListMaxArg(questData.posterior);
psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);

para = qpFit(structVec, questData.qpPF, psiParamsQuest, questData.nOutcomes, ...
    'lowerBounds', [min(this.estDomain) min(this.slopeRange) min(this.guessRate) min(this.lapseRate)], ...
    'upperBounds', [max(this.estDomain) max(this.slopeRange) max(this.guessRate) max(this.lapseRate)]);

threshold = para(1);

% If more complicated return is required
p = inputParser;
p.addParameter('showPlot',  false);
p.addParameter('newFigure', false);
p.addParameter('pointSize', 25);
p.addParameter('returnData',  false);
parse(p, varargin{:});

if ((p.Results.showPlot) || (p.Results.returnData))
    
    stimVal = sort(unique(stimVec));
    numPost = zeros(1,length(stimVal));
    numTotal = zeros(1,length(stimVal));
    
    pCorrect = zeros(1,length(stimVal));
    sePoint = zeros(1, length(stimVal));
    
    for idx = 1:length(stimVal)
        resp = responseVec(stimVec == stimVal(idx));
        numPost(idx) = sum(resp); numTotal(idx) = length(resp);
        
        % compute p correct
        pCorrStim = sum(resp) / length(resp);
        pCorrect(idx) = pCorrStim;
        
        % compute standard error for p of binomial distribution
        sePoint(idx) = sqrt(pCorrStim * (1 - pCorrStim) / length(resp));
    end
    
    % Fit psychometric curve with Palamedes
    % Use the Weibull function
    PF = @PAL_Logistic;
    
    %Threshold and Slope are free parameters, guess and lapse rate are fixed
    paramsFree = [1 1 0 0];
    
    searchGrid.alpha = linspace(min(this.estDomain), max(this.estDomain), 50);
    searchGrid.beta  = linspace(min(this.slopeRange), max(this.slopeRange), 50);
    searchGrid.gamma = 0.5;  %scalar here (since fixed) but may be vector
    searchGrid.lambda = 0.0;  %ditto
    
    para = PAL_PFML_Fit(stimVal, numPost, ...
        numTotal, searchGrid, paramsFree, PF);
    
    nRun = 1e3;
    paraSD = PAL_PFML_BootstrapNonParametric(...
        stimVal, numPost, numTotal, [], paramsFree, nRun, PF,...
        'searchGrid', searchGrid);
    
    if (p.Results.returnData)
        dataOut.paraSD = paraSD;
    end
end

if ((p.Results.showPlot) || (p.Results.returnData))
    if p.Results.newFigure
        figure();
    end
    
    if (p.Results.showPlot)
        for idx = 1:length(stimVal)            
            scatter(stimVal(idx), pCorrect(idx), p.Results.pointSize * 100 / length(stimVec) * numTotal(idx), ...
                'MarkerEdgeColor', zeros(1, 3), 'MarkerFaceColor', ones(1, 3) * 0.5, 'MarkerFaceAlpha', 0.5);
            hold on;
        end
        errorbar(stimVal, pCorrect, sePoint * 2, 'k', 'LineWidth', 1.0, 'LineStyle', 'none', 'Marker', 'none');
    end
    
    stimSpace = this.estDomain;
    fitCurve = PAL_Logistic(para, stimSpace');
    
    if (p.Results.showPlot)
        plot(stimSpace, fitCurve, '-', 'Color', [1.0 0.2 0.0], 'LineWidth', 3);
        
        xlim([min(this.estDomain), max(this.estDomain)]);
        xlabel('Contrast');
        
        ylim([0, 1]);
        ylabel('P Correct');
    end
    
    if (p.Results.returnData)
        dataOut.examinedContrasts = stimVal;
        dataOut.pCorrect = pCorrect;
        dataOut.sePoint = sePoint;
        dataOut.examinedContrastsFit = stimSpace;
        dataOut.pCorrectFit = fitCurve;
    else
        dataOut = [];
    end
    
end
