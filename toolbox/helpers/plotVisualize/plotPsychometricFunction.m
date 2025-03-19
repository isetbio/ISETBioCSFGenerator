function [parameterValuesExamined,pCorrect] = plotPsychometricFunction(questObj, threshold, ...
    fittedPsychometricParams, thresholdParameters, pdfFileName, varargin)

    p = inputParser;
    p.addParameter('xRange', []);
    p.parse(varargin{:});
    xRange = p.Results.xRange;

    % Retrieve the parameter values that Quest examined performance
    [stimVec, responseVec] = questObj.combineData();
    stimVal = unique(stimVec);
    parameterValuesExamined = 10.^(stimVal)*thresholdParameters.maxParamValue;

    % Retrieve the measured performance
    pCorrect = zeros(1,numel(stimVal));
    for idx = 1:numel(stimVal)
        prop = responseVec(stimVec == stimVal(idx));
        pCorrect(idx) = sum(prop) / length(prop);
    end

    % Retrieve the fitted psychometric function
    fittedPsychometricFunction = questObj.qpPF(questObj.estDomain', fittedPsychometricParams);
    examinedParameterAxis = 10.^(questObj.estDomain)*thresholdParameters.maxParamValue;

    hFig = figure(100); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 400 350]);

    % Plot the fitted psychometric function
    plot(examinedParameterAxis, fittedPsychometricFunction(:, 2), '-', 'Color', [1.0 0.2 0.0], 'LineWidth', 3);
    hold on;
    % Plot the measured psychometric curve data points
    scatter(parameterValuesExamined, pCorrect, 14 * 100 / length(stimVec) * length(prop), ...
        'filled', 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [1 0. 0.], 'MarkerFaceAlpha', 0.5, 'LineWidth', 1.5);

    % Plot the computed threshold
    if (~isempty(xRange))
        plot([xRange(1) threshold], thresholdParameters.thresholdCriterion*[1 1], 'k--', 'LineWidth', 1.5);
    else
        plot([parameterValuesExamined(1) threshold], thresholdParameters.thresholdCriterion*[1 1], 'k--', 'LineWidth', 1.5);
    end
    plot(threshold*[1 1], [0 thresholdParameters.thresholdCriterion], 'k--', 'LineWidth', 1.5);
    plot(threshold, 0.03, 'kv', 'MarkerSize', 12, 'LineWidth', 1.5, 'MarkerFaceColor', [1 0.9 0.9], 'MarkerEdgeColor', [0 0. 0.]);
    
    % Finalize plot
    if (length(examinedParameterAxis) > 1)
        xlim([min(examinedParameterAxis), max(examinedParameterAxis)]);
    end
    ylim([0, 1]);
    xlabel('size (degs)');
    ylabel('Pcorrect');

    set(gca, 'FontSize', 16);
    set(gca, 'YTick', 0:0.1:1, 'XScale', 'log');

    if (~isempty(xRange))
        set(gca, 'XLim', xRange);
    end

    grid on; box off
    title(sprintf('threshold: %2.3f degs', threshold));

    if (~isempty(pdfFileName))
        NicePlot.exportFigToPDF(pdfFileName,hFig, 300);
    end
end
