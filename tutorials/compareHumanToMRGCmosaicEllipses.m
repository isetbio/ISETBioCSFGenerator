function compareHumanToMRGCmosaicEllipses()

    load('humanEllipses.mat')

    hFig = generateFigure(1, humanEllipses, 'cones_Poisson_rcePoisson_Summary.mat', 210, ' cones (Poisson noise, ideal observer) ');
    NicePlot.exportFigToPDF('conesPoisson.pdf', hFig, 300);


    %hFig = generateFigure(2, humanEllipses, 'cones_Gaussian_rceTemplateDistance_Summary.mat', 120, ' cones (Gaussian noise, template observer) ');
    %NicePlot.exportFigToPDF('conesGaussian.pdf', hFig, 300);

    hFig = generateFigure(3, humanEllipses, 'mRGCs_Gaussian_rceTemplateDistance_Summary.mat', 0.75, ' mRGCs (Gaussian noise, template observer) ');
    NicePlot.exportFigToPDF('mRGCsGaussian.pdf', hFig, 300);


end

function hFig=generateFigure(figureNo, humanEllipses, mRGCmosaicEllipsesFileName, mRGCmosaicEllipseScalar, mRGCmosaicLabel)

    hFig = figure(figureNo); clf;
    set(hFig, 'Position', [10 10 1400 1400], 'Color', [1 1 1]);
    ax = subplot('Position', [0.107 0.09, 0.89 0.89]);

    % Draw axes
    plot(ax,[-1 1], [0 0], 'k-');
    hold(ax, 'on');
    plot(ax,[0 0], [-1 1],  'k-');

    % Load the mRGCmosaic ellipses 
    load(mRGCmosaicEllipsesFileName, 'allData');
    theReferenceContrastLabels = keys(allData);

    % Draw
    drawIndividualPoints = false;
    
    for iRef = 1:numel(theReferenceContrastLabels)
        d = allData(theReferenceContrastLabels{iRef});
        if (drawIndividualPoints)
            x = d.referenceLMconeContrast(1)+d.thresholdDeltaConeContrasts(1,:);
            y = d.referenceLMconeContrast(2)+d.thresholdDeltaConeContrasts(2,:);
            meanX = mean(x);
            meanY = mean(y);
            x = meanX + (x-meanX) * mRGCmosaicEllipseScalar;
            y = meanY + (y-meanY) * mRGCmosaicEllipseScalar;
            plot(ax,x, y, 'ro', 'MarkerSize', 8, 'LineWidth', 1.0, 'MarkerEdgeColor', 0.5*[0.7 0.5 0.5], 'MarkerFaceColor', [0.7 0.5 0.5]);
        end

        x = d.fittedEllipse(:,1);
        y = d.fittedEllipse(:,2);
        meanX = mean(x);
        meanY = mean(y);
        d.fittedEllipse(:,1) = meanX + (x-meanX) * mRGCmosaicEllipseScalar;
        d.fittedEllipse(:,2) = meanY + (y-meanY) * mRGCmosaicEllipseScalar;

        p1 = patch(...
            'Faces',1:size(d.fittedEllipse,1), ...
            'Vertices',d.fittedEllipse, ...
            'FaceColor', [1 0.5 0.5], ...
            'EdgeColor', [1 0 0], ...
            'LineWidth', 1.5, ...
            'FaceAlpha', 0.5, ...
            'EdgeAlpha', 1.0);
    end
    clear 'd'
    for iRef = 1:size(humanEllipses,2)/2
        xData = humanEllipses(:,(iRef-1)*2+1);
        yData = humanEllipses(:,(iRef-1)*2+2);
        xData(end+1) = xData(1);
        yData(end+1) = yData(1);
        d.fittedEllipse(:,1) = xData(:);
        d.fittedEllipse(:,2) = yData(:);
        p2 = patch(...
            'Faces',1:size(d.fittedEllipse,1), ...
            'Vertices',d.fittedEllipse, ...
            'FaceColor', [0.8 0.8 0.8], ...
            'EdgeColor', [0.5 0.5 0.5], ...
            'LineWidth', 1.5, ...
            'FaceAlpha', 0.5, ...
            'EdgeAlpha', 1.0);
    end

    hL = legend(ax, [p1 p2], {mRGCmosaicLabel, ' human'}, 'Location', 'NorthWest', 'NumColumns',1);
    set(hL, 'EdgeColor', [0.9 0.9 0.9], 'Color', [0.95 0.95 0.95]);

    maxVisualizedContrast = 0.21;
    set(ax, 'FontSize', 40);
    set(ax, 'XLim', maxVisualizedContrast*[-1 1], 'YLim', maxVisualizedContrast*[-1 1], 'XTick', -1:0.05:1, 'YTick', -1:0.05:1);
    set(ax, 'LineWidth', 1.5, 'XColor', [0.3 0.3 0.3], 'YColor', [0.3 0.3 0.3])
    grid(ax, 'on')
    xlabel(ax, 'L cone contrast');
    ylabel(ax, 'M cone contrast');

end

% 
% % Load cone thresholds with Gaussian noise
% load('cones_Gaussian_rceTemplateDistance_Summary.mat', 'allData');
% theReferenceContrastLabels = keys(allData);
% 
% for iRef = 1:numel(theReferenceContrastLabels)
%     d = allData(theReferenceContrastLabels{iRef});
%     plot(ax,d.referenceLMconeContrast(1)+d.thresholdDeltaConeContrasts(1,:), ...
%          d.referenceLMconeContrast(2)+d.thresholdDeltaConeContrasts(2,:), 'ko', 'MarkerSize', 8, 'LineWidth', 1.0, 'MarkerFaceColor',  [1 1 0]);
%     p2 = plot(ax,d.fittedEllipse(:,1), d.fittedEllipse(:,2), 'r-', 'LineWidth', 1.0);
% end
% clear 'allData'
% 
% % Load mRGC thresholds with Gaussian noise
% load('mRGCs_Gaussian_rceTemplateDistance_Summary.mat')
% theReferenceContrastLabels = keys(allData);
% 
% for iRef = 1:numel(theReferenceContrastLabels)
%     d = allData(theReferenceContrastLabels{iRef});
%     plot(ax,d.referenceLMconeContrast(1)+d.thresholdDeltaConeContrasts(1,:), ...
%          d.referenceLMconeContrast(2)+d.thresholdDeltaConeContrasts(2,:), 'ko', 'MarkerSize', 8, 'LineWidth', 1.0, 'MarkerFaceColor',  [1 1 0]);
%     p3 = plot(ax,d.fittedEllipse(:,1), d.fittedEllipse(:,2), 'b-', 'LineWidth', 1.0);
% end
% clear 'allData'
% 
% legend(ax,[p1 p2 p3], {'cones (Poisson)', 'cones (Gaussian)', 'mRGCs (linear)'}, 'Location', 'NorthOutside', 'NumColumns',3);
% set(ax, 'FontSize', 20);
% set(ax, 'XLim', [-0.35 0.35], 'YLim', [-0.35 0.35], 'XTick', -1:0.1:1, 'YTick', -1:0.1:1);
% grid(ax, 'on')
% xlabel(ax, 'L cone contrast');
% ylabel(ax, 'M cone contrast');
% 
