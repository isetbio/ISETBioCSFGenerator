function plotLMthresholds()
% UTTBSkip

    maxThreshold = 0.041;
    maxThreshold = 0.175;
    maxConeVisualizedThreshold = 0.25 * maxThreshold; 

    maxThreshold = [];
    maxConeVisualizedThreshold = [];

    theScriptName = 't_isoResponseLMplaneEllipses';
    
    mosaicEccDegs = [0 0];
    mosaicSizeDegs = [0.59 0.59];
    mRGCsNum = 4628;
    plotCones_vs_mRGCs(maxThreshold, maxConeVisualizedThreshold, theScriptName, mosaicEccDegs, mosaicSizeDegs, mRGCsNum);

    mosaicEccDegs = [-4 0];
    mosaicSizeDegs = [2.1 2.1];
    mRGCsNum = 4633;
    plotCones_vs_mRGCs(maxThreshold, maxConeVisualizedThreshold, theScriptName, mosaicEccDegs, mosaicSizeDegs, mRGCsNum);


    
    mosaicEccDegs = [-7.0 0];
    mosaicSizeDegs = [3.2 3.2];
    mRGCsNum = 4680;
    plotCones_vs_mRGCs(maxThreshold, maxConeVisualizedThreshold, theScriptName, mosaicEccDegs, mosaicSizeDegs, mRGCsNum);


    mosaicEccDegs = [-10.0 0];
    mosaicSizeDegs = [2.9 2.9];
    mRGCsNum = 2155;
    plotCones_vs_mRGCs(maxThreshold, maxConeVisualizedThreshold, theScriptName, mosaicEccDegs, mosaicSizeDegs, mRGCsNum);


    %maxThreshold = 0.041  * 4680/2195;
    mosaicEccDegs = [-14.0 0];
    mosaicSizeDegs = [4.1 4.1];
    mRGCsNum = 2195;
    plotCones_vs_mRGCs(maxThreshold, maxConeVisualizedThreshold, theScriptName, mosaicEccDegs, mosaicSizeDegs, mRGCsNum);



    mosaicEccDegs = [-19.0 0];
    mosaicSizeDegs = [6 6];
    mRGCsNum = 2146;
    plotCones_vs_mRGCs(maxThreshold, maxConeVisualizedThreshold, theScriptName, mosaicEccDegs, mosaicSizeDegs, mRGCsNum);


    mosaicEccDegs = [-25.0 0];
    mosaicSizeDegs = [5.8 5.8];
    mRGCsNum = 1069;
    plotCones_vs_mRGCs(maxThreshold, maxConeVisualizedThreshold, theScriptName, mosaicEccDegs, mosaicSizeDegs, mRGCsNum);



    mosaicEccDegs = [-32.0 0];
    mosaicSizeDegs = [9 9];
    mRGCsNum = 1090;
    plotCones_vs_mRGCs(maxThreshold, maxConeVisualizedThreshold, theScriptName, mosaicEccDegs, mosaicSizeDegs, mRGCsNum);


end


function plotCones_vs_mRGCs(maxVisualizedThreshold, maxConeVisualizedThreshold, theScriptName, mosaicEccDegs, mosaicSizeDegs, mRGCsNum)
    
    % Common params
    opticsType = 'loadComputeReadyRGCMosaic';   % {'loadComputeReadyRGCMosaic', 'adaptiveOptics6MM'} 
    stimulusSpatialFrequency = 0;
    orientationDegs = 90;
    presentationMode = 'counterphasemodulated';
    useMetaContrast = true;
    useConeContrast = true;
    useFixationalEMs = false;
    whichNoiseFreeNre = 'mRGCMosaic';
    whichNoisyInstanceNre = 'Gaussian';
    whichClassifierEngine = 'rceTemplateDistance';

    % Cone-based LM thresholds
    mRGCOutputSignalType = 'cones';             % {'cones', 'mRGCs'}
    [figureFileBaseDir, resultsFileBaseDir] = setupFigureDirectory(theScriptName, ...
        useMetaContrast,useConeContrast,useFixationalEMs,...
        whichNoiseFreeNre,whichNoisyInstanceNre,...
        whichClassifierEngine,mRGCOutputSignalType);

    thresholdsDataFileName = ...
        sprintf('%sCSF_SF_%2.2fCPD_Optics_%s_EccDegs_x%2.1f_%2.1f_SizeDegs_%2.1fx%2.1f_OriDegs_%2.0f_%s.mat', ...
        mRGCOutputSignalType, ...
        stimulusSpatialFrequency, ...
        opticsType, ...
        mosaicEccDegs(1), mosaicEccDegs(2), ...
        mosaicSizeDegs(1),mosaicSizeDegs(2), ...
        orientationDegs, ...
        presentationMode);

    load(fullfile(resultsFileBaseDir,thresholdsDataFileName), 'examinedDirectionsOnLMplane', 'thresholdContrasts');
    coneThresholds = thresholdContrasts;


    % mRGC-based LM thresholds
    mRGCOutputSignalType = 'mRGCs';
    [~, resultsFileBaseDir] = setupFigureDirectory(theScriptName, ...
        useMetaContrast,useConeContrast,useFixationalEMs,...
        whichNoiseFreeNre,whichNoisyInstanceNre,...
        whichClassifierEngine,mRGCOutputSignalType);
    thresholdsDataFileName = ...
        sprintf('%sCSF_SF_%2.2fCPD_Optics_%s_EccDegs_x%2.1f_%2.1f_SizeDegs_%2.1fx%2.1f_OriDegs_%2.0f_%s.mat', ...
        mRGCOutputSignalType, ...
        stimulusSpatialFrequency, ...
        opticsType, ...
        mosaicEccDegs(1), mosaicEccDegs(2), ...
        mosaicSizeDegs(1),mosaicSizeDegs(2), ...
        orientationDegs, ...
        presentationMode);
    fullfile(resultsFileBaseDir,thresholdsDataFileName)
    load(fullfile(resultsFileBaseDir,thresholdsDataFileName), 'examinedDirectionsOnLMplane', 'thresholdContrasts');
    mRGCThresholds = thresholdContrasts;


    % Visualize the LM thresholds to be compared
    visualizeComparedLMplaneThresholds(...
            examinedDirectionsOnLMplane, ...
            coneThresholds, mRGCThresholds, ...
            'cones (L+M+S)', 'mRGCs (L+M+S)', ...
            maxConeVisualizedThreshold, maxVisualizedThreshold, ...
            mosaicSizeDegs, mRGCsNum, ...
            fullfile(figureFileBaseDir, sprintf('mRGCs_vs_cones_ecc_%2.0f_degs.pdf', mosaicEccDegs(1))));

 end


function  visualizeComparedLMplaneThresholds(...
            examinedDirectionsOnLMplane, ...
            thresholds1, thresholds2, ...
            legend1, legend2, ...
            maxVisualizedThreshold, maxVisualizedThreshold2, ...
            mosaicSizeDegs, mRGCsNum, ...
            pdfFileName)

    if (isempty(maxVisualizedThreshold)) && (isempty(maxVisualizedThreshold2))
        maxVisualizedThreshold = 1.05;
        maxVisualizedThreshold2 = 1.05;

        idx = find(examinedDirectionsOnLMplane == 90);
        thresholds1 = thresholds1 / thresholds1(idx);
        thresholds2 = thresholds2 / thresholds2(idx);

        xyTicks =  -2:1:2;
        xyTickLabels = {'-2', '-1', '0', '+1', '+2'};
        xyLims = 1.75*[-1 1];
        coneThresholdMarkerSize = 20;
        mRGCThresholdMarkerSize  = 20;

        thresholdGainToAccountForMaxVisualizedThreshold2 = 1.0;
 
        xAxisLabel = 'normalized L-cone threshold';
        yAxisLabel = 'normalized M-cone threshold';
    else  
        if (isempty(maxVisualizedThreshold))
            maxVisualizedThreshold = 1.05*max([max(thresholds1) max(thresholds2)]);
            if (maxVisualizedThreshold < 0.055)
                maxVisualizedThreshold = 0.055;
            end
        end
    
        if (isempty(maxVisualizedThreshold2))
            maxVisualizedThreshold2 = 1.05*max([max(thresholds1) max(thresholds2)]);
            if (maxVisualizedThreshold2 < 0.055)
                maxVisualizedThreshold2 = 0.055;
            end
        end

        xyTicks =  -.20:0.05:0.20;
        xyTickLabels = {'-.20', '-.15', '-.10', '-.05',  '0', '+.05', '+.10', '+.15', '+.20'};
        xyLims = maxVisualizedThreshold2*[-1 1];
        coneThresholdMarkerSize = 10;
        mRGCThresholdMarkerSize  = 20;
        thresholdGainToAccountForMaxVisualizedThreshold2 = maxVisualizedThreshold2/maxVisualizedThreshold;

        xAxisLabel = 'threshold contrast (L-cone)';
        yAxisLabel = 'threshold contrast (M-cone)';
    end


    
    hFig = figure(2); clf;
    set(hFig, 'Color', [1 1 1]);
    ff = PublicationReadyPlotLib.figureComponents('1x1 standard tall figure');
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);

    % Plot the axes
    plot(theAxes{1,1},  maxVisualizedThreshold2*[-1 1], [0 0], 'k-', 'LineWidth', 1.0);
    hold(theAxes{1,1}, 'on');
    plot(theAxes{1,1},  [0 0], maxVisualizedThreshold2*[-1 1], 'k-', 'LineWidth', 1.0);

 
    x1 = cosd(examinedDirectionsOnLMplane) .* thresholds1 * thresholdGainToAccountForMaxVisualizedThreshold2;
    y1 = sind(examinedDirectionsOnLMplane) .* thresholds1 * thresholdGainToAccountForMaxVisualizedThreshold2;

    x2 = cosd(examinedDirectionsOnLMplane) .* thresholds2;
    y2 = sind(examinedDirectionsOnLMplane) .* thresholds2;
    
    if (numel(examinedDirectionsOnLMplane)>6)
        [z, a, b, rotationRadians] = fitEllipseToXYpoints(...
            [x1(:) y1(:)]', ...
            'nonLinear', false);
        %[z, a, b, rotationRadians] = fitellipse(thresholdConeContrasts, 'linear');
    
        % Plot the fitted ellipse
        % form the parameter vector
        npts = 100;
        t = linspace(0, 2*pi, npts);
    
        % Rotation matrix
        Q = [cos(rotationRadians), -sin(rotationRadians); sin(rotationRadians) cos(rotationRadians)];
        % Ellipse points
        X = Q * [a * cos(t); b * sin(t)] + repmat(z, 1, npts);

        % PLot the ellipse points
        plot(theAxes{1,1}, X(1,:), X(2,:), 'k-', 'LineWidth', 6.0);
        plot(theAxes{1,1}, X(1,:), X(2,:), 'c-', 'Color', [1 0.5 0.7], 'LineWidth',2.0);


        [z, a, b, rotationRadians] = fitEllipseToXYpoints(...
            [x2(:) y2(:)]', ...
            'nonLinear', false);
        %[z, a, b, rotationRadians] = fitellipse(thresholdConeContrasts, 'linear');
    
        % Plot the fitted ellipse
        % form the parameter vector
        npts = 100;
        t = linspace(0, 2*pi, npts);
    
        % Rotation matrix
        Q = [cos(rotationRadians), -sin(rotationRadians); sin(rotationRadians) cos(rotationRadians)];
        % Ellipse points
        X = Q * [a * cos(t); b * sin(t)] + repmat(z, 1, npts);

        % PLot the ellipse points
        plot(theAxes{1,1}, X(1,:), X(2,:), 'k-', 'LineWidth', 6.0);
        plot(theAxes{1,1}, X(1,:), X(2,:), 'c-', 'Color', [0.5 1 0.85], 'LineWidth',2.0);

    end

    % LM plane thresholds-1
    
    p1 = scatter(theAxes{1,1}, x1,y1, coneThresholdMarkerSize^2, ...
        'MarkerFaceColor', [1 0.5 0.7], 'MarkerEdgeColor', [1 0.5 0.7]*0.5, ...
        'MarkerFaceAlpha', 0.6, 'LineWidth', 1.5);
    hold(theAxes{1,1}, 'on');


    % LM plane thresholds-2
    
    p2 = scatter(theAxes{1,1}, x2,y2, mRGCThresholdMarkerSize^2, ...
        'MarkerFaceColor', [0.5 1 0.85], 'MarkerEdgeColor', [0.5 1 0.85]*0.5, ...
        'MarkerFaceAlpha', 0.6, 'LineWidth', 1.5);

   

    hold(theAxes{1,1}, 'off');
    set(theAxes{1,1}, 'XLim', xyLims, 'YLim', xyLims);
    axis(theAxes{1,1}, 'square');
    box(theAxes{1,1}, 'off');
    set(theAxes{1,1}, 'Color', 'none', 'Box', 'off', 'XColor', [0.1 0.1 0.1], 'YColor', [0.1 0.1 0.1]);
    set(theAxes{1,1}, 'XTick', xyTicks, 'XTickLabel', xyTickLabels)
    set(theAxes{1,1}, 'YTick', xyTicks, 'YTickLabel', xyTickLabels)
    legend(theAxes{1,1}, [p1 p2], {'input cone mosaic', sprintf('mRGC mosaic (%2.1f^o\\times%2.1f^o, N=%d)', mosaicSizeDegs(1), mosaicSizeDegs(2), mRGCsNum)}, 'Location', 'NorthOutside', 'NumColumns', 2)
    %PublicationReadyPlotLib.offsetAxes(theAxes{1,1},ff, xLims, yLims);
    PublicationReadyPlotLib.labelAxes(theAxes{1,1},ff, xAxisLabel, yAxisLabel);
    PublicationReadyPlotLib.applyFormat(theAxes{1,1},ff);
    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);

    fprintf('PDF saved in %s\n', pdfFileName);
    %set(theThresholdAxes, 'FontSize', 24);

end

function [figureFileBaseDir, resultsFileBaseDir] = setupFigureDirectory(theScriptName, ...
    useMetaContrast,useConeContrast,useFixationalEMs,...
    whichNoiseFreeNre,whichNoisyInstanceNre,...
    whichClassifierEngine,mRGCOutputSignalType)

    % Make sure local/figures directory exists so we can write out our figures in peace
    projectBaseDir = ISETBioCSFGeneratorRootPath;
    if (~exist(fullfile(projectBaseDir,'local', theScriptName, 'figures'),'dir'))
        mkdir(fullfile(projectBaseDir,'local', theScriptName, 'figures'));
        fprintf('Generated figure directory at %s\n', fullfile(projectBaseDir,'local', theScriptName, 'figures'))
    end

    if (~exist(fullfile(projectBaseDir,'local', theScriptName, 'results'),'dir'))
        mkdir(fullfile(projectBaseDir,'local', theScriptName, 'results'));
        fprintf('Generated results directory at %s\n', fullfile(projectBaseDir,'local', theScriptName, 'figures'))
    end


    figureFileBaseDir = fullfile(projectBaseDir,'local',theScriptName,'figures', ...
        sprintf('%s_Meta_%d_ConeContrast_%d_FEMs_%d_%s_%s_%s_%s', theScriptName, ...
        useMetaContrast,useConeContrast,useFixationalEMs,whichNoiseFreeNre,whichNoisyInstanceNre,...
        whichClassifierEngine,mRGCOutputSignalType));

    if (~exist(figureFileBaseDir, 'dir'))
        mkdir(figureFileBaseDir);
        fprintf('Generated figure sub-directory at %s\n', figureFileBaseDir);
    end

    resultsFileBaseDir = fullfile(projectBaseDir,'local',theScriptName,'results', ...
        sprintf('%s_Meta_%d_ConeContrast_%d_FEMs_%d_%s_%s_%s_%s', theScriptName, ...
        useMetaContrast,useConeContrast,useFixationalEMs,whichNoiseFreeNre,whichNoisyInstanceNre,...
        whichClassifierEngine,mRGCOutputSignalType));
    
    if (~exist(resultsFileBaseDir, 'dir'))
        mkdir(resultsFileBaseDir);
        fprintf('Generated figure sub-directory at %s\n', resultsFileBaseDir);
    end

end
