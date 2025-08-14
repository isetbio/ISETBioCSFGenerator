function plotCSF()

    if (1==1)

        mosaicEccDegs = [0 0];
        mosaicSizeDegs = [0.6 0.6];

        maxCSFvisualized = 100;
        maxSFvisualizedCones = 90;
        maxSFvisualizedAO = 200;

        
        colorCones = [1 0.5 0.7];
        colorMRGCPhysiologicalOpticsLuminance = [0.5 1 0.85];
        plotCones_vs_mRGCs(maxCSFvisualized, maxSFvisualized, maxSFvisualizedCones, mosaicEccDegs, mosaicSizeDegs, colorCones, colorMRGCPhysiologicalOpticsLuminance);

        colorMRGCPhysiologicalOpticsRedGreen = [1 0.0 0.7];
        opticsType = 'loadComputeReadyRGCMosaic'; 
        plotAchromatic_vs_LMopponent(maxCSFvisualized, maxSFvisualized, opticsType, mosaicEccDegs, mosaicSizeDegs, colorMRGCPhysiologicalOpticsLuminance, colorMRGCPhysiologicalOpticsRedGreen);
        
        
        colorMRGCAdaptiveOpticsLuminance = 0.5*colorMRGCPhysiologicalOpticsLuminance;
        colorMRGCAdaptiveOpticsRedGreen = 0.5*colorMRGCPhysiologicalOpticsRedGreen;
        opticsType = 'adaptiveOptics6MM';
        plotAchromatic_vs_LMopponent(maxCSFvisualized, maxSFvisualized, opticsType, mosaicEccDegs, mosaicSizeDegs, colorMRGCAdaptiveOpticsLuminance, colorMRGCAdaptiveOpticsRedGreen);
        


        stimulusChroma = 'luminance';
        plotPhysiologicalOptics_vs_adaptiveOptics(maxCSFvisualized, stimulusChroma, maxSFvisualized, maxSFvisualizedAO, mosaicEccDegs, mosaicSizeDegs, colorMRGCPhysiologicalOpticsLuminance, colorMRGCAdaptiveOpticsLuminance);
        
        color2 = [1 0.0 0.0];
        color1 = [0.5 1 0.85];
        stimulusChroma = 'red-green';
        plotPhysiologicalOptics_vs_adaptiveOptics(maxCSFvisualized, stimulusChroma, maxSFvisualized, maxSFvisualizedAO, mosaicEccDegs, mosaicSizeDegs, colorMRGCPhysiologicalOpticsRedGreen, colorMRGCAdaptiveOpticsRedGreen);
        

        pause;

        mosaicEccDegs = [-4 0];
        mosaicSizeDegs = [2.1 2.1];
        maxSFvisualized = 40;
        maxSFvisualizedCones = 50;
        color1 = [1 0.5 0.7];
        color2 = [0.5 1 0.85];
        plotCones_vs_mRGCs(maxCSFvisualized, maxSFvisualized, maxSFvisualizedCones, mosaicEccDegs, mosaicSizeDegs, color1, color2);
    

        

        mosaicEccDegs = [-7 0];
        mosaicSizeDegs = [3.2 3.2];
        maxSFvisualized = 25;
        maxSFvisualizedCones = 50;
        color1 = [1 0.5 0.7];
        color2 = [0.5 1 0.85];
        plotCones_vs_mRGCs(maxCSFvisualized, maxSFvisualized, maxSFvisualizedCones, mosaicEccDegs, mosaicSizeDegs, color1, color2);
    


    else

        mosaicEccDegs = [0 0];
        mosaicSizeDegs = [2 2];
    
        maxSFvisualized = 100;
        maxSFvisualizedAO = 200;
    
    
        color1 = [1 0.5 0.7];
        color2 = [0.5 1 0.85];
    
        
        plotCones_vs_mRGCs(200, maxSFvisualized, [], mosaicEccDegs, mosaicSizeDegs, color1, color2);
        plotAchromatic_vs_LMopponent(200, maxSFvisualized, mosaicEccDegs, mosaicSizeDegs, color1, color2);
        
        
        stimulusChroma = 'luminance';
        plotPhysiologicalOptics_vs_adaptiveOptics(200, stimulusChroma, maxSFvisualized, maxSFvisualizedAO, mosaicEccDegs, mosaicSizeDegs, color1, color2);
        
        stimulusChroma = 'red-green';
        plotPhysiologicalOptics_vs_adaptiveOptics(200, stimulusChroma, maxSFvisualized, maxSFvisualizedAO, mosaicEccDegs, mosaicSizeDegs, color1, color2);
    end 

    %plotAchromatic_vs_LMopponentAO(70);
    %plotFoveal_vs_Parafoveal2DegPhysioOptics(70);
    %plotFoveal_vs_Parafoveal2DegAdaptiveOptics(70);
end

function plotAchromatic_vs_LMopponent(maxCSF, maxSFvisualized, opticsType, mosaicEccDegs, mosaicSizeDegs, color1, color2)
    mRGCOutputSignalType = 'mRGCs';
     
    stimulusChroma = 'luminance';
    orientationDegs = 120;
    presentationMode = 'drifted';

    theScriptName = 't_spatialCSF';
    useMetaContrast = true;
    useConeContrast = true;
    useFixationalEMs = false;
    whichNoiseFreeNre = 'mRGCMosaic';
    whichNoisyInstanceNre = 'Gaussian';
    whichClassifierEngine = 'rceTemplateDistance';
   [figureFileBaseDir, resultsFileBaseDir] = setupFigureDirectory(theScriptName, ...
        useMetaContrast,useConeContrast,useFixationalEMs,...
        whichNoiseFreeNre,whichNoisyInstanceNre,...
        whichClassifierEngine,mRGCOutputSignalType);

    thresholdsDataFileName = ...
        sprintf('%sCSF_%s_Optics_%s_EccDegs_x%2.1f_%2.1f_SizeDegs_%2.1fx%2.1f_OriDegs_%2.0f_%s.mat', ...
        mRGCOutputSignalType, ...
        stimulusChroma, ...
        opticsType, ...
        mosaicEccDegs(1), mosaicEccDegs(2), ...
        mosaicSizeDegs(1),mosaicSizeDegs(2), ...
        orientationDegs, ...
        presentationMode);
    load(thresholdsDataFileName, 'spatialFreqs', 'thresholdContrasts');
    luminanceThresholds = 1./thresholdContrasts;

    stimulusChroma = 'red-green';
    thresholdsDataFileName = ...
        sprintf('%sCSF_%s_Optics_%s_EccDegs_x%2.1f_%2.1f_SizeDegs_%2.1fx%2.1f_OriDegs_%2.0f_%s.mat', ...
        mRGCOutputSignalType, ...
        stimulusChroma, ...
        opticsType, ...
        mosaicEccDegs(1), mosaicEccDegs(2), ...
        mosaicSizeDegs(1),mosaicSizeDegs(2), ...
        orientationDegs, ...
        presentationMode);
    load(thresholdsDataFileName, 'spatialFreqs', 'thresholdContrasts');
    lmOpponentThresholds = 1./thresholdContrasts;

    lumDir = [1 1 1]';
    redgreenDir = [1 -1 0]';
    normFactor = norm(lumDir)/norm(redgreenDir);
    lmOpponentThresholds = lmOpponentThresholds / normFactor;

    plotComparedDataSets(spatialFreqs, ...
        luminanceThresholds, lmOpponentThresholds, maxCSF, maxSFvisualized, [], ...
        'mRGCs (L+M+S)',  'mRGCs (L-M)', ...
        color1, color2, ...
        fullfile(figureFileBaseDir, sprintf('luma_vs_red_green_eccDegs_%2.1f_%s.pdf', mosaicEccDegs(1), opticsType)));
end


function plotCones_vs_mRGCs(maxCSF, maxSFvisualized, maxSFvisualizedCones, mosaicEccDegs, mosaicSizeDegs, color1, color2)
    mRGCOutputSignalType = 'cones';             % {'cones', 'mRGCs'}
    opticsType = 'loadComputeReadyRGCMosaic';   % {'loadComputeReadyRGCMosaic', 'adaptiveOptics6MM'} 
    stimulusChroma = 'luminance';
    orientationDegs = 120;
    presentationMode = 'drifted';

    theScriptName = 't_spatialCSF';
    useMetaContrast = true;
    useConeContrast = true;
    useFixationalEMs = false;
    whichNoiseFreeNre = 'mRGCMosaic';
    whichNoisyInstanceNre = 'Gaussian';
    whichClassifierEngine = 'rceTemplateDistance';

    [figureFileBaseDir, resultsFileBaseDir] = setupFigureDirectory(theScriptName, ...
        useMetaContrast,useConeContrast,useFixationalEMs,...
        whichNoiseFreeNre,whichNoisyInstanceNre,...
        whichClassifierEngine,mRGCOutputSignalType);

    thresholdsDataFileName = ...
        sprintf('%sCSF_%s_Optics_%s_EccDegs_x%2.1f_%2.1f_SizeDegs_%2.1fx%2.1f_OriDegs_%2.0f_%s.mat', ...
        mRGCOutputSignalType, ...
        stimulusChroma, ...
        opticsType, ...
        mosaicEccDegs(1), mosaicEccDegs(2), ...
        mosaicSizeDegs(1),mosaicSizeDegs(2), ...
        orientationDegs, ...
        presentationMode);

    load(fullfile(resultsFileBaseDir,thresholdsDataFileName), 'spatialFreqs', 'thresholdContrasts');
    coneThresholds = 1./thresholdContrasts;

    mRGCOutputSignalType = 'mRGCs';
    [~, resultsFileBaseDir] = setupFigureDirectory(theScriptName, ...
        useMetaContrast,useConeContrast,useFixationalEMs,...
        whichNoiseFreeNre,whichNoisyInstanceNre,...
        whichClassifierEngine,mRGCOutputSignalType);
    thresholdsDataFileName = ...
        sprintf('%sCSF_%s_Optics_%s_EccDegs_x%2.1f_%2.1f_SizeDegs_%2.1fx%2.1f_OriDegs_%2.0f_%s.mat', ...
        mRGCOutputSignalType, ...
        stimulusChroma, ...
        opticsType, ...
        mosaicEccDegs(1), mosaicEccDegs(2), ...
        mosaicSizeDegs(1),mosaicSizeDegs(2), ...
        orientationDegs, ...
        presentationMode);
    load(fullfile(resultsFileBaseDir,thresholdsDataFileName), 'spatialFreqs', 'thresholdContrasts');
    mRGCThresholds = 1./thresholdContrasts;

    plotComparedDataSets(spatialFreqs, ...
        coneThresholds, mRGCThresholds, maxCSF, maxSFvisualizedCones, maxSFvisualized, ...
        'cones (L+M+S)', 'mRGCs (L+M+S)', ...
        color1, color2, ...
        fullfile(figureFileBaseDir, sprintf('mRGCs_vs_cones_eccDegs_%2.1f.pdf', mosaicEccDegs(1))));
end



function plotPhysiologicalOptics_vs_adaptiveOptics(maxCSF, stimulusChroma, maxSFvisualized, maxSFvisualizedAO, mosaicEccDegs, mosaicSizeDegs, color1, color2)
    mRGCOutputSignalType = 'mRGCs';
    opticsType = 'adaptiveOptics6MM';   % {'loadComputeReadyRGCMosaic', 'adaptiveOptics6MM'} 
   
    orientationDegs = 120;
    presentationMode = 'drifted';

    theScriptName = 't_spatialCSF';
    useMetaContrast = true;
    useConeContrast = true;
    useFixationalEMs = false;
    whichNoiseFreeNre = 'mRGCMosaic';
    whichNoisyInstanceNre = 'Gaussian';
    whichClassifierEngine = 'rceTemplateDistance';

    [figureFileBaseDir, resultsFileBaseDir] = setupFigureDirectory(theScriptName, ...
        useMetaContrast,useConeContrast,useFixationalEMs,...
        whichNoiseFreeNre,whichNoisyInstanceNre,...
        whichClassifierEngine,mRGCOutputSignalType);

    thresholdsDataFileName = ...
        sprintf('%sCSF_%s_Optics_%s_EccDegs_x%2.1f_%2.1f_SizeDegs_%2.1fx%2.1f_OriDegs_%2.0f_%s.mat', ...
        mRGCOutputSignalType, ...
        stimulusChroma, ...
        opticsType, ...
        mosaicEccDegs(1), mosaicEccDegs(2), ...
        mosaicSizeDegs(1),mosaicSizeDegs(2), ...
        orientationDegs, ...
        presentationMode);
    load(thresholdsDataFileName, 'spatialFreqs', 'thresholdContrasts');
    adaptiveOpticsThresholds = 1./thresholdContrasts;
    
    if (strcmp(stimulusChroma, 'red-green'))
        lumDir = [1 1 1]';
        redgreenDir = [1 -1 0]';
        normFactor = norm(lumDir)/norm(redgreenDir);
        adaptiveOpticsThresholds = adaptiveOpticsThresholds / normFactor;
    end


    opticsType = 'loadComputeReadyRGCMosaic'; 
    
    thresholdsDataFileName = ...
        sprintf('%sCSF_%s_Optics_%s_EccDegs_x%2.1f_%2.1f_SizeDegs_%2.1fx%2.1f_OriDegs_%2.0f_%s.mat', ...
        mRGCOutputSignalType, ...
        stimulusChroma, ...
        opticsType, ...
        mosaicEccDegs(1), mosaicEccDegs(2), ...
        mosaicSizeDegs(1),mosaicSizeDegs(2), ...
        orientationDegs, ...
        presentationMode);
    load(thresholdsDataFileName, 'spatialFreqs', 'thresholdContrasts');
    physiologicalOpticsThresholds = 1./thresholdContrasts;

    if (strcmp(stimulusChroma, 'red-green'))
        lumDir = [1 1 1]';
        redgreenDir = [1 -1 0]';
        normFactor = norm(lumDir)/norm(redgreenDir);
        physiologicalOpticsThresholds = physiologicalOpticsThresholds / normFactor;
    end

    plotComparedDataSets(spatialFreqs, ...
        physiologicalOpticsThresholds, adaptiveOpticsThresholds, maxCSF, maxSFvisualized, maxSFvisualizedAO, ...
        'physiological optics', 'adaptive optics', ...
        color1, color2, ...
        fullfile(figureFileBaseDir, sprintf('physio_vs_adaptiveOptics_%s.pdf', stimulusChroma)));
end


function plotFoveal_vs_Parafoveal2DegAdaptiveOptics(maxCSF)
    mRGCOutputSignalType = 'mRGCs';
    opticsType = 'adaptiveOptics6MM';   % {'loadComputeReadyRGCMosaic', 'adaptiveOptics6MM'} 
    stimulusChroma = 'luminance';
    mosaicEccDegs = [0 0];
    mosaicSizeDegs = [1 1];
    orientationDegs = 120;
    presentationMode = 'drifted';

    thresholdsDataFileName = ...
        sprintf('%sCSF_%s_Optics_%s_EccDegs_x%2.1f_%2.1f_SizeDegs_%2.1fx%2.1f_OriDegs_%2.0f_%s.mat', ...
        mRGCOutputSignalType, ...
        stimulusChroma, ...
        opticsType, ...
        mosaicEccDegs(1), mosaicEccDegs(2), ...
        mosaicSizeDegs(1),mosaicSizeDegs(2), ...
        orientationDegs, ...
        presentationMode);


    load(thresholdsDataFileName, 'spatialFreqs', 'thresholdContrasts');
    fovealThresholds = 1./thresholdContrasts;


    mosaicEccDegs = [-2 0];
    thresholdsDataFileName = ...
        sprintf('%sCSF_%s_Optics_%s_EccDegs_x%2.1f_%2.1f_SizeDegs_%2.1fx%2.1f_OriDegs_%2.0f_%s.mat', ...
        mRGCOutputSignalType, ...
        stimulusChroma, ...
        opticsType, ...
        mosaicEccDegs(1), mosaicEccDegs(2), ...
        mosaicSizeDegs(1),mosaicSizeDegs(2), ...
        orientationDegs, ...
        presentationMode);
    load(thresholdsDataFileName, 'spatialFreqs', 'thresholdContrasts');
    parafovealThresholds = 1./thresholdContrasts;

    plotComparedDataSets(spatialFreqs, ...
        fovealThresholds, parafovealThresholds, maxCSF, ...
        'foveal mosaic, AO', 'parafoveal (2 deg) mosaic, AO', 'foveal_vs_parafoveal_AO.pdf');

end


function plotFoveal_vs_Parafoveal2DegPhysioOptics(maxCSF)
    mRGCOutputSignalType = 'mRGCs';
    opticsType = 'loadComputeReadyRGCMosaic';   % {'loadComputeReadyRGCMosaic', 'adaptiveOptics6MM'} 
    stimulusChroma = 'luminance';
    mosaicEccDegs = [0 0];
    mosaicSizeDegs = [1 1];
    orientationDegs = 120;
    presentationMode = 'drifted';

    thresholdsDataFileName = ...
        sprintf('%sCSF_%s_Optics_%s_EccDegs_x%2.1f_%2.1f_SizeDegs_%2.1fx%2.1f_OriDegs_%2.0f_%s.mat', ...
        mRGCOutputSignalType, ...
        stimulusChroma, ...
        opticsType, ...
        mosaicEccDegs(1), mosaicEccDegs(2), ...
        mosaicSizeDegs(1),mosaicSizeDegs(2), ...
        orientationDegs, ...
        presentationMode);
    load(thresholdsDataFileName, 'spatialFreqs', 'thresholdContrasts');
    fovealThresholds = 1./thresholdContrasts;


    mosaicEccDegs = [-2 0];
    thresholdsDataFileName = ...
        sprintf('%sCSF_%s_Optics_%s_EccDegs_x%2.1f_%2.1f_SizeDegs_%2.1fx%2.1f_OriDegs_%2.0f_%s.mat', ...
        mRGCOutputSignalType, ...
        stimulusChroma, ...
        opticsType, ...
        mosaicEccDegs(1), mosaicEccDegs(2), ...
        mosaicSizeDegs(1),mosaicSizeDegs(2), ...
        orientationDegs, ...
        presentationMode);
    load(thresholdsDataFileName, 'spatialFreqs', 'thresholdContrasts');
    parafovealThresholds = 1./thresholdContrasts;

    plotComparedDataSets(spatialFreqs, ...
        fovealThresholds, parafovealThresholds, maxCSF, ...
        'foveal mosaic & optics', 'parafoveal (2 deg) mosaic & optics', 'foveal_vs_parafoveal.pdf');

end


function plotAchromatic_vs_LMopponentAO(maxCSF)
    mRGCOutputSignalType = 'mRGCs';
    opticsType = 'adaptiveOptics6MM';   % {'loadComputeReadyRGCMosaic', 'adaptiveOptics6MM'} 
    stimulusChroma = 'luminance';
    mosaicEccDegs = [0 0];
    mosaicSizeDegs = [1 1];
    orientationDegs = 120;
    presentationMode = 'drifted';

    thresholdsDataFileName = ...
        sprintf('%sCSF_%s_Optics_%s_EccDegs_x%2.1f_%2.1f_SizeDegs_%2.1fx%2.1f_OriDegs_%2.0f_%s.mat', ...
        mRGCOutputSignalType, ...
        stimulusChroma, ...
        opticsType, ...
        mosaicEccDegs(1), mosaicEccDegs(2), ...
        mosaicSizeDegs(1),mosaicSizeDegs(2), ...
        orientationDegs, ...
        presentationMode);
    load(thresholdsDataFileName, 'spatialFreqs', 'thresholdContrasts');
    luminanceThresholds = 1./thresholdContrasts;

    stimulusChroma = 'red-green';
    thresholdsDataFileName = ...
        sprintf('%sCSF_%s_Optics_%s_EccDegs_x%2.1f_%2.1f_SizeDegs_%2.1fx%2.1f_OriDegs_%2.0f_%s.mat', ...
        mRGCOutputSignalType, ...
        stimulusChroma, ...
        opticsType, ...
        mosaicEccDegs(1), mosaicEccDegs(2), ...
        mosaicSizeDegs(1),mosaicSizeDegs(2), ...
        orientationDegs, ...
        presentationMode);
    load(thresholdsDataFileName, 'spatialFreqs', 'thresholdContrasts');
    lmOpponentThresholds = 1./thresholdContrasts;

    lumDir = [1 1 1]';
    redgreenDir = [1 -1 0]';
    normFactor = norm(lumDir)/norm(redgreenDir);
    lmOpponentThresholds = lmOpponentThresholds / normFactor;

    plotComparedDataSets(spatialFreqs, ...
        luminanceThresholds, lmOpponentThresholds, maxCSF, ...
        'physiological optics (L+M+S), AO', 'physiological optics (L-M), AO', 'luma_vs_red_green_AO.pdf');
end






function plotComparedDataSets(spatialFreqs, ...
    thresholds1, thresholds2, maxCSFvisualized, maxSFvisualized, maxSFvisualizedHigh, ...
    legend1, legend2, color1, color2, pdfFileName)

    hFig = figure(1); clf;
    set(hFig, 'Color', [1 1 1]);
    ff = PublicationReadyPlotLib.figureComponents('1x1 standard tall figure');
    theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);

    hold(theAxes{1,1}, 'on');


    % CSF-1
    if (~isempty(maxSFvisualized))
        idx = find(spatialFreqs <= maxSFvisualized);
    else
        idx = 1:numel(spatialFreqs);
    end

    plot(theAxes{1,1}, spatialFreqs(idx), thresholds1(idx), '-', 'LineWidth', 3, 'Color', color1*0.5);
    plot(theAxes{1,1}, spatialFreqs(idx), thresholds1(idx), '-', 'LineWidth', 1.5, 'Color', color1);

    p1 = scatter(theAxes{1,1}, spatialFreqs(idx), thresholds1(idx), 300, ...
        'MarkerFaceColor', color1, 'MarkerEdgeColor', color1*0.5, ...
        'MarkerFaceAlpha', 1.0, 'LineWidth', 1.5);

    idx0 = idx;

    % CSF-2
    if (~isempty(maxSFvisualizedHigh))
        idx = find(spatialFreqs <= maxSFvisualizedHigh);
    else
        idx = 1:numel(spatialFreqs);
    end

    plot(theAxes{1,1}, spatialFreqs(idx), thresholds2(idx), '-', 'LineWidth', 3, 'Color', color2*0.5);
    plot(theAxes{1,1}, spatialFreqs(idx), thresholds2(idx), '-', 'LineWidth', 1.5, 'Color', color2);


    p2 = scatter(theAxes{1,1}, spatialFreqs(idx),  thresholds2(idx), 300, ...
        'MarkerFaceColor', color2, 'MarkerEdgeColor', color2*0.5, ...
        'MarkerFaceAlpha', 1.0, 'LineWidth', 1.5);

    scatter(theAxes{1,1}, spatialFreqs(idx0), thresholds1(idx0), 300, ...
        'MarkerFaceColor', color1, 'MarkerEdgeColor', color1*0.5, ...
        'MarkerFaceAlpha', 0.5, 'LineWidth', 1.5);

    scatter(theAxes{1,1}, spatialFreqs(idx),  thresholds2(idx), 300, ...
        'MarkerFaceColor', color2, 'MarkerEdgeColor', color2*0.5, ...
        'MarkerFaceAlpha', 0.5, 'LineWidth', 1.5);

    legend(theAxes{1,1}, [p1, p2], {legend1, legend2})
    % Axes scaling and ticks
    set(theAxes{1,1}, 'XTick', [0.1 0.3 1 3 10 30 100]);
    set(theAxes{1,1}, 'YTick', 0:10:200);

    % Finalize figure using the Publication-Ready format
    xLims = [0.05 200]; 

    yLims = [1 maxCSFvisualized];
    set(theAxes{1,1}, 'XScale', 'log', 'YScale', 'log', 'YTick', [1 3 10 30 100]);
    set(theAxes{1,1}, 'XLim', xLims, 'YLim', yLims);
    %PublicationReadyPlotLib.offsetAxes(theAxes{1,1},ff, xLims, yLims);
    PublicationReadyPlotLib.labelAxes(theAxes{1,1},ff, 'spatial frequency (c/deg)', 'contrast sensitivity');
    PublicationReadyPlotLib.applyFormat(theAxes{1,1},ff);

    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);

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
