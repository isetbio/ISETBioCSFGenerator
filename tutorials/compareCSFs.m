function compareCSFs()
% UTTBSkip

    % Colors
    backgroundColor = [1 1 1];
    colorCones = [1 0.5 0.7];
    colorMRGCPhysiologicalOpticsLuminance = [0.85 0.85 0.85];
    colorMRGCPhysiologicalOpticsRedGreen = [1 0.0 0.7];
    colorMRGCAdaptiveOpticsLuminance = 0.75*colorMRGCPhysiologicalOpticsLuminance;
    colorMRGCAdaptiveOpticsRedGreen = 0.75*colorMRGCPhysiologicalOpticsRedGreen;

    


    % REVIEWER FIGURES (Gaussian noise effect)

    colorMRGCPhysiologicalOpticsLuminanceDoubleNoise = [0.85 0.85 0.0];
    colorMRGCPhysiologicalOpticsRedGreenDoubleNoise = [1.0 0.0 0.0];
    colorMRGCAdaptiveOpticsLuminanceDoubleNoise = 0.75*colorMRGCPhysiologicalOpticsLuminanceDoubleNoise;
    colorMRGCAdaptiveOpticsRedGreenDoubleNoise = 0.75 * colorMRGCPhysiologicalOpticsRedGreenDoubleNoise;

    % PARAFOVEAL MOSAIC (14 degs)
    maxCSFvisualized = 40;
    psfOrientationDegs = 15;
    mosaicEccDegs = [-14 0];
    mosaicSizeDegs = [4.1 4.1];

    maxSFvisualized = 80;
    maxSFvisualizedCones = 90;
    maxSFvisualizedAO = 120;

    % L+M+S  (physiological optics)
    opticsType = 'loadComputeReadyRGCMosaic'; 
    mRGCOutputSignalType = 'mRGCs';
    finalizePlotWithLegends = '';
    hFig = [];
    ax = [];
    p0 = [];
    p00 = [];
    [hFig, ax, p0, p00] = plotAchromatic_vs_LMopponent(hFig, ax, p0, p00, finalizePlotWithLegends, ...
        maxCSFvisualized, maxSFvisualized, opticsType, mRGCOutputSignalType, ...
        mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, ...
        colorMRGCPhysiologicalOpticsLuminance, colorMRGCPhysiologicalOpticsRedGreen, backgroundColor);

    finalizePlotWithLegends = {' 1 x sigma (L+M+S)',  ' 1 x sigma (L-M)', ' 2 x sigma (L+M+S)',  ' 2 x sigma (L-M)'};
    [hFig, ax, p0, p00] = plotAchromatic_vs_LMopponent(hFig, ax, p0, p00, finalizePlotWithLegends, ...
        maxCSFvisualized, maxSFvisualized, opticsType, mRGCOutputSignalType, ...
        mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs+1, ...
        colorMRGCPhysiologicalOpticsLuminanceDoubleNoise, colorMRGCPhysiologicalOpticsRedGreenDoubleNoise, backgroundColor);


    % L+M+S  (adaptive optics)
    opticsType = 'adaptiveOptics6MM';
    mRGCOutputSignalType = 'mRGCs';
    finalizePlotWithLegends = '';
    hFig = [];
    ax = [];
    p0 = [];
    p00 = [];
    [hFig, ax, p0, p00] = plotAchromatic_vs_LMopponent(hFig, ax, p0, p00, finalizePlotWithLegends, ...
        maxCSFvisualized, maxSFvisualized, opticsType, mRGCOutputSignalType, ...
        mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, ...
        colorMRGCAdaptiveOpticsLuminance, colorMRGCAdaptiveOpticsRedGreen, backgroundColor);

    finalizePlotWithLegends = {' 1 x sigma (L+M+S)',  ' 1 x sigma (L-M)', ' 2 x sigma (L+M+S)',  ' 2 x sigma (L-M)'};
    [hFig, ax, p0, p00] = plotAchromatic_vs_LMopponent(hFig, ax, p0, p00, finalizePlotWithLegends, ...
        maxCSFvisualized, maxSFvisualized, opticsType, mRGCOutputSignalType, ...
        mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs+1, ...
        colorMRGCAdaptiveOpticsLuminanceDoubleNoise, colorMRGCAdaptiveOpticsRedGreenDoubleNoise, backgroundColor);


    pause

    % DEFAULT FIGURES

    % FOVEAL MOSAIC
    psfOrientationDegs = 120;
    mosaicEccDegs = [0 0];
    mosaicSizeDegs = [0.6 0.6];

    maxCSFvisualized = 150;
    maxSFvisualized = 90;
    maxSFvisualizedCones = 90;
    maxSFvisualizedAO = 200;

    % input cone mosaic- vs. mRGC mosaic- CSF
    plotCones_vs_mRGCs(maxCSFvisualized, maxSFvisualized, maxSFvisualizedCones, ...
        mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, ...
        colorCones, colorMRGCPhysiologicalOpticsLuminance, backgroundColor);
    

    % L+M+S vs. L-M (physiological optics)
    opticsType = 'loadComputeReadyRGCMosaic'; 
    mRGCOutputSignalType = 'mRGCs';
    finalizePlotWithLegends = '';
    hFig = [];
    ax = [];
    p0 = [];
    p00 = [];
    [hFig, ax, p0, p00] = plotAchromatic_vs_LMopponent(hFig, ax, p0, p00, finalizePlotWithLegends, ...
        maxCSFvisualized, maxSFvisualized, opticsType, mRGCOutputSignalType, ...
        mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, ...
        colorMRGCPhysiologicalOpticsLuminance, colorMRGCPhysiologicalOpticsRedGreen, backgroundColor);
    
    finalizePlotWithLegends = {' mRGCs (L+M+S)',  ' mRGCs (L-M)', ' cones (L+M+S)', ' cones (L-M)'};
    mRGCOutputSignalType = 'cones';
    plotAchromatic_vs_LMopponent(hFig, ax, p0, p00, finalizePlotWithLegends , ...
        maxCSFvisualized, maxSFvisualized, opticsType, mRGCOutputSignalType, ...
        mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, ...
        colorMRGCPhysiologicalOpticsLuminance, colorMRGCPhysiologicalOpticsRedGreen, backgroundColor);
    



    % L+M+S vs. L-M (adaptive optics)
    opticsType = 'adaptiveOptics6MM';
    mRGCOutputSignalType = 'mRGCs';
    finalizePlotWithLegends = '';
    hFig = [];
    ax = [];
    p0 = [];
    p00 = [];
    [hFig, ax, p0, p00] = plotAchromatic_vs_LMopponent(hFig, ax, p0, p00, finalizePlotWithLegends, ...
        maxCSFvisualized, maxSFvisualizedAO, opticsType, mRGCOutputSignalType, ...
        mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, ...
        colorMRGCAdaptiveOpticsLuminance, colorMRGCAdaptiveOpticsRedGreen, backgroundColor);

    finalizePlotWithLegends = {' mRGCs (L+M+S)',  ' mRGCs (L-M)', ' cones (L+M+S)', ' cones (L-M)'};
    mRGCOutputSignalType = 'cones';
    plotAchromatic_vs_LMopponent(hFig, ax, p0, p00, finalizePlotWithLegends , ...
        maxCSFvisualized, maxSFvisualizedAO, opticsType, mRGCOutputSignalType, ...
        mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, ...
        colorMRGCAdaptiveOpticsLuminance, colorMRGCAdaptiveOpticsRedGreen, backgroundColor);
   
    
    if (1==2)

    % physiological vs. adaptive optics (L+M+S)
    stimulusChroma = 'luminance';
    plotPhysiologicalOptics_vs_adaptiveOptics(maxCSFvisualized, stimulusChroma, maxSFvisualized, maxSFvisualizedAO, mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, colorMRGCPhysiologicalOpticsLuminance, colorMRGCAdaptiveOpticsLuminance, backgroundColor);
    
    % physiological vs. adaptive optics (L-M)
    stimulusChroma = 'red-green';
    plotPhysiologicalOptics_vs_adaptiveOptics(maxCSFvisualized, stimulusChroma, maxSFvisualized, maxSFvisualizedAO, mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, colorMRGCPhysiologicalOpticsRedGreen, colorMRGCAdaptiveOpticsRedGreen, backgroundColor);
    end



    % PARAFOVEAL MOSAIC (4 degs)
    psfOrientationDegs = 45;
    mosaicEccDegs = [-4 0];
    mosaicSizeDegs = [2.1 2.1];
    maxSFvisualized = 80;
    maxSFvisualizedCones = 90;
    maxSFvisualizedAO = 120;

    % input cone mosaic- vs. mRGC mosaic- CSF
    plotCones_vs_mRGCs(maxCSFvisualized, maxSFvisualized, maxSFvisualizedCones, mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, colorCones, colorMRGCPhysiologicalOpticsLuminance, backgroundColor);



    % L+M+S vs. L-M (physiological optics)
    opticsType = 'loadComputeReadyRGCMosaic'; 
    mRGCOutputSignalType = 'mRGCs';
    finalizePlotWithLegends = '';
    hFig = [];
    ax = [];
    p0 = [];
    p00 = [];
    [hFig, ax, p0, p00] = plotAchromatic_vs_LMopponent(hFig, ax, p0, p00, finalizePlotWithLegends, ...
        maxCSFvisualized, maxSFvisualized, opticsType, mRGCOutputSignalType, ...
        mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, ...
        colorMRGCPhysiologicalOpticsLuminance, colorMRGCPhysiologicalOpticsRedGreen, backgroundColor);
    
    finalizePlotWithLegends = {' mRGCs (L+M+S)',  ' mRGCs (L-M)', ' cones (L+M+S)', ' cones (L-M)'};
    mRGCOutputSignalType = 'cones';
    plotAchromatic_vs_LMopponent(hFig, ax, p0, p00, finalizePlotWithLegends , ...
        maxCSFvisualized, maxSFvisualized, opticsType, mRGCOutputSignalType, ...
        mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, ...
        colorMRGCPhysiologicalOpticsLuminance, colorMRGCPhysiologicalOpticsRedGreen, backgroundColor);
    

    % L+M+S vs. L-M (adaptive optics)
    opticsType = 'adaptiveOptics6MM';
    mRGCOutputSignalType = 'mRGCs';
    finalizePlotWithLegends = '';
    hFig = [];
    ax = [];
    p0 = [];
    p00 = [];
    [hFig, ax, p0, p00] = plotAchromatic_vs_LMopponent(hFig, ax, p0, p00, finalizePlotWithLegends, ...
        maxCSFvisualized, maxSFvisualizedAO, opticsType, mRGCOutputSignalType, ...
        mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, ...
        colorMRGCAdaptiveOpticsLuminance, colorMRGCAdaptiveOpticsRedGreen, backgroundColor);

    finalizePlotWithLegends = {' mRGCs (L+M+S)',  ' mRGCs (L-M)', ' cones (L+M+S)', ' cones (L-M)'};
    mRGCOutputSignalType = 'cones';
    plotAchromatic_vs_LMopponent(hFig, ax, p0, p00, finalizePlotWithLegends , ...
        maxCSFvisualized, maxSFvisualizedAO, opticsType, mRGCOutputSignalType, ...
        mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, ...
        colorMRGCAdaptiveOpticsLuminance, colorMRGCAdaptiveOpticsRedGreen, backgroundColor);
   

    

    % physiological vs. adaptive optics (L+M+S)
    stimulusChroma = 'luminance';
    plotPhysiologicalOptics_vs_adaptiveOptics(maxCSFvisualized, stimulusChroma, maxSFvisualized, maxSFvisualizedAO, mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, colorMRGCPhysiologicalOpticsLuminance, colorMRGCAdaptiveOpticsLuminance, backgroundColor);
    
    % physiological vs. adaptive optics (L-M)
    stimulusChroma = 'red-green';
    plotPhysiologicalOptics_vs_adaptiveOptics(maxCSFvisualized, stimulusChroma, maxSFvisualized, maxSFvisualizedAO, mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, colorMRGCPhysiologicalOpticsRedGreen, colorMRGCAdaptiveOpticsRedGreen, backgroundColor);
  

    % PARAFOVEAL MOSAIC (14 degs)
    maxCSFvisualized = 150;
    psfOrientationDegs = 15;
    mosaicEccDegs = [-14 0];
    mosaicSizeDegs = [4.1 4.1];

    % input cone mosaic- vs. mRGC mosaic- CSF
    plotCones_vs_mRGCs(maxCSFvisualized, maxSFvisualized, maxSFvisualizedCones, ...
        mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, ...
        colorCones, colorMRGCPhysiologicalOpticsLuminance, backgroundColor);


    % L+M+S vs. L-M (physiological optics)
    opticsType = 'loadComputeReadyRGCMosaic'; 
    mRGCOutputSignalType = 'mRGCs';
    finalizePlotWithLegends = '';
    hFig = [];
    ax = [];
    p0 = [];
    p00 = [];
    [hFig, ax, p0, p00] = plotAchromatic_vs_LMopponent(hFig, ax, p0, p00, finalizePlotWithLegends, ...
        maxCSFvisualized, maxSFvisualized, opticsType, mRGCOutputSignalType, ...
        mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, ...
        colorMRGCPhysiologicalOpticsLuminance, colorMRGCPhysiologicalOpticsRedGreen, backgroundColor);
    
    finalizePlotWithLegends = {' mRGCs (L+M+S)',  ' mRGCs (L-M)', ' cones (L+M+S)', ' cones (L-M)'};
    mRGCOutputSignalType = 'cones';
    plotAchromatic_vs_LMopponent(hFig, ax, p0, p00, finalizePlotWithLegends , ...
        maxCSFvisualized, maxSFvisualized, opticsType, mRGCOutputSignalType, ...
        mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, ...
        colorMRGCPhysiologicalOpticsLuminance, colorMRGCPhysiologicalOpticsRedGreen, backgroundColor);
  


    % L+M+S vs. L-M (adaptive optics)
    opticsType = 'adaptiveOptics6MM';
    mRGCOutputSignalType = 'mRGCs';
    finalizePlotWithLegends = '';
    hFig = [];
    ax = [];
    p0 = [];
    p00 = [];
    [hFig, ax, p0, p00] = plotAchromatic_vs_LMopponent(hFig, ax, p0, p00, finalizePlotWithLegends, ...
        maxCSFvisualized, maxSFvisualizedAO, opticsType, mRGCOutputSignalType, ...
        mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, ...
        colorMRGCAdaptiveOpticsLuminance, colorMRGCAdaptiveOpticsRedGreen, backgroundColor);

    finalizePlotWithLegends = {' mRGCs (L+M+S)',  ' mRGCs (L-M)', ' cones (L+M+S)', ' cones (L-M)'};
    mRGCOutputSignalType = 'cones';
    plotAchromatic_vs_LMopponent(hFig, ax, p0, p00, finalizePlotWithLegends , ...
        maxCSFvisualized, maxSFvisualizedAO, opticsType, mRGCOutputSignalType, ...
        mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, ...
        colorMRGCAdaptiveOpticsLuminance, colorMRGCAdaptiveOpticsRedGreen, backgroundColor);
   

    % PARAFOVEAL MOSAIC (25 degs)
    maxCSFvisualized = 150;
    psfOrientationDegs = 0;
    mosaicEccDegs = [-25 0];
    mosaicSizeDegs = [7 7];

    % input cone mosaic- vs. mRGC mosaic- CSF
    plotCones_vs_mRGCs(maxCSFvisualized, maxSFvisualized, maxSFvisualizedCones, mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, colorCones, colorMRGCPhysiologicalOpticsLuminance, backgroundColor);

    % L+M+S vs. L-M (physiological optics)
    opticsType = 'loadComputeReadyRGCMosaic'; 
    mRGCOutputSignalType = 'mRGCs';
    finalizePlotWithLegends = '';
    hFig = [];
    ax = [];
    p0 = [];
    p00 = [];
    [hFig, ax, p0, p00] = plotAchromatic_vs_LMopponent(hFig, ax, p0, p00, finalizePlotWithLegends, ...
        maxCSFvisualized, maxSFvisualized, opticsType, mRGCOutputSignalType, ...
        mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, ...
        colorMRGCPhysiologicalOpticsLuminance, colorMRGCPhysiologicalOpticsRedGreen, backgroundColor);
    
    finalizePlotWithLegends = {' mRGCs (L+M+S)',  ' mRGCs (L-M)', ' cones (L+M+S)', ' cones (L-M)'};
    mRGCOutputSignalType = 'cones';
    plotAchromatic_vs_LMopponent(hFig, ax, p0, p00, finalizePlotWithLegends , ...
        maxCSFvisualized, maxSFvisualized, opticsType, mRGCOutputSignalType, ...
        mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, ...
        colorMRGCPhysiologicalOpticsLuminance, colorMRGCPhysiologicalOpticsRedGreen, backgroundColor);
  


    % L+M+S vs. L-M (adaptive optics)
    opticsType = 'adaptiveOptics6MM';
    mRGCOutputSignalType = 'mRGCs';
    finalizePlotWithLegends = '';
    hFig = [];
    ax = [];
    p0 = [];
    p00 = [];
    [hFig, ax, p0, p00] = plotAchromatic_vs_LMopponent(hFig, ax, p0, p00, finalizePlotWithLegends, ...
        maxCSFvisualized, maxSFvisualizedAO, opticsType, mRGCOutputSignalType, ...
        mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, ...
        colorMRGCAdaptiveOpticsLuminance, colorMRGCAdaptiveOpticsRedGreen, backgroundColor);

    finalizePlotWithLegends = {' mRGCs (L+M+S)',  ' mRGCs (L-M)', ' cones (L+M+S)', ' cones (L-M)'};
    mRGCOutputSignalType = 'cones';
    plotAchromatic_vs_LMopponent(hFig, ax, p0, p00, finalizePlotWithLegends , ...
        maxCSFvisualized, maxSFvisualizedAO, opticsType, mRGCOutputSignalType, ...
        mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, ...
        colorMRGCAdaptiveOpticsLuminance, colorMRGCAdaptiveOpticsRedGreen, backgroundColor);
   

    pause
    % PARAFOVEAL MOSAIC (7 degs)
    maxCSFvisualized = 150;
    psfOrientationDegs = 30;
    mosaicEccDegs = [-7 0];
    mosaicSizeDegs = [3.2 3.2];

    % input cone mosaic- vs. mRGC mosaic- CSF
    plotCones_vs_mRGCs(maxCSFvisualized, maxSFvisualized, maxSFvisualizedCones, mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, colorCones, colorMRGCPhysiologicalOpticsLuminance, backgroundColor);

    % L+M+S vs. L-M (physiological optics)
    opticsType = 'loadComputeReadyRGCMosaic';
    mRGCOutputSignalType = 'mRGCs';
    plotAchromatic_vs_LMopponent([], [], [], [], {}, maxCSFvisualized, maxSFvisualized, opticsType, mRGCOutputSignalType, mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, colorMRGCPhysiologicalOpticsLuminance, colorMRGCPhysiologicalOpticsRedGreen, backgroundColor);
    
    
    % physiological vs. adaptive optics (L+M+S)
    stimulusChroma = 'luminance';
    plotPhysiologicalOptics_vs_adaptiveOptics(maxCSFvisualized, stimulusChroma, maxSFvisualized, maxSFvisualizedAO, mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, colorMRGCPhysiologicalOpticsLuminance, colorMRGCAdaptiveOpticsLuminance, backgroundColor);
    
    % physiological vs. adaptive optics (L-M)
    stimulusChroma = 'red-green';
    plotPhysiologicalOptics_vs_adaptiveOptics(maxCSFvisualized, stimulusChroma, maxSFvisualized, maxSFvisualizedAO, mosaicEccDegs, mosaicSizeDegs, psfOrientationDegs, colorMRGCPhysiologicalOpticsRedGreen, colorMRGCAdaptiveOpticsRedGreen, backgroundColor);

end

function [hFig, ax, p0, p00] = plotAchromatic_vs_LMopponent(...
    hFig, ax, p0, p00, finalizeFigureWithLegends, ...
    maxCSF, maxSFvisualized, opticsType, ...
    mRGCOutputSignalType, mosaicEccDegs, mosaicSizeDegs, orientationDegs, ...
    color1, color2, backgroundColor)
    
    stimulusChroma = 'luminance';
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
    luminanceSpatialFreqs = spatialFreqs;


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
    lmOpponentSpatialFreqs = spatialFreqs;

    if (isempty(finalizeFigureWithLegends))
        finalizeFigure = false;
        theLegends = {};
    else
        finalizeFigure = true;
        theLegends = finalizeFigureWithLegends;
    end

    switch (mRGCOutputSignalType)
        case 'cones'
            theSymbols = {'^', '^'};
        case 'mRGCs'
            theSymbols = {'o', 'o'};
        otherwise
            error('Unknown mRGCOutputSignalType: ''%s''.', mRGCOutputSignalType);
    end

    [hFig, ax, p0, p00] = plotComparedDataSets(...
        hFig, ax, p0, p00, finalizeFigure, ...
        luminanceSpatialFreqs, lmOpponentSpatialFreqs, ...
        luminanceThresholds, lmOpponentThresholds, ...
        maxCSF, maxSFvisualized, [], ...
        theSymbols, theLegends, ...
        color1, color2, backgroundColor, ...
        fullfile(figureFileBaseDir, sprintf('%s_luma_vs_red_green_eccDegs_%2.1f_%s.pdf', mRGCOutputSignalType, mosaicEccDegs(1), opticsType)));
end


function plotCones_vs_mRGCs(maxCSF, maxSFvisualized, maxSFvisualizedCones, mosaicEccDegs, mosaicSizeDegs, orientationDegs, color1, color2, backgroundColor)
    mRGCOutputSignalType = 'cones';             % {'cones', 'mRGCs'}
    opticsType = 'loadComputeReadyRGCMosaic';   % {'loadComputeReadyRGCMosaic', 'adaptiveOptics6MM'} 
    stimulusChroma = 'luminance';
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
    coneSpatialFreqs = spatialFreqs;

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
    mRGCspatialFreqs = spatialFreqs;

    plotComparedDataSets(...
        [], [], [], [], true, ...
        coneSpatialFreqs, mRGCspatialFreqs, ...
        coneThresholds, mRGCThresholds, maxCSF, maxSFvisualizedCones, maxSFvisualized, ...
        {'^', 'o'}, {'cones (L+M+S)', 'mRGCs (L+M+S)'}, ...
        color1, color2, backgroundColor, ...
        fullfile(figureFileBaseDir, sprintf('mRGCs_vs_cones_eccDegs_%2.1f.pdf', mosaicEccDegs(1))));
end



function plotPhysiologicalOptics_vs_adaptiveOptics(maxCSF, stimulusChroma, maxSFvisualized, maxSFvisualizedAO, mosaicEccDegs, mosaicSizeDegs, orientationDegs, color1, color2, backgroundColor)
    mRGCOutputSignalType = 'mRGCs';
    opticsType = 'adaptiveOptics6MM';   % {'loadComputeReadyRGCMosaic', 'adaptiveOptics6MM'} 
   
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
    adaptiveOpticsSpatialFreqs = spatialFreqs;

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
    physiologicalOpticsSpatialFreqs = spatialFreqs;

    if (strcmp(stimulusChroma, 'red-green'))
        lumDir = [1 1 1]';
        redgreenDir = [1 -1 0]';
        normFactor = norm(lumDir)/norm(redgreenDir);
        physiologicalOpticsThresholds = physiologicalOpticsThresholds / normFactor;
    end

    plotComparedDataSets(...
        [], [], [], [], true, ...
        physiologicalOpticsSpatialFreqs, adaptiveOpticsSpatialFreqs, ...
        physiologicalOpticsThresholds, adaptiveOpticsThresholds, maxCSF, maxSFvisualized, maxSFvisualizedAO, ...
        {'o', 'o'}, {'physiological optics', 'adaptive optics'}, ...
        color1, color2, backgroundColor, ...
        fullfile(figureFileBaseDir, sprintf('physio_vs_adaptiveOptics_%s_eccDegs_%2.1f.pdf', stimulusChroma, mosaicEccDegs(1))));
end




function [hFig, ax, p0, p00] = plotComparedDataSets(hFig, ax, p0, p00, finalizeFigure, ...
    spatialFreqs1, spatialFreqs2, ...
    thresholds1, thresholds2, maxCSFvisualized, maxSFvisualized, maxSFvisualizedHigh, ...
    theSymbols, theLegends, color1, color2, backgroundColor, pdfFileName)


    ff = PublicationReadyPlotLib.figureComponents('1x1 standard tall figure');
    if (isempty(hFig)) && (isempty(ax))
        hFig = figure(1); clf;
        set(hFig, 'Color', [1 1 1]);
        theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);
        ax = theAxes{1,1};
    end


    hold(ax, 'on');


    % CSF-1
    if (~isempty(maxSFvisualized))
        idx = find(spatialFreqs1 <= maxSFvisualized);
    else
        idx = 1:numel(spatialFreqs1);
    end

    plot(ax, spatialFreqs1(idx), thresholds1(idx), '-', 'LineWidth', 3, 'Color', color1*0.5);
    plot(ax, spatialFreqs1(idx), thresholds1(idx), '-', 'LineWidth', 1.5, 'Color', color1);

    p1 = plot(ax, spatialFreqs1(idx), thresholds1(idx), 'o', 'MarkerSize', sqrt(300), ...
        'Marker', theSymbols{1}, 'MarkerFaceColor', color1, 'MarkerEdgeColor', color1*0.5, ...
        'LineWidth', 1.5);

    idx0 = idx;

    % CSF-2
    if (~isempty(maxSFvisualizedHigh))
        idx = find(spatialFreqs2 <= maxSFvisualizedHigh);
    else
        idx = 1:numel(spatialFreqs2);
    end

    plot(ax, spatialFreqs2(idx), thresholds2(idx), '-', 'LineWidth', 3, 'Color', color2*0.5);
    plot(ax, spatialFreqs2(idx), thresholds2(idx), '-', 'LineWidth', 1.5, 'Color', color2);


    p2 = plot(ax, spatialFreqs2(idx),  thresholds2(idx), 'o', 'MarkerSize', sqrt(300), ...
        'Marker', theSymbols{2}, 'MarkerFaceColor', color2, 'MarkerEdgeColor', color2*0.5, ...
        'LineWidth', 1.5);

    scatter(ax, spatialFreqs1(idx0), thresholds1(idx0), 300, ...
        'Marker', theSymbols{1}, 'MarkerFaceColor', color1, 'MarkerEdgeColor', color1*0.5, ...
        'MarkerFaceAlpha', 0.5, 'LineWidth', 1.5);

    scatter(ax, spatialFreqs2(idx),  thresholds2(idx), 300, ...
        'Marker', theSymbols{2}, 'MarkerFaceColor', color2, 'MarkerEdgeColor', color2*0.5, ...
        'MarkerFaceAlpha', 0.5, 'LineWidth', 1.5);



    if (finalizeFigure)
        if (~isempty(p0) && (~isempty(p00))) && (~isempty(theLegends))
            legend(ax, [p0, p00, p1, p2], theLegends, 'Location', 'SouthWest')
        else
            legend(ax, [p1, p2], theLegends, 'Location', 'SouthWest')
        end
        % Axes scaling and ticks
        set(ax, 'XTick', [0.1 0.3 1 3 10 30 100]);
        set(ax, 'YTick', 0:10:200);
    
        ff.backgroundColor = backgroundColor;
        % Finalize figure using the Publication-Ready format
        xLims = [0.05 200]; 
    
        yLims = [1 maxCSFvisualized];
        set(ax, 'XScale', 'log', 'YScale', 'log', 'YTick', [1 3 10 30 100]);
        set(ax, 'XLim', xLims, 'YLim', yLims);
        
    
        %PublicationReadyPlotLib.offsetAxes(theAxes{1,1},ff, xLims, yLims);
        PublicationReadyPlotLib.labelAxes(ax,ff, 'spatial frequency (c/deg)', 'contrast sensitivity');
        PublicationReadyPlotLib.applyFormat(ax,ff);
    
        set(hFig, 'Color', backgroundColor);
        NicePlot.exportFigToPDF(pdfFileName, hFig, 300);
    end

    p0 = p1;
    p00 = p2;

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
