function plotMultipleIsoresponseLMplaneEllipseThresholds(Pts_coneContrast)
%
% Script to integrate LMellipse data points computes for a set of reference contrasts
% These data are previously computed by running t_isoresponseLMplaneEllipses()
% The integrated data are then loaded by compareHumanToMRGCmosaicEllipses()
% which generates Figure 13 of the 2025 R01


% History:
%    11/01/25  NPC  Wrote it.

% Examples:
%{

% UTTBSkip

    plotMultipleIsoresponseLMplaneEllipseThresholds(Pts_coneContrast);

%}

    % Choose which condition to load data for
    useMetaContrast =~true;
    useConeContrast = true;
    mRGCOutputSignalType = 'mRGCs'; % Choose from {'cones' %'mRGCs'};
    noiseType =  'Gaussian_rceTemplateDistance'; % Choose from {'Gaussian_rceTemplateDistance', 'Poisson_rcePoisson'}


    matlabDir = strrep(isetbioRootPath, '/toolboxes/isetbio', '');
    rootDir = fullfile(matlabDir, 'toolboxes/ISETBioCSFGenerator/local/t_isoresponseLMplaneEllipses/results');


    if (strcmp(mRGCOutputSignalType, 'mRGCs'))
        nonLinearityString = 'ON_OFF_Emulation_nLinearRect_''half''_c50_0.07_n1.2_s1.0';
        nonLinearityString = 'ON_OFF_Emulation_nLinearRect_''half''_c50_0.05_n2.0_s1.0';
        %nonLinearityString = 'LINEAR';
    else
        nonLinearityString = '';
    end

    if (~isempty(nonLinearityString))
        resultsFileBaseDir = fullfile(rootDir, ...
            sprintf('t_isoresponseLMplaneEllipses_Meta_%d_ConeContrast_%d_FEMs_0_mRGCMosaic_%s_%s_%s', ...
                useMetaContrast, useConeContrast, noiseType, mRGCOutputSignalType, nonLinearityString));
    else
        resultsFileBaseDir = fullfile(rootDir, ...
            sprintf('t_isoresponseLMplaneEllipses_Meta_%d_ConeContrast_%d_FEMs_0_mRGCMosaic_%s_%s', ...
                useMetaContrast, useConeContrast, noiseType, mRGCOutputSignalType));
    end


    figureFileBaseDir = strrep(resultsFileBaseDir, 'results', 'figures');

    examinedSpatialFrequencyCPD = 0;
    opticsType = 'loadComputeReadyRGCMosaic';
    mosaicEccDegs = [0 0.0];
    mosaicSizeDegs = [1 1];
    orientationDegs = 90;
    
    presentationMode = 'static';

    initializeFig = true;

    theThresholdAxes = []; 
    maxVisualizedThreshold = 0.35;

    allData = containers.Map();

    for iC = 1:size(Pts_coneContrast,1)

        referenceLMSconeContrast = Pts_coneContrast(iC,:);

        thresholdsDataFileName = ...
                sprintf('%sCSF_SF_%2.2fCPD_Optics_%s_EccDegs_x%2.1f_%2.1f_SizeDegs_%2.1fx%2.1f_OriDegs_%2.0f_RefLMScontrast_%2.1f_%2.1f_%2.1f%s.mat', ...
                mRGCOutputSignalType, ...
                examinedSpatialFrequencyCPD, ...
                opticsType, ...
                mosaicEccDegs(1), mosaicEccDegs(2), ...
                mosaicSizeDegs(1),mosaicSizeDegs(2), ...
                orientationDegs, ...
                referenceLMSconeContrast(1)*100, ...
                referenceLMSconeContrast(2)*100, ...
                referenceLMSconeContrast(3)*100, ...
                presentationMode);


        load(fullfile(resultsFileBaseDir,thresholdsDataFileName), ...
            'options', ...
            'theLMSconeContrastDirections', 'theDeltaLMSconeContrastDirections', ...
            'stimulusRMSLMconeContrast', 'examinedDirectionsOnLMplane', 'thresholdContrasts');


        gratingSceneParams = struct( ...
            'meanLuminanceCdPerM2', options.meanLuminanceCdPerM2, ...
            'meanChromaticityXY', options.meanChromaticityXY, ...
            'backgroundLMSconeExcitations', options.backgroundLMSconeExcitations, ...
            'spectralSupport', 400:20:750, ...
            'fovDegs', options.stimSizeDegs, ...
            'pixelsNum', options.pixelsNum, ...
            'spatialEnvelope', 'rect', ...
            'spatialEnvelopeRadiusDegs', options.stimSizeDegs, ...
            'orientation', options.orientationDegs, ...
            'presentationMode', 'flashed', ... 
            'duration', options.frameDurationSeconds, ...
            'frameDurationSeconds', options.frameDurationSeconds, ...
            'spatialPhase', options.spatialPhaseDegs);

        % Figure for plotting the thresholds on the LM plane
        % along with the stimuli
        hFigStimuliAndThresholds = figure(2346); clf;
        set(hFigStimuliAndThresholds, 'Color', [1 1 1], 'Position', [10 10 1200 1200]);
        set(hFigStimuliAndThresholds, 'HandleVisibility', 'off');

        exportFig = true;
        figName = sprintf('refC_%2.1f_%2.1f_%2.1f%s', 100*referenceLMSconeContrast(1), 100*referenceLMSconeContrast(2), 100*referenceLMSconeContrast(3));
        skippedDirections = 0;
        [thresholdDeltaConeContrasts, theThresholdAxes, theFittedEllipsePoints] = visualizeIsoThresholdEllipsesOnLMplane(...
            theLMSconeContrastDirections, ...
            theDeltaLMSconeContrastDirections, ...
            stimulusRMSLMconeContrast, ...
            examinedSpatialFrequencyCPD, ...
            gratingSceneParams, ...
            examinedDirectionsOnLMplane, ...
            skippedDirections, ...
            thresholdContrasts, ...
            maxVisualizedThreshold, ...
            referenceLMSconeContrast, ...
            figureFileBaseDir, hFigStimuliAndThresholds, ...
            exportFig, initializeFig, theThresholdAxes, figName);

        initializeFig = false;


        theThresholdsDataFile = sprintf('%s_%s_%s.mat', mRGCOutputSignalType, noiseType, figName);

        allData(figName) = struct(...
            'referenceLMconeContrast', referenceLMSconeContrast, ...
            'thresholdDeltaConeContrasts', thresholdDeltaConeContrasts, ...
            'fittedEllipse', theFittedEllipsePoints);

        save(fullfile(resultsFileBaseDir,theThresholdsDataFile), 'referenceLMSconeContrast', 'thresholdDeltaConeContrasts', 'theFittedEllipsePoints');
        fprintf('LM cone contrast thresholds saved in %s\n', theThresholdsDataFile);
    end


    figName = sprintf('Summary%s',mRGCOutputSignalType);
    set(hFigStimuliAndThresholds, 'HandleVisibility', 'on');

    if ((~isempty(figureFileBaseDir)) && (exportFig))
        theFigName = fullfile(figureFileBaseDir,figName);
        NicePlot.exportFigToPNG(theFigName, hFigStimuliAndThresholds, 300);
    end

    theSummaryDataFile = fullfile(resultsFileBaseDir, sprintf('%s_%s_Summary.mat', mRGCOutputSignalType, noiseType));
    save(theSummaryDataFile, 'allData');

    fprintf(sprintf('Summary threshold data for all reference points saved in %s\n', theSummaryDataFile));

end
