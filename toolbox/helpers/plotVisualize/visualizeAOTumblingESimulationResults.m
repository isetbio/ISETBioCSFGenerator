function visualizeAOTumblingESimulationResults(questObj, threshold, fittedPsychometricParams, ...
    thresholdParameters, whichSceneFrame, tumblingEsceneEngines, backgroundSceneEngine, theNeuralEngine, pdfFileName)

    % Choose which sizes to display
    fittedPsychometricFunction = questObj.qpPF(questObj.estDomain', fittedPsychometricParams);
    examinedParameterAxis = 10.^(questObj.estDomain)*thresholdParameters.maxParamValue;

    hFig = figure(3); clf;
    
    % Flag indicating whether to visualize the noise-free cone mosaic
    % excitations or noisy instances (generates a video)
    visualizeNoiseFreeMosaicActivation = true;

    % Define plot range in degrees
    domainSize = 0.2;

    if (visualizeNoiseFreeMosaicActivation == false)
        % If we generate a video of noisy response instances, do so for a high-performance value
        % and also increase the mosaic integration time to 2 seconds
        performanceValuesExamined = 0.95;
        theNeuralEngine.neuralPipeline.noiseFreeResponse.coneMosaic.integrationTime = 2;

        % Figure setup
        set(hFig, 'Position', [10 10 1500 350], 'Color', [1 1 1]);

        % Video setup
        videoOBJ = VideoWriter('NoisyModulations', 'MPEG-4');
        videoOBJ.FrameRate = 10;
        videoOBJ.Quality = 100;
        videoOBJ.open();
    else
        % Performance levels to examine
        performanceValuesExamined = [0.26 0.5 0.8];
        % Figure setup
        set(hFig, 'Position', [10 10 1500 1000], 'Color', [1 1 1]);
    end

    activationColorMap = brewermap(1024, '*RdBu');

    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', numel(performanceValuesExamined), ...
       'colsNum', numel(tumblingEsceneEngines), ...
       'heightMargin',  0.06, ...
       'widthMargin',    0.05, ...
       'leftMargin',     0.07, ...
       'rightMargin',    0.03, ...
       'bottomMargin',   0.05, ...
       'topMargin',      0.05);


    for letterSizeIndex = 1:numel(performanceValuesExamined)
        [~,idx] = min(abs(squeeze(fittedPsychometricFunction(:, 2))-performanceValuesExamined(letterSizeIndex)));
        letterSizeDegs = examinedParameterAxis(idx);
        for letterRotationIndex = 1:numel(tumblingEsceneEngines)
            % Retrieve the sceneEngine
            sceneEngine = tumblingEsceneEngines{letterRotationIndex};

            % Compute the E scene
            theSceneSequence = sceneEngine.compute(letterSizeDegs);
            theTestScene = theSceneSequence{whichSceneFrame};

            if (letterRotationIndex == 1)
                % Compute the background scene
                theSceneSequence = backgroundSceneEngine.compute(letterSizeDegs);
                theBackgroundScene = theSceneSequence{whichSceneFrame};
            end

            % Compute the optical image of the test scene
            theOI = oiCompute(theNeuralEngine.neuralPipeline.noiseFreeResponse.optics,theTestScene);
            
            if (letterRotationIndex == 1)
                % Compute the optical image of the background scene
                theBackgroundOI = oiCompute(theNeuralEngine.neuralPipeline.noiseFreeResponse.optics,theBackgroundScene);
            end

            % Compute cone mosaic activations to the test scene
            [theNoiseFreeConeMosaicActivation, noisyResponseInstances] = ...
                theNeuralEngine.neuralPipeline.noiseFreeResponse.coneMosaic.compute(theOI, 'nTrials', 4);

            if (letterRotationIndex == 1)
                % Compute cone mosaic activations to the background scene
                theNoiseFreeBackgroundConeMosaicActivation = ...
                    theNeuralEngine.neuralPipeline.noiseFreeResponse.coneMosaic.compute(theBackgroundOI);
            end

            % Compute noise-free cone modulations
            theNoiseFreeConeMosaicModulation = 100*VisualizeEExcitationsToContrast(...
                theNoiseFreeConeMosaicActivation, theNoiseFreeBackgroundConeMosaicActivation);

            % Compute noisy cone modulations
            theNoisyConeMosaicModulations = 100*VisualizeEExcitationsToContrast(...
                noisyResponseInstances, theNoiseFreeBackgroundConeMosaicActivation);

            domainVisualizationTicks = struct('x', -domainSize:0.1:domainSize, 'y', -domainSize:0.1:domainSize);
            domainVisualizationTicksForThisPlot = domainVisualizationTicks ;
            if (letterSizeIndex<numel(performanceValuesExamined))
                domainVisualizationTicksForThisPlot.x = [];
            end
            if (letterRotationIndex>1)
                domainVisualizationTicksForThisPlot.y = [];
            end

            d = sqrt(sum(theNeuralEngine.neuralPipeline.noiseFreeResponse.coneMosaic.coneRFpositionsDegs.^2,2));
            roiConeIndices = find(d<domainSize);
            maxModulation = max(abs(theNoiseFreeConeMosaicModulation(roiConeIndices)));

            ax = subplot('Position', subplotPosVectors(letterSizeIndex, letterRotationIndex).v);
            if (visualizeNoiseFreeMosaicActivation)
                theNeuralEngine.neuralPipeline.noiseFreeResponse.coneMosaic.visualize(...
                    'figureHandle', hFig', 'axesHandle', ax, ...
                    'activation', theNoiseFreeConeMosaicModulation, ...
                    'activationRange', maxModulation*[-1 1], ...
                    'verticalActivationColorBar', true, ...
                    'activationColorMap', activationColorMap, ...
                    'colorbarTickLabelColor', [0.3 0.3 0.3],...
                    'domain', 'degrees', ...
                    'domainVisualizationLimits', [-domainSize domainSize -domainSize domainSize], ...
                    'domainVisualizationTicks', domainVisualizationTicksForThisPlot, ...
                    'crossHairsOnMosaicCenter', true, ...
                    'crossHairsColor',[1 0.2 0.2], ...
                    'noXLabel', (letterSizeIndex<numel(performanceValuesExamined)), ...
                    'noYLabel', (letterRotationIndex>1), ...
                    'plotTitle', sprintf('performance level: %2.2f\n(letter size: %2.3f degs)', ...
                            performanceValuesExamined(letterSizeIndex), letterSizeDegs));
            else
                for iTrial = 1:size(noisyResponseInstances,1)
                
                    theNeuralEngine.neuralPipeline.noiseFreeResponse.coneMosaic.visualize(...
                        'figureHandle', hFig', 'axesHandle', ax, ...
                        'activation', theNoisyConeMosaicModulations(iTrial,:,:), ...
                        'activationRange', 2*max(abs(theNoiseFreeConeMosaicModulation(:)))*[-1 1], ...
                        'verticalActivationColorBar', true, ...
                        'activationColorMap', activationColorMap, ...
                        'colorbarTickLabelColor', [0 0 1],...
                        'domain', 'degrees', ...
                        'domainVisualizationLimits', [-domainSize domainSize -domainSize domainSize], ...
                        'domainVisualizationTicks', domainVisualizationTicksForThisPlot, ...
                        'crossHairsOnMosaicCenter', true, ...
                        'crossHairsColor', [1 0.2 0.2], ...
                        'noXLabel', (letterSizeIndex<numel(performanceValuesExamined)), ...
                        'noYLabel', (letterRotationIndex>1), ...
                        'plotTitle', sprintf('performance level: %2.2f\n(letter size: %2.3f degs)', ...
                                performanceValuesExamined(letterSizeIndex), letterSizeDegs));
                     drawnow;
                     videoOBJ.writeVideo(getframe(hFig));
                end

            end

%             visualizeScene(theScenes{letterSizeIndex,letterRotationIndex}, ...
%                 'spatialSupportInDegs', true, ...
%                 'crossHairsAtOrigin', true, ...
%                 'displayRadianceMaps', false, ...
%                 'axesHandle', ax, ...
%                 'noTitle', true);
        end
    end

    if (visualizeNoiseFreeMosaicActivation == false)
        videoOBJ.close();
    end

    if (~isempty(pdfFileName))
        NicePlot.exportFigToPDF(pdfFileName,hFig, 300);
    end
end

% Helper.  
function c = VisualizeEExcitationsToContrast(e,b)
    c = bsxfun(@minus, e, b);
    c = bsxfun(@times, c, 1./b);
end