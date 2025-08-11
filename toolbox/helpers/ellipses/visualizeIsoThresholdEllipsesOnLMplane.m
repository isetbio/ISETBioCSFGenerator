function theThresholdAxes = visualizeIsoThresholdEllipsesOnLMplane(...
    rmsLMconeContrast, examinedSpatialFrequencyCPD, gratingSceneParams, ...
    theChromaticDirections, skippedDirections, theThresholds, ...
    maxThreshold, theVisualizedConeContrastOffset, ...
    figureFileBaseDir, hFig, ...
    exportFig, initializeFig, theThresholdAxes)

    
    %gratingSceneParams.warningInsteadOfErrorOnOutOfGamut = true;

    visualizedLMangles = 0:15:345;
    theLMSconeContrastDirections = computeLMSconeContrastDirections(rmsLMconeContrast, visualizedLMangles);

    % Create the background scene engine
    theNullSceneEngine = createGratingSceneEngine(...
            [0 0 0.4], examinedSpatialFrequencyCPD, ...
            gratingSceneParams);

    % Compute the scene sequence
    theNullSceneSequence = theNullSceneEngine.compute(0.0);

    % Draw
    set(hFig, 'HandleVisibility', 'on');
    figure(hFig);
    
    plotColorCircle = false;

    if (initializeFig)
        clf;
        
        if (plotColorCircle)
            % Fill the figure with the background scene
            theBackgroundAxes = axes('Position', [0 0 1 1]);
            % Fill the figure with the background scene
            theNullSceneEngine.visualizeStaticFrame(...
                    theNullSceneSequence, ...
                    'frameToVisualize', 1, ...
                    'axesHandle', theBackgroundAxes, ...
                    'sRGBforSceneVisualization', true);
            axis(theBackgroundAxes, 'equal');
    
        
            % Superimpose the examined stimuli
            panelWidth = 0.085;
            panelR = 0.5-0.5*panelWidth;
        
            for iDir = 1:(skippedDirections+1):size(theLMSconeContrastDirections,2)
        
                % Create a grating scene engine with the examined chromatic direction
                theSceneEngine = createGratingSceneEngine(...
                    theLMSconeContrastDirections(:,iDir), examinedSpatialFrequencyCPD, ...
                    gratingSceneParams);
        
                % Compute the stimulus scene at max contrast
                testContrast = 1;
                theSceneSequence = theSceneEngine.compute(testContrast);
        
                % Visualize the stimulus scene
                figure(hFig);
                rmsC = norm(theLMSconeContrastDirections(:,iDir));
                theAxes = axes('Position', [panelR+panelR*0.95*theLMSconeContrastDirections(1,iDir)/rmsC panelR + panelR*0.95*theLMSconeContrastDirections(2,iDir)/rmsC panelWidth panelWidth]);
                theSceneEngine.visualizeStaticFrame(...
                    theSceneSequence, ...
                    'frameToVisualize', 1, ...
                    'axesHandle', theAxes, ...
                    'sRGBforSceneVisualization', true);
                axis(theAxes, 'equal');
                set(theAxes, 'Xcolor', 'none', 'YColor', 'none');
                
            end % iDir
        end

        if (isempty(theThresholdAxes))
            thresholdFigureHalfSize = 0.42;
            theThresholdAxes = axes(...
                'Position', [0.5-thresholdFigureHalfSize 0.5-thresholdFigureHalfSize thresholdFigureHalfSize*2 thresholdFigureHalfSize*2]);
            set(theThresholdAxes, 'XLim', maxThreshold*[-1 1], 'YLim', maxThreshold*[-1 1]);
            axis(theThresholdAxes, 'square');
            box(theThresholdAxes, 'off');
            set(theThresholdAxes, 'Color', 'none', 'Box', 'off', 'XColor', 'none', 'YColor', 'none');
        end
        
    end % if (initializeFig)


    thresholdConeContrasts(1,:) = cosd(theChromaticDirections) .* theThresholds;
    thresholdConeContrasts(2,:) = sind(theChromaticDirections) .* theThresholds;

    
    % Plot the axes
    plot(theThresholdAxes,  maxThreshold*[-1 1], [0 0], 'k-', 'LineWidth', 1.0);
    hold(theThresholdAxes, 'on');
    plot(theThresholdAxes,  [0 0], maxThreshold*[-1 1], 'k-', 'LineWidth', 1.0);

    % Plot the data points
    scatter(theThresholdAxes, theVisualizedConeContrastOffset(1) + thresholdConeContrasts(1,:), ...
        theVisualizedConeContrastOffset(2) + thresholdConeContrasts(2,:), 111, ...
        'MarkerEdgeColor', [0.99 0.0 0.0], 'MarkerFaceColor', [0.95 0.5 0.5], ...
        'LineWidth', 2.0, 'MarkerFaceAlpha', 0.7);

    if (numel(theChromaticDirections)>6)
        [z, a, b, rotationRadians] = fitEllipseToXYpoints(...
            thresholdConeContrasts, ...
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
        h = plot(theThresholdAxes, theVisualizedConeContrastOffset(1) + X(1,:), theVisualizedConeContrastOffset(2) + X(2,:), 'k-', 'LineWidth', 6.0);
        h = plot(theThresholdAxes, theVisualizedConeContrastOffset(1) + X(1,:), theVisualizedConeContrastOffset(2) + X(2,:), 'c-', 'LineWidth',2.0);
    end

    
    
    hold(theThresholdAxes, 'off');
    set(theThresholdAxes, 'XLim', maxThreshold*[-1 1], 'YLim', maxThreshold*[-1 1]);
    axis(theThresholdAxes, 'square');
    box(theThresholdAxes, 'off');
    set(theThresholdAxes, 'Color', 'none', 'Box', 'off', 'XColor', [0.1 0.1 0.1], 'YColor', [0.1 0.1 0.1]);
    xlabel(theThresholdAxes, 'threshold contrast (L-cone)');
    ylabel(theThresholdAxes, 'threshold contrast (M-cone)');
    set(theThresholdAxes, 'FontSize', 24);
    drawnow;

    if ((~isempty(figureFileBaseDir)) && (exportFig))
        NicePlot.exportFigToPNG(fullfile(figureFileBaseDir,'/stimuliOnLMplane.png'), hFig, 300);
    end

    %set(hFig, 'HandleVisibility', 'off');
end
