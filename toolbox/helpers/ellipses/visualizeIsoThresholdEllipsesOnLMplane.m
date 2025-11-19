function [thresholdDeltaConeContrasts, theThresholdAxes, theFittedEllipsePoints] = visualizeIsoThresholdEllipsesOnLMplane(...
    theLMSconeContrastDirections, ...
    theDeltaLMSconeContrastDirections, ...
    rmsLMconeContrast, ...
    examinedSpatialFrequencyCPD, ...
    gratingSceneParams, ...
    theChromaticAnglesDegs, ...
    skippedDirections, ...
    theThresholds, ...
    maxThreshold, ...
    referenceLMSconeContrast, ...
    figureFileBaseDir, hFig, ...
    exportFig, initializeFig, theThresholdAxes, figName)


    %gratingSceneParams.warningInsteadOfErrorOnOutOfGamut = true;

    gratingSceneParamsForNullScene = gratingSceneParams;
    if (~isempty(gratingSceneParamsForNullScene.backgroundLMSconeExcitations))
        gratingSceneParamsForNullScene.backgroundLMSconeExcitations = ...
                gratingSceneParamsForNullScene.backgroundLMSconeExcitations .* (1+referenceLMSconeContrast);
    end

    % Create the background scene engine
    theNullSceneEngine = createGratingSceneEngine(...
            [0 0 0], examinedSpatialFrequencyCPD, ...
            gratingSceneParams);

    % Compute the scene sequence
    theNullSceneSequence = theNullSceneEngine.compute(0.0);

    % Draw
    set(hFig, 'HandleVisibility', 'on');
    figure(hFig);
    
    plotColorCircle = ~true;

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
        
                % Compute the stimulus scene at a high contrast
                testContrast = 0.1;
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

    % Validate that they are by comparing their RMS values to what was loaded
    assert(all(sqrt(sum(theDeltaLMSconeContrastDirections.^2,1))-rmsLMconeContrast) == 0, 'not validating');
    

    thresholdDeltaConeContrasts(1,:) = cosd(theChromaticAnglesDegs) .* rmsLMconeContrast .* theThresholds;
    thresholdDeltaConeContrasts(2,:) = sind(theChromaticAnglesDegs) .* rmsLMconeContrast .* theThresholds;

    LconeContrastThresholds = thresholdDeltaConeContrasts(1,:) + referenceLMSconeContrast(1);
    MconeContrastThresholds = thresholdDeltaConeContrasts(2,:) + referenceLMSconeContrast(2);

    % Plot the axes
    hold(theThresholdAxes, 'on');
    plot(theThresholdAxes,  maxThreshold*[-1 1], [0 0], 'k-', 'LineWidth', 1.0);
    hold(theThresholdAxes, 'on');
    plot(theThresholdAxes,  [0 0], maxThreshold*[-1 1], 'k-', 'LineWidth', 1.0);



    % Plot the data points
    plot(theThresholdAxes, LconeContrastThresholds, MconeContrastThresholds, 'k-');
    scatter(theThresholdAxes, LconeContrastThresholds, MconeContrastThresholds, 111, ...
        'MarkerEdgeColor', [0.99 0.0 0.0], 'MarkerFaceColor', [0.95 0.5 0.5], ...
        'LineWidth', 2.0, 'MarkerFaceAlpha', 0.7);

    if (numel(theChromaticAnglesDegs)>6)
        [z, a, b, rotationRadians] = fitEllipseToXYpoints(...
            thresholdDeltaConeContrasts, ...
            'nonLinear', false);
        %[z, a, b, rotationRadians] = fitellipse(thresholdDeltaConeContrasts, 'linear');
    
        % Plot the fitted ellipse
        % form the parameter vector
        npts = 100;
        t = linspace(0, 2*pi, npts);
    
        % Rotation matrix
        Q = [cos(rotationRadians), -sin(rotationRadians); sin(rotationRadians) cos(rotationRadians)];
        % Ellipse points
        X = Q * [a * cos(t); b * sin(t)] + repmat(z, 1, npts);

        theFittedEllipsePoints = (bsxfun(@plus, referenceLMSconeContrast(1:2), X'));

        % PLot the ellipse points
        h = plot(theThresholdAxes, theFittedEllipsePoints(:,1), theFittedEllipsePoints(:,2), 'k-', 'LineWidth', 6.0);
        h = plot(theThresholdAxes, theFittedEllipsePoints(:,1), theFittedEllipsePoints(:,2), 'c-', 'LineWidth',2.0);
    end

    
    
    hold(theThresholdAxes, 'off');
    set(theThresholdAxes, 'XLim', maxThreshold*[-1 1], 'YLim', maxThreshold*[-1 1]);
    axis(theThresholdAxes, 'square');
    box(theThresholdAxes, 'off');
    set(theThresholdAxes, 'Color', 'none', 'Box', 'off', 'XColor', [0.1 0.1 0.1], 'YColor', [0.1 0.1 0.1]);
    xlabel(theThresholdAxes, 'L-cone contrast');
    ylabel(theThresholdAxes, 'M-cone contrast');
    set(theThresholdAxes, 'FontSize', 24);
    drawnow;

    if ((~isempty(figureFileBaseDir)) && (exportFig))
        NicePlot.exportFigToPNG(fullfile(figureFileBaseDir,sprintf('%s.png', figName)), hFig, 300, 'beVerbose');
    end

    set(hFig, 'HandleVisibility', 'off');
end
