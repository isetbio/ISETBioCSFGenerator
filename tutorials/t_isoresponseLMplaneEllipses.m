function t_isoresponseLMplaneEllipses
% Compute iso-response contours on the LM-plane
%
% Syntax:
%   t_isoresponseLMplaneEllipses
%
% Description:
%   Illustrates how to compute threshold iso-response ellipses on the LM
%   contrast plane using either CMosaic- or  mRGCmosaic- based response engines. 
%   In the case that an mRGCMosaic response engine is used, this tutorial also 
%   illustrates two additional options;
%   (a) how to simulate an mRGCmosaic with light increment (ON) and 
%       light-decrement (OFF) neurons, although we currently only have ON-center mRGCmosaics
%   (b) how to simulate a non-linear activation function 
%
%   This tutorial was used to generate the examples ellipses that were
%   included in the computational aim of the 2025 R-01 grant.
%
% History:
%   01/27/2025  NPC   Wrote it.

%% Close any stray figs
hAllFigs = findall(groot,'Type','figure');

% Close all figures
for i = 1:numel(hAllFigs)
    set(hAllFigs(i), 'HandleVisibility', 'on')
end

close all;


%% Freeze rng for replicatbility and validation
rng(1);


debugStimulusConfig = ~true;

verbose = true;
useMetaContrast = true; 
useConeContrast = true;
useFixationalEMs = false;

% Cone mosaic performance params
whichNoiseFreeNre = 'excitationsCmosaic';
whichNoisyInstanceNre = 'Poisson';
whichClassifierEngine = 'rcePoisson';
logThreshLimitLow =  4.0;
logThreshLimitHigh = 1.0;

% mRGC mosaic performance params
mRGCOutputSignalType = 'mRGCs';         % Select between {'cones', 'mRGCs'}

whichNoiseFreeNre = 'mRGCMosaic';
whichNoisyInstanceNre = 'Gaussian';
whichClassifierEngine = 'rceTemplateDistance'; 
logThreshLimitLow =  4.0;
logThreshLimitHigh = 0.0;

% No temporal filter
temporalFilter = [];

% For best results, set the testTrials >= 1024
testTrials = 256;

thresholdParams = struct( ...
        'logThreshLimitLow', logThreshLimitLow , ...
        'logThreshLimitHigh', logThreshLimitHigh, ...
        'logThreshLimitDelta', 0.01, ...
        'slopeRangeLow', 0.1, ...
        'slopeRangeHigh', 3, ...
        'slopeDelta', 0.05, ...
        'thresholdCriterion', 0.81606);

% Optical image padding
oiPadMethod = 'zero';

% Visualization options
visualizeEachScene = false;
visualizeEachCompute = false;
maxVisualizedNoisyResponseInstances = 2;

figureFileBaseDir = setupFigureDirectory(mfilename, ...
    useMetaContrast, useConeContrast, useFixationalEMs, ...
    whichNoiseFreeNre, whichNoisyInstanceNre,...
    whichClassifierEngine, mRGCOutputSignalType);

% Specify stimuli as uniform fields, with a background lum of 
% 100 cd/m2, (xy) chroma of [0.3 0.32], that last for presented for 
% a single frame that lasts for 50 msec. 
examinedSpatialFrequencyCPD = 0;
gratingSceneParams = struct( ...
        'meanLuminanceCdPerM2', 100, ...
        'meanChromaticityXY', [0.30 0.32], ...
        'presentationMode', 'flashed', ... 
        'spatialEnvelope', 'disk', ...
        'duration', 50/1000, ...
        'frameDurationSeconds', 50/1000, ...
        'pixelsNum', 256, ...
        'fovDegs', 1.0, ...
        'spatialEnvelopeRadiusDegs', 0.5, ...
        'spectralSupport', 400:20:750 ...
        );

% Max RMS contrast (so as to keep stimuli within the display gamut)
rmsLMconeContrast = 0.12;

if (debugStimulusConfig)
    % Alter params to maximize the visibility of response
    gratingSceneParams.warningInsteadOfErrorOnOutOfGamut = true;
    gratingSceneParams.spatialEnvelope = 'disk';
    rmsLMconeContrast = 1.0;
    thresholdParams.logThreshLimitLow = 0.1;
end

% Compute the test LMS cone contrast directions

% HiRes RUN
examinedDirectionsOnLMplane = [];

% FAST RUN
examinedDirectionsOnLMplane = 0:45:345;

[theLMSconeContrastDirections, examinedDirectionsOnLMplane] = ...
    computeLMSconeContrastDirections(rmsLMconeContrast, examinedDirectionsOnLMplane);

% Plot all stimuli
skippedDirections = 0;
if (length(examinedDirectionsOnLMplane) > 50)
    % Too many, plot every other stimulus
    skippedDirections = 2;
end



%% Configure our stimulus scene engines 
[theNullSceneEngine, theTestSceneEngines] = configureStimulusSceneEngines(...
    theLMSconeContrastDirections, examinedSpatialFrequencyCPD, gratingSceneParams);


%% Neural response engines
%

% Setup the noise-free neural response engine
switch (whichNoiseFreeNre)
    case 'excitationsCmosaic'
        % Set the compute function
        nreNoiseFreeComputeFunction = @nreNoiseFreeCMosaic;

        % Get default params struct
        nreNoiseFreeParams = nreNoiseFreeCMosaic([],[],[],[], ...
            'oiPadMethod',oiPadMethod);

        % Modify mosaic size (smaller, 85%, that the stimulus)
        nreNoiseFreeParams.coneMosaicParams.sizeDegs = gratingSceneParams.fovDegs * 0.85 * [1 1];

        % Set the cone mosaic integration time to match the stimulus frame duration
        nreNoiseFreeParams.coneMosaicParams.timeIntegrationSeconds = gratingSceneParams.frameDurationSeconds;

        % Handle cone contrast setting
        if (useConeContrast)
            nreNoiseFreeParams.coneMosaicParams.outputSignalType = 'coneContrast';
        else
            nreNoiseFreeParams.coneMosaicParams.outputSignalType = 'coneExcitations';
        end

        % No temporal filter
        nreNoiseFreeParams.temporalFilter = temporalFilter;

        % Custom nre visualization function, or set to []
        customVisualizationFunctionHandle = @nreVisualizeCMosaic;

    case 'mRGCMosaic'
        % Set the compute function
        nreNoiseFreeComputeFunction = @nreNoiseFreeMidgetRGCMosaic;

        % Get default params struct
        nreNoiseFreeParams = nreNoiseFreeMidgetRGCMosaic([],[],[],[], ...
            'oiPadMethod',oiPadMethod);
        

        % Select one of the pre-computed mRGC mosaics by specifying its 
        % eccentricity, size, and type (currently only 'ONcenterMidgetRGC')
        nreNoiseFreeParams.mRGCMosaicParams.eccDegs = [0 0];
        nreNoiseFreeParams.mRGCMosaicParams.sizeDegs = [2.0 2.0];
        nreNoiseFreeParams.mRGCMosaicParams.rgcType = 'ONcenterMidgetRGC';

        % Modify mosaic size (smaller, 85%, that the stimulus) and integration time to match the stimulus
        cropSize = gratingSceneParams.fovDegs * 0.85 * [1 1];
        nreNoiseFreeParams.mRGCMosaicParams.cropParams = struct(...
            'sizeDegs', cropSize, ...
            'eccentricityDegs', [] ...
            );

        % Set the input cone mosaic integration time to match the stimulus frame duration
        nreNoiseFreeParams.mRGCMosaicParams.coneIntegrationTimeSeconds = gratingSceneParams.frameDurationSeconds;

        % Specify that mRGCs will operate on cone-contrast responses
        nreNoiseFreeParams.mRGCMosaicParams.inputSignalType = 'coneContrast';

        % Since we specified that mRGCs will operate on cone-contrast responses,
        % will need a null scene for normalization
        nreNoiseFreeParams.nullStimulusSceneSequence = theNullSceneEngine.compute(0.0);

        % No temporal filter
        nreNoiseFreeParams.temporalFilter = temporalFilter;

        % Simulate incomplete background adaptation - puts mosaic in
        % saturation regime of the activation function
        % This has an effect only if we add an 'activationFunctionParams'
        % struct field in 'nreNoiseInstancesParams', whose 'type' field is set to 'halfwaveSigmoidalRectifier'
        nreNoiseFreeParams.mRGCMosaicParams.responseBias = 0.005;
        theVisualizedConeContrastOffset = [0.0 0.0];

        % Puts mosaic in heavy saturation regime
        %nreNoiseFreeParams.mRGCMosaicParams.responseBias = 0.05;
        %theVisualizedConeContrastOffset = [0.1 0.1];

        % Simulate ON-OFF mosaic
        % All odd-indexes mRGCs will be treated as OFF-center, by inverting their 
        % noise-free responses polarity in nreNoiseFreeMidgetRGCMosaic()
        nreNoiseFreeParams.mRGCMosaicParams.simulateONOFFmosaic = true;

        % Custom nre visualization function, or set to []
        customVisualizationFunctionHandle = @nreVisualizeMRGCmosaic;

    otherwise
        error('Unsupported noise free neural response engine: ''%s''.', whichNoiseFreeNre);
end  % switch (whichNoiseFreeNre)

% Setup the noisy neural response engine

switch (whichNoisyInstanceNre)
    case 'Poisson'
        nreNoisyInstancesComputeFunction = @nreNoisyInstancesPoisson;
        nreNoiseInstancesParams = nreNoisyInstancesPoisson;
       
        % No activation function for cone mosaic excitations
        %nreNoiseInstancesParams.activationFunctionParams = [];

    case 'Gaussian'
        nreNoisyInstancesComputeFunction = @nreNoisyInstancesGaussian;
        nreNoiseInstancesParams = nreNoisyInstancesGaussian;

        % Set Gaussian sigma depending on the noise-free NRE
        switch (whichNoiseFreeNre)
            case 'excitationsCmosaic'
                gaussianSigma = 50;

            case 'mRGCMosaic'
                switch (mRGCOutputSignalType)
                    case 'mRGCs'
                        gaussianSigma = 20 * 1e-3;
                    case 'cones'
                        gaussianSigma = 50;
                    otherwise
                        error('Unknown mRGC output signal type specified');
                end

            otherwise
                error('Unsupported noise free neural response engine: ''%s''.', whichNoiseFreeNre);
        end  % switch (whichNoiseFreeNre)


        % Set the sigma parameter
        nreNoiseInstancesParams.sigma = gaussianSigma;

        % Saturating activation function. If the type is specified as
        % 'halfWaveSigmoidalRectifier', the response bias specified above
        % like so: nreNoiseFreeParams.mRGCMosaicParams.responseBias
        % specifies how the noise-free response is pushed into different
        % regimes of the sigmoidal function
        nreNoiseInstancesParams.activationFunctionParams = struct(...
            'type', 'halfwaveSigmoidalRectifier', ...                 % choose between {'linear', 'halfwaveRectifier', 'halfwaveSigmoidalRectifier'}
            'exponent', 2.0, ...                                      % only relevant for 'halfwaveSigmoidalRectifier'
            'semiSaturationReponseAmplitude', 0.15*gaussianSigma, ... % only relevant for 'halfwaveSigmoidalRectifier'
            'gain', 3*gaussianSigma, ...                              % only relevant for 'halfwaveSigmoidalRectifier'
            'visualize', visualizeEachCompute ...
            );


    otherwise
        error('Unsupported noisy instances neural response engine: ''%s''.', whichNoisyInstanceNre);
end % switch (whichNoisyInstanceNre)



%% If we use cone contrast, we will neeed a null scene for normalization.
if (useConeContrast) 
    nreNoiseFreeParams.nullStimulusSceneSequence = theNullSceneEngine.compute(0.0);
end

%% Create the neural engine
theNeuralEngine = neuralResponseEngine( ...
    nreNoiseFreeComputeFunction, ...
    nreNoisyInstancesComputeFunction, ...
    nreNoiseFreeParams, ...
    nreNoiseInstancesParams ...
    );

% Set the neuralEngine's various visualization properties
theNeuralEngine.visualizeEachCompute = visualizeEachCompute;
theNeuralEngine.maxVisualizedNoisyResponseInstances = maxVisualizedNoisyResponseInstances;
theNeuralEngine.customVisualizationFunctionHandle = customVisualizationFunctionHandle;

%% If using meta contrast, set this up.
if (useMetaContrast)
    metaSceneEngineParams = sceMetaContrast;
    theMetaSceneEngine = sceneEngine(@sceMetaContrast,metaSceneEngineParams);

    % Create nreMetaContrast using the actual scene and neural engines
    metaNeuralResponseEngineNoiseFreeParams = nreNoiseFreeMetaContrast;
    metaNeuralResponseEngineNoiseFreeParams.contrast0 = 0;
    metaNeuralResponseEngineNoiseFreeParams.contrast1 = 1;
    metaNeuralResponseEngineNoiseFreeParams.neuralEngine = theNeuralEngine;

    metaNeuralResponseEngineNoisyInstanceParams = nreNoisyInstancesMetaContrast;
    metaNeuralResponseEngineNoisyInstanceParams.neuralEngine = theNeuralEngine;
end

if (useFixationalEMs) 
    error('Need to set up fEM')
else
    trainFixationalEMObj = [];
    testFixationalEMObj = [];
end

%% Instantiate the responseClassifierEngine
%
% rcePoisson makes decision by performing the Poisson likelihood ratio
% test. This is the ideal observer for the Poisson noise cone excitations
%
% Set up parameters associated with use of this classifier.
switch (whichClassifierEngine)
    case {'rcePoisson'}
        theClassifierEngine = responseClassifierEngine(@rcePoisson);
        classifierEngineParams = struct(...
            'trainFlag', 'none', ...
            'testFlag', 'random', ...
            'nTrain', 1, 'nTest', testTrials);

    case {'rceTemplateDistance'}
        theClassifierEngine = responseClassifierEngine(@rceTemplateDistance);
        classifierEngineParams = struct(...
            'trainFlag', 'none', ...
            'testFlag', 'random', ...
            'nTrain', 1, 'nTest', testTrials);

    otherwise
        error('Unsupported response classifier engine: ''%s''.', whichClassifierEngine)
end


%% Parameters for threshold estimation/quest engine
questEngineParams = struct( ...
    'qpPF',@qpPFWeibullLog, ...
    'minTrial', 1280, ...
    'maxTrial', 1280, ...
    'numEstimator', 1, ...
    'stopCriterion', 0.05);

%% Setup figures
psychometricCurvesFig = figure(2345); clf;
set(psychometricCurvesFig, 'Position', [10 10 1200 1200]);
axLeft = cell(1, length(examinedDirectionsOnLMplane));
axRight = cell(1, length(examinedDirectionsOnLMplane));
axLeftRight = cell(1, length(examinedDirectionsOnLMplane));
for iChromaDirection = 1:length(examinedDirectionsOnLMplane)

    margin = 0.01;
    cols = 5;
    width = (1-(cols+2)*margin) / cols;
    rows = ceil(length(examinedDirectionsOnLMplane)/cols);
    height = (1-(rows+2)*margin) / rows; 
    bottom = 2*margin + floor((iChromaDirection-1)/cols) * (height + margin);
    if (mod(iChromaDirection-1,cols) == 0)
        axLeftRight{iChromaDirection} = axes('Position', [2*margin bottom width height]);
    elseif (mod(iChromaDirection-1,cols) == 1)
        axLeftRight{iChromaDirection} = axes('Position', [3*margin+width bottom width height]);
    elseif (mod(iChromaDirection-1,cols) == 2)
        axLeftRight{iChromaDirection} = axes('Position', [4*margin+2*width bottom width height]);
    elseif (mod(iChromaDirection-1,cols) == 3)
        axLeftRight{iChromaDirection} = axes('Position', [5*margin+3*width bottom width height]);
    else
        axLeftRight{iChromaDirection} = axes('Position', [6*margin+4*width bottom width height]);
    end

end % for iChromaDirection
set(psychometricCurvesFig, 'HandleVisibility', 'off');

% Figure for plotting the thresholds on the LM plane
% along with the stimuli
hFigStimuliAndThresholds = figure(2346); clf;
set(hFigStimuliAndThresholds, 'Color', [1 1 1], 'Position', [10 10 1200 1200]);
set(hFigStimuliAndThresholds, 'HandleVisibility', 'off');

% Max visualized threshold
maxVisualizedThreshold = 0.2; % 0.05;


recomputeThresholds = true;
if (recomputeThresholds)

    logThreshold = zeros(1, length(examinedDirectionsOnLMplane));
    
    %% Compute psychometric curves for each chromatic direction
    for iChromaDirection = 1:length(examinedDirectionsOnLMplane)
        theSceneEngine = theTestSceneEngines{iChromaDirection};
    
        % Set the sceneEngine's visualizeEachCompute property
        theSceneEngine.visualizeEachCompute = visualizeEachScene;
    
        if (useMetaContrast)
            % Instantiate meta contrast neural engine for this spatial
            % frequency and use it to compute threshold
            metaNeuralResponseEngineNoiseFreeParams.sceneEngine = theSceneEngine;
    
            theMetaNeuralEngine = neuralResponseEngine(...
                @nreNoiseFreeMetaContrast, ...
                @nreNoisyInstancesMetaContrast, ...
                metaNeuralResponseEngineNoiseFreeParams, ...
                metaNeuralResponseEngineNoisyInstanceParams);
            
            % Update visualizeEachCompute
            theMetaNeuralEngine.visualizeEachCompute = theNeuralEngine.visualizeEachCompute;
    
            % Compute the threshold for our grating scene with meta scene and
            % and neural response engines. This function does a lot of work,
            % see the function itself, as well as function computePerformance.
            [logThreshold(iChromaDirection), questObj, ~, psychometricCurveParams(iChromaDirection,:)] = ...
                computeThreshold(theMetaSceneEngine, theMetaNeuralEngine, theClassifierEngine, ...
                classifierEngineParams, thresholdParams, questEngineParams, ...
                'TAFC', true, 'useMetaContrast', useMetaContrast, ...
                'trainFixationalEM', trainFixationalEMObj, ...
                'testFixationalEM', testFixationalEMObj, ...
                'verbose', verbose, ...
                'maxVisualizedNoisyResponseInstances', maxVisualizedNoisyResponseInstances);
    
        end % useMetaContrast
    
        % Plot stimulus & psychometric curve
        set(psychometricCurvesFig, 'HandleVisibility', 'on');
        % Plot data and fitted psychometric curve
        questObj.plotMLE(2.5,'para', psychometricCurveParams(iChromaDirection,:), ...
            'axesHandle', axLeftRight{iChromaDirection});  %axRight{iChromaDirection});
        drawnow;
        set(psychometricCurvesFig, 'HandleVisibility', 'off');
    
    
        if (iChromaDirection == length(examinedDirectionsOnLMplane))

            % Visualize thresholds up to this point
            theChromaticDirections = examinedDirectionsOnLMplane(1:iChromaDirection);
            theThresholds = 10 .^ (logThreshold(1:iChromaDirection));
        
            % Initialization
            exportFig = (iChromaDirection == length(examinedDirectionsOnLMplane));
            initializeFig = true;
            if (initializeFig)
                theThresholdAxes = [];
            end
    
            theThresholdAxes = visualizeStimuliAndThresholdsOnLMPlane(...
                rmsLMconeContrast, ...
                examinedSpatialFrequencyCPD, gratingSceneParams, ...
                theChromaticDirections, skippedDirections, theThresholds, maxVisualizedThreshold, theVisualizedConeContrastOffset, ...
                figureFileBaseDir, hFigStimuliAndThresholds, ...
                exportFig, initializeFig, theThresholdAxes);
        end
    end % iChromaDirection

else
    % Load hard-wired, previously computed thresholds
    [theThresholds, theVisualizedConeContrastOffset] = previouslyComputedThresholds();
    exportFig = true;
    initializeFig = true;
    
    visualizeStimuliAndThresholdsOnLMPlane(...
            rmsLMconeContrast, ...
            examinedSpatialFrequencyCPD, gratingSceneParams, ...
            examinedDirectionsOnLMplane, skippedDirections, ...
            theThresholds, maxVisualizedThreshold, theVisualizedConeContrastOffset, ...
            figureFileBaseDir, hFigStimuliAndThresholds, ...
            exportFig, initializeFig, []);
end % recomputeThresholds


end % - function end


%
% SUPPORTING FUNCTIONS
%


function [LMSconeContrasts, examinedDirectionsOnLMplane] = computeLMSconeContrastDirections(...
    rmsLMconeContrast, examinedDirectionsOnLMplane)

    if (isempty(examinedDirectionsOnLMplane))
        examinedDirectionsOnLMplane = [0:10:20 25:5:65 70:10:170];
        examinedDirectionsOnLMplane = cat(2, examinedDirectionsOnLMplane, 180+examinedDirectionsOnLMplane);
    end

    LMSconeContrasts = zeros(3, numel(examinedDirectionsOnLMplane));
    if (numel(rmsLMconeContrast) == numel(examinedDirectionsOnLMplane))
        LMSconeContrasts(1,:) = rmsLMconeContrast .* cosd(examinedDirectionsOnLMplane);
        LMSconeContrasts(2,:) = rmsLMconeContrast .* sind(examinedDirectionsOnLMplane);
    else
        LMSconeContrasts(1,:) = rmsLMconeContrast(1) * cosd(examinedDirectionsOnLMplane);
        LMSconeContrasts(2,:) = rmsLMconeContrast(1) * sind(examinedDirectionsOnLMplane);
    end
end

function theThresholdAxes = visualizeStimuliAndThresholdsOnLMPlane(...
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

    set(hFig, 'HandleVisibility', 'off');
end


function [theNullSceneEngine, theTestSceneEngines] = configureStimulusSceneEngines(...
    theLMSconeContrastDirections, examinedSpatialFrequencyCPD, gratingSceneParams)

    % Create the background scene engine
    theNullSceneEngine = createGratingSceneEngine(...
            [0 0 0.4], examinedSpatialFrequencyCPD, ...
            gratingSceneParams);

    % Store test scene engines
    theTestSceneEngines = cell(1,size(theLMSconeContrastDirections,2));

    for iDir = 1:size(theLMSconeContrastDirections,2)
        % Create a grating scene engine with the examined chromatic direction
        theTestSceneEngines{iDir} = createGratingSceneEngine(...
            theLMSconeContrastDirections(:,iDir), examinedSpatialFrequencyCPD, ...
            gratingSceneParams);
    end % iDir

end




function figureFileBaseDir = setupFigureDirectory(theScriptName, ...
    useMetaContrast,useConeContrast,useFixationalEMs,...
    whichNoiseFreeNre,whichNoisyInstanceNre,...
    whichClassifierEngine,mRGCOutputSignalType)

    % Make sure local/figures directory exists so we can write out our figures in peace
    projectBaseDir = ISETBioCSFGeneratorRootPath;
    if (~exist(fullfile(projectBaseDir,'local', theScriptName, 'figures'),'dir'))
        mkdir(fullfile(projectBaseDir,'local', theScriptName, 'figures'));
        fprintf('Generated figure directory at %s\n', fullfile(projectBaseDir,'local', theScriptName, 'figures'))
    end

    figureFileBaseDir = fullfile(projectBaseDir,'local',mfilename,'figures', ...
        sprintf('%s_Meta_%d_ConeContrast_%d_FEMs_%d_%s_%s_%s_%s', mfilename, ...
        useMetaContrast,useConeContrast,useFixationalEMs,whichNoiseFreeNre,whichNoisyInstanceNre,...
        whichClassifierEngine,mRGCOutputSignalType));
    if (~exist(figureFileBaseDir, 'dir'))
        mkdir(figureFileBaseDir);
        fprintf('Generated figure sub-directory at %s\n', figureFileBaseDir);
    end


end


function [theThresholds, theVisualizedConeContrastOffset] = previouslyComputedThresholds()

theConeMosaicThresholds = [ ...
        0.0092 ...
    0.0089 ...
    0.0088 ...
    0.0101 ...
    0.0085 ...
    0.0096 ...
    0.0099 ...
    0.0107 ...
    0.0109 ...
    0.0128 ...
    0.0116 ...
    0.0126 ...
    0.0133 ...
    0.0132 ...
    0.0160 ...
    0.0108 ...
    0.0137 ...
    0.0102 ...
    0.0115 ...
    0.0083 ...
    0.0089 ...
    0.0085 ...
    0.0083 ...
    0.0085 ...
    0.0090 ...
    0.0080 ...
    0.0080 ...
    0.0092 ...
    0.0100 ...
    0.0098 ...
    0.0105 ...
    0.0123 ...
    0.0116 ...
    0.0126 ...
    0.0122 ...
    0.0125 ...
    0.0148 ...
    0.0135 ...
    0.0137 ...
    0.0128 ...
    0.0121 ...
    0.0109 ...
    0.0094 ...
    0.0101 ...
    0.0099 ...
    0.0091 ...
        ];

    theThresholdsLowSaturation = [ ...
         0.0077  ...
    0.0086 ...
    0.0109 ...
    0.0142 ...
    0.0152 ...
    0.0185 ...
    0.0175 ...
    0.0206 ...
    0.0201 ...
    0.0185 ...
    0.0177 ...
    0.0130 ...
    0.0122 ...
    0.0093 ...
    0.0078 ...
    0.0064 ...
    0.0059 ...
    0.0057 ...
    0.0048 ...
    0.0056 ...
    0.0060 ...
    0.0065 ...
    0.0068 ...
    0.0071 ...
    0.0093 ...
    0.0115 ...
    0.0155 ...
    0.0149 ...
    0.0174 ...
    0.0187 ...
    0.0234 ...
    0.0212 ...
    0.0198 ...
    0.0161 ...
    0.0144 ...
    0.0114 ...
    0.0083 ...
    0.0078 ...
    0.0061 ...
    0.0063 ...
    0.0061 ...
    0.0052 ...
    0.0048 ...
    0.0057 ...
    0.0059 ...
    0.0074 ...
        ];


    theThresholdsHeavySaturation = [...
        0.0208 ...
    0.0250 ...
    0.0335 ...
    0.0374 ...
    0.0430 ...
    0.0487 ...
    0.0540 ...
    0.0593 ...
    0.0543 ...
    0.0515 ...
    0.0426 ...
    0.0345 ...
    0.0341 ...
    0.0283 ...
    0.0215 ...
    0.0177 ...
    0.0171 ...
    0.0164 ...
    0.0166 ...
    0.0175 ...
    0.0138 ...
    0.0172 ...
    0.0195 ...
    0.0202 ...
    0.0257 ...
    0.0291 ...
    0.0348 ...
    0.0436 ...
    0.0475 ...
    0.0532 ...
    0.0597 ...
    0.0500 ...
    0.0430 ...
    0.0403 ...
    0.0364 ...
    0.0311 ...
    0.0256 ...
    0.0185 ...
    0.0174 ...
    0.0175 ...
    0.0138 ...
    0.0151 ...
    0.0143 ...
    0.0169 ...
    0.0173 ...
    0.0186 ...
        ];

    theThresholds = theThresholdsLowSaturation;
    theVisualizedConeContrastOffset = [0.0 0.0];

    %theThresholds = theThresholdsHeavySaturation;
    %theVisualizedConeContrastOffset = [0.15 0.15];

    %theThresholds = theConeMosaicThresholds
    %theVisualizedConeContrastOffset = [0.0 0.0];

end


function [z, a, b, alpha] = fitEllipseToXYpoints(xyPoints, varargin)

	% Parse input
    p = inputParser;
    p.addParameter('maxIterations', 200, @isnumeric);
    p.addParameter('tolerance', 1e-5, @isnumeric);
    p.addParameter('constraint', 'bookstein', @(x)(ismember(x, {'bookstein', 'trace'})));
    p.addParameter('nonLinear', true, @islogical);
    p.parse(varargin{:});

	rows = size(xyPoints,1);
	cols = size(xyPoints,2);
	assert(rows == 2, 'xyPoints must be a 2 x N matrix');
	assert(cols >=6, 'xyPoints must have at least 6 points');

	fitParams = struct(...
		'nonLinear', p.Results.nonLinear, ...
		'constraint', p.Results.constraint, ...
		'maxIterations', p.Results.maxIterations, ...
		'tolerance', p.Results.tolerance);

	% Remove centroid
	centroid = mean(xyPoints, 2);
	xyPoints = bsxfun(@minus, xyPoints, centroid);

	% Obtain a linear estimate
	switch fitParams.constraint
    	case 'bookstein'
        	[z, a, b, alpha] = fitbookstein(xyPoints);
    	case 'trace'
       		[z, a, b, alpha] = fitggk(xyPoints);
	end % switch


	if (fitParams.nonLinear)
		% Initial conditions
	    z0     = z;
	    a0     = a;
	    b0     = b;
	    alpha0 = alpha;

	    % Fit
    	[z, a, b, alpha, converged, isCircle] = fitNonLinear(xyPoints, z0, a0, b0, alpha0, params);

    	% Return linear estimate if GN doesn't converge or if the data points fall on a circle
	    if (~converged) || (isCircle)
	        fprintf('*** FailureToConverge: Gauss-Newton did not converge, returning linear estimate.\n');
	        z = z0;
	        a = a0;
	        b = b0;
	        alpha = alpha0;
	    end
	end % if (fitParams.nonLinear)

	% Add the centroid back on
	z = z + centroid;
end


function [z, a, b, alpha] = fitbookstein(x)
	%FITBOOKSTEIN   Linear ellipse fit using bookstein constraint
	%   lambda_1^2 + lambda_2^2 = 1, where lambda_i are the eigenvalues of A

	% Convenience variables
	m  = size(x, 2);
	x1 = x(1, :)';
	x2 = x(2, :)';

	% Define the coefficient matrix B, such that we solve the system
	% B *[v; w] = 0, with the constraint norm(w) == 1
	B = [x1, x2, ones(m, 1), x1.^2, sqrt(2) * x1 .* x2, x2.^2];

	% To enforce the constraint, we need to take the QR decomposition
	[Q, R] = qr(B);

	% Decompose R into blocks
	R11 = R(1:3, 1:3);
	R12 = R(1:3, 4:6);
	R22 = R(4:6, 4:6);

	% Solve R22 * w = 0 subject to norm(w) == 1
	[U, S, V] = svd(R22);
	w = V(:, 3);

	% Solve for the remaining variables
	v = -R11 \ R12 * w;

	% Fill in the quadratic form
	A        = zeros(2);
	A(1)     = w(1);
	A([2 3]) = 1 / sqrt(2) * w(2);
	A(4)     = w(3);
	bv       = v(1:2);
	c        = v(3);

	% Find the parameters
	[z, a, b, alpha] = conic2parametric(A, bv, c);
end

function [z, a, b, alpha] = fitggk(x)
	% Linear least squares with the Euclidean-invariant constraint Trace(A) = 1
	% Convenience variables
	m  = size(x, 2);
	x1 = x(1, :)';
	x2 = x(2, :)';

	% Coefficient matrix
	B = [2 * x1 .* x2, x2.^2 - x1.^2, x1, x2, ones(m, 1)];
	v = B \ -x1.^2;

	% For clarity, fill in the quadratic form variables
	A        = zeros(2);
	A(1,1)   = 1 - v(2);
	A([2 3]) = v(1);
	A(2,2)   = v(2);
	bv       = v(3:4);
	c        = v(5);

	% find parameters
	[z, a, b, alpha] = conic2parametric(A, bv, c);
end


function [z, a, b, alpha, converged, isCircle] = fitNonLinear(x, z0, a0, b0, alpha0, params)
	% Gauss-Newton least squares ellipse fit minimising geometric distance 

	% Get initial rotation matrix
	Q0 = [cos(alpha0), -sin(alpha0); sin(alpha0) cos(alpha0)];
	m = size(x, 2);

	% Get initial phase estimates
	phi0 = angle( [1 i] * Q0' * (x - repmat(z0, 1, m)) )';
	u = [phi0; alpha0; a0; b0; z0];

	% Iterate using Gauss Newton
	converged = false;

	for nIts = 1:params.maxIterations
	    % Find the function and Jacobian
	    [f, J, isCircle] = computeJacobian(u);
    
    	if (isCircle)
    		fprintf('Ellipse is near-circular - nonlinear fit may not succeed\n.')
    	end

	    % Solve for the step and update u
	    h = -J \ f;
	    u = u + h;
    
	    % Check for convergence
	    delta = norm(h, inf) / norm(u, inf);
	    if delta < params.tolerance
	        converged = true;
	        break
	    end
	end

	alpha = u(end-4);
	a = u(end-3);
	b = u(end-2);
	z = u(end-1:end);

	% ---- Nested function ---
	function [f, J, isCircle] = computeJacobian(u)
        % Define the system of nonlinear equations and Jacobian. 

        % Tolerance for whether it is a circle
        circTol = 1e-5;
        
        % Unpack parameters from u
        phi   = u(1:end-5);
        alpha = u(end-4);
        a     = u(end-3);
        b     = u(end-2);
        z     = u(end-1:end);
        
        % If it is a circle, the Jacobian will be singular, and the
        % Gauss-Newton step won't work. 
        %TODO: This can be fixed by switching to a Levenberg-Marquardt
        %solver
        if (abs(a - b) / (a + b) < circTol)
            isCircle = true;
        else
        	isCircle = false;
        end

        % Convenience trig variables
        c = cos(phi);
        s = sin(phi);
        ca = cos(alpha);
        sa = sin(alpha);
        
        % Rotation matrices
        Q    = [ca, -sa; sa, ca];
        Qdot = [-sa, -ca; ca, -sa];

        % Preallocate function and Jacobian variables
        f = zeros(2 * m, 1);
        J = zeros(2 * m, m + 5);
        for i = 1:m
            rows = (2*i-1):(2*i);
            % Equation system - vector difference between point on ellipse
            % and data point
            f((2*i-1):(2*i)) = x(:, i) - z - Q * [a * cos(phi(i)); b * sin(phi(i))];
            
            % Jacobian
            J(rows, i) = -Q * [-a * s(i); b * c(i)];
            J(rows, (end-4:end)) = ...
                [-Qdot*[a*c(i); b*s(i)], -Q*[c(i); 0], -Q*[0; s(i)], [-1 0; 0 -1]];
        end
    end % ---- Nested function ---
end

function [z, a, b, alpha] = conic2parametric(A, bv, c)
	% Diagonalise A - find Q, D such at A = Q' * D * Q
	[Q, D] = eig(A);
	Q = Q';

	% If the determinant < 0, it's not an ellipse
	if prod(diag(D)) <= 0 
	    error('NotEllipse', 'Linear fit did not produce an ellipse');
	end

	% We have b_h' = 2 * t' * A + b'
	t = -0.5 * (A \ bv);

	c_h = t' * A * t + bv' * t + c;

	z = t;
	a = sqrt(-c_h / D(1,1));
	b = sqrt(-c_h / D(2,2));
	alpha = atan2(Q(1,2), Q(1,1));
end % conic2parametric
