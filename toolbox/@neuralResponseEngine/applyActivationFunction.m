function theNeuralResponses = applyActivationFunction(theNoiseFreeResponses, theNoisyResponseInstances, activationFunctionParams)
% Function for applying an activation function to the noisy response instances
%
% Syntax:
%   theNeuralResponses = applyActivationFunction(theNoiseFreeResponses, ...
%                           theNoisyResponseInstances, activationFunctionParams);
%
% Description:
%   Function for applying an activation function to the noisy response
%   instances. The implicit assumption here is that the activation function
%   is applied after noise is added to the noise-free response. For
%   example when the noisy-response instances represent the noisy membrane
%   potential of ganglion cells, and the activation function represents the
%   action potential generation mechanism. The function is called from 
%   noisy response instancen compute functions, currently reNoisyInstancesGaussian
%   if an activationFunctionParams struct has been specified
%
% Inputs:
%   - theNoiseFreeResponses           - [1 x mNeurons] matrix of noise-free responses (only used for visualizing the activation function)
%   - theNoisyResponseInstances       - [nTrials x mNeurons] matrix of noisy-response instances
%   - activationFunctionParams        - struct with the activation function params
%
% Outputs:
%   - theNeuralResponses              - [nTrials x mNeurons] matrix of noisy-response instances after the activation function is applied
%
% Usage:
%   This function is triggered if an activationFunctionParams struct has
%   been specified. For an example see t_isoresponseLMplaneEllipses
%

% History:
%    01/27/2025  NPC   Wrote it.

    if (isempty(activationFunctionParams))
        theNeuralResponses = theNoisyResponseInstances;
        return;
    end

    % Save for visualization
    if (activationFunctionParams.visualize)
        thePreActivationFunctionNoisyResponseInstances = theNoisyResponseInstances;
    end

    switch (activationFunctionParams.type)
        case 'linear'
            theNeuralResponses = theNoisyResponseInstances;

        case 'halfwaveRectifier'
            theNeuralResponses = theNoisyResponseInstances;

            % 2. half-wave rectify
            theNeuralResponses(theNeuralResponses(:)<0) = 0;

        case 'halfwaveSigmoidalRectifier'
            theNeuralResponses = theNoisyResponseInstances;

            % 2. half-wave rectify
            theNeuralResponses(theNeuralResponses(:)<0) = 0; 

            % 3. Naka-Ruston sigmoidal activation
            theNeuralResponses = (theNeuralResponses.^activationFunctionParams.exponent) ./ ...
                (theNeuralResponses.^activationFunctionParams.exponent + (activationFunctionParams.semiSaturationReponseAmplitude)^activationFunctionParams.exponent);

            % 4. gain
            theNeuralResponses = activationFunctionParams.gain * theNeuralResponses;
            
        otherwise
            error('Unsupported activation function type: ''%s''.', activationFunctionParams.type);
    end % noiseFreeComputeParams.mRGCMosaicParams.activationFunctionType


    if (activationFunctionParams.visualize)
        visualizationActivationFunction(theNoiseFreeResponses, thePreActivationFunctionNoisyResponseInstances, theNeuralResponses, activationFunctionParams);
    end

end

function visualizationActivationFunction(theNoiseFreeResponses, thePreActivationFunctionNeuralResponses, theNeuralResponses, activationFunctionParams)
    if (size(theNoiseFreeResponses,1) == 1)
        theNoiseFreeResponses = repmat(theNoiseFreeResponses, [size(theNeuralResponses,1) 1]);
    end

    maxVisualizedResponse = prctile(abs(thePreActivationFunctionNeuralResponses(:)), 99.9);
    thePreActivationFunctionResponseRange = maxVisualizedResponse*[-1 1];
    % Quantize responses to 256 levels
    responseBins = linspace(thePreActivationFunctionResponseRange(1), thePreActivationFunctionResponseRange(2), 256);
    

    hFig = figure(33);
    clf;
    set(hFig, 'Position', [10 10 1500 600], 'Color', [1 1 1]);
    ax1 = subplot('Position', [0.05 0.1 0.42 0.85]);
    ax2 = subplot('Position', [0.55 0.1 0.42 0.85]);

    % Histogram of responses 
    [N1,edges] = histcounts(theNoiseFreeResponses(:), responseBins);
    [N2,edges] = histcounts(thePreActivationFunctionNeuralResponses(:), responseBins);
    N1 = N1/max(N1(:));
    N2 = N2/max(N2(:));
    b1 = bar(ax1, edges(1:end-1), N1(:), 1);
    hold(ax1, 'on');
    b2 = bar(ax1, edges(1:end-1), N2(:), 1);
    b1.FaceColor = [1 0 0];
    b1.EdgeColor = 'none';
    b2.FaceColor = [.6 .6 .6];
    b2.FaceAlpha = 0.5;
    b2.EdgeColor = 'none';
    hold(ax1, 'off');

    set(ax1, 'XLim', [thePreActivationFunctionResponseRange(1) thePreActivationFunctionResponseRange(2)], ...
            'YLim', [0 1.05], ...
            'YTick', [], ...
            'YColor', 'none');

    yyaxis(ax1, 'right')
    plot(ax1, thePreActivationFunctionNeuralResponses, theNeuralResponses, '.', 'Color', [0 0 0], 'MarkerSize', 16);
    xlabel(ax1,'input responses');
    ylabel(ax1,'output responses');
    set(ax1, 'XLim', [thePreActivationFunctionResponseRange(1) thePreActivationFunctionResponseRange(2)], ...
            'YLim', [thePreActivationFunctionResponseRange(1) thePreActivationFunctionResponseRange(2)], ...
            'YColor', [0 0 0], 'FontSize', 16);
    title(ax1, '\color[rgb]{1 0 0} noise-free responses \color[rgb]{0.5 0.5 0.5} noisy response instances');
    box(ax1, 'off');


    [N,edges] = histcounts(theNeuralResponses(:), responseBins);
    N = N/max(N(:));
    b = bar(ax2, edges(1:end-1), N(:), 1);
    b.FaceColor = [0.45 0.45 0.45];
    b.EdgeColor = 'none';

    set(ax2, 'XLim', [thePreActivationFunctionResponseRange(1) thePreActivationFunctionResponseRange(2)], ...
            'YLim', [0 1.05], ...
            'YTick', [], ...
            'FontSize', 16);
    xlabel(ax2,'output noisy response instances');
    box(ax2, 'off');
    drawnow;
end