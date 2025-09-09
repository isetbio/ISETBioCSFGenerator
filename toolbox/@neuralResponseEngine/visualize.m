function visualize(obj, neuralResponses, temporalSupportSeconds, varargin)
% Gateway for visualization of @neuralResponseEngine responses. 
%
% Syntax:
%    visualize(obj, neuralResponses, temporalSupportSeconds, varargin);
%
% Description:
%   Gateway for visualization of @neuralResponseEngine responses.
%
%   If the 'customVisualizationFunctionHandle' public property is set to
%   some custom visualization function, the custom visualization function
%   is called. If it is not, we defaut to the generic visualization of
%   responses as implemented here in the defaultVisualization() function
%
% The visualize() method is called both from the noiseFree and the noisyInstances compute methods
%
% Inputs:
%    neuralResponses                - a 3D matrix, [kTrials x mNeurons x tFrames], of neural responses
%    temporalSupportSeconds         - the temporal support for the stimulus frames, in seconds
%    
% Optional key/value input arguments:
%    'figureHandle'                 - either [] or a figure handle
%    'axesHandle'                   - either [] or an axes handle
%    'responseLabel'                - string, just a comment about the responses visualized
%    'clearAxesBeforeDrawing'       - boolean, whether to clear the axes before drawing
%
% Outputs:
%   None
%
% Example usage:
%
% noiseFreeResponseParams = nreNoiseFreeMRGCMosaic([],[],[],[]);
% noisyInstancesParams = nreNoisyInstancesGaussian;
%
% theMRGCMosaicNeuralEngine = neuralResponseEngine( ...
%    @nreNoiseFreeMidgetRGCMosaic, ...
%    @nreNoisyInstancesGaussian, ...
%    noiseFreeResponseParams, ...
%    noisyInstancesParams);
%
% theMRGCMosaicNeuralEngine.visualizeEachCompute = true;
% theMRGCMosaicNeuralEngine.maxVisualizedNoisyResponseInstances = 2;
% theMRGCMosaicNeuralEngine.visualize(someNeuralResponses,theTemporalSupportInSecondsOfTheNeuralResponses);
%
% For more usage examples, see t_spatialCSF
%

% History:
%    01/11/2025  NPC  Wrote it


    if (~isempty(obj.customVisualizationFunctionHandle))
        if (~isempty(obj.neuralPipeline))
            % User specified custom visualization function and we have a non-empty neural pipeline
            obj.customVisualizationFunctionHandle(obj.neuralPipeline, neuralResponses, temporalSupportSeconds, ...
                obj.maxVisualizedNoisyResponseInstances, obj.visualizationMetaData, varargin{:});
        end
    else
        % Just the default visualization
        defaultVisualization(neuralResponses, temporalSupportSeconds, ...
            obj.maxVisualizedNoisyResponseInstances, varargin{:});
    end
end

function defaultVisualization(neuralResponses, temporalSupportSeconds, maxVisualizedInstances, varargin)
   
    [figureHandle, axesHandle, clearAxesBeforeDrawing, responseLabel] = ...
        neuralResponseEngine.parseVisualizationOptions(varargin{:});

    %% Set figure size
    if (isempty(figureHandle))
        figureHandle = figure(); clf;
        set(figureHandle, 'Position', [10 10 700 700], 'Color', [1 1 1]);
    end
    if (isempty(axesHandle))
        axesHandle = subplot('Position', [0.09 0.07 0.85 0.88]);
    end
    
    % If we were passed more than one axes handles, just use the first one
    if (iscell(axesHandle))&&(numel(axesHandle) > 1)
        axesHandle = axesHandle{1};
    end

    if (clearAxesBeforeDrawing)
        cla(axesHandle);
    end
    

    [nInstances, mNeurons, tBins] = size(neuralResponses);
    activationRange = prctile(neuralResponses(:), [0 100]);

    % Treat special case of zero activation range
    if (activationRange(1) == activationRange(2))
        activationRange = activationRange(1) + [-10*eps 10*eps];
    end
    
    % Check for consistency in dimensions
    assert(tBins == size(temporalSupportSeconds,2), 'Inconsistency in time dimension');

    % Visualized instances num
    visualizedInstancesNum = min([nInstances maxVisualizedInstances]);

    if (numel(temporalSupportSeconds)>1)
        dt = temporalSupportSeconds(2)-temporalSupportSeconds(1);
        XLim = [temporalSupportSeconds(1)-dt/2 temporalSupportSeconds(end)+dt/2];
        XTick = temporalSupportSeconds(1):dt:temporalSupportSeconds(end);
    else
        XLim = [temporalSupportSeconds(1)-0.01 temporalSupportSeconds(1)+0.01];
        XTick = temporalSupportSeconds(1);
    end

    if (mNeurons == 1)
        plot(axesHandle, temporalSupportSeconds, neuralResponses(1:visualizedInstancesNum,1,:), 'k-', 'LineWidth', 1.5);
        xlabel(axesHandle, 'time (msec)');
        ylabel(axesHandle, 'response amplitude');
        title(axesHandle, responseLabel);
        set(axesHandle, 'XLim', XLim, 'YLim', activationRange, 'XTick', XTick);
        set(axesHandle, 'fontSize', 16);
        drawnow;
    elseif (mNeurons > 1)
        set(axesHandle, 'fontSize', 16);
        xlabel(axesHandle, 'time (msec)');
        ylabel(axesHandle, 'neuron index');
        
        for iTrial = 1:visualizedInstancesNum
            theNeuralResponseMap = squeeze(neuralResponses(iTrial,:,:));
            if (size(theNeuralResponseMap,2) == mNeurons)
                 theNeuralResponseMap =  theNeuralResponseMap';
            end

            imagesc(axesHandle, temporalSupportSeconds, 1:mNeurons, theNeuralResponseMap);
            axis(axesHandle, 'xy');
            set(axesHandle, ...
                'CLim', activationRange, ...
                'XLim', XLim, ...
                'XTick', XTick, ...
                'YLim', [1 mNeurons]);
            set(axesHandle, 'fontSize', 16);
            colorbar(axesHandle, 'NorthOutside')
            title(axesHandle, sprintf('%s (trial: %d)',responseLabel, iTrial));
            
            colormap(axesHandle, gray);
            drawnow;
        end
    end
end
