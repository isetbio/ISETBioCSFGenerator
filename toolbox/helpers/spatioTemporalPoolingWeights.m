function [poolingWeights,quadraturePoolingWeights, ...
    noiseFreeSpatioTemporalResponseNullStimulus, ...
    noiseFreeSpatioTemporalResponseTestStimulus] = spatioTemporalPoolingWeights(theNeuralEngine,theScene,theNullScene, temporalSupport)
% Compute spatiotemporal pooling weights for a given scene 
%
% Syntax:
%   [poolingWeights,quadraturePoolingWeights, ...
%    noiseFreeSpatioTemporalResponseNullStimulus, ...
%    noiseFreeSpatioTemporalResponseTestStimulus] = spatioTemporalPoolingWeights(theNeuralEngine,theScene,theNullScene, temporalSupport)
%
% Description:
%    Given a neural engine, a test scene (or a scene sequence) and the 
%    corresponding null scene (zero contrast), compute spatial or spatiotemporal 
%    weighting kernels (linear and quadrature) based on the noise-free response of the 
%    neural engine to the test scene and the null scene
%
% Inputs:
%   theNeuralEngine      - The neuralEngine to be employed
%   theScene             - The test scene (or scene sequence) for which to
%                          compute matching linear and quadrature spatial pooling kernels
%   theNullScene         - The zero contrast scene
%   temporalSupport      _ Temporal support for the scene sequence
%
%
% Outputs:
%   poolingWeights                               - The linear spatiotemporal pooling weights
%   quadraturePoolingWeights                     - The quadrature spatiotemporal pooling weights
%   noiseFreeSpatioTemporalResponseNullStimulus  - The noise-free spatiotemporal response to the null scene
%   noiseFreeSpatioTemporalResponseTestStimulus  - The noise-free spatiotemporal response to the test scene
%
% Optional key/value pairs:
%   none
%
% See also: t_modulatedGratingsSceneGenerateion, t_spatialCSF
%

% History:
%   8/28/21  npc  Wrote it.

    % Compute the noise-free response to the NULL stimulus
    fprintf('Computing noise-free response to the NULL stimulus\n');
    responseInstancesNum = 1;
    response = theNeuralEngine.compute(...
        theNullScene, ...
        temporalSupport, ...
        responseInstancesNum, ...
        'noiseFlags', {'none'});
    noiseFreeSpatioTemporalResponseNullStimulus = response('none');
    
    fprintf('Computing noise-free response to the TEST stimulus\n');
    % Compute the noise-free response to the TEST stimulus
    response = theNeuralEngine.compute(...
        theScene, ...
        temporalSupport, ...
        responseInstancesNum, ...
        'noiseFlags', {'none'});
    noiseFreeSpatioTemporalResponseTestStimulus = response('none');
    
    
    % Only keep the last time bin - this is a spatial kernel only, not a
    % spatiotemporal kernel
    noiseFreeSpatioTemporalResponseTestStimulus = noiseFreeSpatioTemporalResponseTestStimulus(1,end,:);
    noiseFreeSpatioTemporalResponseNullStimulus = noiseFreeSpatioTemporalResponseNullStimulus(1,end,:);
    
    % Compute the differential response and use this to derive the linear
    % pooling kernel
    diffResponse = noiseFreeSpatioTemporalResponseTestStimulus - noiseFreeSpatioTemporalResponseNullStimulus;
    diffResponse = diffResponse / max(abs(diffResponse(:)));
    
    
    % Generate image
    interpolationMethod = 'linear';
    extrapolationMethod = 'linear';
    rfPos = theNeuralEngine.neuralPipeline.coneMosaic.coneRFpositionsDegs;
    xq = rfPos(:,1);
    xq = xq-mean(xq(:));
    yq = rfPos(:,2);
    yq = yq-mean(yq(:));
    xx = linspace(min(xq(:)), max(xq(:)),256);
    yy = linspace(min(yq(:)), max(yq(:)),256);
    sigma = max(abs(xq(:)))/3;
    
    [X,Y] = meshgrid(xx,yy);
    [X2,Y2] = ndgrid(yy,xx);
    envelope = exp(-0.5*(X/sigma).^2) .* exp(-0.5*(Y/sigma).^2);
    
    poolingWeights = 0*diffResponse;
    quadraturePoolingWeights = 0*diffResponse;
    
    displayPoolingWeights = true;
    cMap = brewermap(1000, '*RdBu');
    
    for tBin = 1:size(poolingWeights,2)
        dR = diffResponse(1,1,:);
        F = scatteredInterpolant(xq(:),yq(:), dR(:), ...
            interpolationMethod, extrapolationMethod);
        poolingWeights2DMap = F(X,Y);
        
         
        quadraturePoolingWeights2DMap = (imag(hilbert(poolingWeights2DMap')))';
        %quadraturePoolingWeights2DMap = imag(hilbert(poolingWeights2DMap));
    
        % Gaussian envelope
        poolingWeights2DMap = envelope .* poolingWeights2DMap;
        quadraturePoolingWeights2DMap = envelope .* quadraturePoolingWeights2DMap;
        

        F = griddedInterpolant(X2,Y2,quadraturePoolingWeights2DMap);
        qW = F(yq(:), xq(:));
        qW = qW / sum(abs(qW(:)));
        quadraturePoolingWeights(1,tBin,:) = reshape(qW, [1 numel(qW)]) / max(abs(qW(:)));

        F = griddedInterpolant(X2,Y2,poolingWeights2DMap);
        pW = F(yq(:), xq(:));
        pW = pW / sum(abs(pW(:)));
        poolingWeights(1,tBin,:) = reshape(pW, [1 numel(pW)]) / max(abs(pW(:)));

        
        if (displayPoolingWeights)
            
%             hFig = figure(11);
%             ax = subplot(1,2,1);
%             imagesc(ax,poolingWeights2DMap, [-1 1]);
%             axis(ax, 'image');
%             title(sprintf('direct, tBin = %d', tBin));
%             
%             ax = subplot(1,2,2);
%             imagesc(ax,quadraturePoolingWeights2DMap,  [-1 1]);
%             axis(ax,'image');
%             title(sprintf('quadrature, tBin = %d', tBin));
            
            hFig = figure(12);
            activationRange = [-1 1] * max([max(abs(squeeze(poolingWeights(1,tBin,:)))) ...
                                            max(abs(squeeze(quadraturePoolingWeights(1,tBin,:))))]);
            ax = subplot(1,2,1);
            theNeuralEngine.neuralPipeline.coneMosaic.visualize(...
                'figureHandle', hFig, 'axesHandle', ax, ...
                'activation', poolingWeights(1,tBin,:), ...
                'activationRange', activationRange, ...
                'crossHairsOnMosaicCenter', true, ...
                'activationColorMap', cMap, ...
                'plotTitle', 'direct');
            
            hFig = figure(12);
            ax = subplot(1,2,2);
            theNeuralEngine.neuralPipeline.coneMosaic.visualize(...
                'figureHandle', hFig, 'axesHandle', ax, ...
                'activation', quadraturePoolingWeights(1,tBin,:), ...
                'activationRange', activationRange, ...
                'crossHairsOnMosaicCenter', true, ...
                'activationColorMap', cMap, ...
                'plotTitle', 'quadrature');
            
            drawnow;
            pause(0.1);
        end
    end
    

end
