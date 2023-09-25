function t_neuralResponseComputeCustomizedPipeline
% Use the @neuralResponseEngine object with an externally generated pipeline
%
% Syntax:
%    t_neuralResponseComputeCustomizedPipeline
%
% Description:
%    Demonstrates how to use a neuralEngine with an externally-generate neural pipeline.
%
% Inputs:
%    None.
%
% Outputs:
%   None.
%
% Optional key/value pairs:
%    None.
%
% See also: t_thresholdEngine, t_sceneGeneration, t_responseClassifier
%

% History:
%    03/30/2021  NPC  Wrote it.

    % Close figures
    close all;

    % Instantiate the scene engine with a compute function.
    % This is a function that the USER has to supply.
    sceneComputeFunction = @sceGrating;
    sceneParams = sceneComputeFunction();
    sceneParams.fovDegs = 3;  
    sceneParams.spatialFrequencyCyclesPerDeg = 10;
    sceneParams.spatialEnvelopeRadiusDegs = sceneParams.fovDegs/2;
    sceneParams.coneContrastModulation = [0.75 0.75 0.75];
    sceneParams.temporalModulationParams.stimOnFrameIndices = [1];
    sceneParams.temporalModulationParams.stimDurationFramesNum = 1;
    sceneParams.frameDurationSeconds = 20/1000;
    
    
    theSceneEngine = sceneEngine(sceneComputeFunction, sceneParams);
    

    % Generate the components of the neural pipeline
    % (1) a @cMosaic object
    theCMosaic = cMosaic('sizeDegs', [0.5 0.5], 'eccentricityDegs', [1 0]);
    
    % (2) optics corresponding to the mosaic's eccentricity
    oiEnsemble = theCMosaic.oiEnsembleGenerate(theCMosaic.eccentricityDegs, ...
                'zernikeDataBase', 'Polans2015', ...
                'subjectID', 10, ...
                'pupilDiameterMM', 3);
    theOptics = oiEnsemble{1};
            
    % Instantiate a neural response engine using a desired neural compute
    % function, here nrePhotopigmentExcitationsCmosaicSingleShot().
    theNeuralEngine = neuralResponseEngine(@nrePhotopigmentExcitationsCmosaicSingleShot);
    
    % Install an  an externally-supplied pipeline (here the @cMosaic and
    % the optics object we generated above.
    theNeuralEngine.customNeuralPipeline(struct(...
                    'coneMosaic', theCMosaic, ...
                    'optics', theOptics));
                
                
    % Specify a pedestal luminance with 70% contrast
    testContrast = 0.7;

    % Compute the scene sequence
    [theSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(testContrast);

    % Compute instances of neural responses to the input scene sequence.
    %
    % We specify the types of noise to be applied to the computed responses.
    % This has to be a cell array.  All nre compute functions need to understand 
    % 'none' and 'random'. Specific nre compute functions are allowed understand additional options
    % as appropriate to the specific model.
    %
    % It is possible to freeze the noise by specifying a seed for the
    % randome number generator through the 'rngSeed' key/value pair.  The
    % compute function should restore the rng to its current state if this
    % is passed, but not otherwise.
    instancesNum = 8;
    noiseFlags = {'random'};
    [theResponses, theResponseTemporalSupportSeconds] = theNeuralEngine.compute(...
            theSceneSequence, ...
            theSceneTemporalSupportSeconds, ...
            instancesNum, ...
            'noiseFlags', noiseFlags, ...
            'rngSeed', [] ...
            );

    noisyInstances = theResponses('random');
    assert(size(noisyInstances,1) == instancesNum);
    assert(size(noisyInstances,2) == length(theResponseTemporalSupportSeconds));
        
    % Visualize all responses
    renderNeuralResponse( theCMosaic, theResponses('random'));

end

function renderNeuralResponse(theCMosaic, theResponses)
    hFig = figure(); clf;
    meanResponse = squeeze(mean(theResponses,1));
    ax = subplot(1,2,1);
    theCMosaic.visualize('figureHandle', hFig, 'axesHandle', ax, ...
        'activation', squeeze(theResponses(1,:,:)), 'visualizedConeAperture', 'geometricArea');
    title('1st response instance');
    
    ax = subplot(1,2,2);
    theCMosaic.visualize('figureHandle', hFig, 'axesHandle', ax, ...
        'activation', meanResponse, 'visualizedConeAperture', 'geometricArea');
    title('mean response (%s)');
    colormap(gray);
end
