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
% See also: t_spatialCSFcMosaic, t_sceneGeneration

% History:
%    03/30/2021  NPC  Wrote it.

    % Close figures
    close all;

    % Instantiate the scene engine with a compute function.
    %
    % This is a function that the user has to supply, but we have a number
    % of various useful ones written already.
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
    %
    % A cMosaic object
    theCMosaic = cMosaic('sizeDegs', [0.5 0.5], 'eccentricityDegs', [1 0], ...
        'integrationTime', sceneParams.frameDurationSeconds);
    
    % Optics corresponding to the mosaic's eccentricity
    oiEnsemble = theCMosaic.oiEnsembleGenerate(theCMosaic.eccentricityDegs, ...
                'zernikeDataBase', 'Polans2015', ...
                'subjectID', 10, ...
                'pupilDiameterMM', 3);
    theOptics = oiEnsemble{1};
            
    % Instantiate a neural response engine using desired neural engine 
    % noiseFree and noisyInstances compute functions.  Here we illustrate
    % a cMosaic-based nre.
    theNeuralEngine = neuralResponseEngine(@nreNoiseFreePhotopigmentExcitationsCMosaic, ...
        @nreNoisyInstancesPoisson);
    
    % Install an  an externally-supplied pipeline (here the @cMosaic and
    % the optics object we generated above.  Note that this must be matched
    % up with the particular compute functions that the neural engine has
    % been initialized with - different compute functions take different
    % pipeline parameters.
    noiseFreeCustomPipeline = struct(...
                    'coneMosaic', theCMosaic, ...
                    'optics', theOptics);
    noisyInstancesCustomPipeline = [];
    theNeuralEngine.customNeuralPipeline(noiseFreeCustomPipeline,noisyInstancesCustomPipeline);
                        
    % Specify a pedestal luminance with 70% contrast
    testContrast = 0.7;

    % Compute the scene sequence
    [sceneSequence, sceneTemporalSupportSeconds] = theSceneEngine.compute(testContrast);

    % Compute instances of neural responses to the input scene sequence.
    %
    % First the noise free response
    instancesNum = 8;
    [noiseFreeResponse, temporalSupportSeconds] = theNeuralEngine.computeNoiseFree( ...
        sceneSequence,sceneTemporalSupportSeconds);

    % And then the noisy instances
    noisyInstances = theNeuralEngine.computeNoisyInstances(noiseFreeResponse,temporalSupportSeconds, ...
        instancesNum,'random');
    assert(size(noisyInstances,1) == instancesNum);
    assert(size(noisyInstances,3) == length(temporalSupportSeconds));
        
    % Visualize all responses
    renderNeuralResponse( theCMosaic, noisyInstances);
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
