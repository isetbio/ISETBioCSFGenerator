 % Meta contrast computation of spatial CSF in different color directions
%
% Description:
%    Use ISETBioCSFGenerator to run out CSFs in different color directions.
%    This example uses an ideal Poisson observer and circularly
%    windowed gratings of constant size, and is set up to illustrate the
%    meta contrast approach.
%
%    The meta contrast approach works by taking advantage of the function that
%    for fixed spatial parameters, we can compute responses to any scene
%    by synthesis from the responses from just two precomputed contrasts.
%
%    Here's the idea, which also works for a scene sequence. Let's suppose
%    that r0 is the noise-free response to a background scene, and r1 is the
%    noise-free response to a scene with a contrast of c1.  We can rewrite 
%       r1 = r0 + c1*[(r1-r0)/c1] = r0 + c1*deltaR
%    Let deltaR = [(r1-r0)/c1].  If I then consider some other contrast c2
%    where I want the noise-free response, I can then compute it as
%       r2 = r0 + c2*deltaR
%    when the nre's noise-free responses are linear in the photons.  The
%    calculation above is much faster than recomputing r2 through the nre,
%    at least in realistic case.
%
%    In addition, adding the noise to the noise-free response is
%    comptuationally fast.
%
%    The name meta contrast here has nothing to do with the concept of meta
%    contrast masking.
%
% See also: t_spatialCSF, t_thresholdEngine, t_modulatedGratingsSceneGeneration,
%           t_chromaticThresholdContour, computeThreshold, computePerformance
%

% History:
%   08/11/24  dhb   Wrote it.

function thresholdRet = t_metaContastCSF(varargin)
    p = inputParser;
    p.addParameter('filter', struct('spectralSupport',[],'transmission',[]), @isstruct);
    p.addParameter('doValidationCheck', true, @islogical);
    parse(p, varargin{:});
    filter = p.Results.filter;
    doValidationCheck = p.Results.doValidationCheck;

    % Close out old figs
    close all;

    % Freeze rng for replicatbility
    rng(0);
    
    % List of spatial frequencies to be tested.
    spatialFreqs = [4, 8, 16, 32];
    if (length(spatialFreqs) ~= 4 | ~all(spatialFreqs == [4, 8, 16, 32]))
        doValidationCheck = false;
    end
    
    % Choose stimulus chromatic direction specified as a 1-by-3 vector
    % of L, M, S cone contrast.  These vectors get normalized below, so only
    % their direction matters in the specification.
    stimType = 'luminance';
    switch (stimType)
        case 'luminance'
            chromaDir = [1.0, 1.0, 1.0]';
        case 'red-green'
            chromaDir = [1.0, -1.0, 0.0]';
            doValidationCheck = false;
        case 'L-isolating'
            chromaDir = [1.0, 0.0, 0.0]';
            doValidationCheck = false;
    end
    
    % Set the RMS cone contrast of the stimulus. Things may go badly if you
    % exceed the gamut of the monitor, so we are conservative and set this at a
    % value that is within gamut of typical monitors and don't worry about it
    % further for this tutorial.  A vector length contrast of 0.1 should be OK.
    rmsContrast = 0.1;
    chromaDir = chromaDir / norm(chromaDir) * rmsContrast;
    assert(abs(norm(chromaDir) - rmsContrast) <= 1e-10);
    
    %% Create actaul neural response engine
    %
    % This calculations isomerizations in a patch of cone mosaic with Poisson
    % noise, and includes optical blur.
    neuralParams = nrePhotopigmentExcitationsCmosaic;
    neuralParams.coneMosaicParams.sizeDegs = [0.5 0.5]; 
    neuralParams.coneMosaicParams.timeIntegrationSeconds = 0.1;
    theNeuralEngine = neuralResponseEngine(@nrePhotopigmentExcitationsCmosaic, neuralParams);
    if (~all(neuralParams.coneMosaicParams.sizeDegs == [0.5 0.5]))
        doValidationCheck = false;
    end
    if (neuralParams.coneMosaicParams.timeIntegrationSeconds ~= 0.1)
        doValidationCheck = false;
    end
    
    %% Instantiate the Poisson responseClassifierEngine
    %
    % rcePoisson makes decision by performing the Poisson likelihood ratio test
    % Also set up parameters associated with use of this classifier.
    classifierEngine = responseClassifierEngine(@rcePoisson);
    classifierPara = struct('trainFlag', 'none', ...
                            'testFlag', 'random', ...
                            'nTrain', 1, 'nTest', 128);
    
    %% Parameters for threshold estimation/quest engine
    % The actual threshold varies enough with the different engines that we
    % need to adjust the contrast range that Quest+ searches over, as well as
    % the range of psychometric function slopes. Threshold limits are computed
    % as 10^-logThreshLimitVal.  The reason it is log units is that below we
    % define the PF for the questEngine as @qpPFWeibullLog. Had we used the
    % default (@qpPFWeibull), the units would have been dB.
    %
    % Also note explicit passing of proportion correct criterion for threshold.
    % The default value of 0.81606 matches the parameterization of mQUESTPlus'
    % Weibull PFs, when lapse rate is 0 and guess rate is 0.5.  But it seems
    % better to pass it explicitly so we know what it is. Keeping 0.81606 for
    % backward compatibilty.
    %
    % There are two separate structures below. The conceptual distinction
    % between them is not entirely clear.  These are interpretted by
    % computeThreshold.
    %
    % See t_thresholdEngine.m for more on options of the two different mode of
    % operation (fixed numer of trials vs. adaptive)
    thresholdPara = struct('logThreshLimitLow', 2.4, ...
                           'logThreshLimitHigh', 0.0, ...
                           'logThreshLimitDelta', 0.02, ...
                           'slopeRangeLow', 1/20, ...
                           'slopeRangeHigh', 50/20, ...
                           'slopeDelta', 2.5/20, ...
                           'thresholdCriterion', 0.81606);
    
    questEnginePara = struct( ...
        'qpPF',@qpPFWeibullLog, ...
        'minTrial', 1280, ...
        'maxTrial', 1280, ...
        'numEstimator', 1, ...
        'stopCriterion', 0.05);
    
    %% Compute threshold for each spatial frequency
    % See toolbox/helpers for functions createGratingScene computeThreshold
    dataFig = figure();
    logThreshold = zeros(1, length(spatialFreqs));
    for idx = 1:length(spatialFreqs)
        % Create a grating scene engine with a particular chromatic direction,
        % spatial frequency, and temporal duration.  Put grating in sine phase
        % becuase that keeps the spatial mean constant across spatial
        % frequencies.
        %
        % Create scene produces square scenes.  We use the min of the mosaic
        % field size to pick a reasonable size
        gratingSceneEngine = createGratingScene(chromaDir, spatialFreqs(idx),...
            'fovDegs', min(neuralParams.coneMosaicParams.sizeDegs), ...
            'duration', neuralParams.coneMosaicParams.timeIntegrationSeconds, ...
            'spatialPhase', 90, ...
            'filter', filter...
            );

        % Create the sceMetaContrast scene engine
        metaSceneEngineParams = sceMetaContrast;
        theMetaSceneEngine = sceneEngine(@sceMetaContrast,metaSceneEngineParams);
        
        % Create nreMetaContrast using the actual scene and neural engines
        metaNeuralResponseEngineParams = nreMetaContrast;
        metaNeuralResponseEngineParams.contrast0 = 0;
        metaNeuralResponseEngineParams.contrast1 = 1;
        metaNeuralResponseEngineParams.sceneEngine = gratingSceneEngine;
        metaNeuralResponseEngineParams.neuralEngine = theNeuralEngine;
        metaNeuralResponseEngineParams.noiseAddMethod = 'nre';
        theMetaNeuralEngine = neuralResponseEngine(@nreMetaContrast,metaNeuralResponseEngineParams);

        % Compute the threshold for our grating scene with meta scene and
        % and neural response engines. This function does a lot of
        % work, see t_tresholdEngine and the function itself, as well as
        % function computePerformance.
        [logThreshold(idx), questObj, ~, para(idx,:)] = ...
            computeThreshold(theMetaSceneEngine, theMetaNeuralEngine, classifierEngine, ...
            classifierPara, thresholdPara, questEnginePara, 'TAFC',true);
        
        % Plot stimulus
        figure(dataFig);
        subplot(length(spatialFreqs), 2, idx * 2 - 1);
        
        visualizationContrast = 1.0;
        [theSceneSequence] = gratingSceneEngine.compute(visualizationContrast);
        gratingSceneEngine.visualizeStaticFrame(theSceneSequence);
        
        % Plot data and psychometric curve 
        % with a marker size of 2.5
        subplot(length(spatialFreqs), 2, idx * 2);
        questObj.plotMLE(2.5,'para',para(idx,:));
    end
    set(dataFig, 'Position',  [0, 0, 800, 800]);
    
    % Convert returned log threshold to linear threshold
    threshold = 10 .^ logThreshold;
    
    %% Plot Contrast Sensitivity Function
    theCsfFig = figure();
    loglog(spatialFreqs, 1 ./ threshold, '-ok', 'LineWidth', 2);
    xticks(spatialFreqs); xlim([spatialFreqs(1), spatialFreqs(end)]);
    yticks([2,5,10,20,50]); ylim([1, 50]);
    xlabel('Spatial Frequency (cyc/deg)');
    ylabel('Sensitivity');
    set(theCsfFig, 'Position',  [800, 0, 600, 800]);
    
    %% Do a check on the answer
    % 
    % So that if we break something in the future we will have
    % a chance of knowing it.
    if (doValidationCheck)
        validationThresholds = [0.0351    0.0827    0.1534    0.5529];
        if (any(abs(threshold-validationThresholds)./validationThresholds > 0.25))
            error('Do not replicate validation thresholds to 25%. Check that parameters match, or for a bug.');
        end
    end

    %% Return a value if it was requested
    if (nargout > 0)
        thresholdRet = threshold;
    end
end
