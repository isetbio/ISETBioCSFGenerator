function dataOut = nreMidgetRGCMosaicSingleShot(...
    neuralEngineOBJ, neuralResponseParamsStruct, sceneSequence, ...
    sceneSequenceTemporalSupport, instancesNum, varargin)
% Compute function for computation of mRGCMosaic activations witout eye
% movements and a static stimulus (1 frame)
%
% Syntax:
%   dataOut = nreMidgetRGCMosaicSingleShot(...
%    neuralEngineOBJ, neuralResponseParamsStruct, sceneSequence, ...
%    sceneSequenceTemporalSupport, instancesNum, varargin);
%
% Description:
%    Function serving as the computeFunctionHandle for a @neuralResponseEngine
%    object. There are 2 ways to use this function.
%
%       [1] If called directly and with no arguments, 
%           dataOut = nreMidgetRGCMosaicSingleShot()
%       it does not compute anything and simply returns a struct with the 
%       defaultParams (optics and midget RGC mosaic params) that define the neural 
%       compute pipeline for this computation.
%
%       [2] If called from a parent @neuralResponseEngine object, 
%       it computes 'instancesNum' of midget RGC response sequences 
%       in response to the passed 'sceneSequence'.
%
%    It is not a good idea to try to call this function with arguments
%    directly - it should be called by the compute method of its parent
%    @neuralResponseEngine.
%
% Inputs:
%    neuralEngineOBJ                - the parent @neuralResponseEngine object that
%                                     is calling this function as its computeFunctionHandle
%    neuralResponseParamsStruct     - a struct containing properties of the
%                                     employed neural chain.
%    sceneSequence                  - a cell array of scenes defining the frames of a stimulus
%    sceneSequenceTemporalSupport   - the temporal support for the stimulus frames, in seconds
%    instancesNum                   - the number of response instances to compute
%
% Optional key/value input arguments:
%    'noiseFlags'                   - Cell array of strings containing labels
%                                     that encode the type of noise to be included
%                                     Valid values are: 
%                                        - 'none' (noise-free responses)
%                                        - 'random' (noisy response instances)
%                                     Default is {'random'}.
%   'rngSeed'                       - Integer.  Set rng seed. Empty (default) means don't touch the
%                                     seed.
%
% Outputs:
%    dataOut  - A struct that depends on the input arguments. 
%
%               If called directly with no input arguments, the returned struct contains
%               the defaultParams (optics and midget RGC mosaic params) that define the neural 
%               compute pipeline for this computation.  This can be useful
%               for a user interested in knowing what needs to be supplied
%               to this.
%
%             - If called from a parent @neuralResponseEngine), the returned
%               struct is organized as follows:
%                .neuralResponses : dictionary of responses indexed with 
%                                   labels corresponding to the entries of
%                                   the 'noiseFlags'  optional argument
%                .temporalSupport : the temporal support of the neural
%                                   responses, in seconds
%                .neuralPipeline  : a struct containing the optics and the mRGC mosaic
%                                   employed in the computation (only returned if 
%                                   the parent @neuralResponseEngine object has 
%                                   an empty neuralPipeline property)
%
%       The computed neural responses can be extracted as:
%           neuralResponses('one of the entries of noiseFlags') 
%       and are arranged in a matrix of:
%           [instancesNum x tTimeBins x mRGCs] 
%
% See Also:
%     t_neuralResponseCompute

% History:
%    05/04/23  NPC  Wrote it
%    11/03/24  FH   Edited it to incorporate metaContrast

% Examples:
%{
    % Usage case #1. Just return the default neural response params
    defaultParams = nreMidgetRGCMosaicSingleShot()

    % Usage case #2. Compute noise free, noisy, and repeatable noisy response instances
    % using a parent @neuralResponseEngine object and the default neural response params

    % Instantiate the parent @neuralResponseEngine object
    theNeuralEngineOBJ = neuralResponseEngine(@nreMidgetRGCMosaicSingleShot);

    % Retrieve the default params
    defaultNeuralEnginePipelineParams = nreMidgetRGCMosaicSingleShot();

    % Modify certain params of interest
    neuralResponsePipelineParams = defaultNeuralEnginePipelineParams;

    % Perhaps we want to use custom optics (not the optics that were used to optimize
    % the mRGCMosaic). We can pass the optics here.
    % neuralResponsePipelineParams.customOpticsToEmploy = oiCreate();

    % Perhaps we want to set the input cone mosaic integration time. 
    % Here, we set it to 200 msec
    neuralResponsePipelineParams.mRGCMosaicParams.coneIntegrationTimeSeconds = 200/1000;

    % Perhaps we want to manipulate noise in the input cone mosaic as well
    % as the mRGC mosaic. Here we set no (Poisson) noise for input cone mosaic
    % and we set random (not frozen) Gaussian noise in the mRGC responses
    % with a standard deviation of 400. The units of the Gaussian noise
    % must be appropriate to the units of the inputConeMosaic response. By
    % default, the inputConeMosaic response is excitations/integration time
    % So here, we set the noise sigma to 400 excitations/integration time
    neuralResponsePipelineParams.noiseParams.inputConeMosaicNoiseFlag = 'random';
    neuralResponsePipelineParams.noiseParams.mRGCMosaicNoiseFlag = 'random';
    neuralResponsePipelineParams.noiseParams.mRGCMosaicVMembraneGaussianNoiseSigma = 500;

    % Perhaps we want to compute mRGCMosaic responses not on the raw cone
    % mosaic excitation responses, but on their modulation with respect to
    % the cone excitation response to a null stimulus. This can be done as follows:
    % neuralResponsePipelineParams.theNullStimulusScene = nullStimulusScene;
    % See t_spatialCSFmRGCMosaic.m to see exactly how to do this

    theNeuralEngineOBJ = neuralResponseEngine(@nreMidgetRGCMosaicSingleShot, neuralResponsePipelineParams);

    % Instantiate a @sceneEngine object
    defaultGratingParams = sceGrating();

    % Modify some scene properties
    gratingParams = defaultGratingParams;

    % Make it an achromatic stimulus
    gratingParams.coneContrastModulation = 0.6*[1 1 1];

    % Make it 2 c/deg
    gratingParams.spatialFrequencyCyclesPerDeg = 2.0;

    % Make it a horizontally oriented grating
    gratingParams.orientationDegs = 0;

    % Make it 5x5 degs wide
    gratingParams.fovDegs = 5
    gratingParams.spatialEnvelopeRadiusDegs = gratingParams.fovDegs/6;

    % Make the grating a single-frame stimulus
    gratingParams.temporalModulationParams.stimDurationFramesNum = 1;
    gratingParams.temporalModulationParams.stimOnFrameIndices = 1;
    theSceneEngineOBJ = sceneEngine(@sceGrating, gratingParams);
    testContrast = 1.0;
    [theTestSceneSequence, theTestSceneTemporalSupportSeconds] = ...
        theSceneEngineOBJ.compute(testContrast);

    % Compute 16 response instances for a number of different noise flags
    instancesNum = 16;
    noiseFlags = {'random', 'none'};
    [theResponses, theResponseTemporalSupportSeconds] = theNeuralEngineOBJ.compute(...
            theTestSceneSequence, ...
            theTestSceneTemporalSupportSeconds, ...
            instancesNum, ...
            'noiseFlags', noiseFlags ...
            );

    % Retrieve the noise-free responses and the noisy response instances
    noiseFreeResponses = theResponses('none');
    noisyResponseInstances = theResponses('random');
    activationRange = [-1 1] * ...
        max([prctile(abs(noiseFreeResponses(:)),90) ...
             prctile(abs(noisyResponseInstances(:)),90)...
            ]);

    % Visualize the noise-free response
    instanceNo = 1;
    theNeuralEngineOBJ.neuralPipeline.mRGCMosaic.visualize(...
        'activation', noiseFreeResponses(instanceNo,:,:), ...
        'activationRange', activationRange, ...
        'verticalActivationColorBarInside', true, ...
        'plotTitle', 'noise-free response');

    % Visualize a noisy response instance
    instanceNo = 1;
    theNeuralEngineOBJ.neuralPipeline.mRGCMosaic.visualize(...
        'activation', noisyResponseInstances(instanceNo,:,:), ...
        'activationRange', activationRange, ...
        'verticalActivationColorBarInside', true, ...
        'plotTitle', 'noisy response instance');
%}

    % Check input arguments. If called with zero input arguments, just return the default params struct
    if (nargin == 0)
        dataOut = generateDefaultParams();
        return;
    end

    % Parse the input arguments
    p = inputParser;
    p.addParameter('noiseFlags', {'random'});
    p.addParameter('rngSeed',[], @(x) (isempty(x) || isnumeric(x)));
    p.addParameter('amputatescenes',true,@islogical)
    p.addParameter('justAddNoise',false,@islogical)
    p.addParameter('theBackgroundRetinalImage', struct('type', 'opticalimage'), @isstruct);
    varargin = ieParamFormat(varargin);
    p.parse(varargin{:});
    
    % Retrieve the response noiseFlag labels and validate them.
    noiseFlags = p.Results.noiseFlags;
    passedRngSeed = p.Results.rngSeed;
    neuralEngineOBJ.validateNoiseFlags(noiseFlags);

    % For each noise flag we generate a corresponing neural response, and all 
    % neural responses are stored in a dictionary indexed by the noiseFlag label.
    % Setup theNeuralResponses dictionary, loading empty responses for now
    theNeuralResponses = containers.Map();
    for idx = 1:length(noiseFlags)
        theNeuralResponses(noiseFlags{idx}) = [];
    end

    if (isempty(neuralEngineOBJ.neuralPipeline))
        % Generate neural pipeline here
        % Step1. Load a compute-ready mRGC mosaic given the passed params
        theMRGCmosaic = mRGCMosaic.loadComputeReadyRGCMosaic(...
            neuralResponseParamsStruct.mRGCMosaicParams, ...
            neuralResponseParamsStruct.opticsParams, ...
            neuralResponseParamsStruct.mRGCMosaicParams.retinalRFmodelParams);

        % Crop it to desired size
        theMRGCmosaic.cropToSizeAtEccentricity(...
            neuralResponseParamsStruct.mRGCMosaicParams.cropParams.sizeDegs, ...
            neuralResponseParamsStruct.mRGCMosaicParams.cropParams.eccentricityDegs);

        % Set input cone mosaic integration time
        theMRGCmosaic.inputConeMosaic.integrationTime = ...
            neuralResponseParamsStruct.mRGCMosaicParams.coneIntegrationTimeSeconds;

        % Set the mRGCMosaic noise flag
        theMRGCmosaic.noiseFlag = neuralResponseParamsStruct.noiseParams.mRGCMosaicNoiseFlag;

        % Set the mRGCMosaic noise
        theMRGCmosaic.vMembraneGaussianNoiseSigma = ...
            neuralResponseParamsStruct.noiseParams.mRGCMosaicVMembraneGaussianNoiseSigma;

        if (~isempty(neuralResponseParamsStruct.customOpticsToEmploy))
            % Employ the passed optics
            theOptics = neuralResponseParamsStruct.customOpticsToEmploy;
        else
            % Retrieve the optics that were used to optimize theMRGCmosaic
            if (~isempty(theMRGCmosaic.theNativeOptics))
                theOptics = theMRGCmosaic.theNativeOptics;
            elseif (~isempty(theMRGCmosaic.theCustomOptics))
                theOptics = theMRGCmosaic.theCustomOptics;
            else
                error('No optics found in the mRGCMosaic object!')
            end
        end

        returnTheNeuralPipeline = true;
    else
        % Load the optics from the previously computed neural pipeline
        theOptics = neuralEngineOBJ.neuralPipeline.optics;
        % Load the mRGC mosaic from the previously computed neural pipeline
        theMRGCmosaic = neuralEngineOBJ.neuralPipeline.mRGCMosaic;
        returnTheNeuralPipeline = false;
    end
    
    % Compute the sequence of optical images corresponding to the sequence 
    % of scenes if we are not just adding noise
    if (~p.Results.justAddNoise)
        % Get the number of scene sequences
        framesNum = numel(sceneSequence);
        theListOfOpticalImages = cell(1, framesNum);
        for frame = 1:framesNum
            theListOfOpticalImages{frame} = oiCompute(theOptics, sceneSequence{frame},'padvalue','mean');
        end
        % Generate an @oiSequence object containing the list of computed optical images
        theOIsequence = oiArbitrarySequence(theListOfOpticalImages, sceneSequenceTemporalSupport);
    end


    % Compute theConeMosaicNullResponse if the inputSignalType is set to
    % 'cone modulations' and if we have theNullStimulusScene
    if (strcmp(neuralResponseParamsStruct.mRGCMosaicParams.inputSignalType, 'cone_modulations')) && ...
       (~isempty(neuralResponseParamsStruct.theNullStimulusScene)) && (~p.Results.justAddNoise)
        
        fprintf('Operating on cone modulations\n');
        % Retrieve the null stimulus scene
        theNullStimulusScene = neuralResponseParamsStruct.theNullStimulusScene;

        % Compute the optical image of the null scene
        theNullSceneOI = oiCompute(theOptics, theNullStimulusScene,'padvalue','mean');

        % Compute theConeMosaicNullResponse, i.e., the input cone mosaic response to the NULL scene
        theConeMosaicNullResponse = theMRGCmosaic.inputConeMosaic.compute(...
            theNullSceneOI, ...
            'nTrials', 1);

        % Compute normalizing response (for computing the  modulated response)
        coneIndicesWithZeroNullResponse = find(theConeMosaicNullResponse == 0);
        coneMosaicNormalizingResponse = 1./theConeMosaicNullResponse;
        coneMosaicNormalizingResponse(coneIndicesWithZeroNullResponse) = 0;
        coneMosaicNormalizingResponse = reshape(coneMosaicNormalizingResponse, [1 1 numel(coneMosaicNormalizingResponse)]);
    elseif (~p.Results.justAddNoise)
        theConeMosaicNullResponse = [];
        fprintf('Operating on cone excitations\n');
    else
        theConeMosaicNullResponse = [];
    end

     % Set rng seed if one was passed. 
    if (~isempty(passedRngSeed))
        oldSeed = rng(passedRngSeed);
    end

    % Compute cone mosaic responses depending on theMRGCmosaic.inputConeMosaic.noiseFlag
    % Save the current coneMosaic noiseFlag so we can restore edit it later
    lastInputConeMosaicNoiseFlag = theMRGCmosaic.inputConeMosaic.noiseFlag;

    % Don't need to compute if we are just adding noise
    if (~p.Results.justAddNoise)
        switch (neuralResponseParamsStruct.noiseParams.inputConeMosaicNoiseFlag)
            case 'none'
                % Compute noise-free cone mosaic response instances
                fprintf('\tComputing noise-free cone mosaic responses\n');
                [theNoiseFreeConeMosaicResponses, ~, ~, ~, coneMosaicTemporalSupportSeconds] = ...
                        theMRGCmosaic.inputConeMosaic.compute(theOIsequence.frameAtIndex(1), ...
                            'nTrials', 1 ...
                );
    
                % Repmat so we have instancesNum identical copies of noise-free
                % cone mosaic responses
                theConeMosaicResponses = repmat(theNoiseFreeConeMosaicResponses, [instancesNum 1 1]);
    
            case 'frozen'
                % Compute input cone mosaic noisy response instances with a specified random noise seed for repeatability
                if (~isempty(passedRngSeed))
                    fprintf('\tComputing noisy cone mosaic responses with a frozen seed (%d)\n', passedRngSeed);
                    [~, theConeMosaicResponses, ~, ~, coneMosaicTemporalSupportSeconds] = ...
                            theMRGCmosaic.inputConeMosaic.compute(theOIsequence.frameAtIndex(1), ...
                                'nTrials', instancesNum, ...
                                'seed', passedRngSeed ...        % the passed random seed
                    );
                else
                   error('The input cone mosaic has a frozen noise specification, but no rngSeed was passed')
                end
    
            case 'random'
                % 1 in a million
                useSeed = randi(1e6,1,1);
                fprintf('\tComputing noisy cone mosaic responses with a random seed (%d)\n', useSeed);
        
                % Compute input cone mosaic noisy response instances with a  random noise seed
                [~, theConeMosaicResponses, ~, ~, coneMosaicTemporalSupportSeconds] = ...
                        theMRGCmosaic.inputConeMosaic.compute(theOIsequence.frameAtIndex(1), ...
                            'nTrials', instancesNum, ...
                            'seed', useSeed ...        % random seed
                    );
            otherwise
                error('Unknown inputConeMosaic noise flag: ''%s''.\n', theMRGCmosaic.inputConeMosaic.noiseFlag);
        end
    end

    % Restore the original cone mosaic noise flag
    theMRGCmosaic.inputConeMosaic.noiseFlag = lastInputConeMosaicNoiseFlag;

    % Transform the cone excitation responses to cone modulation responses
    if (~isempty(theConeMosaicNullResponse)) && (~p.Results.justAddNoise)
        % Transform the noise-free cone mosaic response modulation to a contrast response
        % i.e., relative to the cone mosaic response to the null (zero contrast) stimulus. 
        % This mimics the photocurrent response which is normalized with respect to the 
        % mean cone activation
        theConeMosaicResponses = ...
                    bsxfun(@times, bsxfun(@minus, theConeMosaicResponses, theConeMosaicNullResponse), ...
                    coneMosaicNormalizingResponse);
    end

    % Compute responses for each type of noise flag requested
    for idx = 1:length(noiseFlags)
        % Save the current the mRGC mosaic noiseFlag
        lastMRGCMosaicNoiseFlag = theMRGCmosaic.noiseFlag;
    
        % Set the mRGC mosaic noiseFlag to the tested value
        theMRGCmosaic.noiseFlag = ieParamFormat(noiseFlags{idx});

        % Compute the mRGC mosaic response
        switch (theMRGCmosaic.noiseFlag)
            case 'none'
                fprintf('\tComputing noise-free mRGC responses \n');
                if (p.Results.justAddNoise)
                    % Don't need to compute if we are just adding noise,
                    % because the noise free responses were passed in.           
                    theNeuralResponses(noiseFlags{idx}) = sceneSequence;
                    temporalSupportSeconds = sceneSequenceTemporalSupport;
                else
                    [theNeuralResponses(noiseFlags{idx}), ~, temporalSupportSeconds] = theMRGCmosaic.compute( ...
                        theConeMosaicResponses, coneMosaicTemporalSupportSeconds);
                end
            case 'frozen'
                if (~isempty(passedRngSeed))
                    fprintf('\tComputing noisy mRGC responses with a frozen seed (%d)\n', passedRngSeed);
                    if (p.Results.justAddNoise)
                        theNeuralResponses(noiseFlags{idx}) = theMRGCmosaic.noisyInstances(sceneSequence, 'noiseFlag', 'frozen', 'seed', passedRngSeed); 
                        temporalSupportSeconds = sceneSequenceTemporalSupport;
                    else
                        [~,theNeuralResponses(noiseFlags{idx}), temporalSupportSeconds] = theMRGCmosaic.compute( ...
                            theConeMosaicResponses, coneMosaicTemporalSupportSeconds, ...
                            'seed', passedRngSeed ...  % the passed random seed
                        );
                    end
                else
                     error('The mRGCmosaic has a frozen noise specification, but no rngSeed was passed')
                end

            case 'random'
                if (p.Results.justAddNoise)
                    theNeuralResponses(noiseFlags{idx}) = theMRGCmosaic.noisyInstances(sceneSequence,'noiseFlag','donotset');
                    temporalSupportSeconds = sceneSequenceTemporalSupport;
                else
                    % 1 in a million
                    useSeed = randi(1e6,1,1);
                    fprintf('\tComputing noisy mRGC responses with a random seed (%d)\n', useSeed);
                    % Compute noisy mRGC mosaic response instances
                    [~,theNeuralResponses(noiseFlags{idx}), temporalSupportSeconds] = theMRGCmosaic.compute( ...
                        theConeMosaicResponses, coneMosaicTemporalSupportSeconds, ...
                        'seed', useSeed ...
                     );
                end
            otherwise
                error('Unknown mRGCMosaic noise flag: ''%s''.\n', mRGCMosaicNoiseFlagToSet);
        end

        % Restore the mRGC mosaic noiseFlag
        theMRGCmosaic.noiseFlag = lastMRGCMosaicNoiseFlag;
    end

    % Restore rng seed if we set it
    if (~isempty(passedRngSeed))
        rng(oldSeed);
    end
    
    % Assemble the dataOut struct
    dataOut = struct(...
        'neuralResponses', theNeuralResponses, ...
    	'temporalSupport', temporalSupportSeconds);
    
    if (returnTheNeuralPipeline)
        dataOut.neuralPipeline.optics = theOptics;
        dataOut.neuralPipeline.mRGCMosaic = theMRGCmosaic;
    end

end

function p = generateDefaultParams()
    % Default params for this compute function
    opticsParams = struct(...
            'ZernikeDataBase', 'Polans2015', ...
            'examinedSubjectRankOrder', 6, ...
            'pupilDiameterMM', 3.0, ...
            'analyzedEye', 'right eye', ...
            'refractiveErrorDiopters', 0.0, ...
            'positionDegs', [] ...
        );

    % Neurons of the pre-computed mRGCMosaic have spatial RFs that were 
    % optimized using a double exponential surround model with parameters around
    % those of the 4-th H1 neuron recorded by Packer&Dacey (2002):
    % "Receptive field structure of H1 horizontal cells in macaque monkey
    % retina", (2002) JoV, 2, 272-292
    retinalRFmodelParams = struct(...
        'conePoolingModel', 'arbitraryCenterConeWeights_doubleExpH1cellIndex4SurroundWeights' ...
    );

    % Cropping params. If sizeDegs is empty there is no cropping.
    % We can crop the mRGCmosaic at different eccentricities. 
    % Passing an empty value for eccentricityDegs will crop
    % the mosaic at its center.
    cropParams = struct(...
        'sizeDegs', [], ...
        'eccentricityDegs', []);

    mRGCMosaicParams = struct(...
        'eccDegs', [7 0], ...
        'sizeDegs',  [6 3], ...
        'rgcType', 'ONcenterMidgetRGC', ...
        'cropParams', cropParams, ...
        'retinalRFmodelParams', retinalRFmodelParams, ...
        'inputSignalType', 'cone_modulations', ...
        'coneIntegrationTimeSeconds', 10/1000);

    % If opticsToEmploy is [], we use the optics that were used to optimize
    % the mRGCMosaic
    customOpticsToEmploy = [];

    % Noise params
    noiseParams = struct(...
        'inputConeMosaicNoiseFlag', 'random', ...
        'mRGCMosaicNoiseFlag', 'random', ...
        'mRGCMosaicVMembraneGaussianNoiseSigma', 0.15 ...
        );

    % If the user sets the nullStimulusScene, then mRGC responses are
    % computed on modulated cone responses with respect to the cone
    % response to the null stimulus scene. This is computed ...
    theNullStimulusScene = [];


    % Assemble all params in a struct
    p = struct(...
        'opticsParams', opticsParams, ...
        'mRGCMosaicParams', mRGCMosaicParams, ...
        'noiseParams', noiseParams, ...
        'customOpticsToEmploy', customOpticsToEmploy, ...
        'theNullStimulusScene', theNullStimulusScene ...
    );
end