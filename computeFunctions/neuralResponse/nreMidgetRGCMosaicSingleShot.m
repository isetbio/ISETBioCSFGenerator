function dataOut = nreMidgetRGCMosaicSingleShot(...
    neuralEngineOBJ, neuralResponseParamsStruct, sceneSequence, ...
    sceneSequenceTemporalSupport, instancesNum, varargin)
% Compute function for computation of midgetRGCMosaic activations witout eye movements 
%

    % Check input arguments. If called with zero input arguments, just return the default params struct
    if (nargin == 0)
        dataOut = generateDefaultParams();
        return;
    end

    % Parse the input arguments
    p = inputParser;
    p.addParameter('noiseFlags', {'random'});
    p.addParameter('rngSeed',[], @(x) (isempty(x) || isnumeric(x)));
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


    % Compute the sequence of optical images corresponding to the sequence of scenes
    framesNum = numel(sceneSequence);
    theListOfOpticalImages = cell(1, framesNum);
    for frame = 1:framesNum
        theListOfOpticalImages{frame} = oiCompute(sceneSequence{frame}, theOptics);
    end

    % Generate an @oiSequence object containing the list of computed optical images
    theOIsequence = oiArbitrarySequence(theListOfOpticalImages, sceneSequenceTemporalSupport);

    % Retrieve the null stimulus scene
    theNullStimulusScene = neuralResponseParamsStruct.theNullStimulusScene;

    % If we have a null stimulus scene, compute theConeMosaicNullResponse 
    if (~isempty(theNullStimulusScene))
        % Compute the optical image of the null scene
        theNullSceneOI = oiCompute(theNullStimulusScene, theOptics);

        % Compute theConeMosaicNullResponse, i.e., the input cone mosaic response to the NULL scene
        theConeMosaicNullResponse = theMRGCmosaic.inputConeMosaic.compute(...
            theNullSceneOI, ...
            'nTrials', 1);

        % Compute normalizing response (for computing the  modulated response)
        coneIndicesWithZeroNullResponse = find(theConeMosaicNullResponse == 0);
        coneMosaicNormalizingResponse = 1./theConeMosaicNullResponse;
        coneMosaicNormalizingResponse(coneIndicesWithZeroNullResponse) = 0;
        coneMosaicNormalizingResponse = reshape(coneMosaicNormalizingResponse, [1 1 numel(coneMosaicNormalizingResponse)]);
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

    switch (neuralResponseParamsStruct.noiseParams.inputConeMosaicNoiseFlag)
        case 'none'
            fprintf('\tComputing noise-free cone mosaic responses\n');
            % Compute noise-free cone mosaic response instances
            [theNoiseFreeConeMosaicResponses, ~, ~, ~, coneMosaicRemporalSupportSeconds] = ...
                    theMRGCmosaic.inputConeMosaic.compute(theOIsequence.frameAtIndex(1), ...
                        'nTrials', 1 ...
            );

            % Repmat so we have instancesNum identical copies of noise-free
            % cone mosaic responses
            theConeMosaicResponses = repmat(theNoiseFreeConeMosaicResponses, [instancesNum 1 1]);

        case 'frozen'
            if (~isempty(passedRngSeed))
                fprintf('\tComputing noisy cone mosaic responses with a frozen seed (%d)\n', passedRngSeed);
                % Compute input cone mosaic noisy response instances with a specified random noise seed for repeatability
                [~, theConeMosaicResponses, ~, ~, coneMosaicRemporalSupportSeconds] = ...
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
            [~, theConeMosaicResponses, ~, ~, coneMosaicRemporalSupportSeconds] = ...
                    theMRGCmosaic.inputConeMosaic.compute(theOIsequence.frameAtIndex(1), ...
                        'nTrials', instancesNum, ...
                        'seed', useSeed ...        % random seed
                );
        otherwise
            error('Unknown inputConeMosaic noise flag: ''%s''.\n', theMRGCmosaic.inputConeMosaic.noiseFlag);
    end

    % Restore the original cone mosaic noise flag
    theMRGCmosaic.inputConeMosaic.noiseFlag = lastInputConeMosaicNoiseFlag;

    % Transform the cone excitation responses to cone modulation responses
    if (~isempty(theConeMosaicNullResponse))
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

        mRGCMosaicNoiseFlagToSet = ieParamFormat(noiseFlags{idx});

        if (~strcmp(theMRGCmosaic.noiseFlag,  mRGCMosaicNoiseFlagToSet))
            fprintf(2, 'Overriding mRGCMosaic noiseFlag from ''%s'' to ''%s''.\n', mRGCMosaicNoiseFlagToSet, theMRGCmosaic.noiseFlag)
            mRGCMosaicNoiseFlagToSet = theMRGCmosaic.noiseFlag;
        end

        % Save the current the mRGC mosaic noiseFlag
        lastMRGCMosaicNoiseFlag = theMRGCmosaic.noiseFlag;
    
        % Set the mRGC mosaic noiseFlag to the tested value
        theMRGCmosaic.noiseFlag = mRGCMosaicNoiseFlagToSet;

        % Compute the mRGC mosaic response
        switch (mRGCMosaicNoiseFlagToSet)
            case 'none'
                fprintf('\tComputing noise-free mRGC responses \n');
                [theNeuralResponses(noiseFlags{idx}), ~, temporalSupportSeconds] = theMRGCmosaic.compute( ...
                    theConeMosaicResponses, coneMosaicRemporalSupportSeconds);
            case 'frozen'
                if (~isempty(passedRngSeed))
                    fprintf('\tComputing noisy mRGC responses with a frozen seed (%d)\n', passedRngSeed);
                    [~,theNeuralResponses(noiseFlags{idx}), temporalSupportSeconds] = theMRGCmosaic.compute( ...
                        theConeMosaicResponses, coneMosaicRemporalSupportSeconds, ...
                        'seed', passedRngSeed ...  % the passed random seed
                    );
                else
                     error('The mRGCmosaic has a frozen noise specification, but no rngSeed was passed')
                end

            case 'random'
                % 1 in a million
                useSeed = randi(1e6,1,1);
                fprintf('\tComputing noisy mRGC responses with a random seed (%d)\n', useSeed);
                % Compute noisy mRGC mosaic response instances
                [~,theNeuralResponses(noiseFlags{idx}), temporalSupportSeconds] = theMRGCmosaic.compute( ...
                    theConeMosaicResponses, coneMosaicRemporalSupportSeconds, ...
                    'seed', useSeed ...
                 );

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