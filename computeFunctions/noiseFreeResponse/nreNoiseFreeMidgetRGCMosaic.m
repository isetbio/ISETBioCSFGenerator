function dataOut = nreNoiseFreeMidgetRGCMosaic( ...
    neuralEngine, noiseFreeComputeParams, sceneSequence, ...
    sceneSequenceTemporalSupport, varargin)
% Compute function for computation of mRGCMosaic activations
%
% Syntax:
%   dataOut = nreNoiseFreeMidgetRGCMosaic(...
%    neuralEngine, noiseFreeComputeParams, sceneSequence, ...
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
%    neuralEngine                   - the parent @neuralResponseEngine object that
%                                     is calling this function as its computeFunctionHandle
%    noiseFreeComputeParams         - a struct containing properties of the
%                                     employed neural chain.
%    sceneSequence                  - a cell array of scenes defining the frames of a stimulus
%    sceneSequenceTemporalSupport   - the temporal support for the stimulus frames, in seconds
%
% Outputs:
%    dataOut  - A struct that depends on the input arguments.
%
%               If called directly with no input arguments, the returned struct contains
%               the defaultParams, which here is the empty struct.
%
%             - If called from a parent @neuralResponseEngine), the returned
%               struct is organized as follows:
%
%              .neuralResponses : matrix of neural responses.  Each column
%                                 is the response vector for one frame of the input.
%              .temporalSupport : the temporal support of the neural
%                                   responses, in seconds
%              .noiseFreeResponsePipeline : a struct that the parent
%                                   @neuralResponseEngine
%                                   can tuck away.  Only returned if the
%                                   parent object has an empty parameters
%                                   struct. For this routine, it has
%                                   information about the oi and cMosaic.
%
% Optional key/value input arguments:
%   'verbose'               - Logical (default true). Print things out.
%
% See Also:
%   t_neuralResponseCompute

% History:
%    05/04/23  NPC  Wrote it
%    11/03/24  FH   Edited it to incorporate metaContrast
%    12/23/24  dhb  Rewrite for major architecture redo.

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
    noiseFreeResponsePipelineParams = defaultNeuralEnginePipelineParams;

    % Perhaps we want to use custom optics (not the optics that were used to optimize
    % the mRGCMosaic). We can pass the optics here.
    % noiseFreeResponsePipelineParams.customOpticsToEmploy = oiCreate();

    % Perhaps we want to set the input cone mosaic integration time. 
    % Here, we set it to 200 msec
    noiseFreeResponsePipelineParams.mRGCMosaicParams.coneIntegrationTimeSeconds = 200/1000;

    % Perhaps we want to manipulate noise in the input cone mosaic as well
    % as the mRGC mosaic. Here we set no (Poisson) noise for input cone mosaic
    % and we set random (not frozen) Gaussian noise in the mRGC responses
    % with a standard deviation of 400. The units of the Gaussian noise
    % must be appropriate to the units of the inputConeMosaic response. By
    % default, the inputConeMosaic response is excitations/integration time
    % So here, we set the noise sigma to 400 excitations/integration time
    noiseFreeResponsePipelineParams.noiseParams.inputConeMosaicNoiseFlag = 'random';
    noiseFreeResponsePipelineParams.noiseParams.mRGCMosaicNoiseFlag = 'random';
    noiseFreeResponsePipelineParams.noiseParams.mRGCMosaicVMembraneGaussianNoiseSigma = 500;

    % Perhaps we want to compute mRGCMosaic responses not on the raw cone
    % mosaic excitation responses, but on their modulation with respect to
    % the cone excitation response to a null stimulus. This can be done as follows:
    % noiseFreeResponsePipelineParams.theNullStimulusScene = nullStimulusScene;
    % See t_spatialCSFmRGCMosaic.m to see exactly how to do this

    theNeuralEngineOBJ = neuralResponseEngine(@nreMidgetRGCMosaicSingleShot, noiseFreeResponsePipelineParams);

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
%
% Allow for possibility that other nre's take key/value pairs that we
% can ignore, so set KeepUnmatched to true.
p = inputParser;
p.KeepUnmatched = true;
%p.addParameter('theBackgroundRetinalImage', struct('type', 'opticalimage'), @isstruct);
p.addParameter('verbose',true,@islogical);
varargin = ieParamFormat(varargin);
p.parse(varargin{:});
verbose = p.Results.verbose;

% Create/get key objects on first call
if (isempty(neuralEngine.neuralPipeline) | ~isfield(neuralEngine.neuralPipeline,'noiseFreeResponse'))
    % Step1. Load a compute-ready mRGC mosaic given the passed params.
    %
    % This also produces the cone mosaic on which the RGC mosaic is built.
    theMRGCmosaic = mRGCMosaic.loadComputeReadyRGCMosaic(...
        noiseFreeComputeParams.mRGCMosaicParams, ...
        noiseFreeComputeParams.opticsParams, ...
        noiseFreeComputeParams.mRGCMosaicParams.retinalRFmodelParams);

    % Crop it to desired size
    theMRGCmosaic.cropToSizeAtEccentricity(...
        noiseFreeComputeParams.mRGCMosaicParams.cropParams.sizeDegs, ...
        noiseFreeComputeParams.mRGCMosaicParams.cropParams.eccentricityDegs);

    % Set input cone mosaic integration time
    theMRGCmosaic.inputConeMosaic.integrationTime = ...
        noiseFreeComputeParams.mRGCMosaicParams.coneIntegrationTimeSeconds;

    % No noise added for this pipeline
    theMRGCmosaic.inputConeMosaic.noiseFlag = 'none';
    theMRGCmosaic.noiseFlag = 'none';

    % Set the optics if passed.  This can be used to overrride the default
    % optics that were used to compute the mRGC object.
    if (~isempty(noiseFreeComputeParams.customOpticsToEmploy))
        % Employ the passed optics
        theOptics = noiseFreeComputeParams.customOpticsToEmploy;
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

    % Compute theConeMosaicNullResponse if the inputSignalType is set to
    % 'cone modulations' and if we have theNullStimulusScene
    if ( (strcmp(noiseFreeComputeParams.mRGCMosaicParams.inputSignalType,'cone_modulations')) && ...
            (~isempty(noiseFreeComputeParams.theNullStimulusScene)) )

        if (verbose)
            fprintf('Operating on cone modulations\n');
        end

        % Retrieve the null stimulus scene
        nullStimulusScene = noiseFreeComputeParams.theNullStimulusScene;

        % Compute the optical image of the null scene
        nullSceneOI = oiCompute(theOptics, nullStimulusScene,'padvalue','mean');

        % Compute theConeMosaicNullResponse, i.e., the input cone mosaic response to the NULL scene
        coneMosaicNullResponse = theMRGCmosaic.inputConeMosaic.compute(...
            nullSceneOI, ...
            'nTrials', 1);

        % Compute normalizing response (for computing the  modulated response)
        coneIndicesWithZeroNullResponse = find(coneMosaicNullResponse == 0);
        coneMosaicNormalizingResponse = 1./coneMosaicNullResponse;
        coneMosaicNormalizingResponse(coneIndicesWithZeroNullResponse) = 0;
        coneMosaicNormalizingResponse = reshape(coneMosaicNormalizingResponse, [1 1 numel(coneMosaicNormalizingResponse)]);
    else
        coneMosaicNullResponse = [];
        if (verbose)
            fprintf('Operating on cone excitations\n');
        end
    end

    % Flag that we need to tuck this stuff away on return
    returnTheNoiseFreePipeline = true;
else
    % Get the optics from the previously computed neural pipeline
    theOptics = neuralEngine.neuralPipeline.noiseFreeResponse.optics;

    % Get the mRGC mosaic from the previously computed neural pipeline
    theMRGCmosaic = neuralEngine.neuralPipeline.noiseFreeResponse.mRGCMosaic;

    % Get null response info
    coneMosaicNullResponse = neuralEngine.neuralPipeline.noiseFreeResponse.coneMosaicNullResponse;
    coneMosaicNormalizingResponse = ...
        neuralEngine.neuralPipeline.noiseFreeResponse.coneMosaicNormalizingResponse;

    % No need to store anything
    returnTheNoiseFreePipeline = false;
end

% Compute the sequence of optical images corresponding to the sequence
% of scenes, and convert these into an OISequence
framesNum = numel(sceneSequence);
listOfOpticalImages = cell(1, framesNum);
for frame = 1:framesNum
    listOfOpticalImages{frame} = oiCompute(theOptics, sceneSequence{frame},'padvalue','mean');
end
theOIsequence = oiArbitrarySequence(listOfOpticalImages, sceneSequenceTemporalSupport);

% Compute cone reponses to oiSequence


% Transform the cone excitation responses to cone modulation responses
if (~isempty(neuralEngine.neuralPipeline.noiseFreeResponse.coneMosaicNullResponse))
    % Transform the noise-free cone mosaic response modulation to a contrast response
    % i.e., relative to the cone mosaic response to the null (zero contrast) stimulus.
    % This mimics the photocurrent response which is normalized with respect to the
    % mean cone activation
    theConeMosaicResponses = ...
        bsxfun(@times, bsxfun(@minus, theConeMosaicResponses, ...
        coneMosaicNullResponse), ...
        coneMosaicNormalizingResponse);
end

% Compute the mRGC response
if (verbose)
    fprintf('\tComputing noise-free mRGC responses \n');
end
[theNeuralResponses, ~, temporalSupportSeconds] = theMRGCmosaic.compute( ...
    theConeMosaicResponses, coneMosaicTemporalSupportSeconds);

% Assemble the dataOut struct
dataOut = struct(...
    'neuralResponses', theNeuralResponses, ...
    'temporalSupport', temporalSupportSeconds);

if (returnTheNoiseFreePipeline)
    dataOut.noiseFreeResponsePipeline.optics = theOptics;
    dataOut.noiseFreeResponsePipeline.mRGCMosaic = theMRGCmosaic;
    dataOut.noiseFreeResponsePipeline.mRGCMosaic = theMRGCmosaic;
    dataOut.noiseFreeResponsePipeline.coneMosaicNullResponse = coneMosaicNullResponse;
    dataOut.noiseFreeResponsePipeline.coneMosaicNullResponse.coneMosaicNormalizingResponse ...
        = dataOut.noiseFreeResponsePipeline.coneMosaicNullResponse;
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