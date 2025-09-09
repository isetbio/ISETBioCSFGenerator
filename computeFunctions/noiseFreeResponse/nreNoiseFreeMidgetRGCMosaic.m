function dataOut = nreNoiseFreeMidgetRGCMosaic( ...
    neuralEngineOBJ, noiseFreeComputeParams, sceneSequence, ...
    sceneSequenceTemporalSupport, varargin)
% Compute function for computation of mRGCMosaic activations
%
% Syntax:
%   dataOut = nreNoiseFreeMidgetRGCMosaic(...
%    neuralEngineOBJ, noiseFreeComputeParams, sceneSequence, ...
%    sceneSequenceTemporalSupport, instancesNum, varargin);
%
% Description:
%    Function serving as the computeFunctionHandle for a @neuralResponseEngine
%    object. There are 2 ways to use this function.
%
%       [1] If called directly and with no arguments,
%           dataOut = nreScenePhotonNoise()
%           dataOut = nreScenePhotonNoise([],[],[],[],varargin)

%       it does not compute anything and simply returns a struct with the
%       defaultParams (optics and coneMosaic params) that define the neural
%       compute pipeline for this computation.  In the second usage form,
%       key/value pairs might be used to control the default parameters.
%
%       [2] If called from a parent @neuralResponseEngine object,
%       it computes 'instancesNum' of midget RGC response sequences
%       in response to the passed 'sceneSequence'.
%
%    It is not a good idea to try to call this function with arguments
%    directly - it should be called by the compute method of its parent
%    @neuralResponseEngine.
%
%    In addition to computing, this function checks the `visualizeEachCompute` 
%    flag of the neuralEngineOBJ and, if it set, calls the nreVisualizeMRGCmosaic()
%    visualization function. This causes figures to appear that visualize
%    the noise-free spatiotemporal activation of the mRGC mosaic, 
%    which is helpful for debugging.
%    Note that everything runs much more slowly in this case.

% Inputs:
%    neuralEngineOBJ                   - the parent @neuralResponseEngine object that
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
%             - If called directly with no input arguments or just
%               key/value pairs, the returned struct contains
%               the defaultParams. These describe the optics and cMosaic.
%               Passing key/value pairs allows control of the details of
%               the components generated, etc.
%
%             - If called from a parent @neuralResponseEngine), the returned
%               struct is organized as follows:
%%
%              .neuralResponses : matrix of neural responses.  The first dimension
%                                 is number of instances.  Often there is just one.
%                                 Each "column" of the second dimension 
%                                 is the response vector for one frame of the input.
%                                 The third dimension indexes the frames.
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
%   'fixationalEM'         - Empty (default) or a fixationalEM object
%                            that describes one eye movement path.  If
%                            the latter, this must have one position per
%                            frame of the passed scene sequence, in which
%                            case it is applied.
%   'opticsType'           - String or struct (default 'loadComputeReadyRGCMosaic'). Specify type
%                            of optics to use.
%                               - 'loadComputeReadyRGCMosaic' - Produced by
%                               the loadComputeReadRGCMosaic machinery.
%                               - 'oiEnsembleGenerate' - Reasonable human
%                               wavefront optics.
%                               - 'BerkeleyAO' - Models optics of subject in
%                               some of the Berkeley AO psychophysics
%                               systems.
%                             If it is a struct, it is taken to be the
%                             desired oi to use
%   'verbose'               - Logical (default true). Print things out.
%   'oiPadMethod'           - String (default 'zero'). Passed to the oi
%                             compute method to determine how image is
%                             padded for convolution with PSF.  Options are
%                             'zero' and 'mean'.
%
% See Also:
%   t_neuralResponseCompute

% History:
%    05/04/23  NPC  Wrote it
%    11/03/24  FH   Edited it to incorporate metaContrast
%    12/23/24  dhb  Rewrite for major architecture redo.

% Parse the input arguments
%
% Allow for possibility that other nre's take key/value pairs that we
% can ignore, so set KeepUnmatched to true.
p = inputParser;
p.KeepUnmatched = true;
p.addParameter('fixationalEM', [], @(x)(isempty(x) || (isa(x,'fixationalEM'))));
p.addParameter('opticsType','loadComputeReadyRGCMosaic',@(x)(ischar(x) | isstruct(x)));
p.addParameter('oiPadMethod','zero',@ischar);
p.addParameter('verbose',true,@islogical);
p.addParameter('visualizeActivationFunction', true, @islogical);
varargin = ieParamFormat(varargin);
p.parse(varargin{:});
fixationalEMObj = p.Results.fixationalEM;
opticsType = p.Results.opticsType;
oiPadMethod = p.Results.oiPadMethod;
verbose = p.Results.verbose;
visualizeActivationFunction = p.Results.visualizeActivationFunction;

% Check input arguments. If called with zero input arguments, just return the default params struct
if (nargin == 0 | isempty(neuralEngineOBJ))
    dataOut = generateDefaultParams(opticsType,oiPadMethod);
    return;
end

% Get number of frames
framesNum = numel(sceneSequence);

% Create/get key objects on first call
if (isempty(neuralEngineOBJ.neuralPipeline) | ~isfield(neuralEngineOBJ.neuralPipeline,'noiseFreeResponse'))
    % Step1. Load a compute-ready mRGC mosaic given the passed params.
    %
    % This also produces the cone mosaic on which the RGC mosaic is built.
    % The calling routine generates optics as well.
    [theOptics,theMRGCmosaic] = generateOpticsAndMosaicFromParams(...
        noiseFreeComputeParams.opticsParams, ...
        [], ...
        noiseFreeComputeParams.mRGCMosaicParams);

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

    % Handle contrast versus excitations
    if (strcmp(noiseFreeComputeParams.mRGCMosaicParams.inputSignalType,'coneContrast'))
        % Compute theConeMosaicNullResponse if the inputSignalType is set to
        % 'coneContrast' and if we have nullStimulusSceneSequence.
        if (~isempty(noiseFreeComputeParams.nullStimulusSceneSequence))

            if (verbose)
                fprintf('Operating on cone modulations\n');
            end

            % Retrieve the null stimulus scene sequence
            nullStimulusSceneSequence = noiseFreeComputeParams.nullStimulusSceneSequence;

            % Compute the optical image of the null scene
            listOfNullOpticalImages = cell(1, framesNum);
            for frame = 1:framesNum
                listOfNullOpticalImages{frame} = oiCompute(theOptics, nullStimulusSceneSequence{frame},'padvalue',oiPadMethod);
            end
            nullOIsequence = oiArbitrarySequence(listOfNullOpticalImages, sceneSequenceTemporalSupport);
            clear listOfNullOpticalImages;

            % Compute theConeMosaicNullResponse, i.e., the input cone mosaic response to the NULL scene
            coneMosaicNullResponse = theMRGCmosaic.inputConeMosaic.compute(...
                nullOIsequence, ...
                'nTrials', 1);

            % Compute normalizing response (for computing the  modulated response)
            coneIndicesWithZeroNullResponse = find(coneMosaicNullResponse == 0);
            coneMosaicNormalizingResponse = 1./coneMosaicNullResponse;
            coneMosaicNormalizingResponse(coneIndicesWithZeroNullResponse) = 0;
        else
            error('Input type ''coneContrast'' specified, but no normalizing stimulus provided.');
        end

    elseif (strcmp(noiseFreeComputeParams.mRGCMosaicParams.inputSignalType,'coneExcitations'))
        % Operating on cone excitations
        coneMosaicNullResponse = [];
        coneMosaicNormalizingResponse = [];
        if (verbose)
            fprintf('Operating on cone excitations\n');
        end
    else
        error('The input signal type must be ''coneContrasts'' or ''coneExciations''');
    end

    % Flag that we need to tuck this stuff away on return
    returnTheNoiseFreePipeline = true;
else
    % Get the optics from the previously computed neural pipeline
    theOptics = neuralEngineOBJ.neuralPipeline.noiseFreeResponse.optics;

    % Get the mRGC mosaic from the previously computed neural pipeline
    theMRGCmosaic = neuralEngineOBJ.neuralPipeline.noiseFreeResponse.mRGCMosaic;

    % Get null response info
    coneMosaicNullResponse = neuralEngineOBJ.neuralPipeline.noiseFreeResponse.coneMosaicNullResponse;
    coneMosaicNormalizingResponse = ...
        neuralEngineOBJ.neuralPipeline.noiseFreeResponse.coneMosaicNormalizingResponse;

    % No need to store anything
    returnTheNoiseFreePipeline = false;
end

% Compute the sequence of optical images corresponding to the sequence
% of scenes, and convert these into an OISequence
framesNum = numel(sceneSequence);
listOfOpticalImages = cell(1, framesNum);
for frame = 1:framesNum
    listOfOpticalImages{frame} = oiCompute(theOptics, sceneSequence{frame},'padvalue',oiPadMethod);
end
theOIsequence = oiArbitrarySequence(listOfOpticalImages, sceneSequenceTemporalSupport);
clear listOfOpticalImages;

% Compute cone reponses to oiSequence
if (verbose)
    fprintf('\tComputing noise-free cone mosaic responses\n');
end

if (isempty(fixationalEMObj))
    [noiseFreeConeMosaicResponses, ~, ~, ~, temporalSupportSeconds] = ...
        theMRGCmosaic.inputConeMosaic.compute(theOIsequence, ...
        'nTrials', 1, ...
        'withFixationalEyeMovements', false ...
        );
else
    % We were passed a fixational EM.  Check that it is OK for our
    % purposes, and then comptute with it.

    % Also return the path in units of microns
    fixationalEMObj.emPosMicrons = ...
        theMRGCmosaic.inputConeMosaic.distanceDegreesToDistanceMicronsForCmosaic(fixationalEMObj.emPosArcMin / 60);

    % Compute
    theMRGCmosaic.inputConeMosaic.emSetFixationalEMObj(fixationalEMObj);
    [noiseFreeConeMosaicResponses, ~, ~, ~, temporalSupportSeconds] = ...
        theMRGCmosaic.inputConeMosaic.compute(theOIsequence, ...
        'nTrials', 1, ...
        'withFixationalEyeMovements', true ...
        );
end

% Transform the cone excitation responses to cone modulation responses if
% needed.
if (~isempty(coneMosaicNullResponse))
    % Transform the noise-free cone mosaic response modulation to a contrast response
    % i.e., relative to the cone mosaic response to the null (zero contrast) stimulus.
    % This mimics the photocurrent response which is normalized with respect to the
    % mean cone activation
    
    noiseFreeConeMosaicResponses = ...
            bsxfun(@times, bsxfun(@minus, noiseFreeConeMosaicResponses, ...
            coneMosaicNullResponse), ...
            coneMosaicNormalizingResponse);
end

% Handle what kind of output we are asked for
switch (noiseFreeComputeParams.mRGCMosaicParams.outputSignalType)
    case 'cones'
        if (verbose)
            fprintf('\tReturning cone signals rather than mRGC responses\n');
        end

        % Just use the cone signals.  We have their temporal support set
        % above.
        theNeuralResponses = noiseFreeConeMosaicResponses;
        theNeuralResponses = permute(theNeuralResponses,[1, 3, 2]);
    
    case 'mRGCs'
        % Compute the mRGC response
        if (verbose)
            fprintf('\tComputing noise-free mRGC responses \n');
        end
        [theNeuralResponses, ~, temporalSupportSeconds] = theMRGCmosaic.compute( ...
            noiseFreeConeMosaicResponses, temporalSupportSeconds);
        theNeuralResponses = permute(theNeuralResponses,[1, 3, 2]);
end

%% Apply temporal filter if needed
if (~isempty(noiseFreeComputeParams.temporalFilter))
    filterTemporalSupport = noiseFreeComputeParams.temporalFilter.temporalSupport;
    filterValues = noiseFreeComputeParams.temporalFilter.filterValues;

    % Loop over instances and responses
    nInstances = size(theNeuralResponses,1);
    mResponses = size(theNeuralResponses,2);
    nTimePoints = size(theNeuralResponses,3);

    newNeuralResponses = zeros(nInstances,mResponses,nTimePoints);
    for ii = 1:nInstances
        for jj = 1:mResponses
            % Get the temporal response for this instance/response dim
            theResponse = squeeze(theNeuralResponses(ii,jj,:));

            % Apply temporal filter.  We don't explicitly pad here.  If you
            % want the end of the response generated by the filter, you'll
            % want to pad the input with extra zero frames.
            filteredResponse = conv(theResponse,filterValues(:));

            % Just return the values at time points where we had input.
            % I had thought that this is what calling the conv function
            % with the 'same' argument would do, but it picks out the
            % center of the full return, which is not what we want.
            filteredResponse = filteredResponse(1:nTimePoints);

            % Put it back
            newNeuralResponses(ii,jj,:) = filteredResponse;
        end
    end

    % Get new values into output variable
    theNeuralResponses = newNeuralResponses;
    clear newNeuralResponses;
end

% If 'simulateONOFFmosaic' has been set, flip the response sign of the odd-numbered mRGCs
if (strcmp(noiseFreeComputeParams.mRGCMosaicParams.outputSignalType, 'mRGCs'))
    if (isfield(noiseFreeComputeParams.mRGCMosaicParams,'simulateONOFFmosaic'))
        if (noiseFreeComputeParams.mRGCMosaicParams.simulateONOFFmosaic)
            mResponses = size(theNeuralResponses,2);
            fprintf('Flipping signs of odd-numbered RGCs\n')
            oddNumberedRGCs = 1:2:mResponses;
            theNeuralResponses(:,oddNumberedRGCs,:) = -theNeuralResponses(:,oddNumberedRGCs,:);
        end
    end
end

% If we have a halfwave rectifier nonlinear activation function and a response bias has been specified, 
% apply the response bias to the noise-free responses.
% In the current implementation, the actual activation function is applied to the noisy response instances
% See comments in @neuralResponseEngine.applyActivationFunction
if (isfield(noiseFreeComputeParams.mRGCMosaicParams, 'responseBias'))
    theNeuralResponses = theNeuralResponses + noiseFreeComputeParams.mRGCMosaicParams.responseBias;
end


% Check the visualizeEachCompute flag of the neuralEngineOBJ, and if set to true,
% call the neuralEngineOBJ.visualize() function to visualize the responses

% Empty the visualization meta data
visualizationMetaData = [];

if (neuralEngineOBJ.computeVisualizationMetaData || neuralEngineOBJ.visualizeEachCompute)
    % Save the cone mosaic responses in the visualizationMetaData struct
    visualizationMetaData.noiseFreeConeMosaicResponses = noiseFreeConeMosaicResponses;
end

% Update the stored visualization metadata
neuralEngineOBJ.updataVisualizationMetadata(visualizationMetaData);

% Visualize
if (neuralEngineOBJ.visualizeEachCompute)

    hFig = figure(1000);
    set(hFig, 'Position', [350 700 1650 550]);

    neuralEngineOBJ.visualize(theNeuralResponses, temporalSupportSeconds, ...
        'responseLabel', 'noise-free_mRGC', ...
        'figureHandle', hFig);
end

% Assemble the dataOut struct
dataOut = struct(...
    'neuralResponses', theNeuralResponses, ...
    'temporalSupport', temporalSupportSeconds);


if (returnTheNoiseFreePipeline)
    noiseFreeResponsePipeline = struct();
    noiseFreeResponsePipeline.optics = theOptics;
    noiseFreeResponsePipeline.mRGCMosaic = theMRGCmosaic;
    noiseFreeResponsePipeline.coneMosaicNullResponse = coneMosaicNullResponse;
    noiseFreeResponsePipeline.coneMosaicNormalizingResponse = coneMosaicNormalizingResponse;
    dataOut.noiseFreeResponsePipeline = noiseFreeResponsePipeline;
end

end

function p = generateDefaultParams(opticsType,oiPadMethod)

% Genrate optics params through common bottleneck
opticsParams = generateOpticsParams(opticsType,oiPadMethod);

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

mosaicParams = struct(...
    'type', 'mRGCMosaic', ...
    'spatialCompactnessSpectralPurityTradeoff', 1.0, ...
    'surroundOptimizationSubString', 'PackerDacey2002H1freeLowH1paramsNarrowVisualSTFparamTolerance_vSTF_1.0_1.0', ...
    'eccDegs', [7 0], ...
    'sizeDegs',  [6 3], ...
    'rgcType', 'ONcenterMidgetRGC', ...
    'cropParams', cropParams, ...
    'retinalRFmodelParams', retinalRFmodelParams, ...
    'inputSignalType', 'coneExcitations', ...
    'outputSignalType', 'mRGCs', ...
    'responseBias', 0, ...
    'coneIntegrationTimeSeconds', 10/1000);

% If the user sets the nullStimulusSceneSequence, then mRGC responses are
% computed on modulated cone responses with respect to the cone
% response to the null stimulus scene sequence.
nullStimulusSceneSequence = [];

% Temporal filter
%
% The user can supply a temporal filter.  It is applied
% to each output sequence, on the assumption that its time
% base matches that of the signal, but we do pass the filter
% time base so that this could be generalized in the future.
temporalFilter = [];

% Assemble all params in a struct
p = struct(...
    'opticsParams', opticsParams, ...
    'mRGCMosaicParams', mosaicParams, ...
    'nullStimulusSceneSequence', nullStimulusSceneSequence, ...
    'temporalFilter', temporalFilter ...
    );

end