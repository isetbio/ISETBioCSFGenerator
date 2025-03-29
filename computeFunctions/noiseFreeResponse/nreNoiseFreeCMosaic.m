function dataOut = nreNoiseFreeCMosaic(...
    neuralEngineOBJ, noiseFreeComputeParams, sceneSequence, ...
    sceneSequenceTemporalSupport, varargin)
% Compute function for computation of cone excitations
%
% Syntax:
%   dataOut = nreNoiseFreeCMosaic(...
%    neuralEngineOBJ, noiseFreeComputeParams, sceneSequence, ...
%    sceneSequenceTemporalSupport, varargin);
%
% Description:
%    Function serving as the computeFunctionHandle for a @neuralResponseEngine
%    object.  This version uses the cone exciations as its neural response.
%
%    There are 2 ways to use this function.
%
%       [1] If called directly and with no arguments,
%           dataOut = nreScenePhotonNoise()
%           dataOut = nreScenePhotonNoise([],[],[],[],varargin)
%
%       it does not compute anything and simply returns a struct with the
%       defaultParams (optics and coneMosaic params) that define the neural
%       compute pipeline for this computation.  In the second usage form,
%       key/value pairs might be used to control the default parameters.
%
%       [2] If called from a parent @neuralResponseEngine object,
%       it computes 'instancesNum' of cone photopigment excitation sequences
%       in response to the passed 'sceneSequence'.
%
%    It is not a good idea to try to call this function with arguments
%    directly - it should be called by the compute method of its parent
%    @neuralResponseEngine.
%
%    In addition to computing, this function checks the `visualizeEachCompute`
%    flag of the neuralEngineOBJ and, if it set, calls the nreVisualizeCMosaic()
%    visualization function. This causes figures to appear that visualize
%    the noise-free spatiotemporal activation of the cone mosaic, along with
%    any eye movements that may be in effect, which is helpful for debugging.
%    Note that everything runs much more slowly in this case.
%
% Inputs:
%    neuralEngineOBJ                - the parent @neuralResponseEngine object that
%                                     is calling this function as its computeFunctionHandle
%    noiseFreeComputeParams         - a struct containing properties of the employed neural chain.
%    sceneSequence                  - a cell array of scenes defining the frames of a stimulus
%    sceneSequenceTemporalSupport   - the temporal support for the stimulus frames, in seconds
%
% Outputs:
%    dataOut  - A struct that depends on the input arguments.
%
%             - If called directly with no input arguments or just
%               key/value pairs, the returned struct contains
%               the defaultParams. These describe the optics and cMosaic.
%               Passing key/value pairs allows control of the details of
%               the components generated, etc.
%
%             - If called from a parent @neuralResponseEngine), the returned
%               struct is organized as follows:
%
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
%   Note, to pass key/value pairs you need to pass something for the first
%   four arguments.  If you are just trying to get the default parameters,
%   these can all be empty.
%
%   'fixationalEM'         - Empty (default) or a fixationalEM object
%                            that describes one eye movement path.  If
%                            the latter, this must have one position per
%                            frame of the passed scene sequence, in which
%                            case it is applied.
%   'opticsType'           - String or struct (default 'oiEnsembleGenerate'). Specify type
%                            of optics to use.
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
%     t_neuralResponseCompute, t_spatialCSF

% History:
%    03/29/2021  npc  Wrote it by adapting nrePhotopigmentExcitationsConeMosaicHexWithNoEyeMovements
%    04/10/2023  fh   Edited it so that with correctly set varargin, this
%                       function behaves the same as
%                       nrePhotopigmentExcitationsCMosaicWithNoEyeMovements.m
%    12/19/2024  dhb  Rewrite for major architecture redo.

% Examples:
%{
    % Clear to avoid confusion
    clear;

    % Instantiate a @sceneEngine object and generate a test scene
    theSceneEngine = sceneEngine(@sceUniformFieldTemporalModulation);
    testContrast = 0.1;
    [sceneSequence, sceneTemporalSupportSeconds] = ...
        theSceneEngine.compute(testContrast);

    % Usage case #1. Just return the default neural response params
    nreParams = nreNoiseFreeCMosaic()
    nreParams.coneMosaicParams.timeIntegrationSeconds = ...
        sceneTemporalSupportSeconds(2)-sceneTemporalSupportSeconds(1);
    
    % Usage case #2. Compute noise free responses 
    % using a parent @neuralResponseEngine object and the default neural response params
    %
    % Instantiate the parent @neuralResponseEngine object, with Poisson noise
    theneuralEngineOBJ = neuralResponseEngine(@nreNoiseFreeCMosaic, ...
        @nreNoisyInstancesPoisson, ...
        nreParams, []);
    
    % Compute noise free response
    [noiseFreeResponse, temporalSupportSeconds] = theneuralEngineOBJ.computeNoiseFree(...
            sceneSequence, ...
            sceneTemporalSupportSeconds ...
            );

    % Add Poisson noise
    instancesNum = 16;
    [noisyInstances, temporalSupportSeconds] = theneuralEngineOBJ.computeNoisyInstances(...
            noiseFreeResponse, ...
            temporalSupportSeconds, ...
            instancesNum, ...
            'random' ...
            );

    % Get multiple noise free instances, by not adding noise
    instancesNum = 16;
    [noNoiseInstances, temporalSupportSeconds] = theneuralEngineOBJ.computeNoisyInstances(...
            noiseFreeResponse, ...
            temporalSupportSeconds, ...
            instancesNum, ...
            'none' ...
            );
%}

% Parse the input arguments.
%
% Allow for possibility that other nre's take key/value pairs that we
% can ignore, so set KeepUnmatched to true.
p = inputParser;
p.KeepUnmatched = true;

% The nre knows about these and can call this function with them set
% appropriately.
p.addParameter('fixationalEM', [], @(x)(isempty(x) || (isa(x,'fixationalEM'))));
p.addParameter('verbose',true,@islogical);

% These affect the default parameters returned
p.addParameter('opticsType','oiEnsembleGenerate',@(x)(ischar(x) | isstruct(x)));
p.addParameter('oiPadMethod','zero',@ischar)
varargin = ieParamFormat(varargin);
p.parse(varargin{:});
fixationalEMObj = p.Results.fixationalEM;
opticsType = p.Results.opticsType;
oiPadMethod = p.Results.oiPadMethod;
verbose = p.Results.verbose;


% Check input arguments. If called with zero input arguments, just return the default params struct
if (nargin == 0 | isempty(neuralEngineOBJ))
    dataOut = generateDefaultParams(opticsType,oiPadMethod);
    return;
end


% Get the number of scene sequence frames
framesNum = numel(sceneSequence);

% Create/get key objects
if (isempty(neuralEngineOBJ.neuralPipeline) | ~isfield(neuralEngineOBJ.neuralPipeline,'noiseFreeResponse'))

    % Generate optics and mosaic
    [theOptics,theConeMosaic] = generateOpticsAndMosaicFromParams(...
        noiseFreeComputeParams.opticsParams, ...
        [], ...
        noiseFreeComputeParams.coneMosaicParams);

    % Handle contrast versus excitations
    if (strcmp(noiseFreeComputeParams.coneMosaicParams.outputSignalType,'coneContrast'))
        % Compute theConeMosaicNullResponse if the outputSignalType is set to
        % 'coneContrast' and if we have nullStimulusSceneSequence.
        if (~isempty(noiseFreeComputeParams.nullStimulusSceneSequence))
            if (verbose)
                fprintf('Operating on cone contrast\n');
            end

            % Retrieve the null stimulus scene sequence
            nullStimulusSceneSequence = noiseFreeComputeParams.nullStimulusSceneSequence;

            % Compute the optical image of the null scene
            listOfNullOpticalImages = cell(1, framesNum);
            for frame = 1:framesNum
                listOfNullOpticalImages{frame} = oiCompute(theOptics, nullStimulusSceneSequence{frame},'padvalue',noiseFreeComputeParams.opticsParams.oiPadMethod);
            end
            nullOIsequence = oiArbitrarySequence(listOfNullOpticalImages, sceneSequenceTemporalSupport);
            clear listOfNullOpticalImages;

            % Compute theConeMosaicNullResponse, i.e., the input cone mosaic response to the NULL scene
            coneMosaicNullResponse = theConeMosaic.compute(...
                nullOIsequence, ...
                'nTrials', 1);

            % Compute normalizing response (for computing the  modulated response)
            coneIndicesWithZeroNullResponse = find(coneMosaicNullResponse == 0);
            coneMosaicNormalizingResponse = 1./coneMosaicNullResponse;
            coneMosaicNormalizingResponse(coneIndicesWithZeroNullResponse) = 0;
        else
            error('Input type ''coneContrast'' specified, but no normalizing stimulus provided.');
        end

    elseif (strcmp(noiseFreeComputeParams.coneMosaicParams.outputSignalType,'coneExcitations'))
        % Operating on cone excitations
        coneMosaicNullResponse = [];
        coneMosaicNormalizingResponse = [];
        if (verbose)
            fprintf('Operating on cone excitations\n');
        end
    else
        error('The input signal type must be ''coneContrasts'' or ''coneExciations''');
    end

    returnTheNoiseFreePipeline = true;
else
    % Get the optics from the previously computed neural pipeline
    % stored in the object
    theOptics = neuralEngineOBJ.neuralPipeline.noiseFreeResponse.optics;

    % Get the cone mosaic from the previously computed neural pipeline
    theConeMosaic = neuralEngineOBJ.neuralPipeline.noiseFreeResponse.coneMosaic;

    % Get the info for converting to contrast
    coneMosaicNullResponse =  neuralEngineOBJ.neuralPipeline.noiseFreeResponse.coneMosaicNullResponse;
    coneMosaicNormalizingResponse = ...
        neuralEngineOBJ.neuralPipeline.noiseFreeResponse.coneMosaicNormalizingResponse;

    % No need to store anything
    returnTheNoiseFreePipeline = false;
end

% Compute the sequence of optical images corresponding to the sequence of scenes
if framesNum == 1
    theOIsequence = oiCompute(theOptics, sceneSequence{1},'padvalue',noiseFreeComputeParams.opticsParams.oiPadMethod);
else
    theListOfOpticalImages = cell(1, framesNum);
    for frame = 1:framesNum
        theListOfOpticalImages{frame} = oiCompute(theOptics, sceneSequence{frame},'padvalue',noiseFreeComputeParams.opticsParams.oiPadMethod);
    end

    % Generate an @oiSequence object containing the list of computed optical images
    theOIsequence = oiArbitrarySequence(theListOfOpticalImages, sceneSequenceTemporalSupport);
end

% Compute responses depending on whether there is a fixationalEM or not.
%
% NEED TO CHECK THAT ONLY ONE EM PATH WAS PASSED.  THAT IS THE ONLY CASE WE
% ALLOW AT THIS LEVEL.  ALSO IN OTHER NRE'S.
if (isempty(fixationalEMObj))
    % No fixationalEM passed. So don't move the eyes.
    [theNeuralResponses, ~, ~, ~, temporalSupportSeconds] = ...
        theConeMosaic.compute(theOIsequence, ...
        'nTrials', 1, ...
        'withFixationalEyeMovements', false ...
        );
else
    % We were passed a fixational EM.  Check that it is OK for our
    % purposes, and then comptute with it.

    % Also return the path in units of microns
    fixationalEMObj.emPosMicrons = ...
        theConeMosaic.distanceDegreesToDistanceMicronsForCmosaic(fixationalEMObj.emPosArcMin / 60);

    % Compute
    theConeMosaic.emSetFixationalEMObj(fixationalEMObj);
    [theNeuralResponses, ~, ~, ~, temporalSupportSeconds] = ...
        theConeMosaic.compute(theOIsequence, ...
        'nTrials', 1, ...
        'withFixationalEyeMovements', true ...
        );
end

% Transform the cone excitation responses to cone modulation responses if
% needed.
if (strcmp(noiseFreeComputeParams.coneMosaicParams.outputSignalType,'coneContrast'))
    % Transform the noise-free cone mosaic response modulation to a contrast response
    % i.e., relative to the cone mosaic response to the null (zero contrast) stimulus.
    % This mimics the photocurrent response which is normalized with respect to the
    % mean cone activations.
    theNeuralResponses = ...
        bsxfun(@times, bsxfun(@minus, theNeuralResponses, ...
        coneMosaicNullResponse), ...
        coneMosaicNormalizingResponse);
end

%% Put neural responses into the right format
%
% Convert cMosaic return convention to nre return convention.  We have
% to deal with unfortunate special casing because of the way Matlab
% and/or cMosaic handle singleton dimensions, so that the returned
% responses always have the responses in the column(s).
%
% Be sure to do this before applying temporal filter, so that is applied
% along the expected direction.
theNeuralResponses = permute(theNeuralResponses,[1 3 2]);

%% Apply temporal filter if needed
if (~isempty(noiseFreeComputeParams.temporalFilter))
    filterValues = noiseFreeComputeParams.temporalFilter.filterValues;

    % If the filterValues is a string, we need to compute the filter values
    if (ischar(filterValues))
        switch (filterValues)
            case 'photocurrentImpulseResponseBased'
                if (~exist('coneMosaicNullResponse', 'var'))||(isempty(coneMosaicNullResponse))
                    if (isempty(noiseFreeComputeParams.nullStimulusSceneSequence))
                        error('Specified ''%s'' temporal filter, but did not provide a nullStimulus (background) which is necessary to compute the photocurrent impulse repsonse .', filterValues);
                    else
                        % Retrieve the null stimulus scene sequence
                        nullStimulusSceneSequence = noiseFreeComputeParams.nullStimulusSceneSequence;

                        % Compute the optical image of the null scene
                        listOfNullOpticalImages = cell(1, framesNum);
                        for frame = 1:framesNum
                            listOfNullOpticalImages{frame} = oiCompute(theOptics, nullStimulusSceneSequence{frame},'padvalue',noiseFreeComputeParams.opticsParams.oiPadMethod);
                        end
                        nullOIsequence = oiArbitrarySequence(listOfNullOpticalImages, sceneSequenceTemporalSupport);
                        clear listOfNullOpticalImages;

                        % Compute theConeMosaicNullResponse, i.e., the input cone mosaic response to the NULL scene
                        coneMosaicNullResponse = theConeMosaic.compute(...
                            nullOIsequence, ...
                            'nTrials', 1);
                    end
                end % if (~exist('coneMosaicNullResponse', 'var'))|| (isempty(coneMosaicNullResponse)

                % Call the standalone photocurrent impulse response compute function
                % The contrast of the impulse stimulus
                theConeFlashImpulseContrast = 0.01*[1 1 1];

                % Visualize the computed impulse responses only if the
                % neuralEngineOBJ.visualizeEachCompute flag has been set to true
                visualizePhotocurrentImpulseResponses = neuralEngineOBJ.visualizeEachCompute;

                % Compute the impulse responses
                thePhotocurrentImpulseResponseStruct = CMosaicNrePhotocurrentImpulseResponses(...
                    theConeMosaic, coneMosaicNullResponse, theConeFlashImpulseContrast, framesNum, ...
                    visualizePhotocurrentImpulseResponses);

                % Retrieve the computed filterValues and its temporal support from the returned photocurrentImpulseResponseStruct
                filterTemporalSupport = thePhotocurrentImpulseResponseStruct.temporalSupportSeconds;

                % Ensure temporal support are consistent
                assert(all(size(filterTemporalSupport) == size(temporalSupportSeconds)), 'mismatch in temporal support lengths');
                assert(all(temporalSupportSeconds == filterTemporalSupport), 'mismatch in temporal support values');

                filterValues = thePhotocurrentImpulseResponseStruct.coneDensityWeightedPhotocurrentImpulseResponse;

                % Save the temporal filter for returning it in dataOut
                noiseFreeComputeParams.temporalFilter.temporalSupport = filterTemporalSupport;
                noiseFreeComputeParams.temporalFilter.filterValues = filterValues;

                % Update the noiseFreeComputeParams with the newly computed temporal filter so we can cache it and not have to recompute it again
                neuralEngineOBJ.updateParamsStruct(noiseFreeComputeParams, neuralEngineOBJ.noisyInstancesComputeParams);

                %fprintf(2, 'We are re-computing the filter values\n')
                %pause
            otherwise
                error('Unsupported temporal filter method: ''%s''.', filterValues);
        end % switch

    else % ~ischar(filterValues)
        %fprintf(2, 'We retrieve the CACHED filter values instead of recomputing them !!!')
        %pause
        filterTemporalSupport = noiseFreeComputeParams.temporalFilter.temporalSupport;
    end

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

    % Reaality check, only applies if you have set the filter equal to a
    % temporal delta function, so generally is kept at false.
    CHECK = false;
    if (CHECK)
        if (max(abs(theNeuralResponses(:)-newNeuralResponses(:))) > 1e-10)
            error('Filter did something it should not have');
        end
    end

    % Update responses with filtered version
    theNeuralResponses = newNeuralResponses;
    clear newNeuralResponses;
end

%% Detailed visualization option
%
% Check the visualizeEachCompute flag of the neuralEngineOBJ, and if set to true,
% call the neuralEngineOBJ.visualize() function to visualize the responses
%
% Empty the visualization meta data
visualizationMetaData = [];

if (neuralEngineOBJ.computeVisualizationMetaData || neuralEngineOBJ.visualizeEachCompute)
    % Placeholder here for some future data
    % visualizationMetaData.someField = someData;
end

% Update the stored visualization metadata
neuralEngineOBJ.updataVisualizationMetadata(visualizationMetaData);

% Do the detailed visualization, if requested.
if (neuralEngineOBJ.visualizeEachCompute)
    hFig = figure(1000);
    set(hFig, 'Position', [350 700 1650 550]);
    neuralEngineOBJ.visualize(theNeuralResponses, temporalSupportSeconds, ...
        'responseLabel', 'noise-free cMosaic responses', ...
        'figureHandle', hFig);
end

%% Assemble the dataOut struct
dataOut = struct(...
    'neuralResponses', theNeuralResponses, ...
    'temporalSupport', temporalSupportSeconds);

if (returnTheNoiseFreePipeline)
    dataOut.noiseFreeResponsePipeline.optics = theOptics;
    dataOut.noiseFreeResponsePipeline.coneMosaic = theConeMosaic;
    dataOut.noiseFreeResponsePipeline.temporalFilter = noiseFreeComputeParams.temporalFilter;
    dataOut.noiseFreeResponsePipeline.coneMosaicNullResponse = coneMosaicNullResponse;
    dataOut.noiseFreeResponsePipeline.coneMosaicNormalizingResponse = coneMosaicNormalizingResponse;
end
end

function p = generateDefaultParams(opticsType,oiPadMethod)

% Generate optics parameters
opticsParams = generateOpticsParams(opticsType,oiPadMethod);

% cMosaic params
mosaicParams = struct(...
    'type', 'cMosaic', ...
    'wave', 400:10:700, ...
    'sizeDegs', 0.3*[1 1], ...
    'eccDegs', [0 0], ...
    'timeIntegrationSeconds', 5/1000, ...
    'outputSignalType', 'coneExcitations' ...
    );

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

% Default params for this compute function
p = struct(...
    'opticsParams', opticsParams, ...
    'coneMosaicParams', mosaicParams, ...
    'nullStimulusSceneSequence', nullStimulusSceneSequence, ...
    'temporalFilter', temporalFilter ...
    );
end
