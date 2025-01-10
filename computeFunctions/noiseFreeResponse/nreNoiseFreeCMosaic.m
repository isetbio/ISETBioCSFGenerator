function dataOut = nreNoiseFreeCMosaic(...
    neuralEngine, noiseFreeComputeParams, sceneSequence, ...
    sceneSequenceTemporalSupport, varargin)
% Compute function for computation of cone excitations
%
% Syntax:
%   dataOut = nreNoiseFreeCMosaic(...
%    neuralEngine, noiseFreeComputeParams, sceneSequence, ...
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
% Inputs:
%    neuralEngine                   - the parent @neuralResponseEngine object that
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
    theNeuralEngine = neuralResponseEngine(@nreNoiseFreeCMosaic, ...
        @nreNoisyInstancesPoisson, ...
        nreParams, []);
    
    % Compute noise free response
    [noiseFreeResponse, temporalSupportSeconds] = theNeuralEngine.computeNoiseFree(...
            sceneSequence, ...
            sceneTemporalSupportSeconds ...
            );

    % Add Poisson noise
    instancesNum = 16;
    [noisyInstances, temporalSupportSeconds] = theNeuralEngine.computeNoisyInstances(...
            noiseFreeResponse, ...
            temporalSupportSeconds, ...
            instancesNum, ...
            'random' ...
            );

    % Get multiple noise free instances, by not adding noise
    instancesNum = 16;
    [noNoiseInstances, temporalSupportSeconds] = theNeuralEngine.computeNoisyInstances(...
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
if (nargin == 0 | isempty(neuralEngine))
    dataOut = generateDefaultParams(opticsType,oiPadMethod);
    return;
end

% Get the number of scene sequence frames
framesNum = numel(sceneSequence);

% Create/get key objects
if (isempty(neuralEngine.neuralPipeline) | ~isfield(neuralEngine.neuralPipeline,'noiseFreeResponse'))

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
    theOptics = neuralEngine.neuralPipeline.noiseFreeResponse.optics;

    % Get the cone mosaic from the previously computed neural pipeline
    theConeMosaic = neuralEngine.neuralPipeline.noiseFreeResponse.coneMosaic;
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

% Compute responses depending on whether there is a fixationalEM or
% not.
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
    %
    % Check NEED TO ADD

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
            filteredResponse = conv(theResponse,filterValues,'same');

            % Put it back
            newNeuralResponses(ii,jj,:) = filteredResponse;
        end
    end

    % Get new values into output variable
    theNeuralResponses = newNeuralResponses;
    clear newNeuralResponses;
end


% Visualize responses
if (neuralEngine.visualizeEachCompute)
    visualizeNoiseFreeNeuralResponse(theConeMosaic, theNeuralResponses, temporalSupportSeconds);
end


% Convert cMosaic return convention to nre return convention.  We have
% to deal with unfortunate special casing because of the way Matlab
% and/or cMosaic handle singleton dimensions, so that the returned
% responses always have the responses in the column(s).
theNeuralResponses = permute(theNeuralResponses,[1 3 2]);



% Assemble the dataOut struct
dataOut = struct(...
    'neuralResponses', theNeuralResponses, ...
    'temporalSupport', temporalSupportSeconds);

if (returnTheNoiseFreePipeline)
    dataOut.noiseFreeResponsePipeline.optics = theOptics;
    dataOut.noiseFreeResponsePipeline.coneMosaic = theConeMosaic;
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


% Visualization function for nreNoiseFreeCMosaic compute function
function visualizeNoiseFreeNeuralResponse(theConeMosaic, theNeuralResponses, temporalSupportSeconds)

    hFig = figure(1999);
    clf;
    set(hFig, 'Position', [10 10 1700 500], 'Color', [1 1 1]);
    axCmosaicActivation = subplot(1,2,1);
    axConeResponseTimeResponses = subplot(1,2,2);

    [nTrials, nTimePoints, nCones] = size(theNeuralResponses);
    activationRange = prctile(theNeuralResponses(:),[5 95]);

    theLconeResponses = theNeuralResponses(:,:,theConeMosaic.lConeIndices);
    theMconeResponses = theNeuralResponses(:,:,theConeMosaic.mConeIndices);
    theSconeResponses = theNeuralResponses(:,:,theConeMosaic.sConeIndices);

    iTrial = 1;
    m = max(squeeze(theLconeResponses(iTrial,:,:)),[], 1);
    [~,idx] = sort(m);
    theLconeResponses = theLconeResponses(iTrial,:,idx);
    m = max(squeeze(theMconeResponses(iTrial,:,:)),[], 1);
    [~,idx] = sort(m);
    theMconeResponses = theMconeResponses(iTrial,:,idx);
    m = max(squeeze(theSconeResponses(iTrial,:,:)),[], 1);
    [~,idx] = sort(m);
    theSconeResponses = theSconeResponses(iTrial,:,idx);

    if (~isempty(theConeMosaic.fixEMobj)) 
        emPathsDegs = theConeMosaic.fixEMobj.emPosArcMin/60;
    else
        emPathsDegs = [];
    end


    for iPoint = 1:nTimePoints
        if (isempty(emPathsDegs))
            theConeMosaic.visualize(...
                'figureHandle', hFig, ...
                'axesHandle', axCmosaicActivation, ...
                'activation', theNeuralResponses(iTrial, iPoint,:), ...
                'activationRange', activationRange, ...
                'plotTitle', sprintf('noise free cone excitations (t: %2.1f msec)', temporalSupportSeconds(iPoint)*1e3));
            
        else
            theConeMosaic.visualize(...
                'figureHandle', hFig, ...
                'axesHandle', axCmosaicActivation, ...
                'activation', theNeuralResponses(iTrial, iPoint,:), ...
                'activationRange', activationRange, ...
                'currentEMposition', squeeze(emPathsDegs(1,iPoint,:)), ...
                'displayedEyeMovementData', struct('trial', 1, 'timePoints', 1:iPoint), ...
                'plotTitle', sprintf('noise free cone excitations (t: %2.1f msec)', temporalSupportSeconds(iPoint)*1e3));
 
        end

        drawnow;
    end

    mosaicSpatioTemporalActivation = [];
    mosaicSpatioTemporalActivation(:, 1:numel(theConeMosaic.lConeIndices)) = squeeze(theLconeResponses(iTrial,:,:));
    mosaicSpatioTemporalActivation(:, size(mosaicSpatioTemporalActivation,2)+(1:numel(theConeMosaic.mConeIndices))) = squeeze(theMconeResponses(iTrial,:,:));
    mosaicSpatioTemporalActivation(:, size(mosaicSpatioTemporalActivation,2)+(1:numel(theConeMosaic.sConeIndices))) = squeeze(theSconeResponses(iTrial,:,:));
    
    dt = 0.5*(temporalSupportSeconds(2)-temporalSupportSeconds(1));
    imagesc(axConeResponseTimeResponses, temporalSupportSeconds*1e3, 1:theConeMosaic.conesNum, mosaicSpatioTemporalActivation');
    hold(axConeResponseTimeResponses, 'on');
    LconeRect.x = [temporalSupportSeconds(1)-dt temporalSupportSeconds(end)+dt temporalSupportSeconds(end)+dt temporalSupportSeconds(1)-dt temporalSupportSeconds(1)-dt]*1e3;
    LconeRect.y = [1  1 numel(theConeMosaic.lConeIndices) numel(theConeMosaic.lConeIndices) 1];
    MconeRect.x = [temporalSupportSeconds(1)-dt temporalSupportSeconds(end)+dt temporalSupportSeconds(end)+dt temporalSupportSeconds(1)-dt temporalSupportSeconds(1)-dt]*1e3;
    MconeRect.y = numel(theConeMosaic.lConeIndices) + [1  1 numel(theConeMosaic.mConeIndices) numel(theConeMosaic.mConeIndices) 1];
    SconeRect.x = [temporalSupportSeconds(1)-dt temporalSupportSeconds(end)+dt temporalSupportSeconds(end)+dt temporalSupportSeconds(1)-dt temporalSupportSeconds(1)-dt]*1e3;
    SconeRect.y = numel(theConeMosaic.lConeIndices) + numel(theConeMosaic.mConeIndices) + [1  1 numel(theConeMosaic.sConeIndices) numel(theConeMosaic.sConeIndices) 1];
    
    plot(axConeResponseTimeResponses, LconeRect.x, LconeRect.y, 'r-', 'LineWidth', 2);
    plot(axConeResponseTimeResponses, MconeRect.x, MconeRect.y, 'g-', 'LineWidth', 2);
    plot(axConeResponseTimeResponses, SconeRect.x, SconeRect.y, 'c-', 'LineWidth', 2);
    for i = 1:numel(temporalSupportSeconds)
        plot(axConeResponseTimeResponses, (temporalSupportSeconds(i)-dt)*[1 1]*1e3, [1 theConeMosaic.conesNum], 'k-', 'LineWidth', 1.0);
    end

    colormap(brewermap(1024, '*greys'));

    hold(axConeResponseTimeResponses, 'off');
    set(axConeResponseTimeResponses, ...
        'XLim', [temporalSupportSeconds(1)-dt temporalSupportSeconds(end)+dt]*1e3, ...
        'YLim', [1 theConeMosaic.conesNum], ...
        'CLim', [activationRange(1) activationRange(2)]);
    axis(axConeResponseTimeResponses, 'xy');
    xlabel(axConeResponseTimeResponses, 'time (msec)');
    ylabel(axConeResponseTimeResponses, 'cone index');
    set(axConeResponseTimeResponses, 'FontSize', 16, 'Color', [1 1 1]);
    drawnow;

end

