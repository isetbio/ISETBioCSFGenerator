function dataOut = nreNoiseFreePhotopigmentExcitationsCmosaic(...
    neuralEngine, noiseFreeComputeParams, sceneSequence, ...
    sceneSequenceTemporalSupport, varargin)
% Compute function for computation of cone excitations
%
% Syntax:
%   dataOut = nreNoiseFreePhotopigmentExcitationsCmosaic(...
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
%       it does not compute anything and simply returns a struct with the 
%       defaultParams (optics and coneMosaic params) that define the neural 
%       compute pipeline for this computation.
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
%   'fixationalEM'           - Empty (default) or a fixationalEM object
%                              that describes one eye movement path.  If
%                              the latter, this must have one position per
%                              frame of the passed scene sequence, in which
%                              case it is applied.
% 
% See Also:
%     t_neuralResponseCompute

% History:
%    03/29/2021  npc  Wrote it by adapting nrePhotopigmentExcitationsConeMosaicHexWithNoEyeMovements
%    04/10/2023  fh   Edited it so that with correctly set varargin, this
%                       function behaves the same as
%                       nrePhotopigmentExcitationsCmosaicWithNoEyeMovements.m
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
    nreParams = nreNoiseFreePhotopigmentExcitationsCmosaic()
    nreParams.coneMosaicParams.timeIntegrationSeconds = ...
        sceneTemporalSupportSeconds(2)-sceneTemporalSupportSeconds(1);
    
    % Usage case #2. Compute noise free responses 
    % using a parent @neuralResponseEngine object and the default neural response params
    %
    % Instantiate the parent @neuralResponseEngine object, with Poisson noise
    theNeuralEngine = neuralResponseEngine(@nreNoiseFreePhotopigmentExcitationsCmosaic, ...
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

    % Check input arguments. If called with zero input arguments, just return the default params struct
    if (nargin == 0)
        dataOut = generateDefaultParams();
        return;
    end
    
    % Parse the input arguments.
    %
    % Allow for possibility that other nre's take key/value pairs that we
    % can ignore, so set KeepUnmatched to true.
    p = inputParser;
    p.KeepUnmatched = true;
    p.addParameter('fixationalEM', [], @(x)(isempty(x) || (isa(x,'fixationalEM'))));
    varargin = ieParamFormat(varargin);
    p.parse(varargin{:});
    fixationalEMObj = p.Results.fixationalEM;

    % Get the number of scene sequence frames
    framesNum = numel(sceneSequence);
    
    % Create/get key objects
    if (isempty(neuralEngine.neuralPipeline) | ~isfield(neuralEngine.neuralPipeline,'noiseFreeResponse'))
        % Generate the @cMosaic object
        theConeMosaic = cMosaic(...
            'wave', noiseFreeComputeParams.coneMosaicParams.wave, ...
            'sizeDegs', noiseFreeComputeParams.coneMosaicParams.sizeDegs, ...
            'eccentricityDegs', noiseFreeComputeParams.coneMosaicParams.eccDegs, ...
            'integrationTime', noiseFreeComputeParams.coneMosaicParams.timeIntegrationSeconds, ...
            'noiseFlag', 'none' ...
            );

        % Generate optics appropriate for the mosaic's eccentricity
        oiEnsemble = theConeMosaic.oiEnsembleGenerate(noiseFreeComputeParams.coneMosaicParams.eccDegs, ...
            'zernikeDataBase', 'Polans2015', ...
            'subjectID', noiseFreeComputeParams.opticsParams.PolansSubject, ...
            'pupilDiameterMM', noiseFreeComputeParams.opticsParams.pupilDiameterMM);
        theOptics = oiEnsemble{1};
        clear oiEnsemble
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
        theOIsequence = oiCompute(theOptics, sceneSequence{1},'padvalue','mean');
    else
        theListOfOpticalImages = cell(1, framesNum);
        for frame = 1:framesNum
            theListOfOpticalImages{frame} = oiCompute(theOptics, sceneSequence{frame},'padvalue','mean');
        end

        % Generate an @oiSequence object containing the list of computed optical images
        theOIsequence = oiArbitrarySequence(theListOfOpticalImages, sceneSequenceTemporalSupport);
    end

    % Compute responses depending on whether there is a fixationalEM or
    % not.
    if (isempty(fixationalEMObj))
        % No fixationalEM passed. So don't move the eyes.
        [theNeuralResponsesRaw, ~, ~, ~, temporalSupportSeconds] = ...
            theConeMosaic.compute(theOIsequence, ...
            'nTrials', 1, ...
            'withFixationalEyeMovements', false ...
            );
    else
        % We were passed a fixational EM.  Check that it is OK for our
        % purposes, and then comptute with it.
        %
        % Check

        % Compute
        theConeMosaic.emSetFixationalEMObj(fixationalEMObj);
        [theNeuralResponsesRaw, ~, ~, ~, temporalSupportSeconds] = ...
            theConeMosaic.compute(theOIsequence, ...
            'nTrials', 1, ...
            'withFixationalEyeMovements', true ...
            );
    end

    % Convert cMosaic return convention to nre return convention.  We have
    % to deal with unfortunate special casing because of the way Matlab
    % and/or cMosaic handle singleton dimensions, so that the returned
    % responses always have the responses in the column(s).
    theNeuralResponses(:,:) = theNeuralResponsesRaw(1,:,:);
    clear theNeuralResponsesRaw
    if (length(temporalSupportSeconds) > 1)
           theNeuralResponses = theNeuralResponses';
    end

    % Assemble the dataOut struct
    dataOut = struct(...
        'neuralResponses', theNeuralResponses, ...
    	'temporalSupport', temporalSupportSeconds);
    
    if (returnTheNoiseFreePipeline)
        dataOut.noiseFreeResponsePipeline.optics = theOptics;
        dataOut.noiseFreeResponsePipeline.coneMosaic = theConeMosaic;
    end
end

function p = generateDefaultParams()
    % Default params for this compute function
    p = struct(...
        'opticsParams', struct(...
            'PolansSubject', 10, ...
            'pupilDiameterMM', 3.0 ...
        ), ...
        'coneMosaicParams', struct(...
        'wave', 400:10:700, ...
            'sizeDegs', 0.3*[1 1], ...
            'eccDegs', [0 0], ...
            'timeIntegrationSeconds', 5/1000 ...
        ) ...
    );
end
