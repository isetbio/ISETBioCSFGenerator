function dataOut = nreNoiseFreeAOPhotopigmentExcitationsCMosaic(...
    neuralEngine, neuralResponseParamsStruct, sceneSequence, ...
    sceneSequenceTemporalSupport, varargin)
% Compute function for computation of cone excitations witout eye movements
%
% Syntax:
%   dataOut = nreNoiseFreeAOPhotopigmentExcitations(...
%    neuralEngineOBJ, neuralResponseParamsStruct, sceneSequence, ...
%    sceneSequenceTemporalSupport, varargin);
%
% Description:
%    Function serving as the computeFunctionHandle for a @neuralResponseEngine
%    object.  This is derived from nreNoiseFreePhotopigmentExcitationsCmosaic
%    and is elaborated to describe AO optics.
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
%                                     This should just be empty for this
%                                     compute function.
%    sceneSequence                  - a cell array of scenes defining the frames of a stimulus
%    sceneSequenceTemporalSupport   - the temporal support for the stimulus frames, in seconds
%
% Optional key/value input arguments:
%   'verbose'                       - Boolean, default false. Print out
%                                     diagnostic information.
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
% See Also:
%     nreNoiseFreePhotopigmentExcitationsCmosaic, t_neuralResponseCompute,
%     t_AODisplay, t_BerkeleyAOTumblingEThreshold

% History:
%    03/23/2021  npc  Adapted from nreAOPhotopigmentExcitationsWithNoEyeMovements for the @cMosaic object
%    12/21/2024  dhb  Upated for new architecture

% Check input arguments. If called with zero input arguments, just return the default params struct
if (nargin == 0)
    dataOut = generateDefaultParams();
    return;
end

% Parse the input arguments
p = inputParser;
p.KeepUnmatched = true;
p.addParameter('verbose',false,@islogical);
varargin = ieParamFormat(varargin);
p.parse(varargin{:});

if (isempty(neuralEngine.neuralPipeline) | ~isfield(neuralEngine.neuralPipeline,'noiseFreeResponse'))
    % Generate the optics
    %
    % Some checks
    if (~strcmp(neuralResponseParamsStruct.opticsParams.type,'AO'))
        error('Expecting neuralResponseParamsStruct.opticsParams to be ''AO''');
    end
    
    % Set up wavefront optics object
    %
    % Compute pupil function using 'no lca' key/value pair to turn off LCA.
    % You can turn it back on to compare the effect.
    %
    % Deal with best focus by specifying that the wavefront parameters
    % were measured at the wavelength we want to say is in focus. This
    % is a little bit of a hack but seems OK for the diffraction limited case
    % we're using here.
    wvfP = wvfCreate('calc wavelengths', neuralResponseParamsStruct.opticsParams.wls, ...
        'zcoeffs', neuralResponseParamsStruct.opticsParams.zCoeffs, ...
        'name', sprintf('humanAO-%d', neuralResponseParamsStruct.opticsParams.pupilDiameterMM));
    wvfP = wvfSet(wvfP, 'measured wavelength', neuralResponseParamsStruct.opticsParams.accommodatedWl);
    wvfP = wvfSet(wvfP, 'measured pupil size', neuralResponseParamsStruct.opticsParams.pupilDiameterMM);
    wvfP = wvfSet(wvfP, 'calc pupil size', neuralResponseParamsStruct.opticsParams.pupilDiameterMM);
    wvfP = wvfSet(wvfP,'zcoeffs', neuralResponseParamsStruct.opticsParams.defocusAmount, 'defocus');
    if (~neuralResponseParamsStruct.opticsParams.defeatLCA)
        wvfP = wvfSet(wvfP,'lcaMethod','human');
    end
    
    % Compute pupil function and PSF
    %
    % Whether LCA should be included depends on your apparatus and
    % is controlled by the boolean defeatLCA in the computation of
    % the pupil function.
    wvfP = wvfCompute(wvfP);
    
    % Generate optical image object from the wavefront object
    theOptics = wvf2oi(wvfP,'humanlens',true);
    
    % Set the fNumber to correspond to the pupil size
    focalLengthMM = oiGet(theOptics,'focal length')*1000;
    theOptics = oiSet(theOptics, 'optics fnumber', focalLengthMM/neuralResponseParamsStruct.opticsParams.pupilDiameterMM);
    
    % Generate the cone mosaic
    theConeMosaic = cMosaic('wave', neuralResponseParamsStruct.coneMosaicParams.wave, ...
        'sizeDegs', neuralResponseParamsStruct.coneMosaicParams.fovDegs*[1 1], ...
        'integrationTime', neuralResponseParamsStruct.coneMosaicParams.timeIntegrationSeconds, ...
        'noiseFlag', 'none');
   
    % Set flag
    returnTheNeuralPipeline = true;
else
    % Load the optics from the previously computed neural pipeline
    theOptics = neuralEngine.neuralPipeline.optics;
    
    % Load the cone mosaic from the previously computed neural pipeline
    theConeMosaic = neuralEngine.neuralPipeline.coneMosaic;
    
    % Set flag
    returnTheNeuralPipeline =  false;
end

% Compute the sequence of optical images corresponding to the sequence of scenes
framesNum = numel(sceneSequence);
theListOfOpticalImages = cell(1, framesNum);
for frame = 1:framesNum
    theListOfOpticalImages{frame} = oiCompute(theOptics,sceneSequence{frame});
end

% Generate an @oiSequence object containing the list of computed optical images
theOIsequence = oiArbitrarySequence(theListOfOpticalImages, sceneSequenceTemporalSupport);

% Some diagnosis
if (p.Results.verbose)
    fprintf('Optical image aperuture diameter = %g mm\n',opticsGet(oiGet(theListOfOpticalImages{1},'optics'),'aperture diameter')*1000);
    fprintf('Optical image focal length = %g mm\n',opticsGet(oiGet(theListOfOpticalImages{1},'optics'),'focal length')*1000);
    theWl = 400;
    theFrame = 1;
    index = find(theListOfOpticalImages{theFrame}.spectrum.wave == theWl);
    temp = theListOfOpticalImages{theFrame}.data.photons(:,:,index);
    fprintf('At %d nm, frame %d, oi mean, min, max: %g, %g, %g\n',theWl,theFrame,mean(temp(:)),min(temp(:)),max(temp(:)));
    theWl = 550;
    index = find(theListOfOpticalImages{theFrame}.spectrum.wave == theWl);
    temp = theListOfOpticalImages{theFrame}.data.photons(:,:,index);
    fprintf('At %d nm, frame %d, oi mean, min, max: %g, %g, %g\n',theWl,theFrame,mean(temp(:)),min(temp(:)),max(temp(:)));
end

% Compute responses for each type of noise flag requested
[theNeuralResponsesRaw, ~, ~, ~, temporalSupportSeconds] = ...
    theConeMosaic.compute(theOIsequence, ...
    'nTrials', 1, ...
    'withFixationalEyeMovements', false ...
    );

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

% % Compute responses for each type of noise flag requested
% for idx = 1:length(noiseFlags)
%     if (contains(ieParamFormat(noiseFlags{idx}), 'none'))
%         % Compute the noise-free response
%         % To do so, first save the current mosaic noiseFlag
%         lastConeMosaicNoiseFlag = theConeMosaic.noiseFlag;
% 
%         % Set the coneMosaic.noiseFlag to 'none';
%         theConeMosaic.noiseFlag = 'none';
% 
%         % Compute noise-free response instances
%         %[theNeuralResponses(noiseFlags{idx}), ~, ~, ~, temporalSupportSeconds] = theConeMosaic.compute(theOIsequence.frameAtIndex(1));
%         [theNeuralResponses(noiseFlags{idx}), ~, ~, ~, temporalSupportSeconds] = theConeMosaic.compute(theOIsequence);
% 
%         % Restore the original noise flag
%         theConeMosaic.noiseFlag = lastConeMosaicNoiseFlag;
% 
%     elseif (~isempty(rngSeed))
%         % Compute noisy response instances with a specified random noise seed for repeatability
%         % [~, theNeuralResponses(noiseFlags{idx}), ~, ~, temporalSupportSeconds] = theConeMosaic.compute(theOIsequence.frameAtIndex(1), ...
%         %   'nTrials', instancesNum, 'seed', rngSeed);
%         [~, theNeuralResponses(noiseFlags{idx}), ~, ~, temporalSupportSeconds] = theConeMosaic.compute(theOIsequence, ...
%           'nTrials', instancesNum, 'seed', rngSeed);
% 
%     elseif (contains(ieParamFormat(noiseFlags{idx}), 'random'))
%         % Because computeForOISequence freezes noise, if we want
%         % unfrozen noise (which is the case if we are here),
%         % we have to pass it a randomly chosen seed.
%         useSeed = randi(32000,1,1);
% 
%         % Compute noisy response instances
%         % [~, theNeuralResponses(noiseFlags{idx}), ~, ~, temporalSupportSeconds] = theConeMosaic.compute(theOIsequence.frameAtIndex(1), ...
%         %   'nTrials', instancesNum, 'seed', useSeed);
%         [~, theNeuralResponses(noiseFlags{idx}), ~, ~, temporalSupportSeconds] = theConeMosaic.compute(theOIsequence, ...
%           'nTrials', instancesNum, 'seed', useSeed);
%     end
% end
% 
% % Restore rng seed if we set it
% if (~isempty(rngSeed))
%     rng(oldSeed);
% end
% 
% % Assemble the dataOut struct
% dataOut = struct(...
%     'neuralResponses', theNeuralResponses, ...
%     'temporalSupport', temporalSupportSeconds);
% if (returnTheNeuralPipeline)
%     dataOut.neuralPipeline.optics = theOptics;
%     dataOut.neuralPipeline.coneMosaic = theConeMosaic;
% end
end

function p = generateDefaultParams()
% Default params for this compute function
p = struct(...
    'opticsParams', ...
        struct(...
        'wls', 400:10:700, ...
        'type', 'AO', ...
        'pupilDiameterMM', 6.0, ...
        'defocusAmount', 0.1, ...
        'accommodatedWl', 550, ...
        'zCoeffs', zeros(66,1), ...
        'defeatLCA', false ...
        ), ...
    'coneMosaicParams', ...
        struct(...
        'wave', 400:10:700, ...
        'fovDegs', 0.3, ...
        'timeIntegrationSeconds', 5/1000 ...
        ) ...
    );
end
