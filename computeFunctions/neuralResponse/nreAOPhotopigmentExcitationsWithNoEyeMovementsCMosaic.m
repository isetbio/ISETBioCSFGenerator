function dataOut = nreAOPhotopigmentExcitationsWithNoEyeMovementsCMosaic(...
    neuralEngineOBJ, neuralResponseParamsStruct, sceneSequence, ...
    sceneSequenceTemporalSupport, instancesNum, varargin)
% Compute function for computation of cone excitations witout eye movements
%
% Syntax:
%   dataOut = nreAOPhotopigmentExcitationsWithNoEyeMovements(...
%    neuralEngineOBJ, neuralResponseParamsStruct, sceneSequence, ...
%    sceneSequenceTemporalSupport, instancesNum, varargin);
%
% Description:
%    Function serving as the computeFunctionHandle for a @neuralResponseEngine
%    object.  This is derived from nrePhotopigmentExcitationsWithNoEyeMovements 
%    and describes AO optics.
%
%    There are 2 ways to use this function.
%
%       [1] If called directly and with no arguments,
%           dataOut = nrePhotopigmentExcitationsWithNoEyeMovements()
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
%   'verbose'                       - Boolean. Print out diagnostic stuff?
%                                     Default false.
%
% Outputs:
%    dataOut  - A struct that depends on the input arguments.
%
%               If called directly with no input arguments, the returned struct contains
%               the defaultParams (optics and coneMosaic) that define the neural
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
%                .neuralPipeline  : a struct containing the optics and cone mosaic
%                                   employed in the computation (only returned if
%                                   the parent @neuralResponseEngine object has
%                                   an empty neuralPipeline property)
%
%       The computed neural responses can be extracted as:
%           neuralResponses('one of the entries of noiseFlags')
%       and are arranged in a matrix of:
%           [instancesNum x mCones x tTimeBins]
%
% See Also:
%     t_neuralResponseCompute, t_AODisplay

% History:
%    03/23/2021  npc  Adapted from nreAOPhotopigmentExcitationsWithNoEyeMovements for the @cMosaic object

% Examples:
%{
    % Usage case #1. Just return the default neural response params
    defaultParams = nreAOPhotopigmentExcitationsWithNoEyeMovementsCMosaic()

    % Usage case #2. Compute noise free, noisy, and repeatable (seed: 346) noisy response instances
    % using a parent @neuralResponseEngine object and the default neural response params

    % Instantiate the parent @neuralResponseEngine object
    theNeuralEngineOBJ = neuralResponseEngine(@nreAOPhotopigmentExcitationsWithNoEyeMovementsCMosaic);

    % Instantiate a @sceneEngine object and generate a test scene
    theSceneEngineOBJ = sceneEngine(@sceUniformFieldTemporalModulation);
    testContrast = 0.1;
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

    % Retrieve the different computed responses
    noiseFreeResponses = theResponses('none');
    randomNoiseResponseInstances = theResponses('random');
%}

% Check input arguments. If called with zero input arguments, just return the default params struct
if (nargin == 0)
    dataOut = generateDefaultParams();
    return;
end

% Parse the input arguments
p = inputParser;
p.addParameter('noiseFlags', {'random'});
p.addParameter('rngSeed',[],@(x) (isempty(x) | isnumeric(x)));
p.addParameter('verbose',false,@islogical);
varargin = ieParamFormat(varargin);
p.parse(varargin{:});

% Retrieve the response noiseFlag labels and validate them.
noiseFlags = p.Results.noiseFlags;
rngSeed = p.Results.rngSeed;
neuralEngineOBJ.validateNoiseFlags(noiseFlags);

% For each noise flag we generate a corresponing neural response, and all
% neural responses are stored in a dictionary indexed by the noiseFlag label.
% Setup theNeuralResponses dictionary, loading empty responses for now
theNeuralResponses = containers.Map();
for idx = 1:length(noiseFlags)
    theNeuralResponses(noiseFlags{idx}) = [];
end

if (isempty(neuralEngineOBJ.neuralPipeline))
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
    
    % Compute pupil function and PSF
    %
    % Whether LCA should be included depends on your apparatus and
    % is controlled by the boolean defeatLCA in the computation of
    % the pupil function.
    wvfP = wvfComputePupilFunction(wvfP,false,'no lca',neuralResponseParamsStruct.opticsParams.defeatLCA);
    wvfP = wvfComputePSF(wvfP);
    
    % Generate optical image object from the wavefront object
    theOptics = wvf2oi(wvfP);
    
    % Set the fNumber to correspond to the pupil size
    focalLengthMM = oiGet(theOptics,'focal length')*1000;
    theOptics = oiSet(theOptics, 'optics fnumber', focalLengthMM/neuralResponseParamsStruct.opticsParams.pupilDiameterMM);
    
    % Generate the cone mosaic
    theConeMosaic = cMosaic('sizeDegs', neuralResponseParamsStruct.coneMosaicParams.fovDegs*[1 1], ...
        'integrationTime', neuralResponseParamsStruct.coneMosaicParams.timeIntegrationSeconds);
    
    % Set flag
    returnTheNeuralPipeline = true;
else
    % Load the optics from the previously computed neural pipeline
    theOptics = neuralEngineOBJ.neuralPipeline.optics;
    
    % Load the cone mosaic from the previously computed neural pipeline
    theConeMosaic = neuralEngineOBJ.neuralPipeline.coneMosaic;
    
    % Set flag
    returnTheNeuralPipeline =  false;
end

% Compute the sequence of optical images corresponding to the sequence of scenes
framesNum = numel(sceneSequence);
theListOfOpticalImages = cell(1, framesNum);
for frame = 1:framesNum
    theListOfOpticalImages{frame} = oiCompute(sceneSequence{frame}, theOptics);
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

% Set rng seed if one was passed. Not clear we need to do this because
% all the randomness is in the @coneMosaic compute object, but it
% doesn't hurt to do so, if we ever choose a random number at this
% level.
if (~isempty(rngSeed))
    oldSeed = rng(rngSeed);
end

% Compute responses for each type of noise flag requested
for idx = 1:length(noiseFlags)
    if (contains(ieParamFormat(noiseFlags{idx}), 'none'))
        % Compute the noise-free response
        % To do so, first save the current mosaic noiseFlag
        lastConeMosaicNoiseFlag = theConeMosaic.noiseFlag;
        
        % Set the coneMosaic.noiseFlag to 'none';
        theConeMosaic.noiseFlag = 'none';
        
        % Compute noise-free response instances
        [theNeuralResponses(noiseFlags{idx}), ~, ~, ~, temporalSupportSeconds] = theConeMosaic.compute(theOIsequence.frameAtIndex(1));
        
        
        % Restore the original noise flag
        theConeMosaic.noiseFlag = lastConeMosaicNoiseFlag;
        
    elseif (~isempty(rngSeed))
        % Compute noisy response instances with a specified random noise seed for repeatability
        [~, theNeuralResponses(noiseFlags{idx}), ~, ~, temporalSupportSeconds] = theConeMosaic.compute(theOIsequence.frameAtIndex(1), ...
          'nTrials', instancesNum, 'seed', rngSeed);

    elseif (contains(ieParamFormat(noiseFlags{idx}), 'random'))
        % Because computeForOISequence freezes noise, if we want
        % unfrozen noise (which is the case if we are here),
        % we have to pass it a randomly chosen seed.
        useSeed = randi(32000,1,1);
        
        % Compute noisy response instances
        [~, theNeuralResponses(noiseFlags{idx}), ~, ~, temporalSupportSeconds] = theConeMosaic.compute(theOIsequence.frameAtIndex(1), ...
          'nTrials', instancesNum, 'seed', useSeed);
      
    end
end

% Restore rng seed if we set it
if (~isempty(rngSeed))
    rng(oldSeed);
end

% Assemble the dataOut struct
dataOut = struct(...
    'neuralResponses', theNeuralResponses, ...
    'temporalSupport', temporalSupportSeconds);
if (returnTheNeuralPipeline)
    dataOut.neuralPipeline.optics = theOptics;
    dataOut.neuralPipeline.coneMosaic = theConeMosaic;
end
end

function p = generateDefaultParams()
% Default params for this compute function
p = struct(...
    'opticsParams', struct(...
    'wls', 400:10:700, ...
    'type', 'AO', ...
    'pupilDiameterMM', 6.0, ...
    'defocusAmount', 0.1, ...
    'accommodatedWl', 550, ...
    'zCoeffs', zeros(66,1), ...
    'defeatLCA', false ... 
    ), ...
    'coneMosaicParams', struct(...
    'fovDegs', 0.3, ...
    'timeIntegrationSeconds', 5/1000 ...
    ) ...
    );
end
