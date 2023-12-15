function dataOut = nrePhotopigmentExcitationsCmosaicWithNoEyeMovements(...
    neuralEngineOBJ, neuralResponseParamsStruct, sceneSequence, ...
    sceneSequenceTemporalSupport, instancesNum, varargin)
% Calculates cone excitations without considering fixational eye movements.
% This function allows an input of a sequence of scene, but it only
% computes cone excitations given the first scene and neglects the rest. 
% Since this is a special case, it has been deprecated, and when it's 
% called, it automatically redirects to a more general function 
% nrePhotopigmentExcitationsCmosaic.m.
%
% Syntax:
%   dataOut = nrePhotopigmentExcitationsCmosaicWithNoEyeMovements(...
%    neuralEngineOBJ, neuralResponseParamsStruct, sceneSequence, ...
%    sceneSequenceTemporalSupport, instancesNum, varargin);
%
% See: 
%   nrePhotopigmentExcitationsCmosaic.m
%
% History:
%    09/26/2020  npc  Wrote it.
%    10/05/2020  dhb  Apply ieParamFormat to varargin for all keys.
%    10/05/2020  dhb  Rename. Work on comments.
%    10/05/2020  dhb  Rewrite to use 'rngSeed' key/value pair.
%    10/17/2020  dhb  Use randomly chosen seed for mosaic compute operation
%                     if rngSeed is set to [].  Save/restore rng state when
%                     an explicit seed is passed.
%                dhb  Just return one response instance in the no noise
%                     case.
%    10/19/2020  dhb  Fix comment to reflect fact that we now return
%                     instancesNum instances in noise free case.
%    03/14/2023  mw   Switch from using coneMosaicHex to cMosaic. Add
%                     control of eccentricity.
%    09/28/2023  fh   Moved this function to the deprecated file.  
%

% Check input arguments. If called with zero input arguments, just return 
% the default params struct
if (nargin == 0)
    warning(['This function has been deprecated. Consider using a more ',...
        'general function nrePhotopigmentExcitationsCmosaic.m. If there are ',...
        'more than one scene, and you would like to keep only the first scene,',...
        ' like how this deprecated function does, you could set ',...
        '''amputateScenes'' to true when calling the general function.']);
    dataOut = nrePhotopigmentExcitationsCmosaic();
    return;
end

%if there is only one scene, then we do not need to amputate
if length(sceneSequence) == 1; amputateScenes = false;
else; amputateScenes = true; end
%add the extra key value pair
varargin_extended = [varargin, 'amputateScenes', amputateScenes];
%redirect to the general function
dataOut = nrePhotopigmentExcitationsCmosaic(neuralEngineOBJ,...
    neuralResponseParamsStruct, sceneSequence,  ...
    sceneSequenceTemporalSupport, instancesNum, varargin_extended{:});
        

%     if (nargin == 0)
%         dataOut = generateDefaultParams();
%         return;
%     end
% 
%     % Parse the input arguments
%     p = inputParser;
%     p.addParameter('noiseFlags', {'random'});
%     p.addParameter('rngSeed',[],@(x) (isempty(x) | isnumeric(x)));
%     varargin = ieParamFormat(varargin);
%     p.parse(varargin{:});
% 
%     % Retrieve the response noiseFlag labels and validate them.
%     noiseFlags = p.Results.noiseFlags;
%     rngSeed = p.Results.rngSeed;
%     neuralEngineOBJ.validateNoiseFlags(noiseFlags);
% 
%     % For each noise flag we generate a corresponing neural response, and all 
%     % neural responses are stored in a dictionary indexed by the noiseFlag label.
%     % Setup theNeuralResponses dictionary, loading empty responses for now
%     theNeuralResponses = containers.Map();
%     for idx = 1:length(noiseFlags)
%         theNeuralResponses(noiseFlags{idx}) = [];
%     end
% 
%     if (isempty(neuralEngineOBJ.neuralPipeline))
%         % Generate the optics
% %         theOptics = oiCreate(neuralResponseParamsStruct.opticsParams.type, neuralResponseParamsStruct.opticsParams.pupilDiameterMM);  
% 
%         % Generate the cone mosaic
%         % Modified by Mengxin 2023/03/14: switched from using coneMosaicHex
%         % to cMosaic; eccentricity added
%         theConeMosaic = cMosaic(...
%             'sizeDegs', neuralResponseParamsStruct.coneMosaicParams.sizeDegs, ... 
%             'eccentricityDegs', neuralResponseParamsStruct.coneMosaicParams.eccDegs, ...
%             'integrationTime', neuralResponseParamsStruct.coneMosaicParams.timeIntegrationSeconds ...
%             );
% 
%         % Generate optics appropriate for the mosaic's eccentricity
%         oiEnsemble = theConeMosaic.oiEnsembleGenerate(neuralResponseParamsStruct.coneMosaicParams.eccDegs, ...
%             'zernikeDataBase', 'Polans2015', ...
%             'subjectID', neuralResponseParamsStruct.opticsParams.PolansSubject, ...
%             'pupilDiameterMM', neuralResponseParamsStruct.opticsParams.pupilDiameterMM);
%         theOptics = oiEnsemble{1};
%         returnTheNeuralPipeline = true;
%     else
%         % Load the optics from the previously computed neural pipeline
%         theOptics = neuralEngineOBJ.neuralPipeline.optics;
%         % Load the cone mosaic from the previously computed neural pipeline
%         theConeMosaic = neuralEngineOBJ.neuralPipeline.coneMosaic;
%         returnTheNeuralPipeline =  false;
%     end
% 
%     % Compute the sequence of optical images corresponding to the sequence of scenes
%     framesNum = numel(sceneSequence);
%     theListOfOpticalImages = cell(1, framesNum);
%     for frame = 1:framesNum
%         theListOfOpticalImages{frame} = oiCompute(theOptics, sceneSequence{frame});
%     end
% 
%     % Generate an @oiSequence object containing the list of computed optical images
%     theOIsequence = oiArbitrarySequence(theListOfOpticalImages, sceneSequenceTemporalSupport);
% 
%     % Zero eye movements
%     eyeMovementsNum = theOIsequence.maxEyeMovementsNumGivenIntegrationTime(theConeMosaic.integrationTime);
%     emPaths = zeros(instancesNum, eyeMovementsNum, 2);
% 
%     % Set rng seed if one was passed. Not clear we need to do this because
%     % all the randomness is in the @coneMosaic compute object, but it
%     % doesn't hurt to do so, if we ever choose a random number at this
%     % level.
%     if (~isempty(rngSeed))
%         oldSeed = rng(rngSeed);
%     end
% 
%     timeSamplesNum = eyeMovementsNum;
% 
%     % Compute responses for each type of noise flag requested
%     % Modified by Mengxin 2023/03/14: from computeForOISequence to compute
%     for idx = 1:length(noiseFlags)
%         if (contains(ieParamFormat(noiseFlags{idx}), 'none'))
%             % Compute the noise-free response
%             % To do so, first save the current mosaic noiseFlag
%             lastConeMosaicNoiseFlag = theConeMosaic.noiseFlag;
% 
%             % Set the coneMosaic.noiseFlag to 'none';
%             theConeMosaic.noiseFlag = 'none';
% 
%             % Compute noise-free response instances
%             [theNeuralResponses(noiseFlags{idx}), ~, ~, ~, temporalSupportSeconds] = ...
%                 theConeMosaic.compute(theOIsequence.frameAtIndex(1), ...
%                 'nTrials', instancesNum, ...
%                 'nTimePoints', timeSamplesNum ...
%             );
% 
%             % Restore the original noise flag
%             theConeMosaic.noiseFlag = lastConeMosaicNoiseFlag;
% 
%         elseif (~isempty(rngSeed))
%             % Compute noisy response instances with a specified random noise seed for repeatability
%             [~,theNeuralResponses(noiseFlags{idx}), ~, ~, temporalSupportSeconds] = ...
%                 theConeMosaic.compute(theOIsequence.frameAtIndex(1), ...
%                 'nTrials', instancesNum, ...
%                 'nTimePoints', timeSamplesNum, ...
%                 'withFixationalEyeMovements', flase, ...
%                 'seed', rngSeed ...        % random seed
%             );
% 
%         elseif (contains(ieParamFormat(noiseFlags{idx}), 'random'))
%             % Because computeForOISequence freezes noise, if we want
%             % unfrozen noise (which is the case if we are here), 
%             % we have to pass it a randomly chosen seed.
%             useSeed = randi(32000,1,1);
% 
%             % Compute noisy response instances
%             [~, theNeuralResponses(noiseFlags{idx}), ~, ~, temporalSupportSeconds] = ...
%                 theConeMosaic.compute(theOIsequence.frameAtIndex(1), ...
%                 'nTrials', instancesNum, ...
%                 'nTimePoints', timeSamplesNum, ...
%                 'withFixationalEyeMovements', false, ...
%                 'seed', useSeed ...        % random seed
%             );
%         end
%     end
% 
%     % Restore rng seed if we set it
%     if (~isempty(rngSeed))
%         rng(oldSeed);
%     end
% 
%     % Temporal support for the neural response
% %     temporalSupportSeconds = theConeMosaic.timeAxis; 
% 
%     % Assemble the dataOut struct
%     dataOut = struct(...
%         'neuralResponses', theNeuralResponses, ...
% 	    'temporalSupport', temporalSupportSeconds);
%     if (returnTheNeuralPipeline)
%         dataOut.neuralPipeline.optics = theOptics;
%         dataOut.neuralPipeline.coneMosaic = theConeMosaic;
%     end
% end
% 
% function p = generateDefaultParams()
%     % Default params for this compute function
%     % Modified by Mengxin 2023/03/14 to match the parameters for cMosaic
%     p = struct(...
%         'opticsParams', struct(...
%             'PolansSubject', 10, ...
%             'pupilDiameterMM', 3.0 ...
%         ), ...
%         'coneMosaicParams', struct(...
%             'sizeDegs', 0.3*[1 1], ...
%             'eccDegs', [0 0], ...
%             'timeIntegrationSeconds', 5/1000 ...
%         ) ...
%     );
% end
