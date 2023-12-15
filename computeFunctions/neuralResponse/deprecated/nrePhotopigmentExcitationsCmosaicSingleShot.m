function dataOut = nrePhotopigmentExcitationsCmosaicSingleShot(...
    neuralEngineOBJ, neuralResponseParamsStruct, sceneSequence, ...
    sceneSequenceTemporalSupport, instancesNum, varargin)
% Calculates cone excitations without considering fixational eye movements.
% This function is designed for scenes with a single frame only. It has
% been deprecated, so when called, it automatically redirects to a more
% general function, nrePhotopigmentExcitationsCmosaic.m.
%
% Syntax:
%   dataOut = nrePhotopigmentExcitationsConeMosaicSingleShot(...
%    neuralEngineOBJ, neuralResponseParamsStruct, sceneSequence, ...
%    sceneSequenceTemporalSupport, instancesNum, varargin);
%    
% See:
%    nrePhotopigmentExcitationsCmosaic.m

% History:
%    03/29/2021  npc  Wrote it by adapting nrePhotopigmentExcitationsConeMosaicHexWithNoEyeMovements
%    09/28/2023  fh   Move this function to the deprecated file.  

% Check input arguments. If called with zero input arguments, just return 
% the default params struct
if (nargin == 0)
    warning(['This function has been deprecated. Consider using a more',...
        ' general function nrePhotopigmentExcitationsCmosaic.m']);
    dataOut = nrePhotopigmentExcitationsCmosaic();
    return;
end

%call the master function 
%the total number of frames
total_seq = length(sceneSequence);
if total_seq ~=1
    error(['The input scene is not a single shot! Consider ',...
        'calling nrePhotopigmentExcitationsCmosaic.m']);
else
    dataOut = nrePhotopigmentExcitationsCmosaic(neuralEngineOBJ,...
        neuralResponseParamsStruct, sceneSequence, ...
        sceneSequenceTemporalSupport, instancesNum, varargin{:});
end

    % %the total number of frames
    % total_seq = length(sceneSequence);
    % %if the scene is not single shot, then send a warning
    % if total_seq ~=1
    %     clear str_callingMasterFunc
    %     error(['The input scene is not a single shot! Consider ',...
    %         'calling nrePhotopigmentExcitationsCmosaic.m']);
    % end
    % 
    % % Parse the input arguments
    % p = inputParser;
    % p.addParameter('noiseFlags', {'random'});
    % p.addParameter('rngSeed',[],@(x) (isempty(x) | isnumeric(x)));
    % varargin = ieParamFormat(varargin);
    % p.parse(varargin{:});
    % 
    % % Retrieve the response noiseFlag labels and validate them.
    % noiseFlags = p.Results.noiseFlags;
    % rngSeed = p.Results.rngSeed;
    % neuralEngineOBJ.validateNoiseFlags(noiseFlags);
    % 
    % % For each noise flag we generate a corresponing neural response, and all 
    % % neural responses are stored in a dictionary indexed by the noiseFlag label.
    % % Setup theNeuralResponses dictionary, loading empty responses for now
    % theNeuralResponses = containers.Map();
    % for idx = 1:length(noiseFlags)
    %     theNeuralResponses(noiseFlags{idx}) = [];
    % end
    % 
    % if (isempty(neuralEngineOBJ.neuralPipeline))
    %     % Generate the @cMosaic object
    %     theConeMosaic = cMosaic(...
    %         'sizeDegs', neuralResponseParamsStruct.coneMosaicParams.sizeDegs, ...
    %         'eccentricityDegs', neuralResponseParamsStruct.coneMosaicParams.eccDegs, ...
    %         'integrationTime', neuralResponseParamsStruct.coneMosaicParams.timeIntegrationSeconds ...
    %         );
    % 
    %     % Generate optics appropriate for the mosaic's eccentricity
    %     oiEnsemble = theConeMosaic.oiEnsembleGenerate(neuralResponseParamsStruct.coneMosaicParams.eccDegs, ...
    %         'zernikeDataBase', 'Polans2015', ...
    %         'subjectID', neuralResponseParamsStruct.opticsParams.PolansSubject, ...
    %         'pupilDiameterMM', neuralResponseParamsStruct.opticsParams.pupilDiameterMM);
    %     theOptics = oiEnsemble{1};
    %     returnTheNeuralPipeline = true;
    % else
    %     % Load the optics from the previously computed neural pipeline
    %     theOptics = neuralEngineOBJ.neuralPipeline.optics;
    %     % Load the cone mosaic from the previously computed neural pipeline
    %     theConeMosaic = neuralEngineOBJ.neuralPipeline.coneMosaic;
    %     returnTheNeuralPipeline =  false;
    % end
    % 
    % % Compute the sequence of optical images corresponding to the sequence of scenes
    % theOIsequence = oiCompute(theOptics, sceneSequence{1});
    % 
    % % Set rng seed if one was passed. Not clear we need to do this because
    % % all the randomness is in the @coneMosaic compute object, but it
    % % doesn't hurt to do so, if we ever choose a random number at this
    % % level.
    % if (~isempty(rngSeed))
    %     oldSeed = rng(rngSeed);
    % end
    % 
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
    %         [theNeuralResponses(noiseFlags{idx}), ~, ~, ~, temporalSupportSeconds] = ...
    %             theConeMosaic.compute(theOIsequence, ...
    %             'nTrials', instancesNum ...
    %         );
    % 
    %         % Restore the original noise flag
    %         theConeMosaic.noiseFlag = lastConeMosaicNoiseFlag;
    % 
    %     elseif (~isempty(rngSeed))
    %         % Compute noisy response instances with a specified random noise seed for repeatability
    %         [~,theNeuralResponses(noiseFlags{idx}), ~, ~, temporalSupportSeconds] = ...
    %             theConeMosaic.compute(theOIsequence, ...
    %             'nTrials', instancesNum, ...
    %             'seed', rngSeed ...        % random seed
    %         );
    % 
    %     elseif (contains(ieParamFormat(noiseFlags{idx}), 'random'))
    %         % Because computeForOISequence freezes noise, if we want
    %         % unfrozen noise (which is the case if we are here), 
    %         % we have to pass it a randomly chosen seed.
    %         useSeed = randi(32000,1,1);
    % 
    %         % Compute noisy response instances
    %         [~,theNeuralResponses(noiseFlags{idx}), ~, ~, temporalSupportSeconds] = ...
    %             theConeMosaic.compute(theOIsequence, ...
    %             'nTrials', instancesNum, ...
    %             'seed', useSeed ...        % random seed
    %         );
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
    % 
    % if (returnTheNeuralPipeline)
    %     dataOut.neuralPipeline.optics = theOptics;
    %     dataOut.neuralPipeline.coneMosaic = theConeMosaic;
    % end
end

% function p = generateDefaultParams()
%     % Default params for this compute function
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
