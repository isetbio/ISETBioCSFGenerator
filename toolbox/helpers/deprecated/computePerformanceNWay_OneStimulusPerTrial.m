function [predictions, theClassifierEngine, responses] = computePerformanceNWay_OneStimulusPerTrial(theScenes, ...
    temporalSupport, nTrain, nTest, theNeuralEngine, theClassifierEngine, trainNoiseFlag, testNoiseFlag, ...
    saveResponses, visualizeAllComponents)
% Compute performance of a classifier given multiple scenes, one for each 
% alternative choice, a neural engine, and a classifier engine. This 
% function should be called when the task is N-alternative forced-choice 
% given one stimulus per trial. This function has been merged with
% computePerformanceTAFC.m, and when it's called, it automatically
% redirects to a more general function, computePerformance.m. 
%
% Syntax:
%    [predictions, theClassifierEngine, responses] = ...
%        computePerformanceNWay_OneStimPerTrial(theScenes, temporalSupport, nTrain, nTest, theNeuralEngine, ...
%           theClassifierEngine, trainNoiseFlag, testNoiseFlag, saveResponses, visualizeAllComponents)
%
% See: 
%   computePerformance.m
%
% History:
%   10/23/20  dhb  Comments.
%   05/10/23  fh   Added a warning and moved the file to the deprecated
%                   folder

warning(['This function has been deprecated. Please consider using a more ',...
    'general function computePerformance.m!']);

%direct to the general function
[predictions, theClassifierEngine, responses] = computePerformance(theScenes, ...
    temporalSupport, nTrain, nTest, theNeuralEngine, theClassifierEngine, ...
    trainNoiseFlag, testNoiseFlag, 'TAFC', false, 'saveResponses', saveResponses,...
    'visualizeAllComponents',visualizeAllComponents);


% % Empty responses
% responses = [];
% 
% % Number of alternatives
% nAlternatives = length(theScenes);
% 
% % Train the classifier.
% %
% % If trainNoiseFlag is empty, then the passed classifier has already been trained
% % and training is skipped.  Otherwise trainFlag is passed to the stimulus
% % generation routine to indicate what type of noise (typically 'none' or
% % 'random') should be used in the training.
% %
% % Note use of combineContainers to reformat the individual responses the
% % way we need them.
% if (~isempty(trainNoiseFlag))
%     inSampleStimResponsesCell = cell(1,nAlternatives);
%     % Generate stimuli for training
%     for aa = 1:nAlternatives
%         [inSampleStimResponsesCell{aa}, ~] = theNeuralEngine.compute(...
%             theScenes{aa}, ...
%             temporalSupport, ...
%             nTrain, ...
%             'noiseFlags', {trainNoiseFlag});
%     end
%     inSampleStimResponses = combineContainers(inSampleStimResponsesCell);
% 
%     % Visualization.
%     if visualizeAllComponents
%         visualizeConeResps(theNeuralEngine, inSampleStimResponses,...
%             trainNoiseFlag, nAlternatives);
%     end
% 
%     % Train the classifier. This shows the usage to extact information
%     % from the container retrned as the first return value from the neural
%     % response engine - we index the responses by the string contained in
%     % the variable trainFlag (which was itself passed to the neural
%     % repsonse engine above.)
%     %
%     % Once extracted from the container, the responses are a 3 dimensional
%     % matrix, with the dimensions indexing [instancesNum x mNeuralDim x tTimeBins].
%     %   instancesNum   - number of response instances
%     %   mNeuralDim     - dimension of neural response at one timepoint
%     %   tTimeBins      - number of time points in stimulus sequence.
%     theClassifierEngine.compute('train', inSampleStimResponses(trainNoiseFlag),[]);
% 
%     % Save computed response instances
%     if (saveResponses)
%         responses.inSampleStimResponses = inSampleStimResponses;
%     end
% end
% 
% % Predict using trained classifier.
% %
% % Generate stimulus for prediction, NULL stimulus.  The variable testFlag
% % indicates what type of noise is used to generate the stimuli used for
% % prediction.  Typically 'random'.
% %
% % Note that for compatibility, put all the instances into a single
% % matrix in a single container at the end.
% 
% 
% useOldImplementation = false;
% if (useOldImplementation)
%     whichAlternatives = randi(nAlternatives,1,nTest);
%     outOfSampleStimResponsesCell = cell(1,nTest);
%     tic
%     for tt = 1:nTest
%         % Get responses for scene for this trial
%         [outOfSampleStimResponsesCell{tt}, ~] = theNeuralEngine.compute(...
%             theScenes{whichAlternatives(tt)}, ...
%             temporalSupport, ...
%             1, ...
%             'noiseFlags', {testNoiseFlag});
%     end
%     toc
% 
% else
%     % Alternative implementation here ...
%     assert(mod(nTest, nAlternatives) == 0, 'The number of test trials must be an integer multiple of the number of alternative choices');
%     nTestsPerAlternative = nTest/nAlternatives;
% 
%     whichAlternatives = [];
%     for theAlternative = 1:nAlternatives
%         tmp = repmat(theAlternative, [nTestsPerAlternative 1]);
%         whichAlternatives = cat(1, whichAlternatives, tmp);
%     end
% 
%     outOfSampleStimResponsesCell = cell(1, nAlternatives);
% 
%     eStart = tic;
%     for iAlternative = 1:nAlternatives
%         [outOfSampleStimResponsesCell{iAlternative}, ~] = theNeuralEngine.compute(...
%             theScenes{iAlternative}, ...
%             temporalSupport, ...
%             nTestsPerAlternative, ...
%             'noiseFlags', {testNoiseFlag});
%     end
%     e = toc(eStart);
%     fprintf('computePerformanceNWay_OneStimPerTrial: Took %0.1f secs to generate mean responses for all alternatives\n',e);
% 
% end
% 
% outOfSampleStimResponses = combineContainersMat(outOfSampleStimResponsesCell);
% 
% 
% % Do the prediction
% dataOut = theClassifierEngine.compute('predict', ...
%     outOfSampleStimResponses(testNoiseFlag),whichAlternatives);
% 
% % Save computed response instances
% if (saveResponses)
%     responses.outOfSampleStimResponses = outOfSampleStimResponse;
% end
% 
% % Set return variable.  For each trial 0 means wrong and 1 means right.
% % Taking mean(response) gives fraction correct.
% predictions = dataOut.trialPredictions;
% 
% end
% 
% function visualizeConeResps(theNeuralEngine, inSampleStimResponses,...
%     trainNoiseFlag, nAlternatives)
% if (isfield(theNeuralEngine.neuralPipeline, 'coneMosaic'))
%     diffResponse = inSampleStimResponses(trainNoiseFlag); % - inSampleNullStimResponses(trainNoiseFlag);
%     hFig = figure(998);
%     set(hFig, 'Position', [10 10 1600 400]);
% 
%     for aa = 1:nAlternatives
% 
%         % Contrast response
%         tmp = diffResponse{aa};
%         idx = theNeuralEngine.neuralPipeline.coneMosaic.lConeIndices;
%         meanLconeActivation = mean(tmp(idx(:)));
%         tmp(idx) = (tmp(idx)-meanLconeActivation)/meanLconeActivation;
%         idx = theNeuralEngine.neuralPipeline.coneMosaic.mConeIndices;
%         meanMconeActivation = mean(tmp(idx(:)));
%         tmp(idx) = (tmp(idx)-meanMconeActivation)/meanMconeActivation;
%         idx = theNeuralEngine.neuralPipeline.coneMosaic.sConeIndices;
%         if (~isempty(idx))
%             meanSconeActivation = mean(tmp(idx(:)));
%             tmp(idx) = (tmp(idx)-meanSconeActivation)/meanSconeActivation;
%         end
%         diffResponse{aa} = tmp;
% 
%         if (aa == 1)
%             ax = subplot(1, nAlternatives+1,1);
%             % Visualize the activation
%             theNeuralEngine.neuralPipeline.coneMosaic.visualize(...
%                 'figureHandle', hFig, ...
%                 'axesHandle', ax);
%         end
% 
%         ax = subplot(1, nAlternatives+1,aa+1);
%         % Visualize the activation
%         theNeuralEngine.neuralPipeline.coneMosaic.visualize(...
%             'figureHandle', hFig, ...
%             'axesHandle', ax, ...
%             'activation', squeeze(diffResponse{aa}), ...
%             'activationRange', prctile(squeeze(diffResponse{aa}),[1 99]), ...
%             'verticalActivationColorBarInside', true);
%     end
% 
%     % Also visualize the full absorptions density
%     figNo = 999;
%     theNeuralEngine.neuralPipeline.coneMosaic.visualizeFullAbsorptionsDensity(figNo);
% end
% 
% if (isfield(theNeuralEngine.neuralPipeline, 'mRGCmosaic'))
%     theNeuralEngine.neuralPipeline.mRGCmosaic.visualizeResponses(...
%         responseTemporalSupportSeconds, inSampleTestStimResponses(trainNoiseFlag), ...
%         'stimulusTemporalSupportSeconds', temporalSupport,...
%         'stimulusSceneSequence', testScene);
% end
% end