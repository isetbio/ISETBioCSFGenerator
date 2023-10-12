function [logThreshold, questObj, psychometricFunction, fittedPsychometricParams] = ...
    computeThresholdNWay_OneStimulusPerTrial(theSceneEngine, theNeuralEngine, classifierEngine, ...
    classifierPara, thresholdPara, questEnginePara, varargin)
% Compute contrast threshold for a given scene, neural response engine, and 
% classifier engine This function has been deprecated because there is a
% more general function computeThreshold.m that is applicable for both N-way
% forced choice and TAFC tasks.
%
% See descriptions: computeThreshold.m
%  
% History: 
%  12/07/21  dhb  Wrote NWay version from original TAFC version.
%  05/10/23  fh   Edited the code to direct people to the more general
%                   function computeThreshold.m and moved this to the 
%                   deprecated folder

warning(['This function has been deprecated. Consider using a more ',...
        'general function computeThreshold.m. ']);
varargin_appended = [varargin, 'TAFC', false];
[logThreshold, questObj, psychometricFunction, fittedPsychometricParams] = ...
    computeThreshold(theSceneEngine, theNeuralEngine, classifierEngine, ...
    classifierPara, thresholdPara, questEnginePara, varargin_appended{:});

% % Parse
% p = inputParser;
% p.addParameter('beVerbose',  true, @islogical);
% p.addParameter('extraVerbose',false, @islogical);
% p.addParameter('visualizeStimulus', false, @islogical);
% p.addParameter('visualizeAllComponents', false, @islogical);
% p.addParameter('datasavePara', [], @(x)(isempty(x)||(isstruct(x))));
% 
% parse(p, varargin{:});
% beVerbose = p.Results.beVerbose;
% visualizeStimulus = p.Results.visualizeStimulus;
% visualizeAllComponents = p.Results.visualizeAllComponents;
% datasavePara = p.Results.datasavePara;
% 
% % PF is based on correct/incorrect, so number of response
% % outcomes is 2 independent of number of alternatives.
% nOutcomes = 2;
% 
% % Construct a QUEST threshold estimator estimate threshold
% %
% % The questThreshold estimator is associated with a psychometric function,
% % by default qpPFWeibull.  It's important that the stimulus units be
% % compatable with those expected by the psychometric function.  qpPFWeibull
% % is coded in dB, which is confusing to many.  Probably better to work in
% % log10 units, which can be done by passing qpPFWeibullLog as the PF.
% estDomain  = -thresholdPara.logThreshLimitLow : thresholdPara.logThreshLimitDelta : -thresholdPara.logThreshLimitHigh;
% slopeRange = thresholdPara.slopeRangeLow: thresholdPara.slopeDelta : thresholdPara.slopeRangeHigh;
% 
% % Check whether the quest parameters include a PF, and set default if not.
% if (~isfield(questEnginePara,'qpPF'))
%     qpPF = @qpPFWeibull;
% else
%     qpPF = questEnginePara.qpPF;
% end
% 
% % Set defaults for guess/lapse rates if not passed.
% if (~isfield(thresholdPara,'guessRate'))
%     guessRate = 0.5;
% else
%     guessRate  = thresholdPara.guessRate ;
% end
% if (~isfield(thresholdPara,'lapseRate'))
%     lapseRate = 0;
% else
%     lapseRate = thresholdPara.lapseRate;
% end
% 
% % Handle quest method.
% if (isfield(questEnginePara, 'employMethodOfConstantStimuli'))&&(questEnginePara.employMethodOfConstantStimuli)
%     estimator = questThresholdEngine(...
%         'validation', true, 'nRepeat', questEnginePara.nTest, ...
%         'estDomain', estDomain, 'slopeRange', slopeRange, ...
%         'qpPF', qpPF, 'guessRate', guessRate, 'lapseRate', lapseRate);
% else
%     estimator = questThresholdEngine(...
%         'minTrial', questEnginePara.minTrial, 'maxTrial', questEnginePara.maxTrial, ...
%         'estDomain', estDomain, 'slopeRange', slopeRange, ...
%         'numEstimator', questEnginePara.numEstimator, ...
%         'stopCriterion', questEnginePara.stopCriterion, ...
%         'qpPF', qpPF, 'guessRate', guessRate, 'lapseRate', lapseRate);
% end
% 
% % Threshold estimation with QUEST+
% % Get the initial stimulus contrast from QUEST+
% [logContrast, nextFlag] = estimator.nextStimulus();
% 
% % Loop over trials.
% testedContrasts = [];
% 
% if ((isstruct(datasavePara)) && isfield(datasavePara, 'saveMRGCResponses') && (datasavePara.saveMRGCResponses) && ...
%         (isfield(datasavePara, 'destDir')) && (ischar(datasavePara.destDir)) && ...
%         (isfield(datasavePara, 'condExamined')))
%     neuralEngineSaved = false;
% else
%     datasavePara.saveMRGCResponses = false;
% end
% 
% % Dictionary to store the measured psychometric function which is returned to the user
% psychometricFunction = containers.Map();
% 
% while (nextFlag)
%     % Convert log contrast -> contrast
%     testContrast = 10 ^ logContrast;
% 
%     % Label for pCorrect dictionary
%     stimParamValueLabel = sprintf('paramNormalizedValue = %2.6f', testContrast);
%     if (beVerbose)
%         fprintf('Testing %s\n', stimParamValueLabel);
%     end
% 
%     % Have we already built the classifier for this contrast?
%     testedIndex = find(testContrast == testedContrasts);
%     if (isempty(testedIndex))
%         % No.  Save contrast in list
%         testedContrasts = [testedContrasts testContrast];
%         testedIndex = find(testContrast == testedContrasts);
% 
%         % Generate the scenes for each alternative, at the test contrast
%         for oo = 1:length(theSceneEngines)
%             [theTestSceneSequences{testedIndex}{oo}, theSceneTemporalSupportSeconds] = ...
%                 theSceneEngines{oo}.compute(testContrast);
%         end
% 
% 
%         % Train classifier for this TEST contrast and get predicted
%         % correct/incorrect predictions.  This function also computes the
%         % neural responses needed to train and predict.
%         [predictions, theTrainedClassifierEngines{testedIndex}, responses] = computePerformanceNWay_OneStimPerTrial(...
%             theTestSceneSequences{testedIndex}, ...
%             theSceneTemporalSupportSeconds, classifierPara.nTrain, classifierPara.nTest, ...
%             theNeuralEngine, classifierEngine, classifierPara.trainFlag, classifierPara.testFlag, ...
%             datasavePara.saveMRGCResponses, visualizeAllComponents);
% 
%         % Update the psychometric function with data point for this contrast level
%         psychometricFunction(stimParamValueLabel) = mean(predictions);
% 
%     else
%         % Classifier is already trained, just get predictions
%         [predictions, ~, ~] = computePerformanceNWay_OneStimPerTrial(...
%             theTestSceneSequences{testedIndex}, ...
%             theSceneTemporalSupportSeconds, classifierPara.nTrain, classifierPara.nTest, ...
%             theNeuralEngine, theTrainedClassifierEngines{testedIndex}, [], classifierPara.testFlag, ...
%             false);
% 
%         % Update the psychometric function with data point for this contrast level
%         previousData = psychometricFunction(stimParamValueLabel);
%         currentData = cat(2,previousData,mean(predictions));
%         psychometricFunction(stimParamValueLabel) = currentData;
%     end
% 
%     % Tell QUEST+ what we ran (how many trials at the given contrast) and
%     % get next stimulus contrast to run.
%     [logContrast, nextFlag] = ...
%         estimator.multiTrial(logContrast * ones(1, classifierPara.nTest), predictions);
% 
% end
% 
% % Return threshold value
% if (~isfield(thresholdPara,'thresholdCriterion'))
%     thresholdCriterion = 0.81606;
% else
%     thresholdCriterion = thresholdPara.thresholdCriterion;
% end
% [threshold, para] = estimator.thresholdMLE('showPlot', false, ...
%     'thresholdCriterion', thresholdCriterion);
% if (beVerbose)
%     fprintf('Maximum likelihood fit parameters: %0.2f, %0.2f, %0.2f, %0.2f\n', ...
%         para(1), para(2), para(3), para(4));
%     fprintf('Threshold (criterion proportion correct %0.4f: %0.2f (log10 units)\n', ...
%         thresholdCriterion,threshold);
% end
% 
% % Return the quest+ object wrapper for plotting and/or access to data
% questObj = estimator;
% 
% end