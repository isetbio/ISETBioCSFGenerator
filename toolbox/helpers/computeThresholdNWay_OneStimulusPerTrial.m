function [threshold, questObj, psychometricFunction, para] = computeThresholdNWay_OneStimulusPerTrial(theSceneEngines, theNeuralEngine, classifierEngine, ...
    classifierPara, thresholdPara, questEnginePara, varargin)
% Compute contrast threshold for a given scene, neural response engine, and classifier engine
%
% Syntax:
%    [threshold, questObj, psychometricFunction, para] = ...
%        computeThresholdNWay_OneStimulusPerTrial(theSceneEngines, theNeuralEngine, ...
%        classifierEngine, classifierPara, thresholdPara, questEnginePara)  
%
% Description:
%    Uses Quest+ and the ISETBioCSFGenerator objects to obtain
%    computational observer contrast threshold for a given scene structure.
%
%    There is some art to using this function, in that you need to control
%    thrings such as how many trianing and test instances to use with the
%    classifier, how densely to tell Quest+ to sample the stimulus space
%    and over what range, etc.  The the three passed parameter structs
%    provide this control.  See t_thresholdEngine and t_spatialCSF for what
%    they control and some advice on how to set them. Indeed, understanding
%    those two tutorials should allow you to make effective use of this
%    function.
%
% Inputs:
%   theSceneEngines       - sceneEngine object for stimulus generation.
%                           This is a cell array with one scene per alternative.
%   theNeuralEngine       - neuralResponseEngine object
%   classifierEngine      - responseClassifierEngine
%   classifierPara        - Parameter struct associated with the classifier engine
%   thresholdPara         - Parameter struct associated with threshold estimation
%   questEnginePara       - Parameter struct for running the questThresholdEngine
%
% Outputs:
%   threshold             - Estimated threshold value
%   questObj              - questThresholdEngine object, which
%                           contains information about all the trials run.
%   psychometricFunction  - Dictionary (indexed by contrast level) with the 
%                           psychometric function.
%   para                  - Parameters of psychometric function fit,
%                           matched to PF used in the questThresholdEngine object
%
% Optional key/value pairs:
%   'beVerbose'           - Logical. Provide some printout? Default true.
%   'extraVerbose'        - Logical.  More detailed printout? Default false.
%   'visualizeStimulus'   - Logical. Provide stimulus visualization.
%                           Default false.
%   'visualizeAllComponents' - Logical. All component visualization.
%                           Default false. If set to true, it visualizes
%                           the mosaic responses to all the stimuli (multiple contrasts)
%                           which are computed by the neural engine.
%   'datasaveParameters'  - Parameters related to data saving. Default
%                           empty. When not empty, this has to be a struct
%                           with fields indicating which responses to save.
%                           Right now, the only accepted field is
%                           'saveMRGCResponses' which saved responses of
%                           the mRGC mosaic attached to an MRGC neural engine
%
% See also:
%    t_spatialCSF, t_thresholdEngine,
%    computePerformanceNWay_OneStimulusPerTrial.
%  

% History: 
%  12/07/21  dhb  Wrote NWay version from original TAFC version.

% Parse
p = inputParser;
p.addParameter('beVerbose',  true, @islogical);
p.addParameter('extraVerbose',false, @islogical);
p.addParameter('visualizeStimulus', false, @islogical);
p.addParameter('visualizeAllComponents', false, @islogical);
p.addParameter('datasavePara', [], @(x)(isempty(x)||(isstruct(x))));

parse(p, varargin{:});
beVerbose = p.Results.beVerbose;
visualizeStimulus = p.Results.visualizeStimulus;
visualizeAllComponents = p.Results.visualizeAllComponents;
datasavePara = p.Results.datasavePara;

% PF is based on correct/incorrect, so number of response
% outcomes is 2 independent of number of alternatives.
nOutcomes = 2;

% Construct a QUEST threshold estimator estimate threshold
%
% The questThreshold estimator is associated with a psychometric function,
% by default qpPFWeibull.  It's important that the stimulus units be
% compatable with those expected by the psychometric function.  qpPFWeibull
% is coded in dB, which is confusing to many.  Probably better to work in
% log10 units, which can be done by passing qpPFWeibullLog as the PF.
estDomain  = -thresholdPara.logThreshLimitLow : thresholdPara.logThreshLimitDelta : -thresholdPara.logThreshLimitHigh;
slopeRange = thresholdPara.slopeRangeLow: thresholdPara.slopeDelta : thresholdPara.slopeRangeHigh;

% Check whether the quest parameters include a PF, and set default if not.
if (~isfield(questEnginePara,'qpPF'))
    qpPF = @qpPFWeibull;
else
    qpPF = questEnginePara.qpPF;
end

% Set defaults for guess/lapse rates if not passed.
if (~isfield(thresholdPara,'guessRate'))
    guessRate = 0.5;
else
    guessRate  = thresholdPara.guessRate ;
end
if (~isfield(thresholdPara,'lapseRate'))
    lapseRate = 0;
else
    lapseRate = thresholdPara.lapseRate;
end

% Handle quest method.
if (isfield(questEnginePara, 'employMethodOfConstantStimuli'))&&(questEnginePara.employMethodOfConstantStimuli)
    estimator = questThresholdEngine('blocked',true,...
        'validation', true, 'nRepeat', questEnginePara.nTest, ...
        'estDomain', estDomain, 'slopeRange', slopeRange, ...
        'qpPF', qpPF, 'guessRate', guessRate, 'lapseRate', lapseRate);
else
    if (classifierPara.nTest > 1)
        blockedVal = true;
    else
        blockedVal = false;
    end
    estimator = questThresholdEngine('blocked',blockedVal,...
        'minTrial', questEnginePara.minTrial, 'maxTrial', questEnginePara.maxTrial, ...
        'estDomain', estDomain, 'slopeRange', slopeRange, ...
        'numEstimator', questEnginePara.numEstimator, ...
        'stopCriterion', questEnginePara.stopCriterion, ...
        'qpPF', qpPF, 'guessRate', guessRate, 'lapseRate', lapseRate);
end

% Threshold estimation with QUEST+
% Get the initial stimulus contrast from QUEST+
[logContrast, nextFlag] = estimator.nextStimulus();

% Loop over trials.
testedContrasts = [];

if ((isstruct(datasavePara)) && isfield(datasavePara, 'saveMRGCResponses') && ...
        (datasavePara.saveMRGCResponses) && (isfield(datasavePara, 'destDir')) &&...
        (ischar(datasavePara.destDir)) && (isfield(datasavePara, 'condExamined')))
    neuralEngineSaved = false;
else
    datasavePara.saveMRGCResponses = false;
end

% Dictionary to store the measured psychometric function which is returned to the user
psychometricFunction = containers.Map();

while (nextFlag)
    % Convert log contrast -> contrast
    testContrast = 10 ^ logContrast;
    
    % Label for pCorrect dictionary
    stimParamValueLabel = sprintf('paramNormalizedValue = %2.6f', testContrast);
    if (beVerbose)
        fprintf('Testing %s\n', stimParamValueLabel);
    end
    
    % Have we already built the classifier for this contrast?
    testedIndex = find(testContrast == testedContrasts);
    if (isempty(testedIndex))
        % No.  Save contrast in list
        testedContrasts = [testedContrasts testContrast];
        testedIndex = find(testContrast == testedContrasts);
        
        % Generate the scenes for each alternative, at the test contrast
        for oo = 1:length(theSceneEngines)
            [theTestSceneSequences{testedIndex}{oo}, theSceneTemporalSupportSeconds] = ...
                theSceneEngines{oo}.compute(testContrast);
        end
      
        
        % Train classifier for this TEST contrast and get predicted
        % correct/incorrect predictions.  This function also computes the
        % neural responses needed to train and predict.
        [predictions, theTrainedClassifierEngines{testedIndex}, responses] = ...
            computePerformanceNWay_OneStimulusPerTrial(...
            theTestSceneSequences{testedIndex}, ...
            theSceneTemporalSupportSeconds, classifierPara.nTrain, classifierPara.nTest, ...
            theNeuralEngine, classifierEngine, classifierPara.trainFlag, classifierPara.testFlag, ...
            datasavePara.saveMRGCResponses, visualizeAllComponents);
        
        % Update the psychometric function with data point for this contrast level
        psychometricFunction(stimParamValueLabel) = mean(predictions);
        
    else
        % Classifier is already trained, just get predictions
        [predictions, ~, ~] = computePerformanceNWay_OneStimulusPerTrial(...
            theTestSceneSequences{testedIndex}, ...
            theSceneTemporalSupportSeconds, classifierPara.nTrain, classifierPara.nTest, ...
            theNeuralEngine, theTrainedClassifierEngines{testedIndex}, [], classifierPara.testFlag, ...
            false);
        
        % Update the psychometric function with data point for this contrast level
        previousData = psychometricFunction(stimParamValueLabel);
        currentData = cat(2,previousData,mean(predictions));
        psychometricFunction(stimParamValueLabel) = currentData;
    end
    
    % Tell QUEST+ what we ran (how many trials at the given contrast) and
    % get next stimulus contrast to run.
    if (estimator.validation)
        % Method of constant stimuli
        [logContrast, nextFlag] = ...
            estimator.multiTrial(logContrast * ones(classifierPara.nTest,1), predictions');
    else
        % Quest
        [logContrast, nextFlag] = ...
          estimator.multiTrialQuestBlocked(logContrast * ones(classifierPara.nTest,1), predictions');
    end

end

% Return threshold value
if (~isfield(thresholdPara,'thresholdCriterion'))
    thresholdCriterion = 0.81606;
else
    thresholdCriterion = thresholdPara.thresholdCriterion;
end
[threshold, para] = estimator.thresholdMLE('showPlot', false, ...
    'thresholdCriterion', thresholdCriterion);
if (beVerbose)
    fprintf('Maximum likelihood fit parameters: %0.2f, %0.2f, %0.2f, %0.2f\n', ...
        para(1), para(2), para(3), para(4));
    fprintf('Threshold (criterion proportion correct %0.4f: %0.2f (log10 units)\n', ...
        thresholdCriterion,threshold);
end

% Return the quest+ object wrapper for plotting and/or access to data
questObj = estimator;

end