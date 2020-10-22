function [threshold, questObj] = computeThreshold(theSceneEngine, theNeuralEngine, classifierEngine, classifierPara, thresholdPara, questEnginePara)
%computeThreshold  Compute contrast threshold for a given scene (i.e.,
%stimulus), neural response engine, and classifier engine
%
% Usage:
%  computeThreshold(theSceneEngine, theNeuralEngine, classifierEngine, classifierPara, thresholdPara, questEnginePara);
%   See t_spatialCSF.m for example usage
%   
%
% Inputs:
%   theSceneEngine        - sceneEngine object for stimulus generation
%   theNeuralEngine       - neuralResponseEngine object
%   classifierEngine         - responseClassifierEngine
%   classifierPara             - Parameter associated with the classifier engine
%   thresholdPara           - Parameter associated with threshold estimation
%   questEnginePara      - Parameter for running the questThresholdEngine
%
%   Also see t_thresholdEngine.m and t_spatialCSF.m
%
% Outputs:
%   threshold                    - Estimated threshold value
%   questObj                     - questThresholdEngine object, which
%                                           contains reference to all the simulated 
%                                           stimulus - response data

% Construct a QUEST threshold estimator estimate threshold
estDomain  = -thresholdPara.logThreshLimitLow : thresholdPara.logThreshLimitDelta : -thresholdPara.logThreshLimitHigh;
slopeRange = thresholdPara.slopeRangeLow: thresholdPara.slopeDelta : thresholdPara.slopeRangeHigh;

estimator = questThresholdEngine('minTrial', questEnginePara.minTrial, 'maxTrial', questEnginePara.maxTrial, ...
    'estDomain', estDomain, 'slopeRange', slopeRange, 'numEstimator', questEnginePara.numEstimator);

% Generate the NULL stimulus (zero contrast)
nullContrast = 0.0;
[theNullSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(nullContrast);

% Threshold estimation with QUEST+
% Get the initial stimulus contrast from QUEST+
[logContrast, nextFlag] = estimator.nextStimulus();

% Loop over trials.
testedContrasts = [];
while (nextFlag)
    
    % Convert log contrast -> contrast
    testContrast = 10 ^ logContrast;
    
    % Have we already built the classifier for this contrast?
    testedIndex = find(testContrast == testedContrasts);
    if (isempty(testedIndex))
        % No.  Save contrast in list
        testedContrasts = [testedContrasts testContrast];
        testedIndex = find(testContrast == testedContrasts);
        
        % Generate the TEST scene sequence for the given contrast
        [theTestSceneSequences{testedIndex}, ~] = theSceneEngine.compute(testContrast);
        
        % Train classifier for this TEST contrast and get predicted
        % correct/incorrect predictions.  This function also computes the
        % neural responses needed to train and predict.
        [predictions, theTrainedClassifierEngines{testedIndex}] = computePerformance(...
            theNullSceneSequence, theTestSceneSequences{testedIndex}, ...
            theSceneTemporalSupportSeconds, classifierPara.nTrain, classifierPara.nTest, ...
            theNeuralEngine, classifierEngine, classifierPara.trainFlag, classifierPara.testFlag);
        
    else
        % Classifier is already trained, just get predictions
        predictions = computePerformance(...
            theNullSceneSequence, theTestSceneSequences{testedIndex}, ...
            theSceneTemporalSupportSeconds, classifierPara.nTrain, classifierPara.nTest, ...
            theNeuralEngine, theTrainedClassifierEngines{testedIndex}, [], classifierPara.testFlag);
    end
    
    % Tell QUEST+ what we ran (how many trials at the given contrast) and
    % get next stimulus contrast to run.
    [logContrast, nextFlag] = ...
        estimator.multiTrial(logContrast * ones(1, classifierPara.nTest), predictions);

end

% Return threshold value
[threshold, para] = estimator.thresholdMLE('showPlot', false);
fprintf('Maximum likelihood fit parameters: %0.2f, %0.2f, %0.2f, %0.2f\n', ...
    para(1), para(2), para(3), para(4));

% Return the quest+ object wrapper
% for plotting and/or access to data
questObj = estimator;

end