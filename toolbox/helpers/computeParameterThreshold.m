function [paramValueThreshold, questObj, psychometricFunction, fittedPsychometricParams]  = computeParameterThreshold(...
    theSceneEngines, theNeuralEngine, classifierEngine, classifierPara, ...
    thresholdPara, questEnginePara, varargin)
% Compute threshold for a parameter (not necessarily threshold) in an N-alternative forced choice paradigm
%
% Syntax:
%    [paramValueThreshold, questObj, psychometricFunction, fittedPsychometricParams] = computeParameterThreshold( ...
%        theSceneEngines, theNeuralEngine, classifierEngine, ...
%        classifierPara, thresholdPara, questEnginePara)  
%
% Description:
%     Uses Quest+ and the ISETBioCSFGenerator objects to obtain
%     computational observer threshold for the varied parameter in an N-Way AFC paradigm 
%
%
% Inputs:
%   theSceneEngines       - Cell array of scene engines, each responsible for a different of the N-alternative stimuli.
%   theNeuralEngine       - @neuralResponseEngine object to compute neural responses.
%   classifierEngine      - @responseClassifierEngine object that implements observer decision model.  
%   classifierPara        - Parameter struct associated with the classifier engine
%   thresholdPara         - Parameter struct associated with threshold estimation
%   questEnginePara       - Parameter struct for running the questThresholdEngine
%
% Outputs:
%   paramValueThreshold      - Estimated threshold value
%   questObj                 - questThresholdEngine object, which
%                              contains information about all the trials run.
%   psychometricFunction     - Dictionary (indexed by contrast level) with the 
%                              psychometric function.
%   fittedPsychometricParams - Parameters of psychometric function fit,
%                              matched to PF used in the questThresholdEngine object
%
% Optional key/value pairs:
%   'beVerbose'           - Logical. Provide some printout? Default false.
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
%   'roundEstDomainToIntegers' - Just keep the unique integers in the
%                           linear estimation domain.
%
% See also:
%    t_spatialCSF, t_thresholdEngine,
%    computePerformanceNWay_OneStimulusPerTrial.
%  

% History: 
%  03/01/22  NPC     Wrote it by adapting computeThresholdNWay_OneStimulusPerTrial

    % Parse
    p = inputParser;
    p.addParameter('beVerbose',  false, @islogical);
    p.addParameter('datasavePara', [], @(x)(isempty(x)||(isstruct(x))));
    p.addParameter('visualizeAllComponents', false, @islogical);
    p.addParameter('roundEstDomainToIntegers',false,@islogical);
    p.parse(varargin{:});
    datasavePara = p.Results.datasavePara;
    beVerbose = p.Results.beVerbose;
    visualizeAllComponents = p.Results.visualizeAllComponents;

    % Construct a QUEST threshold estimator
    variedParameterEstimationDomain = -thresholdPara.logThreshLimitLow : thresholdPara.logThreshLimitDelta : -thresholdPara.logThreshLimitHigh;
    if (p.Results.roundEstDomainToIntegers)
        variedParameterEstimationDomain = log10(unique(round(10.^variedParameterEstimationDomain)));
    end
    variedParameterSlopeRange = thresholdPara.slopeRangeLow: thresholdPara.slopeDelta : thresholdPara.slopeRangeHigh;
   
    % Check whether the quest parameters include a PF, and set default if not.
    if (~isfield(questEnginePara,'qpPF'))
        qpPF = @qpPFWeibull;
    else
        qpPF = questEnginePara.qpPF;
    end

    % Set defaults for guess/lapse rates if not passed.
    if (~isfield(thresholdPara,'guessRate'))
        fprintf(2, 'The guessRate threshold parameter was not set by the user. Assuming it is 0.5');
        guessRate = 0.5;
    else
        guessRate  = thresholdPara.guessRate ;
    end
    if (~isfield(thresholdPara,'lapseRate'))
        fprintf(2, 'The lapseRate threshold parameter was not set by the user. Assuming it is 0.0')
        lapseRate = 0;
    else
        lapseRate = thresholdPara.lapseRate;
    end

    % Handle quest method.
    if (isfield(questEnginePara, 'employMethodOfConstantStimuli'))&&(questEnginePara.employMethodOfConstantStimuli)
        estimator = questThresholdEngine(...
            'validation', true, 'blocked', true, 'nRepeat', questEnginePara.nTest, ...
            'estDomain',variedParameterEstimationDomain, 'slopeRange', variedParameterSlopeRange, ...
            'qpPF', qpPF, 'guessRate', guessRate, 'lapseRate', lapseRate);
    else
        if (classifierPara.nTest > 1)
            blockedVal = true;
        else
            blockedVal = false;
        end
        estimator = questThresholdEngine('blocked',blockedVal, ...
            'minTrial', questEnginePara.minTrial, 'maxTrial', questEnginePara.maxTrial, ...
            'estDomain', variedParameterEstimationDomain, 'slopeRange',  variedParameterSlopeRange, ...
            'numEstimator', questEnginePara.numEstimator, ...
            'stopCriterion', questEnginePara.stopCriterion, ...
            'qpPF', qpPF, 'guessRate', guessRate, 'lapseRate', lapseRate);
    end

    % Threshold estimation with QUEST+
    % Get the initial normalized parameter value from QUEST+
    [logNormalizedParamValue, nextFlag] = estimator.nextStimulus();

    % Loop over trials.
    testedNormalizedParamValues = [];

    if ((isstruct(datasavePara)) && isfield(datasavePara, 'saveMRGCResponses') && (datasavePara.saveMRGCResponses) && ...
            (isfield(datasavePara, 'destDir')) && (ischar(datasavePara.destDir)) && ...
            (isfield(datasavePara, 'condExamined')))
        neuralEngineSaved = false;
    else
        datasavePara.saveMRGCResponses = false;
    end

    % Dictionary to store the measured psychometric function which is returned to the user
    psychometricFunction = containers.Map();

    testCounter = 0;
    while (nextFlag)
        % Convert log normalized param value -> normalized param value
        normalizedParamValue = 10 ^ logNormalizedParamValue;
    
        % Convert normalized param value to actual param value
        paramValue = thresholdPara.maxParamValue*normalizedParamValue;

        % Label for pCorrect dictionary
        parameterLabel = sprintf('param value = %2.4f', paramValue);
        if (beVerbose)
            fprintf('computeParameterThreshold: logNormalizedParamValue = %0.2f; linear normalization factor = %0.2f; linear %s\n', ...
                logNormalizedParamValue, thresholdPara.maxParamValue, parameterLabel);
        end
    
        % Have we already built the classifier for this normalized param value?
        testedIndex = find(normalizedParamValue == testedNormalizedParamValues);
        if (isempty(testedIndex))
            % No.  Save this normalized param value in the list of examined
            % normalized param values
            if (isempty(testedIndex))
                % Only add to the list of tested values if we have not
                % tested this value before
                testedNormalizedParamValues(numel(testedNormalizedParamValues)+1) = normalizedParamValue;
            end
            testedIndex = find(normalizedParamValue == testedNormalizedParamValues);
            
            % Generate the scenes for each alternative, at the test param value
            for oo = 1:length(theSceneEngines)
                [theTestSceneSequences{testedIndex}{oo}, theSceneTemporalSupportSeconds] = ...
                    theSceneEngines{oo}.compute(paramValue);
            end
            
            % Train classifier for this parameter value and get predicted
            % correct/incorrect predictions.  This function also computes the
            % neural responses needed to train and predict.
            %
            % Note that because the classifer engine is a handle class, we
            % need to use a copy method to do the caching.  Otherwise the
            % cached pointer will simply continue to point to the same
            % object, and be updated by future training.
            eStart = tic;
            [predictions, tempClassifierEngine, responses] = computePerformanceNWay_OneStimPerTrial(...
                theTestSceneSequences{testedIndex}, ...
                theSceneTemporalSupportSeconds, classifierPara.nTrain, classifierPara.nTest, ...
                theNeuralEngine, classifierEngine, classifierPara.trainFlag, classifierPara.testFlag, ...
                datasavePara.saveMRGCResponses, visualizeAllComponents);
            theTrainedClassifierEngines{testedIndex} = tempClassifierEngine.copy;
            testCounter = testCounter + 1;
            e = toc(eStart);
            if (beVerbose)
                fprintf('computeParameterThreshold: Training and predicting test block %d took %0.1f secs\n',testCounter,e);
            end
            
            % Update the psychometric function with data point for this contrast level
            psychometricFunction(parameterLabel) = mean(predictions);
            if (beVerbose)
                fprintf('computeParamterThreshold: Length of psychometric function %d, test counter %d\n',length(psychometricFunction),testCounter);
            end
        else
            % Reality check
            if (estimator.validation & estimator.blocked)
                error('Should not be repeating any constrast for validation blocked method');
            end
            
            % Classifier is already trained, just get predictions
            eStart = tic;
            [predictions, ~, ~] = computePerformanceNWay_OneStimPerTrial(...
                theTestSceneSequences{testedIndex}, ...
                theSceneTemporalSupportSeconds, classifierPara.nTrain, classifierPara.nTest, ...
                theNeuralEngine, theTrainedClassifierEngines{testedIndex}, [], classifierPara.testFlag, ...
                false);
            testCounter = testCounter + 1;
            e = toc(eStart);
            if (beVerbose)
                fprintf('computeParameterThreshold: Predicting test block %d no training took %0.1f secs\n',testCounter,e);
            end

            % Update the psychometric function with data point for this parameter level
            previousData = psychometricFunction(parameterLabel);
            currentData = cat(2,previousData,mean(predictions));
            psychometricFunction(parameterLabel) = currentData;

            if (beVerbose)
                fprintf('computeParamterThreshold: Length of psychometric function %d, test counter %d\n',length(psychometricFunction),testCounter);
            end
        end
    
        % Tell QUEST+ what we ran (how many trials at the given normalized param value) and
        % get next normalized param value to run.
        eStart = tic;
        checkTwoWays = false;
        if (checkTwoWays & ~estimator.validation)
            questEstimatorSave = estimator.estimators{1};
            [logNormalizedParamValue1, nextFlag1] = ...
                estimator.multiTrial(logNormalizedParamValue * ones(1, classifierPara.nTest), predictions);
            estimator.estimators{1} = questEstimatorSave;
        end
        if (estimator.validation)
            % Method of constant stimuli
            [logNormalizedParamValue, nextFlag] = ...
                estimator.multiTrial(logNormalizedParamValue * ones(classifierPara.nTest,1), predictions);
        else
            % Quest
            [logNormalizedParamValue, nextFlag] = ...
                estimator.multiTrialQuestBlocked(logNormalizedParamValue * ones(classifierPara.nTest,1), predictions);
        end
        if (checkTwoWays & ~estimator.validation)
            if (logNormalizedParamValue ~= logNormalizedParamValue1)
                error('Two next values don''t match');
            end
            if (nextFlag ~= nextFlag1)
                error('Two next flags don''t match');
            end
        end
        e = toc(eStart);
        if (beVerbose)
            fprintf('computeParamterThreshold: Updating took %0.1f secs\n',e);
        end

    end  % (while nextFlag)

    if (beVerbose)
        fprintf('computeParameterThreshold: Ran %d test levels of %d trials per block of tests\n',testCounter,classifierPara.nTest);
        if (estimator.validation)
            fprintf('\tValidation mode, nRepeat set to %d\n',estimator.nRepeat);
        end
        fprintf('\tRecorded number single trials (%d) divided by number of blocks: %0.1f\n',testCounter,estimator.nTrial/testCounter);
    end

    % Return threshold criterion
    if (~isfield(thresholdPara,'thresholdCriterion'))
        thresholdCriterion = 0.81606;
    else
        thresholdCriterion = thresholdPara.thresholdCriterion;
    end

    % Param threshold (log normalized value)
    [logNormalizedParamValueThreshold, fittedPsychometricParams, thresholdDataOut] = estimator.thresholdMLE('showPlot', false, ...
        'thresholdCriterion', thresholdCriterion, 'returnData', true);

    % Convert log normalized param value -> normalized param value
    normalizedParamValueThreshold = 10 ^ logNormalizedParamValueThreshold;
    
    % Convert normalized param value to actual param value
    paramValueThreshold = thresholdPara.maxParamValue*normalizedParamValueThreshold;

    if (beVerbose)
        fprintf('Maximum likelihood fit parameters: %0.2f, %0.2f, %0.2f, %0.2f\n', ...
            fittedPsychometricParams(1), fittedPsychometricParams(2), fittedPsychometricParams(3), fittedPsychometricParams(4));
        fprintf('Threshold (criterion proportion correct %0.4f: %0.2f \n', ...
            thresholdCriterion, paramValueThreshold);
    end

    % Return the quest+ object wrapper for plotting and/or access to data
    questObj = estimator;
end

