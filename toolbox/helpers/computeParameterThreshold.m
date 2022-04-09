function [paramValueThreshold, questObj, psychometricFunction, fittedPsychometricParams]  = computeParameterThreshold(...
    theSceneEngines, theNeuralEngine, classifierEngine, classifierPara, ...
    thresholdPara, questEnginePara, varargin)
    
    % Parse
    p = inputParser;
    p.addParameter('beVerbose',  true, @islogical);
    p.addParameter('datasavePara', [], @(x)(isempty(x)||(isstruct(x))));
    p.addParameter('visualizeAllComponents', false, @islogical);
    p.parse(varargin{:});
    datasavePara = p.Results.datasavePara;
    beVerbose = p.Results.beVerbose;
    visualizeAllComponents = p.Results.visualizeAllComponents;


    % Construct a QUEST threshold estimator
    variedParameterEstimationDomain = -thresholdPara.logThreshLimitLow : thresholdPara.logThreshLimitDelta : -thresholdPara.logThreshLimitHigh;
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
            'validation', true, 'nRepeat', questEnginePara.nTest, ...
            'estDomain',variedParameterEstimationDomain, 'slopeRange', variedParameterSlopeRange, ...
            'qpPF', qpPF, 'guessRate', guessRate, 'lapseRate', lapseRate);
    else
        estimator = questThresholdEngine(...
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

    while (nextFlag)
        % Convert log normalized param value -> normalized param value
        normalizedParamValue = 10 ^ logNormalizedParamValue;
    
        % Convert normalized param value to actual param value
        paramValue = thresholdPara.maxParamValue*normalizedParamValue;

        % Label for pCorrect dictionary
        parameterLabel = sprintf('param value= %2.4f', paramValue);
        if (beVerbose)
            fprintf('Testing %s\n', parameterLabel);
        end
    
        % Have we already built the classifier for this normalized param value?
        testedIndex = find(normalizedParamValue == testedNormalizedParamValues);
        if (isempty(testedIndex))
            % No.  Save this normalized param value in the list of examined
            % normalized param values
            testedNormalizedParamValues(numel(testedNormalizedParamValues)+1) = normalizedParamValue;
            testedIndex = find(normalizedParamValue == testedNormalizedParamValues);
            
            % Generate the scenes for each alternative, at the test param value
            for oo = 1:length(theSceneEngines)
                [theTestSceneSequences{testedIndex}{oo}, theSceneTemporalSupportSeconds] = ...
                    theSceneEngines{oo}.compute(paramValue);
            end
          
            
            % Train classifier for this parameter value and get predicted
            % correct/incorrect predictions.  This function also computes the
            % neural responses needed to train and predict.
            [predictions, theTrainedClassifierEngines{testedIndex}, responses] = computePerformanceNWay_OneStimPerTrial(...
                theTestSceneSequences{testedIndex}, ...
                theSceneTemporalSupportSeconds, classifierPara.nTrain, classifierPara.nTest, ...
                theNeuralEngine, classifierEngine, classifierPara.trainFlag, classifierPara.testFlag, ...
                datasavePara.saveMRGCResponses, visualizeAllComponents);
            
            % Update the psychometric function with data point for this contrast level
            psychometricFunction(parameterLabel) = mean(predictions);
            
    
        else
            % Classifier is already trained, just get predictions
            [predictions, ~, ~] = computePerformanceNWay_OneStimPerTrial(...
                theTestSceneSequences{testedIndex}, ...
                theSceneTemporalSupportSeconds, classifierPara.nTrain, classifierPara.nTest, ...
                theNeuralEngine, theTrainedClassifierEngines{testedIndex}, [], classifierPara.testFlag, ...
                false);
            
            % Update the psychometric function with data point for this parameter level
            previousData = psychometricFunction(parameterLabel);
            currentData = cat(2,previousData,mean(predictions));
            psychometricFunction(parameterLabel) = currentData;
        end
    
        % Tell QUEST+ what we ran (how many trials at the given normalized param value) and
        % get next normalized param value to run.
        [logNormalizedParamValue, nextFlag] = ...
            estimator.multiTrial(logNormalizedParamValue * ones(1, classifierPara.nTest), predictions);

    end  % (while nextFlag)

    % Return threshold criterion
    if (~isfield(thresholdPara,'thresholdCriterion'))
        thresholdCriterion = 0.81606;
    else
        thresholdCriterion = thresholdPara.thresholdCriterion;
    end

    % Param threshold (log normalized value)
    [logNormalizedParamValueThreshold, fittedPsychometricParams] = estimator.thresholdMLE('showPlot', false, ...
        'thresholdCriterion', thresholdCriterion);

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

