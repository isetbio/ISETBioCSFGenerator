function [logThreshold, questObj, psychometricFunction, fittedPsychometricParams] = computeThreshold(task, theSceneEngine, theNeuralEngine, classifierEngine, ...
    classifierPara, thresholdPara, questEnginePara, varargin)
% Compute contrast threshold for a given scene, neural response engine, and classifier engine
%
% Syntax:
%    [logThreshold, questObj, psychometricFunction, fittedPsychometricParams] = computeThresholdTAFC( ...
%        theSceneEngine, theNeuralEngine, classifierEngine, ...
%        classifierPara, thresholdPara, questEnginePara)  
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
%   theSceneEngine        - sceneEngine object for stimulus generation
%   theNeuralEngine       - neuralResponseEngine object
%   classifierEngine      - responseClassifierEngine
%   classifierPara        - Parameter struct associated with the classifier engine
%   thresholdPara         - Parameter struct associated with threshold estimation
%   questEnginePara       - Parameter struct for running the questThresholdEngine
%
% Outputs:
%   logThreshold              - Estimated threshold value (log)
%   questObj                  - questThresholdEngine object, which
%                               contains information about all the trials run.
%   psychometricFunction      - Dictionary (indexed by contrast level) with the 
%                               psychometric function.
%   fittedPsychometricParams  - Parameters of psychometric function fit,
%                               matched to PF used in the questThresholdEngine object
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
%    t_spatialCSF, t_thresholdEngine, computePerformanceTAFC
%  

% History: 
%  10/23/20  dhb  Added commments.
%  12/04/20  npc  Added option to run method of constant stimuli. 
%                 Also added psychometricFunction return argument.
%  05/01/23  npc  Modifications to use the multiTrialQuestBlocked method
%  06/08/23  npc  Modifications to run under the BetterCaching branch
%                 (copyable responseClassifierEngine)

%check if the input task name is valid
if ~strcmp(task, 'TAFC') && ~strcmp(task,'NWay_OneStimulusPerTrial')
    error('Invalid task name! You can either input ''TAFC'' or ''NWay_OneStimulusPerTrial''');
end

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
    estimator = questThresholdEngine(...
        'validation', true, 'blocked', true, 'nRepeat', questEnginePara.nTest, ...
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

% Generate the NULL stimulus (zero contrast)
if strcmp(task,'TAFC')
    nullContrast = 0.0;
    [theNullSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(nullContrast);
end

% Some diagnosis
if (p.Results.extraVerbose)
    theWl = 400;
    theFrame = 1;
    index = find(theNullSceneSequence{theFrame}.spectrum.wave == theWl);
    temp = theNullSceneSequence{theFrame}.data.photons(:,:,index);
    fprintf('At %d nm, frame %d, null scene mean, min, max: %g, %g, %g\n',...
        theWl,theFrame,mean(temp(:)),min(temp(:)),max(temp(:)));
    theWl = 550;
    index = find(theNullSceneSequence{theFrame}.spectrum.wave == theWl);
    temp = theNullSceneSequence{theFrame}.data.photons(:,:,index);
    fprintf('At %d nm, frame %d, null scene mean, min, max: %g, %g, %g\n',...
        theWl,theFrame,mean(temp(:)),min(temp(:)),max(temp(:)));
end


% Threshold estimation with QUEST+
% Get the initial stimulus contrast from QUEST+
[logContrast, nextFlag] = estimator.nextStimulus();

% Loop over trials.
testedContrasts = [];

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
    % Convert log contrast -> contrast
    testContrast = 10 ^ logContrast;
    
    % Label for pCorrect dictionary
    contrastLabel = sprintf('C = %2.4f%%', testContrast*100);
    if (beVerbose)
        fprintf('Testing %s\n', contrastLabel);
    end
    
    % Have we already built the classifier for this contrast?
    testedIndex = find(testContrast == testedContrasts);
    if (isempty(testedIndex))
        % No.  Save contrast in list
        testedContrasts(numel(testedContrasts)+1) = testContrast;
        testedIndex = find(testContrast == testedContrasts);

        % Generate the scenes for each alternative, at the test contrast
        if length(theSceneEngine) >= 2
            for oo = 1:length(theSceneEngine)
                [theSceneSequences{testedIndex}{oo}, theSceneTemporalSupportSeconds] = ...
                    theSceneEngine{oo}.compute(testContrast);
            end
        else
            [theSceneSequence{testedIndex}, theSceneTemporalSupportSeconds] = ...
                    theSceneEngine.compute(testContrast);
            %concatenate theTestSceneSequences and theNullSceneSequences
            %together
            theSceneSequences{testedIndex} = {theSceneSequence{testedIndex},...
                theNullSceneSequence};
        end
        
        % Some diagnosis
        if (p.Results.extraVerbose)
            theWl = 400;
            theFrame = 1;
            index = find(theTestSceneSequences{testedIndex}{theFrame}.spectrum.wave == theWl);
            temp = theTestSceneSequences{testedIndex}{theFrame}.data.photons(:,:,index);
            fprintf('At %d nm, frame %d, test scene %d mean, min, max: %g, %g, %g\n',...
                theWl,theFrame,testedIndex,mean(temp(:)),min(temp(:)),max(temp(:)));
            theWl = 550;
            index = find(theTestSceneSequences{testedIndex}{theFrame}.spectrum.wave == theWl);
            temp = theTestSceneSequences{testedIndex}{theFrame}.data.photons(:,:,index);
            fprintf('At %d nm, frame %d, test scene %d mean, min, max: %g, %g, %g\n',...
                theWl,theFrame,testedIndex,mean(temp(:)),min(temp(:)),max(temp(:)));
        end
        
        % Visualize the drifting sequence
        if (visualizeStimulus)
            theSceneEngine.visualizeSceneSequence(...
                theSceneSequences{testedIndex}{1}, theSceneTemporalSupportSeconds);
        end
        
        
        % Update the classifier engine pooling params for this particular test contrast
        if (isfield(classifierEngine.classifierParams, 'pooling')) && ...
           (~strcmp(classifierEngine.classifierParams.pooling, 'none'))
            
            fprintf('Computing pooling kernels for contrast %f\n', testContrast*100);
       
            % Compute pooling weights
            switch (classifierEngine.classifierParams.pooling.type)
                case 'linear'
                   [poolingWeights.direct, ~, ...
                       noiseFreeNullResponse, noiseFreeTestResponse] = ...
                       spatioTemporalPoolingWeights(theNeuralEngine,...
                       theTestSceneSequences{testedIndex},...
                       theNullSceneSequence, theSceneTemporalSupportSeconds);
                   
                case 'quadratureEnergy'
                   [poolingWeights.direct, poolingWeights.quadrature, ...
                       noiseFreeNullResponse, noiseFreeTestResponse] = ...
                       spatioTemporalPoolingWeights(theNeuralEngine,...
                       theTestSceneSequences{testedIndex},...
                       theNullSceneSequence, theSceneTemporalSupportSeconds);
                   
                otherwise
                    error('Unknown classifier engine pooling type: ''%s''.',...
                        classifierEngine.classifierParams.pooling.type)
            end
            
            % Update the classifier engine's pooling weights
            classifierEngine.updateSpatioTemporalPoolingWeightsAndNoiseFreeResponses(...
                poolingWeights,noiseFreeNullResponse, noiseFreeTestResponse);
        end
        
        % Train classifier for this TEST contrast and get predicted
        % correct/incorrect predictions.  This function also computes the
        % neural responses needed to train and predict.

        %
        % Note that because the classifer engine is a handle class, we
        % need to use a copy method to do the caching.  Otherwise the
        % cached pointer will simply continue to point to the same
        % object, and be updated by future training.
        eStart = tic;
        % [predictions, tempClassifierEngine, responses] = computePerformanceTAFC(...
        %     theSceneSequences{testedIndex}, theSceneTemporalSupportSeconds,...
        %     classifierPara.nTrain, classifierPara.nTest, theNeuralEngine,...
        %     classifierEngine, classifierPara.trainFlag, classifierPara.testFlag, ...
        %     datasavePara.saveMRGCResponses, visualizeAllComponents);
        [predictions, tempClassifierEngine, responses] = computePerformance(...
            task, theSceneSequences{testedIndex}, theSceneTemporalSupportSeconds,...
            classifierPara.nTrain, classifierPara.nTest, theNeuralEngine,...
            classifierEngine, classifierPara.trainFlag, classifierPara.testFlag, ...
            datasavePara.saveMRGCResponses, visualizeAllComponents);
        
        % Copy the trained classifier
        theTrainedClassifierEngines{testedIndex} = tempClassifierEngine.copy;

        testCounter = testCounter + 1;
        e = toc(eStart);
        if (beVerbose)
                fprintf('computeThresholdTAFC: Training and predicting test block %d took %0.1f secs\n',testCounter,e);
        end

        % Update the psychometric function with data point for this contrast level
        psychometricFunction(contrastLabel) = mean(predictions);
        
        if (beVerbose)
            fprintf('computeThresholdTAFC: Length of psychometric function %d, test counter %d\n',...
                length(psychometricFunction),testCounter);
        end

        % Save computed responses only the first time we test this contrast
        if (datasavePara.saveMRGCResponses)
            theMRGCmosaic = theNeuralEngine.neuralPipeline.mRGCmosaic;
            
            % Save neural engine
            if (neuralEngineSaved == false)
                if (~exist(datasavePara.destDir, 'dir'))
                    fprintf('Creating destination directory: ''%s''\n.',datasavePara.destDir);
                    mkdir(datasavePara.destDir);
                end
                mosaicFileName = fullfile(datasavePara.destDir,'mRGCmosaic.mat');
                fprintf('Saving mRGC mosaic to %s.\n', mosaicFileName);
                save(mosaicFileName, 'theMRGCmosaic', '-v7.3');
                neuralEngineSaved = true;
            end
            
            % Save responses
            responseFileName = fullfile(datasavePara.destDir, ...
                sprintf('responses_%s_ContrastLevel_%2.2f.mat', datasavePara.condExamined, testContrast*100));
            fprintf('Saving computed responses to %s.\n', responseFileName);
            save(responseFileName, 'responses',  '-v7.3');
        end
    else

        % Reality check
        if (estimator.validation & estimator.blocked)
             error('Should not be repeating any constrast for validation blocked method');
        end

        % Classifier is already trained, just get predictions
        eStart = tic;
        [predictions, ~, ~] = computePerformance(task, theSceneSequences{testedIndex}, ...
            theSceneTemporalSupportSeconds, classifierPara.nTrain, classifierPara.nTest, ...
            theNeuralEngine, theTrainedClassifierEngines{testedIndex}, [], classifierPara.testFlag, ...
            false, false);
        
        testCounter = testCounter + 1;
        e = toc(eStart);
        if (beVerbose)
            fprintf('computeThresholdTAFC: Predicting test block %d no training took %0.1f secs\n',testCounter,e);
        end

        % Update the psychometric function with data point for this contrast level
        previousData = psychometricFunction(contrastLabel);
        currentData = cat(2,previousData,mean(predictions));
        psychometricFunction(contrastLabel) = currentData;

        if (beVerbose)
           fprintf('computeThresholdTAFC: Length of psychometric function %d, test counter %d\n',...
               length(psychometricFunction),testCounter);
        end
    end
    

    % Tell QUEST+ what we ran (how many trials at the given contrast) and
    % get next stimulus contrast to run.
    eStart = tic;
    if (estimator.validation)
        % Method of constant stimuli
        [logContrast, nextFlag] = ...
            estimator.multiTrial(logContrast * ones( classifierPara.nTest,1), predictions);
    else
        % Quest
        [logContrast, nextFlag] = ...
          estimator.multiTrialQuestBlocked(logContrast * ones(classifierPara.nTest,1), predictions);
    end

    e = toc(eStart);
    if (beVerbose)
            fprintf('computeThresholdTAFC: Updating took %0.1f secs\n',e);
    end

end  % while (nextFlag)


if (beVerbose)
   fprintf('computeThresholdTAFC: Ran %d test levels of %d trials per block of tests\n',...
       testCounter,classifierPara.nTest);
   if (estimator.validation)
            fprintf('\tValidation mode, nRepeat set to %d\n',estimator.nRepeat);
   end
   fprintf('\tRecorded number single trials (%d) divided by number of blocks: %0.1f\n',...
       testCounter,estimator.nTrial/testCounter);
end

% Return threshold value. For the mQUESTPlus Weibull PFs, the first
% parameter of the PF fit is the 0.81606 proportion correct threshold,
% when lapse rate is 0 and guess rate is 0.5.  Better to make this an
% explicit parameter, however.  We use default of 0.81606 if not passed for
% backward compatibility.
if (~isfield(thresholdPara,'thresholdCriterion'))
    thresholdCriterion = 0.81606;
else
    thresholdCriterion = thresholdPara.thresholdCriterion;
end

% Param threshold (log value)
[logThreshold, fittedPsychometricParams, thresholdDataOut] = ...
    estimator.thresholdMLE('showPlot', false, ...
    'thresholdCriterion', thresholdCriterion, 'returnData', true);


if (beVerbose)
    fprintf('Maximum likelihood fit parameters: %0.2f, %0.2f, %0.2f, %0.2f\n', ...
            fittedPsychometricParams(1), fittedPsychometricParams(2),...
            fittedPsychometricParams(3), fittedPsychometricParams(4));
    fprintf('Threshold (criterion proportion correct %0.4f: %0.2f (log10 units)\n', ...
        thresholdCriterion,logThreshold);
end

% Return the quest+ object wrapper for plotting and/or access to data
questObj = estimator;

end