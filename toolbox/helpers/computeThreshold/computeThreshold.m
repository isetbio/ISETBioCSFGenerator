function [logThreshold, questObj, psychometricFunction, estimatorFittedPsychometricParams, ...
    trialByTrialStimulusAlternatives,trialByTrialPerformance,triailByTrialWhichResponses] = ...
    computeThreshold(theSceneEngine, theNeuralEngine, classifierEngine, ...
    classifierPara, thresholdPara, questEnginePara, varargin)
% Compute contrast threshold for different scenes, neural response engine, 
% and classifier engine. This function is suitable for both TAFC and 
% N-alternative forced-choice tasks. 
%
% Syntax:
%    [logThreshold, questObj, psychometricFunction, fittedPsychometricParams, ...
%        trialByTrialStimulusAlternatives,trialByTrialPerformance,triailByTrialWhichResponses] = computeThreshold( ...
%        theSceneEngine, theNeuralEngine, classifierEngine, ...
%        classifierPara, thresholdPara, questEnginePara)  
%
% Description:
%    Uses Quest+ and the ISETBioCSFGenerator objects to obtain
%    computational observer contrast threshold for a given scene structure.
%
%    There is some art to using this function, in that you need to control
%    things such as how many trianing and test instances to use with the
%    classifier, how densely to tell Quest+ to sample the stimulus space
%    and over what range, etc.  The the three passed parameter structs
%    provide this control.  See t_spatialCSF and t_spatialCSF for what
%    they control and some advice on how to set them. Indeed, understanding
%    those two tutorials should allow you to make effective use of this
%    function.
%
%    The thresholdPara structure may contain a field maxParamValue. This is
%    used to scale the linear parameter values returned by the
%    questThresholdEngine, outside of that engine.  Thus the values in
%    questObj are not scaled by this parameter, but those simulated by
%    theSceneEngine are so scaled. See outputs below for a detailed
%    descripiont of where the scalar has been applied and where not. This
%    is probably not the best design, and in general we rather strongly
%    recommend that you not use the maxParamValue field, so that all the
%    stimulus values are consistent across the code. We provide this field
%    for backwards compatibilty with some older code.
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
%   logThreshold              - Estimated threshold value (log). The linear
%                               value that this is the log10 of has been
%                               scaled by thresholdParam.maxParamValue.
%   questObj                  - questThresholdEngine object, which
%                               contains information about all the trials
%                               run. The values in here are not scaled by
%                               thresholdParam.maxParamValue.
%   psychometricFunction      - Dictionary (aka Container, indexed by contrast level) with the 
%                               psychometric function. The dictionary
%                               entries are indexed by a string containing
%                               linear parameter values multiplied by 100
%                               (to make it percent if it is linear
%                               contrast).  This linear value has been
%                               multiplied by thresholdParam.maxParamValue.
%   estimatorFittedPsychometricParams  - Parameters of psychometric function fit,
%                               matched to PF used in the
%                               questThresholdEngine object. These have not
%                               been scaled by thresholdParm.maxParamValue.
%   trialByTrialStimulusAlternatives - Dictionary with same keys as
%                               psychometric fuction that has a column
%                               vector for each contrast level that gives
%                               the trial-by-trial stimulus alternative
%                               (integer 1-N) for each trial at that
%                               contrast level.
%   trialByTrialPerformance   - One more dictionary matched to
%                               trialByTrialStimulusAlternatives that gives correct (1) or incorrect
%                               (0) for each trial.
%   triailByTrialWhichResponses     - One more dictionary matched to
%                               trialByTrialStimulusAlternatives that gives response (same codes as whichAlternatives) for each trial.
%
% Optional key/value pairs:
%   'verbose'           - Logical. Provide some printout? Default true.
%   'extraVerbose'        - Logical.  More detailed printout? Default false.
%   'visualizeStimulus'   - Logical. Provide stimulus visualization.
%                           Default false.https://canvas.upenn.edu/courses/1841993v
%   'visualizeAllComponents' - Logical. All component visualization.
%                           Default false. If set to true, it visualizes
%                           the mosaic responses to all the stimuli (multiple contrasts)
%                           which are computed by the neural engine.
%   'datasavePara'        - Parameters related to data saving. Default
%                           empty. When not empty, this has to be a struct
%                           with fields indicating which responses to save.
%                           Right now, the only accepted field is
%                           'saveMRGCResponses' which saved responses of
%                           the mRGC mosaic attached to an MRGC neural engine
%   'TAFC'                - logical. Whether this is a two-interval forced-choice task
%                           or N-way one-stimulus-per-trial task. Default
%                           false.  When it is true, this routine expects
%                           to receive a single scene engine and figures
%                           out "contrast" threshold for discriinating
%                           that against a contrast of 0.  There are other
%                           sorts of TAFC tasks one could consider, and
%                           this routine would need to be elaborated to
%                           handle them, or else the caller would need to
%                           explicitly formulate the TAFC version as an
%                           2-alternative with one stimulus per trial problem.
%     useMetaContrast     - Passed arguments are for meta contrast setup
%     trainFixationalEM   - EM object with EM paths for training.  This can
%                           be empty, which means no EM on training.ƒthr
%     testFixaitonalEM    - EM object with EM paths for testing. This can
%                           be the same as trainFixationalEM.  Default
%                           empty, which means no EM on testing.
%     maxVisualizedNoisyResponseInstanceStimuli - Maximum number of stimuli
%                            for which the visualizeNoisyResponseInstances
%                            flag stays true.

%
% See also:
%    t_spatialCSF, computePerformance  

% History: 
%  10/23/20  dhb  Added commments.
%  12/04/20  npc  Added option to run method of constant stimuli. 
%                 Also added psychometricFunction return argument.
%  05/01/23  npc  Modifications to use the multiTrialQuestBlocked method
%  06/08/23  npc  Modifications to run under the BetterCaching branch
%                 (copyable responseClassifierEngine)
%  04/10/24  fh   Merged computeThresholdTAFC.m and
%                   computeThresholdNWay_OneStimulusPerTrial.m by adding 
%                   a key/pair pair specifying whether the task is TAFC.
%  12/23/24  dhb  Simplified as part of new architecture push.

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('verbose',  true, @islogical);
p.addParameter('extraVerbose', false, @islogical);
p.addParameter('visualizeStimulus', false, @islogical);
p.addParameter('visualizeAllComponents', false, @islogical);
p.addParameter('datasavePara', [], @(x)(isempty(x) || (isstruct(x))));
p.addParameter('TAFC', false, @islogical);
p.addParameter('useMetaContrast', false, @islogical);
p.addParameter('trainFixationalEM', [], @(x)(isempty(x) || (isa(x,'fixationalEM'))));
p.addParameter('testFixationalEM', [], @(x)(isempty(x) || (isa(x,'fixationalEM'))));
p.addParameter('maxVisualizedNoisyResponseInstanceStimuli',1,@isnumeric);

parse(p, varargin{:});
verbose = p.Results.verbose;
visualizeStimulus = p.Results.visualizeStimulus;
visualizeAllComponents = p.Results.visualizeAllComponents;
datasavePara = p.Results.datasavePara;
isTAFC = p.Results.TAFC;
trainFixationalEMObj = p.Results.trainFixationalEM;
testFixationalEMObj = p.Results.testFixationalEM;
maxVisualizedNoisyResponseInstanceStimuli = p.Results.maxVisualizedNoisyResponseInstanceStimuli;

% Construct a QUEST threshold estimator estimate threshold
%
% The questThreshold estimator is associated with a psychometric function,
% by default qpPFWeibull.  It's important that the stimulus units be
% compatable with those expected by the psychometric function.  qpPFWeibull
% is coded in dB, which is confusing to many.  Probably better to work in
% log10 units, which can be done by passing qpPFWeibullLog as the PF.
if (~isfield(thresholdPara,'maxParamValue'))
    thresholdPara.maxParamValue = 1;
end
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

% Threshold estimation with QUEST+
% Get the initial stimulus contrast from QUEST+
[estimatorLogContrast, nextFlag] = estimator.nextStimulus();

% Loop over trials.
testedContrasts = [];

% Generate the NULL stimulus (zero contrast) if this is a TAFC task.
% Doing this here means we only do it once across the range of test
% contrasts studied.

if isTAFC
    nullContrast = 0.0;
    [theNullSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(nullContrast);
end

% Set up for saving diagnostic information if desired.
if ((isstruct(datasavePara)) && isfield(datasavePara, 'saveMRGCResponses') && (datasavePara.saveMRGCResponses) && ...
        (isfield(datasavePara, 'destDir')) && (ischar(datasavePara.destDir)) && ...
        (isfield(datasavePara, 'condExamined')))
    neuralEngineSaved = false;
else
    datasavePara.saveMRGCResponses = false;
end

% Dictionary to store the measured psychometric function which is returned
% to the user.  We also make and return containers that have the
% trial-by-trial stimulus alternative and trial-by-trial performance for
% each tested stimulus level
psychometricFunction = containers.Map();
trialByTrialStimulusAlternatives = containers.Map();
trialByTrialPerformance = containers.Map();
triailByTrialWhichResponses = containers.Map();

testCounter = 0;
stimCounter = 0;
while (nextFlag)
    % Convert log contrast -> contrast
    testContrast = thresholdPara.maxParamValue*(10 ^ estimatorLogContrast);
    logContrast = log10(testContrast);
    
    % Contrast label for pCorrect dictionary
    contrastLabel = sprintf('C = %2.4f%%', testContrast*100);
    if (verbose)
        fprintf('Testing parameter value %0.4g\n',testContrast)
    end
    
    % Have we already built the classifier for this contrast?
    testedIndex = find(testContrast == testedContrasts);
    if (isempty(testedIndex))
        % No. We need to train and predict.
        % 
        % Save contrast in list
        testedContrasts(numel(testedContrasts)+1) = testContrast;
        testedIndex = find(testContrast == testedContrasts);
        if (p.Results.useMetaContrast)
            if ~isTAFC
                % NWay: Generate the scenes for each alternative, at the test contrast
                % TAFC: Use the null scene and generate the test scene
                [theSceneSequences{testedIndex}, theSceneTemporalSupportSeconds] = ...
                    theSceneEngine.compute(testContrast);
            else
                % TAFC: Generate the test scene; we have the null from outside the
                % trial loop above.
                [theSceneSequence{testedIndex}, theSceneTemporalSupportSeconds] = ...
                    theSceneEngine.compute(testContrast);
                theSceneSequences{testedIndex} = {theNullSceneSequence,...
                    theSceneSequence{testedIndex}};
            end
        else
            if ~isTAFC
                % NWay: Generate the scenes for each alternative, at the test contrast
                % TAFC: Use the null scene and generate the test scene
                for oo = 1:length(theSceneEngine)
                    [theSceneSequences{testedIndex}{oo}, theSceneTemporalSupportSeconds] = ...
                        theSceneEngine{oo}.compute(testContrast);
                end
            else
                % TAFC: Generate the test scene; we have the null from outside the
                % trial loop above.
                [theSceneSequence{testedIndex}, theSceneTemporalSupportSeconds] = ...
                    theSceneEngine.compute(testContrast);
                theSceneSequences{testedIndex} = {theNullSceneSequence,...
                    theSceneSequence{testedIndex}};
            end
        end
        
        % Some diagnosis
        if (p.Results.extraVerbose)
            theWl = 400;
            theStim = 2; %TAFC: 2nd stim is the test
            theFrame = 1;
            Wl_vec = theSceneSequences{testedIndex}{theStim}{theFrame}.spectrum.wave;
            [~, index] = min(abs(Wl_vec - theWl));
            temp = theSceneSequences{testedIndex}{theStim}{theFrame}.data.photons(:,:,index);
            fprintf('At %d nm, frame %d, test scene %d, mean, min, max: %g, %g, %g\n',...
                Wl_vec(index),theFrame,testedIndex,mean(temp(:)),min(temp(:)),max(temp(:)));
            theWl = 550;
            [~, index] = min(abs(Wl_vec - theWl));
            temp = theSceneSequences{testedIndex}{theStim}{theFrame}.data.photons(:,:,index);
            fprintf('At %d nm, frame %d, test scene %d, mean, min, max: %g, %g, %g\n',...
                Wl_vec(index),theFrame,testedIndex,mean(temp(:)),min(temp(:)),max(temp(:)));
        end
        
        % Visualize the stimulus sequence.  This might crash for the
        % 4AFC version.
        if (visualizeStimulus)
            theStim = 2;
            if (length(theSceneEngine) > 1)
                sE = theSceneEngine{theStim};
            else
                sE = theSceneEngine;
            end
            sE.visualizeSceneSequence(theSceneSequences{testedIndex}{theStim},...
                theSceneTemporalSupportSeconds);
        end
        
        % Train classifier and get predicted
        % correct/incorrect predictions.  This function also computes the
        % neural responses needed to train and predict.
        %
        % Note that because the classifer engine is a handle class, we
        % need to use a copy method to do the caching.  Otherwise the
        % cached pointer will simply continue to point to the same
        % object, and be updated by future training.
        eStart = tic;

        % Find number of training and test fixational EM paths.  There are
        % four cases we allow.
        %   a) Both have 1 path, where staying still (empty object passed)
        %   is considered one path, as is one actual passed path.  In this
        %   case, we train on the training path and test on the test path.
        %   b) There is one training path and multiple test paths.  In this
        %   case, we train on the training path and test on each of the
        %   test paths, returning the average performance across all of
        %   them.
        %   c) There are multiple training paths and one test path.  We
        %   train on each of the training paths and test on the single test
        %   path, and return the average performance.
        %   d) There are N training paths and N test paths.  We train on
        %   each of the training paths and test on the corresponding test
        %   path, and return the average performance.
        %
        % Other cases are an error.
        if (~isempty(trainFixationalEMObj))
            nTrainFixationalEMPaths = size(trainFixationalEMObj.emPosArcMin,1);
            allTrainEMPaths = trainFixationalEMObj.emPosArcMin;
            trainFixationalEMObj.emPosMicrons = [];
            trainFixationalEMObj.emPos = [];
        else
            nTrainFixationalEMPaths = 1;
        end
        if (~isempty(testFixationalEMObj))
            nTestFixationalEMPaths = size(testFixationalEMObj.emPosArcMin,1);
            allTestEMPaths = testFixationalEMObj.emPosArcMin;
            testFixationalEMObj.emPosMicrons = [];
            testFixationalEMObj.emPos = [];
        else
            nTestFixationalEMPaths = 1;
        end

        % Check that passed EM conditions are supported.  Don't remove this
        % check unless you carefully fix the places below that count on it
        % passing.
        if (nTrainFixationalEMPaths ~= 1 & nTestFixationalEMPaths ~= 1)
            if (nTrainFixationalEMPaths ~= nTestFixationalEMPaths)
                error('Unsupported combination of traning and test EM paths');
            end
        end

        % Train the classifier for each EM path
        for eeTrain = 1:nTrainFixationalEMPaths
            % Create an EM object with just the eeTrain'th em path
            if (~isempty(trainFixationalEMObj))
                trainFixationalEMObj.emPosArcMin = allTrainEMPaths(eeTrain,:,:);
                trainFixationalEMObj.emPosMicrons = [];
                trainFixationalEMObj.emPos = [];
            end

            % Train the classifier for the eeTrain'th path.  Be sure to
            % copy the trained classifier, because it is a handle class.
            [~, tempTempClassifierEngine, ~, ~] = computePerformance(...
                theSceneSequences{testedIndex}, theSceneTemporalSupportSeconds,...
                classifierPara.nTrain, 0, theNeuralEngine,...
                classifierEngine, classifierPara.trainFlag, classifierPara.testFlag, ...
                'TAFC', isTAFC, 'saveResponses', datasavePara.saveMRGCResponses,...
                'visualizeAllComponents', visualizeAllComponents, ...
                'fixationalEM', trainFixationalEMObj, ...
                'useMetaContrast', p.Results.useMetaContrast ...
                );
            tempClassifierEngine{eeTrain} = tempTempClassifierEngine.copy;
        end

        % Restore training EM object
        if (~isempty(trainFixationalEMObj))
            trainFixationalEMObj.emPosArcMin = allTrainEMPaths;
            trainFixationalEMObj.emPosMicrons = [];
            trainFixationalEMObj.emPos = [];
        end

        % Bump some counters
        testCounter = testCounter + 1;
        stimCounter = stimCounter + 1;
        e = toc(eStart);
        if (verbose)
            fprintf('computeThreshold: Training test block %d took %0.1f secs\n',testCounter,e);
        end

        % Turn off nre detailed visualization if test counter exceeds the
        % specified number that we want it turned on for.
        if (stimCounter > maxVisualizedNoisyResponseInstanceStimuli)
            if (iscell(theNeuralEngine))
                for nn = 1:length(theNeuralEngine)
                    theNeuralEngine{nn}.visualizeEachCompute = false;
                end
            else
                theNeuralEngine.visualizeEachCompute = false;
            end
        end
        
        % Copy the trained classifiers so we hae them later
        eStart = tic;
        for eeTrain = 1:nTrainFixationalEMPaths
            theTrainedClassifierEngines{testedIndex}{eeTrain} = tempClassifierEngine{eeTrain};
        end
        e = toc(eStart);
        if (verbose)
                fprintf('computeThreshold: Copying classifer took %0.1f secs\n',e);
        end
        
        % Now predict. Handle each of the four fixations EM cases
        % separately.
        %
        % DHB: IT IS POSSIBLE WE COULD GET RID OF THIS SECTION AND AT THIS
        % POINT JUST DROP INTO THE PREDICT-ONLY CASE CODE BELOW, NOW THAT
        % WE HAVE SEPARATED TRAINING FROM PREDICTION IN THIS
        % TRAIN-AND-PREDICT BRANCH OF THE CODE.  BUT SAVING THAT FOR A
        % MOMENT WHEN EVERYTHING IS A BIT MORE STABLE.
        if (nTrainFixationalEMPaths == 1 && nTestFixationalEMPaths == 1)
            [predictions, ~, ~, whichAlternatives, whichResponses] = computePerformance(theSceneSequences{testedIndex}, ...
                theSceneTemporalSupportSeconds, 0, classifierPara.nTest, ...
                theNeuralEngine, theTrainedClassifierEngines{testedIndex}{1}, [], classifierPara.testFlag, ...
                'TAFC', isTAFC, ...
                'saveResponses', false, ...
                'visualizeAllComponents', false, ...
                'fixationalEM', testFixationalEMObj, ...
                'useMetaContrast', p.Results.useMetaContrast ...
                );

        elseif (nTrainFixationalEMPaths == 1 && nTestFixationalEMPaths > 1)
            % Set up accumulation variable
            predictions = [];
            whichAlternatives = [];
            responses = [];

             % Make sure numbers work out
            if (rem(classifierPara.nTest,nTestFixationalEMPaths) ~= 0)
                error('nTest must be an integer multiple of number of EM paths')
            end
            nTestPerEMPath = classifierPara.nTest/nTestFixationalEMPaths;

            % Loop over test EM paths getting predictions for each
            for eeTest = 1:nTestFixationalEMPaths
                % Set the test EM path for this time through the loop
                if (~isempty(testFixationalEMObj))
                    testFixationalEMObj.emPosArcMin = allTestEMPaths(eeTest,:,:);
                    testFixationalEMObj.emPosMicrons = [];
                    testFixationalEMObj.emPos = [];
                end

                % Predict
                [predictionsTemp, ~, ~, whichAlternativesTemp, whichResponsesTemp] = computePerformance(theSceneSequences{testedIndex}, ...
                    theSceneTemporalSupportSeconds, 0, nTestPerEMPath, ...
                    theNeuralEngine, theTrainedClassifierEngines{testedIndex}{1}, [], classifierPara.testFlag, ...
                    'TAFC', isTAFC, ...
                    'saveResponses', false, ...
                    'visualizeAllComponents', false, ...
                    'fixationalEM', testFixationalEMObj, ...
                    'useMetaContrast', p.Results.useMetaContrast ...
                );

                % Accumulate predictions over EM paths
                predictions = [predictions ; predictionsTemp];
                whichAlternatives = [whichAlternatives ; whichAlternativesTemp];
                whichResponses = [whichResponses ; whichResponsesTemp];
            end

            % Restore test EM object paths
            if (~isempty(testFixationalEMObj))
                testFixationalEMObj.emPosArcMin = allTestEMPaths;
                testFixationalEMObj.emPosMicrons = [];
                testFixationalEMObj.emPos = [];
            end

        elseif (nTrainFixationalEMPaths > 1 && nTestFixationalEMPaths == 1)
            % Set up accumulation variable
            predictions = [];
            whichAlternatives = [];
            responses = [];

             % Make sure numbers work out
            if (rem(classifierPara.nTest,nTrainFixationalEMPaths) ~= 0)
                error('nTest must be an integer miultiple of number of EM paths')
            end
            nTestPerEMPath = classifierPara.nTest/nTrainFixationalEMPaths;

            % Loop over training EM paths getting predictions for each
            for eeTrain = 1:nTrainFixationalEMPaths
                % Predict
                [predictionsTemp, ~, ~, whichAlternativesTemp, whichResponsesTemp] = computePerformance(theSceneSequences{testedIndex}, ...
                    theSceneTemporalSupportSeconds, 0, nTestPerEMPath, ...
                    theNeuralEngine, theTrainedClassifierEngines{testedIndex}{eeTrain}, [], classifierPara.testFlag, ...
                    'TAFC', isTAFC, ...
                    'saveResponses', false, ...
                    'visualizeAllComponents', false, ...
                    'fixationalEM', testFixationalEMObj, ...
                    'useMetaContrast', p.Results.useMetaContrast ...
                    );

                % Accumulate predictions over EM paths
                predictions = [predictions ; predictionsTemp];
                whichAlternatives = [whichAlternatives ; whichAlternativesTemp];
                whichResponses = [whichResponses ; whichResponsesTemp];

            end

        else
            % Set up accumulation variable
            predictions = [];
            whichAlternatives = [];
            responses = [];

            % Make sure numbers work out
            if (rem(classifierPara.nTest,nTestFixationalEMPaths) ~= 0)
                error('nTest must be an integer miultiple of number of EM paths')
            end
            nTestPerEMPath = classifierPara.nTest/nTestFixationalEMPaths;

            % Loop over paired train and test EM paths getting predictions for each
            for eeTest = 1:nTestFixationalEMPaths
                % Define for clarity of code
                eeTrain = eeTest;

                % Set the test EM path for this time through the loop
                if (~isempty(testFixationalEMObj))
                    testFixationalEMObj.emPosArcMin = allTestEMPaths(eeTest,:,:);
                    testFixationalEMObj.emPosMicrons = [];
                    testFixationalEMObj.emPos = [];
                end

                % Predict
                [predictionsTemp, ~, ~, whichAlternativesTemp, whichResponsesTemp] = computePerformance(theSceneSequences{testedIndex}, ...
                    theSceneTemporalSupportSeconds, 0, classifierPara.nTest, ...
                    theNeuralEngine, theTrainedClassifierEngines{testedIndex}{eeTrain}, [], classifierPara.testFlag, ...
                    'TAFC', isTAFC, ...
                    'saveResponses', false, ...
                    'visualizeAllComponents', false, ...
                    'fixationalEM', testFixationalEMObj, ...
                    'useMetaContrast', p.Results.useMetaContrast ...
                );

                % Accumulate predictions over EM paths
                predictions = [predictions ; predictionsTemp];
                whichAlternatives = [whichAlternatives ; whichAlternativesTemp];
                whichResponses = [whichResponses ; whichResponsesTemp];
            end

            % Restore train and test EM object paths
            if (~isempty(testFixationalEMObj))
                testFixationalEMObj.emPosArcMin = allTestEMPaths;
                testFixationalEMObj.emPosMicrons = [];
                testFixationalEMObj.emPos = [];
            end
        end

        % Update the psychometric function with data point for this
        % contrast level.  Also the trial-by-trial containers.
        psychometricFunction(contrastLabel) = mean(predictions);
        trialByTrialStimulusAlternatives(contrastLabel) = whichAlternatives;
        triailByTrialWhichResponses(contrastLabel) = whichResponses;
        trialByTrialPerformance(contrastLabel) = predictions';    
        if (verbose)
            fprintf('computeThreshold: Length of psychometric function %d, test counter %d\n',...
                length(psychometricFunction),testCounter);
        end

        % Save computed responses only the first time we test this contrast
        %
        % DHB: DOES THIS NEED TO BE UPDATEED FOR EACH EM PATH, OR IS IT GOOD AS
        % IS.  I AM NOT QUITE SURE WHAT THIS IS DOING.
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
        % Classifier is already trained, just get predictions.  The
        % whichAlternatives returned variable is only meaningful at present
        % for the rcePoisson classification engines.

        % Reality check
        if (estimator.validation & estimator.blocked)
             error('Should not be repeating any constrast for validation blocked method');
        end

        % Compute with trained classifier. We have already visualized for
        % this stimulus if we are visualizing, so we don't do it again here.
        eStart = tic;
        savevisualizeEachCompute = theNeuralEngine.visualizeEachCompute;
        theNeuralEngine.visualizeEachCompute = false;

        % Now predict. Handle each of the four fixations EM cases
        % separately.
        if (nTrainFixationalEMPaths == 1 && nTestFixationalEMPaths == 1)
            [predictions, ~, ~, whichAlternatives, whichResponses] = computePerformance(theSceneSequences{testedIndex}, ...
                theSceneTemporalSupportSeconds, 0, classifierPara.nTest, ...
                theNeuralEngine, theTrainedClassifierEngines{testedIndex}{1}, [], classifierPara.testFlag, ...
                'TAFC', isTAFC, ...
                'saveResponses', false, ...
                'visualizeAllComponents', false, ...
                'fixationalEM', testFixationalEMObj, ...
                'useMetaContrast', p.Results.useMetaContrast ...
                );

        elseif (nTrainFixationalEMPaths == 1 && nTestFixationalEMPaths > 1)
            % Set up accumulation variable
            predictions = [];
            whichAlternatives = [];
            responses = [];

             % Make sure numbers work out
            if (rem(classifierPara.nTest,nTestFixationalEMPaths) ~= 0)
                error('nTest must be an integer miultiple of number of EM paths')
            end
            nTestPerEMPath = classifierPara.nTest/nTestFixationalEMPaths;

            % Loop over test EM paths getting predictions for each
            for eeTest = 1:nTestFixationalEMPaths
                % Set the test EM path for this time through the loop
                if (~isempty(testFixationalEMObj))
                    testFixationalEMObj.emPosArcMin = allTestEMPaths(eeTest,:,:);
                    testFixationalEMObj.emPosMicrons = [];
                    testFixationalEMObj.emPos = [];
                end

                % Predict
                [predictionsTemp, ~, ~, whichAlternativesTemp, whichResponsesTemp] = computePerformance(theSceneSequences{testedIndex}, ...
                    theSceneTemporalSupportSeconds, 0, nTestPerEMPath, ...
                    theNeuralEngine, theTrainedClassifierEngines{testedIndex}{1}, [], classifierPara.testFlag, ...
                    'TAFC', isTAFC, ...
                    'saveResponses', false, ...
                    'visualizeAllComponents', false, ...
                    'fixationalEM', testFixationalEMObj, ...
                    'useMetaContrast', p.Results.useMetaContrast ...
                );

                % Accumulate predictions over EM paths
                predictions = [predictions ; predictionsTemp];
                whichAlternatives = [whichAlternatives ; whichAlternativesTemp];
                whichResponses = [whichResponses ; whichResponsesTemp];
            end

            % Restore test EM object paths
            if (~isempty(testFixationalEMObj))
                testFixationalEMObj.emPosArcMin = allTestEMPaths;
                testFixationalEMObj.emPosMicrons = [];
                testFixationalEMObj.emPos = [];
            end

        elseif (nTrainFixationalEMPaths > 1 && nTestFixationalEMPaths == 1)
            % Set up accumulation variable
            predictions = [];
            whichAlternatives = [];
            responses = [];

             % Make sure numbers work out
            if (rem(classifierPara.nTest,nTrainFixationalEMPaths) ~= 0)
                error('nTest must be an integer miultiple of number of EM paths')
            end
            nTestPerEMPath = classifierPara.nTest/nTrainFixationalEMPaths;

            % Loop over training EM paths getting predictions for each
            for eeTrain = 1:nTrainFixationalEMPaths
                % Predict
                [predictionsTemp, ~, ~, whichAlternativesTemp, whichResponsesTemp] = computePerformance(theSceneSequences{testedIndex}, ...
                    theSceneTemporalSupportSeconds, 0, nTestPerEMPath, ...
                    theNeuralEngine, theTrainedClassifierEngines{testedIndex}{eeTrain}, [], classifierPara.testFlag, ...
                    'TAFC', isTAFC, ...
                    'saveResponses', false, ...
                    'visualizeAllComponents', false, ...
                    'fixationalEM', testFixationalEMObj, ...
                    'useMetaContrast', p.Results.useMetaContrast ...
                    );

                % Accumulate predictions over EM paths
                predictions = [predictions ; predictionsTemp];
                whichAlternatives = [whichAlternatives ; whichAlternativesTemp];
                whichResponses = [whichResponses ; whichResponsesTemp];
            end

        else
            % Set up accumulation variable
            predictions = [];
            whichAlternatives = [];
            responses = [];

            % Make sure numbers work out
            if (rem(classifierPara.nTest,nTestFixationalEMPaths) ~= 0)
                error('nTest must be an integer miultiple of number of EM paths')
            end
            nTestPerEMPath = classifierPara.nTest/nTestFixationalEMPaths;

            % Loop over paired train and test EM paths getting predictions for each
            for eeTest = 1:nTestFixationalEMPaths
                % Define for clarity
                eeTrain = eeTest;

                if (~isempty(testFixationalEMObj))
                    testFixationalEMObj.emPosArcMin = allTestEMPaths(eeTest,:,:);
                    testFixationalEMObj.emPosMicrons = [];
                    testFixationalEMObj.emPos = [];
                end

                % Predict
                [predictionsTemp, ~, ~, whichAlternativesTemp, whichResponsesTemp] = computePerformance(theSceneSequences{testedIndex}, ...
                    theSceneTemporalSupportSeconds, 0, classifierPara.nTest, ...
                    theNeuralEngine, theTrainedClassifierEngines{testedIndex}{eeTrain}, [], classifierPara.testFlag, ...
                    'TAFC', isTAFC, ...
                    'saveResponses', false, ...
                    'visualizeAllComponents', false, ...
                    'fixationalEM', testFixationalEMObj, ...
                    'useMetaContrast', p.Results.useMetaContrast ...
                );

                % Accumulate predictions over EM paths
                predictions = [predictions ; predictionsTemp];
                whichAlternatives = [whichAlternatives ; whichAlternativesTemp];
                whichResponses = [whichResponses ; whichResponsesTemp];
            end

            % Restore train EM object paths
            if (~isempty(testFixationalEMObj))
                testFixationalEMObj.emPosArcMin = allTestEMPaths;
                testFixationalEMObj.emPosMicrons = [];
                testFixationalEMObj.emPos = [];
            end
        end
       
        % Restore the visualizeEachCompute flag.
        theNeuralEngine.visualizeEachCompute = savevisualizeEachCompute;

        testCounter = testCounter + 1;
        e = toc(eStart);
        if (verbose)
            fprintf('computeThreshold: Predicting test block %d no training took %0.1f secs\n',testCounter,e);
        end

        % Update the psychometric function with data point for this
        % contrast level. Also update the trial-by-trial containers.
        previousData = psychometricFunction(contrastLabel);
        currentData = cat(2,previousData,mean(predictions));
        psychometricFunction(contrastLabel) = currentData;

        prevTemp = trialByTrialStimulusAlternatives(contrastLabel);
        currentTemp = cat(1,prevTemp,whichAlternatives);
        trialByTrialStimulusAlternatives(contrastLabel) = currentTemp;

        prevTemp = triailByTrialWhichResponses(contrastLabel);
        currentTemp = cat(1,prevTemp,whichResponses);
        triailByTrialWhichResponses(contrastLabel) = currentTemp;

        prevTemp = trialByTrialPerformance(contrastLabel);
        currentTemp = cat(1,prevTemp,predictions');
        trialByTrialPerformance(contrastLabel) = currentTemp;

        if (verbose)
           fprintf('computeThreshold: Length of psychometric function %d, test counter %d\n',...
               length(psychometricFunction),testCounter);
        end
    end
    
    % Tell QUEST+ what we ran (how many trials at the given contrast) and
    % get next stimulus contrast to run.
    eStart = tic;
    if (estimator.validation)
        % Method of constant stimuli
        [estimatorLogContrast, nextFlag] = ...
            estimator.multiTrial(estimatorLogContrast * ones(classifierPara.nTest,1), predictions);
    else
        % Quest
        [estimatorLogContrast, nextFlag] = ...
          estimator.multiTrialQuestBlocked(estimatorLogContrast * ones(classifierPara.nTest,1), predictions);
    end

    e = toc(eStart);
    if (verbose)
            fprintf('computeThreshold: Updating took %0.1f secs\n',e);
    end
end 

if (verbose)
   fprintf('computeThreshold: Ran %d test levels of %d trials per block of tests\n',...
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
[estimatorLogThreshold, estimatorFittedPsychometricParams, estimatorThresholdDataOut] = ...
    estimator.thresholdMLE('showPlot', false, ...
    'thresholdCriterion', thresholdCriterion, 'returnData', true);
thresholdDataOut = estimatorThresholdDataOut;
thresholdDataOut.examinedContrasts = thresholdPara.maxParamValue*thresholdDataOut.examinedContrasts;
thresholdDataOut.examinedContrastsFit = thresholdPara.maxParamValue*thresholdDataOut.examinedContrastsFit;

% Take the scalar into account for the return values
threshold = thresholdPara.maxParamValue*10.^estimatorLogThreshold;
logThreshold = log10(threshold);

if (verbose)
    fprintf('Maximum likelihood fit parameters: %0.2f, %0.2f, %0.2f, %0.2f\n', ...
            estimatorFittedPsychometricParams(1), estimatorFittedPsychometricParams(2),...
            estimatorFittedPsychometricParams(3), estimatorFittedPsychometricParams(4));
    fprintf('Threshold (criterion proportion correct %0.4f): %0.2f (log10 units)\n', ...
        thresholdCriterion,logThreshold);
end

% Return the quest+ object wrapper for plotting and/or access to data
questObj = estimator;

end