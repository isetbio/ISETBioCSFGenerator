function [threshold, questObj] = computeThresholdTAFC(theSceneEngine, theNeuralEngine, classifierEngine, ...
    classifierPara, thresholdPara, questEnginePara, visualizationPara, datasavePara)
% Compute contrast threshold for a given scene, neural response engine, and classifier engine
%
% Syntax:
%    [threshold, questObj] = computeThresholdTAFC(theSceneEngine, theNeuralEngine, classifierEngine, classifierPara, thresholdPara, questEnginePara)  
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
%   visualizationPara     - Parameter struct for visualizing different components
%   datasavePara          - Parameter struct for saving intermediate results
%
% Outputs:
%   threshold             - Estimated threshold value
%   questObj              - questThresholdEngine object, which
%                           contains information about all the trials run.
% Optional key/value pairs:
%
% See also:
%    t_spatialCSF, t_thresholdEngine, computePerformanceTAFC
%  

% History: 
%  10/23/20  dhb  Added commments.

% Construct a QUEST threshold estimator estimate threshold
estDomain  = -thresholdPara.logThreshLimitLow : thresholdPara.logThreshLimitDelta : -thresholdPara.logThreshLimitHigh;
slopeRange = thresholdPara.slopeRangeLow: thresholdPara.slopeDelta : thresholdPara.slopeRangeHigh;
estimator = ...
    questThresholdEngine('minTrial', questEnginePara.minTrial, 'maxTrial', questEnginePara.maxTrial, ...
                                                 'estDomain', estDomain, 'slopeRange', slopeRange, ...
                                                 'numEstimator', questEnginePara.numEstimator, 'stopCriterion', questEnginePara.stopCriterion);

% Generate the NULL stimulus (zero contrast)
nullContrast = 0.0;
[theNullSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(nullContrast);

% Threshold estimation with QUEST+
% Get the initial stimulus contrast from QUEST+
[logContrast, nextFlag] = estimator.nextStimulus();

% Loop over trials.
testedContrasts = [];

if (datasavePara.saveMRGCResponses)
    neuralEngineSaved = false;
end

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
        
        % Visualize the drifting sequence
        if (visualizationPara.visualizeStimulus)
            theSceneEngine.visualizeSceneSequence(theTestSceneSequences{testedIndex}, theSceneTemporalSupportSeconds);
        end
        
        % Train classifier for this TEST contrast and get predicted
        % correct/incorrect predictions.  This function also computes the
        % neural responses needed to train and predict.
        [predictions, theTrainedClassifierEngines{testedIndex}, responses] = computePerformanceTAFC(...
            theNullSceneSequence, theTestSceneSequences{testedIndex}, ...
            theSceneTemporalSupportSeconds, classifierPara.nTrain, classifierPara.nTest, ...
            theNeuralEngine, classifierEngine, classifierPara.trainFlag, classifierPara.testFlag, ...
            datasavePara.saveMRGCResponses, visualizationPara.visualizeAllComponents);
        
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
        % Classifier is already trained, just get predictions
        [predictions, ~, responses] = computePerformanceTAFC(...
            theNullSceneSequence, theTestSceneSequences{testedIndex}, ...
            theSceneTemporalSupportSeconds, classifierPara.nTrain, classifierPara.nTest, ...
            theNeuralEngine, theTrainedClassifierEngines{testedIndex}, [], classifierPara.testFlag, ...
            false);
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

% Return the quest+ object wrapper for plotting and/or access to data
questObj = estimator;

end