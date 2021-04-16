function [threshold, questObj, psychometricFunction] = computeThresholdTAFC(theSceneEngine, theNeuralEngine, classifierEngine, ...
    classifierPara, thresholdPara, questEnginePara, varargin)
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
%
% Outputs:
%   threshold             - Estimated threshold value
%   questObj              - questThresholdEngine object, which
%                           contains information about all the trials run.
%   psychometricFunction  - Dictionary (indexed by contrast level) with the 
%                           psychometric function 
%
% Optional key/value pairs:
%   'beVerbose'           - Logical. Provide some printout? Default true.
%   'extraVerbose'        - Logical.  More detailed printout? Default false.
%   'visualizeStimulus'   - Logical. Provide stimulus visualization.
%                           Default false.
%   'visualizeAllComponents' - Logical. All component visualization.
%                           Default false.
%   'datasaveParameters'  - Parameters related to data saving. Default
%                           empty.
%
% See also:
%    t_spatialCSF, t_thresholdEngine, computePerformanceTAFC
%  

% History: 
%  10/23/20  dhb  Added commments.
%  12/04/20  npc  Added option to run method of constant stimuli. 
%                 Also added psychometricFunction return argument.

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
estDomain  = -thresholdPara.logThreshLimitLow : thresholdPara.logThreshLimitDelta : -thresholdPara.logThreshLimitHigh;
slopeRange = thresholdPara.slopeRangeLow: thresholdPara.slopeDelta : thresholdPara.slopeRangeHigh;

if (isfield(questEnginePara, 'employMethodOfConstantStimuli'))&&(questEnginePara.employMethodOfConstantStimuli)
    estimator = questThresholdEngine(...
        'validation', true, 'nRepeat', questEnginePara.nTest, ...
        'estDomain', estDomain, 'slopeRange', slopeRange);
else
    estimator = questThresholdEngine(...
        'minTrial', questEnginePara.minTrial, 'maxTrial', questEnginePara.maxTrial, ...
        'estDomain', estDomain, 'slopeRange', slopeRange, ...
        'numEstimator', questEnginePara.numEstimator, ...
        'stopCriterion', questEnginePara.stopCriterion);
end

% Generate the NULL stimulus (zero contrast)
nullContrast = 0.0;
[theNullSceneSequence, theSceneTemporalSupportSeconds] = theSceneEngine.compute(nullContrast);

% Some diagnosis
if (p.Results.p.Results.extraVerbose)
    theWl = 400;
    theFrame = 1;
    index = find(theNullSceneSequence{theFrame}.spectrum.wave == theWl);
    temp = theNullSceneSequence{theFrame}.data.photons(:,:,index);
    fprintf('At %d nm, frame %d, null scene mean, min, max: %g, %g, %g\n',theWl,theFrame,mean(temp(:)),min(temp(:)),max(temp(:)));
    theWl = 550;
    index = find(theNullSceneSequence{theFrame}.spectrum.wave == theWl);
    temp = theNullSceneSequence{theFrame}.data.photons(:,:,index);
    fprintf('At %d nm, frame %d, null scene mean, min, max: %g, %g, %g\n',theWl,theFrame,mean(temp(:)),min(temp(:)),max(temp(:)));
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
        testedContrasts = [testedContrasts testContrast];
        testedIndex = find(testContrast == testedContrasts);
        
        % Generate the TEST scene sequence for the given contrast
        [theTestSceneSequences{testedIndex}, ~] = theSceneEngine.compute(testContrast);
        
        % Some diagnosis
        if (p.Results.p.Results.extraVerbose)
            theWl = 400;
            theFrame = 1;
            index = find(theTestSceneSequences{testedIndex}{theFrame}.spectrum.wave == theWl);
            temp = theTestSceneSequences{testedIndex}{theFrame}.data.photons(:,:,index);
            fprintf('At %d nm, frame %d, test scene %d mean, min, max: %g, %g, %g\n',theWl,theFrame,testedIndex,mean(temp(:)),min(temp(:)),max(temp(:)));
            theWl = 550;
            index = find(theTestSceneSequences{testedIndex}{theFrame}.spectrum.wave == theWl);
            temp = theTestSceneSequences{testedIndex}{theFrame}.data.photons(:,:,index);
            fprintf('At %d nm, frame %d, test scene %d mean, min, max: %g, %g, %g\n',theWl,theFrame,testedIndex,mean(temp(:)),min(temp(:)),max(temp(:)));
        end
        
        % Visualize the drifting sequence
        if (visualizeStimulus)
            theSceneEngine.visualizeSceneSequence(theTestSceneSequences{testedIndex}, theSceneTemporalSupportSeconds);
        end
        
        % Train classifier for this TEST contrast and get predicted
        % correct/incorrect predictions.  This function also computes the
        % neural responses needed to train and predict.
        [predictions, theTrainedClassifierEngines{testedIndex}, responses] = computePerformanceTAFC(...
            theNullSceneSequence, theTestSceneSequences{testedIndex}, ...
            theSceneTemporalSupportSeconds, classifierPara.nTrain, classifierPara.nTest, ...
            theNeuralEngine, classifierEngine, classifierPara.trainFlag, classifierPara.testFlag, ...
            datasavePara.saveMRGCResponses, visualizeAllComponents);
        
        % Update the psychometric function with data point for this contrast level
        psychometricFunction(contrastLabel) = mean(predictions);
        
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
        [predictions, ~, ~] = computePerformanceTAFC(...
            theNullSceneSequence, theTestSceneSequences{testedIndex}, ...
            theSceneTemporalSupportSeconds, classifierPara.nTrain, classifierPara.nTest, ...
            theNeuralEngine, theTrainedClassifierEngines{testedIndex}, [], classifierPara.testFlag, ...
            false);
        
        % Update the psychometric function with data point for this contrast level
        previousData = psychometricFunction(contrastLabel);
        currentData = cat(2,previousData,mean(predictions));
        psychometricFunction(contrastLabel) = currentData;
    end
    
    
    % Tell QUEST+ what we ran (how many trials at the given contrast) and
    % get next stimulus contrast to run.
    [logContrast, nextFlag] = ...
        estimator.multiTrial(logContrast * ones(1, classifierPara.nTest), predictions);

end

% Return threshold value
[threshold, para] = estimator.thresholdMLE('showPlot', false);
if (beVerbose)
    fprintf('Maximum likelihood fit parameters: %0.2f, %0.2f, %0.2f, %0.2f\n', ...
        para(1), para(2), para(3), para(4));
end

% Return the quest+ object wrapper for plotting and/or access to data
questObj = estimator;

end