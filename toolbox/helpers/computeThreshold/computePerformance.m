function [predictions, theClassifierEngine, responses, whichAlternatives] = computePerformance(theScenes, ...
    temporalSupport, nTrain, nTest, theNeuralEngine, theClassifierEngine, trainNoiseFlag, testNoiseFlag, ...
    varargin)
% Compute performance of a classifier given different scenes, a neural 
% engine, and a classifier engine. This function is suitable for both TAFC 
% and N-alternative forced-choice tasks. 
%
% Syntax:
%    [predictions, theClassifierEngine, responses, whichAlternatives] = ...
%        computePerformance(theScenes, temporalSupport, nTrain, nTest, theNeuralEngine, ...
%       theClassifierEngine, trainNoiseFlag, testNoiseFlag, varargin)
%
% Description:
%     Train a classifier on a discrimination and report back a vector of
%     1's and 0's indicating correct and incorrect trials respectively.
%
%     This uses the ISETBioCSFGeneratorFramework and works because the uers
%     passes a set of objects with standardized API.  These describe the
%     two scenes to be discriminated, the neural pipeline that processes
%     these scenes, and the classifer.
%
% Inputs:
%     theScenes             - A collection of scenes (type: cell)
%                             If the task is TAFC, then the cell has a size
%                             of 1 x 2, the first for the null stimulus,
%                             and the second for the test stimulus.
%                             If the task is NWay_OneStimulusPerTrial, then
%                             the cell has a size of 1 x #alternative
%                             stimuli.
%     temporalSupport       - Temporal support vector (in seconds) for
%                             scene sequences.
%     nTrain                - Number of null and test response instances
%                             used in classifer training.  The two types of
%                             instances are paired and a nTrain TAFC task is
%                             simulated.
%     nTest                 - Number of null and test response instances
%                             used in classifer training.  The two types of
%                             instances are paired and nTest TAFC
%                             trials are simulated for evaluating
%                             performance.
%     theNeuralEngine       - @neuralResponseEngine object to compute
%                             neural responses.
%     theClassifierEngine   - @responseClassifierEngine object that
%                             implements observer decision model.  This is
%                             assumed untrained if trainedNoiseFlag
%                             contains a string, and trained if
%                             trainedNoiseFlag is empty.
%     trainNoiseFlag        - String.  Type of noise to be used in training
%                             the classifier. This flag are passed to
%                             theNeuralEngine to generate the training
%                             response instances. Typically either 'none' or
%                             'random' depending on whether the desired
%                             classifier is signal known exactly ('none')
%                             or signal known statistically ('random').
%     testNoiseFlag          - String. Type of noise to be used in
%                             evaluating performance. This flag are passed to
%                             theNeuralEngine to generate the test
%                             response instances. Typically 'random'.
%
% Outputs:
%     predictions            - Vector of 1's (correct) and 0's (incorrect)
%                              that gives trial by trial performance of the
%                              tested classifier in the TAFC task.
%                              Contains nTest entries.
%     theClassifierEngine    - Trained version of passed classifier object.
%
%     responses              - Neural responses computed
%     whichAlternatives      - Vector with integers that tells us which
%                              stimulus alternative was presented on each 
%                              of the stimulated trials.  This should have
%                              the same length as the predictions vector
%                              returned. This is currently only meaningful
%                              for classification engines whose names start
%                              with 'rcePoisson'.  
%
% Optional key/value pairs:
%     TAFC                  - logical. Whether this is a two-interval 
%                             forced-choice task or N-way one-stimulus-per-trial
%                             task. Default false.
%     useMetaContrast       - Passed arguments are for meta contrast setup
%     saveResponses         - Logical (default false). Whether to return the computed
%                             response instances
%     visualizeAllComponents - Logical (default false). Whether to visualize or not.
%     verbose               - Logical (default true). Print out stuff.
%     fixationalEM          - Empty (default) or a fixationalEM object
%                             that describes one eye movement path.  If
%                             the latter, this must have one position per
%                             frame of the passed scene sequence.
%
% See also
%   t_spatialCSF,  computeThreshold
%

% History:
%   10/23/20  dhb  Comments.
%   04/10/24  fh   Merged computePerformanceTAFC.m and
%                   computePerformanceNWay_OneStimulusPerTrial.m by adding 
%                   a key/pair pair specifying whether the task is TAFC or 
%                   NWay_OneStimulusPerTrial.

% Parse input
p = inputParser;
p.addParameter('TAFC', false, @islogical);
p.addParameter('useMetaContrast', false, @islogical);
p.addParameter('saveResponses',false, @islogical);
p.addParameter('visualizeAllComponents', false, @islogical);
p.addParameter('verbose', true, @islogical);
p.addParameter('fixationalEM', [], @(x)(isempty(x) || (isa(x,'fixationalEM'))));

parse(p, varargin{:});
isTAFC = p.Results.TAFC;
saveResponses = p.Results.saveResponses;
visualizeAllComponents = p.Results.visualizeAllComponents;
fixationalEMObj = p.Results.fixationalEM;

% Empty responses
responses = [];
if (p.Results.useMetaContrast && ~isTAFC)
    nScenes = length(theNeuralEngine);
else
    nScenes = length(theScenes);
end

% Train the classifier.
%
% If trainNoiseFlag is empty, then the passed classifier has already been trained
% and training is skipped.  Otherwise trainFlag is passed to the stimulus
% generation routine to indicate what type of noise ('none' or
% 'random') should be used in the training. Using 'random' inherits
% whatever noise model is in the nre.
if (~isempty(trainNoiseFlag))
    % Generate responses for training
    %
    % The responses are a 3 dimensional
    % matrix, with the dimensions indexing [instancesNum x mNeuralDim x tTimeBins].
    %   instancesNum   - number of response instances
    %   mNeuralDim     - dimension of neural response at one timepoint
    %   tTimeBins      - number of time points in stimulus sequence.
    % Note that if nTimeBins is 1, the last dimension is implicit,
    % following Matlab conventions.
    inSampleStimResponsesCell = cell(1,nScenes);
    for n = 1:nScenes
        if (p.Results.useMetaContrast && ~isTAFC)
            % Noise free response for nth alternative scene sequence
            [inSampleNoiseFreeStimResponsesCell{n}, ~] = theNeuralEngine{n}.computeNoiseFree(...
                theScenes, ...
                temporalSupport, ...
                'fixationalEM',fixationalEMObj);

            % Add noise (or not) to nth alternative responses
            [inSampleStimResponsesCell{n}, ~] = theNeuralEngine{n}.computeNoisyInstances( ...
                inSampleNoiseFreeStimResponsesCell{n}, ...
                temporalSupport, ...
                nTrain, ...
                trainNoiseFlag);
        else
            % Noise free response for nth alternative scene sequence
            [inSampleNoiseFreeStimResponsesCell{n}, ~] = theNeuralEngine.computeNoiseFree(...
                theScenes{n}, ...
                temporalSupport, ...
                'fixationalEM',fixationalEMObj);

            % Add noise (or not) to nth alternative responses
            [inSampleStimResponsesCell{n}, ~] = theNeuralEngine.computeNoisyInstances( ...
                inSampleNoiseFreeStimResponsesCell{n}, ...
                temporalSupport, ...
                nTrain, ...
                trainNoiseFlag);
        end
    end

    % If it's TAFC, massage the responses to be the concatenation
    % of the responses to the two stimuli that were actually
    % passed. This is because rcePoisson is set up to handle N
    % alternatives, one stimulus per trial and doing this turns the
    % TAFC case work into that format.
    %
    % If rcePoisson and not TAFC, then we just continue on with the cell
    % array we already have because we are explictly doing an N-alternative
    % one stimulus per trial task.
    %
    % We don't overwrite the original inSampleStimResonsesCell
    % array, because for visualization below it is convenient to
    % have the responses to the individual stimuli still available.
    if (isTAFC)
        % Concatenate them along the 2nd dimension (responses) in both
        % orders
        cat1 = cat(2, inSampleStimResponsesCell{1}, ...
            inSampleStimResponsesCell{2}); %[null, test]
        cat2 = cat(2, inSampleStimResponsesCell{2}, ...
            inSampleStimResponsesCell{1}); %[test, null]

        % Nicely put the concatenated cells back to the container
        inSampleStimResponsesMassagedCell = {cat1, cat2};
    else
        inSampleStimResponsesMassagedCell = inSampleStimResponsesCell;
    end

    % Train the classifier.
    theClassifierEngine.compute('train', inSampleStimResponsesMassagedCell,[]);
    clear inSampleStimResponsesMassagedCell

    % Visualization the cone excitation for selected stimulus
    if visualizeAllComponents
        % TAFC: 1st stim is the test; NWay: 1st stim is the correct stim
        theStim = 1; 
        visualizeConeResps(theNeuralEngine, inSampleStimResponsesCell, theStim);
    end

    % Save computed response instances.  These are the responses to the
    % individually passed stimuli, even when the actually classifier is
    % trained on concatenated resoponses
    if (saveResponses)
        responses.inSampleStimResponses = inSampleStimResponsesCell;
    end
end

% Predict using trained classifier.
%
% Generate stimulus for prediction, NULL stimulus.  The variable testFlag
% indicates what type of noise is used to generate the stimuli used for
% prediction.  Typically 'random'.

% Generate responses for prediction
%
% The responses are a 3 dimensional
% matrix, with the dimensions indexing [instancesNum x mNeuralDim x tTimeBins].
%   instancesNum   - number of response instances
%   mNeuralDim     - dimension of neural response at one timepoint
%   tTimeBins      - number of time points in stimulus sequence.
% Note that if nTimeBins is 1, the last dimension is implicit,
% following Matlab conventions.
%
% The way nTest gets passed to the neural engine is a little different for
% TAFC than true N-way.  We handle that here.
outSampleStimResponsesCell = cell(1,nScenes);
if isTAFC
    nTest_eachScene = nTest;
else
    nTest_eachScene = nTest/nScenes;
    assert(mod(nTest, nScenes) == 0, ['The number of test trials must be an',...
        ' integer multiple of the number of alternative choices']);
end
eStart = tic;
for n = 1:nScenes
    if (p.Results.useMetaContrast && ~isTAFC)
        % Noise free response for nth alternative scene sequence
        [outSampleNoiseFreeStimResponsesCell{n}, ~] = theNeuralEngine{n}.computeNoiseFree(...
            theScenes, ...
            temporalSupport, ...
            'fixationalEM',fixationalEMObj);

        % Add noise (or not) to nth alternative responses
        [outSampleStimResponsesCell{n}, ~] = theNeuralEngine{n}.computeNoisyInstances( ...
            outSampleNoiseFreeStimResponsesCell{n}, ...
            temporalSupport, ...
            nTest_eachScene, ...
            testNoiseFlag);
    else
        % Noise free response for nth alternative scene sequence
        [outSampleNoiseFreeStimResponsesCell{n}, ~] = theNeuralEngine.computeNoiseFree(...
            theScenes{n}, ...
            temporalSupport, ...
            'fixationalEM',fixationalEMObj);

        % Add noise (or not) to nth alternative responses
        [outSampleStimResponsesCell{n}, ~] = theNeuralEngine.computeNoisyInstances( ...
            outSampleNoiseFreeStimResponsesCell{n}, ...
            temporalSupport, ...
            nTest_eachScene, ...
            testNoiseFlag);
    end
end
e = toc(eStart);
if (p.Results.verbose)
    fprintf('computePerformance: Took %0.1f secs to generate test responses for all alternatives\n',e);
end

% If it's TAFC , massage the responses to be
% the concatenation of the responses to the two
% stimuli that were actually passed. This is because rcePoisson is set up to
% handle N alternatives, one stimulus per trial and doing this turns
% the TAFC case work into that format.
%
% If not TAFC, then we just continue on with the cell
% array we already have because we are explictly doing an N-alternative
% one stimulus per trial task.
%
% We don't overwrite the original outSampleStimResonsesCell
% array, because for visualization below it is convenient to
% have the responses to the individual stimuli still available.
if (isTAFC)
    % Concatenate them along the 2nd dimension (responses) in
    % null/test order.  We don't need to intermix test/cell as well
    % because the observer is not biased and doesn't care about the
    % order.
    outSampleStimResponsesMassaged = [];
    for nn = 1:nScenes
        outSampleStimResponsesMassaged = ...
            cat(2, outSampleStimResponsesMassaged, outSampleStimResponsesCell{nn});
    end
    whichAlternatives = ones(nTest_eachScene, 1);
else
    % Stack up the responses for each alternative
    outSampleStimResponsesMassaged = [];
    for nn = 1:nScenes
        outSampleStimResponsesMassaged = ...
            cat(1, outSampleStimResponsesMassaged, outSampleStimResponsesCell{nn});
    end
    whichAlternatives = repmat(1:nScenes,[nTest_eachScene, 1]);
    whichAlternatives = whichAlternatives(:); 
end

% Predict
dataOut = theClassifierEngine.compute('predict', outSampleStimResponsesMassaged, whichAlternatives);
clear outSampleStimResponsesMassaged

% Save computed response instances
if (saveResponses)
    responses.outSampleStimResponses = outSampleStimResponsesCell;
end
    
% Set return variable.  For each trial 0 means wrong and 1 means right.
% Taking mean(response) gives fraction correct.
predictions = dataOut.trialPredictions;

end

% Helper function for diagnostic visualizations.
function visualizeConeResps(theNeuralEngine, neuralResponses, selectedScene)

% If we've got a cone mosaic, visualize responses assuming they go with the
% cone mosaic
if (isfield(theNeuralEngine.neuralPipeline, 'coneMosaic'))
    Responses = neuralResponses{selectedScene};
    % Visualize the activation
    theNeuralEngine.neuralPipeline.coneMosaic.visualize('activation', ...
        squeeze(Responses), 'verticalActivationColorBarInside', true);

    % Also visualize the full absorptions density.  This is buggy because
    % it doesn't depend on the passed responses the way it should.
    figNo = 999;
    theNeuralEngine.neuralPipeline.coneMosaic.visualizeFullAbsorptionsDensity(figNo);

% If we don't have a cone mosaic but have an mRGC mosaic, visualize
% responses asumming they go with the mRGC mosaic
elseif (isfield(theNeuralEngine.neuralPipeline, 'mRGCmosaic'))
    theNeuralEngine.neuralPipeline.mRGCmosaic.visualizeResponses(...
        responseTemporalSupportSeconds, neuralResponses{selectedScene}, ...
        'stimulusTemporalSupportSeconds', temporalSupport,...
        'stimulusSceneSequence', testScene);
end

end