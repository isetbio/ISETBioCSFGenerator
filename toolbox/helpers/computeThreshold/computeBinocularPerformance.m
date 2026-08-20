function [predictions, theClassifierEngine, responses, whichAlternatives, whichResponses, whichMetaDataLeftEye, whichMetaDataRightEye] = ...
    computeBinocularPerformance(theLeftEyeScenes, theRightEyeScenes, temporalSupport, nTrain, nTest, ...
    theNeuralEngineLeftEye, theNeuralEngineRightEye, theClassifierEngine, trainNoiseFlag, testNoiseFlag, ...
    varargin)
% Compute performance of a classifier given different scenes, a neural
% engine, and a classifier engine. This function is suitable for both TAFC
% and N-alternative forced-choice tasks.
%
% Syntax:
%    [predictions, theClassifierEngine, responses, whichAlternatives] = ...
%       computeBinocularPerformance(theLeftEyeScenes, theRightEyeScenes, temporalSupport, nTrain, nTest, ...
%       theNeuralEngineLeftEye, theNeuralEngineRightEye,
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
%     theLeftEyeScenes        - A collection of scenes for the left eye(type: cell)
%                               If the task is TAFC, then the cell has a size
%                               of 1 x 2, the first for the null stimulus,
%                               and the second for the test stimulus.
%                               If the task is NWay_OneStimulusPerTrial, then
%                               the cell has a size of 1 x #alternative
%                               stimuli.
%     theRightEyeScenes       - A collection of scenes for the right eye(type: cell)
%                               If the task is TAFC, then the cell has a size
%                               of 1 x 2, the first for the null stimulus,
%                               and the second for the test stimulus.
%                               If the task is NWay_OneStimulusPerTrial, then
%                               the cell has a size of 1 x #alternative
%                               stimuli.
%     temporalSupport         - Temporal support vector (in seconds) for
%                               scene sequences.
%     nTrain                  - Number of null and test response instances
%                               used in classifer training.  The two types of
%                               instances are paired and a nTrain TAFC task is
%                               simulated. Training is skipped if the
%                               trainNoiseFlag (see below) is empty, or if
%                               this is 0.
%     nTest                   - Number of null and test response instances
%                               used in classifer training.  The two types of
%                               instances are paired and nTest TAFC
%                               trials are simulated for evaluating
%                               performance.  If this is passed as 0, then no
%                               testing happens.
%     theNeuralEngineLeftEye  - @neuralResponseEngine object to compute
%                               neural responses (left eye).
%     theNeuralEngineRightEye - @neuralResponseEngine object to compute
%                               neural responses (right eye).
%     theClassifierEngine     - @responseClassifierEngine object that
%                               implements observer decision model.  This is
%                               assumed untrained if trainedNoiseFlag
%                               contains a string, and trained if
%                               trainedNoiseFlag is empty.
%     trainNoiseFlag          - String.  Type of noise to be used in training
%                               the classifier. This flag are passed to
%                               theNeuralEngine to generate the training
%                               response instances. Typically either 'none' or
%                               'random' depending on whether the desired
%                               classifier is signal known exactly ('none')
%                               or signal known statistically ('random').  If
%                               this is empty, the classifier is assumed
%                               trained.
%     testNoiseFlag            - String. Type of noise to be used in
%                               evaluating performance. This flag are passed to
%                               theNeuralEngine to generate the test
%                               response instances. Typically 'random'.
%
% Outputs:
%     predictions            - Vector of 1's (correct) and 0's (incorrect)
%                              that gives trial by trial performance of the
%                              tested classifier in the TAFC task.
%                              Contains nTest entries.  Returned as empty
%                              if nTest == 0.
%     theClassifierEngine    - Trained version of passed classifier object.
%     responses              - Neural responses computed for testing,
%                              Returned as empty if nTest == 0;
%     whichAlternatives      - Vector with integers that tell us which
%                              stimulus alternative was presented on each
%                              of the stimulated test trials.  This should have
%                              the same length as the predictions vector
%                              returned.
%                              If nTest == 0, this is returned as empty.
%     whichResponses        -  Vector with integers that tell us which
%                              alternative the classifier picked on each trial.
%                              Currently only meaningful for rcePoisson and
%                              rceTemplateDistance classifiers.
%                              If nTest == 0, this is returned as empty.
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
p.addParameter('visualizeBinocularSummation', true, @islogical);
p.addParameter('verbose', true, @islogical);
p.addParameter('fixationalEM', [], @(x)(isempty(x) || (isa(x,'fixationalEM'))));
p.addParameter('conditionLabel', '', @ischar);
p.addParameter('binocularSummationPipeline', [], @(x)(isempty(x) || (isstruct(x))));

parse(p, varargin{:});
isTAFC = p.Results.TAFC;
saveResponses = p.Results.saveResponses;
visualizeAllComponents = p.Results.visualizeAllComponents;
visualizeBinocularSummation = p.Results.visualizeBinocularSummation;
fixationalEMObj = p.Results.fixationalEM;
conditionLabel = p.Results.conditionLabel;
binocularSummationPipeline = p.Results.binocularSummationPipeline;


% Empty responses
responses = [];
if (p.Results.useMetaContrast && ~isTAFC)
    assert(length(theNeuralEngineLeftEye) == length(theNeuralEngineRightEye), 'Left and Right eye neural engines must have equal numerosity')
    nScenes = length(theNeuralEngineLeftEye);
else
    assert(length(theLeftEyeScenes) == length(theRightEyeScenes), 'Left and Right eye scenes must have equal numerosity')
    nScenes = length(theLeftEyeScenes);
end




% Train the classifier.
%
% If trainNoiseFlag is empty, then the passed classifier has already been trained
% and training is skipped.  Otherwise trainFlag is passed to the stimulus
% generation routine to indicate what type of noise ('none' or
% 'random') should be used in the training. Using 'random' inherits
% whatever noise model is in the nre.
if (~isempty(trainNoiseFlag) & nTrain ~= 0)
    % Generate responses for training
    %
    % The responses are a 3 dimensional
    % matrix, with the dimensions indexing [instancesNum x mNeuralDim x tTimeBins].
    %   instancesNum   - number of response instances
    %   mNeuralDim     - dimension of neural response at one timepoint
    %   tTimeBins      - number of time points in stimulus sequence.
    % Note that if nTimeBins is 1, the last dimension is implicit,
    % following Matlab conventions.

    % Left eye in-sample responses
    [inSampleStimResponsesCellLeftEye, neuralResponseEngineMetaDataLeftEye, ...
     neuralResponseTemporalSupport, theNeuralEngineLeftEye] = ...
        computeInSampleResponses(nTrain, nScenes, p.Results.useMetaContrast, isTAFC, ...
            theNeuralEngineLeftEye, fixationalEMObj, trainNoiseFlag, ...
            theLeftEyeScenes, temporalSupport, ... 
            visualizeAllComponents);

    % Right eye in-sample responses
    [inSampleStimResponsesCellRightEye, neuralResponseEngineMetaDataRightEye, ...
     neuralResponseTemporalSupport, theNeuralEngineRightEye] = ...
        computeInSampleResponses(nTrain, nScenes, p.Results.useMetaContrast, isTAFC, ...
            theNeuralEngineRightEye, fixationalEMObj, trainNoiseFlag, ...
             theRightEyeScenes, temporalSupport, ... 
             visualizeAllComponents);


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

        % Compute binocular responses to the NULL stimulus
        [inSampleStimResponsesCell{1}, theLinearBinocularSummationNULLresponses] = applyBinocularSummationPipeline(...
                inSampleStimResponsesCellLeftEye{1}, ...
                inSampleStimResponsesCellRightEye{1}, ...
                binocularSummationPipeline);

        % Compute binocular responses to the TEST stimulus
        [inSampleStimResponsesCell{2}, theLinearBinocularSummationTESTresponses] = applyBinocularSummationPipeline(...
                inSampleStimResponsesCellLeftEye{2}, ...
                inSampleStimResponsesCellRightEye{2}, ...
                binocularSummationPipeline);

        if (visualizeBinocularSummation)

            visualizeBinocularSummationKernelsAndSignals(...
                inSampleStimResponsesCellLeftEye{2}, ...
                inSampleStimResponsesCellRightEye{2}, ...
                inSampleStimResponsesCell{2}, ...
                theLinearBinocularSummationTESTresponses, ...
                binocularSummationPipeline, ...
                conditionLabel);
        end

        cat1 = cat(2, inSampleStimResponsesCell{1}, ...
            inSampleStimResponsesCell{2}); %[null, test]
        cat2 = cat(2, inSampleStimResponsesCell{2}, ...
            inSampleStimResponsesCell{1}); %[test, null]

        % Nicely put the concatenated cells back to the container
        inSampleStimResponsesMassagedCell = {cat1, cat2};
    else
        % Concatenate Left and Right eye responses
        inSampleStimResponsesCell = applyBinocularSummationPipeline(...
                inSampleStimResponsesCellLeftEye, ...
                inSampleStimResponsesCellRightEye, ...
                binocularSummationPipeline);

        inSampleStimResponsesMassagedCell = inSampleStimResponsesCell;
    end

    % Train the classifier.
    theClassifierEngine.compute('train', inSampleStimResponsesMassagedCell,[]);
    clear inSampleStimResponsesMassagedCell

    % Visualization the cone excitation for selected stimulus
    if visualizeAllComponents
        % TAFC: 1st stim is the test; NWay: 1st stim is the correct stim
        theStim = 1;
        visualizeConeResps(theNeuralEngineLeftEye, inSampleStimResponsesCellLeftEye, theStim);
        visualizeConeResps(theNeuralEngineRightEye, inSampleStimResponsesCellRightEye, theStim);
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


if (nTest ~= 0)
    
    % Left eye out-of-sample responses
    [outSampleStimResponsesCellLeftEye, neuralResponseEngineMetaDataLeftEye, nTest_eachScene, ...
     neuralResponseTemporalSupport, theNeuralEngineLeftEye] =  computeOutOfSampleResponses(...
        nTest, nScenes, p.Results.useMetaContrast, isTAFC,  ...
        theNeuralEngineLeftEye, fixationalEMObj, testNoiseFlag, ...
        theLeftEyeScenes, temporalSupport, ... 
        visualizeAllComponents, p.Results.verbose);


    % Right eye out-of-sample responses
    [outSampleStimResponsesCellRightEye, neuralResponseEngineMetaDataRightEye, nTest_eachScene,  ...
     neuralResponseTemporalSupport, theNeuralEngineRightEye] = computeOutOfSampleResponses(...
        nTest, nScenes, p.Results.useMetaContrast, isTAFC,  ...
        theNeuralEngineRightEye, fixationalEMObj, testNoiseFlag, ...
        theRightEyeScenes, temporalSupport, ... 
        visualizeAllComponents, p.Results.verbose);


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
        outSampleStimResponsesCell = cell(1, nScenes);
        for nn = 1:nScenes

            % Binocular summation of Left and Right eye responses
            [outSampleStimResponsesCell{nn},  theLinearBinocularSummationTESTresponses]= applyBinocularSummationPipeline(...
                outSampleStimResponsesCellLeftEye{nn}, ...
                outSampleStimResponsesCellRightEye{nn}, ...
                binocularSummationPipeline);

            if (visualizeBinocularSummation)
                visualizeBinocularSummationKernelsAndSignals(...
                    outSampleStimResponsesCellLeftEye{nn}, ...
                    outSampleStimResponsesCellRightEye{nn}, ...
                    outSampleStimResponsesCell{nn}, ...
                    theLinearBinocularSummationTESTresponses, ...
                    binocularSummationPipeline, ...
                    sprintf('%s\n-- scene: %d of %d --', conditionLabel, nn, nScenes));
            end

            outSampleStimResponsesMassaged = ...
                cat(2, outSampleStimResponsesMassaged, outSampleStimResponsesCell{nn});
        end
        whichAlternatives = ones(nTest_eachScene, 1);
    else
        % Stack up the responses for each alternative
        outSampleStimResponsesMassaged = [];
        outSampleStimResponsesCell = cell(1, nScenes);
        for nn = 1:nScenes

            % Concatenate Left and Right eye responses
            outSampleStimResponsesCell{nn} = applyBinocularSummationPipeline(...
                outSampleStimResponsesCellLeftEye{nn}, ...
                outSampleStimResponsesCellRightEye{nn}, ...
                binocularSummationPipeline);

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

    % Set whichResponses return variable
    whichResponses = dataOut.whichAlternativePredicted(:);

    % Set whichMetaData return variable
    whichMetaDataLeftEye = neuralResponseEngineMetaDataLeftEye;
    whichMetaDataRightEye = neuralResponseEngineMetaDataRightEye;

else
    % Not testing, only training
    predictions = [];
    responses = [];
    whichAlternatives = [];
    whichMetaData = [];
end

end



function visualizeBinocularSummationKernelsAndSignals(...
                theLeftEyeConeMosaicResponses, ...
                theRightEyeConeMosaicResponses, ...
                theBinocularSummationResponses, ...
                theLinearBinocularSummationResponses, ...
                theBinocularSummationPipeline, ...
                conditionLabel)

    
    if (~isfield(theBinocularSummationPipeline, 'components')) || ...
       ((isfield(theBinocularSummationPipeline, 'components'))&&(isempty(theBinocularSummationPipeline.components)))
        fprintf('Binocular pipeline of type has no components to visualize\n', theBinocularSummationPipeline.type)
        return;
    end

    
    theLeftConeMosaic = [];
    theRightConeMosaic = [];
    if (isfield(theBinocularSummationPipeline.components, 'leftEyeConeMosaic'))
        theLeftConeMosaic = theBinocularSummationPipeline.components.leftEyeConeMosaic;
    end

    if (isfield(theBinocularSummationPipeline.components, 'rightEyeConeMosaic'))
        theRightConeMosaic = theBinocularSummationPipeline.components.rightEyeConeMosaic;
    end


    [timeBinsNum, theRightMosaicConesNum] = size(theBinocularSummationPipeline.components.rightConeMosaicSummationWeights);
    [nReps, theRightMosaicConesNumCheck] = size(theRightEyeConeMosaicResponses);
    
    [timeBinsNumCheck, theLeftMosaicConesNum] = size(theBinocularSummationPipeline.components.leftConeMosaicSummationWeights);
    [nRepsCheck, theLeftMosaicConesNumCheck] = size(theLeftEyeConeMosaicResponses);

    assert(nReps == nRepsCheck, 'inconsistent cones num between LEFT and RIGHT response reps');
    assert(theRightMosaicConesNum == theRightMosaicConesNumCheck, 'inconsistent cones num between RIGHT cone mosaic weights and responses');
    assert(theLeftMosaicConesNum == theLeftMosaicConesNumCheck, 'inconsistent cones num between LEFT cone mosaic weights and responses');
    assert(timeBinsNum==timeBinsNumCheck, 'inconsistent time bins num between LEFT cone mosaic weights and responses');

    % Reshape
    rightEyeRF = reshape(theBinocularSummationPipeline.components.rightConeMosaicSummationWeights, [1 timeBinsNum theRightMosaicConesNum]);
    rightEyeInputSignal = reshape(theRightEyeConeMosaicResponses, [nReps timeBinsNum theRightMosaicConesNum]);

    leftEyeRF = reshape(theBinocularSummationPipeline.components.leftConeMosaicSummationWeights, [1 timeBinsNum theLeftMosaicConesNum]);
    leftEyeInputSignal = reshape(theLeftEyeConeMosaicResponses, [nReps timeBinsNum theLeftMosaicConesNum]);

    rfColormap = brewermap(256, '*reds');
    rfColormap = cat(1, rfColormap, brewermap(256, 'blues'));
    rfColormap = rfColormap(end:-1:1,:);

    activationColormap = brewermap(256, 'greys');

    maxRF = max([max(abs(rightEyeRF(:))) max(abs(leftEyeRF(:)))]);
    maxActivation = max([prctile(abs(rightEyeInputSignal(:)), 95) prctile(abs(leftEyeInputSignal(:)),95)]);

    [~,repResultingInMinLinearActivation] = min(theLinearBinocularSummationResponses(:));
    [~,repResultingInMaxLinearActivation] = max(theLinearBinocularSummationResponses(:));
    [~,repResultingInZeroLinearActivation] = min(abs(theLinearBinocularSummationResponses(:)));

    visualizeReps = unique([...
        repResultingInMinLinearActivation ...
        repResultingInZeroLinearActivation ...
        repResultingInMaxLinearActivation]);


    mosaicSizeDegsRound = round(max(theLeftConeMosaic.sizeDegs)*10)/10;
    domainVisualizationLimits(1:2) = theLeftConeMosaic.eccentricityDegs(1) + 0.5*[-1 1]*mosaicSizeDegsRound;
    domainVisualizationLimits(3:4) = theLeftConeMosaic.eccentricityDegs(2) + 0.5*[-1 1]*mosaicSizeDegsRound;
    domainVisualizationTicks = struct(...
        'x', round(10*(theLeftConeMosaic.eccentricityDegs(1) + 0.5*mosaicSizeDegsRound*[-1 0 1]))/10, ...
        'y', round(10*(theLeftConeMosaic.eccentricityDegs(2) + 0.5*mosaicSizeDegsRound*[-1 0 1]))/10);

    hFig = figure(99); clf;
    set(hFig, 'Position', [10 10 1500 850]);

    for idx = 1:numel(visualizeReps)

        iRep = visualizeReps(idx);

        if (~isempty(theLeftConeMosaic)) 
        
             ax = subplot(2,3,1);
             theLeftConeMosaic.visualize(...
                 'figureHandle', hFig, ...
                 'axesHandle', ax, ...
                 'activation', leftEyeInputSignal(iRep,:,:), ...
                 'activationColormap', activationColormap, ...
                 'activationRange', maxActivation*[-1 1], ...
                 'domainVisualizationLimits', domainVisualizationLimits, ...
                 'domainVisualizationTicks', domainVisualizationTicks, ...
                 'plotTitle', sprintf('cone mosaic modulation (LE)\n(instance: %d)', iRep));
    
             ax = subplot(2,3,2);
             theLeftConeMosaic.visualize(...
                 'figureHandle', hFig, ...
                 'axesHandle', ax, ...
                 'activation', leftEyeRF, ...
                 'activationColormap', rfColormap, ...
                 'activationRange', maxRF*[-1 1], ...
                 'domainVisualizationLimits', domainVisualizationLimits, ...
                 'domainVisualizationTicks', domainVisualizationTicks, ...
                 'plotTitle', 'V1 cone pooling weights (LE)');
    
         end

         if (~isempty(theRightConeMosaic)) 
    
             ax = subplot(2,3,4);
             theRightConeMosaic.visualize(...
                 'figureHandle', hFig, ...
                 'axesHandle', ax, ...
                 'activation', rightEyeInputSignal(iRep,:,:), ...
                 'activationColormap', activationColormap, ...
                 'activationRange', maxActivation*[-1 1], ...
                 'domainVisualizationLimits', domainVisualizationLimits, ...
                 'domainVisualizationTicks', domainVisualizationTicks, ...
                 'plotTitle', sprintf('cone mosaic modulation (RE)\n(instance: %d)', iRep));

             ax = subplot(2,3,5);
             theRightConeMosaic.visualize(...
                 'figureHandle', hFig, ...
                 'axesHandle', ax, ...
                 'activation', rightEyeRF, ...
                 'activationColormap', rfColormap, ...
                 'activationRange', maxRF*[-1 1], ...
                 'domainVisualizationLimits', domainVisualizationLimits, ...
                 'domainVisualizationTicks', domainVisualizationTicks, ...
                 'plotTitle', 'V1 cone pooling weights (RE)');
         end
     

         % The light weighting functions of the cone pooling

         theRightConeMosaicROI = regionOfInterest(...
            'geometryStruct', struct(...
                'units', 'degs', ...
                'shape', 'rect', ...
                'center', theRightConeMosaic.eccentricityDegs, ...
                'width', theRightConeMosaic.sizeDegs(1), ...
                'height', 0.1, ...
                'rotation', 0.0...
            ));

         theLeftConeMosaicROI = regionOfInterest(...
            'geometryStruct', struct(...
                'units', 'degs', ...
                'shape', 'rect', ...
                'center', theLeftConeMosaic.eccentricityDegs, ...
                'width', theLeftConeMosaic.sizeDegs(1), ...
                'height', 0.1, ...
                'rotation', 0.0...
            ));


         visualizedConeIndices = theRightConeMosaicROI.indicesOfPointsInside(theRightConeMosaic.coneRFpositionsDegs);
         theRightEyeXcoords = squeeze(theRightConeMosaic.coneRFpositionsDegs(visualizedConeIndices,1));
         theRightEyeConePoolingWeights = squeeze(rightEyeRF(1,1,visualizedConeIndices));

         visualizedConeIndices = theLeftConeMosaicROI.indicesOfPointsInside(theLeftConeMosaic.coneRFpositionsDegs);
         theLeftEyeXcoords = squeeze(theLeftConeMosaic.coneRFpositionsDegs(visualizedConeIndices,1));
         theLeftEyeConePoolingWeights = squeeze(leftEyeRF(1,1,visualizedConeIndices));

         ax = subplot(2,3,3);
         plot(ax, theRightEyeXcoords, theRightEyeConePoolingWeights, 'r.');
         hold (ax, 'on');
         plot(ax, theLeftEyeXcoords, theLeftEyeConePoolingWeights, 'b.');
         axis(ax, 'square')
         set(ax, 'YLim', [-1 1], 'XLim',  domainVisualizationLimits(1:2), 'XTick', domainVisualizationTicks.x);
         xlabel(ax, 'space, x (degs)');
         ylabel(ax, 'pooling amplitude');
         

         ax = subplot(2,3,6);
         plot(ax, theLinearBinocularSummationResponses, theBinocularSummationResponses, 'b.');
         hold (ax, 'on');
         plot(ax, theLinearBinocularSummationResponses(iRep), theBinocularSummationResponses(iRep), 'ro');

         axis(ax, 'equal');
         axis(ax, 'square');
         xlabel(ax, 'linear binocular summation response');
         ylabel(ax, 'non-linear binocular summation response');
         title(ax, sprintf('activation function\n%s\n(%d instances, %d time bins)', ...
             conditionLabel, ...
             size(theBinocularSummationResponses,1), size(theBinocularSummationResponses,2)));

         drawnow;
    end % iRep

end


function [theBinocularResponse, theLinearBinocularSummationResponse] = applyBinocularSummationPipeline(...
    leftEyeResponses, rightEyeResponses, ...
    binocularSummationPipeline)

    switch (binocularSummationPipeline.type)
        case 'left + right eye response concatenation'
            theBinocularResponse = cat(2, leftEyeResponses, rightEyeResponses);
            theLinearBinocularSummationResponse = [];

        case 'zeroDisparityTunedSimpleCell'

            [timeBinsNum, theRightMosaicConesNum] = size(binocularSummationPipeline.components.rightConeMosaicSummationWeights);
            [nReps, theRightMosaicConesNumCheck] = size(rightEyeResponses);
            
            [timeBinsNumCheck, theLeftMosaicConesNum] = size(binocularSummationPipeline.components.leftConeMosaicSummationWeights);
            [nRepsCheck, theLeftMosaicConesNumCheck] = size(leftEyeResponses);

            assert(nReps == nRepsCheck, 'inconsistent cones num between LEFT and RIGHT response reps');
            assert(theRightMosaicConesNum == theRightMosaicConesNumCheck, 'inconsistent cones num between RIGHT cone mosaic weights and responses');
            assert(theLeftMosaicConesNum == theLeftMosaicConesNumCheck, 'inconsistent cones num between LEFT cone mosaic weights and responses');
            assert(timeBinsNum==timeBinsNumCheck, 'inconsistent time bins num between LEFT cone mosaic weights and responses');

            % Reshape
            rightEyeRF = reshape(binocularSummationPipeline.components.rightConeMosaicSummationWeights, [1 timeBinsNum theRightMosaicConesNum]);
            rightEyeInputSignal = reshape(rightEyeResponses, [nReps timeBinsNum theRightMosaicConesNum]);

            leftEyeRF = reshape(binocularSummationPipeline.components.leftConeMosaicSummationWeights, [1 timeBinsNum theLeftMosaicConesNum]);
            leftEyeInputSignal = reshape(leftEyeResponses, [nReps timeBinsNum theLeftMosaicConesNum]);

            % Preallocate memory
            theLinearBinocularSummationResponse = zeros(nReps, timeBinsNum, 1);

            % Compute the limear binocular summation of LE and RE signals weighted by
            % the LE and RE receptive fields
            parfor iRep = 1:nReps
                theLinearBinocularSummationResponse(iRep,:) = ...
                    dot(leftEyeInputSignal(iRep, 1:timeBinsNum, :),  leftEyeRF(1, 1:timeBinsNum,:),  3) + ...
                    dot(rightEyeInputSignal(iRep, 1:timeBinsNum, :), rightEyeRF(1, 1:timeBinsNum,:), 3);
            end
        
            % Apply the instantaneous static non-linearity
            theBinocularResponse = binocularSummationPipeline.components.activationFunction(...
                theLinearBinocularSummationResponse, ...
                binocularSummationPipeline.components.activationFunctionParams);

        otherwise
            error('Unknown binocular summation pipeline type: ''%s''.', binocularSummationPipeline.type)
    end

end



function [outSampleStimResponsesCell, neuralResponseEngineMetaData, nTest_eachScene, ...
        neuralResponseTemporalSupport, theNeuralEngine] = computeOutOfSampleResponses(nTest, nScenes, useMetaContrast, is2AFC, ...
            theNeuralEngine, fixationalEMObj, testNoiseFlag, ...
            theScenes, temporalSupport, ... 
            visualizeAllComponents, beVerbose)

    outSampleNoiseFreeStimResponsesCell = cell(1, nScenes);
    outSampleStimResponsesCell = cell(1,nScenes);
    neuralResponseEngineMetaData = cell(1, nScenes);

    if is2AFC
        nTest_eachScene = nTest;
    else
        nTest_eachScene = nTest/nScenes;
        assert(mod(nTest, nScenes) == 0, ['The number of test trials must be an',...
            ' integer multiple of the number of alternative choices']);
    end
    eStart = tic;

    % Get the video filename base and then append for each alternative.
    % Get basename for response visualization.
    %
    % Below we append a counter for each stimulus onto this
    if (iscell(theNeuralEngine))
        for nn = 1:length(theNeuralEngine)
            responseAltVideoFileNameBase{nn} = theNeuralEngine{nn}.responseVideoFileName;
        end
    else
        responseAltVideoFileNameBase = theNeuralEngine.responseVideoFileName;
    end


    for n = 1:nScenes
        % Adjust visualization video filename for alternative
        if (iscell(theNeuralEngine))
                theNeuralEngine{n}.responseVideoFileName = [responseAltVideoFileNameBase{n} sprintf('_scene%d',n)];
        else
            theNeuralEngine.responseVideoFileName = [responseAltVideoFileNameBase sprintf('_scene%d',n)];
        end

        if (useMetaContrast && ~is2AFC)
            % Noise free response for nth alternative scene sequence
            [outSampleNoiseFreeStimResponsesCell{n}, neuralResponseTemporalSupport] = theNeuralEngine{n}.computeNoiseFree(...
                theScenes, ...
                temporalSupport, ...
                'fixationalEM',fixationalEMObj);

            if (visualizeAllComponents)
                % Call the noise-free compute function directly to get the dataOut.metaData
                dataOut = theNeuralEngine{n}.noiseFreeComputeFunction(...
                    theNeuralEngine{n}, ...
                    theNeuralEngine{n}.noiseFreeComputeParams, ...
                    theScenes, ...
                    temporalSupport, ...
                    'fixationalEM',fixationalEMObj, ...
                    'visualizeActivationFunction', true);

                if (isfield(dataOut, 'metaData'))
                    dataOut.metaData.conditionLabel = conditionLabel;
                    neuralResponseEngineMetaData{n} = dataOut.metaData;
                end
            end

            % Add noise (or not) to nth alternative responses
            [outSampleStimResponsesCell{n}, ~] = theNeuralEngine{n}.computeNoisyInstances( ...
                outSampleNoiseFreeStimResponsesCell{n}, ...
                neuralResponseTemporalSupport, ...
                nTest_eachScene, ...
                testNoiseFlag);
        else
            % Noise free response for nth alternative scene sequence
            [outSampleNoiseFreeStimResponsesCell{n}, neuralResponseTemporalSupport] = theNeuralEngine.computeNoiseFree(...
                theScenes{n}, ...
                temporalSupport, ...
                'fixationalEM',fixationalEMObj);

            if (visualizeAllComponents)
                % Call the noise-free compute function directly to get the dataOut.metaData
                dataOut = theNeuralEngine.noiseFreeComputeFunction(...
                    theNeuralEngine, ...
                    theNeuralEngine.noiseFreeComputeParams, ...
                    theScenes{n}, ...
                    temporalSupport, ...
                    'fixationalEM',fixationalEMObj, ...
                    'visualizeActivationFunction', true);

                if (isfield(dataOut, 'metaData'))
                    dataOut.metaData.conditionLabel = conditionLabel;
                    neuralResponseEngineMetaData{n} = dataOut.metaData;
                end
            end

            % Add noise (or not) to nth alternative responses
            [outSampleStimResponsesCell{n}, ~] = theNeuralEngine.computeNoisyInstances( ...
                outSampleNoiseFreeStimResponsesCell{n}, ...
                neuralResponseTemporalSupport, ...
                nTest_eachScene, ...
                testNoiseFlag);
        end
    end

    % Restore responseVideoFileName to its pre-loop value so that a second
    % call to computePerformance (e.g. for the test phase after training)
    % does not accumulate an extra _scene<N> suffix.
    if (iscell(theNeuralEngine))
        for nn = 1:length(theNeuralEngine)
            theNeuralEngine{nn}.responseVideoFileName = responseAltVideoFileNameBase{nn};
        end
    else
        theNeuralEngine.responseVideoFileName = responseAltVideoFileNameBase;
    end

    e = toc(eStart);
    if (beVerbose)
        fprintf('computePerformance: Took %0.1f secs to generate test responses for all alternatives\n',e);
    end
end



function [inSampleStimResponsesCell, neuralResponseEngineMetaData, neuralResponseTemporalSupport, theNeuralEngine] = ...
    computeInSampleResponses(nTrain, nScenes, useMetaContrast, is2AFC, ...
    theNeuralEngine, fixationalEMObj, trainNoiseFlag, ...
    theScenes, temporalSupport, ...
    visualizeAllComponents)

    inSampleNoiseFreeStimResponsesCell = cell(1, nScenes);
    inSampleStimResponsesCell = cell(1,nScenes);
    neuralResponseEngineMetaData = cell(1, nScenes);
    
    for n = 1:nScenes
        if (useMetaContrast && ~is2AFC)
            % We will only call the nre visualization function of we are testing.
            % That simplifies the amount of stuff that gets written out.
            saveVisualizeEachCompute = theNeuralEngine{n}.visualizeEachCompute;
            theNeuralEngine{n}.visualizeEachCompute = false;
    
            % Noise free response for nth alternative scene sequence
            [inSampleNoiseFreeStimResponsesCell{n}, neuralResponseTemporalSupport] = theNeuralEngine{n}.computeNoiseFree(...
                theScenes, ...
                temporalSupport, ...
                'fixationalEM',fixationalEMObj);
    
            if (visualizeAllComponents)
                % Call the noise-free compute function directly to get the dataOut.metaData
                dataOut = theNeuralEngine{n}.noiseFreeComputeFunction(...
                    theNeuralEngine{n}, ...
                    theNeuralEngine{n}.noiseFreeComputeParams, ...
                    theScenes, ...
                    temporalSupport, ...
                    'fixationalEM',fixationalEMObj, ...
                    'visualizeActivationFunction', true);
    
                if (isfield(dataOut, 'metaData'))
                    dataOut.metaData.conditionLabel = conditionLabel;
                    neuralResponseEngineMetaData{n} = dataOut.metaData;
                end
            end
    
            % Add noise (or not) to nth alternative responses
            [inSampleStimResponsesCell{n}, ~] = theNeuralEngine{n}.computeNoisyInstances( ...
                inSampleNoiseFreeStimResponsesCell{n}, ...
                neuralResponseTemporalSupport, ...
                nTrain, ...
                trainNoiseFlag);
    
            % Put back visualization flag
            theNeuralEngine{n}.visualizeEachCompute = saveVisualizeEachCompute;
        else
            % We will only call the nre visualization function of we are testing.
            % That simplifies the amount of stuff that gets written out.
            saveVisualizeEachCompute = theNeuralEngine.visualizeEachCompute;
            theNeuralEngine.visualizeEachCompute = false;
    
            % Noise free response for nth alternative scene sequence
            [inSampleNoiseFreeStimResponsesCell{n}, neuralResponseTemporalSupport] = theNeuralEngine.computeNoiseFree(...
                theScenes{n}, ...
                temporalSupport, ...
                'fixationalEM',fixationalEMObj);
    
            if (visualizeAllComponents)
                % Call the noise-free compute function directly to get the dataOut.metaData
                dataOut = theNeuralEngine.noiseFreeComputeFunction(...
                    theNeuralEngine, ...
                    theNeuralEngine.noiseFreeComputeParams, ...
                    theScenes{n}, ...
                    temporalSupport, ...
                    'fixationalEM',fixationalEMObj, ...
                    'visualizeActivationFunction', true);
    
                if (isfield(dataOut, 'metaData'))
                    dataOut.metaData.conditionLabel = conditionLabel;
                    neuralResponseEngineMetaData{n} = dataOut.metaData;
                end
            end
    
            % Add noise (or not) to nth alternative responses
            [inSampleStimResponsesCell{n}, ~] = theNeuralEngine.computeNoisyInstances( ...
                inSampleNoiseFreeStimResponsesCell{n}, ...
                neuralResponseTemporalSupport, ...
                nTrain, ...
                trainNoiseFlag);
    
            % Put back visualization flag
            theNeuralEngine.visualizeEachCompute = saveVisualizeEachCompute;
        end
    end
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