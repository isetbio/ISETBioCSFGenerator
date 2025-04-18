function dataOut = nreNoisyInstancesGaussian(...
    neuralEngineOBJ, noisyInstancesComputeParams, noiseFreeResponses, ...
    temporalSupportSeconds, instancesNum, noiseFlag, varargin)
% Compute function for adding Poisson noise to neural responses
%
% Syntax:
%   dataOut = nreNoisyInstancesGaussian(...
%    neuralEngineOBJ, noisyInstancesComputeParams, sceneSequence, ...
%    sceneSequenceTemporalSupport, instancesNum, varargin);
%
% Description:
%    Function serving as the computeFunctionHandle for a @neuralResponseEngine
%    object.  This version adds Gaussian noise.
%
%    Currently, iid Gaussian noise is implemented, with specified standard
%    deviation.
%
%    There are 2 ways to use this function.
%
%       [1] If called directly and with no arguments, 
%           dataOut = nreNoisyInstancesPoisson()
%       it does not compute anything and simply returns a struct with the 
%       default parameter structure for this computation.
%
%       [2] If called from a parent @neuralResponseEngine object, 
%       it computes 'instancesNum' of noisy response instances using
%       the Gaussian noise model, based on the passed noise free responses.
%
%       It is not a terribly good idea to try to call this function with arguments
%       directly - it should be called by the compute method of its parent
%       @neuralResponseEngine.
%
% In addition to computing, this function checks the `visualizeEachCompute` 
%    flag of the neuralEngineOBJ and, if it set, calls the appropriate
%    visualization function. This causes figures to appear that visualize
%    the noisy spatiotemporal response instances, which is helpful for debugging.
%    Note that everything runs much more slowly in this case.
%
% Inputs:
%    neuralEngineOBJ   - the parent @neuralResponseEngine object that
%                                     is calling this function as its computeFunctionHandle
%    noisyInstanceComputeParams - a struct containing properties of the employed neural chain.
%                                     This should just be empty for this function.                                
%    noiseFreeResponses  - a matrix of noise free neural
%                                     responses, with each column providing the responses for one frame.
%    temporalSupportSeconds  - the temporal support for the stimulus frames, in seconds
%    instancesNum        - the number of response instances to compute
%    noiseFlag               - string: 'random' (default) or 'none'.
%                                     The defalt value of 'random' means
%                                     actually add the noise.  The value of
%                                     'none' means return instancesNum
%                                     instances of noise versions of the
%                                     passed responses. This option may
%                                     seem silly, but it simplifies some
%                                     cases of calling code by avoiding
%                                     special cases at that level.
%
% Outputs:
%    dataOut  - A struct that depends on the input arguments. 
%
%               If called directly with no input arguments, the returned struct contains
%               the defaultParams, which here is the empty struct.
%
%             - If called from a parent @neuralResponseEngine, the returned
%               struct is organized as follows:
%
%              .noisyResponseInstances : matrix of noisy response
%                                   instances. Indexing is
%                                   (instance,:,frame) where the colon
%                                   indicates the vector of responses for
%                                   that instance/frame pair.  So the
%                                   dimension of the 3D matrix is
%                                   [instancesNum x mResponses x tTimeBins].
%              .temporalSupport : the temporal support of the neural
%                                   responses, in seconds.  This is just
%                                   returned as what was passed.
%              .noisyInstancesPipeline : a struct that the parent
%                                   @neuralResponseEngine
%                                   can tuck away.  Only returned if the
%                                   parent object has an empty parameters
%                                   struct. For this routine, it is just a
%                                   dummy struct as it isn't needed by this
%                                   routine.
%
%       The computed noisy response instances are a 3D matrix that has the form:
%           [instancesNum x mResponses x nFrames] 
%       where the mResponses is the dimension of the passed response vector for each
%       frame.
% 
% Optional key/value input arguments:
%   'rngSeed'                 - Integer or string.  Set rng seed. Empty (default) means don't touch the
%                                     seed. An integer uses that as the
%                                     seed.  A string (e.g. 'shuffle') sets
%                                     the seed by passing the string into
%                                     the rng() function. When the rng seed
%                                     is set, the old seed is saved and
%                                     restored upon return. Note that
%                                     setting the seed is slow, so you want
%                                     to invoke this option with some
%                                     trepidation.
%
% See Also:
%     nreNoiseFreeSceneAsResponse, t_neuralResponseCompute

% History:
%    09/26/2020  dhb  Wrote it.
%    10/19/2020  dhb  Fix comment to reflect fact that we now return
%                     instancesNum instances in noise free case.

    % Check input arguments. If called with zero input arguments, just return the default params struct
    if (nargin == 0)
        dataOut = generateDefaultParams();
        return;
    end
    
    % Parse the input arguments
    p = inputParser;
    p.addParameter('rngSeed',[],@(x) (isempty(x) | isnumeric(x) | ischar(x)));
    varargin = ieParamFormat(varargin);
    p.parse(varargin{:});

    % Return a neural pipeline struct to keep parent object happy if
    % it hasn't yet stored one.
    if (isempty(neuralEngineOBJ.neuralPipeline) | ~isfield(neuralEngineOBJ.neuralPipeline,'noisyInstances'))
        % At present, sigma must be a scalar.  Maybe in the future it could be
        % a covariance matrix.
        sigma = noisyInstancesComputeParams.sigma;
        if (~isscalar(sigma))
            error('The standard deviation parameter sigma must be a scalar');
        end

        returnTheNoisyInstancesPipeline = true;
    else
        sigma = neuralEngineOBJ.neuralPipeline.noisyInstances.sigma;

        returnTheNoisyInstancesPipeline =  false;
    end

    % Set rng seed
    if (~isempty(p.Results.rngSeed))
        oldSeed = rng(p.Results.rngSeed);
    end

    % Check whether the noise-free responses are cone contrast based, for
    % visualizaition. This is used so that the visualization routine
    % can more sensibly display things.  Each nre has this information in
    % different fields or not at all, so this is a little fussy. Next time
    % this breaks we should encapsulate this in one place so that can be
    % firxed once rather than fixing in every nreNoisyInstances routine.
    if (isfield(neuralEngineOBJ.noiseFreeComputeParams,'coneMosaicParams'))
        noiseFreeResponsesAreContrastResponses = strcmp(neuralEngineOBJ.noiseFreeComputeParams.coneMosaicParams.outputSignalType,'coneContrast');
    elseif (isfield(neuralEngineOBJ.noiseFreeComputeParams,'mRGCMosaicParams'))
        noiseFreeResponsesAreContrastResponses = strcmp(neuralEngineOBJ.noiseFreeComputeParams.mRGCMosaicParams.inputSignalType,'coneContrast');
    else
        noiseFreeResponsesAreContrastResponses = false;
    end


    % Get dimensions of noiseFreeResponses
    [noiseFreeInstancesNum, responseDim, framesNum] = size(noiseFreeResponses);
    if (noiseFreeInstancesNum ~= 1)
        error('Need to generalize beyond one noise-free instance');
    end

    noisyResponseInstances = zeros(instancesNum,responseDim,framesNum);

    switch (noiseFlag)
        case {'random'}
            for ii = 1:instancesNum
                % Get and add Gaussian noise
                noise = sigma .* randn(size(noiseFreeResponses));
                noisyResponseInstances(ii,:,:) = noiseFreeResponses + noise;
            end
        case {'none'}
            for ii = 1:instancesNum
                noisyResponseInstances(ii,:,:) = noiseFreeResponses;
            end
        otherwise
            error('noiseFlag must be ''random'' or ''none''');
    end

    % If an activation function has been specified, apply it to the noisy repsonse instances here
    if (isfield(noisyInstancesComputeParams, 'activationFunctionParams'))
       noisyResponseInstances = neuralResponseEngine.applyActivationFunction(...
            noiseFreeResponses, noisyResponseInstances, noisyInstancesComputeParams.activationFunctionParams);
    end

    % Restore
    if (~isempty(p.Results.rngSeed))
        rng(oldSeed);
    end
    
    % Check the visualizeEachCompute flag of the neuralEngineOBJ , and if set to true,
    % call the appropriate visualization function to visualize the generated 
    % spatiotemporal noisy response instances.
    if (neuralEngineOBJ.visualizeEachCompute)
        % Visualize computed data
        hFig = figure(1001);
        set(hFig, 'Position', [350 25 1650 550]);
        neuralEngineOBJ.visualize(noisyResponseInstances, temporalSupportSeconds, ...
            'figureHandle', hFig, ...
            'responseLabel', 'noisy response instances (Gaussian)', ...
            'responseVideoFileName', neuralEngineOBJ.responseVideoFileName, ...
            'neuralPipelineID', neuralEngineOBJ.ID, ...
            'visualizeResponsesAsModulations', noiseFreeResponsesAreContrastResponses);
    end

    % Assemble the dataOut struct
    dataOut = struct(...
        'neuralResponses', noisyResponseInstances, ...
    	'temporalSupport', temporalSupportSeconds);

    if (returnTheNoisyInstancesPipeline)
        noisyInstancesPipeline.sigma = sigma;
        dataOut.noisyInstancesPipeline = noisyInstancesPipeline;
    end
end

function p = generateDefaultParams()
    % Default params for this compute function
    p = struct( ...
        'sigma',1 ...
        );
end
