function dataOut = nreNoisyInstancesPoisson(...
    neuralEngineOBJ, noisyInstancesComputeParams, noiseFreeResponses, ...
    temporalSupportSeconds, instancesNum, noiseFlag, varargin)
% Compute function for adding Poisson noise to neural responses
%
% Syntax:
%   dataOut = nreNoisyInstancesPoisson(...
%    neuralEngineOBJ, noisyInstancesComputeParams, sceneSequence, ...
%    sceneSequenceTemporalSupport, instancesNum, varargin);
%
% Description:
%    Function serving as the computeFunctionHandle for a @neuralResponseEngine
%    object.  This version adds Poisson noise
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
%       a Poisson noise model, based on the passed noise free responses.
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


% Inputs:
%    neuralEngineOBJ                - the parent @neuralResponseEngine object that
%                                     is calling this function as its computeFunctionHandle
%    noisyInstanceComputeParams     - a struct containing properties of the employed neural chain.
%                                     This should just be empty for this function.                                
%    noiseFreeResponses             - a matrix of noise free neural
%                                     responses, with each column providing the responses for one frame.
%    temporalSupportSeconds         - the temporal support for the stimulus frames, in seconds
%    instancesNum                   - the number of response instances to compute
%    noiseFlag                      - string: 'random' (default) or 'none'.
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
%   'rngSeed'                       - Integer or string.  Set rng seed. Empty (default) means don't touch the
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
    
    % Return a dummy neural pipeline struct to keep parent object happy if
    % it hasn't yet stored one.
    if (isempty(neuralEngineOBJ.neuralPipeline) | ~isfield(neuralEngineOBJ.neuralPipeline,'noisyInstances'))
        noisyInstancesPipeline.dummyfield = 1;
        returnTheNoisyInstancesPipeline = true;
    else
        returnTheNoisyInstancesPipeline =  false;
    end

    % Set rng seed
    if (~isempty(p.Results.rngSeed))
        oldSeed = rng(p.Results.rngSeed);
    end

    % Check whether the noise-free response is a cone contrast based
    % response and raise an error
    noiseFreeResponsesAreContrastResponses = strcmp(neuralEngineOBJ.noiseFreeComputeParams.coneMosaicParams.outputSignalType,'coneContrast');
    if (noiseFreeResponsesAreContrastResponses)
        error('You cannot have Poisson noisy instances when the mean response is a contrast-based response\n');
    end

    % Compute noisy response instances
    [noiseFreeInstancesNum, responseDim, framesNum] = size(noiseFreeResponses);
    if (noiseFreeInstancesNum ~= 1)
        error('Need to generalize beyond one noise-free instance');
    end

    noisyResponseInstances = zeros(instancesNum,responseDim,framesNum);
    switch (noiseFlag)
        case {'random'}
            % Find index to low mean noise free responses
            poissonCriterion = 25;
            idx = find(noiseFreeResponses(:) < poissonCriterion);
            
            for ii = 1:instancesNum
                % You might think that you could just use Matlab's
                % poissnrnd as in this commented out code.  But it is very
                % slow. So instead we use the Gaussian appoximation once
                % the mean gets reasonably large.  The Gaussian
                % approximation still keeps variance proportional to mean,
                % and we round so will still always get an integer.  It is
                % considerably faster.
                %
                % noisyResponseInstances(ii,:,:) = poissrnd(noiseFreeResponses);
                
                % Gaussian approximation for large mean; for small counts (less
                % than poissonCriterion) create a true Poisson sample.
                noise = sqrt(noiseFreeResponses) .* randn(size(noiseFreeResponses));
                noisyResponseInstances(ii,:,:) = round(noiseFreeResponses + noise);
                if ~isempty(idx)
                    noisyResponseInstances(ii,idx) = poissrnd(noiseFreeResponses(idx));
                end
            end
        case {'none'}
            for ii = 1:instancesNum
                noisyResponseInstances(ii,:,:) = noiseFreeResponses;
            end
        otherwise
            error('noiseFlag must be ''random'' or ''none''');
    end

    % Restore
    if (~isempty(p.Results.rngSeed))
        rng(oldSeed);
    end
    
    % Check the visualizeEachCompute flag of the neuralEngineOBJ, and if set to true,
    % visualize the generated noisy response instances.
    if (neuralEngineOBJ.visualizeEachCompute)
        % Visualize computed data
        hFig = figure(1001);
        set(hFig, 'Position', [350 25 1650 550]);
        neuralEngineOBJ.visualize(noisyResponseInstances, temporalSupportSeconds, ...
            'figureHandle', hFig, ...
            'responseLabel', 'Poisson response instances', ...
            'responseVideoFileName', neuralEngineOBJ.responseVideoFileName, ...
            'neuralPipelineID', neuralEngineOBJ.ID, ...
            'visualizeResponsesAsModulations', noiseFreeResponsesAreContrastResponses);

    end

    % Assemble the dataOut struct
    dataOut = struct(...
        'neuralResponses', noisyResponseInstances, ...
    	'temporalSupport', temporalSupportSeconds);

    if (returnTheNoisyInstancesPipeline)
        dataOut.noisyInstancesPipeline = noisyInstancesPipeline;
    end
end

function p = generateDefaultParams()
    % Default params for this compute function
    p = struct([]);
end
