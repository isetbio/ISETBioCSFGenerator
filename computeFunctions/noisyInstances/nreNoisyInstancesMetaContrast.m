function dataOut = nreNoisyInstancesMetaContrast(...
    neuralEngine, metaNoisyInstancesParams, noiseFreeResponses, ...
    temporalSupportSeconds, instancesNum, noiseFlag, varargin)
% Meta contrast nre noisy instances compute function
%
% Syntax:
%   dataOut = nreNoisyInstancesMetaContrast(...
%      neuralEngine, metaNoisyInstancesParams, noiseFreeResponses, ...
%      temporalSupportSeconds, instancesNum, noiseFlag)
%
% Description:
%    Function serving as the computeFunctionHandle for a @neuralResponseEngine
%    object.  This version adds noise using the compute method of the
%    underlying nre.
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
% Inputs:
%    neuralEngine                   - the parent @neuralResponseEngine object that
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
%   'rngSeed'                       - Integer or string.  Set rng seed.
%                                     Empty (default) means don't touch the
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
%     nreNoiseFreeMetaContrast,  t_metaContrastCSF

% History:
%    08/11/24    dhb  Wrote it.
%    12/20/2024  dhb  Create for major architecture redo.

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
if (isempty(neuralEngine.neuralPipeline) | ~isfield(neuralEngine.neuralPipeline,'noisyInstances'))
    noisyInstancesPipeline.dummyfield = 1;
    returnTheNoisyInstancesPipeline = true;
else
    returnTheNoisyInstancesPipeline =  false;
end

% Set rng seed
if (~isempty(p.Results.rngSeed))
    oldSeed = rng(p.Results.rngSeed);
end

% Compute noisy response instances
[noiseFreeInstancesNum, responseDim, framesNum] = size(noiseFreeResponses);
if (noiseFreeInstancesNum ~= 1)
    error('Need to generalize beyond one noise-free instance');
end
noisyResponseInstances = zeros(instancesNum,responseDim,framesNum);
switch (noiseFlag)
    case {'random'}
        switch (metaNoisyInstancesParams.noiseAddMethod)
            case 'nre'
                noisyResponseInstances = metaNoisyInstancesParams.neuralEngine.computeNoisyInstances(...
                    noiseFreeResponses, ...
                    temporalSupportSeconds, ...
                    instancesNum, ...
                    'random');
            otherwise
                error('Unknown noise add method specified');
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

% Assemble the dataOut struct
dataOut = struct(...
    'neuralResponses', noisyResponseInstances, ...
    'temporalSupport', temporalSupportSeconds);

if (returnTheNoisyInstancesPipeline)
    dataOut.noisyInstancesPipeline = noisyInstancesPipeline;
end

end

% Default params for this compute function
function p = generateDefaultParams()
p = struct( ...
    'noiseAddMethod', 'nre' ...
    );
end
