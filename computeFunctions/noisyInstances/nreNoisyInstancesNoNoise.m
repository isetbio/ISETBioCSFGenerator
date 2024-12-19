function dataOut = nreNoisyInstancesNoNoise(...
    neuralEngineOBJ, noisyInstancesComputeParams, noiseFreeResponses, ...
    temporalSupportSeconds, instancesNum, noiseFlag, varargin)
% Compute function for making instances without actually adding any noise
%
% Syntax:
%   dataOut = nreNoisyInstancesNoNoise(...
%    neuralEngineOBJ, noisyInstancesComputeParams, sceneSequence, ...
%    sceneSequenceTemporalSupport, instancesNum, varargin);
%
% Description:
%    Function serving as the computeFunctionHandle for a @neuralResponseEngine
%    object.  This version does not add noise, which might be useful for
%    testing or for creating multiple response instances in a manner
%    compatible with other nre noise adding functions.
%
%    There are 2 ways to use this function.
%
%       [1] If called directly and with no arguments, 
%           dataOut = nreNoisyInstancesNoNoise()
%       it does not compute anything and simply returns a struct with the 
%       default parameter structure for this computation.
%
%       [2] If called from a parent @neuralResponseEngine object, 
%       it computes 'instancesNum' of response instances but doesn't
%       add any noise.
%
%       It is not a terribly good idea to try to call this function with arguments
%       directly - it should be called by the compute method of its parent
%       @neuralResponseEngine.
%
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
%                                     special cases at that level. In this
%                                     particular routine, we can ignore the
%                                     noise flag since it does the same
%                                     thing in either case.
%
% Optional key/value input arguments:
%   'rngSeed'                       - Integer or string. Default is empty. This does not do anything
%                                     but is here for compatibility with other nre noise computation functions.
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
%                                   Here mResponses is the number of pixels times the
%                                   number of frames.
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
% See Also:
%     t_neuralResponseCompute

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
    p.KeepUnmatched = true;
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
    
    % Compute noisy response instances. Here we aren't adding noise,
    % so the noisy response is just the noise free response.
    [responseDim, framesNum] = size(noiseFreeResponses);
    noisyResponseInstances = zeros(instancesNum,responseDim,framesNum);
    switch (noiseFlag)
        case {'random', 'none'}
            for ii = 1:instancesNum
                noisyResponseInstances(ii,:,:) = noiseFreeResponses;
            end
        otherwise
            error('noiseFlag must be ''random'' or ''none''');
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
