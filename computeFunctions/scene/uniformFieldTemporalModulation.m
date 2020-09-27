function dataOut = uniformFieldTemporalModulation(testContrast,sceneParamsStruct)
% Compute function for computation of cone excitations witout eye movements 
%
% Syntax:
%   dataOut = uniformFieldTemporalModulation(testContrast,sceneParamsStruct);
%
% Description:
%    Compute function to be used as a computeFunctionHandle for a @sceneEngine
%    object. There are 2 ways to use this function.
%
%       [1] If called directly and with no arguments, 
%
%               dataOut = uniformFieldTemporalModulation()
%
%           it does not compute anything and simply returns a struct with the 
%           defaultParams that define the scene.
%
%       [2] If called from a parent @sceneEngine object, it computes a cell array 
%           of scenes defining the frames of a stimulus and the temporal
%           support of the frames
%
% Inputs:
%    testContrast                - the contrast for the scene to be generated
%                               
%    sceneParamsStruct            - a struct containing properties of the scene
%
%
% Optional key/value input arguments: none 
%
% Outputs:
%    dataOut  - A struct that depends on the input arguments. 
%
%               If called directly with no input arguments, the returned struct contains
%               the defaultParams that define the scene
%
%             - If called from a parent @sceneEngine, the returned
%               struct is organized as follows:
%
%              .sceneSequence : a cell array of scenes defining the frames of the generated scene sequence
%                               
%              .temporalSupport : the temporal support of the frames of the generated scene sequence, in seconds            
%
%
%
% See Also:
%     t_sceneGeneration

% History:
%    09/26/2020  NPC  Wrote it.
%
%   Examples:
%{
    % Usage case #1. Just return the default scene params
    defaultParams = uniformFieldTemporalModulation()

    % Usage case #2. Compute a stimulus sequence using a parent @sceneEngine 
    % object and the default neural response params

    % Instantiate the parent @sceneEngine object
    theSceneEngineOBJ = sceneEngine(@uniformFieldTemporalModulation);

    % Generate a test scene
    testContrast = 0.1;
    [theTestSceneSequence, theTestSceneTemporalSupportSeconds] = ...
        theSceneEngineOBJ.compute(testContrast);
    
%}
`
    % Check input arguments. If called with zero input arguments, just return the default params struct
    if (nargin == 0)
        dataOut = generateDefaultParams();
        return;
    end
    
    % Create an equal photon uniform field
    uniformScene = sceneCreate('uniform equal photon', sceneParamsStruct.sizePixels);

    % Set the scene width in degrees
    uniformScene = sceneSet(uniformScene, 'wAngular', sceneParamsStruct.fovDegs);

    % 1 meter away
    uniformScene = sceneSet(uniformScene, 'distance', 1.0);

    % background: adjust radiance according to desired  mean luminance
    meanLuminance = sceneParamsStruct.meanLuminanceCdPerM2;
    backgroundScene = sceneAdjustLuminance(uniformScene, meanLuminance);

    % pedestal: adjust radiance according to desired  mean luminance and test contrast
    testLuminance = meanLuminance * (1.0 + testContrast);
    pedestalScene = sceneAdjustLuminance(uniformScene, testLuminance);
    
    % Generate temporal support for the scene sequence
    temporalSupportSeconds = (0:(sceneParamsStruct.stimDurationFramesNum-1))*(sceneParamsStruct.frameDurationSeconds);
    
    % Generate the scene sequence
    theSceneSequence = cell(1, sceneParamsStruct.stimDurationFramesNum);
    for frameIndex = 1:sceneParamsStruct.stimDurationFramesNum
        if (ismember(frameIndex, sceneParamsStruct.stimOnsetFramesIndices))
            theSceneSequence{frameIndex} = pedestalScene;
        else
            theSceneSequence{frameIndex} = backgroundScene;
        end
    end
    
    % Assemble dataOut struct
    dataOut.sceneSequence = theSceneSequence;
    dataOut.temporalSupport = temporalSupportSeconds;
end

function p = generateDefaultParams()
    p = struct(...
        'fovDegs', 0.25, ...                        % 0.25 degs across
        'meanLuminanceCdPerM2', 100, ...            % 100 cd/m2 mean luminance
        'frameDurationSeconds', 50/1000, ...        % 50 msec frame duration
        'stimDurationFramesNum', 4, ...             % total time: 200 msec
        'stimOnsetFramesIndices', [2 3], ...        % modulate luminance at frames 1 and 2, so between 50 and 150 msec
        'sizePixels', 64 ...                        % 64 x 64 pixels
    );
end
