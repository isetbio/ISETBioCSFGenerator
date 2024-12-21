function dataOut = sceMetaContrast(sceneEngineOBJ,testContrast,sceneParamsStruct)
% Compute function for generating a meta contrast scenes. These just return
% the numerical value of the contrast in a one pixel, one wavelength
% scene,one frame. The actual full temporal sequence is handled at the meta
% contrast nre stage.
%
% Syntax:
%   dataOut = sceMetaContrast(obj,testContrast,sceneParamsStruct);
%
% Description:
%    Compute function to be used as a computeFunctionHandle for a @sceneEngine
%    object. There are 2 ways to use this function.
%
%       [1] If called directly and with no arguments, 
%               dataOut = sceUniformFieldTemporalModulation();
%           it does not compute anything and simply returns a struct with the 
%           defaultParams that define the scene.
%
%       [2] If called with arguments, as it is from a parent @sceneEngine object,
%           it computes a cell array of scenes defining the frames of a
%           stimulus and the temporal support of the frames. These are
%           returned as named fields of the returned dataOut struct.
%
%    All scene functions used with the sceneEngine class must conform to
%    this API.
%
% Inputs:
%    sceneEngineOBJ             - Calling @sceneEngine object.  This is
%                                  currently unused, but passing it allows us
%                                  flexibility in the future and matches
%                                  conventions for the other classes in
%                                  this toolbox.
%    testContrast                - Scalar providing the contrast for the
%                                  scene to be generated.                           
%    sceneParamsStruct           - Struct containing properties of the
%                                  scene understood by this function.
%                                  As noted above, execute
%                                  sceMetaContrast at the
%                                  command line to see the structure's
%                                  fields and default values.
%
% Outputs:
%    dataOut  - A struct that depends on the input arguments. 
%
%               If called directly with no input arguments, the returned struct contains
%               the defaultParams that define the scene
%
%             - If called from a parent @sceneEngine, the returned
%               struct is organized as follows:
%                 .sceneSequence : a cell array of scenes defining the frames of the generated grating scene sequence                            
%                 .temporalSupport : the temporal support of the frames of the generated grating scene sequence, in seconds                   
%
% Optional key/value input arguments:
%    None.
%
% Examples:
%    The source code contains examples.
%
% See Also:
%     t_sceneGeneration, t_spatialCSF

% History:
%    08/09/2024  dhb  Wrote it.
%
%   Examples:
%{
    % Usage case #1. Just return the default scene params
    defaultParams = sceMetaContrast()

    % Usage case #2. Compute a stimulus sequence using a parent @sceneEngine 
    % object and the default scene params

    % Instantiate the parent @sceneEngine object
    theSceneEngineOBJ = sceneEngine(@sceMetaContrast);

    % Generate a test scene
    testContrast = 0.1;
    [theTestSceneSequence, temporalSupportSeconds] = ...
        theSceneEngineOBJ.compute(testContrast);  

    % Plot the encoded contrasts
    figure; clf; hold on;
    for ii = 1:length(temporalSupportSeconds)
        plot(temporalSupportSeconds(ii),sceneGet(theTestSceneSequence{ii},'photons'),'ro','MarkerFaceColor','r','MarkerSize',12);
        xlabel('Time (secs)');
        ylabel('Encoded Contrast')
    end
%}

    % Check input arguments. If called with zero input arguments, just return the default params struct
    if (nargin == 0)
        dataOut = generateDefaultParams();
        return;
    end
    
    % Create an equal photon uniform field
    uniformScene = sceneCreate('uniform monochromatic',550,1);

    % Stash contrast in only scene in scene sequence
    theSceneSequence{1} = sceneSet(uniformScene,'photons',testContrast); 

    % Generate temporal support for the scene sequence
    temporalSupportSeconds = 0; 
    
    % Assemble dataOut struct
    dataOut.sceneSequence = theSceneSequence;
    dataOut.temporalSupport = temporalSupportSeconds;
    dataOut.statusReport = struct();
end

function p = generateDefaultParams()
    p = struct( ...
        );
end
