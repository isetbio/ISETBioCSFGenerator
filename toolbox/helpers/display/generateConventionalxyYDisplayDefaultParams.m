function displayParams = generateConventionalxyYDisplayDefaultParams
% Setup default parameters for a conventional display to be used in a scene
%
% Syntax:
%    displayParams = generateConventionalxyYDisplayDefaultParams
%
% Description:
%    This generates a reasonable parameter set for a conventional monitor
%    presentation display. These are accepted by function
%    generateConventionalxyYDisplay.
%
% Inputs:
%    None.
%
% Outputs:
%    displayParams               - A display parameter structure.  See code
%                                  for list of fields and what they are.
% Optional key/value pairs:
%     None.
%
% See also: generateConventionalxyYDisplay, sceGrating,
%           createGratingSceneEngine

% History:
%   01/12/2024  dhb  Refactored this into its own function.

displayParams = struct( ...
    'displayType', 'conventionalxyY', ...           % display type
    'whichDisplay', 'LCD-Apple',...                 % display: the monitor type
    'meanLuminanceCdPerM2', 40, ...                 % display: desired mean luminance. Determines display scale
    'luminanceHeadroom', 0.05, ...                  % display: how much headroom to leave re 2*meanLuminance in display
    'viewingDistanceMeters', 0.57, ...              % display: viewing distance
    'bitDepth', 20, ...                             % display: length of LUT
    'gammaTableExponent', 2.0, ...                  % display: shape of LUT, 1 = Linear
    'spectralSupport', [] ...                       % display: spectral support. Will get filled in from the parent params
    );
end