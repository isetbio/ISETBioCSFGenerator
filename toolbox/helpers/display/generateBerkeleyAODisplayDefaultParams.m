function displayParams = generateBerkeleyAODisplayDefaultParams
% Setup default parameters for a Berkeley AO display to be used in a scene
%
% Syntax:
%    displayParams = generateBerkeleyAODisplayDefaultParams
%
% Description:
%    This generates a reasonable parameter set for a model of a Berkeley AO
%    presentation display. These are accepted by function
%    generateBerkeleyAODisplay.
%
% Inputs:
%    None.
%
% Outputs:
%    displayParams               - A display parameter structure.  See code
%                                  for list of fields and what they are.
%
% Optional key/value pairs:
%     None.
%
% See also: generateBerkelyAODisplay, generateConventionalxyYDisplay

% History:
%   01/12/2024  dhb  Refactored this into its own function.

displayParams = struct( ...
    'displayType', 'berkeleyAO', ...        % Display type
    'displayPixelSize', 512, ...            % Linear display pixel size
    'displayFOVDeg', 1.413, ...             % Linear field size in degrees
    'AOPrimaryWls', [840 650 540], ...      % Display spd center wavelengths
    'AOPrimaryFWHM', [10 10 10], ...        % Display spd FWHM in nm
    'AOAOCornealPowersUW', [141.4 0 0], ... % Display spd power full on
    'viewingDistanceMeters', 3, ...         % Far enough away to be in good focus
    'spectralSupport', [], ...              % Spectral support. Will get filled in from the parent params
    'bitDepth', 10, ...                     % Bit depth length of LUT
    'gammaTableExponent', 1.0 ...           % Shape of LUT as gamma exponent, 1 = Linear
    );
end