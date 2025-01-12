function presentationDisplay = generateConventionalxyYDisplay(displayParams)
% Setup conventional display to be used for scene generation
%
% Syntax:
%    presentationDisplay = generateConventionalxyYDisplay(displayParams)
%
% Description:
%    This generates a display onto which we can describe stimuli. This
%    produces a conventional display with max luminance adjusted so that
%    the passed mean luminance is near the middle of the lumiannce range.
%
% Inputs:
%    displayParams              - Parameter struct describing the display
%                                 properties. See
%                                 generateConventionalxyYDisplayDefaultParams.
%
% Outputs:
%    presentationDisplay        - A ISETBio presentation display object.

% Optional key/value pairs:
%     None.
%
% See also: generateConventionalxyYDisplayDefaultParams, sceGrating,
%           createGratingSceneEngine

% History:
%   01/12/2024  dhb  Refactored this into its own function.

% Check display parameter type
if (~strcmp(displayParams.displayType,'conventionalxyY'))
    error('Incorrect display parameters type passed');
end

% Generic LCD display
presentationDisplay = displayCreate(displayParams.whichDisplay, ...
    'viewing distance', displayParams.viewingDistanceMeters, ...
    'wave', displayParams.spectralSupport ...        % Custom spectral support
    );

% Custom LUT and bit depth
if (displayParams.gammaTableExponent == 1)
    % Linear LUT
    N = 2^displayParams.bitDepth;
    gTable = repmat(linspace(0, 1, N), 3, 1)';
else
    x = linspace(0,1,2^displayParams.bitDepth);
    gTable = x(:).^displayParams.gammaTableExponent;
    gTable = repmat(gTable, [1,3]);
end

% Set the gamma table
presentationDisplay = displaySet(presentationDisplay, 'gTable', gTable);

% Adjust display SPDs so we can generate desired mean luminance,
% with a little headroom as specified.
desiredPeakLuminance = displayParams.meanLuminanceCdPerM2 * 2*(1+displayParams.luminanceHeadroom);
peakLuminance = displayGet(presentationDisplay, 'peak luminance');
scalingFactor = desiredPeakLuminance / peakLuminance;
spds = displayGet(presentationDisplay, 'spd primaries');
presentationDisplay = displaySet(presentationDisplay, 'spd primaries', spds*scalingFactor);

end