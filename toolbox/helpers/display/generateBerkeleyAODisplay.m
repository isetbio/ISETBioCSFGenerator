function presentationDisplay = generateBerkeleyAODisplay(displayParams)
% Setup a display that mimics the Berkeley AO system
%
% Syntax:
%    presentationDisplay = generateBerkeleyAODisplay(displayParams)
%
% Description:
%    This generates a display onto which we can describe stimuli. This
%    produces a display that can mimic the AOSLO system at Berkeley at
%    the frame rate level of temporal precision.
%
%    Probably this could be fairly easily adapted to other AO display
%    systems.
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
% See also: generateConventionalxyYDisplayDefaultParams,
%           createGratingSceneEngine

% History:
%   01/12/2024  dhb  Refactored this into its own function.

% Check display parameter type
if (~strcmp(displayParams.displayType,'berkeleyAO'))
    error('Incorrect display parameters type passed');
end

% We know the dispaly FOV and number of pixels.  We choose a distance
% far enough way and compute DPI to make it work out.
inchesPerMeter = 39.3701;
displaySizeMeters = 2*displayParams.viewingDistanceMeters*tand(displayParams.displayFOVDeg/2);
displayDotsPerMeter = displayParams.displayPixelSize/displaySizeMeters;
displayDPI = displayDotsPerMeter/inchesPerMeter;

% Now build the display given all the wonderful things we know about it.
presentationDisplay = generateCustomDisplay(...
    'viewingDistanceMeters', displayParams.viewingDistanceMeters, ...
    'dotsPerInch', displayDPI, ...
    'wavelengthSupportNanoMeters', displayParams.spectralSupport, ...
    'spectralPowerDistributionWattsPerSteradianM2NanoMeter', displayParams.spd, ...
    'ambientSPDWattsPerSteradianM2NanoMeter', displayParams.ambientSpd, ...
    'gammaTable', repmat((linspace(0,1,2^displayParams.bitDepth)').^displayParams.gammaTableExponent, [1 3]), ...
    'plotCharacteristics', displayParams.plotDisplayCharacteristics);
end