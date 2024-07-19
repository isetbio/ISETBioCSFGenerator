function spdRadiance = generateAOPrimarySpd(wls,centerWl,FWHM,cornealPowerUW,stimulusSizeDeg,pupilDiameterMm)
% spdRadiance = generateAOPrimarySpd(wls,centerWl,FWHM,cornealPowerUW,stimulusSizeDeg)
%
% Description
%    Produce the spectral power distribution on wavelength support given by
%    wls, given primary center wavelength, FWHM, power at the cornea, and
%    stimulus size in degrees.  The returned spectral power distribution
%    has units of radiance, suitable for creating a display object.
%
%    Stimulus power is pecified here as the power entering the eye (going through
%    the pupil).  If the beam just fills the pupil, then this is just the power
%    measured at the cornea.
%
% If you have instead corneal irradiance in UW/mm^2 and this fills or overfills
% the pupil, then multiply by pupil area in mm^2 to get the power entering the eye.
%
%
%    All wls in nm.

% Create relative spd of the primary by modeling it as a normal 
relSpd = normpdf(wls,centerWl,FWHMToStd(FWHM));
unitSpd = relSpd/trapz(wls,relSpd);

% % Make spds that give desired full power.
cornealSpdUW = cornealPowerUW*unitSpd;

% % Get equivalent spectral radiance of background and test increment as a
% % function of the wavelength support.
% %
% % The routine here finds the radiance on an external conventional display that
% % produces the same retinal illuminance as the corneal power specified
% % above.  This is purely geometric calculation; attenuation of light by
% % occular media is not taken into account at this stage.  Note that this
% % conversion routine expects power per wavelength band, not power per nm,
% % as its input, but returns power per nm as its output.  A little
% % confusing, but that is why the spd being passed in is multiplied by the
% % wavelength spacing.
deltaWl = wls(2)-wls(1);
pupilAreaMm2 = pi*((pupilDiameterMm/2)^2);

% This is what we want to return as spd
spdRadiance = AOMonochromaticCornealPowerToRadiance(wls,wls,cornealSpdUW*deltaWl,pupilDiameterMm,stimulusSizeDeg^2);

% Make sure our computed radiance yields the desired corneal
% irradiance when we go in the other direction.
spdCornealIrradianceUWMm2Check = RadianceAndDegrees2ToCornIrradiance(spdRadiance,stimulusSizeDeg^2)*(1e6)*((1e-3)^2);
cornealPowerUWCheck = trapz(wls,spdCornealIrradianceUWMm2Check)*pupilAreaMm2;
if (abs(cornealPowerUWCheck-cornealPowerUW)/cornealPowerUW > 1e-4)
    error('Do not get right cornal power back from computed radiance');
end
