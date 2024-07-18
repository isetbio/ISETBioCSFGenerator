function spd = generateAOPrimarySpd(wls,centerWl,FWHM,uWCornea,stimulusSizeDeg,pupilDiameterMm)
% spd = generateAOPrimarySpd(wls,centerWl,FWHM,uWCornea,stimulusSizeDeg)
%
% Description
%    Produce the spectral power distribution on wavelength support given by
%    wls, given primary center wavelength, FWHM, power at the cornea, and
%    stimulus size in degrees.
%
%    All wls in nm.

% Create relative spd of the primary by modeling it as a normal 
relSpd = normpdf(wls,centerWl,FWHMToStd(FWHM));
unitSpd = relSpd/trapz(twoSpotParams.wls,relSpd);

% % Stimulus power
% %
% % Specified here as the power entering the eye (going through
% % the pupil).  If the beam just fills the pupil, then this is just the power
% % measured at the cornea.
% %
% % If you have instead corneal irradiance in UW/mm^2 and this fills or overfills
% % the pupil, then multiply by pupil area in mm^2 to get the power entering the eye.
% %
% % This code assumes that the power for both background and test is spread
% % out over the background area. This is because typically for an AOSLO we
% % measure power for the whole raster, and then in an experiment switch it
% % off for part of the time to make the spatial pattern we want.
% %
% % The test should be specified as the amount of power added to
% % the background (its incremental power).  This is what you'll get if
% % background and test are in separate channels of the system and you
% % measure them separately.
% imagingBgCornealPowerUW = twoSpotParams.imagingBgPowerUW;
% 
% % Make spds that give desired full power.
% backgroundSpdCornealPowerUW = imagingBgCornealPowerUW*unitSpd;
% 
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
% deltaWl = twoSpotParams.wls(2)-twoSpotParams.wls(1);
% backgroundSizeDegs = twoSpotParams.spotBgDegs;
% pupilDiameterMm = twoSpotParams.pupilDiameterMm;
% pupilAreaMm2 = pi*((pupilDiameterMm/2)^2);

% This is what we want to return as spd
% imagingBgSpdRadiance = AOMonochromaticCornealPowerToRadiance(twoSpotParams.wls,twoSpotParams.wls,backgroundSpdCornealPowerUW*deltaWl,pupilDiameterMm,backgroundSizeDegs^2);
% 
% % Make sure our computed radiance yields the desired corneal
% % irradiance when we go in the other direction.
% imagingBgSpdCornealIrradianceUWMm2Check = RadianceAndDegrees2ToCornIrradiance(imagingBgSpdRadiance,backgroundSizeDegs^2)*(1e6)*((1e-3)^2);
% imagingCornealPowerUWCheck = trapz(twoSpotParams.wls,imagingBgSpdCornealIrradianceUWMm2Check)*pupilAreaMm2;
% if (abs(imagingCornealPowerUWCheck-imagingBgCornealPowerUW)/imagingBgCornealPowerUW > 1e-4)
%     error('Do not get right cornal power back from computed radiance');
% end
