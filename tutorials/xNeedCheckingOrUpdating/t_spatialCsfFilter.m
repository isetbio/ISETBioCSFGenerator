% Compute spatial CSF using different filters
%
% Description:
%    This script computes the contrast sensitivity function (CSF) under
%    different filtering conditions. It evaluates the effect of filters with
%    varying levels of light transmission on contrast sensitivity.
%
% See also: t_spatialCSF, t_thresholdEngine, t_modulatedGratingsSceneGeneration,
%           t_chromaticThresholdContour, computeThreshold, computePerformance
%

% History:
%   07/09/24  fh   Wrote it.

clear; close all;

% Compute the baseline contrast sensitivity function without any filter
threshold0 = t_spatialCSF('doValidationCheck', false); %no input filter
% Convert threshold to sensitivity (sensitivity is the inverse of threshold)
sensitivity0 = 1./threshold0;

%% filter 1
% Define the spectral range for the filter
wave = 400:20:740; % Wavelength range from 400 to 740 nm, stepping by 20 nm

% Create a perfect transmission filter (100% transmittance across all specified wavelengths)
% This could be a row or column vector. It doesn't matter.
transmission_func1 = ones(length(wave),1);

% Define filter properties in a struct
filter1 = struct('spectralSupport', wave, 'transmission', transmission_func1);

% Compute CSF threshold using the perfect transmission filter
threshold1 = t_spatialCSF('filter', filter1,'doValidationCheck', false);
sensitivity1 = 1./threshold1;

%% filter 2
% Create a neutral density filter with 25% transmittance
transmission_func2 = 0.25.*ones(length(wave),1);

% Define filter properties in a struct
filter2 = struct('spectralSupport', wave, 'transmission', transmission_func2);

% Compute CSF threshold using the neutral density filter
threshold2 = t_spatialCSF('filter', filter2,'doValidationCheck', false);
sensitivity2 = 1./threshold2;

%% Plot the results using logarithmic scales for both axes
% Expected Results:
% (1) The contrast sensitivity simulated using the perfect filter (100% transmittance) 
%     is expected to be nearly identical to that obtained without any filter. 
%     Note that there may be slight stochastic variability in the results.
% (2) The contrast sensitivity simulated using the neutral density filter (25% transmittance,
%     which attenuates the light level by a factor of 4) should be approximately half 
%     that of the sensitivity when no filter is applied.

figure; hold on;
% Plot identity line (baseline sensitivity)
h1 = loglog(sensitivity0, sensitivity0, 'k--', 'lineWidth', 2);
% Plot sensitivity with perfect filter
h2 = loglog(sensitivity0, sensitivity1, 'go-', 'lineWidth', 2); 
% Plot half the sensitivity of perfect filter (illustrative)
h3 = loglog(sensitivity0, sensitivity0./2, 'k:', 'lineWidth', 2);
% Plot sensitivity with 25% transmittance
h4 = loglog(sensitivity0, sensitivity2, 'bo-', 'lineWidth', 2); 
grid on; axis square; 
xticks([1,2,5,10,20,50]); yticks([1,2,5,10,20,50])
xlabel('Contrast sensitivity with no input filter');
ylabel('Contrast sensitivity with different filter types');
lgd2 = legend([h2, h4, h1, h3], {'Perfect filter (100% transmittance)',...
    'Neutral density filter (25% transmittance)','Identity line', ...
    'Half of the contrast sensitivity'}, 'Location', 'northwest');
title(lgd2, 'Filter type');
set(gca,'FontSize',12)

%%
% This block handles spectral mismatch between the scene's sampled wavelengths 
% and the filter by applying linear interpolation to adjust the filter's data.

% Define the sampled wavelengths for the scene, using a fine resolution
wave_scene  = 400:5:740; 

% Define the sampled wavelengths for the filter, using a coarser resolution
% than the scene for demonstration of interpolation necessity
wave_filter = 400:20:740;

% Create a Gaussian filter with a peak at 550 nm and a standard deviation of 50 nm
% The Gaussian filter values are normalized to have a maximum of 1
filter_Gaussian = normpdf(wave_filter, 550, 50); 
filter_Gaussian = filter_Gaussian./max(filter_Gaussian);

% Interpolate the Gaussian filter values to match the scene's finer wavelength sampling
% This ensures the filter's transmittance values are defined at every wavelength used by the scene
transmission_wvl = interp1(wave_filter, filter_Gaussian, wave_scene, 'linear');

% Visualizing the original and interpolated filter functions
figure; hold on
plot(wave_filter,filter_Gaussian, 'r-*', 'LineWidth',2); 
plot(wave_scene, transmission_wvl, 'ko:');
xlabel('Sampled wavelength (nm)');
ylabel('Transmittance');
legend({'Original filter function', 'Interpolated filter function'})
set(gca,'FontSize', 12);

