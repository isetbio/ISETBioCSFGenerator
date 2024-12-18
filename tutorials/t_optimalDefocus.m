function [defocus_adj_diopters, defocus_adj_microns] = t_optimalDefocus(...
    Subj, varargin)
%{
The goal of this script is to find the amount of defocus that can lead
to the most compact PSF at the FOVEA. Specifically,
1. We pick a in-focus wavelength (e.g., 550 nm) and set it to be the 
    calculate wavelength 
2. Vary the amount of defocus and convert it to microns
3. Loop through all the defocus
    Add it to the 5th Zernike polynomial
    Compute the PSF
    Find the peak of the PSF / compute the Strehl ratio
4. Find the defocus that corresponds to the highest peak of the PSF / the 
    largest Strehl ratio
5. Validation
    Add the optimal defocus to the 5th Zernike polynomial 
    Repeat everything above 
    The PSF should be the most compact when added defocus = 0

Output: this script will generate four figures
Fig. 1: A heatmap showing the peak of the PSF's with varying in-focus wvl
and added defocus; a slice of the heatmap for in-focus wvl = 550 nm. This
plot also marks at what amount of defocus the peak of the PSF is the
highest or the Strehl ratio is the greatest.

Fig. 2: the point spread function with 0 added defocus at three selected
in-focus wavelength 

Fig. 3-4: same as Figs.1-2 except that we added the optimal amount of defocus so
that the PSF can be the most compact given the in-focus wavelength = 550nm
%}

%% Parse input
p = inputParser;
p.addParameter('measuredPupilMM', 4,@isscalar);  %measured pupil size (mm)
p.addParameter('calcPupilMM', 3, @isscalar);     %pupil size for calculation (mm)
p.addParameter('measuredWvl', 550, @(x)(isscalar(x) && floor(x)==ceil(x) && x <= 800 && x >= 400)); %measured wavelength (nm)
p.addParameter('infocusWvl', 550, @(x)(isscalar(x) && floor(x)==ceil(x) && x <= 650 && x >= 450));  %select one in-focus wavelength for visualization
p.addParameter('optics_dataset', 'Polans2015', @ischar);
p.addParameter('defocus_lb', -1.5, @isscalar);   %the lower bound of defocus
p.addParameter('defocus_ub', 1.5, @isscalar);    %the upper bound of defocus
p.addParameter('num_defocus', 101, @(x)(isscalar(x) && x >= 51));   %the number of defocus within the range
p.addParameter('whichEye', 'right', @(x)(ischar(x) && (ismember(x, {'left', 'right'}))));  
p.addParameter('zernike_jIndex', 3:14, @(x)(isnumeric(x) && (max(x) <= 15) && (min(x) >= 1)));
p.addParameter('verbose', true, @islogical);
p.addParameter('visualization',true, @islogical);
p.addParameter('saveFigs', false, @islogical);

parse(p, varargin{:});
measuredPupilMM = p.Results.measuredPupilMM;
calcPupilMM     = p.Results.calcPupilMM;
measuredWvl     = p.Results.measuredWvl;
infocusWvl_slc  = p.Results.infocusWvl;
optics_dataset  = p.Results.optics_dataset;
addedDefocus_lb = p.Results.defocus_lb;
addedDefocus_ub = p.Results.defocus_ub;
lenDefocus      = p.Results.num_defocus;
whichEye        = p.Results.whichEye;
jIndex          = p.Results.zernike_jIndex;
verbose         = p.Results.verbose; 
visualization   = p.Results.visualization; 
saveFigs        = p.Results.saveFigs;

%% Do some basic checks for the inputs
assert(calcPupilMM <= measuredPupilMM, ['The pupil size for calculation ',...
    'has to be smaller than the measured pupil size!']);
assert(addedDefocus_ub > addedDefocus_lb, ['The upper bound has to be ',...
    'smaller than the lower bound for the range of added defocus!']);

%% We want human lca for this tutorial
humanlca = true;

%% Set some parameters for the focus search
%
% Dave a range of in-focus wavelengths and set them to calc wvl
infocusWvl   = 450:10:650;
calcWvl      = infocusWvl; 
lenCalcWvl   = length(calcWvl);

% Find the index corresponding to the selected in-focus wavelength
idx_wvl      = find(calcWvl == infocusWvl_slc);

% Have a range of defocus value (diopters)
addedDefocus = linspace(addedDefocus_lb,addedDefocus_ub,lenDefocus);

%% Load the data from the selected subject
wvf_0 = wvfLoadWavefrontOpticsData(...
    'wvfZcoefsSource', optics_dataset, 'jIndex', jIndex, ...
    'whichEye', whichEye, 'eccentricity', [0,0], ...
    'whichGroup', Subj, 'verbose', false);
%grab the zernike polynomial corresponding to defocus
defocus_wvf0 = wvfGet(wvf_0, 'zcoeffs', {'defocus'});

% Get the diffraction-limited optics
wvf_diffractionLimited = wvfCreate(...
    'zcoeffs', zeros(1,15), ...
    'measured pupil', measuredPupilMM,...
    'measured wavelength', measuredWvl, ...
    'calc pupil size', calcPupilMM,...
    'calc wavelengths', infocusWvl_slc);
if (humanlca)
    wvf_diffractionLimited = wvfSet(wvf_diffractionLimited,'lcaMethod','human');
else
    wvf_diffractionLimited = wvfSet(wvf_diffractionLimited,'lcaMethod','none');
end
wvf_diffractionLimited = wvfCompute(wvf_diffractionLimited);

% Get the psf given the diffraction-limited optics
psf_diffractionLimited = wvf_diffractionLimited.psf;

%% Main code starts here
%
% We loop through the code twice. In the 1st time, we find the
% optimal amount of defocus using both grid search and fmincon. In the
% 2n time, we add the optimal amount of defocus we just calculated and
% repeat the same thing. This serves as a sanity check. If the code is doing
% what it is supposed to do, then the optimal defocus we find in the 2nd
% time should be 0.

% Initialize the variable that stores the optimal defocus
[optDefocus_diopters_gridSearch, optDefocus_diopters_fmincon] = deal(NaN(1,2));
for counter = 1:2
    % save all the point spread functions with varying in-focus wavelength
    % and added amount of defocus
    psf = cell(lenCalcWvl, lenDefocus);  
    for l = 1:lenCalcWvl
        % copy the wavefront parameter structure
        wvf_l = wvf_0;

        % Change the in-focus wavelength
        wvf_l.wls = calcWvl(l);

        % If it's 550 nm, save it for later use
        if l == idx_wvl; wvf_550 = wvf_l; end

        % loop through each defocus
        for m = 1:lenDefocus
            % If we haven't found the optimal defocus that leads to the
            % sharpest PSF, we make no adjustment because
            % optDefocus_diopters_gridSearch(1) = NaN
            % if we have found that, optDefocus_diopters_gridSearch(2) will
            % not be NaN.
            if counter == 1
                addedDefocus_updated = addedDefocus(m);
            else
                addedDefocus_updated = addedDefocus(m) + ...
                    optDefocus_diopters_gridSearch(counter-1);
            end

            % Compute the psf with added defocus
            [~, psf{l,m}] = psf_addedDefocus(addedDefocus_updated, ...
                defocus_wvf0, wvf_l, measuredPupilMM, humanlca);
        end
    end 

    %% Compute the optimal defocus 
    % Call func compactness_maxVal to get the peak of all the PSF's
    maxPSF = compactness_maxVal(calcWvl, addedDefocus, psf);
    
    % Alternatively, we can call compactness_StrehlRatio, which essentially
    % returns the same results as compactness_StehlRatio but in different
    % units. Strehl ratio is between 0 and 1. 1 means the person has a
    % perfect vision.
    % strehlRatio = compactness_StrehlRatio(calcWvl, addedDefocus, psf, ...
    %     psf_diffractionLimited);
    
    % We do not vary in-focus wavelength here, but instead just fix it at
    % 550nm since it's a reasonable choice
    [maxPSF_wvl550, ~, optDefocus_diopters_gridSearch(counter), maxPSFVal_wvl550] = ...
        compactness_maxVal([calcWvl(idx_wvl)], addedDefocus, psf(idx_wvl,:));

    [strehlRatio_wvl550, ~, ~, maxStrehlRatio_wvl550] = ...
        compactness_StrehlRatio([calcWvl(idx_wvl)], addedDefocus, ...
        psf(idx_wvl,:), psf_diffractionLimited);

    % Or we could call fmincon and have it to find the defocus that
    % leads to the most compact PSF
    if counter == 1; defocus_adj_temp = 0;
    else; defocus_adj_temp = optDefocus_diopters_fmincon(counter-1); end
    [optDefocus_diopters_fmincon(counter), ~] = findOptDefocus_fmincon(...
        defocus_wvf0, defocus_adj_temp, wvf_550, measuredPupilMM, ...
        addedDefocus_lb, addedDefocus_ub);

    %% Display the results and compare
    if verbose
        if counter == 1; disp('Before adjustment:'); 
        else; disp('After adjustment:'); end
        fprintf(['The amount of defocus that leads to the most compact PSF is:\n',...
            '%.3f diopters (by grid search); %.3f diopters (by fmincon)\n'],...
            optDefocus_diopters_gridSearch(counter), ...
            optDefocus_diopters_fmincon(counter));
    end
    
    %% VISUALIZATION
    if visualization
        % Visualize the peak of PSF
        plotMaxPSF(calcWvl, addedDefocus, maxPSF, maxPSF_wvl550,...
            strehlRatio_wvl550, optDefocus_diopters_gridSearch(counter),...
            maxPSFVal_wvl550, maxStrehlRatio_wvl550, measuredPupilMM,...
            saveFigs);
        
        % Fix defocus to be 0, and visualize the PSF with selected calc wavelengths 
        contourPSF(calcWvl, addedDefocus, psf, saveFigs);
    end
end

%do some simple checks
%the optimal amount of defocus (found either by grid search or fmincon)
% after an adjustment is implemented should be 0
assert(abs(optDefocus_diopters_gridSearch(end)) <= 1e-2 && ...
    abs(optDefocus_diopters_fmincon(end)) <= 1e-2, ['After the adjustment, ',...
    'the psf function should be the most compact when 0 diopter is added.']);
%the results returned by grid search and fmincon should be very similar
assert(abs(optDefocus_diopters_gridSearch(end-1) - ...
    optDefocus_diopters_fmincon(end-1)) < diff(addedDefocus(1:2)),...
    ['The two approaches for finding the optimal defocus (grid search and ',...
    'fmincon) should return almost the same result!']);
%return the optimal defocus
defocus_adj_diopters = optDefocus_diopters_fmincon(1); 
defocus_adj_microns = wvfDefocusDioptersToMicrons(-defocus_adj_diopters, measuredPupilMM);

end

%% HELPING FUNCTIONS
function [maxPSF, optCalcW, optDefocus, maxVal] = compactness_maxVal(...
    calcW, defocus, PSF)
    %{
    This function finds the optimal amount of added defocus that would lead
    to the most compact PSF using grid search.
    ---- Inputs ----
    calcW: a range of in-focus wavelengths
    defocus: a range of added defocus
    PSF: the point spread function 
    ---- Outputs ----
    maxPSF: the peak psf for each given calcW and defocus
    optCalcW & optDefocus: the optimal combination of calcW and defocus
        that together leads to the most compact PSF 
    maxVal: the maximum value of maxPSF (i.e., the read-out value
        corresponding to optCalcW & optDefocus)
    %}

    %initialize 
    maxPSF = NaN(length(calcW), length(defocus));
    %loop through all the combinations of the in-focus wvl and added
    %defocus
    for l = 1:length(calcW)
        for d = 1:length(defocus)
            %get the peak
            maxPSF(l,d) = max(PSF{l,d}(:));
        end
    end
    %find the in-focus wvl and added defocus that correspond to the
    %greatest peak 
    [maxVal, idx] = max(maxPSF(:));
    [row, col]    = ind2sub([length(calcW), length(defocus)], idx);
    optCalcW      = calcW(row);
    optDefocus    = defocus(col);
end

function [strehlRatio, optCalcW, optDefocus, maxStrehlRatio] = ...
    compactness_StrehlRatio(calcW, defocus, PSF, PSF_diffractionLimited)
    %{
    This function finds the optimal amount of added defocus that would lead
    to the biggest strehl ratio using grid search.
    ---- Inputs ----
    calcW: a range of in-focus wavelengths
    defocus: a range of added defocus
    PSF: the point spread function 
    PSF_diffractionLimited: the psf without aberrations
    ---- Outputs ----
    maxPSF: the peak psf for each given calcW and defocus
    optCalcW & optDefocus: the optimal combination of calcW and defocus
        that together leads to the biggest strehl ratio
    maxStrehlRatio: the maximum value of maxPSF (i.e., the read-out value
        corresponding to optCalcW & optDefocus)
    %}

    %initialize
    strehlRatio = NaN(length(calcW), length(defocus));
    %get the peak of the PSF given diffraction-limited optics
    PSF_peak = max(PSF_diffractionLimited{1}(:));
    %loop through all the combinations of the in-focus wvl and added
    %defocus
    for l = 1:length(calcW)
        for d = 1:length(defocus)
            %compute the strehl ratio
            strehlRatio(l,d) = max(PSF{l,d}(:))./PSF_peak;
        end
    end
    %find the in-focus wvl and added defocus that correspond to the
    %greatest peak 
    [maxStrehlRatio, idx] = max(strehlRatio(:));
    [row, col]            = ind2sub([length(calcW), length(defocus)], idx);
    optCalcW              = calcW(row);
    optDefocus            = defocus(col);
end

function [peakPSF, psf] = psf_addedDefocus(defocus_added, defocus_wvf,...
    wvf, measuredPupilMM, humanlca)
    %{
    This function returns the peak of the psf and the whole psf given
    different defocus. 
    ---- Inputs ----
    defocus_added: the amount of defocus (in diopters) we would like to add
    defocus_wvf0: the original amount of defocus (the 5th zernike
        polynomial)
    wvf: the input wavefront parameters
    measuredPupilMM: The pupil size for the measured wavefront aberration
    ---- Outputs ----
    psf: computed psf with added defocus
    peakPSF: the peak of the psf
    %}

    %convert diopter to micron
    lcaMicrons   = wvfDefocusDioptersToMicrons(-defocus_added, measuredPupilMM);
    %then we just add the defocus (microns)
    defocus_wvf  = defocus_wvf + lcaMicrons;
    %replace the defocus value with the new one
    wvf          = wvfSet(wvf, 'zcoeffs', defocus_wvf, {'defocus'});
    %compute the point spread function
    if (humanlca)
        wvf          = wvfSet(wvf,'lcaMethod','human');
    else
        wvf          = wvfSet(wvf,'lcaMethod','none');
    end
    wvf          = wvfCompute(wvf);
    %get the psf
    psf          = wvf.psf{1};
    %get the peak
    peakPSF      = max(psf(:));
end

function [optDefocus, maxPSF] = findOptDefocus_fmincon(defocus_wvf, ...
    defocus_adjustment, wvf, measuredPupilMM, lb, ub)
    %{
    This function finds the optimal amount of added defocus that would lead
    to the most compact PSF using fmincon.
    ---- Inputs ----
    defocus_wvf: the original defocus (without any adjustment and any added
        defocus)
    defocus_adjustment: the optimal amount of defocus, which is 0 if it has
        not been found yet.
    wvf: the wavefront parameters
    measuredPupilMM: measured pupil size (mm)
    lb: the lower bound of added defocus
    ub: the upper bound of added defocus
    ---- Outputs ----
    optDefocus: the optimal amount of added defocus 
    maxPSF: the maximum peak PSF that corresponds to the optiaml amount of
        added defocus
    %}

    %call psf_addedDefocus to get the peak psf
    %a minus sign is added because fmincon finds the minimum
    humanlca = true;
    invertedPeakPSF = @(d) -psf_addedDefocus(d + defocus_adjustment, ...
                defocus_wvf, wvf, measuredPupilMM,humanlca);
    %have different initial points to avoid fmincon from getting stuck at
    %some places
    N_runs  = 20;
    init    = rand(1,N_runs).*(ub-lb) + lb;
    options = optimoptions(@fmincon, 'MaxIterations', 1e5, 'Display','off');
    [optDefocus_n, minNegPSF_n] = deal(NaN(1, N_runs));
    for n = 1:N_runs
        %use fmincon to search for the optimal defocus
        [optDefocus_n(n), minNegPSF_n(n)] = fmincon(invertedPeakPSF, init(n), ...
            [],[],[],[],lb,ub,[],options);
    end
    %find the index that corresponds to the minimum value
    [~,idx_min] = min(minNegPSF_n);
    %find the corresponding optimal focus that leads to the highest peak of
    %the psf's
    maxPSF      = -minNegPSF_n(idx_min);
    optDefocus  = optDefocus_n(idx_min);
end

%% PLOTTING FUNCTIONS
function plotMaxPSF(calcWvl, addedDefocus, maxPSF, maxPSF_slc, StrehlR_slc,...
    optDefocus_wvl550,maxPSFVal_slc, maxStrehlR_slc, measuredPupilMM, saveFigs)
    %convert it to microns
    lcaMicrons = wvfDefocusDioptersToMicrons(-addedDefocus, measuredPupilMM);
    c = [143,188,143]./255;
    cmap = [linspace(c(1), 1, 255)', linspace(c(2), 1, 255)',...
        linspace(c(3), 1,255)'];
    
    figure
    subplot(3,1,[1,2])
    imagesc(addedDefocus, calcWvl, maxPSF); 
    % c= colormap(sky); colormap(flipud(c)); colorbar
    colormap(cmap); colorbar;%colorbar('location','northoutside');
    xticks(addedDefocus(1:25:end));
    xticksRow1 = string(round(addedDefocus(1:25:end),1));
    xticksRow2 = string(round(lcaMicrons(1:25:end),1));
    xticksArray = [xticksRow1; xticksRow2];
    xticklabels(strtrim(sprintf('%s\\newline%s\n', xticksArray{:})));
    % xlabel(sprintf('The added defocus (1st row: diopters; 2nd row: microns)'));
    ylabel('In-focus wavelength (nm)'); axis square;
    title('The peak of PSF');
    set(gca,'FontSize',15);
    
    subplot(3,1,3)
    yyaxis left
    plot(addedDefocus, maxPSF_slc,'-','LineWidth',3); hold on
    plot([optDefocus_wvl550,optDefocus_wvl550],[0, maxPSFVal_slc],'k--','LineWidth',2);
    plot(optDefocus_wvl550, maxPSFVal_slc, 'k*','MarkerSize',10);
    text(optDefocus_wvl550+0.1, maxPSFVal_slc, ...
        sprintf('Defocus = %.2f diopters', optDefocus_wvl550),'FontSize',12);  
    box off
    xlim([addedDefocus(1), addedDefocus(end)]);
    xticks(addedDefocus(1:25:end)); yticks(0:0.005:0.1);
    xticklabels(strtrim(sprintf('%s\\newline%s\n', xticksArray{:})));
    ylabel('The peak of PSF'); ylim([0, maxPSFVal_slc*1.3]);
    yticks(round([0, maxPSFVal_slc*0.6, maxPSFVal_slc*1.2],3));
    xlabel(sprintf('The added defocus \n1st row: diopters; 2nd row: microns'));

    yyaxis right
    plot(addedDefocus-0.01, StrehlR_slc,'-','LineWidth',2);   
    yticks(round([0, maxStrehlR_slc*0.6, maxStrehlR_slc*1.2],2));
    ylim([0, maxStrehlR_slc*1.3]);
    ylabel('Maximum Strehl ratio');

    set(gca,'FontSize',15);
    set(gcf,'Units','normalized','Position',[0,0,0.2,0.65])
    if saveFigs
        set(gcf,'PaperUnits','centimeters','PaperSize',[20 28]);
        formattedDateTime = string(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));
        saveas(gcf, strcat('peakPSF_optimizingDefocus_', formattedDateTime,'.pdf'))
    end
end

function contourPSF(calcWvl, addedDefocus, psf, saveFigs)
    lb_idx          = 77;
    ub_idx          = 127;
    lb_adjusted     = ceil((ub_idx - lb_idx)/2);
    ub_adjusted     = ub_idx - lb_idx + 1;
    slc_defocus     = 0;
    slc_calcWvl     = [500,550,600];
    idx_slc_calcWvl = arrayfun(@(idx) find(slc_calcWvl(idx) == calcWvl),...
                        1:length(slc_calcWvl));
    idx_slc_defocus = find(addedDefocus == slc_defocus);
    
    figure
    for i = 1:length(slc_calcWvl)
        subplot(length(slc_calcWvl),1,i)
        contour(psf{idx_slc_calcWvl(i), idx_slc_defocus}(lb_idx:ub_idx,...
            lb_idx: ub_idx)); hold on
        plot([lb_adjusted, lb_adjusted], [0, ub_adjusted],'r-');
        plot([0, ub_adjusted], [lb_adjusted, lb_adjusted],'r-');
        grid on; axis square; 
        xticks([]); yticks([]);
        title(sprintf('Calc wvl = %.d (nm)', slc_calcWvl(i)));
        set(gca,'FontSize',15);
    end
    set(gcf,'Units','normalized','Position',[0,0, 0.15, 0.65]);
    set(gcf,'PaperUnits','centimeters','PaperSize',[20 28]);
    if saveFigs
        set(gcf,'PaperUnits','centimeters','PaperSize',[20 28]);
        formattedDateTime = string(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));
        saveas(gcf, strcat('contourPSF_optimizingDefocus_', formattedDateTime,'.pdf'))
    end
end

