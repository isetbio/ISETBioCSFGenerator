function [defocus_adj, defocus_adj_micons] = t_optimalDefocus(emmetropes_subj,...
    varargin)
%{
The goal of this script is to find the amount of defocus that can lead
to the most compact PSF. Specifically,
1. We pick a in-focus wavelength (550 nm) and set it to be the calculate 
    wavelength 
2. Vary the amount of defocus and convert it to microns
3. Loop through all the defocus
    Add it to the 5th Zernike polynomial
    Compute the PSF
    Find the peak of the PSF / compute the Strehl ratio
4. Find the LCA that corresponds to the highest peak of the PSF / the 
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
p.addParameter('measuredWvl', 550, @(x)(isscalar(x) && floor(x)==ceil(x))); %measured wavelength (nm)
p.addParameter('infocusWvl', 550, @(x)(isscalar(x) && floor(x)==ceil(x)));  %select one in-focus wavelength for visualization
p.addParameter('defocus_lb', -1.5, @isnumber);   %the lower bound of defocus
p.addParameter('defocus_ub', 1.5, @isnumber);    %the upper bound of defocus
p.addParameter('whichEye', 'right', @(x)(ischar(x) && (ismember(x, {'left', 'right'}))));  
p.addParameter('visualization',true, @islogical);

parse(p, varargin{:});
measuredPupilMM = p.Results.measuredPupilMM;
calcPupilMM     = p.Results.calcPupilMM;
measuredWvl     = p.Results.measuredWvl;
infocusWvl_slc  = p.Results.infocusWvl;
addedDefocus_lb = p.Results.defocus_lb;
addedDefocus_ub = p.Results.defocus_ub;
whichEye        = p.Results.whichEye;
visualization   = p.Results.visualization; 

%do some basic checks for the inputs
assert(calcPupilMM < measuredPupilMM, ['The pupil size for calculation ',...
    'has to be smaller than the measured pupil size!']);
assert(addedDefocus_ub > addedDefocus_lb, ['The upper bound has to be ',...
    'smaller than the lower bound for the range of added defocus!']);


%% first pick a subject's optics
%have a range of in-focus wavelength and set them to calc wvl
infocusWvl   = 450:1:650;
calcWvl      = infocusWvl; 
lenCalcWvl   = length(calcWvl);
%select one in-focus wavelength for visualization
idx_wvl      = find(calcWvl == infocusWvl_slc);
%have a range of defocus value (diopters)
lenDefocus   = 101;
addedDefocus = linspace(addedDefocus_lb,addedDefocus_ub,lenDefocus);
%select one defocus for visualizing PSF's with different calc wvls
defocus_slc  = 0;
%convert it to microns
lcaMicrons   = wvfDefocusDioptersToMicrons(-addedDefocus, measuredPupilMM);

%the goal is to find the amount of defocus that makes the PSF as compact as
%possible. Initialize the variable that stores it.
optDefocus_wvl550 = NaN; %for in-focus wavelength = 550 nm

%if we'd like to make an adjustment (i.e., add the amount of defocus to the
%5th zernike polynomial) after we update the variable optDefocus_wvl550
flag_adjustment   = true;

%% Load the data from the selected subject
% wvf_0 = wvfLoadWavefrontOpticsData(...
%     'wvfZcoefsSource', 'JaekenArtal2012', 'jIndex', 3:14, ...
%     'whichEye', 'right', 'eccentricity', [0,0], ...
%     'whichGroup', emmetropes_subj, 'verbose', false);
wvf_0 = wvfLoadWavefrontOpticsData(...
    'wvfZcoefsSource', 'Polans2015', 'jIndex', 3:14, ...
    'whichEye', whichEye, 'eccentricity', [0,0], ...
    'whichGroup', emmetropes_subj, 'verbose', false);
%grab the zernike polynomial corresponding to defocus
defocus_wvf0 = wvfGet(wvf_0, 'zcoeffs', {'defocus'});

% get the diffraction-limited optics
wvf_diffractionLimited = wvfCreate(...
    'zcoeffs',zeros(1,15), ...
    'measured pupil', measuredPupilMM,...
    'measured wavelength', measuredWvl, ...
    'calc pupil size', calcPupilMM,...
    'calc wavelengths', infocusWvl_slc);
wvf_diffractionLimited = wvfComputePSF(wvf_diffractionLimited,[],'nolca',false);
% get the psf given the diffraction-limited optics
psf_diffractionLimited = wvf_diffractionLimited.psf;

%% Main code starts here
while flag_adjustment     
    % save all the point spread functions with varying in-focus wavelength
    % and added amount of defocus
    psf = cell(lenCalcWvl, lenDefocus);    
    for l = 1:lenCalcWvl
        %copy the wavefront parameter structure
        wvf_l = wvf_0;
        %change the in-focus wavelength
        wvf_l.wls = calcWvl(l);
        %if it's 550nm, save it for later use
        if l == idx_wvl; wvf_550 = wvf_l; end

        %loop through each defocus
        for m = 1:lenDefocus
            %if we haven't found the optimal defocus that leads to the
            %sharpest PSF
            if ~isnan(optDefocus_wvl550)
                %then we do the adjustment by adding the defocus
                defocus_added = addedDefocus(m) + optDefocus_wvl550;
                flag_adjustment_temp = false;
            %if we have found the optimal defocus
            else
                %then we just add the defocus (microns)
                defocus_added = addedDefocus(m);
                flag_adjustment_temp = true;
            end
            %compute the psf with added defocus
            [~, psf{l,m}] = psf_addedDefocus(defocus_added, ...
                defocus_wvf0, wvf_l, measuredPupilMM);
        end
    end 
    %% Compute the optimal defocus 
    %Call func compactness_maxVal to get the peak of all the PSF's
    maxPSF = compactness_maxVal(calcWvl, addedDefocus, psf);
    
    %Alternatively, we can call compactness_StrehlRatio, which essentially
    %returns the same results as compactness_StehlRatio but in different
    %units. Strehl ratio is between 0 and 1. 1 means the person has a
    %perfect vision.
    % strehlRatio = compactness_StrehlRatio(calcWvl, addedDefocus, psf, ...
    %     psf_diffractionLimited);
    
    %We do not vary in-focus wavelength here, but instead just fix it at
    %550nm since it's a reasonable choice
    [maxPSF_wvl550, ~, optDefocus_wvl550, maxPSFVal_wvl550] = ...
        compactness_maxVal([calcWvl(idx_wvl)], addedDefocus, psf(idx_wvl,:));
    
    [strehlRatio_wvl550, ~, ~, maxStrehlRatio_wvl550] = ...
        compactness_StrehlRatio([calcWvl(idx_wvl)], addedDefocus, ...
        psf(idx_wvl,:), psf_diffractionLimited);

    %update the flag 
    %if we have already done adjustment by adding the defocus, then we can
    %exit out of the loop by switching flag_adjustment to false
    flag_adjustment = flag_adjustment_temp;

    %Or we could call fmincon and have it to find the defocus that
    %leads to the most compact PSF
    if flag_adjustment; defocus_adj = 0; 
    else; defocus_adj = optDefocus_fmincon; 
        defocus_adj_micons = wvfDefocusDioptersToMicrons(-defocus_adj, measuredPupilMM);
    end
    [optDefocus_fmincon, ~] = findOptDefocus_fmincon(...
        defocus_wvf0, defocus_adj, wvf_550, measuredPupilMM, ...
        addedDefocus_lb, addedDefocus_ub);

    %% display the results and compare
    if flag_adjustment; disp('Before adjustment:'); 
    else; disp('After adjustment:'); end
    fprintf(['The amount of defocus that leads to the most compact PSF is:\n',...
        '%.3f diopters (by grid search); %.3f diopters (by fmincon)\n'],...
        optDefocus_wvl550, optDefocus_fmincon);
    
    %% VISUALIZATION
    if visualization
        % Visualize the peak of PSF
        plotMaxPSF(calcWvl, addedDefocus, lcaMicrons,maxPSF, maxPSF_wvl550,...
            strehlRatio_wvl550, optDefocus_wvl550,maxPSFVal_wvl550,...
            maxStrehlRatio_wvl550, emmetropes_subj, flag_adjustment)
        
        %Fix defocus to be 0, and visualize the PSF with selected calc wavelengths 
        contourPSF(calcWvl, addedDefocus, [500,550,600], defocus_slc, psf,...
            emmetropes_subj, flag_adjustment)
    end
end
end

%% HELPING FUNCTIONS
function [maxPSF, optCalcW, optDefocus, maxVal] = compactness_maxVal(...
    calcW, defocus, PSF)
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

function [peakPSF, psf] = psf_addedDefocus(defocus_added, defocus_wvf0,...
    wvf_m, measuredPupilMM)
    lcaMicrons = wvfDefocusDioptersToMicrons(-defocus_added, measuredPupilMM);
    %then we just add the defocus (microns)
    defocus_wvf1 = defocus_wvf0 + lcaMicrons;
    %replace the defocus value with the new one
    wvf_m = wvfSet(wvf_m, 'zcoeffs', defocus_wvf1, {'defocus'});
    %compute the point spread function
    wvf_m = wvfComputePSF(wvf_m,[],'nolca',false);
    %get the psf
    psf = wvf_m.psf{1};
    %get the peak
    peakPSF = max(psf(:));
end

function [optDefocus, maxPSF] = findOptDefocus_fmincon(defocus_wvf0, ...
    defocus_adjustment, wvf_m, measuredPupilMM, lb, ub)
    peakPSF = @(d) -psf_addedDefocus(d + defocus_adjustment, ...
                defocus_wvf0, wvf_m, measuredPupilMM);
    %have different initial points to avoid fmincon from getting stuck at
    %some places
    N_runs  = 10;
    init    = rand(1,N_runs).*(ub-lb) + lb;
    options = optimoptions(@fmincon, 'MaxIterations', 1e5, 'Display','off');
    [optDefocus_n, minNegPSF_n] = deal(NaN(1, N_runs));
    for n = 1:N_runs
        [optDefocus_n(n), minNegPSF_n(n)] = fmincon(peakPSF, init(n), ...
            [],[],[],[],lb,ub,[],options);
    end
    [~,idx_min] = min(minNegPSF_n);
    maxPSF      = -minNegPSF_n(idx_min);
    optDefocus  = optDefocus_n(idx_min);
end

%% PLOTTING FUNCTIONS
function plotMaxPSF(calcWvl, addedDefocus, lcaMicrons, maxPSF, ...
    maxPSF_slc,StrehlR_slc, optDefocus_wvl550,maxPSFVal_slc,...
    maxStrehlR_slc, sN, flag_adjustment)
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
    title(sprintf('Jaeken & Artal (Subj #%d)',sN));
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
    set(gcf,'PaperUnits','centimeters','PaperSize',[20 28]);
    % if flag_adjustment %before adjustment
    %     saveas(gcf, ['OptimalDefocus_Polans_subj',num2str(sN),'.pdf']);
    % else %after adjustment
    %     saveas(gcf, ['OptimalDefocus_Polans_subj',num2str(sN),'_wAdjustment.pdf']);
    % end
end

function contourPSF(calcWvl, addedDefocus, slc_calcWvl, slc_defocus,...
    psf, sN, flag_adjustment)
    lb_idx          = 77;
    ub_idx          = 127;
    lb_adjusted     = ceil((ub_idx - lb_idx)/2);
    ub_adjusted     = ub_idx - lb_idx + 1;
    idx_slc_calcWvl = arrayfun(@(idx) find(slc_calcWvl(idx) == calcWvl), 1:length(slc_calcWvl));
    idx_slc_defocus = find(addedDefocus == slc_defocus);
    
    figure
    for i = 1:length(slc_calcWvl)
        subplot(length(slc_calcWvl),1,i)
        contour(psf{idx_slc_calcWvl(i), idx_slc_defocus}(lb_idx:ub_idx, lb_idx: ub_idx)); hold on
        plot([lb_adjusted, lb_adjusted], [0, ub_adjusted],'r-');
        plot([0, ub_adjusted], [lb_adjusted, lb_adjusted],'r-');
        grid on; axis square; 
        xticks([]); yticks([]);
        title(sprintf('Calc wvl = %.d (nm)', slc_calcWvl(i)));
        set(gca,'FontSize',15);
    end
    set(gcf,'Units','normalized','Position',[0,0, 0.15, 0.65]);
    set(gcf,'PaperUnits','centimeters','PaperSize',[20 28]);
    % if ~flag_adjustment
    %     saveas(gcf, ['Optics_Polans_sub',num2str(sN),'.pdf']);
    % else
    %     saveas(gcf, ['Optics_Polans_sub',num2str(sN),'_wAdjustment.pdf']);
    % end
end

