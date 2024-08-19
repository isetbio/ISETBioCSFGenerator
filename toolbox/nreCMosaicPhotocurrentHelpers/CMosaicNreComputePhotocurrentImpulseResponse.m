function impulseResponse = cMosaicNreComputePhotocurrentImpulseResponse(timeAxis)

    % Temporary solution until we program photocurrent reponse computation
    % in @cMosaic
    load(fullfile(csfGeneratorRootPath,'osLinearFilters25CDM2.mat'), 'osLinearFilters', 'dt');
    t = (1:size(osLinearFilters,1))*dt;
    
    impulseResponse = 2*interp1(t, squeeze(osLinearFilters(:,1)), timeAxis);
end