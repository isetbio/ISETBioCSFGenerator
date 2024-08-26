function impulseResponse = CMosaicNreComputePhotocurrentImpulseResponse(timeAxis, cutoffFlag)

    % Temporary solution until we program photocurrent reponse computation
    % in @cMosaic
    load(fullfile(csfGeneratorRootPath,'osLinearFilters25CDM2.mat'), 'osLinearFilters', 'dt');
    t = (1:size(osLinearFilters,1))*dt;
    
    % impulseResponse = 2*interp1(t, squeeze(osLinearFilters(:,1)), timeAxis);

    if notDefined('cutoffFlag'), cutoffFlag = 'false'; end
    switch cutoffFlag
        case 'true'
            impulseResponse = 2*interp1(t, squeeze(osLinearFilters(:,1)), timeAxis);
        case 'false'
            impulseResponse = 2*interp1(t, squeeze(osLinearFilters(:,1)), t);
    end
end