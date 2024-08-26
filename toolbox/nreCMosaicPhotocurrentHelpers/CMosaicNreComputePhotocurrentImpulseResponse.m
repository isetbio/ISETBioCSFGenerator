function impulseResponse = CMosaicNreComputePhotocurrentImpulseResponse(timeAxis, interpolateImpulseResponseFlag)

    % Temporary solution until we program photocurrent reponse computation
    % in @cMosaic
    load(fullfile(csfGeneratorRootPath,'osLinearFilters25CDM2.mat'), 'osLinearFilters', 'dt');
    t = (1:size(osLinearFilters,1))*dt;
    
    % impulseResponse = 2*interp1(t, squeeze(osLinearFilters(:,1)), timeAxis);

    if notDefined('interpolateImpulseResponseFlag')
        interpolateImpulseResponseFlag = 'false';
    end

    % QUESTION: What is the factor of 2 below about?
    switch interpolateImpulseResponseFlag
        case 'true'
            impulseResponse = 2*interp1(t, squeeze(osLinearFilters(:,1)), timeAxis);
        case 'false'
            impulseResponse = 2*interp1(t, squeeze(osLinearFilters(:,1)), t);
    end
end