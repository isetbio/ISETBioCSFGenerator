function theOptics = generateOpticsFromParams(opticsParams,theMosaic,mosaicParams)

if (isstruct(opticsParams))
    type = opticsParams.type;
else
    type = opticsParams;
end

switch (type)
    case 'opticalimage'
        % It is already and optical image, we just pass it on through.
        if (~strcmp(opticsParams.type,'opticalimage'))
            error('Something went wrong with optics logic');
        end
        theOptics = opticsParams;

    case 'oiEnsembleGenerate'
        % Check that this will work
        if (~isa(theMosaic,'cMosaic'))
            error('Generating optics with ''oiEnsembleGenerate'' requires theMosaic be a cMosaic');
        end
        
        % Generate wavefront optics appropriate for the mosaic's eccentricity
        oiEnsemble = theMosaic.oiEnsembleGenerate(mosaicParams.eccDegs, ...
            'zernikeDataBase', 'Polans2015', ...
            'subjectID', opticsParams.PolansSubject, ...
            'pupilDiameterMM', opticsParams.pupilDiameterMM);
        theOptics = oiEnsemble{1};
        clear oiEnsemble

    case 'BerkeleyAO'
        % Set up wavefront optics object directly
        %
        % Compute pupil function using 'no lca' key/value pair to turn off LCA.
        % You can turn it back on to compare the effect.
        %
        % Deal with best focus by specifying that the wavefront parameters
        % were measured at the wavelength we want to say is in focus. This
        % is a little bit of a hack but seems OK for the diffraction limited case
        % we're using here.
        wvfP = wvfCreate('calc wavelengths', opticsParams.wls, ...
            'zcoeffs', opticsParams.zCoeffs, ...
            'name', sprintf('humanAO-%d', opticsParams.pupilDiameterMM));
        wvfP = wvfSet(wvfP, 'measured wavelength', opticsParams.accommodatedWl);
        wvfP = wvfSet(wvfP, 'measured pupil size', opticsParams.pupilDiameterMM);
        wvfP = wvfSet(wvfP, 'calc pupil size', opticsParams.pupilDiameterMM);
        wvfP = wvfSet(wvfP,'zcoeffs', opticsParams.defocusAmount, 'defocus');
        if (~opticsParams.defeatLCA)
            wvfP = wvfSet(wvfP,'lcaMethod','human');
        end

        % Compute pupil function and PSF
        %
        % Whether LCA should be included depends on your apparatus and
        % is controlled by the boolean defeatLCA in the computation of
        % the pupil function.
        wvfP = wvfCompute(wvfP);

        % Generate optical image object from the wavefront object
        theOptics = wvf2oi(wvfP,'humanlens',true);

        % Set the fNumber to correspond to the pupil size
        focalLengthMM = oiGet(theOptics,'focal length')*1000;
        theOptics = oiSet(theOptics, 'optics fnumber', focalLengthMM/opticsParams.pupilDiameterMM);

    case 'loadComputeReadyRGCMosaic'
      
    otherwise
        error('Unknown opticsType specified');
end

end