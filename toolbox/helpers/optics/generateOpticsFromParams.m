function theOptics = generateOpticsFromParams(opticsParams,theMosaic)

switch (opticsParams.type)
    case 'opticalimage'
        % It is already and optical image, we just pass it on through.
        theOptics = opticsParams;

    case 'oiEnsembleGenerate'
        % Generate wavefront optics appropriate for the mosaic's eccentricity
        oiEnsemble = theConeMosaic.oiEnsembleGenerate(noiseFreeComputeParams.coneMosaicParams.eccDegs, ...
            'zernikeDataBase', 'Polans2015', ...
            'subjectID', noiseFreeComputeParams.opticsParams.PolansSubject, ...
            'pupilDiameterMM', noiseFreeComputeParams.opticsParams.pupilDiameterMM);
        theOptics = oiEnsemble{1};
        clear oiEnsemble

    case 'BerkeleyAO'
        % Set up wavefront optics object
        %
        % Compute pupil function using 'no lca' key/value pair to turn off LCA.
        % You can turn it back on to compare the effect.
        %
        % Deal with best focus by specifying that the wavefront parameters
        % were measured at the wavelength we want to say is in focus. This
        % is a little bit of a hack but seems OK for the diffraction limited case
        % we're using here.
        wvfP = wvfCreate('calc wavelengths', neuralResponseParamsStruct.opticsParams.wls, ...
            'zcoeffs', neuralResponseParamsStruct.opticsParams.zCoeffs, ...
            'name', sprintf('humanAO-%d', neuralResponseParamsStruct.opticsParams.pupilDiameterMM));
        wvfP = wvfSet(wvfP, 'measured wavelength', neuralResponseParamsStruct.opticsParams.accommodatedWl);
        wvfP = wvfSet(wvfP, 'measured pupil size', neuralResponseParamsStruct.opticsParams.pupilDiameterMM);
        wvfP = wvfSet(wvfP, 'calc pupil size', neuralResponseParamsStruct.opticsParams.pupilDiameterMM);
        wvfP = wvfSet(wvfP,'zcoeffs', neuralResponseParamsStruct.opticsParams.defocusAmount, 'defocus');
        if (~neuralResponseParamsStruct.opticsParams.defeatLCA)
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
        theOptics = oiSet(theOptics, 'optics fnumber', focalLengthMM/neuralResponseParamsStruct.opticsParams.pupilDiameterMM);

    case
    otherwise
        error('Unknown opticsType specified');
end

end