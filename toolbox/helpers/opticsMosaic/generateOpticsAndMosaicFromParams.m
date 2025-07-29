function [theOptics,theMosaic] = generateOpticsAndMosaicFromParams(opticsParams,theMosaic,mosaicParams)
% Generate optics and mosaic from the parameters, based on parameter types.
%
% Syntax:
%    [theOptics,theMosaic] = generateOpticsAndMosaicFromParams(opticsParams,theMosaic,mosaicParams)
%
% Description:
%    Generate optics and mosaic objects from parameters.
%    You can also pass the mosaic, in which case 

% You can also just pass the mosaic in which case its parameters are
% probably not needed.

% Generate theOptics
if (isempty(theMosaic))
    switch (mosaicParams.type)
        case 'cMosaic'
            % Generate the cMosaic object
            theMosaic = cMosaic(...
                'wave', mosaicParams.wave, ...
                'sizeDegs', mosaicParams.sizeDegs, ...
                'eccentricityDegs', mosaicParams.eccDegs, ...
                'integrationTime', mosaicParams.timeIntegrationSeconds, ...
                'noiseFlag', 'none' ...
                );

        case 'mRGCMosaic'

            mosaicParams
            opticsParams
            
            if (strcmp(opticsParams.type,'loadComputeReadyRGCMosaic'))
                opticsParams.visualizePSFonTopOfConeMosaic = true;
                [theMosaic, theOptics] = mRGCMosaic.loadPrebakedMosaic(mosaicParams, opticsParams);
                theMosaic.visualize();
            else
                dummyOpticsParams = generateOpticsParams('loadComputeReadyRGCMosaic');
                theMosaic = mRGCMosaic.loadPrebakedMosaic(mosaicParams, dummyOpticsParams);

                opticsParams.whichOptics = 'nativeOptics'; % choose one from mRGCMosaic.validOpticsModifications
                opticsParams.customRefractionDiopters = 0;

                theOptics = theMosaic.nativeOI(...
                    'opticsModification', opticsParams.whichOptics, ...
                    'customRefractionDiopters', opticsParams.customRefractionDiopters, ...
                    'visualizePSF', true, ...
                    'visualizedWavelengths', 450:20:650);
                theMosaic.visualize();
            end

            

            pause
            %{
            % - - - - - - -- - OLD - - - - - -- - -
            % Generate the mRGC object
            if (strcmp(opticsParams.type,'loadComputeReadyRGCMosaic'))
                % We were passed optics parameters that the load method knows
                % about.  Use them.
                theMosaic = mRGCMosaic.loadComputeReadyRGCMosaic(...
                    mosaicParams, ...
                    opticsParams, ...
                    mosaicParams.retinalRFmodelParams);

                % Retrieve the optics that were used to optimize
                % theMRGCmosaic we loaded, as the parameters say this is
                % what we want here.
                if (~isempty(theMosaic.theNativeOptics))
                    theOptics = theMosaic.theNativeOptics;
                elseif (~isempty(theMosaic.theCustomOptics))
                    error('Not expecting to use theCustomOptics field here');
                else
                    error('No optics found in the mRGCMosaic object!')
                end

            else
                % We will need to override the optics.  So first we generate
                % dummy optics parameters, generate the mRGCMosiac, then
                % generate the optics we want below.
                dummyOpticsParams = generateOpticsParams('loadComputeReadyRGCMosaic');
                theMosaic = mRGCMosaic.loadComputeReadyRGCMosaic(...
                    mosaicParams, ...
                    dummyOpticsParams, ...
                    mosaicParams.retinalRFmodelParams);

                % We pass a cone mosaic so we can use it if needed. The
                % only parameters we need from the mosaicParams are the
                % eccDegs field, and these exist both for cMosaic and mRGC
                % types so we're good.  
                %
                % This line might be a bit fragile.
                theOptics = generateOpticsAndMosaicFromParams(opticsParams,theMosaic.theConeMosaic,mosaicParams);
            end

            %  - - - - - - - - - END OF OLD - - - - -- -
            %}

        otherwise
            error('Unknown mosaic type pased in mosaicParams');
    end
end

% Generate optics
switch (opticsParams.type)
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