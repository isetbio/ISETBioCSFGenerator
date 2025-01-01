function opticsParams = generateOpticsParams(opticsType)
% Generate parameters structure to produce the optics
%
% Syntax:
%    opticsParams = generateOpticsParams(opticsType)
%
% Description;
%    Set up some default optics parameters, based on various possible
%    presets.  If a structure is passed, it is assumed to be an actual oi
%    structure, and is just passed through;

if (isstruct(opticsType))
    opticsParams = opticsType;
else
    switch (opticsType)
        case 'BerkeleyAO'
            opticsParams = struct(...
                'wls', 400:10:700, ...
                'type', 'BerkeleyAO', ...
                'pupilDiameterMM', 6.0, ...
                'defocusAmount', 0.1, ...
                'accommodatedWl', 550, ...
                'zCoeffs', zeros(66,1), ...
                'defeatLCA', false ...
                );

        case 'oiEnsembleGenerate'
            % These parameters are understood by cMosaic's
            % oiEnsembleGenerate method.
            opticsParams = struct(...
                'type', 'oiEnsembleGenerate', ...
                'PolansSubject', 10, ...
                'pupilDiameterMM', 3.0 ...
                );

        case 'loadComputeReadyRGCMosaic'
            % These parameters are understood by the RGCMosaic code that
            % loads precomputed things.
            opticsParams = struct(...
                'ZernikeDataBase', 'Polans2015', ...
                'examinedSubjectRankOrder', 6, ...
                'pupilDiameterMM', 3.0, ...
                'analyzedEye', 'right eye', ...
                'refractiveErrorDiopters', 0.0, ...
                'positionDegs', [] ...
                );

        otherwise
            error('Unknown optics type specified');
    end
end

end