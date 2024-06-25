function presentationDisplay = generatePresentationDisplay(...
    letterSizeDegs, letterSizePixels, spdDataFile, ambientSPDDataFile, plotCharacteristics)

    % Load the ambient SPD
    projectBaseDir = ISETBioJandJRootPath();

    fprintf('Loading ambient SPD from %s\n', fullfile(getpref('ISETBioJandJ','dataDir'),ambientSPDDataFile));
    load(fullfile(getpref('ISETBioJandJ','dataDir'),ambientSPDDataFile), 'spd');
    ambientSPD = spd;
    clear 'spd'
    
    % Check data consistency
    assert(size(ambientSPD, 2) == 2, 'The ambient SPD matrix must be an N x 2 matrix, with the first column being the spectral support');
    if (size(ambientSPD,2) > 2)
        fprintf(2,'\nThe ambient SPD matrix must be an N x 2 matrix, with the first column being the spectral support and the second column being the ambient energy.\n');
        fprintf(2,'The data retrieved from ''%s'', contain %d columns. Ignoring all but the first 2 columns.\n\n', ambientSPDDataFile, size(ambientSPD,2));
    end
    
    ambientWave = ambientSPD(:,1);
    ambientSPD = ambientSPD(:,2);
    ambientSPD = ambientSPD/ (ambientWave(2)-ambientWave(1));
    
    % Load the RGB SPDs
    fprintf('Loading SPDs from %s\n', fullfile(getpref('ISETBioJandJ','dataDir'),spdDataFile));
    load(fullfile(getpref('ISETBioJandJ','dataDir'),spdDataFile), 'spd');
    
    % Check data consistency
    assert(size(spd, 2) == 4, 'The SPD matrix must be an N x 4 matrix, with the first column being the spectral support');
    wave = spd(:,1);
    spd = spd(:,2:4);
    spd = spd / (wave(2)-wave(1));
    assert(size(spd,1) == size(ambientSPD,1), 'The ambient SPD must have the same wavelength entries as the display SPD');
    assert(all(ambientWave == wave), 'The ambient wavelength support must match the wavelength support of the display SPD');
    
    
    presentationDisplay = generateCustomDisplay(...
           'dotsPerInch', 220, ...
           'wavelengthSupportNanoMeters', wave, ...
           'spectralPowerDistributionWattsPerSteradianM2NanoMeter', spd, ...
           'ambientSPDWattsPerSteradianM2NanoMeter', ambientSPD, ...
           'gammaTable', repmat((linspace(0,1,1024)').^2, [1 3]), ...
           'plotCharacteristics', plotCharacteristics);
    
    
    pixelSizeMeters = displayGet(presentationDisplay, 'meters per dot');
    letterSizeMeters = letterSizePixels*pixelSizeMeters;
    desiredViewingDistance = 0.5*letterSizeMeters/(tand(letterSizeDegs/2));
    presentationDisplay = displaySet(presentationDisplay, 'viewing distance', desiredViewingDistance);
    
end