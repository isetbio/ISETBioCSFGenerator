function fixationalEMObj = generateFixationalEyeMovements(stimDurationSeconds, frameDurationSeconds, nTrials, theConeMosaic)
    % Initialize
    fixationalEMObj = fixationalEM;              % Instantiate a fixationalEM object
    fixationalEMObj.microSaccadeType = 'none';   % No microsaccades, just drift
    trainFixationalEMObj.randomSeed = -1;
    
    % Compute number of eye movements
    eyeMovementsPerTrial = stimDurationSeconds/frameDurationSeconds;

    % Generate the em sequence for the passed cone mosaic,
    % which results in a time step equal to the integration time of theConeMosaic
    fixationalEMObj.computeForCmosaic(...
        theConeMosaic, eyeMovementsPerTrial,...
        'nTrials' , nTrials);

    % Set the fixational eye movements into the cone mosaic
    theConeMosaic.emSetFixationalEMObj(fixationalEMObj);
end