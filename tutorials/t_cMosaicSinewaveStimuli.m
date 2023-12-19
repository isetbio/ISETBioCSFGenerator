% Demo computing off-axis responses to sinewave stimuli
%
% Description:
% Demonstrate usage of @cMosaic, +PolansOptics to compute cone excitations 
%    to sinewave stimuli and display the cone mosaic, the PSF and the mosaic
%    cone excitations and modulations
%
% See Also:
%   t_cMosaicOffAxisDistortion
%   t_cMosaicRankedSubjectsOptics

% History:
%    07/20/21  NPC  ISETBIO Team, Copyright 2021 Wrote it.

function t_cMosaicSinewaveStimuli

    %% Control saving of figures.
    %
    % We don't want tutorials saving things into the isetbio source tree
    % willy-nilly
    saveFigures = false;
    figureDir = fullfile(csfGeneratorRootPath,'local',mfilename);
    if (saveFigures)
        if (~exist(figureDir,'dir'))
            mkdir(figureDir);
        end
        fprintf('Will save figures into %s\n',figureDir)
    else
        fprintf('Not saving figures. Set saveFigures to true in the source to save\n');
    end

    if (~exist('csfGeneratorRootPath','file'))
        fprintf('This tutorial requires the ISETBioCSFGenerator repo on your path\n');
        fprintf('Returning without doing anything.\n');
        return;
    end

    whichEye = 'right eye';             % choose between {'right eye', 'left eye'}
    opticsDataBase = 'Artal2012';       % choose between {'Polans2015', 'Artal2012'}
    subjectRankOrder = 5;
    
    % Mosaic size in degrees
    mosaicSizeDegs = [1.5 1.0]; 
    
    % Mosaic eccentricities
    horizontalEccDegs = [0 10];
    
    for iEcc = 1:numel(horizontalEccDegs)
        % Mosaic ecc
        mosaicEccDegs = [horizontalEccDegs(iEcc) 0];

        if (iEcc == 1)
            % Cone excitations activation range
            activationRange = [0 250];
        else
            % Cone excitations activation range
            activationRange = [0 1000];
        end
        
        % Run simulation (low frequency)
        theConeMosaic = [];
        stimulusSpatialFrequencyCPD = 6;
        [theConeMosaic, thePSFData, ...
         noiseFreeConeMosaicTestStimulusActivation, ...
         noisyConeMosaicActivationInstances, ...
         coneMosaicNullStimulusActivation] = simulateCondition(whichEye, ...
              mosaicEccDegs, mosaicSizeDegs, opticsDataBase, subjectRankOrder, ...
              stimulusSpatialFrequencyCPD, theConeMosaic);

        % Visualize results
        visualizeResults(theConeMosaic, thePSFData, ...
            noiseFreeConeMosaicTestStimulusActivation, coneMosaicNullStimulusActivation, ...
            stimulusSpatialFrequencyCPD, activationRange, saveFigures, figureDir);

        % Run simulation (high frequency)
        stimulusSpatialFrequencyCPD = 12;
        [theConeMosaic, thePSFData, ...
         noiseFreeConeMosaicTestStimulusActivation, ...
         noisyConeMosaicActivationInstances, ...
         coneMosaicNullStimulusActivation] = simulateCondition(whichEye, ...
              mosaicEccDegs, mosaicSizeDegs, opticsDataBase, subjectRankOrder, ...
              stimulusSpatialFrequencyCPD, theConeMosaic);

        % Visualize results
        visualizeResults(theConeMosaic, thePSFData, ...
            noiseFreeConeMosaicTestStimulusActivation, coneMosaicNullStimulusActivation, ...
            stimulusSpatialFrequencyCPD, activationRange, saveFigures, figureDir);
    end
end

function visualizeResults(cm, thePSFData, noiseFreeConeMosaicTestStimulusActivation, ...
    coneMosaicNullStimulusActivation, stimulusSpatialFrequencyCPD, activationRange, saveFigures, figureDir)
    
    % Retrieve the PSF at 550 nm
    targetWavelength = 550;
    [~,idx] = min(abs(thePSFData.supportWavelength-targetWavelength));
    psf = squeeze(thePSFData.data(:,:,idx));
    
    % Normalize to unit amplitude
    psf = psf/max(psf(:));
    
    % Get spatial support, adjusting for the mosaic's eccentricity 
    % (so that the PSF is plotted at the center of the mosaic)
    psfSupportDegs = thePSFData.supportX/60;
    psfSupportDegsX = psfSupportDegs + cm.eccentricityDegs(1);
    psfSupportDegsY = psfSupportDegs + cm.eccentricityDegs(2);
    
    fontSize = 20;
    
    % Setup figure
    hFig = figure(2); clf;
    set(hFig, 'Position', [10 10 620 900], 'Color', [1 1 1]);
    
    % Percentage of mosaic to be visualized within the inset
    visualizedFraction = 0.07;
    w = cm.sizeDegs(1);
    % Visualize part of the mosaic
    domainVisualizationLimsInset(1:2) = cm.eccentricityDegs(1) + visualizedFraction*w/2*[-1 1];
    domainVisualizationLimsInset(3:4) = cm.eccentricityDegs(2) + visualizedFraction*w/2*[-1 1];
    domainVisualizationTicks = struct('x',  [nan], 'y', [nan]);
    

    % Visualize the entire mosaic
    ax = axes('Position',[0.08 0.57 0.92 0.4]);
    cm.visualize('figureHandle', hFig, ...
        'axesHandle', ax, ...
        'domain', 'degrees', ...
        'labelCones', true, ...
        'visualizedConeAperture', 'geometricArea', ...
        'fontSize', fontSize, ...
        'noYLabel', ~true, ...
        'noXLabel', true, ...
        'plotTitle', 'cone mosaic');
    
    % Inset showing a high-res version of the mosaic and the PSF
    ax = axes('Position', [0.705 0.75 0.25 0.25], 'Color',[1 1 1]);
    cm.visualize('figureHandle', hFig, ...
        'axesHandle', ax, ...
        'domain', 'degrees', ...
        'domainVisualizationLimits', domainVisualizationLimsInset, ...
        'domainVisualizationTicks', domainVisualizationTicks, ...
        'labelCones', ~true, ...
        'visualizedConeAperture', 'geometricArea', ...
        'fontSize', fontSize, ...
        'noYLabel', true, ...
        'noXLabel', true, ...
        'backgroundColor', [0 0 0], ...
        'plotTitle', ' ');
    
    % Co-visualize the PSF
    hold(ax, 'on');
    cmap = brewermap(1024,'reds');
    alpha = 0.5;
    
    contourLineColor = [0.0 0.0 0.0];
    cMosaic.semiTransparentContourPlot(ax, psfSupportDegsX, psfSupportDegsY, psf, 0.05:0.2:0.95, cmap, alpha, contourLineColor);
  

    % Plot the mosaic's activation by the stimulus
    ax = axes('Position',[0.08 0.07 0.92 0.4]);
    cm.visualize('figureHandle', hFig, ...
              'axesHandle', ax, ...
              'activation', noiseFreeConeMosaicTestStimulusActivation, ... 
              'activationRange', activationRange, ...
              'activationColorMap', brewermap(1024, '*greys'), ...
              'labelCones', ~true, ...
              'crossHairsOnFovea', false, ...
              'visualizedConeAperture', 'geometricArea', ...
              'verticalActivationColorbarInside', true, ...
              'colorbarTickLabelColor', [1 0.5 0], ...
              'domain','degrees', ...
              'noYLabel', ~true, ...
              'noXLabel', ~true, ...
              'fontSize', fontSize, ...
              'backgroundColor', [0 0 0], ...  
              'plotTitle',  sprintf('mosaic activation (spatial frequency = %2.0f c/deg)',stimulusSpatialFrequencyCPD));
    
     if (saveFigures)
        NicePlot.exportFigToPDF(fullfile(sprintf('ECC_%2.0fDEGS_SF_%2.0fCPD.pdf', cm.eccentricityDegs(1), stimulusSpatialFrequencyCPD)),hFig, 300);
     end
          
end

function [cm, thePSFData, noiseFreeConeMosaicTestStimulusActivation, ...
          noisyConeMosaicActivationInstances, ...
          coneMosaicNullStimulusActivation] =  simulateCondition(whichEye, ...
          mosaicEccDegs, mosaicSizeDegs, opticsDataBase, subjectRankOrder, stimulusSpatialFrequencyCPD, cm)
      
    % Setup stimulus params
    sParams = struct(...
                'contrast', 0.7, ...
                'chromaDir', [1 1 1], ...
                'meanLuminanceCdPerM2', 100, ...
                'sf', stimulusSpatialFrequencyCPD, ...
                'fovDegs', max(mosaicSizeDegs), ...
                'orientation', 90, ...
                'spatialEnvelope', 'rect', ...  % choose between {'rect', 'disk', 'soft'}
                'spatialEnvelopeRadiusDegs', max(mosaicSizeDegs)/2, ...
                'spatialPhase', 0, ...
                'pixelsNum', 1024, ...
                'minPixelsNumPerCycle', 11, ...
                'spectralSupport', [450:25:700], ...
                'presentationMode', 'flashed', ...
                'duration', 1000/1000, ...
                'warningInsteadOfErrorOnOutOfGamut', false ...
     );
 
    % Compute the test and null stimulus scenes
    [stimulusSceneSequence, theNullStimulusScene, statusReport] = ...
        CSFGeneratorApp.generate.gratingSceneEngine(sParams, []);
    theTestStimulusScene = stimulusSceneSequence{1};
    
    if (isempty(cm))
        % Generate the cone mosaic
        cm = cMosaic(...
            'whichEye', whichEye, ...         
            'sizeDegs', mosaicSizeDegs, ...    
            'eccentricityDegs', mosaicEccDegs, ...  
            'opticalImagePositionDegs', 'mosaic-centered' ...
            );
    end
    
    % Get optics params for given database and subject ranking
    switch (opticsDataBase)
        case 'Polans2015'
            % Obtain subject IDs ranking in decreasing foveal resolution
            rankedSujectIDs = PolansOptics.constants.subjectRanking;
            testSubjectID = rankedSujectIDs(subjectRankOrder);

            % Determine if we need to subtract the subject's central refraction
            subtractCentralRefraction = PolansOptics.constants.subjectRequiresCentralRefractionCorrection(testSubjectID);

        case 'Artal2012'
            % Obtain subject IDs ranking in decreasing foveal resolution
            rankedSujectIDs = ArtalOptics.constants.subjectRanking(whichEye);
            testSubjectID = rankedSujectIDs(subjectRankOrder);
            
            % Determine if we need to subtract the subject's central refraction
            subtractCentralRefraction = ArtalOptics.constants.subjectRequiresCentralRefractionCorrection(whichEye, testSubjectID);
    end
 
    % Generate optics appropriate for the mosaic's eccentricity  
    [oiEnsemble, psfEnsemble] = cm.oiEnsembleGenerate(cm.eccentricityDegs, ...
                    'zernikeDataBase', opticsDataBase, ...
                    'subjectID', testSubjectID, ...
                    'pupilDiameterMM', 3.0, ...
                    'zeroCenterPSF', true, ...
                    'subtractCentralRefraction', subtractCentralRefraction, ...
                    'wavefrontSpatialSamples', 501);
                
    % Retrieve the OI - used for computation
    theOI = oiEnsemble{1};   
    
    % Retrieve the PSF - used for visualization
    thePSFData = psfEnsemble{1};  
 
    % Compute the optical image of the stimulus
    theStimulusOI = oiCompute(theOI,theTestStimulusScene,'pad value','mean');
 
    % Compute the optical image of the null stimulus
    theNullStimulusOI = oiCompute(theOI,theNullStimulusScene,'pad value','mean');
 
    % Compute the mosaic response to the test stimulus
    nInstances = 2;
    [noiseFreeConeMosaicTestStimulusActivation, ...
     noisyConeMosaicActivationInstances] = cm.compute(theStimulusOI, 'nTrials', nInstances);

    % Compute the mosaic response to the null stimulus
    coneMosaicNullStimulusActivation = cm.compute(theNullStimulusOI);
end


    