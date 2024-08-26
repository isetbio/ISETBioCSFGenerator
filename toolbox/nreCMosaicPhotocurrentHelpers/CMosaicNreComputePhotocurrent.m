function photocurrentResponses = CMosaicNreComputePhotocurrent(coneExcitationResponses, temporalSupportSeconds, noiseFlag, responseSecondsKept, cutoffFlag)
    % % Compute photocurrent responses from the cone excitation responses
    % photocurrentResponses = 0*coneExcitationResponses;
    % 
    % impulseResponse = CMosaicNreComputePhotocurrentImpulseResponse(temporalSupportSeconds);
    
    if notDefined('cutoffFlag'), cutoffFlag = 'false'; end
    % Compute photocurrent responses from the cone excitation responses
    photocurrentResponses = 0*coneExcitationResponses;
    switch cutoffFlag
        case 'true'
            impulseResponse = CMosaicNreComputePhotocurrentImpulseResponse(temporalSupportSeconds, 'true');
        case 'false'
            impulseResponse = CMosaicNreComputePhotocurrentImpulseResponse(temporalSupportSeconds, 'false');
    end
    impulseResponse(1) = 0.0;

    dt = temporalSupportSeconds(2)-temporalSupportSeconds(1);
    
    for iTrial = 1:size(coneExcitationResponses,1)
        parfor iCone = 1:size(coneExcitationResponses,3)
            coneExcitation = squeeze(coneExcitationResponses(iTrial,:,iCone));
            pCurrent = conv(coneExcitation, impulseResponse, 'same');
            if (~strcmp(noiseFlag, 'none'))
                photocurrentResponses(iTrial,:,iCone) = osAddNoise(pCurrent, 'sampTime', dt);
            else
                photocurrentResponses(iTrial,:,iCone) = pCurrent;
            end
            
%             figure(33);
%             subplot(2,1,1);
%             plot(temporalSupportSeconds, coneExcitation);
%             subplot(2,1,2);
%             plot(temporalSupportSeconds, squeeze(photocurrentResponses(iTrial,:,iCone)));
%             pause
        end
    end
    
    % Keep the last specified msec of the photocurrent response
    dt = temporalSupportSeconds(2)-temporalSupportSeconds(1);
    idx = find(temporalSupportSeconds > temporalSupportSeconds(end) +dt - responseSecondsKept);
    fprintf('Computed photocurrent response over a time duration of %d msec\n', numel(temporalSupportSeconds)*dt*1000);
    fprintf('Keeping the last %d mseconds of the response to allow for pCurrent stabilization.\n',  responseSecondsKept*1000);
    photocurrentResponses = photocurrentResponses(:,idx,:);

end