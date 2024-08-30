function photocurrentResponses = CMosaicNreComputePhotocurrent(coneExcitationResponses, temporalSupportSeconds, noiseFlag, responseSecondsKept, interpolateImpulseResponseFlag)
    % % Compute photocurrent responses from the cone excitation responses
    % photocurrentResponses = 0*coneExcitationResponses;
    % 
    % impulseResponse = CMosaicNreComputePhotocurrentImpulseResponse(temporalSupportSeconds);
    
    % Set default on the flag that determines whether we interpate the
    % impulse response to the timebase of the excitation sequence.
    if notDefined('interpolateImpulseResponseFlag')
        interpolateImpulseResponseFlag = 'true';
    
    end
    % Compute photocurrent responses from the cone excitation responses
    photocurrentResponses = 0*coneExcitationResponses;
    switch interpolateImpulseResponseFlag
        case 'true'
            % This case reads the impuluse response and interpolates it to
            % the time base of the excitations.
            %
            % QUESTIONS:
            % 1) Why is there a factor of 2 in the returned impuluse response?
            % 2) How much error in the calculation do we introduce by doing
            % the calculations at the frame rate of the excitations, as
            % opposed to interpolating them to a finer time scale?
            % 3) How should we automatically compute the impulse response
            % that we want from the excitations, rather than reading them
            % from a file where they have been precomputed for some other
            % purpose?
            impulseResponse = CMosaicNreComputePhotocurrentImpulseResponse(temporalSupportSeconds, 'true');
        case 'false'
            % I DON'T THINK THIS CASE WILL CURRENTLY RUN CORRECTLY,
            % ALTHOUGH IT MIGHT RUN WITHOUT CRASHING.
            impulseResponse = CMosaicNreComputePhotocurrentImpulseResponse(temporalSupportSeconds, 'false');
    end

    % WHY IS THIS HAPPENING. You'd think that would be set correctly when
    % we snagged the impuluse response.
    % impulseResponse(1) = 0.0;

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
    
    % for now, we keep all photocurrent response
    fprintf('Computed photocurrent response over a time duration of %d msec\n', numel(temporalSupportSeconds)*dt*1000);

    % Keep the last specified msec of the photocurrent response
    % dt = temporalSupportSeconds(2)-temporalSupportSeconds(1);
    % idx = find(temporalSupportSeconds > temporalSupportSeconds(end) +dt - responseSecondsKept);
    % fprintf('Keeping the last %d mseconds of the response to allow for pCurrent stabilization.\n',  responseSecondsKept*1000);
    % photocurrentResponses = photocurrentResponses(:,idx,:);

end