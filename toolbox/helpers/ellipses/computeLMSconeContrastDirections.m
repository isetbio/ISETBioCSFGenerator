function [LMSconeContrasts, examinedDirectionsOnLMplane] = computeLMSconeContrastDirections(...
    rmsLMconeContrast, examinedDirectionsOnLMplane)

    if (isempty(examinedDirectionsOnLMplane))
        examinedDirectionsOnLMplane = [0:10:20 25:5:65 70:10:170];
        examinedDirectionsOnLMplane = cat(2, examinedDirectionsOnLMplane, 180+examinedDirectionsOnLMplane);
    end

    LMSconeContrasts = zeros(3, numel(examinedDirectionsOnLMplane));
    if (numel(rmsLMconeContrast) == numel(examinedDirectionsOnLMplane))
        LMSconeContrasts(1,:) = rmsLMconeContrast .* cosd(examinedDirectionsOnLMplane);
        LMSconeContrasts(2,:) = rmsLMconeContrast .* sind(examinedDirectionsOnLMplane);
    else
        LMSconeContrasts(1,:) = rmsLMconeContrast(1) * cosd(examinedDirectionsOnLMplane);
        LMSconeContrasts(2,:) = rmsLMconeContrast(1) * sind(examinedDirectionsOnLMplane);
    end
end