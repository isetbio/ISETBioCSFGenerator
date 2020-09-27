function validateNoiseFlags(obj,noiseFlags)

    assert(iscell(noiseFlags), 'noiseFlags must be a cell array');
    
    for idx = 1:length(noiseFlags)
        if ismember(ieParamFormat(noiseFlags{idx}), obj.validNoiseFlags)
            continue;
        elseif contains(ieParamFormat(noiseFlags{idx}), 'rngseed')
            continue;
        else
            fprintf('Valid noise flags are: ''rngSeedSomeInt'', ''%s'' and ''%s''\n', ...
                obj.validNoiseFlags{1}, obj.validNoiseFlags{2});
            error('noiseFlag ''%s'' is not valid', noiseFlags{idx});
        end
        
   end
end

