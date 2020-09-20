% ISETbioAdaPEELocalHook
%
% Template for setting preferences for the ISETbio Adaptive Performance
% Estimator Engine toolbox

% 9/20/2020  npc    Wrote it.

%% Clear prefs
% 
% We use these, clear before setting below.
if (ispref('ISETbioAdaPEE'))
    rmpref('ISETbioAdaPEE');
end


% Root dir
setpref('ISETbioAdaPEE','recipesDir',fullfile(tbLocateToolbox('ISETbioAdaPEE'),'recipes'));
