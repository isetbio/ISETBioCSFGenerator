% ISETCSFGeneratorLocalHook
%
% Template for setting preferences for the ISETBio CSF Generator Toolbox
%
% 9/20/2020  npc    Wrote it.

%% Clear prefs
% 
% We use these, clear before setting below.
if (ispref('ISETbioAdaPEE'))
    rmpref('ISETbioAdaPEE');
end
if (ispref('ISETBioCSFGenerator'))
    rmpref('ISETBioCSFGenerator');
end

% Root dir
setpref('ISETBioBioCSFGenerator','recipesDir',fullfile(tbLocateToolbox('ISETBioCSFGenerator'),'recipes'));
