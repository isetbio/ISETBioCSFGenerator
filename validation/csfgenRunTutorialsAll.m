function status = csfgenRunTutorialsAll
% Run all ISETBioCSFGenerator tutorials and collect up status
%
% Syntax
%    status = csfgenRunTutorialsAll
%
% Description
%    Run all of the ISETBioCSFGenerator tutorials that we think should work, and print
%    out a report at the end as to whether they threw errors, or not.
%    Scripts inside of isetbioRootPath/tutorials are run, except that
%    scripts within the directory 'underDevelopment' are skipped.
%
% Inputs:
%    None
%
% Outputs:
%    status    - Set to true if all tutorials ran OK, false otherwise.
% 
% See also:
%

% History:
%   10/17/20  dhb  Wrote this, because we care.

% User/project specific preferences
p = struct(...
    'rootDirectory',            fileparts(which(mfilename())), ...
    'tutorialsSourceDir',       fullfile(tbLocateToolbox('ISETBioCSFGenerator'), 'tutorials') ...                % local directory where tutorial scripts are located
    );

%% List of scripts to be skipped
%
% Anything with this in its path name is skipped.
scriptsToSkip = {...
    'underDevelopment' ...
    'support' ...
    };

%% Use UnitTestToolbox method to do this.
status = UnitTest.runProjectTutorials(p, scriptsToSkip, 'All');

end