function status = csfgenRunExamplesAll
%% Run all the examples in the ISETBioCSFGenerator tree
%
% Syntax:
%     csfgenRunExamplesAll
%
% Description:
%     Run all the examples in the ISETBioCSFGenerator tree,
%     excepthose that contain a line of the form
%     "% ETTBSkip"
%
% Inputs:
%    None.
%
% Outputs:
%    status    - 1 if all examples run OK, 0 otherwise.
%
% Optional key/value pairs:
%    None.
%
% See also:
%   csfgenRunTutorialsAll

% History:
%   10/17/20 dhb  Wrote it.

[~, functionStatus] = ExecuteExamplesInDirectory(tbLocateToolbox('ISETBioCSFGenerator'),'verbose',false);
status = all(functionStatus ~= -1);