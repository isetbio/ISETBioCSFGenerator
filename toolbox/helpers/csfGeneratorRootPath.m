function rootPath = csfGeneratorRootPath()
% Return the path to the root ISETBioCSFGenerator directory
%
% Syntax:
%   rootPath = csfGeneratorRootPath;
%
% Description:
%    This points at the top level of the ISETBioCSFGenerator tree on the Matlab path.
%
%    Examples are included within the code.
%
%    This function works by using the function mfilename to
%    find itself, and then walks back up the result to the top level of
%    isetbio. Thus, you can't move this function within the ISETBioCSFGenerator tree
%    without also adjusting the number of levels in the walk to match
%    where you move it to.
%
% Inputs:
%    None.
%
% Outputs:
%    rootPath - The root directory for ISETBioCSFGenerator
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    isetbioDataPath, isetbioRootPath

% History:
%    12/15/23  dhb, fh  ISETBioCSFGenerator version.

% Examples:
%{
    fullfile(csfGeneratorRootPath, 'tutorials')
%}

%% Get path to this function and then walk back up to the isetbio root.
pathToMe = mfilename('fullpath');

%% Walk back up the chain
rootPath = fileparts(fileparts(fileparts(pathToMe)));

end
