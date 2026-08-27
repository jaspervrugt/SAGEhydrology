function must_exist(fpth,dispname,caller)
%MUST_EXIST Check that a required file exists.
%
% SYNOPSIS:
%   must_exist(fpth,dispname)
%   must_exist(fpth,dispname,caller)
%
% INPUT:
%   fpth       Full path to required file
%   dispname   Display name shown in error message
%   caller     Optional caller name, e.g. 'read_attr_DE'
%
% OUTPUT:
%   None. The function returns silently if the file exists and throws a
%   region-specific error if the file is missing.
%
% NOTES:
%   This helper is used by regional CAMELS read_attr_XX readers to check
%   whether required input files are present. Passing the reader name as
%   caller preserves reader-specific error identifiers, for example:
%
%       read_attr_DE:missingFile
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, June 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 2 ...
            || isempty(dispname)
        dispname = fpth;
    end

    if nargin < 3 ...
            || isempty(caller)
        caller = mfilename;
    end

    caller = char(string(caller));

    if ~isfile(fpth)
        error([caller ':missingFile'], ...
            ['      Error:%s: ' ...
             'Required file not found: %s\n' ...
             '      Expected path: %s'], ...
            caller,dispname,fpth);
    end

end