function mode = get_assess_mode(prd,bas)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GET_ASSESS_MODE Determine assessment design from period and basin split
%
% SYNOPSIS:
%  mode = get_assess_mode(prd,bas)
%
%   prd         structure with temporal settings
%    .des        first day evaluation period [OPTIONAL]
%    .dee        last day evaluation period  [OPTIONAL]
%   bas         structure with basin split information
%    .K_e        # evaluation basins [OPTIONAL]
%   mode        OUTPUT: assessment design code
%                1 = train basins | train period/mask
%                2 = train basins | train & evaluation period/mask
%                3 = train + eval basins | train period/mask
%                4 = train + eval basins | train & evaluation period/mask
%
% LOGIC:
%   mode = 1 if there are no evaluation basins and no separate evaluation
%            period/mask is defined
%   mode = 2 if there are no evaluation basins but a separate evaluation
%            period/mask is defined
%   mode = 3 if there are evaluation basins but no separate evaluation
%            period/mask is defined [or eval prd/mask = train prd/mask]
%   mode = 4 if there are evaluation basins and a separate evaluation
%            period/mask is defined
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025                                 %
% University of California Irvine                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hasEvalPeriod = false;
if isfield(prd,'des') ...
        && ~isempty(prd.des) ...
    && isfield(prd,'dee') ...
    && ~isempty(prd.dee)
    % A separate evaluation period only exists if:
    %   1) training dates also exist, and
    %   2) eval dates differ from train dates
    if isfield(prd,'dts') ...
            && ~isempty(prd.dts) ...
            && isfield(prd,'dte') ...
            && ~isempty(prd.dte)
        dts = double(prd.dts(:).');
        dte = double(prd.dte(:).');
        des = double(prd.des(:).');
        dee = double(prd.dee(:).');
        hasEvalPeriod = ~(isequal(des,dts) ...
            && isequal(dee,dte));
    else
        % fallback: if training dates are not available, keep old behavior
        hasEvalPeriod = true;
    end

end

K_e = 0;
if nargin >= 2 ...
        && isstruct(bas)
    if isfield(bas,'K_e') ...
            && ~isempty(bas.K_e)
        K_e = double(bas.K_e);
    end
end

if K_e > 0 ...
        && hasEvalPeriod
    mode = 4;
elseif K_e > 0
    mode = 3;
elseif hasEvalPeriod
    mode = 2;
else
    mode = 1;
end

end
