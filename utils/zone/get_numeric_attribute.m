function x = get_numeric_attribute(v,vname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GET_NUMERIC_ATTRIBUTE Convert an attribute variable to numeric values.
%
% SYNOPSIS:
%   x = get_numeric_attribute(v,vname)
%
% INPUT:
%   v       Attribute variable from a table column
%   vname   Attribute variable name, used in warning/error messages
%
% OUTPUT:
%   x       Numeric column vector with one value per basin/row. Numeric and
%           logical inputs are converted directly. Text-valued numeric data
%           are converted with str2double. Non-numeric categorical/text
%           values are encoded as integer category IDs, with missing values
%           retained as NaN.
%
% DESCRIPTION:
%   This helper standardizes regional CAMELS attribute columns to numeric
%   vectors for use in SAGE. It supports numeric, logical, string, cellstr,
%   categorical, and mixed text representations. Missing values such as
%   empty strings, NA, NaN, and NULL are converted to NaN. Text attributes
%   that cannot be interpreted numerically are encoded using sorted 
%   category labels and reported with a warning.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, June 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if isnumeric(v) ...
            || islogical(v)
        x = double(v);
        x = x(:);
        return
    end
    
    s = string(v);
    s = strip(s);
    s(ismissing(s) ...
        | strcmpi(s,"nan") ...
        | strcmpi(s,"na") ...
        | strcmpi(s,"null") ...
        | strcmpi(s,"")) = missing;
    
    x = str2double(s);
    
    nonmiss = ~ismissing(s);
    numericOK = ~nonmiss ...
        | isfinite(x);
    if all(numericOK)
        x = x(:);
        return
    end
    
    % Categorical/timing attributes are encoded deterministically by sorted
    % category. They remain selectable for inspection, but are not defaults
    [cats,~,~] = unique(s(nonmiss));
    cats = sort(cats);
    [~,ic2] = ismember(s(nonmiss),cats);
    x = nan(numel(s),1);
    x(nonmiss) = double(ic2);
    
    warning('read_attr_IND:encodedCategorical', ...
        ['      Warning:read_attr_IND: ' ...
        'Encoded categorical attribute ' ...
        '%s using %d categories.'], ...
        vname,numel(cats));

end
