function value = sage_unit_text(value)
%SAGE_UNIT_TEXT Normalize legacy mojibake in displayed unit metadata.
% This changes display text only, including units loaded from old sessions.
    if iscell(value)
        value = cellfun(@sage_unit_text,value,'UniformOutput',false);
        return
    end
    wasChar = ischar(value);
    value = string(value);
    degree = string(char(176));
    value = replace(value,string(char([195 130 194 176])),degree);
    value = replace(value,string(char([194 176])),degree);
    value = replace(value,string(char([194 186])),degree);
    value = replace(value,string(char(186)),degree);
    if wasChar, value = char(value); end
end
