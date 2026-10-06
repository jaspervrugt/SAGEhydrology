function lines = sage_log_detail_lines(message,width)
%SAGE_LOG_DETAIL_LINES Wrap setup details with the standard five-space indent.
    if nargin<2, width=80; end
    remaining=strtrim(char(message)); lines={};
    while numel(remaining)>width
        window=remaining(1:width);
        cut=find(window==';',1,'last');
        if isempty(cut), cut=find(isspace(window),1,'last'); end
        if isempty(cut), cut=width; end
        lines{end+1}=['     ' strtrim(remaining(1:cut))]; %#ok<AGROW>
        remaining=strtrim(remaining(cut+1:end));
    end
    if ~isempty(remaining), lines{end+1}=['     ' remaining]; end
end
