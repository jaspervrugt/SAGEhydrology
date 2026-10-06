function write_header(mdl,part,details)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WRITE_HEADER Print the selected model and training banner.
%
%  Writes the command-window heading for SITE or SAGE training.
%
% SYNOPSIS:
%   write_header(mdl)
%   write_header(mdl,part)
%   write_header(mdl,part,details)
%
% INPUT ARGUMENTS:
%   mdl             model-selection settings
%    .model          hydrologic model identifier
%    .variant        optional model variant
%   part            optional 'site' or 'sage'; default 'sage'
%   details         optional N-by-2 cell array of metadata label/value pairs
%
% OUTPUT ARGUMENTS:
%   none            text is written to the command window
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Dec. 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 2 ...
            || isempty(part)
        part = 'sage';
    end
    if nargin < 3 || isempty(details)
        details = cell(0,2);
    end
    if ~iscell(details) || size(details,2) ~= 2
        error('write_header:InvalidDetails', ...
            'details must be an N-by-2 cell array of label/value pairs.');
    end
    
    model = mdl.model;
    modelRoot = '';
    if isfield(mdl,'root') && ~isempty(mdl.root)
        modelRoot = mdl.root;
    end
    mname = sage_model_name(model,modelRoot);
    if isfield(mdl,'variant') ...
            && strcmpi(string(mdl.variant),'gchm_ode')
        mname = 'gchm_ode';
    end
    
    if strcmpi(part,'sage')
        modelRow = {'Hydrologic model',char(string(mname))};
        projectMask = ismember(lower(string(details(:,1))), ...
            ["project","title","summary"]);
        if any(strcmpi(string(details(:,1)),'Project'))
            rows = [details(projectMask,:); modelRow; ...
                details(~projectMask,:)];
        else
            rows = [modelRow; details];
        end
    elseif strcmpi(part,'site')
        rows = [{'Individual training',char(string(mname))}; details];
    else
        error('write_header:UnknownPart', ...
            'Unknown entry of variable part: %s',char(string(part)));
    end

    labelWidth = max(15,max(cellfun( ...
        @(x) numel(char(string(x))),rows(:,1))));
    % Keep the timestamped GUI banner compact on smaller screens. Metadata
    % wraps within the fixed box instead of expanding it for long values.
    w = 60;
    lines = local_format_rows(rows,labelWidth,w-4);
    subtitle = 'Process-Based Hydrologic Models across Large Basin Samples';
    indent = '           ';
    % W is the width between the two vertical box edges.  The corner and
    % separator characters replace those edges; they must not add another
    % two columns or the horizontal rules extend beyond the right-hand |.
    border = repmat('─',1,w);
    
    % -----------------------
    % SAGE header
    % -----------------------
    fprintf('\n%s┌%s┐\n',indent,border);
    local_box_line(indent,'S A G E',w,true);
    local_box_line(indent,'',w);
    local_box_line(indent, ...
        'Sensitivity-Aware Learning for Rapid Training of',w,true);
    local_box_line(indent,subtitle,w,true);
    local_box_line(indent,'',w);
    for k = 1:numel(lines)
        local_box_line(indent,lines{k},w);
    end
    local_box_line(indent,'',w);
    local_box_line(indent,'© 2026 Jasper A. Vrugt',w);
    local_box_line(indent,'  All rights reserved',w);
    local_box_line(indent,'',w);
    local_box_line(indent,'✉ jasper@uci.edu',w);
    fprintf('%s└%s┘\n\n\n',indent,border);

end

function lines = local_format_rows(rows,labelWidth,maxWidth)
%LOCAL_FORMAT_ROWS Format and wrap metadata without widening the banner.

    lines = cell(0,1);
    valueWidth = max(1,maxWidth-labelWidth-3);
    for k = 1:size(rows,1)
        label = char(string(rows{k,1}));
        wrapped = local_wrap_text(char(string(rows{k,2})),valueWidth);
        for j = 1:numel(wrapped)
            if j == 1
                lines{end+1,1} = sprintf('%-*s : %s', ...
                    labelWidth,label,wrapped{j}); %#ok<AGROW>
            else
                lines{end+1,1} = sprintf('%-*s   %s', ...
                    labelWidth,'',wrapped{j}); %#ok<AGROW>
            end
        end
    end
end

function lines = local_wrap_text(text,width)
%LOCAL_WRAP_TEXT Wrap plain metadata at word boundaries.

    words = regexp(strtrim(text),'\s+','split');
    if isempty(words) || (isscalar(words) && isempty(words{1}))
        lines = {''};
        return
    end
    lines = cell(0,1);
    current = '';
    for k = 1:numel(words)
        word = words{k};
        % Long dataset identifiers must not push the right border outward.
        if numel(word) > width
            if ~isempty(current)
                lines{end+1,1} = current; %#ok<AGROW>
                current = '';
            end
            while numel(word) > width
                lines{end+1,1} = word(1:width); %#ok<AGROW>
                word = word(width+1:end);
            end
        end
        if isempty(current)
            candidate = word;
        else
            candidate = [current ' ' word];
        end
        if numel(candidate) <= width || isempty(current)
            current = candidate;
        else
            lines{end+1,1} = current; %#ok<AGROW>
            current = word;
        end
    end
    lines{end+1,1} = current;
end

function local_box_line(indent,text,w,centered)
    if nargin < 4, centered = false; end
    text = char(string(text));
    if centered
        textWidth = local_display_width(text);
        padL = floor((w-textWidth)/2);
        padR = w-textWidth-padL;
        text = [repmat(' ',1,padL) text repmat(' ',1,padR)];
    else
        text = ['  ' text];
        padR = max(0,w-local_display_width(text));
        text = [text repmat(' ',1,padR)];
    end
    fprintf('%s|%s|\n',indent,text);
end

function width = local_display_width(text)
%LOCAL_DISPLAY_WIDTH Approximate monospace cells used by MATLAB text areas.
% The envelope is commonly supplied by a fallback font as a two-cell glyph.
    text = char(string(text));
    width = numel(text) + count(text,'✉');
end
