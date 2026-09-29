function latlon = get_latlon(dirD,bas)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GET_LATLON Return basin latitude and longitude coordinates.
%
% SYNOPSIS:
%   latlon = get_latlon(dirD,bas)
%   dirD, bas regional data directory and prepared basin structure. If
%             bas.zone contains nonempty lat/lon arrays, those values
%             are returned without reading a metadata file.
%   bas.id_gauge  basin or gauge identifiers used for metadata matching
%
% OUTPUT:
%   latlon    Kx2 matrix with latitude and longitude of requested basins
%             column 1: latitude
%             column 2: longitude
%
% DESCRIPTION:
%   This helper function is used by exported stand-alone SAGE run scripts
%   to recover basin coordinates from gauge_information.txt. The function
%   reads the gauge metadata file in dirD, normalizes the requested basin
%   identifiers and metadata identifiers, matches the requested basins, and
%   returns their latitude and longitude coordinates.
%
%   For CAMELS-FR gauge_information.txt, the expected column layout is:
%      column 1: gauge identifier
%      column 2: gauge name
%      column 3: latitude
%      column 4: longitude
%      column 5: drainage area
%
%   The returned coordinates are used for exported diagnostics and regional
%   plotting outside the graphical SAGE user interface.
%
% EXAMPLE:
%   latlon = get_latlon(dirD,bas);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written by Jasper A. Vrugt, Feb. 2026                                 %
% University of California, Irvine                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isfield(bas,'zone') && isfield(bas.zone,'lat') ...
            && isfield(bas.zone,'lon') ...
            && ~isempty(bas.zone.lat) && ~isempty(bas.zone.lon)
        latlon = [double(bas.zone.lat(:)) ...
            double(bas.zone.lon(:))];
        return
    end

    if isempty(dirD)
        error('get_latlon: dirD is missing or empty.');
    end
    if ~isstruct(bas) ...
            || ~isfield(bas,'id_gauge') ...
            || isempty(bas.id_gauge)
        error('get_latlon: bas.id_gauge is missing or empty.');
    end
    dirD = char(string(dirD));
    id_USGS = bas.id_gauge;

    f1 = fullfile(dirD,'gauge_information.txt');
    f2 = fullfile(dirD,'gauge_information_COL.CSV');
    
    if isfile(f1)
        fname_gauge = f1;
    elseif isfile(f2)
        fname_gauge = f2;
    else
        error(['Cannot find gauge_information ' ...
            'file in: %s'],dirD);
    end
    
    opts = detectImportOptions(fname_gauge, ...
        'VariableNamingRule','preserve');
    
    T = readtable(fname_gauge,opts);
    
    id_req = normalize_usgs(id_USGS);

    names0 = string(T.Properties.VariableNames);
    names = lower(names0);
    names = regexprep(names, ...
        '[^a-z0-9]+','_');
    names = regexprep(names, ...
        '^x_','');
    names = regexprep(names, ...
        '_$','');

    idCol = find(ismember(names, ...
        ["gauge_id","gage_id"]),1);

    latCol = find(ismember(names, ...
        ["lat","latitude", ...
        "gauge_lat"]),1);

    lonCol = find(ismember(names, ...
        ["lon","long","longitude", ...
        "gauge_lon"]),1);

    if isempty(idCol) ...
            || isempty(latCol) ...
            || isempty(lonCol)

        try
            id1 = normalize_usgs(string(T{:,1}));
            [tf1,~] = ismember(id_req,id1);
        catch
            tf1 = false(size(id_req));
        end

        try
            id2 = normalize_usgs(string(T{:,2}));
            [tf2,~] = ismember(id_req,id2);
        catch
            tf2 = false(size(id_req));
        end

        if all(tf1) ...
                && width(T) >= 4
            idCol = 1;
            latCol = 3;
            lonCol = 4;
        elseif all(tf2) ...
                && width(T) >= 5
            idCol = 2;
            latCol = 4;
            lonCol = 5;
        else
            disp(T.Properties.VariableNames)
            error(['Could not identify gauge ID, ' ...
                'latitude, or longitude columns.']);
        end
    end

    id_meta = normalize_usgs(string(T{:,idCol}));
    [tf,loc] = ismember(id_req,id_meta);

    if ~all(tf)
        miss = id_req(~tf);
        error(['Missing gauge metadata ' ...
            'for basin %s.'],char(miss(1)));
    end

    latlon = double(T{loc,[latCol lonCol]});
end
