function [amr,t] = readamrdata_forestclaw_vtu(dim,Frame,dir)

% Read ForestClaw VTU output and return AMR patch data.
%
%   [amr,t] = readamrdata_forestclaw_vtu(dim,Frame,dir)
%
%   Input:
%     dim   - Patch dimension (2 or 3).
%     Frame - Output frame number.
%     dir   - Directory containing fort_frame_####.vtu (optional).
%
%   Output:
%     amr   - Struct array with fields compatible with legacy ForestClaw
%             MATLAB readers: gridno, level, blockno, mpirank,
%             mx/my/(mz), xlow/ylow/(zlow), dx/dy/(dz), data.
%     t     - Time from fort.t#### when available; otherwise Frame.
%
%   Notes:
%     - XML header is parsed into a struct.
%     - Local helper read_by_name(h, arrays, name, startIdx, endIdx)
%       accepts an optional 1-based inclusive tuple range.

if nargin < 3
    dir = '';
end

if ~isempty(dir)
    lastch = dir(end);
    if ~(lastch == "/" || lastch == "\\")
        % Normalize directory input so filename concatenation is robust.
        dir = [dir filesep];
    end
end

filename = [dir, sprintf('fort_frame_%04d.vtu', Frame)];
if ~exist(filename, 'file')
    amr = [];
    t = [];
    disp(' ');
    disp(['Frame ',num2str(Frame),' (',filename,') does not exist ***']);
    disp(' ');
    return
end

disp(['Reading data from ',filename]);

h = parse_vtu_header(filename);
fid_cleanup = onCleanup(@() fclose(h.fid));

have_field_data = isfield(h.VTKFile.UnstructuredGrid, 'FieldData');
if ~have_field_data
    error('FieldData section not found in VTU file %s; unable to read time and patch metadata. Possible old Forestclaw VTU output.', filename);
end

% Global Variables
num_cells  = h.VTKFile.UnstructuredGrid.Piece.NumberOfCells;

t = find_by_name(h.VTKFile.UnstructuredGrid.FieldData.DataArray, 'TimeValue').data;
patch_dimension = find_by_name(h.VTKFile.UnstructuredGrid.FieldData.DataArray, 'patch_dimension').data;

patch_starts = read_by_name(h, h.VTKFile.UnstructuredGrid.FieldData.DataArray, 'patch_starts');
patch_spacings = read_by_name(h, h.VTKFile.UnstructuredGrid.FieldData.DataArray, 'patch_spacings');
levels = read_by_name(h, h.VTKFile.UnstructuredGrid.FieldData.DataArray, 'levels');

num_patches = numel(levels);
mx = patch_dimension(1);
my = patch_dimension(2);
if dim == 3 
    mz = patch_dimension(3);
end

% All patches have the same number of cells; compute patch map by arithmetic.
ncells_per_patch = num_cells / num_patches;
amr = struct('gridno', {}, ...
             'level', {}, ...
             'blockno', {}, ...
             'mpirank', {}, ...
             'mx', {}, ...
             'my', {}, ...
             'mz', {}, ...
             'xlow', {}, ...
             'ylow', {}, ...
             'zlow', {}, ...
             'dx', {}, ...
             'dy', {}, ...
             'dz', {}, ...
             'data', {});

for ng = 1:num_patches
    amrdata = struct('gridno', [], ...
                     'level', [], ...
                     'blockno', [], ...
                     'mpirank', [], ...
                     'mx', [], ...
                     'my', [], ...
                     'mz', [], ...
                     'xlow', [], ...
                     'ylow', [], ...
                     'zlow', [], ...
                     'dx', [], ...
                     'dy', [], ...
                     'dz', [], ...
                     'data', []);
    % Compute this patch's 1-based cell index range directly from uniform size.
    cell_start= (ng-1) * ncells_per_patch + 1;  % 0-based offset into flat arrays

    % Fill AMR metadata
    amrdata.gridno = read_by_name(h, h.VTKFile.UnstructuredGrid.Piece.CellData.DataArray, 'patchno', cell_start, cell_start);
    amrdata.level = double(levels(ng));
    amrdata.blockno = read_by_name(h, h.VTKFile.UnstructuredGrid.Piece.CellData.DataArray, 'blockno', cell_start, cell_start);
    amrdata.mpirank = read_by_name(h, h.VTKFile.UnstructuredGrid.Piece.CellData.DataArray, 'mpirank', cell_start, cell_start);
    amrdata.mx = mx;
    amrdata.my = my;
    if dim > 2
        amrdata.mz = mz;
    else
        amrdata.mz = [];
    end

    amrdata.xlow = patch_starts(1, ng);
    amrdata.ylow = patch_starts(2, ng);
    if dim > 2
        amrdata.zlow = patch_starts(3, ng);
    else
        amrdata.zlow = [];
    end

    amrdata.dx = patch_spacings(1, ng);
    amrdata.dy = patch_spacings(2, ng);
    if dim > 2
        amrdata.dz = patch_spacings(3, ng);
    else
        amrdata.dz = [];
    end

    % read meqn
    meqn = read_by_name(h, h.VTKFile.UnstructuredGrid.Piece.CellData.DataArray, 'meqn', cell_start, cell_start + ncells_per_patch - 1);
    amrdata.data = meqn;

    amr(ng) = amrdata;
end

end

% Return the DataArray struct whose Name matches, or [] if not found.
function da = find_by_name(arrays, name)
da = [];
for i = 1:numel(arrays)
    if strcmp(arrays(i).Name, name)
        da = arrays(i);
        return;
    end
end
end

% Read a named DataArray by Name from either FieldData or Piece sections.
% Accepts raw XML DataArray structs or normalized entries from collect_appended_arrays.
% Optional startIdx/endIdx select an inclusive 1-based tuple range to minimize disk I/O.
function values = read_by_name(h, arrays, name, startIdx, endIdx)
values = [];
da = find_by_name(arrays, name);
if isempty(da)
    return;
end

ncomp = data_array_num_components(da);

have_range = (nargin >= 4 && ~isempty(startIdx)) || (nargin >= 5 && ~isempty(endIdx));
if nargin < 4 || isempty(startIdx)
    startIdx = 1;
end
if nargin < 5 || isempty(endIdx)
    endIdx = [];
end

da_format = data_array_format(da);
if strcmp(da_format, 'appended')
    if have_range
        % Read only the requested tuple range to minimize disk I/O
        if isempty(endIdx)
            % Need total tuple count; read full length prefix
            if fseek(h.fid, h.payload_start + da.offset, 'bof') ~= 0
                error('Seek failed for DataArray "%s".', da.Name);
            end
            nbytes = double(read_uint64_le_from_fid(h.fid));
            elem_size = vtk_type_size(da.type);
            ntuples = nbytes / (elem_size * ncomp);
            endIdx = ntuples;
        end
        % Validate range
        validate_range_indices(da.Name, startIdx, endIdx);
        % Seek to first tuple and read only requested range
        raw = read_data_array_range(h, da, startIdx, endIdx);
    else
        % Read full array
        if fseek(h.fid, h.payload_start + da.offset, 'bof') ~= 0
            error('Seek failed for DataArray "%s".', da.Name);
        end
        nbytes = double(read_uint64_le_from_fid(h.fid));
        raw = fread(h.fid, nbytes, '*uint8');
        if numel(raw) < nbytes
            error('Unexpected end of file reading DataArray "%s".', da.Name);
        end
    end
    values = cast_appended(raw, da.type, ncomp);
elseif strcmp(da_format, 'ascii')
    values = da.data;
    if have_range
        if isempty(endIdx)
            if ncomp > 1 && ismatrix(values) && size(values,1) == ncomp
                endIdx = size(values,2);
            else
                endIdx = numel(values);
            end
        end
        validate_range_indices(da.Name, startIdx, endIdx);
        if ncomp > 1 && ismatrix(values) && size(values,1) == ncomp
            values = values(:, startIdx:endIdx);
        else
            values = values(startIdx:endIdx);
        end
    end
end
end

% Read only a range of tuples from an appended DataArray by seeking.
% startIdx/endIdx are 1-based inclusive tuple indices.
function raw = read_data_array_range(h, da, startIdx, endIdx)
elem_size = vtk_type_size(da.type);
bytes_per_tuple = elem_size * data_array_num_components(da);
% Seek past 8-byte length prefix, then to startIdx-th tuple (1-based to 0-based).
data_offset = h.payload_start + da.offset + 8 + (startIdx - 1) * bytes_per_tuple;
if fseek(h.fid, data_offset, 'bof') ~= 0
    error('Seek failed reading range from DataArray "%s".', da.Name);
end
nbytes_to_read = (endIdx - startIdx + 1) * bytes_per_tuple;
raw = fread(h.fid, nbytes_to_read, '*uint8');
if numel(raw) < nbytes_to_read
    error('Unexpected end of file reading range from DataArray "%s".', da.Name);
end
end

% Return DataArray tuple width with compatibility for alternate field names.
function ncomp = data_array_num_components(da)
if isfield(da, 'NumberOfComponents') && ~isempty(da.NumberOfComponents)
    ncomp = da.NumberOfComponents;
else
    ncomp = 1;
end
end

% Return DataArray format; normalized appended arrays may omit this field.
function fmt = data_array_format(da)
if isfield(da, 'format') && ~isempty(da.format)
    fmt = da.format;
else
    fmt = '';
end
end

% Validate a requested inclusive tuple range for a named DataArray.
function validate_range_indices(name, startIdx, endIdx)
if ~isscalar(startIdx) || ~isnumeric(startIdx) || ~isfinite(startIdx) || startIdx ~= floor(startIdx)
    error('Invalid startIdx for DataArray "%s": expected a finite integer scalar.', name);
end
if ~isscalar(endIdx) || ~isnumeric(endIdx) || ~isfinite(endIdx) || endIdx ~= floor(endIdx)
    error('Invalid endIdx for DataArray "%s": expected a finite integer scalar.', name);
end
if startIdx < 1
    error('Invalid range for DataArray "%s": startIdx=%d must be >= 1.', name, startIdx);
end
if endIdx < startIdx
    error('Invalid range for DataArray "%s": startIdx=%d must be <= endIdx=%d.', name, startIdx, endIdx);
end
end

% Return the byte width of one scalar element for a VTK type string.
function n = vtk_type_size(vtk_type)
switch vtk_type
    case {'Float64', 'Int64'}
        n = 8;
    case {'Float32', 'Int32'}
        n = 4;
    case 'UInt8'
        n = 1;
    otherwise
        error('Unknown VTK type "%s".', vtk_type);
end
end

% Read a little-endian UInt64 length prefix from the current file position.
function u = read_uint64_le_from_fid(fid)
raw = fread(fid, 8, '*uint8');
if numel(raw) < 8
    error('Unexpected end of file reading 8-byte length prefix.');
end
u = uint64(0);
for k = 0:7
    u = bitor(u, bitshift(uint64(raw(k+1)), 8*k));
end
end

% Convert raw byte payload into MATLAB numeric arrays based on VTK type.
function out = cast_appended(raw, vtk_type, ncomp)
switch vtk_type
    case 'Float64'
        out = typecast(raw, 'double');
    case 'Float32'
        out = typecast(raw, 'single');
    case 'Int32'
        out = typecast(raw, 'int32');
    case 'Int64'
        out = typecast(raw, 'int64');
    case 'UInt8'
        out = uint8(raw);
    otherwise
        error('Unsupported VTK type: %s', vtk_type);
end

if ncomp > 1
    out = reshape(out, ncomp, []);
end
end

% ============================ XML Header Parsing ============================

% Parse a VTU file header into a struct h with fields:
%   h.fid           - open file handle (caller must close)
%   h.payload_start - 0-based file offset of the first payload byte
%   h.VTKFile       - struct mirroring the XML header hierarchy:
%       .type, .version, .byte_order, .header_type
%       .UnstructuredGrid.Piece.NumberOfPoints  (numeric)
%       .UnstructuredGrid.Piece.NumberOfCells   (numeric)
%       .UnstructuredGrid.Piece.Points.DataArray(...)
%       .UnstructuredGrid.Piece.Cells.DataArray(...)
%       .UnstructuredGrid.Piece.CellData.DataArray(...)
%       .UnstructuredGrid.FieldData.DataArray(...)  [if present]
%       .AppendedData.encoding
function h = parse_vtu_header(filename)
fid = fopen(filename, 'r');
if fid < 0
    error('Unable to open %s', filename);
end

chunk_size = 65536;
buf = uint8([]);
needle = uint8('<AppendedData');

while true
    chunk = fread(fid, chunk_size, '*uint8');
    buf = [buf; chunk]; %#ok<AGROW>

    ad_idx = strfind(buf.', needle);
    if ~isempty(ad_idx)
        search_from = ad_idx(1) + numel(needle);
        u_rel = find(buf(search_from:end) == uint8('_'), 1, 'first');
        if ~isempty(u_rel)
            % underscore_pos is the 1-based index of '_' in buf.
            underscore_pos = search_from + u_rel - 1;
            % h.payload_start is the 0-based file offset of the byte after '_'.
            h.fid = fid;
            h.payload_start = underscore_pos;  % equals 0-based offset because MATLAB is 1-based
            h.VTKFile = build_vtu_struct(char(buf(1:underscore_pos-1)).');
            return;
        end
    end

    if isempty(chunk)
        fclose(fid);
        error('No AppendedData payload marker "_" found in %s.', filename);
    end
end
end

% Build a MATLAB struct mirroring the VTU XML header hierarchy.
function vtk = build_vtu_struct(header_text)
vtk = tag_attrs(header_text, 'VTKFile');

piece = tag_attrs(header_text, 'Piece');
piece.NumberOfPoints = str2double(piece.NumberOfPoints);
piece.NumberOfCells  = str2double(piece.NumberOfCells);

section_names = {'Points', 'Cells', 'CellData', 'PointData'};
for k = 1:numel(section_names)
    sname = section_names{k};
    content = extract_between_tags(header_text, sname);
    if ~isempty(content)
        sect = tag_attrs(header_text, sname);
        sect.DataArray = section_data_arrays(content);
        piece.(sname) = sect;
    end
end

vtk.UnstructuredGrid.Piece = piece;

fd_content = extract_between_tags(header_text, 'FieldData');
if ~isempty(fd_content)
    fd = tag_attrs(header_text, 'FieldData');
    fd.DataArray = section_data_arrays(fd_content);
    % Read inline values for ASCII DataArrays immediately during header parse.
    for i = 1:numel(fd.DataArray)
        if strcmp(fd.DataArray(i).format, 'ascii')
            fd.DataArray(i).data = read_ascii_da_content(fd_content, fd.DataArray(i).Name);
        end
    end
    vtk.UnstructuredGrid.FieldData = fd;
end

vtk.AppendedData = tag_attrs(header_text, 'AppendedData');
end

% Extract all XML attributes from the first <tagname ...> tag as a struct.
function s = tag_attrs(text, tagname)
pat = ['<', tagname, '(?=[\s>\/])(\s[^>]*)?>'];
tok = regexp(text, pat, 'tokens', 'once');
if isempty(tok) || isempty(tok{1})
    s = struct();
    return;
end
kv = regexp(tok{1}, '(\w+)="([^"]*)"', 'tokens');
s = struct();
for i = 1:numel(kv)
    s.(kv{i}{1}) = kv{i}{2};
end
end

% Return the content between <tagname ...> and </tagname>, or '' if absent.
function content = extract_between_tags(text, tagname)
open_pat = ['<', tagname, '(?=[\s>])([^>]*)>'];
[~, e] = regexp(text, open_pat, 'start', 'end', 'once');
if isempty(e)
    content = '';
    return;
end
close_pat = ['</', tagname, '>'];
c_rel = regexp(text(e+1:end), close_pat, 'start', 'once');
if isempty(c_rel)
    content = '';
    return;
end
content = text(e+1 : e + c_rel - 1);
end

% Parse all <DataArray ...> elements in text into a struct array.
% Fields: Name, type, NumberOfComponents (numeric), format, offset (numeric), data ([] or numeric).
function das = section_data_arrays(text)
das = struct('Name', {}, 'type', {}, 'NumberOfComponents', {}, 'format', {}, 'offset', {}, 'data', {});
toks = regexp(text, '<DataArray\s+([^>]*)>', 'tokens');
for i = 1:numel(toks)
    attr_str = toks{i}{1};
    da.Name               = read_attr(attr_str, 'Name');
    da.type               = read_attr(attr_str, 'type');
    da.NumberOfComponents = str2double_or(read_attr(attr_str, 'NumberOfComponents'), 1);
    da.format             = read_attr(attr_str, 'format');
    da.offset             = str2double_or(read_attr(attr_str, 'offset'), NaN);
    da.data               = [];
    das(end+1) = da; %#ok<AGROW>
end
end

% Extract numeric values from an inline ASCII DataArray element, matched by Name.
function vals = read_ascii_da_content(text, name)
pat = ['<DataArray\s+[^>]*Name="', name, '"[^>]*>([\s\S]*?)<\/DataArray>'];
tok = regexp(text, pat, 'tokens', 'once');
if isempty(tok) || isempty(tok{1})
    vals = [];
    return;
end
vals = sscanf(strtrim(tok{1}), '%g');
end

% Return str2double(s) when s is non-empty, otherwise return default.
function v = str2double_or(s, default)
if isempty(s)
    v = default;
else
    v = str2double(s);
end
end

% Read one XML attribute value from an attribute string.
function value = read_attr(attrs, key)
pat = [key, '="([^"]+)"'];
tok = regexp(attrs, pat, 'tokens', 'once');
if isempty(tok)
    value = '';
else
    value = tok{1};
end
end



