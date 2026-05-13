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
%     - VTU data is read from AppendedData encoding="raw" sections.
%     - Data columns in amr(ng).data are patch cells; rows are concatenated
%       in order: meqn, aux, rhs, soln, error (when present).

if nargin < 3
    dir = '';
end

if ~isempty(dir)
    lastch = dir(end);
    if ~(lastch == '/' || lastch == '\\')
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
fid_cleanup = onCleanup(@() fclose(h.fid)); %#ok<NASGU>

have_field_data = isfield(h.VTKFile.UnstructuredGrid, 'FieldData');

% Set time
if have_field_data
    t = find_by_name(h.VTKFile.UnstructuredGrid.FieldData.DataArray, 'TimeValue').data;
else
    t = Frame;
end

num_points = h.VTKFile.UnstructuredGrid.Piece.NumberOfPoints;
num_cells  = h.VTKFile.UnstructuredGrid.Piece.NumberOfCells;
arrays = collect_appended_arrays(h.VTKFile);

% Phase 1: read only topology arrays (geometry + metadata, not field data).
% Include explicit metadata arrays written by newer versions of the VTU writer.
topo_names = {'Position', 'connectivity', 'types', 'mpirank', 'blockno', 'patchno', ...
              'patch_dimension', 'levels', 'patch_starts', 'patch_spacings', 'TimeValue'};
topo = decode_specific_arrays(h.fid, h.payload_start, arrays, topo_names);

required_fields = {'Position','connectivity','types','mpirank','blockno','patchno'};
for ireq = 1:numel(required_fields)
    if ~isfield(topo, required_fields{ireq})
        error('Required VTU DataArray "%s" is missing in %s.', required_fields{ireq}, filename);
    end
end

% Convert VTK connectivity to 1-based indexing for MATLAB array access.
points = reshape(topo.Position, 3, []).';
connectivity = double(topo.connectivity) + 1;
connectivity = reshape(connectivity, [], num_cells);
types = double(topo.types(:));

if dim == 2
    expected_type = 9;
    verts_per_cell = 4;
else
    expected_type = 12;
    verts_per_cell = 8;
end

if ~all(types == expected_type)
    error('Unexpected VTK cell type in %s.', filename);
end

if size(connectivity,1) ~= verts_per_cell
    error('Connectivity width does not match expected %d vertices per cell.', verts_per_cell);
end

if size(points,1) ~= piece.num_points
    error('Point count mismatch while reading %s.', filename);
end

if size(connectivity,2) ~= piece.num_cells
    error('Cell count mismatch while reading %s.', filename);
end

% All patches have the same number of cells; compute patch map by arithmetic.
num_patches = numel(unique(double(topo.patchno(:)), 'sorted'));
ncells_per_patch = piece.num_cells / num_patches;
if ncells_per_patch ~= floor(ncells_per_patch)
    error('num_cells (%d) is not evenly divisible by num_patches (%d) in %s.', ...
          piece.num_cells, num_patches, filename);
end
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

field_order = {'meqn', 'aux', 'rhs', 'soln', 'error'};

% Extract global patch_dimension metadata (same for all patches)
if isfield(topo, 'patch_dimension')
    mxmymz_global = reshape(double(topo.patch_dimension), 3, []);
    mx_global = mxmymz_global(1, 1);
    my_global = mxmymz_global(2, 1);
    if dim > 2
        mz_global = mxmymz_global(3, 1);
    else
        mz_global = [];
    end
end

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
    cell_start_0 = (ng-1) * ncells_per_patch;  % 0-based offset into flat arrays
    cell_ids = (cell_start_0+1 : cell_start_0+ncells_per_patch);  % 1-based MATLAB indices

    patch_conn = connectivity(:, cell_ids);
    patch_pts_ids = unique(patch_conn(:));
    patch_pts = points(patch_pts_ids, :);

    % Patch shape (mx, my, mz): use global metadata (same for all patches) or fall back to inference.
    if isfield(topo, 'patch_dimension')
        mx = mx_global;
        my = my_global;
        mz = mz_global;
    else
        if dim == 2
            [mx,my] = infer_2d_shape(patch_conn, points);
            mz = [];
        else
            [mx,my,mz] = infer_3d_shape(patch_conn, points);
        end
    end

    % Patch origin and spacing: prefer explicit FieldData; fall back to inference.
    if isfield(topo, 'patch_starts')
        xyz_low_arr = reshape(double(topo.patch_starts), 3, []);
        xlow = xyz_low_arr(1, ng);
        ylow = xyz_low_arr(2, ng);
        if dim > 2
            zlow = xyz_low_arr(3, ng);
        else
            zlow = [];
        end
    else
        [xlow,~] = infer_axis_spacing(patch_pts(:,1));
        [ylow,~] = infer_axis_spacing(patch_pts(:,2));
        if dim > 2
            [zlow,~] = infer_axis_spacing(patch_pts(:,3));
        end
    end
    
    if isfield(topo, 'patch_spacings')
        dxdydz_arr = reshape(double(topo.patch_spacings), 3, []);
        dx = dxdydz_arr(1, ng);
        dy = dxdydz_arr(2, ng);
        if dim > 2
            dz = dxdydz_arr(3, ng);
        else
            dz = [];
        end
    else
        [~,dx] = infer_axis_spacing(patch_pts(:,1));
        [~,dy] = infer_axis_spacing(patch_pts(:,2));
        if dim > 2
            [~,dz] = infer_axis_spacing(patch_pts(:,3));
        else
            dz = [];
        end
    end

    % Fill legacy AMR metadata expected by plotting/post-processing scripts.
    amrdata.gridno = double(topo.patchno(cell_ids(1))) + 1;
    amrdata.blockno = double(topo.blockno(cell_ids(1)));
    amrdata.mpirank = double(topo.mpirank(cell_ids(1)));
    amrdata.mx = mx;
    amrdata.my = my;
    if dim > 2
        amrdata.mz = mz;
    else
        amrdata.mz = [];
    end

    amrdata.xlow = xlow;
    amrdata.ylow = ylow;
    if dim > 2
        amrdata.zlow = zlow;
    else
        amrdata.zlow = [];
    end

    amrdata.dx = dx;
    amrdata.dy = dy;
    if dim > 2
        amrdata.dz = dz;
    else
        amrdata.dz = [];
    end

    % Phase 2: seek to this patch's cell slice in each field array and read it.
    patch_data = [];
    for k = 1:numel(field_order)
        fname_k = field_order{k};
        fa = find_array(arrays, fname_k);
        if ~isempty(fa)
            raw = read_array_cells(h.fid, h.payload_start, fa, cell_ids(1)-1, numel(cell_ids));
            field_data = reshape(float_array(cast_appended(raw, fa.type, fa.num_components)), [], numel(cell_ids));
            patch_data = [patch_data; field_data]; %#ok<AGROW>
        end
    end
    amrdata.data = patch_data;

    amr(ng) = amrdata; %#ok<AGROW>
end

% If explicit level data was written by the VTU writer, use it directly;
% otherwise infer AMR levels from relative cell spacing (legacy files).
if ~isfield(topo, 'levels')
    amr = assign_levels_from_spacing(amr, dim);
else
    for ng = 1:numel(amr)
        amr(ng).level = double(topo.levels(ng));
    end
end

if dim == 2
    for ng = 1:numel(amr)
        if ~isfield(amr(ng), 'mz')
            amr(ng).mz   = [];
        end
        if ~isfield(amr(ng), 'zlow')
            amr(ng).zlow = [];
        end
        if ~isfield(amr(ng), 'dz')
            amr(ng).dz   = [];
        end
    end
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

% Read only the named DataArrays from the appended section; skip all others.
function values = decode_specific_arrays(fid, payload_start, arrays, names)
values = struct();
for i = 1:numel(arrays)
    if any(strcmp(arrays(i).Name, names))
        a = arrays(i);
        if fseek(fid, payload_start + a.offset, 'bof') ~= 0
            error('Seek failed for DataArray "%s".', a.Name);
        end
        nbytes = double(read_uint64_le_from_fid(fid));
        raw = fread(fid, nbytes, '*uint8');
        if numel(raw) < nbytes
            error('Unexpected end of file reading DataArray "%s".', a.Name);
        end
        values.(a.Name) = cast_appended(raw, a.type, a.num_components);
    end
end
end

% Return the arrays entry whose Name matches, or [] if not found.
function a = find_array(arrays, name)
a = [];
for i = 1:numel(arrays)
    if strcmp(arrays(i).Name, name)
        a = arrays(i);
        return;
    end
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

% Read a contiguous range of cells from a flat appended DataArray by seeking.
% cell_start is 0-based; ncells is the count to read.
function raw = read_array_cells(fid, payload_start, a, cell_start, ncells)
elem_size = vtk_type_size(a.type);
bytes_per_cell = a.num_components * elem_size;
% Skip 8-byte length prefix then jump to cell_start within the data.
data_offset = payload_start + a.offset + 8 + cell_start * bytes_per_cell;
if fseek(fid, data_offset, 'bof') ~= 0
    error('Seek failed reading cells from DataArray "%s".', a.Name);
end
nbytes_to_read = ncells * bytes_per_cell;
raw = fread(fid, nbytes_to_read, '*uint8');
if numel(raw) < nbytes_to_read
    error('Unexpected end of file reading cells from DataArray "%s".', a.Name);
end
end

% Seek to each DataArray by offset and read only its bytes from the open file.
function values = decode_appended_arrays_streaming(fid, payload_start, arrays)
values = struct();
for i = 1:numel(arrays)
    a = arrays(i);
    % Offsets in the VTU header are measured from the byte immediately after '_'.
    if fseek(fid, payload_start + a.offset, 'bof') ~= 0
        error('Seek failed for DataArray "%s".', a.Name);
    end
    nbytes = double(read_uint64_le_from_fid(fid));
    raw = fread(fid, nbytes, '*uint8');
    if numel(raw) < nbytes
        error('Unexpected end of file reading DataArray "%s".', a.Name);
    end
    values.(a.Name) = cast_appended(raw, a.type, a.num_components);
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

% Build connected components of cells that share vertices.
function labels = connected_cell_components(connectivity, npoints)
ncells = size(connectivity, 2);
parent = 1:ncells;
owner = zeros(npoints, 1);

for c = 1:ncells
    verts = connectivity(:, c);
    for iv = 1:numel(verts)
        v = verts(iv);
        if owner(v) == 0
            owner(v) = c;
        else
            % Union-Find merge when two cells touch the same point.
            parent = unite(parent, c, owner(v));
        end
    end
end

labels = zeros(ncells,1);
for c = 1:ncells
    labels(c) = find_root(parent, c);
end
end

% Union operation for cell component labels.
function parent = unite(parent, a, b)
ra = find_root(parent, a);
rb = find_root(parent, b);
if ra ~= rb
    parent(rb) = ra;
end
end

% Find operation for cell component labels.
function r = find_root(parent, x)
r = x;
while parent(r) ~= r
    r = parent(r);
end
end

% Infer (mx,my) assuming patch points lie on a Cartesian x/y grid.
function [mx,my] = infer_2d_shape(patch_conn, points)
patch_point_ids = unique(patch_conn(:));
patch_xy = double(points(patch_point_ids, 1:2));

nx = numel(unique(patch_xy(:,1)));
ny = numel(unique(patch_xy(:,2)));

if nx * ny ~= numel(patch_point_ids)
    error('2D patch points do not form a rectangular Cartesian lattice.');
end

ncells = size(patch_conn, 2);
mx = nx - 1;
my = ny - 1;
if mx * my ~= ncells
    error('Cartesian-grid inference mismatch: mx*my does not equal cell count.');
end
end

% Infer (mx,my,mz) assuming patch points lie on a Cartesian x/y/z grid.
function [mx,my,mz] = infer_3d_shape(patch_conn, points)
patch_point_ids = unique(patch_conn(:));
patch_xyz = double(points(patch_point_ids, 1:3));

nx = numel(unique(patch_xyz(:,1)));
ny = numel(unique(patch_xyz(:,2)));
nz = numel(unique(patch_xyz(:,3)));

if nx * ny * nz ~= numel(patch_point_ids)
    error('3D patch points do not form a rectangular Cartesian lattice.');
end

mx = nx - 1;
my = ny - 1;
mz = nz - 1;

ncells = size(patch_conn, 2);
if mx * my * mz ~= ncells
    error('Cartesian-grid inference mismatch: mx*my*mz does not equal cell count.');
end
end

% Estimate axis origin and spacing from unique coordinate values.
function [x0,dx] = infer_axis_spacing(vals)
u = unique(sort(double(vals(:))));
if numel(u) < 2
    x0 = u(1);
    dx = 0;
    return
end
d = diff(u);
tol = max(1e-12, max(abs(u)) * 1e-10);
d = d(d > tol);
if isempty(d)
    x0 = u(1);
    dx = 0;
else
    x0 = u(1);
    dx = min(d);
end
end

% Convert numeric arrays to double and normalize orientation for reshape.
function arr = float_array(v)
arr = double(v);
if isvector(arr)
    arr = arr(:).';
end
end

% Assign AMR levels from relative spacing tiers.
function amr = assign_levels_from_spacing(amr, dim)
if isempty(amr)
    return
end

all_dx = zeros(numel(amr),1);
for i = 1:numel(amr)
    all_dx(i) = amr(i).dx;
end

levels = spacing_to_levels(all_dx);
for i = 1:numel(amr)
    amr(i).level = levels(i);
    if dim == 2
        if ~isfield(amr(i), 'mz')
            amr(i).mz = [];
        end
        if ~isfield(amr(i), 'zlow')
            amr(i).zlow = [];
        end
        if ~isfield(amr(i), 'dz')
            amr(i).dz = [];
        end
    end
end
end

% Convert unique spacing values to level indices (coarsest -> level 1).
function levels = spacing_to_levels(dx)
udx = unique(sort(dx, 'descend'));
levels = zeros(size(dx));
for i = 1:numel(dx)
    [~,idx] = min(abs(udx - dx(i)));
    levels(i) = idx;
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

% Collect all DataArrays with format="appended" from a VTKFile struct.
% Returns a struct array with fields Name, type, offset, num_components
% compatible with decode_specific_arrays and find_array.
function arrays = collect_appended_arrays(vtk)
arrays = struct('Name', {}, 'type', {}, 'offset', {}, 'num_components', {});
section_names = {'Points', 'Cells', 'CellData', 'PointData'};
piece = vtk.UnstructuredGrid.Piece;
for k = 1:numel(section_names)
    sname = section_names{k};
    if ~isfield(piece, sname) || ~isfield(piece.(sname), 'DataArray')
        continue;
    end
    for i = 1:numel(piece.(sname).DataArray)
        da = piece.(sname).DataArray(i);
        if ~strcmp(da.format, 'appended')
            continue;
        end
        a.Name           = da.Name;
        a.type           = da.type;
        a.offset         = da.offset;
        a.num_components = da.NumberOfComponents;
        arrays(end+1) = a; %#ok<AGROW>
    end
end

% FieldData lives at UnstructuredGrid level, not inside Piece.
if isfield(vtk.UnstructuredGrid, 'FieldData') && ...
        isfield(vtk.UnstructuredGrid.FieldData, 'DataArray')
    for i = 1:numel(vtk.UnstructuredGrid.FieldData.DataArray)
        da = vtk.UnstructuredGrid.FieldData.DataArray(i);
        if ~strcmp(da.format, 'appended')
            continue;
        end
        a.Name           = da.Name;
        a.type           = da.type;
        a.offset         = da.offset;
        a.num_components = da.NumberOfComponents;
        arrays(end+1) = a; %#ok<AGROW>
    end
end
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



