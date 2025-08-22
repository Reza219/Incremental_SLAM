function graph = read_graph_g2o(filename)
% READ_GRAPH_G2O  Parse a 2D SLAM g2o file (SE2 poses, XY landmarks).
%
% SYNTAX
%   graph = read_graph_g2o(filename)
%
% DESCRIPTION
%   Supports the following g2o records:
%     - VERTEX_SE2     id  x  y  theta
%     - VERTEX_XY      id  x  y
%     - EDGE_SE2       id_i id_j  dx dy dth   Ixx Ixy Ixt  Iyy Iyt  Itt
%     - EDGE_SE2_XY    id_i id_l  zx zy       Ixx Ixy Iyy
%
% OUTPUT
%   graph : struct with fields
%       .x         (nx1 double)           stacked state [poses; landmarks]
%       .edges     (struct array)         fields: type∈{'P','L'}, from,to,
%                                         measurement, information, fromIdx,toIdx
%       .idLookup  (struct)               id%d -> .offset (0‑based), .dimension
%       .var2node  (nx1 double)           map: variable index -> node id
%
% NOTES
%   - Offsets in .idLookup are 0‑based (for consistency with your code).
%   - Indices .fromIdx/.toIdx are 1‑based (MATLAB/Octave array indexing).
%   - Information matrices are symmetrized defensively.

% ---- Open file -----------------------------------------------------------
fid = fopen(filename, 'r');
if fid == -1, error('read_graph_g2o:OpenFailed','Cannot open file: %s', filename); end
cleanupObj = onCleanup(@() fclose(fid));

graph = struct('x', [], 'edges', [], 'idLookup', struct(), 'var2node', []);
edges = [];  % will be grown as a struct array

% ---- Parse line by line --------------------------------------------------
while true
  ln = fgetl(fid);
  if ~ischar(ln), break; end
  % Strip comments (# ...) and trim
  hash = find(ln=='#', 1, 'first');
  if ~isempty(hash), ln = ln(1:hash-1); end
  ln = strtrim(ln);
  if isempty(ln), continue; end

  tokens = strsplit(ln, ' ', 'CollapseDelimiters', true);
  tag = tokens{1};
  % Convert numeric tail once (NaNs will be ignored by tag switch)
  vals = str2double(tokens(2:end));

  switch tag
    case 'VERTEX_SE2'
      if numel(vals) < 4, continue; end
      id = vals(1); xytheta = vals(2:4).';
      key = sprintf('id%d', int32(id));
      graph.idLookup.(key) = struct('offset', numel(graph.x), 'dimension', 3);
      graph.x = [graph.x; xytheta];

    case 'VERTEX_XY'
      if numel(vals) < 3, continue; end
      id = vals(1); xy = vals(2:3).';
      key = sprintf('id%d', int32(id));
      graph.idLookup.(key) = struct('offset', numel(graph.x), 'dimension', 2);
      graph.x = [graph.x; xy];

    case 'EDGE_SE2'
      % from to  dx dy dth  Ixx Ixy Ixt  Iyy Iyt  Itt
      if numel(vals) < 11, continue; end
      fromId = int32(vals(1)); toId = int32(vals(2));
      meas   = vals(3:5).';
      ut     = vals(6:11).';  % upper-triangular order
      Info   = [ut(1) ut(2) ut(3);
                ut(2) ut(4) ut(5);
                ut(3) ut(5) ut(6)];
      Info   = (Info + Info.')/2;  % symmetrize defensively
      edges = append_edge(edges, 'P', fromId, toId, meas, Info);

    case 'EDGE_SE2_XY'
      % from pose id, to landmark id; z in robot frame
      if numel(vals) < 5, continue; end
      fromId = int32(vals(1)); toId = int32(vals(2));
      meas   = vals(3:4).';
      ut     = vals(5:7).';
      Info   = [ut(1) ut(2);
                ut(2) ut(3)];
      Info   = (Info + Info.')/2;
      edges = append_edge(edges, 'L', fromId, toId, meas, Info);

    otherwise
      % ignore other tags
  end
end

graph.edges = edges;

% ---- Prepare helper structures ------------------------------------------
% (1) Fill fromIdx/toIdx for edges (1-based indices into x)
for k = 1:numel(graph.edges)
  e = graph.edges(k);
  fk = sprintf('id%d', e.from);
  tk = sprintf('id%d', e.to);
  assert(isfield(graph.idLookup,fk),'read_graph_g2o:MissingFrom','Missing vertex id%d', e.from);
  assert(isfield(graph.idLookup,tk),'read_graph_g2o:MissingTo','Missing vertex id%d', e.to);
  graph.edges(k).fromIdx = graph.idLookup.(fk).offset + 1;
  graph.edges(k).toIdx   = graph.idLookup.(tk).offset + 1;
end

% (2) Build var2node (size = numel(x))
n = numel(graph.x);
graph.var2node = zeros(n,1);
id_keys = fieldnames(graph.idLookup);
for i = 1:numel(id_keys)
  v = graph.idLookup.(id_keys{i});
  nid = sscanf(id_keys{i}, 'id%d');
  if isempty(nid), continue; end
  idx = (v.offset + 1):(v.offset + v.dimension);
  graph.var2node(idx) = nid;
end
end

% ==================== Local helper =======================================
function edges = append_edge(edges, type, fromId, toId, measurement, information)
% Append an edge struct; keeps fields consistent.
E = struct('type', type, 'from', fromId, 'to', toId, ...
           'measurement', measurement, 'information', information, ...
           'fromIdx', [], 'toIdx', []);
if isempty(edges)
  edges = E;
else
  edges(end+1) = E; %#ok<AGROW>
end
end
