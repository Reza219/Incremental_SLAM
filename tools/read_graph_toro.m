function graph = read_graph_toro(filename)
% READ_GRAPH_TORO  Parse a TORO 2D SLAM file (VERTEX2/EDGE2 only).
%
% SYNTAX
%   graph = read_graph_toro(filename)
%
% DESCRIPTION
%   Supports:
%     - VERTEX2    <id>  x  y  theta
%     - EDGE2      <i> <j>  dx dy dth   Ixx Ixy Iyy Itt Ixth Iy th (TORO-style order)
%       (We read the 6 unique entries; any extra tokens on the line are ignored.)
%
% OUTPUT
%   graph : struct with fields
%       .x         (nx1 double)           stacked state [x1;y1;th1; x2;y2;th2; ...]
%       .edges     (struct array)         fields: type='P', from,to,measurement,information,fromIdx,toIdx
%       .idLookup  (struct)               id%d -> .offset (0‑based), .dimension (=3)
%       .var2node  (nx1 double)           map: variable index -> node id
%
% NOTES
%   - Offsets in .idLookup are 0‑based (matches the rest of your tools).
%   - Indices .fromIdx/.toIdx are 1‑based (MATLAB/Octave arrays).
%   - Information matrices are symmetrized defensively.

% ---- Open file -----------------------------------------------------------
fid = fopen(filename, 'r');
if fid == -1, error('read_graph_toro:OpenFailed','Cannot open file: %s', filename); end
cleanupObj = onCleanup(@() fclose(fid));

graph = struct('x', [], 'edges', [], 'idLookup', struct(), 'var2node', []);

edges = [];  % will grow as struct array

% ---- Parse line by line --------------------------------------------------
while true
  ln = fgetl(fid);
  if ~ischar(ln), break; end

  % Strip comments and trim
  hash = find(ln=='#', 1, 'first');
  if ~isempty(hash), ln = ln(1:hash-1); end
  ln = strtrim(ln);
  if isempty(ln), continue; end

  tok = strsplit(ln, ' ', 'CollapseDelimiters', true);
  tag = tok{1};

  switch tag
    case 'VERTEX2'
      if numel(tok) < 5, continue; end
      id = str2double(tok{2});
      vals = [str2double(tok{3}); str2double(tok{4}); str2double(tok{5})];

      key = sprintf('id%d', int32(id));
      % 0-based offset is current length of x
      graph.idLookup.(key) = struct('offset', numel(graph.x), 'dimension', 3);
      graph.x = [graph.x; vals];

    case 'EDGE2'
      if numel(tok) < 13
        % Some TORO files have at least 12 tokens after tag; skip if too short
        continue;
      end
      fromId = int32(str2double(tok{2}));
      toId   = int32(str2double(tok{3}));
      meas   = [str2double(tok{4}); str2double(tok{5}); str2double(tok{6})];

      % TORO provides 6 unique info entries among tokens 7..12 (we take the common layout)
      % Here we follow the same positions as your original implementation.
      Info = zeros(3,3);
      Info(1,1) = str2double(tok{7});
      Info(1,2) = str2double(tok{8});
      Info(2,1) = Info(1,2);
      Info(2,2) = str2double(tok{9});
      Info(3,3) = str2double(tok{10});
      Info(1,3) = str2double(tok{11});
      Info(3,1) = Info(1,3);
      Info(2,3) = str2double(tok{12});
      Info(3,2) = Info(2,3);
      Info = (Info + Info.')/2; % symmetrize defensively

      edges = append_edge(edges, 'P', fromId, toId, meas, Info);

    otherwise
      % ignore other tags
  end
end

graph.edges = edges;

% ---- Fill fromIdx/toIdx (1-based) ---------------------------------------
for k = 1:numel(graph.edges)
  e = graph.edges(k);
  fk = sprintf('id%d', e.from);
  tk = sprintf('id%d', e.to);
  assert(isfield(graph.idLookup,fk), 'read_graph_toro:MissingFrom', 'Missing vertex id%d', e.from);
  assert(isfield(graph.idLookup,tk), 'read_graph_toro:MissingTo',   'Missing vertex id%d', e.to);
  graph.edges(k).fromIdx = graph.idLookup.(fk).offset + 1;
  graph.edges(k).toIdx   = graph.idLookup.(tk).offset + 1;
end

% ---- Build var2node ------------------------------------------------------
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
E = struct('type', type, 'from', fromId, 'to', toId, ...
           'measurement', measurement, 'information', information, ...
           'fromIdx', [], 'toIdx', []);
if isempty(edges)
  edges = E;
else
  edges(end+1) = E; %#ok<AGROW>
end
end
