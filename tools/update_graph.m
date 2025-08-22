function [g_current, edge_nodes, new_nodes, edge_vars, loop_closure] = update_graph(edge_ids, g, g_current, lc_gap)
% UPDATE_GRAPH  Update the incremental graph with new edges (and any new nodes).
%
% SYNTAX
%   [g_current, edge_nodes, new_nodes, edge_vars, loop_closure] = ...
%       update_graph(edge_ids, g, g_current, lc_gap)
%
% DESCRIPTION
%   - Appends edges g.edges(edge_ids) to g_current.edges in a safe, incremental way.
%   - Adds any previously unseen nodes referenced by those edges (copies state from g.x).
%   - Maintains a per-call persistent state for lightweight loop-closure detection.
%   - Treats GNSS ('G') edges as unary: uses .from only and sets toIdx = 0.
%
% INPUTS
%   edge_ids  : vector of global edge indices into g.edges
%   g         : full graph (fields: .x, .edges, .idLookup[id%d].offset/.dimension)
%   g_current : current graph to be extended (same structure as g)
%   lc_gap    : numeric parameter for loop-closure detection (passed through)
%
% OUTPUTS
%   g_current   : updated graph
%   edge_nodes  : unique node ids referenced by the new edges
%   new_nodes   : node ids that were added to g_current in this call
%   edge_vars   : variable indices (in global g.x) belonging to the newly added nodes
%   loop_closure: result of detect_loop_closure_unordered on these edges
%
% NOTES
%   - MATLAB R2020b+ and GNU Octave 7+ compatible (no string class).
%   - Expects a helper: detect_loop_closure_unordered(edges, state, lc_gap).

% ---- Validate ------------------------------------------------------------
assert(isstruct(g) && isfield(g,'edges') && isfield(g,'idLookup'), ...
  'update_graph:InvalidG', 'g must have .edges and .idLookup');
if ~exist('g_current','var') || isempty(g_current)
  g_current = struct('x',zeros(0,1), 'edges', g.edges([]), 'idLookup', struct, 'var2node', zeros(0,1));
else
  if ~isfield(g_current,'x'),        g_current.x        = zeros(0,1); end
  if ~isfield(g_current,'edges'),    g_current.edges    = g.edges([]); end
  if ~isfield(g_current,'idLookup'), g_current.idLookup = struct;      end
  if ~isfield(g_current,'var2node'), g_current.var2node = zeros(0,1);  end
end

edge_ids = edge_ids(:).';
assert(all(edge_ids >= 1 & edge_ids <= numel(g.edges)), 'update_graph:BadEdgeIds');

% ---- Step 1: Identify involved nodes ------------------------------------
edges = g.edges(edge_ids);

% Collect node ids robustly (support unary GNSS edges with only .from)
from_ids = arrayfun(@(e) getfield_safe(e,'from',NaN), edges); %#ok<GFLD>
to_ids   = arrayfun(@(e) getfield_safe(e,'to',NaN),   edges); %#ok<GFLD>
edge_nodes = unique([from_ids(~isnan(from_ids)), to_ids(~isnan(to_ids) & to_ids~=0)]);

% Parse existing node ids in g_current
id_fields = fieldnames(g_current.idLookup);
seen_nodes = zeros(1, numel(id_fields));
for k = 1:numel(id_fields)
  nk = sscanf(id_fields{k}, 'id%d');
  if ~isempty(nk), seen_nodes(k) = nk; else, seen_nodes(k) = NaN; end
end
seen_nodes = seen_nodes(~isnan(seen_nodes));

% New nodes are those referenced by edges but not yet in g_current
new_nodes = setdiff(edge_nodes, seen_nodes);
edge_vars = []; % will be filled after adding nodes

% ---- Step 2: Loop-closure detection (persistent) ------------------------
persistent lc_state
if isempty(lc_state), lc_state = []; end
loop_closure = detect_loop_closure_unordered(edges, lc_state, lc_gap);
lc_state = loop_closure; % let detector carry its own state forward as needed

% ---- Step 3: Add new nodes ----------------------------------------------
if ~isempty(new_nodes)
  % Filter to nodes that exist in g.idLookup (defensive)
  new_keys = arrayfun(@(n) sprintf('id%d', n), new_nodes, 'UniformOutput', false);
  valid = cellfun(@(k) isfield(g.idLookup, k), new_keys);
  new_nodes = new_nodes(valid);
  new_keys  = new_keys(valid);

  if ~isempty(new_nodes)
    dims = arrayfun(@(n) g.idLookup.(sprintf('id%d', n)).dimension, new_nodes);
    total_dim = sum(dims);

    newX        = zeros(total_dim, 1);
    newVar2node = zeros(total_dim, 1);
    newEdgeVars = zeros(total_dim, 1);

    pos = 1;
    for i = 1:numel(new_nodes)
      nid = new_nodes(i);
      key = new_keys{i};
      info = g.idLookup.(key);
      dim  = info.dimension;

      % Copy idLookup entry and assign a new offset (0-based)
      g_current.idLookup.(key) = info;
      g_current.idLookup.(key).offset = numel(g_current.x) + pos - 1;

      % Copy initial state from g.x (global) to g_current.x
      src = (info.offset + 1) : (info.offset + dim);
      dst = pos:(pos+dim-1);

      newX(dst)        = g.x(src);
      newVar2node(dst) = nid;
      newEdgeVars(dst) = src;

      pos = pos + dim;
    end

    % Append to current graph
    g_current.x        = [g_current.x; newX];
    g_current.var2node = [g_current.var2node; newVar2node];
    edge_vars          = unique(newEdgeVars);
  end
end

% ---- Step 4: Set edge index fields (fromIdx/toIdx) -----------------------
n_existing_edges = numel(g_current.edges);
n_new_edges      = numel(edges);
new_edge_indices = (n_existing_edges + 1) : (n_existing_edges + n_new_edges);

for e = 1:n_new_edges
  ei = edges(e);

  % Type (default to non-GNSS if missing)
  etype = 'X';
  if isfield(ei,'type') && ~isempty(ei.type), etype = ei.type; end

  % from index must exist
  assert(isfield(ei,'from'), 'update_graph:MissingFrom', 'Edge missing .from');
  from_key = sprintf('id%d', ei.from);
  assert(isfield(g_current.idLookup, from_key), 'update_graph:UnknownFrom', ...
    'Missing from-node id%d in g_current', ei.from);
  ei.fromIdx = g_current.idLookup.(from_key).offset + 1;

  if etype == 'G'
    % GNSS unary: toIdx = 0
    ei.toIdx = 0;
  else
    % binary edge requires .to
    assert(isfield(ei,'to'), 'update_graph:MissingTo', 'Non-GNSS edge missing .to');
    to_key = sprintf('id%d', ei.to);
    assert(isfield(g_current.idLookup, to_key), 'update_graph:UnknownTo', ...
      'Missing to-node id%d in g_current', ei.to);
    ei.toIdx = g_current.idLookup.(to_key).offset + 1;
  end

  edges(e) = ei; % write back
end

% ---- Step 5: Update per-node edge lists in idLookup ----------------------
for e = 1:n_new_edges
  ei = edges(e);
  % from node
  fk = sprintf('id%d', ei.from);
  if ~isfield(g_current.idLookup.(fk), 'edges') || isempty(g_current.idLookup.(fk).edges)
    g_current.idLookup.(fk).edges = new_edge_indices(e);
  else
    g_current.idLookup.(fk).edges(end+1) = new_edge_indices(e);
  end

  % to node (if binary / non-GNSS)
  if ~(isfield(ei,'type') && ~isempty(ei.type) && ei.type == 'G')
    tk = sprintf('id%d', ei.to);
    if ~isfield(g_current.idLookup.(tk), 'edges') || isempty(g_current.idLookup.(tk).edges)
      g_current.idLookup.(tk).edges = new_edge_indices(e);
    else
      g_current.idLookup.(tk).edges(end+1) = new_edge_indices(e);
    end
  end
end

% ---- Step 6: Append edges to g_current ----------------------------------
if isempty(g_current.edges)
  g_current.edges = edges;
else
  g_current.edges(end+1:end+n_new_edges) = edges;
end

end

% ========== Local utility =================================================
function v = getfield_safe(s, f, default_val)
%GETFIELD_SAFE  Return s.(f) if present, else default_val.
if isstruct(s) && isfield(s,f)
  v = s.(f);
else
  v = default_val;
end
end
