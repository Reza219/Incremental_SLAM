function [affected_vars, affected_edge_ids, affected_nodes] = find_affected(g, dxa, dx_th, affected_vars)
% FIND_AFFECTED  Identify variables, nodes, and edges impacted by significant updates.
%
% SYNTAX
%   [affected_vars, affected_edge_ids, affected_nodes] = ...
%       find_affected(g, dxa, dx_th, affected_vars)
%
% DESCRIPTION
%   Thresholds per-variable updates (dxa) and expands to all variables of the
%   corresponding nodes, then returns the set of incident edge indices.
%
% INPUTS
%   g             : graph struct with fields
%                   .var2node (Nx1 double)        map var index -> node id
%                   .idLookup struct: id%d -> .offset (0-based), .dimension, .edges (optional)
%   dxa           : vector of delta-x values corresponding to affected_vars
%   dx_th         : scalar threshold for significance (e.g., 1e-3)
%   affected_vars : candidate variable indices (global indices into g.x)
%
% OUTPUTS
%   affected_vars     : expanded & deduplicated variable indices (all vars of affected nodes)
%   affected_edge_ids : unique incident edge indices touching affected nodes
%   affected_nodes    : node ids whose variables exceeded threshold
%
% NOTES
%   - MATLAB R2020b+ / GNU Octave 7+ compatible.
%   - If a node lacks `.edges`, it contributes no edges (handled gracefully).

% ---- Validate ------------------------------------------------------------
assert(isvector(dxa) && isvector(affected_vars) && numel(dxa)==numel(affected_vars), ...
  'find_affected:SizeMismatch', 'dxa and affected_vars must be vectors of equal length.');
assert(isfield(g,'var2node') && numel(g.var2node) >= max(affected_vars(:)), ...
  'find_affected:Var2NodeMissing', 'g.var2node missing or too small.');

% ---- 1) Threshold variables ---------------------------------------------
significant = abs(dxa(:)) > dx_th;
affected_vars = affected_vars(:);
affected_vars = affected_vars(significant);

if isempty(affected_vars)
  affected_nodes    = [];
  affected_edge_ids = [];
  return;
end

% ---- 2) Vars -> node ids ------------------------------------------------
affected_nodes = unique(g.var2node(affected_vars));

% ---- 3) Collect all vars & edges per affected node ----------------------
nNodes = numel(affected_nodes);
var_cells  = cell(nNodes,1);
edge_cells = cell(nNodes,1);

for i = 1:nNodes
  nid = affected_nodes(i);
  key = sprintf('id%d', nid);

  if ~isfield(g.idLookup, key)
    % Node not known in lookup; skip gracefully
    var_cells{i}  = []; %#ok<*AGROW>
    edge_cells{i} = [];
    continue;
  end

  node = g.idLookup.(key);

  % Variables belonging to this node (convert 0-based offset to 1-based)
  if isfield(node,'offset') && isfield(node,'dimension') && node.dimension > 0
    vstart = node.offset + 1;
    var_cells{i} = vstart:(vstart + node.dimension - 1);
  else
    var_cells{i} = [];
  end

  % Incident edges (optional field)
  if isfield(node,'edges') && ~isempty(node.edges)
    edge_cells{i} = node.edges(:);
  else
    edge_cells{i} = [];
  end
end

% ---- 4) Flatten & deduplicate -------------------------------------------
if any(cellfun(@isempty, var_cells))
  % If some cells are empty, vertcat([]) is fine; guard for all-empty
  if all(cellfun(@isempty, var_cells)), affected_vars = []; else, affected_vars = unique(vertcat(var_cells{:})); end
else
  affected_vars = unique(vertcat(var_cells{:}));
end

if isempty(affected_vars)
  % No variables after expansion: no edges either
  affected_edge_ids = [];
  return;
end

if all(cellfun(@isempty, edge_cells))
  affected_edge_ids = [];
else
  affected_edge_ids = unique(vertcat(edge_cells{:}));
end
end
