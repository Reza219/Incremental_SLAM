function [R, b, J, r, edge_J_index, p, parent] = linearize_new(J, r, edge_J_index, g, edge_ids, g_full)
% LINEARIZE_NEW  Add new edges to (J,r) and refactor (symbolic + numeric).
%
% SYNTAX
%   [R, b, J, r, edge_J_index, p, parent] = ...
%       linearize_new(J, r, edge_J_index, g, edge_ids, g_full)
%
% INPUTS
%   J, r            : existing global Jacobian and residual
%   edge_J_index    : map from global edge id -> starting row in J/r
%                     (containers.Map OR numeric vector; if numeric, 0 means unset)
%   g               : current graph (current linearization point in g.x)
%   edge_ids        : vector (or cell) of global edge indices to append
%   g_full          : full graph (for edge definitions)
%
% OUTPUTS
%   R        : upper-triangular Cholesky factor of A(p,p), A = J'*J
%   b        : permuted RHS, b = J(:,p)' * r
%   J, r     : updated Jacobian and residual (with new edge rows appended)
%   edge_J_index : updated edge->row mapping (start row of each edge block)
%   p        : fill-reducing column permutation (symamd on A by default)
%   parent   : elimination tree (etree) of A(p,p)
%
% NOTES
%   - Uses Jacobian_edge_Jr for per-edge blocks (which whitens via sqrtm).
%   - Assumes each edge contributes a fixed number of rows equal to numel(measurement).
%   - Expands J’s column count to match numel(g.x) if needed.

  % ---- Normalize and validate -------------------------------------------
  if iscell(edge_ids), edge_ids = cell2mat(edge_ids(:)).'; else, edge_ids = edge_ids(:).'; end
  state_size = numel(g.x);
  [rows_J, cols_J] = size(J);

  % Expand Jacobian width if the state grew
  if state_size > cols_J
    J(:, end+1:state_size) = 0;  %#ok<AGROW>
  end

  % Prepare to collect new blocks
  num_new = numel(edge_ids);
  if num_new == 0
    % No new edges: just (re)factorize with a fresh permutation
    p = symamd(J' * J);
    Jp = J(:, p);
    A  = Jp' * Jp;
    b  = Jp' * r;
    R  = chol(A);
    parent = etree(spones(A));
    return;
  end

  J_cells = cell(num_new, 1);
  r_cells = cell(num_new, 1);

  % Next free row
  next_row = rows_J + 1;

  % ---- Build and append each new edge block -----------------------------
  for i = 1:num_new
    eid = edge_ids(i);
    if eid < 1 || eid > numel(g_full.edges)
      error('linearize_new:EdgeOutOfRange','Edge id %d out of range.', eid);
    end
    e = g_full.edges(eid);

    % Per-edge Jacobian/residual at current linearization
    [Je, re] = Jacobian_edge_Jr(e, g, state_size);
    m = size(Je,1);

    J_cells{i} = Je;
    r_cells{i} = re;

    % Record starting row for this edge
    if isa(edge_J_index, 'containers.Map')
      edge_J_index(eid) = next_row;
    else
      % numeric vector case
      nmap = numel(edge_J_index);
      if eid > nmap
        edge_J_index(nmap+1:eid) = 0; %#ok<AGROW>
      end
      edge_J_index(eid) = next_row;
    end

    next_row = next_row + m;
  end

  % Append to J and r
  if ~isempty(J_cells)
    J_new = vertcat(J_cells{:});
    r_new = vertcat(r_cells{:});
    J = [J; J_new]; %#ok<AGROW>
    r = [r; r_new]; %#ok<AGROW>
  end

  % ---- Symbolic reordering & numeric factorization ----------------------
  % Use symamd on the normal matrix for stability across MATLAB/Octave
  A = J' * J;
  p = symamd(A);          % fill-reducing permutation
  Jp = J(:, p);
  A  = Jp' * Jp;
  b  = Jp' * r;

  % Cholesky (upper) and etree of sparsity
  R = chol(A);
  parent = etree(spones(A));
end
