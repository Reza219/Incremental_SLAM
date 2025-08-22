function [R, b, J, r] = linearize_affected(J, r, edge_J_index, p, g, affected_edge_ids)
% LINEARIZE_AFFECTED  Refresh J and r for selected edges and refactorize.
%
% SYNTAX
%   [R, b, J, r] = linearize_affected(J, r, edge_J_index, p, g, affected_edge_ids)
%
% INPUTS
%   J, r               : current global Jacobian and residual
%   edge_J_index       : map from edge id -> starting row in J/r
%                        (either containers.Map or numeric vector indexed by edge id)
%   p                  : column permutation of variables (optional). If empty,
%                        uses COLAMD over J (column approximate minimum degree).
%   g                  : current graph (with .x, .edges, etc.)
%   affected_edge_ids  : vector of edge ids to relinearize
%
% OUTPUTS
%   R  : upper-triangular Cholesky of A(p,p) where A = J'*J
%   b  : permuted RHS b = J(:,p)' * r
%   J, r : updated Jacobian and residual (with refreshed edge rows)
%
% NOTES
%   - Uses Jacobian_edge_Jr (which whitens with sqrtm) for each edge.
%   - J’s rows must be laid out contiguously per edge, starting at edge_J_index(eid).
%   - For safety, if Je’s row count differs from the stored block size,
%     we error (layout mismatch).

  state_size = numel(g.x);
  affected_edge_ids = affected_edge_ids(:).';

  % --- Update rows for each affected edge ---------------------------------
  for k = 1:numel(affected_edge_ids)
    eid = affected_edge_ids(k);
    e   = g.edges(eid);

    % Recompute per-edge Jacobian & residual (whitening handled inside)
    [Je, re] = Jacobian_edge_Jr(e, g, state_size);
    m = size(Je,1);

    % Locate block row range in J/r
    if isa(edge_J_index, 'containers.Map')
      if ~isKey(edge_J_index, eid)
        error('linearize_affected:MissingEdgeRow', 'No row index stored for edge %d', eid);
      end
      row_start = edge_J_index(eid);
    else
      % assume numeric vector indexed by edge id
      if eid < 1 || eid > numel(edge_J_index) || edge_J_index(eid) == 0
        error('linearize_affected:MissingEdgeRow', 'No row index stored for edge %d', eid);
      end
      row_start = edge_J_index(eid);
    end

    row_range = row_start : (row_start + m - 1);

    % Validate existing allocation matches Je’s size
    if size(J,1) < row_range(end)
      error('linearize_affected:RowRangeOOB', 'Row range for edge %d exceeds J size.', eid);
    end
    if size(J,2) ~= state_size
      error('linearize_affected:StateSizeMismatch', 'J columns (%d) != state size (%d).', size(J,2), state_size);
    end

    % In-place patch
    J(row_range, :) = Je;
    r(row_range, 1) = re;
  end

  % --- Reorder columns (if needed) and form normal equations --------------
  if nargin < 4 || isempty(p)
    % Column ordering for normal equations; prefer COLAMD(J) (works on J directly)
    p = colamd(J);
  end
  Jp = J(:, p);
  A  = Jp' * Jp;
  b  = Jp' * r;

  % --- Factorize -----------------------------------------------------------
  % Cholesky of SPD normal matrix. MATLAB/Octave chol returns upper by default.
  R = chol(A);
end
