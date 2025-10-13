function x_subset = solve_affected(R, b_perm, affected, p, parent, alpha)
% SOLVE_AFFECTED  Partial solve for selected variables using a cached block–Schur strategy.
%
% SYNTAX
%   x_subset = solve_affected(R, b_perm, affected, p, parent, alpha)
%
% DESCRIPTION
%   Solves A*x = b but returns only x(affected). Assumes R is an upper
%   Cholesky factor of A(p,p) (i.e., A_perm = R' * R), where p is the
%   column permutation used to form A(p,p).
%   Steps:
%     1) Map requested variables to the permuted domain.
%     2) Compute elimination-tree reach (dynamic set).
%     3) Reorder to [static; reach]; reuse cached solves on the static block.
%     4) Solve reduced system and scatter back to original order.
%
% INPUTS
%   R        : (n×n) upper-triangular Cholesky of A(p,p)
%   b_perm   : (n×1) permuted RHS (b_perm = J(:,p)' * r)
%   affected : vector of variable indices in the ORIGINAL order
%   p        : permutation vector (size n) such that A_perm = A(p,p)
%   parent   : elimination tree over the permuted ordering (n×1), parent(i)=0/<=0 for root
%   alpha    : scalar in (0,1]; if numel(affected)/n > alpha, do full solve (default 0.5)
%
% OUTPUT
%   x_subset : x(affected) exactly as from a full solve
%
% NOTES
%   - MATLAB/Octave compatible.
%   - Handles empty static or dynamic blocks robustly (Octave-safe empties).

  n = numel(b_perm);
  if nargin < 4 || isempty(p),       p = (1:n)'; end
  if nargin < 5 || isempty(parent),  parent = zeros(n,1); end
  if nargin < 6 || isempty(alpha),   alpha = 0.5; end

  % --- Basic validation
  assert(isvector(b_perm) && numel(b_perm)==n, 'solve_affected:BadRHS');
  assert(isequal(size(R), [n,n]), 'solve_affected:BadR');
  assert(isvector(p) && numel(p)==n && all(p>=1 & p<=n), 'solve_affected:BadPerm');
  affected = affected(:);
  if isempty(affected), x_subset = []; return; end
  assert(all(affected>=1), 'solve_affected:BadAffectedIdx');

  % --- Full solve fallback if request is large
  if numel(affected)/n > alpha
    y = R' \ b_perm;
    dx_perm = R \ y;
    x_full = zeros(n,1);
    x_full(p) = dx_perm;         % scatter to original order
    x_subset = x_full(affected);
    return;
  end

  % === PARTIAL SOLVE ===

  % Inverse permutation: orig -> perm
  inv_p = zeros(n,1); inv_p(p) = (1:n)';
  affected_p = inv_p(affected);
  assert(all(affected_p>=1 & affected_p<=n), 'solve_affected:AffectedOutOfRange');

  % Elimination-tree reach from affected_p
  visited = false(n,1);
  stack = affected_p(:);
  while ~isempty(stack)
    i = stack(end); stack(end) = [];
    if ~visited(i)
      visited(i) = true;
      pi = parent(i);
      if ~isempty(pi) && pi > 0
        stack(end+1) = pi; %#ok<AGROW>
      end
    end
  end
  reach  = find(visited);           % dynamic block
  static = setdiff((1:n)', reach);  % static block

  % Local reordering: [static; reach]
  idx_order = [static; reach];
  R2 = R(idx_order, idx_order);
  b2 = b_perm(idx_order);
  nU = numel(static);
  nR = numel(reach);

  % === CACHE STATIC BLOCK SOLVE ===
  persistent U_prev Ruu Tsolve Usolve
  if isempty(U_prev) || ~isequal(U_prev, static)
    if nU > 0
      Ruu = R2(1:nU, 1:nU);
      Tsolve = @(rhs) Ruu' \ rhs;   % forward (since Ruu is upper)
      Usolve = @(y)    Ruu  \ y;    % back-substitution
    else
      % Empty static block
      Ruu = [];
      Tsolve = @(rhs) rhs;
      Usolve = @(y)    y;
    end
    U_prev = static;
  end

  % --- Solve static block (if present)
  if nU > 0
    y_u = Tsolve(b2(1:nU));
    x_u = Usolve(y_u);
  else
    y_u = []; x_u = [];
  end

  % --- Dynamic RHS and solve (Octave-safe guards for empty blocks)
  if nR > 0
    Rur   = R2(1:nU, nU+1:end);     % size: nU × nR (possibly 0×nR)
    Rrr   = R2(nU+1:end, nU+1:end); % size: nR × nR
    if nU > 0
      % Safe: only form product if static block exists
      b_eff = b2(nU+1:end) - (Rur.' * y_u);
    else
      b_eff = b2(nU+1:end);
    end
    % Solve dynamic block
    y_r = Rrr' \ b_eff;
    x_r = Rrr  \ y_r;
  else
    x_r = [];
  end

  % --- Reconstruct permuted solution
  x_perm_local = [x_u; x_r];        % matches idx_order
  x_perm = zeros(n,1);
  if ~isempty(idx_order)
    x_perm(idx_order) = x_perm_local;
  end

  % --- Scatter to original ordering and return requested entries
  x_full = zeros(n,1);
  x_full(p) = x_perm;
  x_subset = x_full(affected);
end
