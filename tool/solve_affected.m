function x_subset = solve_affected(R, b_perm, affected, p, parent, alpha)
% SOLVE_AFFECTED  Partial solve for selected variables using a cached block–Schur strategy.
%
% SYNTAX
%   x_subset = solve_affected(R, b_perm, affected, p, parent, alpha)
%
% DESCRIPTION
%   Solves A*x = b but returns only x(affected). Assumes R is an upper
%   Cholesky factor of A(p,p) (i.e., A_perm = R' * R), where p is the
%   permutation used to form A(p,p). The method:
%     1) Maps the requested variables to the permuted domain.
%     2) Computes the elimination-tree reach (dynamic set).
%     3) Reorders to [static; reach], reuses cached solves on the static block.
%     4) Solves the reduced system and scatters back to original order.
%
% INPUTS
%   R        : (n×n) upper-triangular Cholesky of A(p,p)
%   b_perm   : (n×1) permuted RHS (b_perm = J(:,p)' * r)
%   affected : vector of variable indices in the ORIGINAL (unpermuted) order
%   p        : permutation vector mapping perm→orig (A_perm = A(p,p)); default 1:n
%   parent   : elimination tree parent array over the permuted ordering (n×1), >=0
%   alpha    : scalar in (0,1]; if numel(affected)/n > alpha, do full solve (default 0.5)
%
% OUTPUT
%   x_subset : x(affected) exactly as from a full solve
%
% NOTES
%   - MATLAB/Octave compatible. Uses persistent cache keyed by the static set.
%   - If affected is empty, returns [] immediately.

  n = numel(b_perm);
  if nargin < 4 || isempty(p),       p = (1:n)'; end
  if nargin < 5 || isempty(parent),  parent = zeros(n,1); end
  if nargin < 6 || isempty(alpha),   alpha = 0.5; end

  % --- Basic validation
  assert(isvector(b_perm) && numel(p)==n && all(p>=1 & p<=n), 'solve_affected:BadInputs');
  assert(isequal(size(R), [n,n]), 'solve_affected:BadR');
  affected = affected(:);
  if isempty(affected), x_subset = []; return; end
  assert(all(affected>=1), 'solve_affected:BadAffected');

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
      if pi > 0, stack(end+1) = pi; end %#ok<AGROW>
    end
  end
  reach  = find(visited);                % dynamic block (R)
  static = setdiff((1:n)', reach);       % static block (U), stays unchanged

  % Local reordering: [static; reach]
  idx_order = [static; reach];
  R2 = R(idx_order, idx_order);
  b2 = b_perm(idx_order);
  nU = numel(static);

  % === CACHE STATIC BLOCK SOLVE ===
  persistent U_prev Ruu Tsolve Usolve
  if isempty(U_prev) || ~isequal(U_prev, static)
    if nU > 0
      Ruu = R2(1:nU, 1:nU);
      Tsolve = @(rhs) Ruu' \ rhs;   % forward (since Ruu is upper)
      Usolve = @(y)    Ruu  \ y;    % back-substitution
    else
      % Empty static part
      Ruu = []; Tsolve = @(rhs) rhs; Usolve = @(y) y;
    end
    U_prev = static;
  end

  % Static block solve
  if nU > 0
    y_u = Tsolve(b2(1:nU));
    x_u = Usolve(y_u);
  else
    y_u = []; x_u = [];
  end

  % Dynamic RHS and solve
  Rur = R2(1:nU, nU+1:end);
  b_eff = b2(nU+1:end) - (Rur.' * y_u);
  Rrr = R2(nU+1:end, nU+1:end);
  y_r = Rrr' \ b_eff;
  x_r = Rrr  \ y_r;

  % Reconstruct permuted solution on requested support
  x_local = [x_u; x_r];
  x_perm = zeros(n,1);
  x_perm(idx_order) = x_local;

  % Scatter to original ordering and return requested entries
  x_full = zeros(n,1);
  x_full(p) = x_perm;
  x_subset = x_full(affected);
end
