function [R, t] = kabsch2d_local(P, Q)
% KABSCH2D_LOCAL  Rigid alignment in 2D (no scale).
%
% SYNTAX
%   [R, t] = kabsch2d_local(P, Q)
%
% INPUTS
%   P : N×2 source points (rows are points)
%   Q : N×2 target points (rows are points)
%
% OUTPUTS
%   R : 2×2 rotation matrix (det(R)=+1)
%   t : 2×1 translation vector so that Q ≈ P*R' + t'
%
% DESCRIPTION
%   Minimizes ‖Q − (P*R' + t)‖_F over rotations R and translation t.
%   Implementation follows the Kabsch algorithm:
%     1) Center P and Q
%     2) H = X' * Y
%     3) SVD(H) = U Σ V'
%     4) R = V * U' (with reflection fix)
%     5) t = μ_Q' − R * μ_P'
%
% NOTES
%   - Requires N ≥ 2 non-collinear points for a well-conditioned estimate.
%   - Works in both MATLAB and GNU Octave.

  % ---- Validate inputs ---------------------------------------------------
  if ~ismatrix(P) || ~ismatrix(Q) || size(P,2)~=2 || size(Q,2)~=2
    error('kabsch2d_local:BadSize','P and Q must be N×2.');
  end
  if size(P,1) ~= size(Q,1) || size(P,1) < 2
    error('kabsch2d_local:BadN','P and Q must have the same N ≥ 2.');
  end
  if any(~isfinite(P(:))) || any(~isfinite(Q(:)))
    error('kabsch2d_local:NaNInf','Inputs contain NaN/Inf.');
  end

  % ---- Center the points -------------------------------------------------
  muP = mean(P, 1);       % 1×2
  muQ = mean(Q, 1);       % 1×2
  X = P - muP;            % centered P
  Y = Q - muQ;            % centered Q

  % ---- Cross-covariance and SVD -----------------------------------------
  H = X' * Y;             % 2×2
  [U, ~, V] = svd(H);     % H = U Σ V'

  % ---- Rotation with reflection fix -------------------------------------
  R = V * U';             % 2×2
  if det(R) < 0
    % Ensure a proper rotation (det=+1)
    V(:,end) = -V(:,end);
    R = V * U';
  end

  % ---- Translation -------------------------------------------------------
  % Q ≈ P*R' + t  ⇒  t' = μ_Q - μ_P * R'
  t = muQ' - R * muP';

  % (Optional) tiny numerical cleanup
  % R = real(R); t = real(t);
end
