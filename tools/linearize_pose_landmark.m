function [e, A, B] = linearize_pose_landmark(x, l, z)
% LINEARIZE_POSE_LANDMARK  Linearize a 2D pose–landmark bearing/range (xy) factor.
%
% SYNTAX
%   [e, A, B] = linearize_pose_landmark(x, l, z)
%
% INPUTS
%   x : 3x1 robot pose [x_i; y_i; theta_i] in world frame
%   l : 2x1 landmark position [l_x; l_y] in world frame
%   z : 2x1 measured landmark position in the *robot frame* at pose x
%       (i.e., z = R_i' * (l - t_i) + noise)
%
% OUTPUTS
%   e : 2x1 residual        = R_i' * (l - t_i) - z
%   A : 2x3 Jacobian w.r.t. x = [∂e/∂x_i, ∂e/∂y_i, ∂e/∂θ_i]
%   B : 2x2 Jacobian w.r.t. l = ∂e/∂l
%
% NOTES
%   - Uses R_i = [c -s; s c], with c = cos(theta_i), s = sin(theta_i).
%   - Whitening is applied outside this function.

% ---- Validate ------------------------------------------------------------
if numel(x) ~= 3, error('linearize_pose_landmark:BadX','x must be 3x1 [x;y;theta].'); end
if numel(l) ~= 2, error('linearize_pose_landmark:BadL','l must be 2x1 [l_x; l_y].'); end
if numel(z) ~= 2, error('linearize_pose_landmark:BadZ','z must be 2x1 [z_x; z_y] (robot frame).'); end

x = x(:); l = l(:); z = z(:);

xi = x(1); yi = x(2); th = x(3);
lx = l(1); ly = l(2);

c = cos(th); s = sin(th);

% Rotation (world->robot): R' = [ c  s; -s  c ]
Rt = [ c,  s;
      -s,  c];

% ---- Residual ------------------------------------------------------------
e = Rt * ([lx; ly] - [xi; yi]) - z;

% ---- Jacobians -----------------------------------------------------------
% A = ∂e/∂x = [∂e/∂x_i, ∂e/∂y_i, ∂e/∂θ_i]
A = [ -c,  -s,  (xi - lx) * s + (ly - yi) * c;
       s,  -c,  (xi - lx) * c + (yi - ly) * s ];

% B = ∂e/∂l = R' (since e = R'*(l - t) - z)
B = [  c,  s;
      -s,  c];
end
