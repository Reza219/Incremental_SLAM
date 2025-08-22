function [e, A] = linearize_pose_gnss(xf, z)
% LINEARIZE_POSE_GNSS  Linearize a 2D GNSS (x,y) prior w.r.t. a pose.
%
% SYNTAX
%   [e, A] = linearize_pose_gnss(xf, z)
%
% INPUTS
%   xf : 3x1 pose state [x; y; theta]
%   z  : 2x1 GNSS measurement [z_x; z_y]  (same units/frame as state)
%
% OUTPUTS
%   e : 2x1 residual  = [x - z_x; y - z_y]
%   A : 2x3 Jacobian  = [1 0 0; 0 1 0]
%
% NOTES
%   - This is a position-only unary factor; heading is not observed.
%   - Whitening is handled outside (e.g., L = sqrtm(W); use Je = L*A, re = L*e).

% ---- Validate (lightweight) ---------------------------------------------
if numel(xf) ~= 3
  error('linearize_pose_gnss:BadX', 'xf must be 3x1 [x;y;theta].');
end
if numel(z) ~= 2
  error('linearize_pose_gnss:BadZ', 'z must be 2x1 [z_x; z_y].');
end

xf = xf(:); z = z(:);

% ---- Residual and Jacobian ----------------------------------------------
e = [xf(1) - z(1);
     xf(2) - z(2)];

A = [1 0 0;
     0 1 0];
end
