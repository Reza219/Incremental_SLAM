function [e, A, B] = linearize_pose_pose(x1, x2, z)
% LINEARIZE_POSE_POSE  Linearize a 2D relative pose constraint.
%
% SYNTAX
%   [e, A, B] = linearize_pose_pose(x1, x2, z)
%
% INPUTS
%   x1 : 3x1 [x1; y1; th1]   first pose in world frame
%   x2 : 3x1 [x2; y2; th2]   second pose in world frame
%   z  : 3x1 [zx; zy; zth]   measured relative pose of x2 w.r.t. x1 (in x1 frame)
%
% OUTPUTS
%   e : 3x1 residual
%       e(1:2) = Rz' * ( R1' * (t2 - t1) - z(1:2) )
%       e(3)   = wrap( th2 - th1 - z(3) )
%   A : 3x3 Jacobian wrt x1
%   B : 3x3 Jacobian wrt x2
%
% NOTES
%   - R1 = [c1 -s1; s1 c1], Rz = [cz -sz; sz cz].
%   - Whitening is done outside this function.

% ---- Validate ------------------------------------------------------------
if numel(x1)~=3 || numel(x2)~=3 || numel(z)~=3
  error('linearize_pose_pose:BadSize','x1,x2,z must be 3x1.');
end
x1 = x1(:); x2 = x2(:); z = z(:);

% ---- Shorthands ----------------------------------------------------------
x_1 = x1(1); y_1 = x1(2); th1 = x1(3);
x_2 = x2(1); y_2 = x2(2); th2 = x2(3);
cz = cos(z(3)); sz = sin(z(3));
c1 = cos(th1);  s1 = sin(th1);

% Rotations
R1  = [ c1, -s1;  s1,  c1];
Rz  = [ cz, -sz;  sz,  cz];

% ---- Residual ------------------------------------------------------------
delta_t  = [x_2 - x_1; y_2 - y_1];        % world-frame translation difference
pred_xy  = R1' * delta_t;                  % predicted relative translation in x1 frame
e_xy     = Rz' * (pred_xy - z(1:2));       % error expressed in measurement frame
e_th     = normalize_angle(th2 - th1 - z(3));
e        = [e_xy; e_th];

% ---- Jacobians -----------------------------------------------------------
% Q1 and Q2 per standard derivation (same as your original)
Q_1 = [ -c1, -s1,  (x_1 - x_2)*s1 + (y_2 - y_1)*c1;
         s1, -c1,  (x_1 - x_2)*c1 + (y_1 - y_2)*s1 ];

Q_2 = [  c1,  s1,  0;
        -s1,  c1,  0 ];

A = [ Rz' * Q_1;  0, 0, -1 ];
B = [ Rz' * Q_2;  0, 0,  1 ];
end
