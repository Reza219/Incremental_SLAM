function A = invt(m)
% INVT  Invert a 2D homogeneous transform (SE2).
%
%   A = invt(m)
%
% INPUT
%   m : 3x3 homogeneous transform [R t; 0 0 1]
%       where R ∈ SO(2), t ∈ ℝ²
%
% OUTPUT
%   A : 3x3 inverse homogeneous transform = [R'  -R'*t; 0 0 1]

  R = m(1:2,1:2);
  t = m(1:2,3);

  Rt = R.';               % since R is rotation
  A = [Rt, -Rt*t; 0 0 1];
end
