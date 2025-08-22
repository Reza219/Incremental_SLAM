function phi = normalize_angle(phi)
% NORMALIZE_ANGLE  Wrap angles to the interval (-pi, pi].
%
% SYNTAX
%   phi = normalize_angle(phi)
%
% INPUT
%   phi : scalar or vector of angles (radians)
%
% OUTPUT
%   phi : same size as input, with all values wrapped into (-pi, pi]
%
% NOTES
%   - Works elementwise for vectors/matrices of angles.
%   - Equivalent to atan2(sin(phi), cos(phi)), but faster.

phi = mod(phi + pi, 2*pi) - pi;
end
