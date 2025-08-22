function A = v2t(v)
%V2T Convert pose vector [x; y; theta] into 2D homogeneous transform.
%
%   A = V2T(v)
%   v : 3x1 vector [x; y; theta]
%   A : 3x3 homogeneous transform matrix

    assert(numel(v) == 3, 'v2t expects a 3x1 pose vector [x; y; theta].');

    c = cos(v(3));
    s = sin(v(3));

    A = [c, -s, v(1);
         s,  c, v(2);
         0,  0, 1];
end
