function v = t2v(A)
%T2V Convert 2D homogeneous transform(s) to pose vector(s) [x; y; theta].
%
%   v = T2V(A)
%   A : 3x3 or 3x3xN homogeneous transform(s)
%   v : 3x1 or 3xN pose vector(s) [x; y; theta]

    if ndims(A) == 2
        v = [A(1:2,3); atan2(A(2,1), A(1,1))];
    elseif ndims(A) == 3
        N = size(A,3);
        v = zeros(3,N);
        for k = 1:N
            v(:,k) = [A(1:2,3,k); atan2(A(2,1,k), A(1,1,k))];
        end
    else
        error('Input must be 3x3 or 3x3xN homogeneous transform(s).');
    end
end
