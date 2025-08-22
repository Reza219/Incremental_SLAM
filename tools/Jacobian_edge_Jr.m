function [Je, re] = Jacobian_edge_Jr(edge, g, state_size)
% JACOBIAN_EDGE_JR  Linearize a single edge into global Jacobian/residual.
%
% SYNTAX
%   [Je, re] = Jacobian_edge_Jr(edge, g, state_size)
%
% DESCRIPTION
%   Computes the whitened contribution (J_e, r_e) for one edge, supporting:
%     - 'P' : pose–pose
%     - 'L' : pose–landmark
%     - 'G' : GNSS/prior on pose (unary; uses .from only, .toIdx may be 0)
%
%   Whitening uses L = sqrtm(W) with W = edge.information (symmetrized).
%
% INPUTS
%   edge       : struct with fields .type in {'P','L','G'}, .fromIdx, .toIdx,
%                .measurement, .information (m×m SPD)
%   g          : graph struct with field .x (current linearization point)
%   state_size : total number of state variables (length of full x)
%
% OUTPUTS
%   Je : sparse(m, state_size) – whitened Jacobian block for this edge
%   re : (m×1) whitened residual for this edge
%
% NOTES
%   - Pose block size = 3, landmark block size = 2.
%   - Indexing is 1‑based; unary edges should have edge.toIdx = 0 or empty.

% ---- Constants -----------------------------------------------------------
dim_pose     = 3;  % [x y theta]
dim_landmark = 2;  % [x y]

% ---- Validate (lightweight) ---------------------------------------------
if ~isfield(edge,'type') || isempty(edge.type)
  error('Jacobian_edge_Jr:NoType','Edge .type is missing.');
end
if ~isfield(edge,'fromIdx') || isempty(edge.fromIdx)
  error('Jacobian_edge_Jr:NoFromIdx','Edge .fromIdx is missing.');
end
if ~isfield(edge,'information') || isempty(edge.information)
  error('Jacobian_edge_Jr:NoInfo','Edge .information is missing.');
end

% ---- Dispatch by type ----------------------------------------------------
switch edge.type
  case 'P'  % pose–pose
    dim_from = dim_pose;  dim_to = dim_pose;
    xf = g.x(edge.fromIdx : edge.fromIdx + dim_from - 1);
    xt = g.x(edge.toIdx   : edge.toIdx   + dim_to   - 1);
    [e, A, B] = linearize_pose_pose(xf, xt, edge.measurement);

  case 'L'  % pose–landmark
    dim_from = dim_pose;  dim_to = dim_landmark;
    xf = g.x(edge.fromIdx : edge.fromIdx + dim_from - 1);
    xt = g.x(edge.toIdx   : edge.toIdx   + dim_to   - 1);
    [e, A, B] = linearize_pose_landmark(xf, xt, edge.measurement);

  case 'G'  % GNSS/prior on pose (unary)
    dim_from = dim_pose;  dim_to = 0;
    xf = g.x(edge.fromIdx : edge.fromIdx + dim_from - 1);
    z  = edge.measurement(:);
    if numel(z) ~= 2
      error('Jacobian_edge_Jr:BadGNSS','GNSS measurement must be length 2.');
    end
    [e, A] = linearize_pose_gnss(xf, z);
    B = [];

  otherwise
    error('Jacobian_edge_Jr:BadType','Unsupported edge type: %s', edge.type);
end

% ---- Whitening with sqrtm ------------------------------------------------
m = numel(e);
W = edge.information;
W = (W + W.')/2;         % symmetrize defensively
L = real(sqrtm(W));      % ensure real in case of tiny imag parts

% ---- Assemble sparse global Jacobian ------------------------------------
Je = sparse(m, state_size);

% Bounds checks (defensive)
if edge.fromIdx < 1 || edge.fromIdx + dim_from - 1 > state_size
  error('Jacobian_edge_Jr:FromRange','fromIdx out of range.');
end
Je(:, edge.fromIdx : edge.fromIdx + dim_from - 1) = L * A;

if dim_to > 0 && ~isempty(edge.toIdx) && edge.toIdx > 0
  if edge.toIdx < 1 || edge.toIdx + dim_to - 1 > state_size
    error('Jacobian_edge_Jr:ToRange','toIdx out of range.');
  end
  Je(:, edge.toIdx : edge.toIdx + dim_to - 1) = L * B;
end

% ---- Whitened residual ---------------------------------------------------
re = L * e;
end
