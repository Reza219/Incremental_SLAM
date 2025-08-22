function [rmse, ids, E_xy_aligned, R, t] = compute_ate_rmse(g_est, g_gt)
% COMPUTE_ATE_RMSE  Absolute Trajectory Error (RMSE) after 2D rigid alignment.
%
% SYNTAX
%   [rmse, ids, E_xy_aligned, R, t] = compute_ate_rmse(g_est, g_gt)
%
% INPUTS
%   g_est : current/estimated graph (fields: .x, .idLookup)
%   g_gt  : ground-truth graph (fields: .x, .idLookup)
%
% OUTPUTS
%   rmse          : scalar ATE RMSE (meters) on translation only; NaN if <2 poses overlap
%   ids           : intersected pose ids used for the metric (in stable order)
%   E_xy_aligned  : aligned estimated XY positions (N×2), for plotting
%   R, t          : alignment rotation (2×2) and translation (2×1), such that
%                   G ≈ E*R' + t'
%
% NOTES
%   - A "pose" is any node with dimension==3 in idLookup ([x y theta]).
%   - No scale change (Kabsch 2D).
%   - If <2 shared pose ids, returns rmse = NaN.

  % --- collect pose ids
  ids_est = pose_ids_3dof(g_est);
  ids_gt  = pose_ids_3dof(g_gt);

  ids = intersect(ids_est(:), ids_gt(:), 'stable');
  if numel(ids) < 2
    rmse = NaN; E_xy_aligned = []; R = eye(2); t = [0;0];
    return;
  end

  % --- gather [x y theta] rows for these ids
  E = poses_from_ids_2d(g_est, ids);   % N×3
  G = poses_from_ids_2d(g_gt,  ids);   % N×3

  % --- rigid align EST->GT on XY (no scale)
  [R, t] = kabsch2d_local(E(:,1:2), G(:,1:2));

  % --- apply alignment and compute RMSE
  E_xy_aligned = (E(:,1:2) * R.') + repmat(t.', size(E,1), 1);
  dxy = G(:,1:2) - E_xy_aligned;
  rmse = sqrt(mean(sum(dxy.^2, 2)));
end

% ===== helpers (kept local to avoid cluttering tools/) ====================
function ids = pose_ids_3dof(g)
  keys = fieldnames(g.idLookup);
  ids = zeros(numel(keys),1); k = 0;
  for i=1:numel(keys)
    v = g.idLookup.(keys{i});
    if isfield(v,'dimension') && v.dimension==3
      nid = sscanf(keys{i}, 'id%d'); if isempty(nid), continue; end
      k=k+1; ids(k)=nid;
    end
  end
  ids = ids(1:k);
end

function X = poses_from_ids_2d(g, ids)
  % returns N×3 [x y theta] in the same order as ids
  X = zeros(numel(ids),3);
  for i=1:numel(ids)
    key = sprintf('id%d', ids(i));
    off = g.idLookup.(key).offset + 1;   % 1-based
    X(i,:) = g.x(off:off+2).';
  end
end
