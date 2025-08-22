function [poses, landmarks] = get_poses_landmarks(g)
% GET_POSES_LANDMARKS  Collect 0-based offsets for pose and landmark blocks.
%
% SYNTAX
%   [poses, landmarks] = get_poses_landmarks(g)
%
% DESCRIPTION
%   Scans g.idLookup and returns the 0-based offsets of pose and landmark
%   variables based on their block dimensions:
%     - Pose blocks:    dimension == 3
%     - Landmark blocks: dimension == 2
%
% INPUT
%   g : graph struct with field
%       .idLookup : struct, keys 'id%d' -> .offset (0-based), .dimension
%
% OUTPUT
%   poses     : column vector of 0-based offsets for pose blocks (sorted)
%   landmarks : column vector of 0-based offsets for landmark blocks (sorted)
%
% NOTES
%   - Offsets are 0-based to match the convention used in idLookup.
%   - To convert to 1-based starting indices, use poses+1 or landmarks+1.

% ---- Validate ------------------------------------------------------------
assert(isstruct(g) && isfield(g,'idLookup'), ...
  'get_poses_landmarks:InvalidGraph', 'g.idLookup is required.');

keys = fieldnames(g.idLookup);
posesC = cell(numel(keys),1);
landC = cell(numel(keys),1);
np = 0; nl = 0;

for i = 1:numel(keys)
  v = g.idLookup.(keys{i});
  if ~isfield(v,'dimension') || ~isfield(v,'offset'), continue; end

  if v.dimension == 3
    np = np + 1; posesC{np} = v.offset;
  elseif v.dimension == 2
    nl = nl + 1; landC{nl} = v.offset;
  end
end

% Trim and concatenate
if np > 0, poses = sort(cell2mat(posesC(1:np)).'); else, poses = []; end
if nl > 0, landmarks = sort(cell2mat(landC(1:nl)).'); else, landmarks = []; end
end
