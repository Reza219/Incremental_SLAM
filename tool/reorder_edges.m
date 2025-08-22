function g = reorder_edges(g)
% REORDER_EDGES  Reorder graph edges to simulate incremental SLAM acquisition.
%
% SYNTAX
%   g = reorder_edges(g)
%
% DESCRIPTION
%   Produces an ordering where each iteration (i) grows the initialized set
%   by exactly one new node via a single *growth* edge (never GNSS), and then
%   (ii) appends all *closure* edges whose endpoints are now initialized.
%
% CONTRACT
%   - Starts from known initial pose id0 (must exist in g.idLookup).
%   - GNSS ('G') edges are never used for growth; they are appended only
%     after the pose exists and kept unary (toIdx = 0).
%   - Pose/landmark edges ('P','L') can be growth or closure depending on
%     whether they connect initialized↔uninitialized or initialized↔initialized.
%
% INPUT / OUTPUT
%   g (struct) with fields:
%     .x        (nx1 double)         : state vector (preserved)
%     .edges    (struct array)       : fields .type, .from, .to, .measurement,
%                                      .information, .fromIdx, .toIdx (indices will be set)
%     .idLookup (struct of structs)  : fields id%d -> .offset (0-based), .dimension
%
% NOTES
%   - MATLAB R2020b+ and GNU Octave 7+ compatible (no `string` class used).
%   - Only reorders g.edges; numeric values are untouched.
%
% See also READ_GRAPH_G2O, OPTIMIZE_INCREMENTAL, LINEARIZE_NEW_EDGES.

% ---- Validate inputs -----------------------------------------------------
assert(isstruct(g) && isfield(g,'edges') && isfield(g,'idLookup'), ...
  'reorder_edges:InvalidGraph', 'g must have .edges and .idLookup');

% Require id0 as the seed
if ~isfield(g.idLookup,'id0')
  error('reorder_edges:MissingId0', 'Vertex id0 not found in g.idLookup.');
end

E = g.edges;
if isempty(E)
  g.edges = E; % nothing to do
  return;
end

% ---- Helpers -------------------------------------------------------------
% Fast lookup from numeric id -> 1-based state offset
get_off = @(nid) g.idLookup.(sprintf('id%d', nid)).offset + 1;

% Track which vertex ids are initialized
initialized = containers.Map('KeyType','double','ValueType','logical');
initialized(0) = true;  % start from id0

% Bookkeeping
N = numel(E);
inserted      = false(N,1);
ordered_edges = E([]);  % empty struct array with same fields

% ---- Main loop: grow by 1, then add all closures ------------------------
while true
  changed   = false;
  new_batch = E([]);    % collect this iteration's edges in order

  % --- (A) Choose exactly one growth edge (skip GNSS) --------------------
  growth_idx = 0;
  for i = 1:N
    if inserted(i), continue; end
    ei = E(i);

    % Treat missing .type as non-GNSS by default
    etype = 'X';
    if isfield(ei,'type') && ~isempty(ei.type)
      etype = ei.type;
    end

    if etype == 'G'
      continue; % never grow with GNSS
    end

    % We require both .from and .to to exist for P/L edges
    if ~isfield(ei,'from') || ~isfield(ei,'to')
      continue;
    end

    from_ok = isKey(initialized, ei.from);
    to_ok   = isKey(initialized, ei.to);

    if xor(from_ok, to_ok)
      % This edge connects initialized ↔ uninitialized: growth candidate
      ei.fromIdx = get_off(ei.from);
      ei.toIdx   = get_off(ei.to);
      new_batch  = [new_batch; ei]; %#ok<AGROW>
      inserted(i) = true;
      changed = true;
      growth_idx = i;

      % Mark both endpoints as initialized after growth
      initialized(ei.from) = true;
      initialized(ei.to)   = true;
      break; % exactly one growth per iteration
    end
  end

  % --- (B) Append all closure edges now possible -------------------------
  % (Includes GNSS unary priors for newly-initialized poses)
  for i = 1:N
    if inserted(i), continue; end
    ei = E(i);

    etype = 'X';
    if isfield(ei,'type') && ~isempty(ei.type)
      etype = ei.type;
    end

    if etype == 'G'
      % GNSS/prior is unary on .from; add once pose exists
      if isfield(ei,'from') && isKey(initialized, ei.from)
        ei.fromIdx = get_off(ei.from);
        ei.toIdx   = 0;  % keep unary
        new_batch  = [new_batch; ei]; %#ok<AGROW>
        inserted(i) = true;
        changed = true;
      end
      continue;
    end

    % P/L closures when both endpoints are initialized
    if isfield(ei,'from') && isfield(ei,'to') && ...
       isKey(initialized, ei.from) && isKey(initialized, ei.to)

      ei.fromIdx = get_off(ei.from);
      ei.toIdx   = get_off(ei.to);
      new_batch  = [new_batch; ei]; %#ok<AGROW>
      inserted(i) = true;
      changed = true;
    end
  end

  if ~changed
    break; % no growth and no closures: we're done
  end

  ordered_edges = [ordered_edges; new_batch]; %#ok<AGROW>

  % Optional: if we failed to find a growth edge but added closures,
  % we’ll loop once more; if we *still* cannot grow, remaining edges will
  % be dropped below (e.g., disconnected components without id0 path).
  if growth_idx == 0
    % no-op: the loop will either find a growth in the next pass or exit
  end
end

% ---- Finalize ------------------------------------------------------------
dropped = sum(~inserted);
if dropped > 0
  warning('reorder_edges:DroppedEdges', ...
    'Could not safely reorder %d edge(s); leaving them out (likely disconnected).', dropped);
end

g.edges = ordered_edges;
end
