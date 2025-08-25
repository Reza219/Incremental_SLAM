function [loop_closure, state] = detect_loop_closure_unordered(edges, state, order_gap)
% DETECT_LOOP_CLOSURE_UNORDERED  Loop-closure flag without assuming pose-ID order.
%
% SYNTAX
%   [loop_closure, state] = detect_loop_closure_unordered(edges, state, order_gap)
%
% DESCRIPTION
%   Returns true if any edge in 'edges' indicates a non-local revisit:
%   - 'P' (pose–pose): both poses seen before and their first-seen orders
%     differ by more than order_gap.
%   - 'L' (pose–landmark): a landmark is re-observed by a pose whose first-seen
%     order differs by more than order_gap from any previous observing pose.
%   - 'G' (GNSS unary) never triggers a loop closure.
%
% INPUTS
%   edges     : struct array with fields .type in {'P','L','G'}, .from, .to
%   state     : running state (struct; may be empty on first call)
%   order_gap : positive scalar threshold on order difference (default 5)
%
% OUTPUTS
%   loop_closure : logical scalar
%   state        : updated state with fields:
%                    .pose_order : Map (double->double) pose_id -> first-seen order
%                    .lm2orders  : Map (double->any)    landmark_id -> vector of pose orders
%                    .next_order : next order index (double)
%
% NOTES
% - MATLAB R2020b+ and GNU Octave 7+ compatible (no string class).
% - This is a lightweight heuristic; it does not inspect geometry, only visit order.

% ---- Defaults & validation ----------------------------------------------
if nargin < 3 || isempty(order_gap) || ~isscalar(order_gap) || ~isfinite(order_gap)
    order_gap = 5;
end
loop_closure = false;

% Initialize/repair 'state' if missing or wrong type/shape
needs_reset = (nargin < 2) || isempty(state) || islogical(state) || ~isstruct(state) ...
    || ~isfield(state,'pose_order') || ~isa(state.pose_order,'containers.Map') ...
    || ~isfield(state,'lm2orders')  || ~isa(state.lm2orders,'containers.Map') ...
    || ~isfield(state,'next_order') || ~isscalar(state.next_order);

if needs_reset
    state = struct( ...
        'pose_order', containers.Map('KeyType','double','ValueType','double'), ...
        'lm2orders',  containers.Map('KeyType','double','ValueType','any'), ...
        'next_order', 1);
end

if isempty(edges), return; end

% ---- Main scan -----------------------------------------------------------
for e = 1:numel(edges)
    ei = edges(e);

    % Defensive: tolerate missing fields/types
    et = 'X';
    if isfield(ei,'type') && ~isempty(ei.type), et = ei.type; end

    switch et
        case 'P'  % pose–pose (binary)
            if ~isfield(ei,'from') || ~isfield(ei,'to'), continue; end
            i = ei.from; j = ei.to;

            i_in = isKey(state.pose_order, i);
            j_in = isKey(state.pose_order, j);

            if i_in && j_in
                oi = state.pose_order(i);
                oj = state.pose_order(j);
                if abs(oi - oj) > order_gap
                    loop_closure = true; return;
                end
            end

            % Assign first-seen orders AFTER checking
            if ~i_in
                state.pose_order(i) = state.next_order; state.next_order = state.next_order + 1;
            end
            if ~j_in
                state.pose_order(j) = state.next_order; state.next_order = state.next_order + 1;
            end

        case 'L'  % pose–landmark (binary: from=pose, to=landmark)
            if ~isfield(ei,'from') || ~isfield(ei,'to'), continue; end
            pid = ei.from; lid = ei.to;

            % Ensure pose has an order
            if ~isKey(state.pose_order, pid)
                state.pose_order(pid) = state.next_order; state.next_order = state.next_order + 1;
            end
            po = state.pose_order(pid);

            % Check non-local re-observation
            if isKey(state.lm2orders, lid)
                prev = state.lm2orders(lid);
                if any(abs(prev(:) - po) > order_gap)
                    loop_closure = true; return;
                end
                state.lm2orders(lid) = [prev(:).' po];  % append
            else
                state.lm2orders(lid) = po;
            end

        case 'G'  % GNSS unary (pose only)
            if ~isfield(ei,'from'), continue; end
            pid = ei.from;
            if ~isKey(state.pose_order, pid)
                state.pose_order(pid) = state.next_order; state.next_order = state.next_order + 1;
            end

        otherwise
            % Unknown type: ignore
    end
end
