function plot_graph(g, affected_nodes, edge_ids, draw_edges, c, xlimits, ylimits)
% PLOT_GRAPH  Visualize a 2D SLAM graph: poses, landmarks, affected nodes, and edges.
%
% SYNTAX
%   plot_graph(g, affected_nodes, edge_ids, draw_edges, c, xlimits, ylimits)
%
% INPUTS
%   g              : graph struct (.x, .edges, .idLookup)
%   affected_nodes : vector of node ids to highlight (optional, can be [])
%   edge_ids       : vector or cell array of edge indices to draw (optional)
%   draw_edges     : logical, whether to draw edges (default true)
%   c              : color char for affected nodes (e.g., 'k'); optional
%   xlimits        : [xmin xmax] or []; optional
%   ylimits        : [ymin ymax] or []; optional
%
% NOTES
%   - GNSS ('G') stars are plotted only once, when their edges first appear.
%   - Uses 0-based offsets stored in g.idLookup; converts to 1-based for indexing.

    if nargin < 4 || isempty(draw_edges), draw_edges = true; end
    if nargin < 5 || isempty(c), c = 'k'; end
    if nargin < 6, xlimits = []; end
    if nargin < 7, ylimits = []; end
    if nargin < 2 || isempty(affected_nodes), affected_nodes = []; end
    if nargin < 3, edge_ids = []; end

    % Normalize edge_ids to a row numeric vector
    if iscell(edge_ids), edge_ids = cell2mat(edge_ids(:)).'; else, edge_ids = edge_ids(:).'; end

    % Persistent figure/axes + GNSS watermark
    persistent fig ax last_plotted_gnss_edge_idx
    if isempty(fig) || ~ishandle(fig)
        fig = figure('Name','SLAM Graph','NumberTitle','off','Color','w',...
                     'Visible','on','CloseRequestFcn',@(src,~) delete(src));
        ax = axes('Parent',fig,'Color','w'); axis(ax,'equal'); grid(ax,'on');
        xlabel(ax,'x'); ylabel(ax,'y');
    end
    if isempty(last_plotted_gnss_edge_idx), last_plotted_gnss_edge_idx = 0; end

    cla(ax); hold(ax,'on');

    % ---- Poses & landmarks (scatter) -------------------------------------
    [pose_off0, lm_off0] = get_poses_landmarks(g);   % 0-based offsets
    if ~isempty(lm_off0)
        % landmarks are 2D blocks [x;y]
        lm_x = g.x(lm_off0+1);
        lm_y = g.x(lm_off0+2);
        plot(ax, lm_x, lm_y, 'ro', 'MarkerSize', 4, 'DisplayName','Landmarks');
    end
    if ~isempty(pose_off0)
        % poses are 3D blocks [x;y;theta]
        p_x = g.x(pose_off0+1);
        p_y = g.x(pose_off0+2);
        plot(ax, p_x, p_y, 'b.', 'MarkerSize', 6, 'DisplayName','Poses');
    end

    % ---- GNSS stars: plot only for *new* GNSS edges ----------------------
    if ~isempty(edge_ids)
        gnss_xy = []; new_max = last_plotted_gnss_edge_idx;
        for ii = edge_ids
            e = g.edges(ii);
            if isfield(e,'type') && e.type=='G' && isfield(e,'measurement') && numel(e.measurement)>=2 ...
                    && ii > last_plotted_gnss_edge_idx
                gnss_xy(:, end+1) = e.measurement(1:2); %#ok<AGROW>
                if ii > new_max, new_max = ii; end
            end
        end
        if ~isempty(gnss_xy)
            % De-duplicate within this batch
            [~, ia] = unique(gnss_xy.', 'rows', 'stable');
            gnss_xy = gnss_xy(:, ia);
            plot(ax, gnss_xy(1,:), gnss_xy(2,:), 'r*', 'MarkerSize', 6, 'DisplayName','GNSS priors');
            last_plotted_gnss_edge_idx = new_max;
        end
    end

    % ---- Affected nodes (highlight) --------------------------------------
    for nid = affected_nodes(:).'
        key = sprintf('id%d', nid);
        if isfield(g.idLookup, key)
            off1 = g.idLookup.(key).offset + 1;  % 1-based
            plot(ax, g.x(off1), g.x(off1+1), [c 'o'], 'MarkerSize', 6, 'LineWidth', 1.0);
        end
    end

    % ---- Edges (optional) ------------------------------------------------
    if draw_edges && ~isempty(edge_ids)
        pose_segs = {}; landmark_segs = {}; gnss_segs = {};
        for ii = edge_ids
            e = g.edges(ii);
            if ~isfield(e,'type') || isempty(e.type), continue; end
            switch e.type
                case 'P'
                    pt1 = g.x(e.fromIdx : e.fromIdx+1);
                    pt2 = g.x(e.toIdx   : e.toIdx+1);
                    pose_segs{end+1} = [pt1, pt2, [NaN;NaN]]; %#ok<AGROW>
                case 'L'
                    pt1 = g.x(e.fromIdx : e.fromIdx+1);
                    pt2 = g.x(e.toIdx   : e.toIdx+1);
                    landmark_segs{end+1} = [pt1, pt2, [NaN;NaN]]; %#ok<AGROW>
                case 'G'
                    pt1 = g.x(e.fromIdx : e.fromIdx+1);
                    z   = e.measurement(1:2);
                    gnss_segs{end+1} = [pt1, z(:), [NaN;NaN]]; %#ok<AGROW>
            end
        end

        if ~isempty(pose_segs)
            all_pose = horzcat(pose_segs{:});
            plot(ax, all_pose(1,:), all_pose(2,:), ':b', 'LineWidth', 0.5, 'HandleVisibility','off');
        end
        if ~isempty(landmark_segs)
            all_lm = horzcat(landmark_segs{:});
            plot(ax, all_lm(1,:), all_lm(2,:), ':r', 'LineWidth', 0.5, 'HandleVisibility','off');
        end
        if ~isempty(gnss_segs)
            all_g = horzcat(gnss_segs{:});
            plot(ax, all_g(1,:), all_g(2,:), ':r', 'LineWidth', 0.5, 'HandleVisibility','off');
        end
    end

    % ---- Axes & finish ----------------------------------------------------
    if ~isempty(xlimits) && ~isempty(ylimits)
        xlim(ax, xlimits); ylim(ax, ylimits);
    else
        axis(ax,'tight'); % keeps equal aspect from earlier
    end
    hold(ax,'off');
end
