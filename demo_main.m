% DEMO MAIN — Incremental 2D SLAM with dataset switch, gating, and SPO toggle
clear functions; clearvars; clc;
addpath('tools');

% =========================
% Config
% =========================
dataset   = 'MIT';          % {'MIT','Intel','CSAIL','FR079','FRH','MITP'}
gating    = 'IGG';          % {'IGG','LCG','none'}  -> info-gain / loop-closure / no gating
use_spo   = true;           % true = selective partial optimization (SPO); false = classic full GN
max_GNi   = 10;
batchSize = 1;

% =========================
% Dataset settings
% =========================
switch dataset
  case 'MIT'      % G2O
    data_file = fullfile('data','MITb_g2o.g2o');
    gt_file   = fullfile('data','MIT_ground_truth.mat');
    xlm = [-230, 40]; ylm = [-80, 150];
    ent_th = 1; dx_th = 1e-3; lc_gap = 4;

  case 'Intel'    % G2O
    data_file = fullfile('data','INTEL_g2o.g2o');
    gt_file   = fullfile('data','Intel_ground_truth.mat');
    xlm = [-16, 18]; ylm = [-27, 5];
    ent_th = 0.72; dx_th = 1e-6; lc_gap = 4;

  case 'CSAIL'    % TORO
    data_file = fullfile('data','CSAIL_P_toro.graph');
    gt_file   = fullfile('data','CSAIL_ground_truth.mat');
    xlm = [-10, 35]; ylm = [-35, 35];
    ent_th = 0.95; dx_th = 1e-5; lc_gap = 4;

  case 'FR079'    % TORO
    data_file = fullfile('data','FR079_P_toro.graph');
    gt_file   = fullfile('data','FR079_ground_truth.mat');
    xlm = [-25, 12]; ylm = [-8, 8];
    ent_th = 0.6; dx_th = 1e-4; lc_gap = 4;

  case 'FRH'      % TORO
    data_file = fullfile('data','FRH_P_toro.graph');
    gt_file   = fullfile('data','FRH_ground_truth.mat');
    xlm = [-10, 400]; ylm = [-140, 90];
    ent_th = 0.45; dx_th = 1e-7; lc_gap = 1;

  case 'MITP'     % MAT (graph stored as S.g)
    data_file = fullfile('data','MITP.mat');
    gt_file   = fullfile('data','MITP_ground_truth.mat');
    xlm = [-230, 40]; ylm = [-80, 150];
    ent_th = 1; dx_th = 1e-3; lc_gap = 4;

  otherwise
    error('Unknown dataset: %s', dataset);
end

% =========================
% Load graphs
% =========================
[~,~,ext] = fileparts(data_file);
if strcmpi(ext,'.mat')
  S = load(data_file); g = S.g;
elseif strcmpi(ext,'.g2o')
  g = read_graph_g2o(data_file);
elseif strcmpi(ext,'.graph')
  g = read_graph_toro(data_file);
else
  error('Unsupported dataset format: %s', ext);
end
gt = load(gt_file);                 % expects gt.g_current
g  = reorder_edges(g);              % GNSS never grows

% =========================
% Parameters & init
% =========================
total_edges     = numel(g.edges);
numIterations   = ceil(total_edges / batchSize);

err             = zeros(numIterations, 1);
gn              = zeros(numIterations, 1);
lc              = zeros(numIterations, 1);
di              = zeros(numIterations, 1);
ate_rmse        = nan(numIterations, 1);
flops_update    = [];
flops_solve     = [];
flops_update_inc= zeros(numIterations, 1);

g_current = struct('x', [], 'edges', [], 'idLookup', struct(), 'var2node', []);
current_batch_start = 1;

J = speye(3);
r = zeros(3, 1);
edge_J_index = [];             % numeric vector; grows in linearize_new
ifo = NaN; lo = 1;               % info-gain proxy bookkeeping

% Plot colors
color_SPO_IGG = 'g';
color_SPO_LCG = 'b';
color_SPO     = 'r';
color_IGG     = 'm';
color_LCG     = 'c';
color_none    = 'k';

% =========================
% Main loop
% =========================
tic
for it = 1:numIterations
  fprintf('Iteration %d\n', it);

  % ---- batch edges ----
  batch_end = min(current_batch_start + batchSize - 1, total_edges);
  edge_ids  = current_batch_start:batch_end;
  current_batch_start = batch_end + 1;

  % ---- update graph with new edges/nodes ----
  [g_current, edge_nodes, new_nodes, new_vars, loop_closure] = update_graph(edge_ids, g, g_current, lc_gap);
  lc(it) = loop_closure;

  % ---- linearize new edges & factorize ----
  [R, b, J, r, edge_J_index, p, parent] = linearize_new(J, r, edge_J_index, g_current, edge_ids, g);

  % flop estimate for new cols
  if ~isempty(new_vars)
    nnz_cols = full(sum(R~=0,1));
    upd0 = sum(nnz_cols(new_vars).^2);
  else
    upd0 = 0;
  endif
  flops_update = [flops_update; upd0]; %#ok<AGROW>
  flops_update_inc(it) = flops_update_inc(it) + upd0;

  % ---- info gain ----
  ifn = sum(log(abs(diag(R)))); ln = length(R);
  if isnan(ifo), ifo = ifn; lo = ln; end
  di(it) = ifn - ifo*(ln/lo); ifo = ifn; lo = ln;

  % ---- SPO vs Full GN branch -------------------------------------------
  if use_spo
    % ========== SPO: partial solver with gating ==========
    nx = numel(g_current.x);
    switch gating
      case 'none'
        c = color_SPO;
        affected_vars = 1:nx;
      case 'IGG'
        c = color_SPO_IGG;
        if di(it) > ent_th
          affected_vars = 1:nx;
        else
          affected_vars = unique(new_vars(:)).';
        end
      case 'LCG'
        c = color_SPO_LCG;
        if lc(it)
          affected_vars = 1:nx;
        else
          affected_vars = unique(new_vars(:)).';
        end
      otherwise
        error('Unknown gating: %s (use IGG | LCG | none)', gating);
    end

    % ---- GN refinement on affected vars ----
    if isempty(affected_vars)
      continue;                        % nothing to do this increment
    end

    dx = zeros(nx,1);
    for k = 1:max_GNi
      % solve on affected vars
      dx(affected_vars) = solve_affected(R, b, affected_vars, p, parent, .3);

      flops_solve = [flops_solve; 2*sum(nnz_cols(affected_vars))];

      % threshold: truly affected vars/edges/nodes
      [affected_vars, affected_edge_ids, ~] = find_affected(g_current, dx(affected_vars), dx_th, affected_vars);

      % stop if converged
      if isempty(affected_vars)
        break;
      end

      % state update
      g_current.x(affected_vars) = g_current.x(affected_vars) - dx(affected_vars);

      % relinearize affected edges & refactorize
      [R, b, J, r] = linearize_affected(J, r, edge_J_index, p, g_current, affected_edge_ids);

      nnz_cols = full(sum(R~=0,1));
      updK = min(2*sum(nnz_cols(affected_vars).^2), sum(nnz_cols.^2));
      flops_update = [flops_update; updK];
      flops_update_inc(it) = flops_update_inc(it) + updK;

      gn(it) = k;
    end

  else
    % ========== Full GN (GNi-like): always full solve ==========
    % decide whether to run the full update this increment
    switch gating
      case 'none'
        c = color_none;
        update = true;
      case 'IGG'
        c = color_IGG;
        if di(it) > ent_th
          update = true;
        else
          update = false;
        end
      case 'LCG'
        c = color_LCG;
        if lc(it)
          update = true;
        else
          update = false;
        end
      otherwise
        error('Unknown gating: %s (use IGG | LCG | none)', gating);
    end

    if update
      for k = 1:max_GNi
        % solve
        dx_p  = R \ (R' \ b);
        dx    = zeros(length(dx_p), 1);
        dx(p) = dx_p;

        flops_solve = [flops_solve; 2*nnz(R)];

        % stop if converged
        if max(abs(dx)) <= dx_th
          break;
        end

        % state update
        g_current.x = g_current.x - dx;

        % relinearize ALL edges and refactorize for next iteration
        all_edge_ids = 1:numel(g_current.edges);
        [R, b, J, r] = linearize_affected(J, r, edge_J_index, p, g_current, all_edge_ids);

        nnz_cols = full(sum(R~=0,1));
        updk = sum(nnz_cols.^2);
        flops_update = [flops_update; updk];
        flops_update_inc(it) = flops_update_inc(it) + updk;

        gn(it) = k;
      end
    end
  end
  % -----------------------------------------------------------------------

  % ---- ATE (aligned) vs ground truth ----
  [ate_rmse(it), ~] = compute_ate_rmse(g_current, gt.g_current);

  % ---- normalized error ----
  total_meas = sum(arrayfun(@(e) numel(e.measurement), g_current.edges(:)));
  err(it) = compute_global_error(g_current) / max(total_meas, 1);
end
fprintf('Elapsed time (s)  %f\n', toc);

% =========================
% plots & stats
% =========================

% normalized error
figure(20); hold on;
semilogy(err, c, 'LineWidth', 3);
xlabel('increment'); ylabel('normalized error');
legend(sprintf('final = %.3g', err(end)));

% Gauss–Newton iterations
figure(5); hold on;
bar(gn, c);
title('GN iterations');
xlabel('increment'); ylabel('#iterations');

% update FLOPs
figure(11); hold on;
bar(flops_update);
title('update FLOPs'); ylabel('\approx \Sigma (nnz_j)^2');
grid on;

% solve FLOPs
figure(12); hold on;
bar(flops_solve);
title('solve FLOPs'); ylabel('\approx \Sigma (nnz_j)^2');
grid on;

% ATE over increments
figure(41); hold on;
plot(ate_rmse, c, 'LineWidth', 3);
grid on;
xlabel('increment'); ylabel('ATE RMSE');

% info gain with loop-closure markers
idxLC = find(lc(:) > 0);
figure(33); hold on;
plot(1:numel(di), di, '-', 'Color', 'g', 'LineWidth', 3);
yl = ylim;                                % for vertical lines
for k = 1:numel(idxLC)
    line([idxLC(k) idxLC(k)], yl, 'Color', 'b', 'LineStyle', ':', 'LineWidth', 3);
end
grid on;
xlabel('increment'); ylabel('information gain');

% cumulative update FLOPs (per increment)
figure(54); hold on;
stairs(1:numIterations, cumsum(flops_update_inc), c, 'LineWidth', 3);
grid on;
xlabel('increment'); ylabel('cumulative update FLOPs');

% print summary stats
fprintf('\nFinal normalized error  %1.6g\n', err(end));
fprintf('Avg   normalized error  %1.6g\n', mean(err));
fprintf('Final ATE (RMSE)        %1.6g\n', ate_rmse(end));
fprintf('Avg   ATE (RMSE)        %1.6g\n', mean(ate_rmse));
fprintf('Avg   solve FLOPs       %1.6g\n', sum(flops_solve)/numIterations);
fprintf('Avg   update FLOPs      %1.8g\n', sum(flops_update)/numIterations);
