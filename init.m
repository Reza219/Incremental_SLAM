% INIT  One-time environment setup for the SLAM repo.
% - Adds the tools/ folder to the path
% - Sets a reproducible RNG seed
% - Prints MATLAB/Octave version info
% - Applies a few sane plotting defaults

function init(seed)
  if nargin < 1 || isempty(seed), seed = 123; end

  % ---- Paths -------------------------------------------------------------
  here = fileparts(mfilename('fullpath'));
  addpath(genpath(fullfile(here, 'tools')));

  % ---- Version / RNG -----------------------------------------------------
  if exist('OCTAVE_VERSION','builtin')
    fprintf('[INIT] GNU Octave %s\n', OCTAVE_VERSION);
    rand('state', seed); randn('state', seed);
  else
    v = ver('MATLAB');
    fprintf('[INIT] MATLAB %s\n', v.Version);
    rng(seed, 'twister');
  end

  % ---- Console / formatting ---------------------------------------------
  format short g;

  % ---- Plot defaults (light, non-invasive) ------------------------------
  set(0, 'DefaultFigureColor', 'w');
  set(0, 'DefaultAxesFontName', 'Helvetica');
  set(0, 'DefaultAxesFontSize', 11);
  set(0, 'DefaultLineLineWidth', 1.25);
  set(0, 'DefaultAxesBox', 'on');

  % ---- Sanity checks -----------------------------------------------------
  if ~exist(fullfile(here, 'data'), 'dir')
    warning('[INIT] data/ folder not found. Create data/ and place datasets inside.');
  end

  if ~exist(fullfile(here, 'tools'), 'dir')
    warning('[INIT] tools/ folder not found. Make sure repository is complete.');
  end

  fprintf('[INIT] Path set. Seed=%d. Ready.\n', seed);
end
