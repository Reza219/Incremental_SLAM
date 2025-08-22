function Fx = compute_global_error(g)
% COMPUTE_GLOBAL_ERROR  Total weighted squared error over all constraints.
%
%   Fx = sum_e e' * W * e, for edge types:
%     'P' : pose–pose
%     'L' : pose–landmark (z in robot frame)
%     'G' : GNSS/prior on pose (position-only)
%
% Inputs
%   g.x        (nx1)              state vector
%   g.edges    (struct array)     .type, .fromIdx, .toIdx, .measurement, .information
%
% Notes
%   - Uses invt() to invert SE(2) transforms.
%   - Information matrices are symmetrized defensively.

  Fx = 0;
  if ~isfield(g,'edges') || isempty(g.edges), return; end
  if ~isfield(g,'x') || isempty(g.x), error('compute_global_error:EmptyState','g.x is empty.'); end

  for eid = 1:numel(g.edges)
    e = g.edges(eid);
    W = (e.information + e.information.')/2;  % ensure symmetry

    switch e.type
      case 'P'  % pose–pose
        xi = g.x(e.fromIdx : e.fromIdx+2);
        xj = g.x(e.toIdx   : e.toIdx+2);

        Xi = v2t(xi);
        Xj = v2t(xj);
        Z  = v2t(e.measurement);

        e_ij = t2v(invt(Z) * (invt(Xi) * Xj));   % use invt()
        Fx   = Fx + e_ij' * W * e_ij;

      case 'L'  % pose–landmark
        x = g.x(e.fromIdx : e.fromIdx+2);   % [x;y;th]
        l = g.x(e.toIdx   : e.toIdx+1);     % [x;y]

        c = cos(x(3)); s = sin(x(3));
        Rt = [ c,  s;    % world->robot rotation (R')
              -s,  c];
        e_il = Rt * (l - x(1:2)) - e.measurement(:);

        Fx = Fx + e_il' * W * e_il;

      case 'G'  % GNSS/prior (unary, position only)
        xi = g.x(e.fromIdx : e.fromIdx+2);
        z  = e.measurement(:);
        if numel(z) ~= 2
          error('compute_global_error:BadGNSS','GNSS measurement must be length 2.');
        end
        e_g = [xi(1) - z(1); xi(2) - z(2)];
        Fx  = Fx + e_g' * W * e_g;

      otherwise
        % ignore unknown types
    end
  end
end
