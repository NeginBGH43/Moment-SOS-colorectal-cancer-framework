%% momentsos_cvx.m
% Moment-SOS lower/upper bounds on P(tau <= T) for the killed multistage CTMC
% (APC -> KRAS -> TP53), solved with CVX, plus an exact baseline for small K.
%
% State of the alive process: x = (x0,x1,x2), x3 = 0, sum(x) <= K.
%   birth    b_i x_i (1 - N/K)^+      death  d_i x_i
%   mutation mu_i x_i  (i = 1,2 move to the next stage and stay alive)
%   killing  mu_3 x_3-stage (stage-2 -> malignant) at rate mu(3)*x2   (first passage)
%
% Requirements: MATLAB >= R2016b (local functions in scripts), CVX, and an SDP solver
% (SDPT3/SeDuMi are bundled with CVX; MOSEK is recommended for d >= 4).
% No Symbolic Toolbox needed.
%
% Sanity check: with the DEFAULT parameters below the exact value must be 0.1558242167.
% The constraint generator was verified in Python (exact_vs_momentsos.py, mirror of this
% algorithm): exact moments satisfy all Liouville rows to ~1e-10 and give PSD matrices.
% This MATLAB file itself has NOT been executed by the author of the file.
% FIX (v2): an earlier version imposed  yt_0 = T  (mass of the occupation measure). That is WRONG:
% the mass of the killed occupation measure is E[min(T,tau)] < T and follows from the Liouville
% rows. With it, d=3 was infeasible and d=2 results were not valid bounds. Re-run everything.

clear; clc;

%% ======================= USER PARAMETERS (edit here) =======================
beta   = [1.0 1.2 1.5];      % proliferation rates   beta_0..beta_2   (1/yr)
delta  = [0.5 0.5 0.5];      % elimination rates     delta_0..delta_2 (1/yr)
mu     = [0.30 0.25 0.20];   % mutation rates  mu_APC, mu_KRAS, mu_TP53 (1/yr)
K      = 8;                  % carrying capacity
xini   = [5 2 0];            % initial state (x0,x1,x2); x3 = 0 and sum(xini) <= K
Tlist  = 2.0;                % horizon(s) in years, e.g. 0.5:0.5:10 for a risk curve
dlist  = [2 3];              % relaxation orders to run
useScaling = true;           % solve in scaled variables  x/K, t/T  (recommended)
solverName = '';             % '' = CVX default;  'mosek', 'sdpt3', 'sedumi', ...
precisionName = 'default';   % 'low','medium','default','high','best'
Kexact_max = 60;             % compute the exact baseline only if K <= this
alpha_tol  = 0.05;           % risk threshold for T_safe (needs several T values)
%% ===========================================================================

assert(sum(xini) <= K, 'sum(xini) must be <= K');
results = struct('T',[],'d',[],'L',[],'U',[],'exact',[],'time',[]);
ir = 0;
for T = Tlist(:)'
    if K <= Kexact_max
        ex = exact_ptau(K, T, xini, beta, delta, mu);
    else
        ex = NaN;
    end
    for d = dlist
        t0 = tic;
        prob = build_problem(d, beta, delta, mu, K, T, xini, useScaling);
        [L, stL] = solve_bound(prob, +1, solverName, precisionName);
        [U, stU] = solve_bound(prob, -1, solverName, precisionName);
        tt = toc(t0);
        ir = ir + 1;
        results(ir) = struct('T',T,'d',d,'L',L,'U',U,'exact',ex,'time',tt);
        ok = '';
        if ~isnan(ex), ok = sprintf('  enclosed=%d', (L <= ex + 1e-6) && (ex <= U + 1e-6)); end
        fprintf('T=%6.3f d=%d  L=%.6f  U=%.6f  exact=%.6f  gap=%.2e  (%s/%s, %.1fs)%s\n', ...
                T, d, L, U, ex, U - L, stL, stU, tt, ok);
    end
end

% ---- conservative risk horizon  T_safe = inf{T : U_F(T) >= alpha_tol}  (largest d) ----
if numel(Tlist) > 1
    dmax = max(dlist);
    sel  = [results.d] == dmax;
    Ts = [results(sel).T]; Us = [results(sel).U];
    k = find(Us >= alpha_tol, 1);
    if isempty(k), fprintf('U_F stays below %.3f on the grid.\n', alpha_tol);
    elseif k == 1, fprintf('T_safe <= %.3f (already above threshold at first T).\n', Ts(1));
    else
        Tsafe = interp1(Us(k-1:k), Ts(k-1:k), alpha_tol);
        fprintf('T_safe(alpha=%.3f) ~ %.3f yr (linear interpolation, d=%d)\n', alpha_tol, Tsafe, dmax);
    end
end

%% =============================== FUNCTIONS =================================
function prob = build_problem(d, beta, delta, mu, K, T, xini, useScaling)
% Variables: yt = moments of the occupation measure in (t,x0,x1,x2), degree <= 2d
%            zt = moments of the terminal measure in (x0,x1,x2),     degree <= 2d
% With scaling: yt_{j,a} = int (t/T)^j (x/K)^a d mu_occ ,  zt_a = E[(x/K)^a ; tau > T]
    Y = monos(4, 2*d);  Z = monos(3, 2*d);
    ny = size(Y,1);  nz = size(Z,1);
    iy = make_index(Y, 4, 2*d);  iz = make_index(Z, 3, 2*d);
    if useScaling
        Dy = T.^Y(:,1) .* K.^sum(Y(:,2:4),2);   Dz = K.^sum(Z,2);
    else
        Dy = ones(ny,1);  Dz = ones(nz,1);
    end

    % ---------- Liouville equalities:  T^j z_a - sum c*y = [j==0] xini^a ----------
    ri = []; ci = []; vv = []; b = []; r = 0;
    for deg = 0:2*d-1
        for j = 0:deg
            for a0 = 0:deg-j
                for a1 = 0:deg-j-a0
                    a2 = deg - j - a0 - a1;  a = [a0 a1 a2];
                    P = A0_x(a, beta, delta, mu, K);          % rows [coef e0 e1 e2]
                    rowc = []; rowv = [];
                    for q = 1:size(P,1)
                        col = iy(lin_key([j P(q,2:4)], 4, 2*d));
                        rowc(end+1) = col;            %#ok<AGROW>
                        rowv(end+1) = -P(q,1) * Dy(col);  %#ok<AGROW>
                    end
                    if j > 0
                        col = iy(lin_key([j-1 a], 4, 2*d));
                        rowc(end+1) = col;            %#ok<AGROW>
                        rowv(end+1) = -j * Dy(col);   %#ok<AGROW>
                    end
                    colz = ny + iz(lin_key(a, 3, 2*d));
                    rowc(end+1) = colz;               %#ok<AGROW>
                    rowv(end+1) = T^j * Dz(colz - ny);%#ok<AGROW>
                    rhs = 0;  if j == 0, rhs = prod(xini .^ a); end
                    s = max(abs(rowv));  if s == 0, s = 1; end      % row normalisation
                    r = r + 1;
                    ri = [ri, r*ones(1,numel(rowc))]; %#ok<AGROW>
                    ci = [ci, rowc];                  %#ok<AGROW>
                    vv = [vv, rowv / s];              %#ok<AGROW>
                    b(r,1) = rhs / s;                 %#ok<AGROW>
                end
            end
        end
    end
    A = sparse(ri, ci, vv, r, ny + nz);

    % ---------- localizing / moment matrices ----------
    if useScaling, Kg = 1; Tg = 1; else, Kg = K; Tg = T; end
    gocc = { [1 0 0 0 0], ...                                  % moment matrix
             [1 0 1 0 0], [1 0 0 1 0], [1 0 0 0 1], ...        % x_i >= 0
             [Kg 0 0 0 0; -1 0 1 0 0; -1 0 0 1 0; -1 0 0 0 1], ...   % K - N >= 0
             [Kg^2 0 0 0 0; -1 0 2 0 0; -1 0 0 2 0; -1 0 0 0 2], ... % K^2 - |x|^2 >= 0
             [1 1 0 0 0], ...                                  % t >= 0
             [Tg 0 0 0 0; -1 1 0 0 0] };                       % T - t >= 0
    gter = { [1 0 0 0], ...
             [1 1 0 0], [1 0 1 0], [1 0 0 1], ...
             [Kg 0 0 0; -1 1 0 0; -1 0 1 0; -1 0 0 1], ...
             [Kg^2 0 0 0; -1 2 0 0; -1 0 2 0; -1 0 0 2] };
    prob.Socc = build_S(Y, iy, 4, 2*d, gocc, d);
    prob.Sterm = build_S(Z, iz, 3, 2*d, gter, d);

    prob.A = A;  prob.b = b;  prob.ny = ny;  prob.nz = nz;
    prob.iobj = iy(lin_key([0 0 0 1], 4, 2*d));          % moment of x2
    prob.objcoef = mu(3) * Dy(prob.iobj);                % P(tau<=T) = mu_TP53 * int x2 d mu_occ
    prob.T = T;
end

function [val, status] = solve_bound(prob, sense, solverName, precisionName)
% sense = +1 : minimise (lower bound);  sense = -1 : maximise (upper bound)
    ny = prob.ny; nz = prob.nz; A = prob.A; b = prob.b;
    cvx_begin quiet
        if ~isempty(solverName), cvx_solver(solverName); end
        cvx_precision(precisionName);
        variable v(ny + nz)
        yt = v(1:ny);  zt = v(ny+1:end);
        if sense > 0
            minimize( prob.objcoef * yt(prob.iobj) )
        else
            maximize( prob.objcoef * yt(prob.iobj) )
        end
        subject to
            A * v == b;
            for k = 1:numel(prob.Socc)
                n = prob.Socc{k}.n;
                M = reshape(prob.Socc{k}.S * yt, n, n);
                (M + M')/2 == semidefinite(n);
            end
            for k = 1:numel(prob.Sterm)
                n = prob.Sterm{k}.n;
                M = reshape(prob.Sterm{k}.S * zt, n, n);
                (M + M')/2 == semidefinite(n);
            end
    cvx_end
    val = cvx_optval;  status = cvx_status;
end

% ---------------------------- exact baseline --------------------------------
function Ptau = exact_ptau(K, T, xini, beta, delta, mu)
    n = nchoosek(K+3, 3);
    idx = zeros(K+1, K+1, K+1);  S = zeros(n, 3);  m = 0;
    for a = 0:K
        for b = 0:K-a
            for c = 0:K-a-b
                m = m + 1;  idx(a+1,b+1,c+1) = m;  S(m,:) = [a b c];
            end
        end
    end
    rr = []; cc = []; vv = []; out = zeros(n,1);
    for s = 1:n
        x = S(s,:);  N = sum(x);
        for i = 1:3
            if x(i) == 0, continue; end
            e = zeros(1,3); e(i) = 1;
            if N < K
                y = x + e;  rate = beta(i)*x(i)*(1 - N/K);
                rr(end+1) = idx(y(1)+1,y(2)+1,y(3)+1); cc(end+1) = s; vv(end+1) = rate; out(s) = out(s) + rate; %#ok<AGROW>
            end
            y = x - e;  rate = delta(i)*x(i);
            rr(end+1) = idx(y(1)+1,y(2)+1,y(3)+1); cc(end+1) = s; vv(end+1) = rate; out(s) = out(s) + rate; %#ok<AGROW>
            if i < 3
                y = x; y(i) = y(i) - 1; y(i+1) = y(i+1) + 1;  rate = mu(i)*x(i);
                rr(end+1) = idx(y(1)+1,y(2)+1,y(3)+1); cc(end+1) = s; vv(end+1) = rate; out(s) = out(s) + rate; %#ok<AGROW>
            end
        end
        out(s) = out(s) + mu(3)*x(3);                  % killing (first passage)
    end
    Q = sparse(rr, cc, vv, n, n) - spdiags(out, 0, n, n);
    p0 = zeros(n,1);  p0(idx(xini(1)+1, xini(2)+1, xini(3)+1)) = 1;
    opts = odeset('Jacobian', Q, 'RelTol', 1e-10, 'AbsTol', 1e-13);
    [~, Y] = ode15s(@(t,p) Q*p, [0 T], p0, opts);
    Ptau = 1 - sum(Y(end,:));
end

% ------------------------- polynomial helpers --------------------------------
% A polynomial is a matrix with rows [coef, e_1, ..., e_n].
function P = combine(P)
    if isempty(P), return; end
    [U, ~, ic] = unique(P(:,2:end), 'rows');
    c = accumarray(ic, P(:,1));
    P = [c, U];
    P(abs(P(:,1)) == 0, :) = [];
end

function C = pmul(A, B)
    n = size(A,1);  m = size(B,1);
    C = [kron(A(:,1), B(:,1)), repelem(A(:,2:end), m, 1) + repmat(B(:,2:end), n, 1)];
    C = combine(C);
end

function P = scal(P, c)          % scale COEFFICIENTS only (never the exponents)
    P(:,1) = c * P(:,1);
end

function P = shiftpoly(a, s)     % prod_k (x_k + s_k)^{a_k}, 3 variables
    P = [1 0 0 0];
    for k = 1:3
        rows = [];
        for m = 0:a(k)
            pw = a(k) - m;
            if pw > 0, c = nchoosek(a(k), m) * s(k)^pw; else, c = nchoosek(a(k), m); end
            if c ~= 0
                e = zeros(1,3); e(k) = m;  rows = [rows; c, e]; %#ok<AGROW>
            end
        end
        P = pmul(P, rows);
    end
end

function P = A0_x(a, beta, delta, mu, K)
% A0 x^a for the killed generator, as a polynomial in x0,x1,x2 (rows [coef e0 e1 e2])
    xa = [1 a];  res = zeros(0,4);
    N = [1 1 0 0; 1 0 1 0; 1 0 0 1];
    for i = 1:3
        ei = zeros(1,3); ei(i) = 1;  xi = [1 ei];
        xiN = pmul(xi, N);  xiN(:,1) = (-beta(i)/K) * xiN(:,1);
        rate = combine([scal(xi, beta(i)); xiN]);                 % beta_i x_i (1 - N/K)
        dif = combine([shiftpoly(a, ei); scal(xa, -1)]);
        if ~isempty(dif), res = [res; pmul(rate, dif)]; end       %#ok<AGROW>
        dif = combine([shiftpoly(a, -ei); scal(xa, -1)]);
        if ~isempty(dif), res = [res; pmul(scal(xi, delta(i)), dif)]; end %#ok<AGROW>
        if i < 3
            f = zeros(1,3); f(i) = -1; f(i+1) = 1;
            dif = combine([shiftpoly(a, f); scal(xa, -1)]);
            if ~isempty(dif), res = [res; pmul(scal(xi, mu(i)), dif)]; end %#ok<AGROW>
        end
    end
    x2 = [1 0 0 1];
    res = [res; pmul(scal(x2, -mu(3)), xa)];                      % killing term
    P = combine(res);
end

% -------------------------- indexing helpers ---------------------------------
function M = monos(nv, deg)      % all exponent vectors with sum <= deg, sorted by degree
    g = cell(1, nv);  [g{:}] = ndgrid(0:deg);
    M = zeros(numel(g{1}), nv);
    for k = 1:nv, M(:,k) = g{k}(:); end
    M = M(sum(M,2) <= deg, :);
    [~, o] = sortrows([sum(M,2), M]);  M = M(o,:);
end

function key = lin_key(E, nv, maxdeg)   % linear key of exponent rows
    base = maxdeg + 1;
    key = 1 + E * (base .^ (0:nv-1))';
end

function idx = make_index(M, nv, maxdeg)  % idx(key) = position of the monomial in M
    idx = zeros((maxdeg+1)^nv, 1);
    idx(lin_key(M, nv, maxdeg)) = 1:size(M,1);
end

function out = build_S(Mon, idx, nv, maxdeg, glist, d)
% For each polynomial g (rows [coef exps]) build S with  vec(Mloc) = S*y,  Mloc(u,w) = sum_e c y_{u+w+e}
    out = cell(1, numel(glist));
    for q = 1:numel(glist)
        g = glist{q};
        gdeg = max(sum(g(:,2:end), 2));
        dd = d - ceil(gdeg/2);
        B = monos(nv, dd);  n = size(B,1);
        [I, J] = ndgrid(1:n, 1:n);  I = I(:);  J = J(:);
        Esum = B(I,:) + B(J,:);
        rows = I + n*(J-1);
        rr = []; cc = []; vv = [];
        for t = 1:size(g,1)
            cols = idx(lin_key(Esum + g(t,2:end), nv, maxdeg));
            rr = [rr; rows];              %#ok<AGROW>
            cc = [cc; cols];              %#ok<AGROW>
            vv = [vv; g(t,1)*ones(n*n,1)];%#ok<AGROW>
        end
        out{q}.S = sparse(rr, cc, vv, n*n, size(Mon,1));
        out{q}.n = n;
    end
end