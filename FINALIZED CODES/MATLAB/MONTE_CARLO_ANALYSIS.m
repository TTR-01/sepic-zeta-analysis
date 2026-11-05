%% Monte Carlo efficiency comparison: SEPIC vs ZETA, all nine cases.
% Draws N_samples Latin-hypercube perturbations of Rds_on, Vd, DCR, ESR,
% L, C, Rg around their nominal values, runs the electro-thermal solver
% for both topologies at each of the nine (Vin, Io) points, and exports
% histogram and summary figures in IEEE and report styles.
%
% Needs dpt_lookup.csv in the working folder.

clc; clear; close all;

%% Nominal values and tolerances --------------------------------------------

p_nom.Rds_on_25 = 25e-3;
p_nom.Vd_25     = 0.5;
p_nom.DCR       = 50e-3;
p_nom.ESR       = 10e-3;
p_nom.L         = 180e-6;
p_nom.C         = 15e-6;
p_nom.Rg        = 3.4;

% +/- fraction of nominal
tol.Rds_on = 0.20; tol.Vd  = 0.10; tol.DCR = 0.10;
tol.ESR    = 0.20; tol.L   = 0.10; tol.C   = 0.10;
tol.Rg     = 0.05;

%% Sampling settings --------------------------------------------------------

N_samples = 1000;
Vin_vec = [9, 12, 15];
Io_vec  = [0.5, 1.0, 1.5];
N_cases = 9;

%% Latin hypercube ---------------------------------------------------------
% Fixed seed so every figure and table reproduces. Set to [] for fresh runs.

rng_seed = 1;
if ~isempty(rng_seed), rng(rng_seed, 'twister'); end
lhs_raw = lhsdesign(N_samples, 7);
samples = zeros(N_samples, 7);
samples(:,1) = (lhs_raw(:,1) * 2 - 1) * tol.Rds_on;
samples(:,2) = (lhs_raw(:,2) * 2 - 1) * tol.Vd;
samples(:,3) = (lhs_raw(:,3) * 2 - 1) * tol.DCR;
samples(:,4) = (lhs_raw(:,4) * 2 - 1) * tol.ESR;
samples(:,5) = (lhs_raw(:,5) * 2 - 1) * tol.L;
samples(:,6) = (lhs_raw(:,6) * 2 - 1) * tol.C;
samples(:,7) = (lhs_raw(:,7) * 2 - 1) * tol.Rg;

eff_SEPIC_all   = zeros(N_samples, N_cases);
eff_ZETA_all    = zeros(N_samples, N_cases);
prob_sepic_wins = zeros(N_cases, 1);
mean_SEPIC = zeros(N_cases, 1);
mean_ZETA  = zeros(N_cases, 1);
std_SEPIC  = zeros(N_cases, 1);
std_ZETA   = zeros(N_cases, 1);

fprintf('Running Monte Carlo (%d samples x %d cases)...\n', N_samples, N_cases);
tic;

%% Sweep ----------------------------------------------------------------
% Both topologies see the same perturbed parameter set per sample, so the
% win probability below is a fair paired comparison, not two independent
% distributions that happen to overlap.

for case_idx = 1:N_cases
    [V_idx, I_idx] = ind2sub([3, 3], case_idx);
    Vin = Vin_vec(V_idx);
    Io  = Io_vec(I_idx);

    fprintf('  Case %d/9: Vin=%dV, Io=%.1fA...\n', case_idx, Vin, Io);

    for k = 1:N_samples
        p = p_nom;
        p.Rds_on_25 = p_nom.Rds_on_25 * (1 + samples(k,1));
        p.Vd_25     = p_nom.Vd_25     * (1 + samples(k,2));
        p.DCR       = p_nom.DCR       * (1 + samples(k,3));
        p.ESR       = p_nom.ESR       * (1 + samples(k,4));
        p.L         = p_nom.L         * (1 + samples(k,5));
        p.C         = p_nom.C         * (1 + samples(k,6));
        p.Rg        = p_nom.Rg        * (1 + samples(k,7));

        res_SEPIC = solve_point(Vin, Io, p, 'SEPIC');
        res_ZETA  = solve_point(Vin, Io, p, 'ZETA');

        eff_SEPIC_all(k, case_idx) = res_SEPIC.eff_system;
        eff_ZETA_all(k, case_idx)  = res_ZETA.eff_system;
    end

    mean_SEPIC(case_idx) = mean(eff_SEPIC_all(:, case_idx));
    mean_ZETA(case_idx)  = mean(eff_ZETA_all(:, case_idx));
    std_SEPIC(case_idx)  = std(eff_SEPIC_all(:, case_idx));
    std_ZETA(case_idx)   = std(eff_ZETA_all(:, case_idx));
    prob_sepic_wins(case_idx) = mean(eff_SEPIC_all(:, case_idx) > eff_ZETA_all(:, case_idx)) * 100;
end
toc;

%% Console summary ----------------------------------------------------------

fprintf('\n========================================================================\n');
fprintf('MONTE CARLO SUMMARY (All 9 Cases)\n');
fprintf('========================================================================\n');
fprintf('%-6s | %-4s | %-4s | %-12s | %-12s | %-10s\n', ...
    'Case', 'Vin', 'Io', 'SEPIC Eff (%)', 'ZETA Eff (%)', 'P(SEPIC>ZETA)');
fprintf('%s\n', repmat('-', 1, 70));
for case_idx = 1:N_cases
    [V_idx, I_idx] = ind2sub([3, 3], case_idx);
    fprintf('%-6d | %-4d | %-4.1f | %.2f ± %.2f | %.2f ± %.2f | %.1f%%\n', ...
        case_idx, Vin_vec(V_idx), Io_vec(I_idx), ...
        mean_SEPIC(case_idx), std_SEPIC(case_idx), ...
        mean_ZETA(case_idx), std_ZETA(case_idx), ...
        prob_sepic_wins(case_idx));
end
fprintf('========================================================================\n');

%% Figures ------------------------------------------------------------------

outDir    = 'report_figures';
styleList = {'ieee', 'report'};

% Cases are ordered Vin-major, so row boundaries in the 3x3 grid and the
% group separators in the summary bar chart both fall where Vin changes.
caseVin = zeros(N_cases, 1);  caseIo = zeros(N_cases, 1);
for cs = 1:N_cases
    [V_idx, I_idx] = ind2sub([3, 3], cs);
    caseVin(cs) = Vin_vec(V_idx);
    caseIo(cs)  = Io_vec(I_idx);
end
caseLabels = arrayfun(@(v, i) sprintf('%dV/%.1fA', v, i), caseVin, caseIo, 'UniformOutput', false);
sepAt = find(diff(caseVin) ~= 0) + 0.5;

clr_sepic = [0.000 0.447 0.698];     % Okabe-Ito blue
clr_zeta  = [0.835 0.369 0.000];     % Okabe-Ito vermillion

% Both topologies share the same bin edges within a case so the two
% histograms overlay honestly, and a single common y-range across all
% nine tiles keeps the visual scale comparable panel to panel.
nBins    = 25;
edges    = cell(N_cases, 1);
maxCount = 0;
for cs = 1:N_cases
    both = [eff_SEPIC_all(:, cs); eff_ZETA_all(:, cs)];
    lo = min(both);  hi = max(both);
    if hi - lo < 1e-6, hi = lo + 1e-3; end
    edges{cs} = linspace(lo, hi, nBins + 1);
    maxCount  = max([maxCount, ...
        max(histcounts(eff_SEPIC_all(:, cs), edges{cs})), ...
        max(histcounts(eff_ZETA_all(:, cs),  edges{cs}))]);
end

for sIdx = 1:numel(styleList)
close all;
ieee_cfg(styleList{sIdx});
apply_ieee_defaults();
IEEE   = ieee_cfg();
figDir = fullfile(outDir, styleList{sIdx});
fprintf('\n--- Exporting "%s" figure set -> %s ---\n', styleList{sIdx}, figDir);

% Fig 11 - 3x3 histogram grid, rows = Vin, columns = Io
fig11 = new_fig('MonteCarlo_All9Cases', IEEE.col2, 6.0);
tl11  = tiledlayout(fig11, 3, 3, 'TileSpacing','compact', 'Padding','compact');
for cs = 1:N_cases
    ax = nexttile(tl11); hold(ax, 'on');
    span = edges{cs}(end) - edges{cs}(1);

    % Filled histograms carry the legend entries; the outlined stairs on
    % top give a second encoding (solid vs dashed) so the plot survives
    % greyscale printing.
    histogram(ax, eff_SEPIC_all(:, cs), 'BinEdges', edges{cs}, 'FaceColor', clr_sepic, ...
        'FaceAlpha', 0.55, 'EdgeColor', 'w', 'LineWidth', 0.3, 'DisplayName', 'SEPIC');
    histogram(ax, eff_ZETA_all(:, cs),  'BinEdges', edges{cs}, 'FaceColor', clr_zeta, ...
        'FaceAlpha', 0.55, 'EdgeColor', 'w', 'LineWidth', 0.3, 'DisplayName', 'ZETA');
    histogram(ax, eff_SEPIC_all(:, cs), 'BinEdges', edges{cs}, 'DisplayStyle', 'stairs', ...
        'EdgeColor', clr_sepic, 'LineStyle', '-',  'LineWidth', 0.9, 'HandleVisibility', 'off');
    histogram(ax, eff_ZETA_all(:, cs),  'BinEdges', edges{cs}, 'DisplayStyle', 'stairs', ...
        'EdgeColor', clr_zeta,  'LineStyle', '--', 'LineWidth', 0.9, 'HandleVisibility', 'off');

    % Mean of each distribution, styled to match the corresponding outline
    xline(ax, mean_SEPIC(cs), '--', 'Color', clr_sepic, 'LineWidth', IEEE.lw, 'HandleVisibility', 'off');
    xline(ax, mean_ZETA(cs),  ':',  'Color', clr_zeta,  'LineWidth', IEEE.lw, 'HandleVisibility', 'off');

    ylim(ax, [0, 1.5 * maxCount]);
    xlim(ax, [edges{cs}(1) - 0.03*span, edges{cs}(end) + 0.03*span]);
    xtickformat(ax, '%.1f');
    if mod(cs - 1, 3) ~= 0, set(ax, 'YTickLabel', []); end
    grid(ax, 'on');

    % Panel tag and win probability, both top-left
    panel_tag(ax, sprintf('(%s) %d V, %.1f A', char('a' + cs - 1), caseVin(cs), caseIo(cs)));
    text(ax, 0.02, 0.84, sprintf('SEPIC wins %.0f%%', prob_sepic_wins(cs)), 'Units', 'normalized', ...
        'FontSize', IEEE.fontSmall, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    if cs == 1, axFirst = ax; end
end
xlabel(tl11, lx('System efficiency \eta_{sys} (%)'), 'Interpreter', IEEE.interp);
ylabel(tl11, 'Count', 'Interpreter', IEEE.interp);
lgd = legend(axFirst, 'Orientation', 'horizontal', 'FontSize', IEEE.fontSmall);
lgd.Layout.Tile = 'north';
ttl(tl11, sprintf('Monte Carlo efficiency distribution (N = %d samples)', N_samples));
save_fig(fig11, figDir, 'Fig11_MonteCarlo_All9Cases');

% Fig 12 - mean +/- 1 sigma and win probability per operating point
fig12 = new_fig('MonteCarlo_Summary', IEEE.col2, 2.8);
tl12  = tiledlayout(fig12, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl12); hold(ax, 'on');
xs = (1:N_cases)';
% Offset the two markers horizontally so overlapping error bars stay
% distinguishable at this column width.
errorbar(ax, xs - 0.12, mean_SEPIC, std_SEPIC, 'o', 'Color', clr_sepic, 'MarkerFaceColor', clr_sepic, ...
    'CapSize', 3, 'LineStyle', 'none', 'DisplayName', 'SEPIC');
errorbar(ax, xs + 0.12, mean_ZETA,  std_ZETA,  's', 'Color', clr_zeta,  'MarkerFaceColor', clr_zeta, ...
    'CapSize', 3, 'LineStyle', 'none', 'DisplayName', 'ZETA');
yLo = min([mean_SEPIC - std_SEPIC; mean_ZETA - std_ZETA]);
yHi = max([mean_SEPIC + std_SEPIC; mean_ZETA + std_ZETA]);
ySp = max(yHi - yLo, 1e-3);
ylim(ax, [yLo - 0.10*ySp, yHi + 0.35*ySp]);
ylabel(ax, lx('System efficiency \eta_{sys} (%)'));
legend(ax, 'Location', 'north', 'Orientation', 'horizontal', 'FontSize', IEEE.fontSmall);
ax_a = ax;

ax = nexttile(tl12); hold(ax, 'on');
bar(ax, xs, prob_sepic_wins, 0.6, 'FaceColor', [0.45 0.45 0.45], 'EdgeColor', 'none');
for k = 1:N_cases
    text(ax, k, prob_sepic_wins(k) + 2, lxf('%.0f%%', prob_sepic_wins(k)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', IEEE.fontTiny);
end
% 50 % is the natural "coin-flip" line: bars above it mean SEPIC wins
% more often than not at that operating point.
yline(ax, 50, '--', lx('50% (tie)'), 'Color', [0.75 0 0], 'LineWidth', 1.0, ...
    'Interpreter', IEEE.interp, 'LabelHorizontalAlignment', 'left', ...
    'LabelVerticalAlignment', 'bottom', 'FontSize', IEEE.fontSmall, 'HandleVisibility', 'off');
ylim(ax, [0, 118]);
ylabel(ax, lx('SEPIC win probability (%)'));
ax_b = ax;

for ax = [ax_a ax_b]
    set(ax, 'XTick', 1:N_cases, 'XTickLabel', caseLabels, 'XTickLabelRotation', 45, ...
        'TickLabelInterpreter', 'none');
    xlim(ax, [0.4, N_cases + 0.6]);
    yl = ylim(ax);
    for s = sepAt'
        plot(ax, [s s], yl, ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.6, 'HandleVisibility', 'off');
    end
    ylim(ax, yl);
    xlabel(ax, lx('Operating point (V_{in} / I_o)'));
    grid(ax, 'on'); ax.XGrid = 'off';
end
panel_tag(ax_a, '(a)');  panel_tag(ax_b, '(b)');
ttl(tl12, lx('Monte Carlo summary: mean \pm 1\sigma and win probability'));
save_fig(fig12, figDir, 'Fig12_MonteCarlo_Summary');
end   % style loop

fprintf('\nFigures exported to ./%s/<style>/ (ieee: PDF+EPS+600 dpi PNG; report: PDF+300 dpi PNG).\n', outDir);


%% ========================================================================
%  Local functions
%% ========================================================================

function res = solve_point(Vin, Io, p, topology)
% Electro-thermal fixed point for one operating point. Loops on eta and
% the two junction temperatures until all three settle, with
% under-relaxation so the coupled loop doesn't oscillate.
    eta_est = 0.90; Tj_M_est = 40; Tj_D_est = 40;
    relax_eta = 0.7; relax_T = 0.5;

    cfg.Vo = 12; cfg.fsw = 340e3;
    cfg.Rds_on_25 = p.Rds_on_25; cfg.Rds_on = p.Rds_on_25; cfg.alpha_Rds = 0.0055;
    cfg.Coss = 400e-12; cfg.Qgd = 25e-9; cfg.Qg = 48e-9;
    cfg.Vgate = 5; cfg.Vgs3 = 3.0; cfg.Rg = p.Rg;
    cfg.t_sw_rise = 84e-9; cfg.t_sw_fall = 15e-9;
    cfg.switchingLossMethod = 'dpt';
    cfg.Vds_max = 55; cfg.Id_max = 47;
    cfg.RJA_Q1 = 62; cfg.Tj_max_Q1 = 175;
    cfg.Cj_diode = 250e-12; cfg.Vrrm_diode = 40; cfg.Ifavg_diode = 3.0;
    cfg.RJA_D1 = 55; cfg.Tj_max_D1 = 150; cfg.T_amb = 40;
    cfg.kE_on = 0.0015; cfg.kE_off = 0.0010;
    cfg.temp_coeff_diode = -0.002; cfg.Vd_25 = p.Vd_25; cfg.Vd = p.Vd_25;
    cfg.DCR_L1 = p.DCR; cfg.DCR_L2 = p.DCR;
    cfg.ripple_factor = 0.3; cfg.rippleMethod = 'simplified';

    cfg.dVC1_target = 0.24; cfg.dVc_target = 0.24; cfg.dVo_target = 0.24; cfg.dVin_target = 0.15;

    % Field names differ between the two solvers only in the cap
    % identifiers; everything else is common.
    if strcmp(topology, 'SEPIC')
        cfg.C1_val = p.C; cfg.C2_val = p.C; cfg.Cin_val = p.C;
        cfg.ESR_C1 = p.ESR; cfg.ESR_C2 = p.ESR; cfg.ESR_Cin = p.ESR;
    else
        cfg.Cc_val = p.C; cfg.Cout_val = p.C; cfg.Cin_val = p.C;
        cfg.ESR_Cc = p.ESR; cfg.ESR_Cout = p.ESR; cfg.ESR_Cin = p.ESR;
        cfg.useCoupled = false; cfg.k_c = 0;
    end

    for iter = 1:20
        cfg.Rds_on   = cfg.Rds_on_25 * (1 + cfg.alpha_Rds * (Tj_M_est - 25));
        cfg.Vd       = cfg.Vd_25 + cfg.temp_coeff_diode * (Tj_D_est - 25);
        cfg.Tj_M_est = Tj_M_est; cfg.Tj_D_est = Tj_D_est;

        if strcmp(topology, 'SEPIC')
            res = calcSEPIC(Vin, Io, eta_est, cfg, p.L, p.L);
        else
            res = calcZETA(Vin, Io, eta_est, cfg, p.L, p.L);
        end

        eta_new  = res.eff_system / 100; Tj_M_new = res.Tj_Q1; Tj_D_new = res.Tj_D1;
        delta_eta = abs(eta_new - eta_est); delta_Tq = abs(Tj_M_new - Tj_M_est); delta_Td = abs(Tj_D_new - Tj_D_est);

        eta_est  = relax_eta * eta_new  + (1 - relax_eta) * eta_est;
        Tj_M_est = relax_T   * Tj_M_new + (1 - relax_T)   * Tj_M_est;
        Tj_D_est = relax_T   * Tj_D_new + (1 - relax_T)   * Tj_D_est;

        if delta_eta < 1e-4 && delta_Tq < 0.1 && delta_Td < 0.1, break; end
    end
end


function res = calcSEPIC(Vin, Io, eta_est, cfg, L1_in, L2_in)
% SEPIC steady-state at one point. Same model as the standalone SEPIC
% script, inlined here so the Monte Carlo loop has no external deps.
% Only eff_system and the two junction temperatures are consumed by the
% caller, so the return struct is trimmed to those three fields.
    Vo = cfg.Vo; fsw = cfg.fsw; Vd = cfg.Vd; Rds_on = cfg.Rds_on;
    Vgate = cfg.Vgate; Vgs3 = cfg.Vgs3; Rg = cfg.Rg; Qgd = cfg.Qgd; Qg = cfg.Qg;
    t_sw_rise = cfg.t_sw_rise; t_sw_fall = cfg.t_sw_fall; Vds_max = cfg.Vds_max; Id_max = cfg.Id_max;
    RJA_Q1 = cfg.RJA_Q1; Tj_max_Q1 = cfg.Tj_max_Q1; RJA_D1 = cfg.RJA_D1; Tj_max_D1 = cfg.Tj_max_D1; T_amb = cfg.T_amb;
    DCR1 = cfg.DCR_L1; DCR2 = cfg.DCR_L2; C1_val = cfg.C1_val; C2_val = cfg.C2_val; Cin_val = cfg.Cin_val;
    ESR_C1 = cfg.ESR_C1; ESR_C2 = cfg.ESR_C2; ESR_Cin = cfg.ESR_Cin;
    dVC1_target = cfg.dVC1_target; dVo_target = cfg.dVo_target; dVin_target = cfg.dVin_target;
    ripple_factor = cfg.ripple_factor; rippleMethod = cfg.rippleMethod; kE_on = cfg.kE_on; kE_off = cfg.kE_off;
    if isfield(cfg, 'Tj_M_est'), Tj_M_est = cfg.Tj_M_est; else, Tj_M_est = T_amb; end
    Pout = Vo * Io; Iin_avg = Io * Vo / (Vin * eta_est);
    relax_D = 0.5; D = (Vo + Vd) / (Vin + Vo + Vd);
    for iter = 1:30
        Iin_est = Io * Vo / (Vin * eta_est); Isw_avg = Iin_est + Io;
        D_new = (Vo + Vd + Io*DCR2) / (Vin + Vo + Vd - Iin_est*DCR1 - Isw_avg*Rds_on);
        D_next = relax_D * D_new + (1 - relax_D) * D;
        if abs(D_next - D) < 1e-7, D = D_next; break; end
        D = D_next;
    end
    L_req  = (Vin * D) / (Iin_avg * ripple_factor * fsw); L2_req = (Vin * D) / (Io * ripple_factor * fsw);
    L1_use = L_req; if L1_in > 0, L1_use = L1_in; end; L2_use = L2_req; if L2_in > 0, L2_use = L2_in; end
    delta_IL1 = Vin * D / (L1_use * fsw); delta_IL2 = Vin * D / (L2_use * fsw);
    IL1_peak = Iin_avg + delta_IL1/2; IL2_peak = Io + delta_IL2/2; IL1_min = Iin_avg - delta_IL1/2; IL2_min = Io - delta_IL2/2;
    IL1_rms = sqrt(Iin_avg^2 + delta_IL1^2/12); IL2_rms = sqrt(Io^2 + delta_IL2^2/12);
    Isw_base = Iin_avg + Io; Isw_peak = Isw_base + (delta_IL1 + delta_IL2)/2; Isw_valley = max(Isw_base - (delta_IL1 + delta_IL2)/2, 0);
    Isw_rms = sqrt(D) * sqrt(Isw_base^2 + (delta_IL1 + delta_IL2)^2/12);
    Vsw_theory = Vin + Vo + Vd; Isw_step = Iin_avg + Io;
    dVC1_actual = (Io * D) / (C1_val * fsw) + (IL1_peak + IL2_peak) * ESR_C1;
    if strcmp(rippleMethod, 'simplified'), dVo_cap = (Io * D) / (C2_val * fsw); else, dVo_cap = (Vin * D) / (8 * L2_use * fsw^2 * C2_val); end
    dVo_actual = dVo_cap + (IL1_peak + IL2_peak) * ESR_C2; dVin_actual = delta_IL1 / (8 * Cin_val * fsw) + delta_IL1 * ESR_Cin;
    ICin_rms = delta_IL1 / sqrt(12); IC1_rms = sqrt(D*(Io^2 + delta_IL2^2/12) + (1-D)*(Iin_avg^2 + delta_IL1^2/12));
    I_L2_off_avg = Io * D / (1 - D); IC2_rms = sqrt((1-D)*(I_L2_off_avg^2 + (delta_IL1+delta_IL2)^2/12) + D*Io^2);
    Pmos_cond = Isw_rms^2 * Rds_on; Pmos_sw_overlap_transition = 0.5 * Vsw_theory * (Isw_valley + Isw_peak) * (t_sw_rise + t_sw_fall) * fsw;
    Ig_on = max((Vgate - Vgs3)/Rg, 0.01); Ig_off = max(Vgs3/Rg, 0.01);
    t_miller_on = Qgd / Ig_on; t_miller_off = Qgd / Ig_off; Pmos_sw_overlap_miller = 0.5 * Vsw_theory * Isw_peak * (t_miller_on + t_miller_off) * fsw;
    dpt = dpt_params();
    if numel(dpt.V_ref) > 1
        E_on_ref_v = interp1(dpt.V_ref, dpt.E_on_ref, Vsw_theory, 'linear','extrap'); E_off_ref_v = interp1(dpt.V_ref, dpt.E_off_ref, Vsw_theory, 'linear','extrap');
        I_on_ref_v = interp1(dpt.V_ref, dpt.I_on_ref, Vsw_theory, 'linear','extrap'); I_off_ref_v = interp1(dpt.V_ref, dpt.I_off_ref, Vsw_theory, 'linear','extrap');
    else
        E_on_ref_v = dpt.E_on_ref * (Vsw_theory / dpt.V_ref); E_off_ref_v = dpt.E_off_ref * (Vsw_theory / dpt.V_ref);
        I_on_ref_v = dpt.I_on_ref(1); I_off_ref_v = dpt.I_off_ref(1);
    end
    I_on_ref_v = max(I_on_ref_v, 0.01); I_off_ref_v = max(I_off_ref_v, 0.01);
    E_on = E_on_ref_v * (Isw_valley / I_on_ref_v); E_off = E_off_ref_v * (Isw_peak / I_off_ref_v);
    E_on = E_on * (1 + kE_on * (Tj_M_est - 25)); E_off = E_off * (1 + kE_off * (Tj_M_est - 25));
    Pmos_sw_overlap_dpt = (E_on + E_off) * fsw; Pmos_sw_gate = Vgate * Qg * fsw;
    Pmos_sw = Pmos_sw_overlap_dpt + Pmos_sw_gate + 0; Pmos_total = Pmos_cond + Pmos_sw;
    Pdiode_cond = Vd * Io; Pdiode_cap = 0.5 * cfg.Cj_diode * Vsw_theory^2 * fsw; Pdiode = Pdiode_cond + Pdiode_cap;
    PL1_copper = IL1_rms^2 * DCR1; PL2_copper = IL2_rms^2 * DCR2; Pind_total = PL1_copper + PL2_copper;
    Pcap_C1 = IC1_rms^2 * ESR_C1; Pcap_C2 = IC2_rms^2 * ESR_C2; Pcap_Cin = ICin_rms^2 * ESR_Cin; Pcap_total = Pcap_C1 + Pcap_C2 + Pcap_Cin;
    Pmos_stage = Pmos_cond + Pmos_sw_overlap_dpt + 0; Ploss_stage = Pmos_stage + Pdiode + Pind_total + Pcap_total;
    Ploss_system = Ploss_stage + Pmos_sw_gate; Pin_system_est = Pout + Ploss_system;
    eff_system = 100 * Pout / Pin_system_est; res.eff_system = eff_system;
    res.Tj_Q1 = T_amb + Pmos_stage * RJA_Q1; res.Tj_D1 = T_amb + Pdiode * RJA_D1;
end


function res = calcZETA(Vin, Io, eta_est, cfg, L1_in, L2_in)
% ZETA steady-state at one point. Same as the standalone ZETA script.
% Only the three fields the caller needs are filled in.
    Vo = cfg.Vo; fsw = cfg.fsw; Vd = cfg.Vd; Rds_on = cfg.Rds_on;
    ripple_factor = cfg.ripple_factor; useCoupled = cfg.useCoupled; k_c = cfg.k_c;
    Vgate = cfg.Vgate; Vgs3 = cfg.Vgs3; Rg = cfg.Rg; Qgd = cfg.Qgd; Qg = cfg.Qg;
    t_sw_rise = cfg.t_sw_rise; t_sw_fall = cfg.t_sw_fall; Vds_max = cfg.Vds_max; Id_max = cfg.Id_max;
    RJA_Q1 = cfg.RJA_Q1; Tj_max_Q1 = cfg.Tj_max_Q1; RJA_D1 = cfg.RJA_D1; Tj_max_D1 = cfg.Tj_max_D1; T_amb = cfg.T_amb;
    DCR1 = cfg.DCR_L1; DCR2 = cfg.DCR_L2; Cc_val = cfg.Cc_val; Cout_val = cfg.Cout_val; Cin_val = cfg.Cin_val;
    ESR_Cc = cfg.ESR_Cc; ESR_Cout = cfg.ESR_Cout; ESR_Cin = cfg.ESR_Cin;
    dVc_target = cfg.dVc_target; dVo_target = cfg.dVo_target; dVin_target = cfg.dVin_target;
    kE_on = cfg.kE_on; kE_off = cfg.kE_off;
    if isfield(cfg, 'Tj_M_est'), Tj_M_est = cfg.Tj_M_est; else, Tj_M_est = T_amb; end
    Pout = Vo * Io; Iin_ideal = Io * Vo / Vin; Iin_adj = Iin_ideal / eta_est;
    relax_D = 0.5; D = (Vo + Vd) / (Vin + Vo + Vd);
    for iter = 1:30
        Iin_est = Io * Vo / (Vin * eta_est); Isw_avg = Iin_est + Io;
        D_new = (Vo + Vd + Io*DCR2) / (Vin + Vo + Vd - Iin_est*DCR1 - Isw_avg*Rds_on);
        D_next = relax_D * D_new + (1 - relax_D) * D;
        if abs(D_next - D) < 1e-7, D = D_next; break; end
        D = D_next;
    end
    L_req = (Vin * D) / (Iin_adj * ripple_factor * fsw); L2_req = (Vin * D) / (Io * ripple_factor * fsw);
    if L1_in > 0, L1_use = L1_in; else, L1_use = L_req; end; if L2_in > 0, L2_use = L2_in; else, L2_use = L2_req; end
    if useCoupled
        L_eff = L1_use * (1 - k_c); delta_IL1 = (Vin * D) / (L_eff * fsw); delta_IL2 = delta_IL1;
    else
        delta_IL1 = (Vin * D) / (L1_use * fsw); delta_IL2 = (Vin * D) / (L2_use * fsw);
    end
    IL1_peak = Iin_adj + delta_IL1/2; IL2_peak = Io + delta_IL2/2; IL1_min = Iin_adj - delta_IL1/2; IL2_min = Io - delta_IL2/2;
    IL1_rms = sqrt(Iin_adj^2 + delta_IL1^2/12); IL2_rms = sqrt(Io^2 + delta_IL2^2/12);
    Isw_base = Iin_adj + Io; Isw_peak = Isw_base + (delta_IL1 + delta_IL2)/2; Isw_valley = max(Isw_base - (delta_IL1 + delta_IL2)/2, 0);
    Isw_rms = sqrt(D) * sqrt(Isw_base^2 + (delta_IL1 + delta_IL2)^2/12);
    Vsw_theory = Vin + Vo + Vd; Isw_step = Iin_adj + Io;
    ICc_rms = sqrt(D * (Io^2 + delta_IL2^2/12) + (1-D) * (Iin_adj^2 + delta_IL1^2/12));
    ICout_rms = delta_IL2 / sqrt(12); delta_Itotal = delta_IL1 + delta_IL2;
    ICin_rms = sqrt(D * (Isw_base^2 + delta_Itotal^2/12) + (1-D) * Iin_adj^2);
    Pmos_cond = Isw_rms^2 * Rds_on; Pmos_sw_overlap_transition = 0.5 * Vsw_theory * (Isw_valley + Isw_peak) * (t_sw_rise + t_sw_fall) * fsw;
    Ig_on = max((Vgate - Vgs3)/Rg, 0.01); Ig_off = max(Vgs3/Rg, 0.01);
    t_miller_on = Qgd / Ig_on; t_miller_off = Qgd / Ig_off; Pmos_sw_overlap_miller = 0.5 * Vsw_theory * Isw_peak * (t_miller_on + t_miller_off) * fsw;
    dpt = dpt_params();
    if numel(dpt.V_ref) > 1
        E_on_ref_v = interp1(dpt.V_ref, dpt.E_on_ref, Vsw_theory, 'linear','extrap'); E_off_ref_v = interp1(dpt.V_ref, dpt.E_off_ref, Vsw_theory, 'linear','extrap');
        I_on_ref_v = interp1(dpt.V_ref, dpt.I_on_ref, Vsw_theory, 'linear','extrap'); I_off_ref_v = interp1(dpt.V_ref, dpt.I_off_ref, Vsw_theory, 'linear','extrap');
    else
        E_on_ref_v = dpt.E_on_ref * (Vsw_theory / dpt.V_ref); E_off_ref_v = dpt.E_off_ref * (Vsw_theory / dpt.V_ref);
        I_on_ref_v = dpt.I_on_ref(1); I_off_ref_v = dpt.I_off_ref(1);
    end
    I_on_ref_v = max(I_on_ref_v, 0.01); I_off_ref_v = max(I_off_ref_v, 0.01);
    E_on = E_on_ref_v * (Isw_valley / I_on_ref_v); E_off = E_off_ref_v * (Isw_peak / I_off_ref_v);
    E_on = E_on * (1 + kE_on * (Tj_M_est - 25)); E_off = E_off * (1 + kE_off * (Tj_M_est - 25));
    Pmos_sw_overlap_dpt = (E_on + E_off) * fsw; Pmos_sw_gate = Vgate * Qg * fsw;
    Pmos_sw = Pmos_sw_overlap_dpt + Pmos_sw_gate + 0; Pmos_total = Pmos_cond + Pmos_sw;
    Pdiode_cond = Vd * Io; Pdiode_cap = 0.5 * cfg.Cj_diode * Vsw_theory^2 * fsw; Pdiode = Pdiode_cond + Pdiode_cap;
    PL1_copper = IL1_rms^2 * DCR1; PL2_copper = IL2_rms^2 * DCR2; Pind_total = PL1_copper + PL2_copper;
    Pcap_Cc = ICc_rms^2 * ESR_Cc; Pcap_Cout = ICout_rms^2 * ESR_Cout; Pcap_Cin = ICin_rms^2 * ESR_Cin; Pcap_total = Pcap_Cc + Pcap_Cout + Pcap_Cin;
    Pmos_stage = Pmos_cond + Pmos_sw_overlap_dpt + 0; Ploss_stage = Pmos_stage + Pdiode + Pind_total + Pcap_total;
    Ploss_system = Ploss_stage + Pmos_sw_gate; Pin_system_est = Pout + Ploss_system;
    eff_system = 100 * Pout / Pin_system_est; res.eff_system = eff_system;
    res.Tj_Q1 = T_amb + Pmos_stage * RJA_Q1; res.Tj_D1 = T_amb + Pdiode * RJA_D1;
end


function d = dpt_params()
% Load the DPT lookup CSV once and cache it. Signature includes size and
% mtime, so editing the CSV forces a reload even in a long-running
% Monte Carlo where this gets called hundreds of thousands of times.
    persistent cached; persistent cached_signature;
    candidates = { fullfile(pwd, 'dpt_lookup.csv') };
    csv_path = '';
    for k = 1:numel(candidates)
        if isfile(candidates{k}), csv_path = candidates{k}; break; end
    end
    if isempty(csv_path), signature = '__NO_FILE__';
    else, info = dir(csv_path); signature = sprintf('%s|%d|%f', csv_path, info.bytes, info.datenum); end
    if isempty(cached) || isempty(cached_signature) || ~strcmp(signature, cached_signature)
        if ~isempty(csv_path)
            T = readtable(csv_path); T = sortrows(T, 'V_bus_V');
            d.V_ref = T.V_bus_V(:)'; d.E_on_ref = T.E_on_uJ(:)' * 1e-6; d.E_off_ref = T.E_off_uJ(:)' * 1e-6;
            if all(ismember({'I_turnoff_A','I_turnon_A'}, T.Properties.VariableNames))
                d.I_off_ref = T.I_turnoff_A(:)'; d.I_on_ref = T.I_turnon_A(:)';
            else
                d.I_off_ref = 5.386 * ones(size(d.V_ref)); d.I_on_ref = 5.386 * ones(size(d.V_ref));
            end
        else
            error('dpt_lookup.csv not found!');
        end
        cached = d; cached_signature = signature;
    else
        d = cached;
    end
end


%% ========================================================================
%  Figure-style helpers
%  Same engine used by the SEPIC and ZETA scripts: one persistent preset,
%  pushed into MATLAB's defaults on every style switch.
%% ========================================================================

function S = ieee_cfg(style)
    persistent cur
    if nargin > 0, cur = lower(strtrim(style)); end
    if isempty(cur), cur = 'ieee'; end

    S.name = cur;
    switch cur
        case 'ieee'
            S.col1 = 3.5;   S.col2 = 7.16;
            S.font = 8;     S.fontSmall = 8;    S.fontTiny = 7;
            S.lw = 1.3;     S.ms = 6;           S.axLW = 0.5;
            S.dpi = 600;    S.titles = false;   S.panelBold = false;
            S.exts = {'pdf','eps','png'};
            S.interp = 'latex';
            S.fontPref = {'Times New Roman','Times','Liberation Serif','Nimbus Roman','DejaVu Serif'};
        case 'report'
            S.col1 = 5.0;   S.col2 = 6.8;
            S.font = 10;    S.fontSmall = 9;    S.fontTiny = 8;
            S.lw = 1.6;     S.ms = 7;           S.axLW = 0.7;
            S.dpi = 300;    S.titles = true;    S.panelBold = true;
            S.exts = {'pdf','png'};
            S.interp = 'tex';
            S.fontPref = {'Calibri','Arial','Helvetica','Liberation Sans','DejaVu Sans'};
        otherwise
            error('ieee_cfg:style', 'Unknown figure style "%s" (use ''ieee'' or ''report'').', cur);
    end
    S.fontName = pick_font(S.fontPref);
    S.ink      = [0.15 0.15 0.15];
end


function name = pick_font(prefs)
    name = get(0, 'DefaultAxesFontName');
    try
        avail = listfonts;
        for k = 1:numel(prefs)
            if any(strcmpi(avail, prefs{k})), name = prefs{k}; return; end
        end
    catch
        % listfonts unavailable -- keep MATLAB's default
    end
end


function apply_ieee_defaults()
    S = ieee_cfg();
    set(0, 'DefaultAxesFontName',               S.fontName);
    set(0, 'DefaultTextFontName',               S.fontName);
    set(0, 'DefaultLegendFontName',             S.fontName);
    set(0, 'DefaultColorbarFontName',           S.fontName);
    set(0, 'DefaultAxesFontSize',               S.font);
    set(0, 'DefaultTextFontSize',               S.font);
    set(0, 'DefaultLegendFontSize',             S.fontSmall);
    set(0, 'DefaultColorbarFontSize',           S.fontSmall);
    set(0, 'DefaultAxesLabelFontSizeMultiplier',1);
    set(0, 'DefaultAxesTitleFontSizeMultiplier',1);
    set(0, 'DefaultAxesTitleFontWeight',        'normal');
    set(0, 'DefaultTextInterpreter',            S.interp);
    set(0, 'DefaultAxesTickLabelInterpreter',   S.interp);
    set(0, 'DefaultColorbarTickLabelInterpreter', S.interp);
    set(0, 'DefaultLegendInterpreter',          S.interp);
    set(0, 'DefaultLineLineWidth',              S.lw);
    set(0, 'DefaultLineMarkerSize',             S.ms);
    set(0, 'DefaultAxesLineWidth',              S.axLW);
    set(0, 'DefaultAxesBox',                    'on');
    set(0, 'DefaultAxesTickDir',                'in');
    set(0, 'DefaultAxesTickLength',             [0.012 0.025]);
    set(0, 'DefaultAxesXColor',                 S.ink);
    set(0, 'DefaultAxesYColor',                 S.ink);
    set(0, 'DefaultAxesGridColor',              [0.5 0.5 0.5]);
    set(0, 'DefaultAxesGridAlpha',              0.35);
    set(0, 'DefaultAxesGridLineStyle',          ':');
    set(0, 'DefaultAxesMinorGridAlpha',         0.12);
    set(0, 'DefaultFigureColor',                'w');
    set(0, 'DefaultLegendBox',                  'off');
end


function fig = new_fig(name, w, h)
% Figure at exact printed size. Height is given for the IEEE width and
% scaled for other presets.
    S = ieee_cfg();
    if abs(w - S.col1) < 1e-9, h = h * S.col1 / 3.5;
    else,                      h = h * S.col2 / 7.16; end
    fig = figure('Name', name, 'NumberTitle','off', 'Color','w', 'Units','inches', ...
        'Position', [0.5 0.5 w h], 'PaperUnits','inches', 'PaperSize', [w h], ...
        'PaperPosition', [0 0 w h], 'PaperPositionMode','manual');
end


function ttl(h, str)
    S = ieee_cfg();
    if S.titles, title(h, str, 'FontWeight','normal', 'FontSize', S.font); end
end


function save_fig(fig, outDir, name)
% Write PDF (+EPS for IEEE) and PNG. If the primary path is locked, write
% a timestamped fallback instead of throwing.
    S = ieee_cfg();
    if ~isfolder(outDir), mkdir(outDir); end
    drawnow;
    useEG = exist('exportgraphics') ~= 0; %#ok<EXIST>
    if ~useEG
        set(fig, 'Units','inches');
        p = get(fig, 'Position');
        set(fig, 'PaperUnits','inches', 'PaperSize',p(3:4), ...
                 'PaperPosition',[0 0 p(3:4)], 'PaperPositionMode','manual');
    end
    exts = S.exts;
    for e = 1:numel(exts)
        target = fullfile(outDir, [name '.' exts{e}]);
        try
            write_one(fig, target, exts{e}, useEG);
        catch err
            alt = fullfile(outDir, sprintf('%s_%s.%s', name, datestr(now,'HHMMSS'), exts{e}));
            warning('save_fig:locked', ...
                'Could not write %s (%s).\n  Saved as %s instead - close the file in your viewer.', ...
                target, strtrim(err.message), alt);
            try
                write_one(fig, alt, exts{e}, useEG);
            catch err2
                warning('save_fig:failed', 'Skipping %s: %s', alt, strtrim(err2.message));
            end
        end
    end
end


function write_one(fig, target, ext, useEG)
% Single export call. Uses exportgraphics when available, print() otherwise.
    S = ieee_cfg();
    if useEG
        if strcmp(ext, 'png')
            args = {'Resolution', S.dpi, 'BackgroundColor','white'};
        else
            args = {'ContentType','vector', 'Resolution', S.dpi, 'BackgroundColor','white'};
        end
        try
            exportgraphics(fig, target, args{:}, 'Padding','figure');
        catch
            exportgraphics(fig, target, args{:});
        end
    else
        switch ext
            case 'pdf', print(fig, target, '-dpdf',  '-painters');
            case 'eps', print(fig, target, '-depsc', '-painters');
            otherwise,  print(fig, target, '-dpng', sprintf('-r%d', S.dpi));
        end
    end
end


function panel_tag(ax, str)
    S = ieee_cfg();
    if S.panelBold, fw = 'bold'; else, fw = 'normal'; end
    text(ax, 0.02, 0.97, str, 'Units','normalized', 'FontSize', S.font, ...
        'FontWeight', fw, 'VerticalAlignment','top', 'HorizontalAlignment','left');
end


function out = lx(s)
% Convert TeX-style markup to LaTeX when interp is 'latex'; pass-through
% otherwise. Multi-line strings become cell arrays (LaTeX needs that).
    S = ieee_cfg();
    if iscell(s), out = cellfun(@lx, s, 'UniformOutput', false); return; end
    if ~strcmp(S.interp, 'latex'), out = s; return; end
    s = char(s);
    if contains(s, newline)
        out = cellfun(@lx, strsplit(s, newline), 'UniformOutput', false);
        return;
    end
    expr = '\\[A-Za-z]+(?:_\{[^}]*\}|_[A-Za-z0-9])*|[A-Za-z]+[0-9]*(?:_\{[^}]*\}|_[A-Za-z0-9])+';
    [m, parts] = regexp(s, expr, 'match', 'split');
    out = tex_escape(parts{1});
    for k = 1:numel(m)
        out = [out, '$', tex_sub(m{k}), '$', tex_escape(parts{k+1})]; %#ok<AGROW>
    end
end


function out = lxf(fmt, varargin)
    out = lx(sprintf(fmt, varargin{:}));
end


function t = tex_escape(t)
    t = strrep(t, '%', '\%');
    t = strrep(t, '&', '\&');
    t = strrep(t, '#', '\#');
    t = strrep(t, char(176), '$^{\circ}$');
end


function t = tex_sub(t)
% Upright roman subscripts: V_{in} -> V_{\mathrm{in}}, I_o -> I_{\mathrm{o}}.
    t = regexprep(t, '_\{([^}]*)\}', '_{\mathrm{$1}}');
    t = regexprep(t, '_([A-Za-z0-9])', '_{\mathrm{$1}}');
end