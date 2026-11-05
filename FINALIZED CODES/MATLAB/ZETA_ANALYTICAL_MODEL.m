%% ZETA converter analytical model
% Steady-state efficiency, component losses, and junction temperatures over
% a (Vin, Io) grid using three switching-loss methods (DPT, transition,
% Miller). Optional coupled inductor via useCoupled / k_c.
%
% Requires MATLAB R2020a+ (tiledlayout, exportgraphics).

clc; clear; close all; clear functions;

%% Style and output location -------------------------------------------------

outDir    = 'report_figures';
styleList = {'ieee', 'report'};

%% Design inputs -------------------------------------------------------------

Vin_vec = [9, 12, 15];
Io_vec  = [0.5, 1.0, 1.5];

Vo            = 12;
fsw           = 340e3;
ripple_factor = 0.30;
rippleMethod  = 'simplified';

% Coupled-inductor option (off by default)
useCoupled = false;
k_c        = -0.95;

%% MOSFET --------------------------------------------------------------------

Rds_on_25 = 25e-3;
alpha_Rds = 0.0055;

Coss = 400e-12;
Qgd  = 25e-9;
Qg   = 48e-9;

Vgate = 5;
Vgs3  = 3.0;
Rg    = 10;

Vds_max = 55;
Id_max  = 47;

RJA_Q1    = 62;
Tj_max_Q1 = 175;

%% Diode ---------------------------------------------------------------------

Vd_25            = 0.5;
Cj_diode         = 250e-12;
temp_coeff_diode = -0.002;

Vrrm_diode  = 40;
Ifavg_diode = 3.0;
RJA_D1      = 55;
Tj_max_D1   = 150;

T_amb = 40;

%% Switching times and temp coefficients -------------------------------------

t_sw_rise = 84e-9;
t_sw_fall = 15e-9;

kE_on  = 0.0015;
kE_off = 0.0010;

%% Loss method(s) ------------------------------------------------------------

switchingLossMethod = 'compare';

if strcmpi(switchingLossMethod, 'compare')
    METHODS_TO_RUN = {'dpt', 'transition', 'miller'};
    PRIMARY_METHOD = 'dpt';
    RUN_COMPARISON = true;
else
    METHODS_TO_RUN = {switchingLossMethod};
    PRIMARY_METHOD = switchingLossMethod;
    RUN_COMPARISON = false;
end

%% Passive parasitics --------------------------------------------------------

DCR_L1 = 50e-3;   DCR_L2 = 50e-3;
ESR_Cc = 10e-3;   ESR_Cout = 10e-3;   ESR_Cin = 10e-3;

%% Ripple targets ------------------------------------------------------------

dVc_target  = 0.02 * Vo;
dVo_target  = 0.02 * Vo;
dVin_target = 0.01 * max(Vin_vec);

%% Inductor sizing sweep ----------------------------------------------------
% Same pattern as the SEPIC script: sweep, take the max requirement, round
% up to a standard value, use L1 = L2 for simplicity.

Cc_temp   = 100e-6;
Cout_temp = 100e-6;
Cin_temp  = 100e-6;

cfg_temp = build_cfg(Vo, fsw, Vd_25, temp_coeff_diode, ...
    ripple_factor, rippleMethod, useCoupled, k_c, ...
    Rds_on_25, alpha_Rds, Coss, Qgd, Qg, ...
    Vgate, Vgs3, Rg, t_sw_rise, t_sw_fall, PRIMARY_METHOD, ...
    Vds_max, Id_max, kE_on, kE_off, ...
    RJA_Q1, Tj_max_Q1, Cj_diode, Vrrm_diode, Ifavg_diode, RJA_D1, Tj_max_D1, T_amb, ...
    DCR_L1, DCR_L2, Cc_temp, Cout_temp, Cin_temp, ESR_Cc, ESR_Cout, ESR_Cin, ...
    dVc_target, dVo_target, dVin_target);

relax_eta = 0.7;
relax_T   = 0.5;

fprintf('=== ZETA Inductor Sizing Sweep ===\n\n');
fprintf('%-6s %-6s %-10s %-10s %-8s\n', 'Vin', 'Io', 'L1_req(uH)', 'L2_req(uH)', 'Conv');
fprintf('%s\n', repmat('-', 1, 44));

InductorSizing = struct([]);

for i = 1:numel(Vin_vec)
    for j = 1:numel(Io_vec)
        Vin_s = Vin_vec(i);
        Io_s  = Io_vec(j);

        eta_est  = 0.90;
        Tj_M_est = T_amb;
        Tj_D_est = T_amb;
        conv     = false;

        for iter = 1:20
            dT = Tj_M_est - 25;
            cfg_temp.Rds_on   = cfg_temp.Rds_on_25 * (1 + cfg_temp.alpha_Rds * dT);
            cfg_temp.Vd       = cfg_temp.Vd_25 + cfg_temp.temp_coeff_diode * (Tj_D_est - 25);
            cfg_temp.Tj_M_est = Tj_M_est;
            cfg_temp.Tj_D_est = Tj_D_est;

            res_s = calcZETA(Vin_s, Io_s, eta_est, cfg_temp, 0, 0);

            eta_new  = res_s.eff_system / 100;
            Tj_M_new = res_s.Tj_Q1;
            Tj_D_new = res_s.Tj_D1;

            d_eta = abs(eta_new  - eta_est);
            d_Tq  = abs(Tj_M_new - Tj_M_est);
            d_Td  = abs(Tj_D_new - Tj_D_est);

            eta_est  = relax_eta * eta_new  + (1 - relax_eta) * eta_est;
            Tj_M_est = relax_T   * Tj_M_new + (1 - relax_T)   * Tj_M_est;
            Tj_D_est = relax_T   * Tj_D_new + (1 - relax_T)   * Tj_D_est;

            if d_eta < 1e-4 && d_Tq < 0.1 && d_Td < 0.1
                conv = true; break;
            end
        end

        dT = Tj_M_est - 25;
        cfg_temp.Rds_on   = cfg_temp.Rds_on_25 * (1 + cfg_temp.alpha_Rds * dT);
        cfg_temp.Vd       = cfg_temp.Vd_25 + cfg_temp.temp_coeff_diode * (Tj_D_est - 25);
        cfg_temp.Tj_M_est = Tj_M_est;
        cfg_temp.Tj_D_est = Tj_D_est;

        res_s = calcZETA(Vin_s, Io_s, eta_est, cfg_temp, 0, 0);
        res_s.converged    = conv;
        res_s.Rds_on_final = cfg_temp.Rds_on;

        if isempty(InductorSizing), InductorSizing = res_s;
        else,                       InductorSizing = [InductorSizing; res_s]; end

        if conv, conv_str = 'OK'; else, conv_str = 'WARN'; end
        fprintf('%-6.0f %-6.2f %-10.3f %-10.3f %-8s\n', Vin_s, Io_s, ...
            res_s.L_req * 1e6, res_s.L2_req * 1e6, conv_str);
    end
end

L1_req_all   = [InductorSizing.L_req];
L2_req_all   = [InductorSizing.L2_req];
L_fixed_calc = max([max(L1_req_all), max(L2_req_all)]);

std_inductors = [100, 120, 150, 180, 220, 270, 330] * 1e-6;
idxL = find(std_inductors >= L_fixed_calc, 1);
if isempty(idxL), error('No standard inductor large enough.'); end
L1_fixed = std_inductors(idxL);
L2_fixed = L1_fixed;

fprintf('\n=== Selected Fixed Inductors ===\n');
fprintf('  Max required L1 = %.2f uH, L2 = %.2f uH\n', ...
    max(L1_req_all) * 1e6, max(L2_req_all) * 1e6);
fprintf('  Practical chosen (matched pair) L1 = L2 = %.0f uH\n\n', L1_fixed * 1e6);

%% Capacitor sizing (ESR-aware) ---------------------------------------------
% Same idea as SEPIC but with the ZETA coupling capacitor Cc (step current
% Isw_step) and the buck-style output cap Cout on the L2 ripple.

Cc_req_list   = [];
Cout_req_list = [];
Cin_req_list  = [];

for idx = 1:numel(InductorSizing)
    r   = InductorSizing(idx);
    Vin = r.Vin;
    Io  = r.Io;

    cfg_temp2 = build_cfg(Vo, fsw, Vd_25, temp_coeff_diode, ...
        ripple_factor, rippleMethod, useCoupled, k_c, ...
        Rds_on_25, alpha_Rds, Coss, Qgd, Qg, ...
        Vgate, Vgs3, Rg, t_sw_rise, t_sw_fall, PRIMARY_METHOD, ...
        Vds_max, Id_max, kE_on, kE_off, ...
        RJA_Q1, Tj_max_Q1, Cj_diode, Vrrm_diode, Ifavg_diode, RJA_D1, Tj_max_D1, T_amb, ...
        DCR_L1, DCR_L2, 100e-6, 100e-6, 100e-6, ESR_Cc, ESR_Cout, ESR_Cin, ...
        dVc_target, dVo_target, dVin_target);

    dT = r.Tj_Q1 - 25;
    cfg_temp2.Rds_on   = cfg_temp2.Rds_on_25 * (1 + cfg_temp2.alpha_Rds * dT);
    cfg_temp2.Vd       = cfg_temp2.Vd_25 + cfg_temp2.temp_coeff_diode * (r.Tj_D1 - 25);
    cfg_temp2.Tj_M_est = r.Tj_Q1;
    cfg_temp2.Tj_D_est = r.Tj_D1;

    res_fixed = calcZETA(Vin, Io, 0.9, cfg_temp2, L1_fixed, L2_fixed);

    D         = res_fixed.D;
    delta_IL1 = res_fixed.delta_IL1;
    delta_IL2 = res_fixed.delta_IL2;
    Iin_avg   = res_fixed.Iin_adj;
    Isw_step  = Iin_avg + Io;

    dVc_esr = Isw_step * ESR_Cc;
    if dVc_target <= dVc_esr, error('Cc ESR alone exceeds target at Vin=%.1f Io=%.1f', Vin, Io); end
    Cc_req = (D * Io) / ((dVc_target - dVc_esr) * fsw);

    dVo_esr = delta_IL2 * ESR_Cout;
    if dVo_target <= dVo_esr, error('Cout ESR alone exceeds target at Vin=%.1f Io=%.1f', Vin, Io); end
    Cout_req = delta_IL2 / (8 * fsw * (dVo_target - dVo_esr));

    dVin_esr = Isw_step * ESR_Cin;
    if dVin_target <= dVin_esr, error('Cin ESR alone exceeds target at Vin=%.1f Io=%.1f', Vin, Io); end
    Cin_req = (D * Io) / ((dVin_target - dVin_esr) * fsw);

    Cc_req_list   = [Cc_req_list;   Cc_req];
    Cout_req_list = [Cout_req_list; Cout_req];
    Cin_req_list  = [Cin_req_list;  Cin_req];
end

std_caps = [1.0,1.2,1.5,1.8,2.2,2.7,3.3,3.9,4.7,5.6,6.8,8.2,10,12,15, ...
            18,22,27,33,39,47,56,68,82,100] * 1e-6;

idxCc   = find(std_caps >= max(Cc_req_list),   1);
idxCout = find(std_caps >= max(Cout_req_list), 1);
idxCin  = find(std_caps >= max(Cin_req_list),  1);

if isempty(idxCc) || isempty(idxCout) || isempty(idxCin)
    error('Capacitor sizing failed: no standard value large enough.');
end

Cc_val   = std_caps(idxCc);
Cout_val = std_caps(idxCout);
Cin_val  = std_caps(idxCin);

fprintf('=== Selected Fixed Capacitors (ESR-aware) ===\n');
fprintf('  Cc=%.1f uF, Cout=%.1f uF, Cin=%.1f uF\n\n', Cc_val * 1e6, Cout_val * 1e6, Cin_val * 1e6);

%% Final configuration -------------------------------------------------------

cfg = build_cfg(Vo, fsw, Vd_25, temp_coeff_diode, ...
    ripple_factor, rippleMethod, useCoupled, k_c, ...
    Rds_on_25, alpha_Rds, Coss, Qgd, Qg, ...
    Vgate, Vgs3, Rg, t_sw_rise, t_sw_fall, PRIMARY_METHOD, ...
    Vds_max, Id_max, kE_on, kE_off, ...
    RJA_Q1, Tj_max_Q1, Cj_diode, Vrrm_diode, Ifavg_diode, RJA_D1, Tj_max_D1, T_amb, ...
    DCR_L1, DCR_L2, Cc_val, Cout_val, Cin_val, ESR_Cc, ESR_Cout, ESR_Cin, ...
    dVc_target, dVo_target, dVin_target);

%% Main evaluation loop -----------------------------------------------------

ResultsByMethod = struct();

for m_idx = 1:numel(METHODS_TO_RUN)

    current_method = METHODS_TO_RUN{m_idx};
    cfg.switchingLossMethod = current_method;

    fprintf('\n=== Evaluating with method: %s ===\n', current_method);

    Results_m = struct([]);

    for i = 1:numel(Vin_vec)
        for j = 1:numel(Io_vec)
            Vin = Vin_vec(i);
            Io  = Io_vec(j);

            eta_est   = 0.90;
            Tj_M_est  = T_amb;
            Tj_D_est  = T_amb;
            converged = false;

            for iter = 1:20
                dT = Tj_M_est - 25;
                cfg.Rds_on   = cfg.Rds_on_25 * (1 + cfg.alpha_Rds * dT);
                cfg.Vd       = cfg.Vd_25 + cfg.temp_coeff_diode * (Tj_D_est - 25);
                cfg.Tj_M_est = Tj_M_est;
                cfg.Tj_D_est = Tj_D_est;

                res = calcZETA(Vin, Io, eta_est, cfg, L1_fixed, L2_fixed);

                eta_new  = res.eff_system / 100;
                Tj_M_new = res.Tj_Q1;
                Tj_D_new = res.Tj_D1;

                d_eta = abs(eta_new  - eta_est);
                d_Tq  = abs(Tj_M_new - Tj_M_est);
                d_Td  = abs(Tj_D_new - Tj_D_est);

                eta_est  = relax_eta * eta_new  + (1 - relax_eta) * eta_est;
                Tj_M_est = relax_T   * Tj_M_new + (1 - relax_T)   * Tj_M_est;
                Tj_D_est = relax_T   * Tj_D_new + (1 - relax_T)   * Tj_D_est;

                if d_eta < 1e-4 && d_Tq < 0.1 && d_Td < 0.1
                    converged = true; break;
                end
            end

            dT = Tj_M_est - 25;
            cfg.Rds_on   = cfg.Rds_on_25 * (1 + cfg.alpha_Rds * dT);
            cfg.Vd       = cfg.Vd_25 + cfg.temp_coeff_diode * (Tj_D_est - 25);
            cfg.Tj_M_est = Tj_M_est;
            cfg.Tj_D_est = Tj_D_est;

            res = calcZETA(Vin, Io, eta_est, cfg, L1_fixed, L2_fixed);
            res.converged    = converged;
            res.Rds_on_final = cfg.Rds_on;
            res.switching_loss_method = current_method;

            res.Isw_peak_status = ternary(res.Isw_peak <= Id_max, 'PASS', 'FAIL Id');
            rev_margin          = 1.30 * (Vin + Vo);
            res.Vdiode_status   = ternary(rev_margin <= Vrrm_diode, 'PASS', 'FAIL Vrrm');
            res.Idiode_status   = ternary(Io <= Ifavg_diode, 'PASS', 'FAIL If');

            if isempty(Results_m), Results_m = res;
            else,                  Results_m = [Results_m; res]; end
        end
    end

    ResultsByMethod.(current_method) = Results_m;
end

Results = ResultsByMethod.(PRIMARY_METHOD);
N       = numel(Results);

%% Summary table ------------------------------------------------------------

T = struct2table(Results);

keepCols = {'Vin','Io','Pout','D', ...
    'eff_stage','eff_system', ...
    'eff_stage_miller','eff_system_miller', ...
    'Ploss_stage','Ploss_system', ...
    'Ploss_stage_miller','Ploss_system_miller', ...
    'Pmos_sw_overlap','Pmos_sw_overlap_miller','Pmos_sw_overlap_transition','Pmos_sw_overlap_dpt', ...
    'Pmos_sw_gate','dVo_actual','Tj_Q1','Tj_D1', ...
    'Vsw_status','Isw_peak_status','Vdiode_status'};

SummaryTable = T(:, keepCols);
SummaryTable.Properties.VariableNames = {'Vin_V','Io_A','Pout_W','Duty', ...
    'StageEfficiency_pct','SystemEfficiency_pct', ...
    'StageEfficiency_Miller_pct','SystemEfficiency_Miller_pct', ...
    'StageLoss_W','SystemLoss_W', ...
    'StageLoss_Miller_W','SystemLoss_Miller_W', ...
    'SelectedMOSFETSwitchLoss_W','MillerSwitchLoss_W','TransitionSwitchLoss_W','DPTSwitchLoss_W', ...
    'GateDriveLoss_W','Ripple_V','Tj_MOSFET_C','Tj_Diode_C', ...
    'MOSFET_V','MOSFET_I','Diode_VDC'};

dpt_info = dpt_params();
tax      = model_taxonomy(PRIMARY_METHOD);

fprintf('\n========================================================\n');
fprintf('ZETA v10.7 FINAL SUMMARY TABLE\n');
fprintf('========================================================\n');
fprintf('  MOSFET: IRLZ44NPbF, RJA = 62 K/W, tr = %.0f ns, tf = %.0f ns\n', ...
    t_sw_rise * 1e9, t_sw_fall * 1e9);
fprintf('  L1 = L2 = %.0f uH, Cc = %.1f uF, Cout = %.1f uF, Cin = %.1f uF\n', ...
    L1_fixed * 1e6, Cc_val * 1e6, Cout_val * 1e6, Cin_val * 1e6);
fprintf('--------------------------------------------------------\n');
fprintf('  Primary method        : %s\n', PRIMARY_METHOD);
fprintf('  Category              : %s\n', tax.category);
fprintf('  Review mapping        : %s\n', tax.review);
if isnumeric(tax.generation)
    fprintf('  Generation            : %d\n', tax.generation);
else
    fprintf('  Generation            : %s\n', tax.generation);
end
fprintf('  Included              : %s\n', strjoin(tax.included, '; '));
fprintf('  Omitted               : %s\n', strjoin(tax.omitted,  '; '));
fprintf('--------------------------------------------------------\n');
fprintf('  DPT source            : %s\n', dpt_info.source);
fprintf('  DPT CSV path          : %s\n', dpt_info.csv_path);
if numel(dpt_info.V_ref) > 1
    fprintf('  DPT V sweep           : %.1f - %.1f V (%d points, interpolated)\n', ...
        min(dpt_info.V_ref), max(dpt_info.V_ref), numel(dpt_info.V_ref));
else
    fprintf('  DPT V reference       : %.0f V (single-point fallback)\n', dpt_info.V_ref);
end
if numel(dpt_info.I_on_ref) > 1
    fprintf('  DPT measurement note  : per-V current scaling (interpolated refs)\n');
else
    fprintf('  DPT measurement note  : scalar I_ref = %.3f A (single reference)\n', ...
        dpt_info.I_ref);
end
fprintf('  E_on/E_off T scaling  : kE_on = %.4f, kE_off = %.4f (per deg C)\n', ...
    kE_on, kE_off);
fprintf('  Thermal model         : steady-state limit of review Eq. (5)-(6)\n');
fprintf('                          (single R_JA per device, no transient Foster)\n');
fprintf('========================================================\n');
disp(SummaryTable);

if RUN_COMPARISON
    fprintf('[compare mode] Component breakdown below uses %s as primary.\n', PRIMARY_METHOD);
    fprintf('[compare mode] Cross-method comparison printed after validation.\n');
end

ComparisonReady = SummaryTable(:, { ...
    'Vin_V','Io_A','Pout_W','Duty', ...
    'StageEfficiency_pct','SystemEfficiency_pct', ...
    'StageEfficiency_Miller_pct','SystemEfficiency_Miller_pct', ...
    'StageLoss_W','SystemLoss_W', ...
    'StageLoss_Miller_W','SystemLoss_Miller_W', ...
    'SelectedMOSFETSwitchLoss_W','MillerSwitchLoss_W','TransitionSwitchLoss_W','DPTSwitchLoss_W', ...
    'GateDriveLoss_W','Ripple_V'});

%% Component-level loss breakdown --------------------------------------------

fprintf('\n========================================================\n');
fprintf('ZETA COMPONENT-LEVEL LOSS BREAKDOWN (Watts, stage vs. system)\n');
fprintf('========================================================\n');

hdr = sprintf('%-18s | %-8s %-8s | %-9s %-9s %-9s %-9s | %-10s %-10s %-10s | %-9s | %-9s %-9s | %-8s %-8s', ...
    'Vin / Io', 'Q1_cond','Q1_sw', 'D1_cond','D1_cap', 'L1_cu','L2_cu', ...
    'Cc_esr','Cout_esr','Cin_esr', 'Gate_drv', 'Stage_W','System_W', 'EtaStg','EtaSys');
fprintf('%s\n', hdr);
fprintf('%s\n', repmat('-', 1, length(hdr)));

compLossStruct = struct('Label',{}, 'Vin',{}, 'Io',{}, ...
    'Q1_cond_W',{}, 'Q1_sw_selected_W',{}, 'Q1_sw_miller_W',{}, ...
    'Q1_sw_transition_W',{}, 'Q1_sw_dpt_W',{}, 'Q1_sw_gate_W',{}, 'Q1_Coss_W',{}, ...
    'D1_cond_W',{}, 'D1_cap_W',{}, 'L1_copper_W',{}, 'L2_copper_W',{}, ...
    'Cc_ESR_W',{}, 'Cout_ESR_W',{}, 'Cin_ESR_W',{}, ...
    'StageLoss_W',{}, 'SystemLoss_W',{}, ...
    'StageEfficiency_pct',{}, 'SystemEfficiency_pct',{});

for k = 1:N
    r       = Results(k);
    lbl     = sprintf('Vin=%dV Io=%.1fA', r.Vin, r.Io);
    Pstage  = r.Ploss_stage;
    Psystem = r.Ploss_system;

    fprintf('%-18s | %8.6fW %8.6fW | %8.6fW %8.6fW | %8.6fW %8.6fW | %8.6fW %8.6fW %8.6fW | %8.6fW | %8.6fW %8.6fW | %5.2f%% %5.2f%%\n', ...
        lbl, r.Pmos_cond, r.Pmos_sw_overlap + r.Pmos_Coss, ...
        r.Pdiode_cond, r.Pdiode_cap, ...
        r.PL1_copper, r.PL2_copper, ...
        r.Pcap_Cc, r.Pcap_Cout, r.Pcap_Cin, ...
        r.Pmos_sw_gate, Pstage, Psystem, r.eff_stage, r.eff_system);

    compLossStruct(k).Label              = lbl;
    compLossStruct(k).Vin                = r.Vin;
    compLossStruct(k).Io                 = r.Io;
    compLossStruct(k).Q1_cond_W          = r.Pmos_cond;
    compLossStruct(k).Q1_sw_selected_W   = r.Pmos_sw_overlap + r.Pmos_Coss;
    compLossStruct(k).Q1_sw_miller_W     = r.Pmos_sw_overlap_miller;
    compLossStruct(k).Q1_sw_transition_W = r.Pmos_sw_overlap_transition;
    compLossStruct(k).Q1_sw_dpt_W        = r.Pmos_sw_overlap_dpt;
    compLossStruct(k).Q1_sw_gate_W       = r.Pmos_sw_gate;
    compLossStruct(k).Q1_Coss_W          = r.Pmos_Coss;
    compLossStruct(k).D1_cond_W          = r.Pdiode_cond;
    compLossStruct(k).D1_cap_W           = r.Pdiode_cap;
    compLossStruct(k).L1_copper_W        = r.PL1_copper;
    compLossStruct(k).L2_copper_W        = r.PL2_copper;
    compLossStruct(k).Cc_ESR_W           = r.Pcap_Cc;
    compLossStruct(k).Cout_ESR_W         = r.Pcap_Cout;
    compLossStruct(k).Cin_ESR_W          = r.Pcap_Cin;
    compLossStruct(k).StageLoss_W        = Pstage;
    compLossStruct(k).SystemLoss_W       = Psystem;
    compLossStruct(k).StageEfficiency_pct  = r.eff_stage;
    compLossStruct(k).SystemEfficiency_pct = r.eff_system;
end

ComponentLoss = struct2table(compLossStruct);
ComponentLoss.Q1_cond_pct  = 100 * ComponentLoss.Q1_cond_W          ./ ComponentLoss.SystemLoss_W;
ComponentLoss.Q1_sw_pct    = 100 * ComponentLoss.Q1_sw_selected_W   ./ ComponentLoss.SystemLoss_W;
ComponentLoss.Q1_gate_pct  = 100 * ComponentLoss.Q1_sw_gate_W       ./ ComponentLoss.SystemLoss_W;
ComponentLoss.D1_cond_pct  = 100 * ComponentLoss.D1_cond_W          ./ ComponentLoss.SystemLoss_W;
ComponentLoss.D1_cap_pct   = 100 * ComponentLoss.D1_cap_W           ./ ComponentLoss.SystemLoss_W;
ComponentLoss.L1_cu_pct    = 100 * ComponentLoss.L1_copper_W        ./ ComponentLoss.SystemLoss_W;
ComponentLoss.L2_cu_pct    = 100 * ComponentLoss.L2_copper_W        ./ ComponentLoss.SystemLoss_W;
ComponentLoss.Cc_ESR_pct   = 100 * ComponentLoss.Cc_ESR_W           ./ ComponentLoss.SystemLoss_W;
ComponentLoss.Cout_ESR_pct = 100 * ComponentLoss.Cout_ESR_W         ./ ComponentLoss.SystemLoss_W;
ComponentLoss.Cin_ESR_pct  = 100 * ComponentLoss.Cin_ESR_W          ./ ComponentLoss.SystemLoss_W;

fprintf('========================================================\n');

%% Validation ---------------------------------------------------------------

all_converged = all([Results.converged]);
all_ccm       = all(strcmp({Results.CCM_status}, 'CCM OK'));
all_ripple    = all(strcmp({Results.Cc_status},   'PASS')) && ...
                all(strcmp({Results.Cout_status}, 'PASS')) && ...
                all(strcmp({Results.Cin_status},  'PASS'));
all_voltage   = all(strcmp({Results.Vsw_status},    'PASS')) && ...
                all(strcmp({Results.Vdiode_status}, 'PASS'));
all_current   = all(strcmp({Results.Isw_peak_status}, 'PASS')) && ...
                all(strcmp({Results.Idiode_status},   'PASS'));
all_thermal   = all(strcmp({Results.MOSFET_thermal_status}, 'OK')) && ...
                all(strcmp({Results.Diode_thermal_status},  'OK'));

fprintf('\n========================================================\n');
fprintf('FINAL VALIDATION SUMMARY\n');
fprintf('========================================================\n');
fprintf('  Electro-thermal convergence : %s\n', ternary(all_converged, 'PASS','FAIL'));
fprintf('  CCM operation               : %s\n', ternary(all_ccm,       'PASS','FAIL'));
fprintf('  Ripple requirements         : %s\n', ternary(all_ripple,    'PASS','FAIL'));
fprintf('  Voltage-rating screening    : %s\n', ternary(all_voltage,   'PASS','FAIL'));
fprintf('  Current-rating screening    : %s\n', ternary(all_current,   'PASS','FAIL'));
fprintf('  Thermal screening           : %s\n', ternary(all_thermal,   'OK','WARNING'));
fprintf('========================================================\n');

%% Cross-method comparison ---------------------------------------------------

if RUN_COMPARISON

    R_dpt   = ResultsByMethod.dpt;
    R_trans = ResultsByMethod.transition;
    R_mill  = ResultsByMethod.miller;

    nCases = numel(R_dpt);

    fprintf('\n========================================================\n');
    fprintf('ZETA CROSS-METHOD COMPARISON — Switching Loss (Coss-normalised)\n');
    fprintf('  dpt = E_on+E_off (includes Coss recovery)\n');
    fprintf('  trans/mill = overlap + Coss (normalised to match dpt)\n');
    fprintf('========================================================\n');
    fprintf('%-12s | %12s %12s %12s | %9s %9s %9s\n', ...
        'Case', 'Psw_dpt', 'Psw_trans+Coss', 'Psw_mill+Coss', ...
        'Eff_dpt', 'Eff_trans', 'Eff_mill');
    fprintf('%s\n', repmat('-', 1, 100));

    for k = 1:nCases
        fprintf('%-12s | %12.5f %12.5f %12.5f | %8.3f%% %8.3f%% %8.3f%%\n', ...
            sprintf('%dV/%.1fA', R_dpt(k).Vin, R_dpt(k).Io), ...
            R_dpt(k).Pmos_sw_comparable_dpt, ...
            R_trans(k).Pmos_sw_comparable_transition, ...
            R_mill(k).Pmos_sw_comparable_miller, ...
            R_dpt(k).eff_system, ...
            R_trans(k).eff_system, ...
            R_mill(k).eff_system);
    end

    fprintf('\n========================================================\n');
    fprintf('METHOD SPREAD (max - min across the three comparable fields)\n');
    fprintf('========================================================\n');
    fprintf('%-12s | %12s | %12s | %10s\n', ...
        'Case', 'Delta_Psw_W', 'Delta_Eff_pp', 'Delta_Tj_C');
    fprintf('%s\n', repmat('-', 1, 56));

    for k = 1:nCases
        Psw = [R_dpt(k).Pmos_sw_comparable_dpt, ...
               R_trans(k).Pmos_sw_comparable_transition, ...
               R_mill(k).Pmos_sw_comparable_miller];
        Eff = [R_dpt(k).eff_system, ...
               R_trans(k).eff_system, ...
               R_mill(k).eff_system];
        Tj  = [R_dpt(k).Tj_Q1, R_trans(k).Tj_Q1, R_mill(k).Tj_Q1];

        fprintf('%-12s | %12.5f | %12.3f | %10.2f\n', ...
            sprintf('%dV/%.1fA', R_dpt(k).Vin, R_dpt(k).Io), ...
            max(Psw) - min(Psw), max(Eff) - min(Eff), max(Tj) - min(Tj));
    end
    fprintf('========================================================\n');

    cmpRows = {};
    for m_idx = 1:numel(METHODS_TO_RUN)
        mm = METHODS_TO_RUN{m_idx};
        Rm = ResultsByMethod.(mm);
        for k = 1:numel(Rm)
            cmpRows{end+1} = struct( ...
                'Vin_V',          Rm(k).Vin, ...
                'Io_A',           Rm(k).Io, ...
                'Method',         mm, ...
                'P_sw_comparable_W', Rm(k).Pmos_sw_comparable_dpt * strcmp(mm,'dpt') + ...
                                     Rm(k).Pmos_sw_comparable_transition * strcmp(mm,'transition') + ...
                                     Rm(k).Pmos_sw_comparable_miller * strcmp(mm,'miller'), ...
                'P_cond_W',       Rm(k).Pmos_cond, ...
                'P_mos_total_W',  Rm(k).Pmos_total, ...
                'P_stage_W',      Rm(k).Ploss_stage, ...
                'P_system_W',     Rm(k).Ploss_system, ...
                'Eff_stage_pct',  Rm(k).eff_stage, ...
                'Eff_system_pct', Rm(k).eff_system, ...
                'Tj_Q1_C',        Rm(k).Tj_Q1, ...
                'Tj_D1_C',        Rm(k).Tj_D1); %#ok<SAGROW>
        end
    end
    CmpTable = struct2table([cmpRows{:}]);

    writetable(CmpTable, 'ZETA_v10_7_method_comparison.csv');
    fprintf('\nComparison CSV exported: ZETA_v10_7_method_comparison.csv\n');
end

%% Figures ------------------------------------------------------------------
% Same structure as the SEPIC script: build the figure once per style
% preset. Data (lossMatrix, shortLabels, sepAt, ...) is prepared once.

clr_vin = [0.000 0.447 0.698; 0.835 0.369 0.000; 0.000 0.620 0.451];
ls_vin  = {'-', '--', '-.'};
mk_vin  = {'o', 's', '^'};
clr_io  = clr_vin;  ls_io = ls_vin;  mk_io = mk_vin;

clr_comp = [ ...
    0.00 0.25 0.45;  0.00 0.45 0.70;  0.45 0.70 0.90;  0.80 0.30 0.00; ...
    0.95 0.75 0.30;  0.00 0.45 0.35;  0.65 0.85 0.75;  0.50 0.30 0.55];

comp_names = {'Q1 cond.','Q1 sw.','Q1 gate','D1 cond.','D1 cap', ...
              'L1 Cu','L2 Cu','Cap ESR'};

lossMatrix = zeros(N, 8);
for k = 1:N
    r = Results(k);
    lossMatrix(k,1) = r.Pmos_cond;
    lossMatrix(k,2) = r.Pmos_sw_overlap + r.Pmos_Coss;
    lossMatrix(k,3) = r.Pmos_sw_gate;
    lossMatrix(k,4) = r.Pdiode_cond;
    lossMatrix(k,5) = r.Pdiode_cap;
    lossMatrix(k,6) = r.PL1_copper;
    lossMatrix(k,7) = r.PL2_copper;
    lossMatrix(k,8) = r.Pcap_Cc + r.Pcap_Cout + r.Pcap_Cin;
end

shortLabels = arrayfun(@(r) sprintf('%dV/%.1fA', r.Vin, r.Io), Results, ...
    'UniformOutput', false);

Vin_all  = [Results.Vin];
sepAt    = find(diff(Vin_all) ~= 0) + 0.5;
Vin_uniq = unique(Vin_all, 'stable');
grpCtr   = arrayfun(@(v) mean(find(Vin_all == v)), Vin_uniq);

for sIdx = 1:numel(styleList)
close all;
ieee_cfg(styleList{sIdx});
apply_ieee_defaults();
IEEE   = ieee_cfg();
figDir = fullfile(outDir, styleList{sIdx});
fprintf('\n--- Exporting "%s" figure set -> %s ---\n', styleList{sIdx}, figDir);

% Fig 1 — stage vs system efficiency, two panels, shared legend
fig1 = new_fig('Efficiency_StageSystem', IEEE.col2, 2.7);
tl1  = tiledlayout(fig1, 1, 2, 'TileSpacing','compact','Padding','compact');
effFields = {'StageEfficiency_pct','SystemEfficiency_pct'};
effTag    = {'(a) Stage','(b) System (incl. gate drive)'};
effLim    = pad_lim([SummaryTable.StageEfficiency_pct; SummaryTable.SystemEfficiency_pct], ...
                    0.15, false, 100);
for p = 1:2
    ax = nexttile(tl1); hold(ax,'on');
    for j = 1:numel(Io_vec)
        idx = abs(SummaryTable.Io_A - Io_vec(j)) < 1e-9;
        plot(ax, SummaryTable.Vin_V(idx), SummaryTable.(effFields{p})(idx), ...
            [ls_io{j} mk_io{j}], 'Color', clr_io(j,:), 'MarkerFaceColor', clr_io(j,:), ...
            'DisplayName', lxf('I_o = %.1f A', Io_vec(j)));
    end
    if p == 2
        [bestEff, bestIdx] = max(SummaryTable.SystemEfficiency_pct);
        plot(ax, SummaryTable.Vin_V(bestIdx), bestEff, 'p', 'MarkerSize', 9, ...
            'MarkerEdgeColor','k','MarkerFaceColor',[0.95 0.85 0.20],'HandleVisibility','off');
        text(ax, SummaryTable.Vin_V(bestIdx), bestEff, lxf('  max %.2f%%', bestEff), ...
            'VerticalAlignment','bottom','FontSize',IEEE.fontSmall);
    end
    xlabel(ax,lx('Input voltage V_{in} (V)'));
    if p == 1, ylabel(ax,lx('Efficiency \eta (%)')); else, set(ax,'YTickLabel',[]); end
    set(ax, 'XTick', Vin_vec);
    xlim(ax, pad_lim(Vin_vec, 0.08)); ylim(ax, effLim);
    grid(ax,'on'); panel_tag(ax, effTag{p});
    if p == 1, axFirst = ax; end
end
lgd = legend(axFirst, 'Orientation','horizontal', 'FontSize', IEEE.fontSmall);
lgd.Layout.Tile = 'north';
ttl(tl1, 'ZETA efficiency - stage vs. system');
save_fig(fig1, figDir, 'Fig1_ZETA_Efficiency_vs_Vin');

% Fig 2 — system efficiency vs load
fig2 = new_fig('Efficiency_vs_Io', IEEE.col1, 2.6);
ax2 = axes(fig2); hold(ax2,'on');
for i = 1:numel(Vin_vec)
    idx = abs(SummaryTable.Vin_V - Vin_vec(i)) < 1e-9;
    plot(ax2, SummaryTable.Io_A(idx), SummaryTable.SystemEfficiency_pct(idx), ...
        [ls_vin{i} mk_vin{i}], 'Color', clr_vin(i,:), 'MarkerFaceColor', clr_vin(i,:), ...
        'DisplayName', lxf('V_{in} = %d V', Vin_vec(i)));
end
xlabel(ax2,lx('Load current I_o (A)'));
ylabel(ax2,lx('System efficiency \eta_{sys} (%)'));
ttl(ax2,'System efficiency vs. load current');
legend(ax2, 'Location','best', 'FontSize', IEEE.fontSmall);
grid(ax2,'on');
set(ax2, 'XTick', Io_vec);
xlim(ax2, pad_lim(Io_vec, 0.08));
ylim(ax2, pad_lim(SummaryTable.SystemEfficiency_pct, 0.15, false, 100));
save_fig(fig2, figDir, 'Fig2_ZETA_Efficiency_vs_Io');

% Fig 3 — MOSFET and diode Tj
fig3 = new_fig('Thermal', IEEE.col2, 2.7);
tl3  = tiledlayout(fig3, 1, 2, 'TileSpacing','compact','Padding','compact');
thCol = {'Tj_MOSFET_C','Tj_Diode_C'};
thMax = [Tj_max_Q1, Tj_max_D1];
thTag = {'(a) MOSFET Q1','(b) Diode D1'};
for p = 1:2
    ax = nexttile(tl3); hold(ax,'on');
    for i = 1:numel(Vin_vec)
        idx = abs(SummaryTable.Vin_V - Vin_vec(i)) < 1e-9;
        plot(ax, SummaryTable.Io_A(idx), SummaryTable.(thCol{p})(idx), ...
            [ls_vin{i} mk_vin{i}], 'Color', clr_vin(i,:), 'MarkerFaceColor', clr_vin(i,:), ...
            'DisplayName', lxf('V_{in} = %d V', Vin_vec(i)));
    end
    Tall = SummaryTable.(thCol{p});
    if max(Tall) > 0.6 * thMax(p)
        yline(ax, thMax(p), '--', lxf('T_{j,max} = %d °C', thMax(p)), ...
            'Color',[0.75 0 0], 'LineWidth',1.0, 'Interpreter', IEEE.interp, 'LabelHorizontalAlignment','left', ...
            'LabelVerticalAlignment','bottom', 'FontSize',IEEE.fontSmall, 'HandleVisibility','off');
        ylim(ax, [min(Tall) - 5, thMax(p) * 1.05]);
    else
        ylim(ax, pad_lim([Tall; T_amb], 0.25, false));
        text(ax, 0.98, 0.04, lxf('T_{j,max} = %d °C (min. margin %.0f K)', ...
            thMax(p), thMax(p) - max(Tall)), 'Units','normalized', ...
            'HorizontalAlignment','right','VerticalAlignment','bottom', ...
            'FontSize',IEEE.fontSmall,'Color',[0.75 0 0]);
    end
    xlabel(ax,lx('Load current I_o (A)'));
    ylabel(ax,lx('Junction temperature T_j (°C)'));
    set(ax, 'XTick', Io_vec); xlim(ax, pad_lim(Io_vec, 0.08));
    grid(ax,'on'); panel_tag(ax, thTag{p});
    if p == 1, axFirst = ax; end
end
lgd = legend(axFirst, 'Orientation','horizontal', 'FontSize', IEEE.fontSmall);
lgd.Layout.Tile = 'north';
ttl(tl3, lxf('Junction temperature (T_{amb} = %d °C)', T_amb));
save_fig(fig3, figDir, 'Fig3_ZETA_Thermal');

% Fig 4 — component loss stack
fig4 = new_fig('ComponentLoss_Bar', IEEE.col2, 3.3);
ax4 = axes(fig4); hold(ax4,'on');
bh = bar(ax4, lossMatrix, 'stacked', 'BarWidth', 0.68);
for m = 1:8
    bh(m).FaceColor = clr_comp(m,:); bh(m).EdgeColor = 'w'; bh(m).LineWidth = 0.4;
end
totals = sum(lossMatrix, 2);
ymax4  = 1.30 * max(totals);
ylim(ax4, [0, ymax4]); xlim(ax4, [0.4, N + 0.6]);
for s = sepAt
    plot(ax4, [s s], [0 ymax4], ':', 'Color',[0.45 0.45 0.45], 'LineWidth',0.6, 'HandleVisibility','off');
end
for g = 1:numel(Vin_uniq)
    text(ax4, grpCtr(g), 0.985*ymax4, lxf('V_{in} = %d V', Vin_uniq(g)), ...
        'HorizontalAlignment','center','VerticalAlignment','top','FontSize',IEEE.font);
end
for k = 1:N
    text(ax4, k, totals(k) + 0.012*ymax4, ...
        lxf('%.3f W\n\\eta = %.1f%%', totals(k), Results(k).eff_system), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',IEEE.fontTiny);
end
set(ax4, 'XTick', 1:N, 'XTickLabel', shortLabels, 'XTickLabelRotation', 45, ...
    'TickLabelInterpreter','none');
xlabel(ax4,lx('Operating point (V_{in} / I_o)')); ylabel(ax4,'Power loss (W)');
ttl(ax4,'Component-level power loss breakdown');
legend(ax4, comp_names, 'Location','northoutside', 'NumColumns', 4, 'FontSize', IEEE.fontSmall);
grid(ax4,'on'); ax4.XGrid = 'off';
save_fig(fig4, figDir, 'Fig4_ZETA_ComponentLoss_Bar');

% Fig 5 — worst-case loss budget
[~, worst] = max(totals); r_worst = Results(worst); vals = lossMatrix(worst,:);
[~, order] = sort(vals, 'descend');
fig5 = new_fig('WorstCase_Budget', IEEE.col1, 2.5);
ax5 = axes(fig5);
bh5 = barh(ax5, vals(order), 'FaceColor','flat', 'EdgeColor','none', 'BarWidth', 0.72);
bh5.CData = clr_comp(order,:);
set(ax5, 'YTick', 1:8, 'YTickLabel', comp_names(order), 'YDir', 'reverse');
xlabel(ax5,'Power loss (W)');
ttl(ax5, lxf('Loss budget, worst case (V_{in} = %d V, I_o = %.1f A)', r_worst.Vin, r_worst.Io));
grid(ax5,'on'); ax5.YGrid = 'off';
xlim(ax5, [0, 1.65 * max(vals)]);
xlims = xlim(ax5);
for m = 1:8
    w = vals(order(m));
    text(ax5, w + 0.015*diff(xlims), m, lxf('%.3f W (%.1f%%)', w, 100*w/sum(vals)), ...
        'VerticalAlignment','middle','FontSize',IEEE.fontSmall);
end
text(ax5, 0.98, 0.05, lxf('V_{in} = %d V, I_o = %.1f A\n\\eta = %.1f%%, total %.3f W', ...
    r_worst.Vin, r_worst.Io, r_worst.eff_system, sum(vals)), 'Units','normalized', ...
    'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',IEEE.fontSmall);
save_fig(fig5, figDir, 'Fig5_ZETA_WorstCase_Budget');

% Fig 6 — fsw sensitivity
fsw_sweep = [100e3, 200e3, 340e3, 500e3, 700e3, 1000e3];
Vin_s = 12; Io_s = 1.0;
nPts = numel(fsw_sweep);
P_q1c_sw = zeros(1,nPts); P_q1sw_sw = zeros(1,nPts); P_d1_sw = zeros(1,nPts);
P_L_sw = zeros(1,nPts); P_cap_sw = zeros(1,nPts); eff_sw = zeros(1,nPts);
for fi = 1:nPts
    cfg_s = build_cfg(Vo, fsw_sweep(fi), Vd_25, temp_coeff_diode, ...
        ripple_factor, rippleMethod, useCoupled, k_c, ...
        Rds_on_25, alpha_Rds, Coss, Qgd, Qg, ...
        Vgate, Vgs3, Rg, t_sw_rise, t_sw_fall, PRIMARY_METHOD, ...
        Vds_max, Id_max, kE_on, kE_off, ...
        RJA_Q1, Tj_max_Q1, Cj_diode, Vrrm_diode, Ifavg_diode, RJA_D1, Tj_max_D1, T_amb, ...
        DCR_L1, DCR_L2, Cc_val, Cout_val, Cin_val, ESR_Cc, ESR_Cout, ESR_Cin, ...
        dVc_target, dVo_target, dVin_target);
    cfg_s.Rds_on = cfg_s.Rds_on_25; cfg_s.Vd = cfg_s.Vd_25;
    cfg_s.Tj_M_est = T_amb; cfg_s.Tj_D_est = T_amb;
    res_s = calcZETA(Vin_s, Io_s, 0.90, cfg_s, L1_fixed, L2_fixed);
    P_q1c_sw(fi)  = res_s.Pmos_cond;
    P_q1sw_sw(fi) = res_s.Pmos_sw;
    P_d1_sw(fi)   = res_s.Pdiode;
    P_L_sw(fi)    = res_s.Pind_total;
    P_cap_sw(fi)  = res_s.Pcap_total;
    eff_sw(fi)    = res_s.eff_system;
end
fsw_kHz = fsw_sweep/1e3; design_kHz = fsw/1e3;
all_losses = [P_q1c_sw, P_q1sw_sw, P_d1_sw, P_L_sw, P_cap_sw];

fig6 = new_fig('Fsw_Sensitivity', IEEE.col1, 4.2);
tl6  = tiledlayout(fig6, 2, 1, 'TileSpacing','compact','Padding','compact');

axA = nexttile(tl6); hold(axA,'on');
lossSeries = {P_q1c_sw, P_q1sw_sw, P_d1_sw, P_L_sw, P_cap_sw};
lossNames  = {'P_{Q1,cond}','P_{Q1,sw}','P_{D1}','P_{L}','P_{cap}'};
lossClr    = clr_comp([1 2 4 6 8], :);
lossMk     = {'o','s','^','d','v'};
lossLs     = {'-','--','-.',':','-'};
for q = 1:5
    plot(axA, fsw_kHz, lossSeries{q}, [lossLs{q} lossMk{q}], 'Color', lossClr(q,:), ...
        'MarkerFaceColor', lossClr(q,:), 'DisplayName', lx(lossNames{q}));
end
ylabel(axA,'Power loss (W)'); ylim(axA, [0, 1.75*max(all_losses)]);
legend(axA, 'Location','north', 'NumColumns', 2, 'FontSize', IEEE.fontSmall);
panel_tag(axA, '(a)');
ttl(axA, lxf('Loss and efficiency vs. f_{sw} (V_{in} = %d V, I_o = %.1f A)', Vin_s, Io_s));

axB = nexttile(tl6); hold(axB,'on');
plot(axB, fsw_kHz, eff_sw, '-o', 'Color',[0.1 0.1 0.1], 'MarkerFaceColor',[0.1 0.1 0.1]);
[effPk, pkI] = max(eff_sw);
plot(axB, fsw_kHz(pkI), effPk, 'p', 'MarkerSize',9, 'MarkerEdgeColor','k', ...
    'MarkerFaceColor',[0.95 0.85 0.20], 'HandleVisibility','off');
text(axB, fsw_kHz(pkI), effPk, lxf('  peak %.2f%% @ %g kHz', effPk, fsw_kHz(pkI)), ...
    'VerticalAlignment','bottom','FontSize',IEEE.fontSmall);
ylabel(axB,lx('System efficiency (%)')); ylim(axB, pad_lim(eff_sw, 0.35, false, 100));
xlabel(axB,lx('Switching frequency f_{sw} (kHz)'));
panel_tag(axB, '(b)');

for ax = [axA axB]
    set(ax, 'XScale','log', 'XTick', fsw_kHz, ...
        'XTickLabel', arrayfun(@(x) sprintf('%g', x), fsw_kHz, 'UniformOutput', false), ...
        'XMinorTick','off');
    xlim(ax, [min(fsw_kHz)*0.9, max(fsw_kHz)*1.1]);
    xline(ax, design_kHz, ':', 'Color',[0.35 0.35 0.35], 'LineWidth',0.8, 'HandleVisibility','off');
    grid(ax,'on');
end
text(axB, design_kHz, min(ylim(axB)), sprintf(' design: %g kHz', design_kHz), ...
    'Color',[0.35 0.35 0.35], 'VerticalAlignment','bottom', 'FontSize',IEEE.fontSmall);
set(axA, 'XTickLabel', []);
linkaxes([axA axB], 'x');
save_fig(fig6, figDir, 'Fig6_ZETA_Fsw_Sensitivity');

% Fig 7 — MOSFET loss split
mos_labels = {'Conduction','Switching overlap','C_{oss}','Gate drive'};
if strcmpi(PRIMARY_METHOD, 'dpt')
    mos_labels{3} = 'C_{oss} (included in overlap)';
end
mos_colors = [0.80 0.30 0.00; 0.95 0.75 0.30; 0.00 0.45 0.35; 0.00 0.25 0.45];

mosM = zeros(N, 4);
for k = 1:N
    r = Results(k);
    mosM(k,:) = [r.Pmos_cond, r.Pmos_sw_overlap, r.Pmos_Coss, r.Pmos_sw_gate];
end
mosShare = 100 * mosM ./ sum(mosM, 2);

fig7 = new_fig('MOSFET_SubLoss', IEEE.col2, 3.0);
tl7  = tiledlayout(fig7, 1, 2, 'TileSpacing','compact','Padding','compact');
for p = 1:2
    ax = nexttile(tl7); hold(ax,'on');
    if p == 1, Y = mosM; else, Y = mosShare; end
    b7 = bar(ax, Y, 'stacked', 'BarWidth', 0.68);
    for m = 1:4
        b7(m).FaceColor = mos_colors(m,:); b7(m).EdgeColor = 'w'; b7(m).LineWidth = 0.4;
    end
    if p == 1
        yTop = 1.12 * max(sum(Y,2));
        ylabel(ax,'MOSFET power loss (W)');
        for k = 1:N
            text(ax, k, sum(Y(k,:)) + 0.01*yTop, sprintf('%.3f', sum(Y(k,:))), ...
                'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',IEEE.fontTiny);
        end
        panel_tag(ax, '(a) Absolute');
    else
        yTop = 100;
        ylabel(ax,lx('Share of MOSFET loss (%)'));
        cum = cumsum(Y, 2);
        for k = 1:N
            for m = 1:4
                if Y(k,m) >= 6
                    c = cum(k,m) - Y(k,m)/2;
                    lum = [0.299 0.587 0.114] * mos_colors(m,:)';
                    if lum < 0.5, tc = 'w'; else, tc = 'k'; end
                    text(ax, k, c, lxf('%.0f%%', Y(k,m)), 'HorizontalAlignment','center', ...
                        'VerticalAlignment','middle','FontSize',IEEE.fontTiny,'Color',tc);
                end
            end
        end
        panel_tag(ax, '(b) Relative');
    end
    ylim(ax, [0, yTop]); xlim(ax, [0.4, N + 0.6]);
    for s = sepAt
        plot(ax, [s s], [0 yTop], ':', 'Color',[0.45 0.45 0.45], 'LineWidth',0.6, 'HandleVisibility','off');
    end
    set(ax, 'XTick', 1:N, 'XTickLabel', shortLabels, 'XTickLabelRotation', 45, ...
        'TickLabelInterpreter','none');
    xlabel(ax,lx('Operating point (V_{in} / I_o)'));
    grid(ax,'on'); ax.XGrid = 'off';
end
lgd = legend(ax, lx(mos_labels), 'Orientation','horizontal','Location','southoutside', ...
    'FontSize', IEEE.fontSmall);
lgd.Layout.Tile = 'south';
ttl(tl7, 'MOSFET loss decomposition');
save_fig(fig7, figDir, 'Fig7_ZETA_MOSFET_SubLoss');

% Fig 9 — heat-maps
effMat = nan(numel(Io_vec), numel(Vin_vec));
tjMat  = nan(numel(Io_vec), numel(Vin_vec));
for k = 1:N
    ii = find(abs(Vin_vec - SummaryTable.Vin_V(k)) < 1e-9, 1);
    jj = find(abs(Io_vec  - SummaryTable.Io_A(k))  < 1e-9, 1);
    effMat(jj, ii) = SummaryTable.SystemEfficiency_pct(k);
    tjMat(jj, ii)  = SummaryTable.Tj_MOSFET_C(k);
end
xt = arrayfun(@(v) sprintf('%d', v),   Vin_vec, 'UniformOutput', false);
yt = arrayfun(@(v) sprintf('%.1f', v), Io_vec,  'UniformOutput', false);

fig9 = new_fig('Heatmaps', IEEE.col2, 2.5);
tl9  = tiledlayout(fig9, 1, 2, 'TileSpacing','compact','Padding','compact');
ax = nexttile(tl9);
plot_heatmap(ax, effMat, xt, yt, '%.2f', '\eta_{sys} (%)', seq_cmap('viridis'));
panel_tag(ax, '(a) System efficiency');
ax = nexttile(tl9);
plot_heatmap(ax, tjMat, xt, yt, '%.1f', 'T_j (°C)', seq_cmap('thermal'));
panel_tag(ax, '(b) MOSFET junction temp.');
save_fig(fig9, figDir, 'Fig9_ZETA_Heatmaps');

% Fig 8 — cross-method comparison
if RUN_COMPARISON
    plot_method_comparison(ResultsByMethod.dpt, ResultsByMethod.transition, ...
                           ResultsByMethod.miller, figDir, 'ZETA');
end
end   % style loop

%% CSV exports --------------------------------------------------------------

SummaryTable.SwitchingLossMethod = repmat({PRIMARY_METHOD}, height(SummaryTable), 1);
writetable(SummaryTable,    'ZETA_v10_7_summary.csv');
writetable(ComponentLoss,   'ZETA_v10_7_component_losses.csv');
writetable(ComparisonReady, 'ZETA_v10_7_comparison_ready.csv');

fprintf('\nCSV files exported:\n');
fprintf('  ZETA_v10_7_summary.csv\n');
fprintf('  ZETA_v10_7_component_losses.csv\n');
fprintf('  ZETA_v10_7_comparison_ready.csv\n');
if RUN_COMPARISON
    fprintf('  ZETA_v10_7_method_comparison.csv\n');
end
fprintf('\nFigures exported to ./%s/<style>/ (ieee: PDF+EPS+600 dpi PNG; report: PDF+300 dpi PNG).\n', outDir);
fprintf('All figures and CSV files exported. Script completed.\n');

%% ========================================================================
%  Helper functions
%% ========================================================================

% The following are identical to the SEPIC file — copy straight across:
%   ieee_cfg, pick_font, apply_ieee_defaults, new_fig, ttl, save_fig,
%   write_one, panel_tag, seq_cmap, plot_heatmap, lx, lxf,
%   tex_escape, tex_sub, pad_lim, style_cases, dpt_params,
%   model_taxonomy, ternary
%
% Two functions differ for ZETA:
%   plot_method_comparison  (takes a 'tag' argument to name the figure)
%   calcZETA                (the ZETA-specific model, below)

function plot_method_comparison(R_dpt, R_trans, R_mill, outDir, tag)
% Same as SEPIC version except the output filename carries a tag.
    S = ieee_cfg();
    n = numel(R_dpt);
    labels = arrayfun(@(r) sprintf('%dV/%.1fA', r.Vin, r.Io), R_dpt, 'UniformOutput', false);

    c1 = [0.000 0.447 0.698]; c2 = [0.835 0.369 0.000]; c3 = [0.000 0.620 0.451];
    names = {'DPT (E_{on}+E_{off})', 'Transition + C_{oss}', 'Miller + C_{oss}'};

    Psw = [[R_dpt.Pmos_sw_comparable_dpt]', ...
           [R_trans.Pmos_sw_comparable_transition]', ...
           [R_mill.Pmos_sw_comparable_miller]'];
    Eff = [[R_dpt.eff_system]', [R_trans.eff_system]', [R_mill.eff_system]'];
    Tj  = [[R_dpt.Tj_Q1]',      [R_trans.Tj_Q1]',      [R_mill.Tj_Q1]'];
    clr = [c1; c2; c3];
    mk  = {'o','s','^'};

    Vin_all = [R_dpt.Vin];
    sepAt   = find(diff(Vin_all) ~= 0) + 0.5;

    fig = new_fig([tag '_MethodComparison'], S.col2, 2.9);
    tl  = tiledlayout(fig, 1, 3, 'TileSpacing','compact','Padding','compact');

    ax = nexttile(tl); hold(ax,'on');
    b = bar(ax, Psw, 'grouped', 'BarWidth', 0.85);
    for m = 1:3, b(m).FaceColor = clr(m,:); b(m).EdgeColor = 'none'; end
    ylabel(ax,lx('P_{sw} (W, C_{oss}-normalised)'));
    ylim(ax, [0, 1.60*max(Psw(:))]);
    legend(ax, lx(names), 'Location','north', 'FontSize', S.fontSmall);
    style_cases(ax, n, labels, sepAt); panel_tag(ax, '(a)');
    ax.XGrid = 'off';

    ax = nexttile(tl); hold(ax,'on');
    for m = 1:3
        plot(ax, 1:n, Eff(:,m), mk{m}, 'Color', clr(m,:), 'MarkerFaceColor', clr(m,:), ...
            'MarkerSize', S.ms + 1, 'LineStyle','none');
    end
    ylabel(ax,lx('System efficiency (%)'));
    ylim(ax, pad_lim(Eff(:), 0.20, false, 100));
    style_cases(ax, n, labels, sepAt); panel_tag(ax, '(b)');

    ax = nexttile(tl); hold(ax,'on');
    for m = 1:3
        plot(ax, 1:n, Tj(:,m), mk{m}, 'Color', clr(m,:), 'MarkerFaceColor', clr(m,:), ...
            'MarkerSize', S.ms + 1, 'LineStyle','none');
    end
    ylabel(ax,lx('MOSFET T_j (°C)'));
    ylim(ax, pad_lim(Tj(:), 0.20, false));
    style_cases(ax, n, labels, sepAt); panel_tag(ax, '(c)');

    ttl(tl, lx('Cross-method comparison (C_{oss}-normalised)'));
    save_fig(fig, outDir, ['Fig8_' tag '_MethodComparison']);
end


function res = calcZETA(Vin, Io, eta_est, cfg, L1_in, L2_in)
% ZETA steady-state at one operating point.
% Key differences from SEPIC:
%   - Coupling cap Cc sees a step current (Isw_step) instead of a fixed
%     charge per cycle -> its ripple calc changes.
%   - Output cap Cout sits on the L2 ripple, so dVo uses delta_IL2.
%   - Input cap Cin sees the full switch step (Isw_base + ripple).
%   - Optional coupled inductor (useCoupled) reduces L_eff by (1 - k_c).
%   - The input current is Iin_ideal / eta_est, not Iin * eta_est.

    assert(Vin > 0 && Io > 0);
    assert(eta_est > 0 && eta_est <= 1);

    Vo          = cfg.Vo;         fsw         = cfg.fsw;
    Vd          = cfg.Vd;         Rds_on      = cfg.Rds_on;
    ripple_factor = cfg.ripple_factor;
    rippleMethod  = cfg.rippleMethod;
    useCoupled    = cfg.useCoupled;   k_c = cfg.k_c;
    Vgate = cfg.Vgate;   Vgs3 = cfg.Vgs3;   Rg = cfg.Rg;
    Qgd   = cfg.Qgd;     Qg   = cfg.Qg;
    t_sw_rise = cfg.t_sw_rise;   t_sw_fall = cfg.t_sw_fall;
    switchingLossMethod = cfg.switchingLossMethod;
    Vds_max = cfg.Vds_max;   Id_max = cfg.Id_max;
    RJA_Q1 = cfg.RJA_Q1;     Tj_max_Q1 = cfg.Tj_max_Q1;
    RJA_D1 = cfg.RJA_D1;     Tj_max_D1 = cfg.Tj_max_D1;   T_amb = cfg.T_amb;
    DCR1 = cfg.DCR_L1;   DCR2 = cfg.DCR_L2;
    Cc_val = cfg.Cc_val;   Cout_val = cfg.Cout_val;   Cin_val = cfg.Cin_val;
    ESR_Cc = cfg.ESR_Cc;   ESR_Cout = cfg.ESR_Cout;   ESR_Cin = cfg.ESR_Cin;
    dVc_target  = cfg.dVc_target;
    dVo_target  = cfg.dVo_target;
    dVin_target = cfg.dVin_target;
    kE_on  = cfg.kE_on;
    kE_off = cfg.kE_off;

    if isfield(cfg, 'Tj_M_est'), Tj_M_est = cfg.Tj_M_est; else, Tj_M_est = T_amb; end

    Pout      = Vo * Io;
    Iin_ideal = Io * Vo / Vin;
    Iin_adj   = Iin_ideal / eta_est;

    % Duty-cycle iteration
    relax_D = 0.5;
    D = (Vo + Vd) / (Vin + Vo + Vd);
    duty_converged = false;
    for iter = 1:30
        Iin_est  = Io * Vo / (Vin * eta_est);
        Isw_avg  = Iin_est + Io;
        D_new    = (Vo + Vd + Io*DCR2) / ...
                   (Vin + Vo + Vd - Iin_est*DCR1 - Isw_avg*Rds_on);
        D_next   = relax_D * D_new + (1 - relax_D) * D;
        if abs(D_next - D) < 1e-7
            D = D_next;
            duty_converged = true;
            break;
        end
        D = D_next;
    end
    if ~duty_converged
        warning('Duty-cycle iteration did not converge at Vin=%.2f V, Io=%.2f A', Vin, Io);
    end
    if ~(D > 0 && D < 1), error('Duty cycle out of range: D=%.4f', D); end

    L_req  = (Vin * D) / (Iin_adj * ripple_factor * fsw);
    L2_req = (Vin * D) / (Io      * ripple_factor * fsw);

    if L1_in > 0, L1_use = L1_in; else, L1_use = L_req;  end
    if L2_in > 0, L2_use = L2_in; else, L2_use = L2_req; end

    if useCoupled
        L_eff     = L1_use * (1 - k_c);
        delta_IL1 = (Vin * D) / (L_eff * fsw);
        delta_IL2 = delta_IL1;
    else
        delta_IL1 = (Vin * D) / (L1_use * fsw);
        delta_IL2 = (Vin * D) / (L2_use * fsw);
    end

    IL1_peak = Iin_adj + delta_IL1/2;
    IL2_peak = Io      + delta_IL2/2;
    IL1_min  = Iin_adj - delta_IL1/2;
    IL2_min  = Io      - delta_IL2/2;
    IL1_rms  = sqrt(Iin_adj^2 + delta_IL1^2/12);
    IL2_rms  = sqrt(Io^2      + delta_IL2^2/12);

    Isw_base   = Iin_adj + Io;
    Isw_peak   = Isw_base + (delta_IL1 + delta_IL2) / 2;
    Isw_valley = max(Isw_base - (delta_IL1 + delta_IL2) / 2, 0);
    Isw_rms    = sqrt(D) * sqrt(Isw_base^2 + (delta_IL1 + delta_IL2)^2/12);

    Vsw_theory = Vin + Vo + Vd;
    Vsw_rating = 1.30 * Vsw_theory;
    Vsw_status = ternary(Vsw_rating <= Vds_max, 'PASS', 'FAIL Vds');

    Isw_step = Iin_adj + Io;

    % Ripple
    dVc_cap    = (D * Io) / (Cc_val * fsw);
    dVc_esr    = Isw_step * ESR_Cc;
    dVc_actual = dVc_cap + dVc_esr;

    dVo_cap    = delta_IL2 / (8 * Cout_val * fsw);
    dVo_esr    = delta_IL2 * ESR_Cout;
    dVo_actual = dVo_cap + dVo_esr;

    dVin_cap    = (D * Io) / (Cin_val * fsw);
    dVin_esr    = Isw_step * ESR_Cin;
    dVin_actual = dVin_cap + dVin_esr;

    Cc_status   = ternary(dVc_actual  <= dVc_target,  'PASS', 'FAIL');
    Cout_status = ternary(dVo_actual  <= dVo_target,  'PASS', 'FAIL');
    Cin_status  = ternary(dVin_actual <= dVin_target, 'PASS', 'FAIL');

    % RMS currents
    ICc_rms    = sqrt( D * (Io^2 + delta_IL2^2/12) + (1-D) * (Iin_adj^2 + delta_IL1^2/12) );
    ICout_rms  = delta_IL2 / sqrt(12);
    delta_Itotal = delta_IL1 + delta_IL2;
    ICin_rms   = sqrt( D * (Isw_base^2 + delta_Itotal^2/12) + (1-D) * Iin_adj^2 );

    % Switching-loss calculations
    Pmos_cond = Isw_rms^2 * Rds_on;

    Pmos_sw_overlap_transition = 0.5 * Vsw_theory * ...
        (Isw_valley + Isw_peak) * (t_sw_rise + t_sw_fall) * fsw;

    Ig_on  = max((Vgate - Vgs3)/Rg, 0.01);
    Ig_off = max( Vgs3        /Rg, 0.01);
    t_miller_on  = Qgd / Ig_on;
    t_miller_off = Qgd / Ig_off;
    Pmos_sw_overlap_miller = 0.5 * Vsw_theory * Isw_peak * ...
        (t_miller_on + t_miller_off) * fsw;

    dpt = dpt_params();
    if numel(dpt.V_ref) > 1
        E_on_ref_v  = interp1(dpt.V_ref, dpt.E_on_ref,  Vsw_theory, 'linear','extrap');
        E_off_ref_v = interp1(dpt.V_ref, dpt.E_off_ref, Vsw_theory, 'linear','extrap');
        I_on_ref_v  = interp1(dpt.V_ref, dpt.I_on_ref,  Vsw_theory, 'linear','extrap');
        I_off_ref_v = interp1(dpt.V_ref, dpt.I_off_ref, Vsw_theory, 'linear','extrap');
    else
        E_on_ref_v  = dpt.E_on_ref  * (Vsw_theory / dpt.V_ref);
        E_off_ref_v = dpt.E_off_ref * (Vsw_theory / dpt.V_ref);
        I_on_ref_v  = dpt.I_on_ref(1);
        I_off_ref_v = dpt.I_off_ref(1);
    end

    I_on_ref_v  = max(I_on_ref_v,  0.01);
    I_off_ref_v = max(I_off_ref_v, 0.01);

    E_on  = E_on_ref_v  * (Isw_valley / I_on_ref_v);
    E_off = E_off_ref_v * (Isw_peak   / I_off_ref_v);

    E_on  = E_on  * (1 + kE_on  * (Tj_M_est - 25));
    E_off = E_off * (1 + kE_off * (Tj_M_est - 25));

    Pmos_sw_overlap_dpt = (E_on + E_off) * fsw;

    switch lower(strtrim(switchingLossMethod))
        case 'transition', Pmos_sw_overlap = Pmos_sw_overlap_transition;
        case 'miller',     Pmos_sw_overlap = Pmos_sw_overlap_miller;
        case 'dpt',        Pmos_sw_overlap = Pmos_sw_overlap_dpt;
        otherwise
            error('Invalid switchingLossMethod. Use ''dpt'', ''transition'' or ''miller''.');
    end

    % Comparable-on-equal-footing variants
    Pmos_sw_gate     = Vgate * Qg * fsw;
    Pmos_Coss_actual = 0.5 * cfg.Coss * Vsw_theory^2 * fsw;
    Pmos_Coss        = Pmos_Coss_actual;
    if strcmpi(strtrim(switchingLossMethod), 'dpt')
        Pmos_Coss = 0;
    end

    Pmos_sw_comparable_dpt        = Pmos_sw_overlap_dpt;
    Pmos_sw_comparable_transition = Pmos_sw_overlap_transition + Pmos_Coss_actual;
    Pmos_sw_comparable_miller     = Pmos_sw_overlap_miller     + Pmos_Coss_actual;

    Pmos_sw    = Pmos_sw_overlap + Pmos_sw_gate + Pmos_Coss;
    Pmos_total = Pmos_cond + Pmos_sw;

    % Diode
    Pdiode_cond = Vd * Io;
    Pdiode_cap  = 0.5 * cfg.Cj_diode * Vsw_theory^2 * fsw;
    Pdiode      = Pdiode_cond + Pdiode_cap;

    % Inductors
    PL1_copper = IL1_rms^2 * DCR1;
    PL2_copper = IL2_rms^2 * DCR2;
    Pind_total = PL1_copper + PL2_copper;

    % Capacitors
    Pcap_Cc    = ICc_rms^2   * ESR_Cc;
    Pcap_Cout  = ICout_rms^2 * ESR_Cout;
    Pcap_Cin   = ICin_rms^2  * ESR_Cin;
    Pcap_total = Pcap_Cc + Pcap_Cout + Pcap_Cin;

    % Roll-ups
    Pmos_stage        = Pmos_cond + Pmos_sw_overlap + Pmos_Coss;
    Pmos_Coss_miller  = 0.5 * cfg.Coss * Vsw_theory^2 * fsw;
    Pmos_stage_miller = Pmos_cond + Pmos_sw_overlap_miller + Pmos_Coss_miller;

    Ploss_stage         = Pmos_stage         + Pdiode + Pind_total + Pcap_total;
    Ploss_stage_miller  = Pmos_stage_miller  + Pdiode + Pind_total + Pcap_total;
    Ploss_system        = Ploss_stage        + Pmos_sw_gate;
    Ploss_system_miller = Ploss_stage_miller + Pmos_sw_gate;

    Pin_stage_est         = Pout + Ploss_stage;
    Pin_system_est        = Pout + Ploss_system;
    Pin_stage_est_miller  = Pout + Ploss_stage_miller;
    Pin_system_est_miller = Pout + Ploss_system_miller;

    eff_stage         = 100 * Pout / Pin_stage_est;
    eff_system        = 100 * Pout / Pin_system_est;
    eff_stage_miller  = 100 * Pout / Pin_stage_est_miller;
    eff_system_miller = 100 * Pout / Pin_system_est_miller;

    % Thermal
    Tj_Q1 = T_amb + Pmos_stage * RJA_Q1;
    Tj_D1 = T_amb + Pdiode     * RJA_D1;
    MOSFET_thermal_status = ternary(Tj_Q1 < Tj_max_Q1, 'OK', 'WARNING');
    Diode_thermal_status  = ternary(Tj_D1 < Tj_max_D1, 'OK', 'WARNING');

    if Ploss_system > 0
        Pmos_pct   = 100 * Pmos_total  / Ploss_system;
        Pdiode_pct = 100 * Pdiode      / Ploss_system;
        Pind_pct   = 100 * Pind_total  / Ploss_system;
        Pcap_pct   = 100 * Pcap_total  / Ploss_system;
    else
        Pmos_pct = 0; Pdiode_pct = 0; Pind_pct = 0; Pcap_pct = 0;
    end

    L1_crit = (Vin * D) / (2 * fsw * Iin_adj);
    L2_crit = (Vin * D) / (2 * fsw * Io);

    if IL1_min > 0 && IL2_min > 0, CCM_status = 'CCM OK';
    else, CCM_status = 'DCM possible'; end

    res.Vin = Vin;   res.Io = Io;   res.Pout = Pout;   res.D = D;
    res.Iin_ideal = Iin_ideal;   res.Iin_adj = Iin_adj;
    res.delta_IL1 = delta_IL1;   res.delta_IL2 = delta_IL2;
    res.L_req = L_req;   res.L2_req = L2_req;
    res.L1_val = L1_use; res.L2_val = L2_use;
    res.L1_crit = L1_crit;   res.L2_crit = L2_crit;
    res.IL1_peak = IL1_peak; res.IL2_peak = IL2_peak;
    res.IL1_min  = IL1_min;  res.IL2_min  = IL2_min;
    res.IL1_rms  = IL1_rms;  res.IL2_rms  = IL2_rms;
    res.Isw_peak = Isw_peak; res.Isw_valley = Isw_valley; res.Isw_rms  = Isw_rms;
    res.E_off_J  = E_off;    res.E_on_J = E_on;
    res.Vsw_theory = Vsw_theory;   res.Vsw_rating = Vsw_rating;
    res.Vsw_status = Vsw_status;
    res.dVc_actual  = dVc_actual;
    res.dVo_actual  = dVo_actual;
    res.dVin_actual = dVin_actual;
    res.Cc_status   = Cc_status;
    res.Cout_status = Cout_status;
    res.Cin_status  = Cin_status;
    res.ICc_rms     = ICc_rms;
    res.ICout_rms   = ICout_rms;
    res.ICin_rms    = ICin_rms;
    res.PL1_copper  = PL1_copper;   res.PL2_copper  = PL2_copper;
    res.Pind_total  = Pind_total;
    res.Pmos_cond   = Pmos_cond;
    res.Pmos_sw_overlap            = Pmos_sw_overlap;
    res.Pmos_sw_overlap_miller     = Pmos_sw_overlap_miller;
    res.Pmos_sw_overlap_transition = Pmos_sw_overlap_transition;
    res.Pmos_sw_overlap_dpt        = Pmos_sw_overlap_dpt;
    res.Pmos_sw_gate  = Pmos_sw_gate;
    res.Pmos_Coss     = Pmos_Coss;
    res.Pmos_Coss_actual = Pmos_Coss_actual;
    res.Pmos_sw_comparable_dpt        = Pmos_sw_comparable_dpt;
    res.Pmos_sw_comparable_transition = Pmos_sw_comparable_transition;
    res.Pmos_sw_comparable_miller     = Pmos_sw_comparable_miller;
    res.Pmos_sw       = Pmos_sw;
    res.Pmos_total    = Pmos_total;
    res.Pdiode_cond   = Pdiode_cond;
    res.Pdiode_cap    = Pdiode_cap;
    res.Pdiode        = Pdiode;
    res.Pcap_Cc       = Pcap_Cc;
    res.Pcap_Cout     = Pcap_Cout;
    res.Pcap_Cin      = Pcap_Cin;
    res.Pcap_total    = Pcap_total;
    res.Pmos_stage        = Pmos_stage;
    res.Pmos_stage_miller = Pmos_stage_miller;
    res.Ploss_stage         = Ploss_stage;
    res.Ploss_system        = Ploss_system;
    res.Ploss_stage_miller  = Ploss_stage_miller;
    res.Ploss_system_miller = Ploss_system_miller;
    res.Pin_stage_est         = Pin_stage_est;
    res.Pin_system_est        = Pin_system_est;
    res.Pin_stage_est_miller  = Pin_stage_est_miller;
    res.Pin_system_est_miller = Pin_system_est_miller;
    res.eff_stage         = eff_stage;
    res.eff_system        = eff_system;
    res.eff_stage_miller  = eff_stage_miller;
    res.eff_system_miller = eff_system_miller;
    res.Ploss_total = Ploss_system;
    res.Pin_est     = Pin_system_est;
    res.eff         = eff_system;
    res.Pmos_pct    = Pmos_pct;
    res.Pdiode_pct  = Pdiode_pct;
    res.Pind_pct    = Pind_pct;
    res.Pcap_pct    = Pcap_pct;
    res.CCM_status  = CCM_status;
    res.Tj_Q1       = Tj_Q1;   res.Tj_D1 = Tj_D1;
    res.MOSFET_thermal_status = MOSFET_thermal_status;
    res.Diode_thermal_status  = Diode_thermal_status;
    res.Rds_on_final   = Rds_on;
    res.duty_converged = duty_converged;
    res.t_miller_on_s  = t_miller_on;
    res.t_miller_off_s = t_miller_off;
    res.Ig_on_A        = Ig_on;
    res.Ig_off_A       = Ig_off;
    res.Tj_M_est_used  = Tj_M_est;
    res.I_on_ref_v     = I_on_ref_v;
    res.I_off_ref_v    = I_off_ref_v;
end


function cfg = build_cfg(Vo, fsw, Vd_25, temp_coeff_diode, ...
        ripple_factor, rippleMethod, useCoupled, k_c, ...
        Rds_on_25, alpha_Rds, Coss, Qgd, Qg, ...
        Vgate, Vgs3, Rg, t_sw_rise, t_sw_fall, switchingLossMethod, ...
        Vds_max, Id_max, kE_on, kE_off, ...
        RJA_Q1, Tj_max_Q1, Cj_diode, Vrrm_diode, Ifavg_diode, RJA_D1, Tj_max_D1, T_amb, ...
        DCR_L1, DCR_L2, Cc_val, Cout_val, Cin_val, ESR_Cc, ESR_Cout, ESR_Cin, ...
        dVc_target, dVo_target, dVin_target)
% Same as SEPIC's build_cfg but with the ZETA-specific coupling fields.
    cfg.Vo = Vo;   cfg.fsw = fsw;
    cfg.Vd_25 = Vd_25; cfg.Vd = Vd_25;
    cfg.temp_coeff_diode = temp_coeff_diode;
    cfg.ripple_factor = ripple_factor; cfg.rippleMethod = rippleMethod;
    cfg.useCoupled = useCoupled; cfg.k_c = k_c;
    cfg.Rds_on_25 = Rds_on_25; cfg.Rds_on = Rds_on_25; cfg.alpha_Rds = alpha_Rds;
    cfg.Coss = Coss; cfg.Qgd = Qgd; cfg.Qg = Qg;
    cfg.Vgate = Vgate; cfg.Vgs3 = Vgs3; cfg.Rg = Rg;
    cfg.t_sw_rise = t_sw_rise; cfg.t_sw_fall = t_sw_fall;
    cfg.switchingLossMethod = switchingLossMethod;
    cfg.Vds_max = Vds_max; cfg.Id_max = Id_max;
    cfg.kE_on = kE_on; cfg.kE_off = kE_off;
    cfg.RJA_Q1 = RJA_Q1; cfg.Tj_max_Q1 = Tj_max_Q1;
    cfg.Cj_diode = Cj_diode; cfg.Vrrm_diode = Vrrm_diode; cfg.Ifavg_diode = Ifavg_diode;
    cfg.RJA_D1 = RJA_D1; cfg.Tj_max_D1 = Tj_max_D1; cfg.T_amb = T_amb;
    cfg.DCR_L1 = DCR_L1; cfg.DCR_L2 = DCR_L2;
    cfg.Cc_val = Cc_val; cfg.Cout_val = Cout_val; cfg.Cin_val = Cin_val;
    cfg.ESR_Cc = ESR_Cc; cfg.ESR_Cout = ESR_Cout; cfg.ESR_Cin = ESR_Cin;
    cfg.dVc_target = dVc_target; cfg.dVo_target = dVo_target; cfg.dVin_target = dVin_target;
    cfg.Tj_M_est = T_amb; cfg.Tj_D_est = T_amb;
end