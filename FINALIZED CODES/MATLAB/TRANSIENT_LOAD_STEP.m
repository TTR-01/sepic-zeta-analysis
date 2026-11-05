%% transientloadstep.m
% MOSFET thermal response to a load step, plotted on twin axes.
%   left  axis: total junction temperature Tj (set by Rth_CA, the slow path)
%   right axis: die-level Foster rise dT_die (fast thermal inertia)
%
% Runs a 500 ms window with the load stepping at 100 ms, then writes
% Fig9 as both vector PDF and 300 dpi PNG.

t_end  = 0.5;          % s
dt     = 1e-4;         % 0.1 ms sample
t      = (0:dt:t_end)';
t_step = 0.1;          % load changes at 100 ms

% MOSFET die losses only (conduction + switching overlap + Coss).
% Case 1 = 9 V / 0.5 A, Case 3 = 9 V / 1.5 A.
P_low  = 0.233;        % W
P_high = 0.945;        % W

P_loss = P_low * (t < t_step) + P_high * (t >= t_step);

% Foster die-to-case network plus a single case-to-ambient resistance.
[Tj, dT_die] = thermal_foster(t, P_loss, 40, 61);

%% Plot ----------------------------------------------------------------------

fig = figure('Name','Foster_Thermal_Transient','NumberTitle','off', ...
    'Position', [200 200 950 480]);

% Total Tj on the left. Limits are fixed so the step is visible against
% the 40 C ambient baseline.
yyaxis left
plot(t*1e3, Tj, 'LineWidth', 2.0, 'Color', [0.00 0.45 0.70]);
ylabel('Junction temperature T_j (°C)');
ylim([50 105]);
yticks(50:5:105);

% Die-level transient on the right. Auto-scaled to the data because the
% amplitude depends on the loss step and the Foster coefficients.
yyaxis right
plot(t*1e3, dT_die, '--', 'LineWidth', 2.0, 'Color', [0.84 0.37 0.00]);
ylabel('Die-level rise \DeltaT_{Foster} (°C)');
ylim([0 max(dT_die)*1.15 + 0.05]);

xline(t_step*1e3, ':', 'Load step 0.5 → 1.5 A at 100 ms', ...
    'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');

xlabel('Time (ms)');
title('MOSFET thermal response: total T_j (left), die-level \DeltaT (right)');
grid on; box on;
xlim([0 500]);

%% Export --------------------------------------------------------------------

outDir = 'report_figures';
if ~isfolder(outDir)
    mkdir(outDir);
end

pdf_path = fullfile(outDir, 'Fig9_Foster_LoadStep.pdf');
png_path = fullfile(outDir, 'Fig9_Foster_LoadStep.png');

try
    % R2020a+: keeps the PDF as vectors and the PNG at 300 dpi.
    exportgraphics(fig, pdf_path, 'ContentType', 'vector');
    exportgraphics(fig, png_path, 'Resolution', 300);
    fprintf('\nPDF saved to: %s\n', pdf_path);
catch
    % Older releases fall back to print().
    print(fig, pdf_path, '-dpdf', '-bestfit');
    print(fig, png_path, '-dpng', '-r300');
    fprintf('\nPDF saved to: %s (via print fallback)\n', pdf_path);
end

%% Console summary -----------------------------------------------------------

% The step happens at 200 ms (100 ms in, then settle to end of window).
% Grab the instantaneous value at the step and the settled mean afterwards
% so the printout shows both the peak and the long-term tail.
dT_at_step  = interp1(t, dT_die, 0.200);
dT_settled  = mean(dT_die(t >= 0.400));

fprintf('\n--- Foster thermal step response ---\n');
fprintf('  Tj initial         = %.2f °C\n', Tj(1));
fprintf('  Tj final           = %.2f °C\n', Tj(end));
fprintf('  dT_die initial     = %.4f °C\n', dT_die(1));
fprintf('  dT_die at 200 ms   = %.4f °C\n', dT_at_step);
fprintf('  dT_die settled     = %.4f °C\n', dT_settled);
fprintf('  dT_die final       = %.4f °C\n', dT_die(end));


%% ------------------------------------------------------------------------
function [Tj, dT_die] = thermal_foster(t, P_loss, T_amb, Rth_CA)
% Four-stage Foster network for die-to-case, plus one lumped case-to-ambient
% resistance. Each stage is a first-order low-pass in thermal domain; the
% trapezoidal input samples the loss at the mid-point between steps so a
% step change doesn't get charged as a full-step impulse.

    if nargin < 3 || isempty(T_amb),  T_amb  = 40; end
    if nargin < 4 || isempty(Rth_CA), Rth_CA = 61; end

    % Datasheet Foster coefficients for the IRLZ44N package
    R   = [206.89, 395.98, 227.36, 170.35] * 1e-3;   % °C/W
    C   = [222.07,  62.26, 371.75,   5.477] * 1e-3;  % J/°C
    tau = R .* C;                                    % s

    t      = t(:);
    P_loss = P_loss(:);
    n      = numel(t);

    dT_die = zeros(n, 1);
    for i = 1:4
        dT_i = zeros(n, 1);
        for k = 2:n
            dt_k  = t(k) - t(k-1);
            P_mid = 0.5 * (P_loss(k) + P_loss(k-1));
            a     = exp(-dt_k / tau(i));
            dT_i(k) = dT_i(k-1) * a + R(i) * P_mid * (1 - a);
        end
        dT_die = dT_die + dT_i;
    end

    % Steady-state contribution from the slow path plus the die-level
    % transient. This is why the left axis barely moves at the step while
    % the right axis spikes and decays.
    Tj = T_amb + Rth_CA * P_loss + dT_die;
end