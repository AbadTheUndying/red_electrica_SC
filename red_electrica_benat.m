clear; 
clc;
close all;

%% CONSTRUIR GRAFO %%

% Nodos
nodeNames = {
'B01 Galicia-Norte'
'B02 Galicia-Sur'
'B03 Asturias-Centro'
'B04 Cantabria-Mataporquera'
'B05 Pais Vasco-Gatika'
'B06 Navarra-Muruarte'
'B07 Aragon-Ebro'
'B08 Aragon-Pirineos'
'B09 CastillaLeon-Oeste'
'B10 CastillaLeon-Este'
'B11 Extremadura-Almaraz'
'B12 Extremadura-Badajoz'
'B13 CLM-Oeste'
'B14 CLM-Este-Slack'
'B15 Madrid-400kV'
'B16 Andalucia-Oeste'
'B17 Andalucia-Este'
'B18 Cataluna-Tarragona'
'B19 Cataluna-Norte'
'B20 CValenciana-Saguntum'
};

% Aristas propuestas
s = [
    1 1 2 3 4 5 6 7 8 ...
    5 6 9 9 9 10 ...
    11 11 13 14 ...
    12 13 13 16 ...
    7 18 18 15 14 17 ...
];

t = [ ...
    2 3 9 4 5 6 7 8 19 ...
    10 10 10 11 13 14 ...
    12 13 14 15 ...
    16 16 17 17 ...
    18 19 20 20 20 20 ...
];

G = graph(s,t,[],nodeNames);


% Coordenadas
lon = [ ...
   -7.5  ... B01 Galicia-Norte
   -8.4  ... B02 Galicia-Sur
   -5.8  ... B03 Asturias-Centro
   -4.0  ... B04 Cantabria-Mataporquera
   -2.9  ... B05 Pais Vasco-Gatika
   -1.6  ... B06 Navarra-Muruarte
   -0.9  ... B07 Aragon-Ebro
    0.6  ... B08 Aragon-Pirineos
   -5.7  ... B09 CastillaLeon-Oeste
   -4.7  ... B10 CastillaLeon-Este
   -5.7  ... B11 Extremadura-Almaraz
   -6.8  ... B12 Extremadura-Badajoz
   -4.0  ... B13 CLM-Oeste
   -3.0  ... B14 CLM-Este-Slack
   -3.7  ... B15 Madrid-400kV
   -6.0  ... B16 Andalucia-Oeste
   -3.2  ... B17 Andalucia-Este
    1.2  ... B18 Cataluna-Tarragona
    2.2  ... B19 Cataluna-Norte
   -0.2  ... B20 CValenciana-Saguntum
];

lat = [ ...
   43.2  ... B01 Galicia-Norte
   42.3  ... B02 Galicia-Sur
   43.4  ... B03 Asturias-Centro
   43.0  ... B04 Cantabria-Mataporquera
   43.3  ... B05 Pais Vasco-Gatika
   42.7  ... B06 Navarra-Muruarte
   41.6  ... B07 Aragon-Ebro
   42.7  ... B08 Aragon-Pirineos
   41.8  ... B09 CastillaLeon-Oeste
   41.8  ... B10 CastillaLeon-Este
   39.8  ... B11 Extremadura-Almaraz
   38.9  ... B12 Extremadura-Badajoz
   38.9  ... B13 CLM-Oeste
   39.8  ... B14 CLM-Este-Slack
   40.4  ... B15 Madrid-400kV
   37.5  ... B16 Andalucia-Oeste
   37.2  ... B17 Andalucia-Este
   41.1  ... B18 Cataluna-Tarragona
   42.1  ... B19 Cataluna-Norte
   39.7  ... B20 CValenciana-Saguntum
];
  
% Representación
figure();
p = plot(G, ...
    'XData', lon, ...
    'YData', lat, ...
    'NodeLabel', nodeNames, ...
    'LineWidth', 1.2, ...
    'MarkerSize', 5);

title('Backbone red peninsular con coordenadas geograficas aproximadas');
xlabel('Longitud');
ylabel('Latitud');
axis equal;
grid on;

% Matrices características del grafo
A = full(adjacency(G));
L = full(laplacian(G));
B = full(incidence(G));


%% MODELO SWING EQUATION %%
% =========================================================================
% Modelos implementados (Dörfler & Bullo, 2012 / Swing Equation notes):
%
%  [SWING]    M_i·θ̈_i = -D_i·θ̇_i + P_i - Σ_j P_ij·sin(θ_i - θ_j)
%
%  [KURAMOTO] D_i·θ̇_i =              P_i - Σ_j P_ij·sin(θ_i - θ_j)
%
%  Unidades:  P [GW],  M_i [GW·s²],  D_i [GW·s],  P_ij [GW],  θ [rad]
%  Red sin pérdidas (lossless): φ_ij = 0  →  solo susceptancias
% =========================================================================

%% -- PARÁMETROS FÍSICOS --------------------------------------------------

N   = numnodes(G);
f0  = 50;               % Hz  (frecuencia nominal UCTE)
w0  = 2*pi*f0;          % rad/s

% -- POTENCIAS Pi [GW] (inyección neta media)
% positivo = nodo exportador / generador neto
% negativo = nodo importador / consumidor neto
Pi = [ ...
   0.623  % B01 Galicia-Norte
   0.265  % B02 Galicia-Sur
  -0.093  % B03 Asturias-Centro
  -0.367  % B04 Cantabria-Mataporquera
  -1.096  % B05 Pais Vasco-Gatika
   0.166  % B06 Navarra-Muruarte
   1.007  % B07 Aragon-Ebro
   0.140  % B08 Aragon-Pirineos
   0.752  % B09 CastillaLeon-Oeste
   0.589  % B10 CastillaLeon-Este
   2.186  % B11 Extremadura-Almaraz
   0.350  % B12 Extremadura-Badajoz
   0.179  % B13 CLM-Oeste
   1.622  % B14 CLM-Este-Slack
  -2.999  % B15 Madrid-400kV
   0.098  % B16 Andalucia-Oeste
  -1.038  % B17 Andalucia-Este
   1.811  % B18 Cataluna-Tarragona
  -2.893  % B19 Cataluna-Norte
  -1.302  % B20 CValenciana-Saguntum
];

sum(Pi)   % ~ 0

% -- INERCIAS Mi [GW·s²]
%    Nodos grandes tienen más inercia
Mi = [
   0.12   % B01 Galicia-Norte       ~eólica + hidro
   0.10   % B02 Galicia-Sur         ~eólica + hidro
   0.08   % B03 Asturias-Centro     ~carbón
   0.07   % B04 Cantabria-Matap.    ~hidro
   0.15   % B05 Pais Vasco-Gatika   ~ciclos combinados
   0.09   % B06 Navarra-Muruarte    ~eólica + ciclos
   0.13   % B07 Aragon-Ebro         ~eólica + hidro
   0.07   % B08 Aragon-Pirineos     ~hidro
   0.16   % B09 CastillaLeon-Oeste  ~eólica + ciclos
   0.14   % B10 CastillaLeon-Este   ~ciclos combinados
   0.28   % B11 Extremadura-Almaraz ~NUCLEAR (alta inercia)
   0.09   % B12 Extremadura-Badajoz ~renovables
   0.08   % B13 CLM-Oeste           ~solar
   0.20   % B14 CLM-Este-Slack      ~nuclear + ciclos
   0.30   % B15 Madrid-400kV        ~gran carga + ciclos
   0.10   % B16 Andalucia-Oeste     ~solar + eólica
   0.13   % B17 Andalucia-Este      ~solar + ciclos
   0.22   % B18 Cataluna-Tarragona  ~NUCLEAR + CC
   0.16   % B19 Cataluna-Norte      ~hidro + ciclos
   0.14   % B20 CValenciana-Sagunt. ~ciclos + nuclear
];

% -- Amortiguamiento Di [GW·s]
%    Ratio ε = M_max/D_min << 1 valida la aproximación Kuramoto
Di = 3 * Mi;

fprintf('Ratio ε = Mmax/Dmin = %.3f  (< 1 valida aprox. Kuramoto)\n', ...
        max(Mi)/min(Di));

% ── Capacidades máximas de transferencia Pij_max [GW]
%    (P_ij = E_i·E_j·|Y_ij|, con V=1 pu, 400 kV doble circuito típico)
%    Mismo orden que las aristas s, t del grafo
Pij_max = [
   1.5   % B01-B02  Galicia Norte-Sur
   1.2   % B01-B03  Galicia-Asturias
   1.8   % B02-B09  Galicia Sur-CyL Oeste
   1.0   % B03-B04  Asturias-Cantabria
   1.2   % B04-B05  Cantabria-P.Vasco
   1.5   % B05-B06  P.Vasco-Navarra
   1.5   % B06-B07  Navarra-Aragon Ebro
   1.0   % B07-B08  Aragon Ebro-Pirineos
   1.2   % B08-B19  Pirineos-Cat.Norte
   2.0   % B05-B10  P.Vasco-CyL Este
   1.8   % B06-B10  Navarra-CyL Este
   2.5   % B09-B10  CyL Oeste-Este  (corredor N-S clave)
   2.0   % B09-B11  CyL Oeste-Extremadura Almaraz
   1.8   % B09-B13  CyL Oeste-CLM Oeste
   2.5   % B10-B14  CyL Este-CLM Este  (corredor clave)
   1.5   % B11-B12  Extr.Almaraz-Badajoz
   2.0   % B11-B13  Extr.Almaraz-CLM Oeste
   2.0   % B13-B14  CLM Oeste-CLM Este
   3.0   % B14-B15  CLM Este-Madrid  (corredor principal)
   1.5   % B12-B16  Extr.Badajoz-And.Oeste
   2.0   % B13-B16  CLM Oeste-And.Oeste
   1.5   % B13-B17  CLM Oeste-And.Este
   1.5   % B16-B17  And.Oeste-Este
   2.0   % B07-B18  Aragon-Cat.Tarragona
   1.5   % B18-B19  Cat.Tarragona-Norte
   1.5   % B18-B20  Cat.Tarragona-C.Valenciana
   3.0   % B15-B20  Madrid-C.Valenciana  (corredor principal)
   2.5   % B14-B20  CLM Este-C.Valenciana
   1.5   % B17-B20  And.Este-C.Valenciana
];

% ── Matriz de acoplamiento K [N×N]  (K_ij = Pij_max para aristas, 0 resto)
edges_mat = [s(:), t(:)];
Ne = size(edges_mat, 1);
K  = zeros(N, N);
for e = 1:Ne
    ii = edges_mat(e,1);  jj = edges_mat(e,2);
    K(ii,jj) = Pij_max(e);
    K(jj,ii) = Pij_max(e);   % red sin pérdidas → simétrica
end

% ── Laplaciano ponderado por capacidades (para análisis espectral)
Lw = diag(sum(K,2)) - K;
ev = eig(Lw);
ev_sorted = sort(real(ev));
fprintf('Conectividad algebraica λ₂(L_w) = %.4f GW\n', ev_sorted(2));
fprintf('(λ₂ > 0 → grafo conectado; mayor valor → red más robusta)\n\n');


%% -- PUNTO DE OPERACIÓN ESTACIONARIO ────────────────────────────────
%  Resolver flujo de potencia DC: P_i = Σ_j K_ij · sin(δ_i - δ_j)
%  B14 (nodo 14) como referencia de ángulo: δ_14 = 0

slack     = 14; % índice del nodo de referencia 
non_slack = setdiff(1:N, slack);

% Punto inicial: ángulos pequeños proporcionales a la potencia
delta_guess           = zeros(N,1); % vector de angulos iniciales
delta_guess(non_slack) = Pi(non_slack) * 0.05;  % escalar para aproximar por angulos pequeños

opts_fs = optimoptions('fsolve', 'Display','off', ...
                        'TolFun',1e-10, 'TolX',1e-10, ...
                        'MaxIter',2000, 'Algorithm','levenberg-marquardt'); % definicion de opciones para fsolve

[delta_free, ~, exitflag] = fsolve( ...
    @(x) pf_mismatch(x, Pi, K, slack, N), ...
    delta_guess(non_slack), opts_fs);   % resolucion del sistema no lineal

if exitflag <= 0
    warning('fsolve: convergencia dudosa. Comprobar factibilidad del flujo.');
end

delta_op          = zeros(N,1);
delta_op(non_slack) = delta_free;   % insertar las soluciones

% Verificar residuos
P_calc    = K * sin(delta_op - delta_op.');  % suma sobre j
P_net     = sum(P_calc, 2);
residuals = P_net - Pi;
fprintf('─── Flujo de potencia estacionario ───\n');
fprintf('Residuo máximo: %.2e GW  (< 1e-6 → convergido)\n\n', max(abs(residuals)));

fprintf('%-28s  %8s  %8s\n','Nodo','δ [°]','P_calc [GW]');
fprintf('%-28s  %8s  %8s\n','----','------','----------');
for i = 1:N
    fprintf('%-28s  %+8.3f  %+8.3f\n', nodeNames{i}, ...
            rad2deg(delta_op(i)), P_net(i));
end
fprintf('\n');

%% ── 3. CONFIGURACIÓN DE PERTURBACIÓN (N-1) ────────────────────────────
%  Escenario: pérdida repentina de B11 Extremadura-Almaraz (nuclear, 2.186 GW)
%  en t = t_fault. Efecto más severo porque es el mayor generador neto.

t_fault = 2.0;    % s  (fallo en t=2 s)
t_end   = 40.0;   % s  (simulación hasta 40 s)

Pi_fault     = Pi;
Pi_fault(11) = 0;   % pérdida de generación en B11
% El desequilibrio (2.186 GW) se redistribuye dinámicamente por la red

fprintf('─── Contingencia simulada ───\n');
fprintf('Fallo en: %s\n', nodeNames{11});
fprintf('Potencia perdida: %.3f GW  en t = %.1f s\n\n', Pi(11), t_fault);

%% ── 4. SIMULACIÓN: SWING EQUATION ────────────────────────────────────
%  Estado: x = [θ₁…θ_N, θ̇₁…θ̇_N]'   (2N variables)
%  dθ_i/dt   = ω_i
%  dω_i/dt   = (P_i - D_i·ω_i - Σ_j K_ij·sin(θ_i-θ_j)) / M_i

x0_sw = [delta_op; zeros(N,1)];   % estado inicial: punto estacionario
ode_opts = odeset('RelTol',1e-7, 'AbsTol',1e-9, 'MaxStep',0.05);

fprintf('Simulando Swing Equation... ');
% Fase pre-fallo
[t_sw1, X_sw1] = ode15s(@(t,x) swing_eq(t,x,Pi,    K,Mi,Di,N), ...
                          [0, t_fault], x0_sw, ode_opts);
% Fase post-fallo (condiciones iniciales desde el final del pre-fallo)
[t_sw2, X_sw2] = ode15s(@(t,x) swing_eq(t,x,Pi_fault,K,Mi,Di,N), ...
                          [t_fault, t_end], X_sw1(end,:)', ode_opts);

t_sw  = [t_sw1; t_sw2(2:end)];
X_sw  = [X_sw1; X_sw2(2:end,:)];
th_sw = X_sw(:, 1:N);           % ángulos [rad]
om_sw = X_sw(:, N+1:2*N);       % desviación de freq. [rad/s]
fq_sw = (w0 + om_sw) / (2*pi);  % frecuencia absoluta [Hz]
fprintf('OK  (%d puntos)\n', length(t_sw));

%% ── 5. SIMULACIÓN: MODELO DE KURAMOTO ────────────────────────────────
%  Estado: θ = [θ₁…θ_N]'   (N variables)
%  D_i·dθ_i/dt = P_i - Σ_j K_ij·sin(θ_i - θ_j)
%  (límite sobreAmortiguado de la Swing Equation, Dörfler & Bullo 2012)

th0_ku = delta_op;

fprintf('Simulando Kuramoto...       ');
[t_ku1, TH_ku1] = ode15s(@(t,th) kuramoto(t,th,Pi,      K,Di,N), ...
                           [0, t_fault], th0_ku, ode_opts);
[t_ku2, TH_ku2] = ode15s(@(t,th) kuramoto(t,th,Pi_fault,K,Di,N), ...
                           [t_fault, t_end], TH_ku1(end,:)', ode_opts);

t_ku  = [t_ku1; t_ku2(2:end)];
TH_ku = [TH_ku1; TH_ku2(2:end,:)];
fprintf('OK  (%d puntos)\n', length(t_ku));

%% ── 6. PARÁMETRO DE ORDEN DE KURAMOTO r(t) ───────────────────────────
%  r(t) = |1/N · Σ_j exp(i·θ_j(t))|
%  r → 1: sincronización perfecta,  r → 0: pérdida total de sincronismo

r_sw = abs(mean(exp(1i * th_sw), 2));
r_ku = abs(mean(exp(1i * TH_ku), 2));

%% ── 7. FIGURA 1: Ángulos de rotor — Swing Equation ───────────────────
colors = lines(N);
figure('Name','Swing Equation – Ángulos de rotor','NumberTitle','off');
hold on;
for i = 1:N
    plot(t_sw, rad2deg(th_sw(:,i)), 'Color',colors(i,:), 'LineWidth',1.2);
end
xline(t_fault,'k--','LineWidth',1.5,'Label','Fallo B11');
legend(nodeNames,'Location','eastoutside','FontSize',7,'Interpreter','none');
xlabel('Tiempo (s)'); ylabel('Ángulo de rotor δ_i (°)');
title('Swing Equation – Ángulos de rotor');
grid on; box on;

%% ── 8. FIGURA 2: Desviación de frecuencia — Swing Equation ───────────
figure('Name','Swing Equation – Frecuencia de nodos','NumberTitle','off');
hold on;
for i = 1:N
    plot(t_sw, fq_sw(:,i), 'Color',colors(i,:), 'LineWidth',1.2);
end
xline(t_fault,'k--','LineWidth',1.5,'Label','Fallo B11');
yline(50,'k:','LineWidth',1,'Label','50 Hz');
legend(nodeNames,'Location','eastoutside','FontSize',7,'Interpreter','none');
xlabel('Tiempo (s)'); ylabel('Frecuencia (Hz)');
title('Swing Equation – Frecuencia de nodos');
grid on; box on;

%% ── 9. FIGURA 3: Comparación Swing vs Kuramoto ────────────────────────
figure('Name','Comparación Swing vs Kuramoto','NumberTitle','off');

subplot(2,2,1);
hold on;
for i = 1:N
    plot(t_sw, rad2deg(th_sw(:,i)), 'Color',colors(i,:), 'LineWidth',1);
end
xline(t_fault,'k--'); title('Swing – Ángulos (°)');
xlabel('t (s)'); ylabel('δ_i (°)'); grid on;

subplot(2,2,2);
hold on;
for i = 1:N
    plot(t_ku, rad2deg(TH_ku(:,i)), 'Color',colors(i,:), 'LineWidth',1);
end
xline(t_fault,'k--'); title('Kuramoto – Ángulos (°)');
xlabel('t (s)'); ylabel('θ_i (°)'); grid on;

subplot(2,2,3);
hold on;
plot(t_sw, r_sw, 'b-', 'LineWidth',2, 'DisplayName','Swing');
plot(t_ku, r_ku, 'r--','LineWidth',2, 'DisplayName','Kuramoto');
xline(t_fault,'k--','LineWidth',1.5);
yline(1,'k:'); yline(0,'k:');
legend; xlabel('t (s)'); ylabel('r(t)');
title('Parámetro de orden (sincronismo r→1)'); grid on;

subplot(2,2,4);
% Diferencia de ángulo entre B11 (afectado) y el slack B14
i_ref = 14; i_flt = 11;
hold on;
plot(t_sw, rad2deg(th_sw(:,i_flt) - th_sw(:,i_ref)), ...
     'b-','LineWidth',1.5,'DisplayName','Swing');
plot(t_ku, rad2deg(TH_ku(:,i_flt) - TH_ku(:,i_ref)), ...
     'r--','LineWidth',1.5,'DisplayName','Kuramoto');
xline(t_fault,'k--','LineWidth',1.5);
legend; xlabel('t (s)'); ylabel('δ_{11} - δ_{14} (°)');
title(sprintf('Ángulo relativo: %s – %s', ...
    strtrim(nodeNames{i_flt}), strtrim(nodeNames{i_ref})));
grid on;
sgtitle(sprintf('Pérdida de generación en %s (t = %.0fs)', ...
    strtrim(nodeNames{11}), t_fault), 'FontWeight','bold');

%% ── 10. FIGURA 4: Red coloreada por ángulo post-fallo ────────────────
t_snap  = t_end;              % instante para el snapshot
[~,idx] = min(abs(t_sw - t_snap));
delta_snap = th_sw(idx, :);   % ángulos en ese instante

figure('Name','Red – ángulos post-fallo','NumberTitle','off');
p2 = plot(G, 'XData',lon, 'YData',lat, 'NodeLabel',nodeNames, ...
           'LineWidth',1.2, 'MarkerSize',10);
colormap(gca, parula);
p2.NodeCData = rad2deg(delta_snap);
colorbar; clim([min(rad2deg(delta_snap))-1, max(rad2deg(delta_snap))+1]);
title(sprintf('Ángulos de rotor en t = %.0f s (post-fallo B11)', t_snap));
xlabel('Longitud'); ylabel('Latitud'); grid on; axis equal;

%% ── 11. ANÁLISIS ESPECTRAL DEL LAPLACIANO ─────────────────────────────
figure('Name','Espectro del Laplaciano ponderado','NumberTitle','off');

subplot(1,2,1);
bar(sort(real(eig(Lw))));
xlabel('Índice'); ylabel('Valor propio λ_k (GW)');
title('Espectro de L_w (Laplaciano ponderado)');
grid on;

subplot(1,2,2);
% Modo de sincronización: segundo vector propio (Fiedler vector)
[V_eig, D_eig] = eig(Lw);
[ev_sort, idx_sort] = sort(diag(D_eig));
fiedler_vec = V_eig(:, idx_sort(2));   % vector de Fiedler
scatter(lon, lat, 200, fiedler_vec, 'filled');
colormap(gca, coolwarm_cmap());
colorbar;
for i = 1:N
    text(lon(i)+0.1, lat(i)+0.1, sprintf('B%02d',i), 'FontSize',7);
end
title(sprintf('Vector de Fiedler  (λ₂ = %.4f GW)',ev_sort(2)));
xlabel('Longitud'); ylabel('Latitud');
grid on; axis equal;

%% ── 12. RESUMEN NUMÉRICO ─────────────────────────────────────────────
fprintf('\n────────────────────────────────────────\n');
fprintf('RESUMEN DE LA SIMULACIÓN\n');
fprintf('────────────────────────────────────────\n');

% Frecuencia mínima post-fallo (Swing)
[f_min, i_min] = min(fq_sw(:));
fprintf('Frecuencia mínima (Swing):  %.4f Hz  en nodo %s  (t = %.2f s)\n', ...
        f_min, nodeNames{mod(i_min-1,N)+1}, t_sw(ceil(i_min/N)));

% ¿Se recupera el sistema?
f_final_sw = mean(fq_sw(end,:));
fprintf('Frecuencia media final:     %.4f Hz\n', f_final_sw);
if abs(f_final_sw - 50) < 0.1
    fprintf('Estado final: SINCRONIZADO (|Δf| < 0.1 Hz)\n');
else
    fprintf('Estado final: DESINCRONIZADO (|Δf| = %.3f Hz)\n', abs(f_final_sw-50));
end

fprintf('Parámetro de orden final r: %.4f (Swing), %.4f (Kuramoto)\n', ...
        r_sw(end), r_ku(end));
fprintf('────────────────────────────────────────\n');

%% ═══════════════════════════════════════════════════════════════════════
%                        FUNCIONES LOCALES
% ═══════════════════════════════════════════════════════════════════════

function dxdt = swing_eq(~, x, Pi_vec, K, Mi, Di, N)
% SWING_EQ  ODE de la Swing Equation para la red completa.
%   Estado: x = [θ₁…θ_N, ω₁…ω_N]'
%   M_i·ω̇_i = P_i - D_i·ω_i - Σ_j K_ij·sin(θ_i - θ_j)
    theta = x(1:N);
    omega = x(N+1:2*N);

    % Potencia eléctrica: P_e,i = Σ_j K_ij·sin(θ_i - θ_j)
    % Forma vectorizada: P_e = sum_j K_ij * sin(theta_i - theta_j)
    dth = theta - theta.';            % matriz de diferencias de ángulo [N×N]
    Pe  = sum(K .* sin(dth), 2);      % suma por filas → [N×1]

    dxdt        = zeros(2*N, 1);
    dxdt(1:N)   = omega;                                      % dθ/dt = ω
    dxdt(N+1:2*N) = (Pi_vec - Di.*omega - Pe) ./ Mi;         % dω/dt
end

function dth_dt = kuramoto(~, theta, Pi_vec, K, Di, N)
% KURAMOTO  Modelo no-uniforme de Kuramoto para red eléctrica.
%   D_i·dθ_i/dt = P_i - Σ_j K_ij·sin(θ_i - θ_j)
%   (límite sobreAmortiguado de la Swing Equation)
    dth    = theta - theta.';         % diferencias angulares [N×N]
    Pe     = sum(K .* sin(dth), 2);   % flujos de potencia [N×1]
    dth_dt = (Pi_vec - Pe) ./ Di;     % [N×1]
end

function mismatch = pf_mismatch(delta_free, Pi_vec, K, slack, N)
% PF_MISMATCH  Residuo de las ecuaciones de flujo de potencia DC.
%   Para usar con fsolve: busca δ* tal que P_i = Σ_j K_ij·sin(δ_i - δ_j)
    non_slack = setdiff(1:N, slack);
    delta     = zeros(N,1);
    delta(non_slack) = delta_free;
    dth       = delta - delta.';
    P_calc    = sum(K .* sin(dth), 2);
    mismatch  = P_calc(non_slack) - Pi_vec(non_slack);
end

function cmap = coolwarm_cmap()
% Mapa de colores azul–blanco–rojo (divergente) para el vector de Fiedler
    n = 64;
    r = [linspace(0.23,1,n/2), linspace(1,0.71,n/2)]';
    g = [linspace(0.30,1,n/2), linspace(1,0.00,n/2)]';
    b = [linspace(0.75,1,n/2), linspace(1,0.15,n/2)]';
    cmap = [r, g, b];
end