clear; clc; close all;

%% ===== FUNCIONES AUXILIARES =====
% Función para calcular distancia real entre nodos
function d = haversine_km(lat1, lon1, lat2, lon2)
    R = 6371.0; % radio medio terrestre [km]
    p1 = deg2rad(lat1);
    p2 = deg2rad(lat2);
    dp = deg2rad(lat2-lat1);
    dl = deg2rad(lon2-lon1);

    a = sin(dp/2)^2 + cos(p1)*cos(p2)*sin(dl/2)^2;
    d = 2*R*asin(sqrt(a));
end

% Funcion swing equation
% Le ofreces el estado x=[theta,thetapunto] y te devuelve dx/dt
function dxdt = swing_eq(~, x, Pi_vec, Kmat, Mi_vec, Di_vec, N)
    theta = x(1:N);
    omega = x(N+1:2*N); %thetapunto

    dth = theta - theta.'; %matriz diferencias de desfase
    Pe  = sum(Kmat .* sin(dth), 2);

    dxdt = zeros(2*N,1);
    dxdt(1:N) = omega;
    dxdt(N+1:2*N) = (Pi_vec - Di_vec.*omega - Pe) ./ Mi_vec;
end

% Funcion para la busqueda del estado estacionario 
% Esto implica thetapunto=0 llegando a Pi=Kij*sen(...)
function mismatch = pf_mismatch(delta_free, Pi_vec, Kmat, slack, N)
    non_slack = setdiff(1:N, slack);
    delta = zeros(N,1);
    delta(non_slack) = delta_free; %vector desfases construido 

    dth = delta - delta.';
    Pcalc = sum(Kmat .* sin(dth), 2);

    mismatch = Pcalc(non_slack) - Pi_vec(non_slack); %debemos de buscar que el mismatch sea 0
end

%% ===== DEFINICION DEL GRAFO =====
nodeNames = {'B01 Galicia-Norte','B02 Galicia-Sur','B03 Asturias-Centro','B04 Cantabria-Mataporquera', ...
            'B05 Pais Vasco-Gatika','B06 Navarra-Muruarte','B07 Aragon-Ebro','B08 Aragon-Pirineos', ...
            'B09 CastillaLeon-Oeste','B10 CastillaLeon-Este','B11 Extremadura-Almaraz','B12 Extremadura-Badajoz', ...
            'B13 CLM-Oeste','B14 CLM-Este-Slack','B15 Madrid-400kV','B16 Andalucia-Oeste', ...
            'B17 Andalucia-Este','B18 Cataluna-Tarragona','B19 Cataluna-Norte','B20 CValenciana-Saguntum'}';

% Aristas propuestas (nodos iniciales y destino)
sEdge = [1 1 2 3 4 5 6 7 8 5 6 9 9 9 10 11 11 13 14 12 13 13 16 7 18 18 15 14 17]; 
tEdge = [2 3 9 4 5 6 7 8 19 10 10 10 11 13 14 12 13 14 15 16 16 17 17 18 19 20 20 20 20];

G = graph(sEdge, tEdge, [], nodeNames);

% Coordenadas geográficas
lon = [-7.5 -8.4 -5.8 -4.0 -2.9 -1.6 -0.9 0.6 -5.7 -4.7 -5.7 -6.8 -4.0 -3.0 -3.7 -6.0 -3.2 1.2 2.2 -0.2]';
lat = [43.2 42.3 43.4 43.0 43.3 42.7 41.6 42.7 41.8 41.8 39.8 38.9 38.9 39.8 40.4 37.5 37.2 41.1 42.1 39.7]';

figure;
plot(G, 'XData', lon, 'YData', lat, 'NodeLabel', nodeNames, 'LineWidth', 1.2, 'MarkerSize', 5);
title('Backbone red peninsular con coordenadas geograficas aproximadas');
xlabel('Longitud'); ylabel('Latitud');
axis equal; grid on;

A = full(adjacency(G));
L = full(laplacian(G));
B = full(incidence(G));

%% ===== PARAMETROS SWING EQUATION =====
N  = numnodes(G);
f0 = 50;
w0 = 2*pi*f0;

% Pi [GW] - Potencia neta media que sale o entra de un nodo
Pi = [0.623 0.265 -0.093 -0.367 -1.096 0.166 1.007 0.140 0.752 0.589 2.186 0.350 0.179 1.622 -2.999 0.098 -1.038 1.811 -2.893 -1.302]';
fprintf('Suma Pi = %.4e GW\n', sum(Pi));

% Ssyn [GW] - Medida de la masa rotante útil (equipo)
Ssyn = [1.10 0.55 0.70 0.30 1.00 0.80 1.10 0.45 1.10 0.80 2.05 0.25 0.55 1.55 0.35 1.80 1.00 3.75 1.20 2.10]';

% Hbar [s] - Cuánta energía cinética por unidad de potencia base tiene ese equipo
Hbar = [3.6 3.6 4.0 3.0 4.7 4.2 3.8 3.4 3.5 3.5 5.8 3.2 4.0 5.0 3.5 4.2 4.2 5.6 4.2 4.0]';

% Juntando ambas para encontrar la Inercia (Mstar y Mi)
Mstar = 2.*Hbar.*Ssyn;   % [GW*s]
% La inercia que nos interesa en las ecuaciones va dividida de w0
Mi = Mstar / w0;         % [GW*s^2/rad]

% Amortiguamiento Di
Di = [0.080 0.060 0.050 0.035 0.075 0.065 0.085 0.055 0.090 0.085 0.130 0.065 0.075 0.110 0.110 0.100 0.085 0.120 0.125 0.100]';

%% ===== MATRIZ DE ACOPLAMIENTO K =====
% MATRIZ DE ACOPLAMIENTO K (Factor de cesion de potencia)
% K_ij ~ n_ij * V^2 / (xkm * L_ij)
% V = 400 kV | xkm = 0.40 ohm/km | L_ij en km | K_ij en GW

V_kV = 400;   % tension backbone
xkm  = 0.40;  % reactancia serie tipica [ohm/km]

edgeNames = {'B01-B02','B01-B03','B02-B09','B03-B04','B04-B05','B05-B06','B06-B07','B07-B08', ...
            'B08-B19','B05-B10','B06-B10','B09-B10','B09-B11','B09-B13','B10-B14','B11-B12', ...
            'B11-B13','B13-B14','B14-B15','B12-B16','B13-B16','B13-B17','B16-B17','B07-B18', ...
            'B18-B19','B18-B20','B15-B20','B14-B20','B17-B20'}';

Nedges = length(sEdge);

% 1) Distancia geografica recta entre nodos [km]
d_geo = zeros(Nedges,1);
for e = 1:Nedges
    d_geo(e) = haversine_km(lat(sEdge(e)), lon(sEdge(e)), lat(tEdge(e)), lon(tEdge(e)));
end

% 2) Factor de correccion del trazado real alpha_ij
%    1.10 = directo | 1.15 = normal | 1.20 = sinuoso | 1.30 = montañoso
alpha = [1.15 1.20 1.20 1.15 1.15 1.10 1.10 1.20 1.30 1.20 1.15 1.10 1.15 1.15 1.10 1.10 1.10 1.10 1.10 1.15 1.15 1.15 1.15 1.15 1.10 1.10 1.10 1.10 1.15]';

% Longitud electrica aproximada del corredor [km]
Lij = alpha .* d_geo;

% 3) Numero de circuitos equivalentes n_ij (2 troncales, 1 normales)
nij = [2 1 1 1 1 1 1 1 1 1 1 2 1 1 2 1 2 2 2 1 1 1 1 1 2 2 2 2 1]';

% 4) Kedge en GW: K = n * V^2 / (xkm * L) -> dividido entre 1000 para GW
Kedge = (nij .* V_kV^2) ./ (xkm .* Lij) / 1000;   % [GW]

% 5) Construir matriz K simetrica
K = zeros(N,N);
for e = 1:Nedges
    K(sEdge(e), tEdge(e)) = Kedge(e);
    K(tEdge(e), sEdge(e)) = Kedge(e);
end

%% ===== EQUILIBRIO ESTACIONARIO =====
slack = 14; % índice del nodo de referencia desfase 0
non_slack = setdiff(1:N, slack);

% Punto inicial: ángulos pequeños proporcionales a la potencia
delta_guess = zeros(N,1); % vector de angulos iniciales
delta_guess(non_slack) = 0.02 * Pi(non_slack); % escalar para aproximar por angulos pequeños

opts_fs = optimoptions('fsolve', 'Display','iter', 'TolFun',1e-10, 'TolX',1e-10, 'MaxIter',2000, 'Algorithm','levenberg-marquardt');

% resolucion del sistema no lineal buscando pf_mismatch=0
[delta_free, ~, exitflag] = fsolve(@(x) pf_mismatch(x, Pi, K, slack, N), delta_guess(non_slack), opts_fs); 

if exitflag <= 0
    warning('No se ha encontrado equilibrio claro.');
end

delta_op = zeros(N,1);
delta_op(non_slack) = delta_free;

resid = sum(K .* sin(delta_op - delta_op.'), 2) - Pi;
fprintf('\nResiduo maximo equilibrio = %.3e GW\n', max(abs(resid)));

%vector definitivo de desfases iniciales en estado estacionario: delta_op
%% ===== CONTINGENCIA =====
t_fault = 2.0;
t_end   = 40.0;

% Escenario A: perdida de generacion en B11
Pi_fault = Pi;
Pi_fault(11) = 0;
K_fault = K;

%% ===== SIMULACION SWING =====
x0 = [delta_op; zeros(N,1)]; %estado inicial del estado estacionario
ode_opts = odeset('RelTol',1e-7, 'AbsTol',1e-9, 'MaxStep',0.02);

[t1, X1] = ode15s(@(tt,x) swing_eq(tt, x, Pi, K, Mi, Di, N), [0 t_fault], x0, ode_opts);
[t2, X2] = ode15s(@(tt,x) swing_eq(tt, x, Pi_fault, K_fault, Mi, Di, N), [t_fault t_end], X1(end,:)', ode_opts);

tSim = [t1; t2(2:end)];
XSim = [X1; X2(2:end,:)];

theta = XSim(:,1:N);
omega = XSim(:,N+1:2*N);
freq  = 50 + omega/(2*pi);

%% ===== FLUJOS DE LINEA =====
Ne = length(sEdge);
Fline = zeros(length(tSim), Ne);

for k = 1:length(tSim)
    for e = 1:Ne
        Fline(k,e) = K(sEdge(e),tEdge(e)) * sin(theta(k,sEdge(e)) - theta(k,tEdge(e)));
    end
end

%% ===== GRAFICAS =====
figure;
plot(tSim, rad2deg(theta), 'LineWidth',1.1);
xline(t_fault,'k--','LineWidth',1.4);
grid on; xlabel('t (s)'); ylabel('\theta_i (grados)'); title('Angulos nodales');

figure;
plot(tSim, freq, 'LineWidth',1.1);
xline(t_fault,'k--','LineWidth',1.4); yline(50,'k:');
grid on; xlabel('t (s)'); ylabel('Frecuencia (Hz)'); title('Frecuencia nodal');

figure;
plot(tSim, Fline, 'LineWidth',1.0);
xline(t_fault,'k--','LineWidth',1.4);
grid on; xlabel('t (s)'); ylabel('Flujo por linea (GW)'); title('Flujos de linea');

%% ===== VOLTAJES AC RECONSTRUIDOS MEJORADOS =====
nodesPlot = [11 14 15];   % B11, B14, B15
Vamp = 1.0;               % amplitud en p.u.
tFine = linspace(1.95, 2.10, 5000); % Malla fina para ver bien 50 Hz

thetaFine = zeros(length(tFine), N);
for i = 1:N
    thetaFine(:,i) = interp1(tSim, theta(:,i), tFine, 'pchip');
end

vAbs = zeros(length(tFine), length(nodesPlot));
vRel = zeros(length(tFine), length(nodesPlot));
dv   = zeros(length(tFine), length(nodesPlot));
vSlack = Vamp * cos(w0*tFine(:) + thetaFine(:,slack));

for m = 1:length(nodesPlot)
    ii = nodesPlot(m);
    vAbs(:,m) = Vamp * cos(w0*tFine(:) + thetaFine(:,ii));
    vRel(:,m) = Vamp * cos(w0*tFine(:) + (thetaFine(:,ii) - thetaFine(:,slack)));
    dv(:,m)   = vAbs(:,m) - vSlack;
end

figure('Name','Voltajes reconstruidos mejorados','NumberTitle','off');
subplot(3,1,1); plot(tFine, vAbs, 'LineWidth', 1.1); xline(t_fault,'k--','LineWidth',1.4);
grid on; xlabel('t (s)'); ylabel('v_i(t) [p.u.]'); title('Voltajes AC absolutos reconstruidos'); legend(nodeNames(nodesPlot), 'Interpreter','none', 'Location','best');
subplot(3,1,2); plot(tFine, vRel, 'LineWidth', 1.1); xline(t_fault,'k--','LineWidth',1.4);
grid on; xlabel('t (s)'); ylabel('v_i^{rel}(t) [p.u.]'); title('Voltajes reconstruidos referidos al slack');
subplot(3,1,3); plot(tFine, dv, 'LineWidth', 1.1); xline(t_fault,'k--','LineWidth',1.4);
grid on; xlabel('t (s)'); ylabel('\Deltav_i'); title('Diferencia de voltaje respecto al slack');

%% ===== ANGULOS RELATIVOS AL SLACK =====
figure('Name','Angulos relativos al slack','NumberTitle','off');
plot(tSim, rad2deg(theta - theta(:,slack)), 'LineWidth', 1.2);
xline(t_fault,'k--','LineWidth',1.4);
grid on; xlabel('t (s)'); ylabel('\theta_i - \theta_{slack} [grados]');
title('Angulos relativos al slack'); legend(nodeNames, 'Interpreter','none', 'Location','eastoutside');

%% ===== TABLA DE ARISTAS =====
Tedges = table(edgeNames, sEdge(:), tEdge(:), d_geo, alpha, Lij, nij, Kedge, ...
    'VariableNames', {'Edge','From','To','d_geo_km','alpha','Lij_km','nij','Kij_GW'});
disp(Tedges);
