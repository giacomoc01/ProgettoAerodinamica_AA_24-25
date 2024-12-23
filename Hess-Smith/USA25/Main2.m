%% traccia Hess Smith (2024)

clc
close all
clear 

%addpath mat_functions

%% Input

U_inf = 1;  % Velocità all'infinito [m/s]
alpha = 2;   % Angolo di incidenza [°]
U_inf_x = U_inf * cos(deg2rad(alpha));
U_inf_y = U_inf * sin(deg2rad(alpha));

U_inf = [U_inf_x; U_inf_y];
U_inf_normal = [-U_inf(2); U_inf(1)]; %versore normale 
U_inf_normal = U_inf_normal ./ norm(U_inf_normal); %vettore calcolato come versore/modulo

TestCase = 0;

LE_X_Position = 0;
LE_Y_Position = 0;

%% Creazione profilo

% numero profilo:
chord=100;
Npannelli=32;
%[x,y]=createProfile('0012',Npannelli,chord);

Corpo = importXfoilProfile('usa25.dat',2,34);

% Prima flippa i vettori
x = flipud(Corpo.x);
y = flipud(Corpo.y);

Corpo.x = x.*chord;
Corpo.y = y.*chord;

x = x - chord/2;

figure(1);
plot(x/chord, y/chord, 'o-', 'LineWidth', 1)
xlabel('x/c')
ylabel('y/c')
title('Airfoil Profile')
axis equal
grid on
legend('Profile points')

%% Creazione di una struttura di pannelli

[Centro, Normale, Tangente, Estremo_1, Estremo_2, alpha, lunghezza, L2G_TransfMatrix, G2L_TransfMatrix] = CreaStrutturaPannelli(Corpo);
        
%% Inizializzazione matrici e vettori

% Ora che ho i pannelli, posso inizializzare la matrice ed i vettori

NCols = sum(Npannelli) + 1; %numero di colonne, penso sia numero dei pannelli+1 perchè c'è la condizione di non penetrazione per ogni pannello e poi una condizione di kutta
NRows = NCols; %stessa cosa di prima
matriceA = zeros(NRows, NCols); %matrice dei coefficienti
TermineNoto = zeros(NRows, 1); %vettore da una colonna e righe quante quelle della matrice

%% Creazione della matrice quadrata As

%questa è la parte che sugli appunti il prof ha chiamato parti A e C della
%matrice A
for i = 1:Npannelli
    index_i = i; % riga

    Centro_qui = Centro(i, :)'; %righe vettori centri dei pannelli
    Normale_qui = Normale(i, :)'; %righe vettori normale dei pannelli

    indexStart_colonna = 0;

        for j = 1:Npannelli
            index_j = indexStart_colonna + j;  % Colonna

            Estremo_1_qui = Estremo_1(j, :)'; %righe estremi
            Estremo_2_qui = Estremo_2(j, :)';

            L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :)); %matrice rotazione local to global
            G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :)); %matrice rotazione global to local

            matriceA(index_i, index_j) = dot(ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Normale_qui);

            matriceA(index_i, sum(Npannelli)+1) = matriceA(index_i, sum(Npannelli)+1) + dot(ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Normale_qui);


        end

end


% Creazione delle componenti dei vettori a_v, c_s e c_v, che sono le parti
% della matrice A che il prof ha chiamato B e D


Centro_Start = Centro(1, :)';
Tangente_Start = Tangente(1, :)';

Centro_End = Centro(end, :)';
Tangente_End = Tangente(end, :)';


b = 0;
for j = 1:Npannelli(1)

    index_j = j;

    Estremo_1_qui = Estremo_1(j, :)';
    Estremo_2_qui = Estremo_2(j, :)';
    L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :)); 
    G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

    a = dot(ViSorgente(Centro_Start, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_Start);
    b = b + dot(ViVortice(Centro_Start, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_Start);

    a = a + dot(ViSorgente(Centro_End, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_End);
    b = b + dot(ViVortice(Centro_End, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_End);


    matriceA(sum(Npannelli) + 1, index_j) = a;

end

matriceA(sum(Npannelli) + 1, sum(Npannelli) + 1) = b;



%% Creazione del termine noto

for j = 1:Npannelli

    Normale_qui = Normale(j, :)';

    index = j;

    TermineNoto(index) = - dot(U_inf, Normale_qui);
end

Tangente_1 = Tangente(1, :)';
Tangente_end = Tangente(end, :)';
TermineNoto(sum(Npannelli) + 1) = - dot(U_inf, (Tangente_1 + Tangente_end)); %formule sugli appunti

%% Risoluzione sistema lineare
Soluzione = linsolve(matriceA,TermineNoto);
q=Soluzione(1:Npannelli); %sorgenti sugli N pannelli
gamma=Soluzione(Npannelli+1); %intensità del vortice
GAMMA=[];
for i=1:Npannelli
  lunghezza(i) = norm(Estremo_2(i, :) - Estremo_1(i, :));
end


%% Calcolo dell'angolo di Theodorsen
M = [x,y];

Zup = M(18:end,:);
Zdown = M(1:17,:);

Zup = [[-50 , 0]; Zup];

a = -50;
b = 50;
n = 1000;
h = (b-a)/(n-1);




%%

k = linspace(a, b, n); % Creiamo una griglia comune lungo x

z_u_interp = @ (k) spline(Zup(:,1), Zup(:,2),k); % Interpolazione lineare per z_u
z_l_interp = @ (k) spline(Zdown(:,1), Zdown(:,2),k); % Interpolazione lineare per z_l

% Calcolo della linea media
z_m = @ (k) (z_u_interp(k) + z_l_interp(k))./2;

Dz_m = @ (k) gradient(z_m(k),k);

f = @ (x) Dz_m(x)./(sqrt((chord.^2)/4 - x.^2));

AlphaThRad = (integral(f,-50,50))./pi;

AlphaThDeg = (AlphaThRad*180)/pi;

fprintf('L''angolo di Theodorsen vale:%f\n',AlphaThDeg)

figure(5);
h1 = plot(x/chord, y/chord, 'o-', 'LineWidth', 1, 'DisplayName', 'Airfoil Profile');
hold on
h2 = plot(x/chord, z_m(x)/chord, '-', 'LineWidth', 1, 'DisplayName', 'Mean Camber Line');
h3 = plot(x/chord, Dz_m(x)/chord, '-', 'LineWidth', 1, 'DisplayName', 'Slope Distribution');

% Theodorsen angle visualization
x_arrow = 0;
y_arrow = 0;
h4 = quiver(x_arrow, y_arrow, cos(AlphaThRad), sin(AlphaThRad), 0.2, 'k', 'LineWidth', 1, 'DisplayName', 'Theodorsen Angle');

% Angle arc
th = linspace(0, AlphaThRad, 50);
r = 0.1;
xarc = r*cos(th);
yarc = r*sin(th);
h5 = plot(xarc, yarc, 'r--', 'LineWidth', 1, 'DisplayName', 'Angle Arc');

% Text label
text(r/2, r/2, ['\alpha_{Th} = ' num2str(AlphaThDeg,'%.2f') '°'], ...
    'HorizontalAlignment', 'left', 'Color', 'k')

xlabel('x/c')
ylabel('y/c')
title('Theodorsen Angle Analysis')
axis equal
grid on
legend([h1 h2 h3 h4 h5], 'Location', 'best')





% Parametri per l'analisi
h_values = linspace(1e-4, 1e-1, 50); % Valori di h (da 10^-4 a 10^-1)
AlphaThDeg_results = zeros(size(h_values)); % Array per salvare i risultati

for i = 1:length(h_values)
    h = h_values(i); % Aggiorna il regolarizzatore
    % Funzione con regolarizzatore
    f = @(x) Dz_m(x) ./ (sqrt((chord.^2)/4 - x.^2 + h));
    % Calcolo dell'angolo di Theodorsen
    AlphaThRad = integral(f, -50, 50) / pi;
    AlphaThDeg_results(i) = AlphaThRad * 180 / pi;
end

% Plot dei risultati (grafico 2D)
figure;
plot(h_values, AlphaThDeg_results, 'ro-', 'LineWidth', 0.5, 'MarkerSize', 5);
xlabel('Regolarizzatore h');
ylabel('Angolo di Theodorsen (\alpha_{Th}) [°]');
title('Sensibilità dell’angolo di Theodorsen rispetto al regolarizzatore h');
grid on;

% Parametri per l'analisi
Npannelli_values = 10:5:100; % Valori del numero di pannelli
AlphaThDeg_results = zeros(size(Npannelli_values)); % Array per i risultati

for i = 1:length(Npannelli_values)
    Npannelli = Npannelli_values(i); % Aggiorna il numero di pannelli

    % Creazione del profilo con il nuovo numero di pannelli
    [x, y] = createProfile('0012', Npannelli, chord); % Funzione per creare il profilo
    x = x - chord / 2; % Centra il profilo sulla corda

    % Separazione parte superiore e inferiore del profilo
    x_u = x(ceil(Npannelli/2):end);
    x_l = x(1:ceil(Npannelli/2));
    z_u = y(ceil(Npannelli/2):end); % Parte superiore
    z_l = y(1:ceil(Npannelli/2)); % Parte inferiore

    % Rimuovere i duplicati (ordinare e mantenere solo i valori unici)
    [x_u, idx_u] = unique(x_u, 'stable'); % 'stable' mantiene l'ordine originale
    z_u = z_u(idx_u); % Riorganizza anche i valori di z_u
    [x_l, idx_l] = unique(x_l, 'stable');
    z_l = z_l(idx_l);

    % Creazione di una griglia comune
    k_grid = linspace(min(x), max(x), 1000); % Griglia uniforme lungo x

    % Interpolazione spline per linea media
    z_u_interp = interp1(x_u, z_u, k_grid, 'spline'); % Parte superiore
    z_l_interp = interp1(x_l, z_l, k_grid, 'spline'); % Parte inferiore
    z_m = (z_u_interp + z_l_interp) / 2; % Linea media

    % Calcolo derivata della linea media
    Dz_m = gradient(z_m, k_grid); % Derivata lungo la griglia uniforme

    % Funzione integranda con regolarizzatore
    f = @(x) interp1(k_grid, Dz_m, x, 'linear', 'extrap') ./ ...
             sqrt((chord.^2) / 4 - x.^2 + 0.01);

    % Calcolo dell'angolo di Theodorsen
    AlphaThRad = integral(f, -chord/2, chord/2) / pi;
    AlphaThDeg_results(i) = AlphaThRad * 180 / pi;
end

% Plot dei risultati
figure;
plot(Npannelli_values, AlphaThDeg_results, '-o', 'LineWidth', 0.5, 'MarkerSize', 5);
ylim([-5 4])
xlabel('Numero di pannelli (N_{pannelli})');
ylabel('Angolo di Theodorsen (\alpha_{Th}) [°]');
title('Sensibilità dell’angolo di Theodorsen rispetto al numero di pannelli');
grid on;



