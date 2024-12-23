%% traccia Hess Smith (2024)

clc
close all
clear 

addpath mat_functions

%% Input
chord=1;
Npannelli=101;
U_inf = 1;  % Velocità all'infinito [m/s]
alpha = 2;   % Angolo di incidenza [°]
U_inf_x = U_inf * cos(deg2rad(alpha));
U_inf_y = U_inf * sin(deg2rad(alpha));

U_inf = [U_inf_x; U_inf_y];
U_inf_normal = [-U_inf(2); U_inf(1)]; %versore normale 
U_inf_normal = U_inf_normal./norm(U_inf_normal); %vettore calcolato come versore/modulo

TestCase = 0;

LE_X_Position = 0;
LE_Y_Position = 0;

%% Creazione profilo

% numero profilo:

[x,y]=createProfile('0012',101,1);

Corpo = importXfoilProfile('NACA_0012.dat',2,103);
% Prima flippa i vettori
x = flipud(Corpo.x);
y = flipud(Corpo.y);
Corpo.x = x.*chord;
Corpo.y = y.*chord;

figure;
plot(x, y, 'o-')
title('Visualizzazione del profilo')
xlabel('x/c')
ylabel('y')
axis equal
grid on;
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
GAMMA=sum(lunghezza,'all')*gamma;
C_l=2*GAMMA*(U_inf(1)*chord);
fprintf('Il coefficiente di portanza per un profilo NACA0012 è: %16f\n',C_l)



%% Calcolo del cp e della velocità sui pannelli
%calcoliamo il cp tramite Bernoulli, con la formula cp,i=1-(componente
%tangenziale della velocità^2/velocità indisturbata^2)

Vt=zeros(Npannelli,1);

for i=1:Npannelli
    Estremo_1_qui = Estremo_1(i, :)';
    Estremo_2_qui = Estremo_2(i, :)';
    L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(i, :, :)); 
    G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(i, :, :));

    Velocita_indotta = zeros(2, 1);
    for j = 1:Npannelli
        Estremo_1_j = Estremo_1(j, :)';
        Estremo_2_j = Estremo_2(j, :)';
        L2G_TransfMatrix_j = squeeze(L2G_TransfMatrix(j, :, :));
        G2L_TransfMatrix_j = squeeze(G2L_TransfMatrix(j, :, :));
        
        % Calcolo del contributo delle sorgenti e vortici
        Velocita_indotta = Velocita_indotta + ...
            Soluzione(j) * ViSorgente(Centro(i, :)', Estremo_1_j, Estremo_2_j, L2G_TransfMatrix_j, G2L_TransfMatrix_j) + ...
            Soluzione(end) * ViVortice(Centro(i, :)', Estremo_1_j, Estremo_2_j, L2G_TransfMatrix_j, G2L_TransfMatrix_j);
    end
    
    % Velocità totale tangenziale
    Velocita_totale = U_inf + Velocita_indotta;
    V_t(i) = dot(Velocita_totale, Tangente(i, :)');
end

% Coefficiente di pressione sui pannelli
Cp = 1 - (V_t / norm(U_inf)).^2;



%% Calcolo dei coefficienti di momento

figure;
plot(Centro(:, 1), Cp, 'o-');
set(gca, 'YDir', 'reverse'); % Cp è negativo in alto
xlabel('x/c');
ylabel('C_p');
title('Distribuzione del coefficiente di pressione');
grid on;
% Lunghezza del pannello e posizione del centro rispetto a LE
%% Calcolo del coefficiente di momento rispetto al Leading Edge
Moment_LE = 0;
r_c=zeros(Npannelli,3);
n_c=zeros(Npannelli,3);
for i = 1:Npannelli
    % Calcolo del momento elementare
    r_c(i,1) = norm(Centro(i, 1)-Estremo_1(i,1)); 
    r_c(i,2)= norm(Centro(i, 2)-Estremo_1(i,2));
    n_c(i,1) = Normale(i, 1)';                                % Vettore normale al pannello
    n_c(i,2)= Normale(i,2);                             
    cross_product(i,:) = cross(r_c(i,:),n_c(i,:));   % Componente lungo z del prodotto vettoriale
    Moment_LE = Moment_LE + Cp(i) * (lunghezza(i) / chord^2) * cross_product(i,3);
end

% Coefficiente di momento normalizzato
C_m_LE = Moment_LE;

fprintf('C_m_LE (coeff. momento rispetto al LE): %f\n', C_m_LE);


%% Calcolo dell'angolo di Theodorsen
M = [x,y];

Zup = M(M(:,2) >=0 , :);
Zdown = M(M(:,2) <0 , :);

Zup = Zup(2:end-1,:);

k = linspace(0, 1, 50); % Creiamo una griglia comune lungo x
z_u_interp = @(k) spline(Zup(:,1), Zup(:,2), k); % Interpolazione lineare per z_u
z_l_interp = @ (k) spline(Zdown(:,1), Zdown(:,2), k); % Interpolazione lineare per z_l

% Calcolo della linea media
z_m = @ (k) (z_u_interp(k) + z_l_interp(k)) / 2;

Dz_m = @ (k) gradient(z_m(k),k);

f = @ (k) (Dz_m(k)./(sqrt((chord.^2)/4 - k.^2)))./pi;

AlphaTh = integral(f,0,1);

fprintf('Angolo di Theodorsen per NACA 0012: %f\n',AlphaTh)



