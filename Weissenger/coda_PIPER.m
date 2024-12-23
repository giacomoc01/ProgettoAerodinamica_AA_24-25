%% SECONDA PARTE PIPER
% CON coda
close all
clear all
clc

%% Configuration
U_Inf_Mag = 1;
beta = 0;
alpha = 5;
U_Inf = U_Inf_Mag * [cosd(alpha)*cosd(beta), sind(beta), sind(alpha)*cosd(beta)];
rho = 1.225;

% Wing and tail configuration
config.NCorpi = 2;
% Wing parameters (index 1)
config.RootChord = [1.5962, 0.55];  % [wing, tail]
config.DihedralAngle = [0, 0];
config.SweepAngle = [0.6, 0.9];
config.TaperRatio = [1, 0.83];
config.AspectRatio = [6.944, 6.9];
config.Span = [10.73, 3.2];
config.LEPosition_X = [0, 6];
config.LEPosition_Y = [0, 0];
config.LEPosition_Z = [0, 1];
config.RotationAngle_X = [0, 0];
config.RotationAngle_Y = [0, 0];
config.RotationAngle_Z = [0, 0];

% Discretization
config.SemiSpanwiseDiscr = [20, 10];  % Less panels for tail
config.ChordwiseDiscr = [20, 10];

%% Preliminary computations
config.SemiSpan = config.Span./2;
config.Surface = 2 * (config.SemiSpan .* config.RootChord .* (1+config.TaperRatio)./2);
config.SurfaceProjected = config.Surface .* cosd(config.DihedralAngle);
config.TipChord = config.RootChord .* config.TaperRatio;
config.MAC = (2/3) .* config.RootChord .* ((1 + config.TaperRatio + config.TaperRatio.^2)./(1 + config.TaperRatio));

%% Create geometry structures for both wing and tail
ControlPoints = cell(config.NCorpi, 1);
InducedPoints = cell(config.NCorpi, 1);
Normals = cell(config.NCorpi, 1);
InfiniteVortices = cell(config.NCorpi, 1);
Vortices = cell(config.NCorpi, 1);
internalMesh = cell(config.NCorpi, 1);
WingExtremes = cell(config.NCorpi, 1);

for iCorpo = 1:config.NCorpi
    [ControlPoints{iCorpo}, InducedPoints{iCorpo}, Normals{iCorpo}, ...
        InfiniteVortices{iCorpo}, Vortices{iCorpo}, internalMesh{iCorpo}, ...
        WingExtremes{iCorpo}] = createStructure(config, iCorpo);
end

%% Matrix initialization
NPanelsTot = 0;
for iCorpo = 1:config.NCorpi
    NPanelsTot = NPanelsTot + 2*config.SemiSpanwiseDiscr(iCorpo)*config.ChordwiseDiscr(iCorpo);
end

matriceA = zeros(NPanelsTot, NPanelsTot);
TermineNoto = zeros(NPanelsTot, 1);

%% Construction of the matrix
rowIndex = 0;
for iCorpo = 1:config.NCorpi

    % Cycle on all of its chordwise panels
    for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
        % Cycle on all of its spanwise panels
        for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)

            % Update row index
            rowIndex = rowIndex + 1;

            columnIndex = 0;

            ControlPointHere = ControlPoints{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
            NormalHere = Normals{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;


            for jCorpo = 1:config.NCorpi

                % Cycle on all of its chordwise panels
                for ChordPanel_j = 1:config.ChordwiseDiscr(jCorpo)
                    % Cycle on all of its spanwise panels
                    for SpanPanel_j = 1:2*config.SemiSpanwiseDiscr(jCorpo)

                        % Update column index
                        columnIndex = columnIndex + 1;

                        % Compute the influence induced by first
                        % semi-infinite vortex
                        Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.toInfty;
                        Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.onWing;
                        U = vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);

                        % Compute the influence induced by finite vortex
                        Extreme_1 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root;
                        Extreme_2 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip;
                        U = U + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);

                        % Compute the influence induced by second
                        % semi-infinite vortex
                        Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.onWing;
                        Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.toInfty;
                        U = U + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);

                        matriceA(rowIndex, columnIndex) = dot(U, NormalHere);


                    end
                end
            end



        end
    end
end


%% Construct right-hand side
rowIndex = 0;
for iCorpo = 1:config.NCorpi
    for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
        for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
            rowIndex = rowIndex + 1;
            NormalHere = Normals{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
            TermineNoto(rowIndex) = -dot(U_Inf, NormalHere);
        end
    end
end

%% Solve the linear system
Solution = linsolve(matriceA, TermineNoto);
Gamma = cell(config.NCorpi, 1);
rowIndex = 0;
for iCorpo = 1:config.NCorpi

    Gamma{iCorpo} = zeros(config.ChordwiseDiscr(iCorpo), config.SemiSpanwiseDiscr(iCorpo)*2 );

    % Cycle on all of its chordwise panels
    for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
        % Cycle on all of its spanwise panels
        for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)

            % Update row index
            rowIndex = rowIndex + 1;

            Gamma{iCorpo}(ChordPanel_i, SpanPanel_i) = Solution(rowIndex);
        end

    end

end


%% Compute the 3D Lift
Lift3D = [0, 0];
for iCorpo = 1:config.NCorpi
    dSpan = config.SemiSpan(iCorpo) / config.SemiSpanwiseDiscr(iCorpo);
    for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
        for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
            GammaHere = Gamma{iCorpo}(ChordPanel_i, SpanPanel_i);
            Lift3D(iCorpo) = Lift3D(iCorpo) + rho * U_Inf_Mag * GammaHere * dSpan;
        end
    end
end
disp(['3D Lift: ' num2str(Lift3D) ' N']);


%% Compute 2D Lift distribution
Lift2D = cell(config.NCorpi, 1);
for iCorpo = 1:config.NCorpi
    Lift2D{iCorpo} = rho * U_Inf_Mag * sum(Gamma{iCorpo}, 1) .* config.MAC(iCorpo);
end
% Plot 2D Lift distribution
figure;
for iCorpo = 1:config.NCorpi
    y = linspace(-config.SemiSpan(iCorpo), config.SemiSpan(iCorpo), 2 * config.SemiSpanwiseDiscr(iCorpo));
    plot(y, Lift2D{iCorpo}, 'LineWidth', 1.5);
    hold on;
end
title('Distribuzione della Portanza (2D)');
xlabel('Spanwise Position [m]');
ylabel('Lift per unit span [N/m]');
grid on;

%% Calcolo dei coefficienti aerodinamici
% Dynamic pressure
q_inf = 0.5 * rho * U_Inf_Mag^2;

% Lift coefficient
for iCorpo = 1:config.NCorpi
    CL_3D(iCorpo) = Lift3D(iCorpo) / (q_inf * config.Surface(iCorpo));
end

%% Calcolo della resistenza indotta
% Fattore di Oswald
e = 0.85;

% Calcolo del coefficiente di resistenza indotta
CDi = [0,0];
for iCorpo = 1:config.NCorpi
    CDi(iCorpo) = CL_3D(iCorpo)^2 / (pi * config.AspectRatio(iCorpo) * e);
end

% Calcolo della resistenza indotta
for iCorpo = 1:config.NCorpi
    Di = CDi(iCorpo) * q_inf * config.Surface(iCorpo);
end

% Stampa dei risultati
fprintf('Coefficiente di portanza: CL = %.4f\n', CL_3D);
fprintf('Coefficiente di resistenza indotta: CDi = %.4f\n', CDi);
fprintf('Resistenza indotta: Di = %.4f N\n', Di);

% Calcolo dell'efficienza aerodinamica
E = CL_3D/CDi;
fprintf('Efficienza aerodinamica: E = %.4f\n', E);

%% Linea CL_alfa
alpha_range = linspace(-10, 10, 21);
CL_values = zeros(config.NCorpi, length(alpha_range));

for i = 1:length(alpha_range)
    alpha = alpha_range(i);
    U_Inf = U_Inf_Mag * [cosd(alpha)*cosd(beta), sind(beta), sind(alpha)*cosd(beta)];

    rowIndex = 0;
    for iCorpo = 1:config.NCorpi
        for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
            for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
                rowIndex = rowIndex + 1;
                NormalHere = Normals{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
                TermineNoto(rowIndex) = -dot(U_Inf, NormalHere);
            end
        end
    end

    Solution = linsolve(matriceA, TermineNoto);

    for iCorpo = 1:config.NCorpi
        Lift = 0;
        dSpan = config.SemiSpan(iCorpo) / config.SemiSpanwiseDiscr(iCorpo);
        GammaHere = reshape(Solution((iCorpo-1)*config.ChordwiseDiscr(iCorpo)*2*config.SemiSpanwiseDiscr(iCorpo)+1:iCorpo*config.ChordwiseDiscr(iCorpo)*2*config.SemiSpanwiseDiscr(iCorpo)), config.ChordwiseDiscr(iCorpo), []);
        Lift = Lift + sum(sum(GammaHere)) * rho * U_Inf_Mag * dSpan;
        CL_values(iCorpo, i) = Lift / (0.5 * rho * U_Inf_Mag^2 * config.Surface(iCorpo));
    end
end

% Plot CL vs alpha for each body
figure;
hold on;
for iCorpo = 1:config.NCorpi
    plot(alpha_range, CL_values(iCorpo, :), 'o-', 'DisplayName', ['Body ' num2str(iCorpo)]);
end
xlabel('Angle of Attack (degrees)');
ylabel('Lift Coefficient (C_L)');
title('Lift Coefficient vs Angle of Attack');
legend('Profilo Alare','Piani di coda');
grid on;
hold off;



%% Plot CL vs CD for wing and tail
% Define alpha range
alpha_range = linspace(-10, 10, 21);
CL_values = zeros(config.NCorpi, length(alpha_range));
CD_values = zeros(config.NCorpi, length(alpha_range));

for i = 1:length(alpha_range)
    alpha = alpha_range(i);
    U_Inf = U_Inf_Mag * [cosd(alpha)*cosd(beta), sind(beta), sind(alpha)*cosd(beta)];

    rowIndex = 0;
    for iCorpo = 1:config.NCorpi
        for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
            for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
                rowIndex = rowIndex + 1;
                NormalHere = Normals{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
                TermineNoto(rowIndex) = -dot(U_Inf, NormalHere);
            end
        end
    end

    Solution = linsolve(matriceA, TermineNoto);

    for iCorpo = 1:config.NCorpi
        Lift = 0;
        dSpan = config.SemiSpan(iCorpo) / config.SemiSpanwiseDiscr(iCorpo);
        GammaHere = reshape(Solution((iCorpo-1)*config.ChordwiseDiscr(iCorpo)*2*config.SemiSpanwiseDiscr(iCorpo)+1:iCorpo*config.ChordwiseDiscr(iCorpo)*2*config.SemiSpanwiseDiscr(iCorpo)), config.ChordwiseDiscr(iCorpo), []);
        Lift = Lift + sum(sum(GammaHere)) * rho * U_Inf_Mag * dSpan;
        CL_values(iCorpo, i) = Lift / (0.5 * rho * U_Inf_Mag^2 * config.Surface(iCorpo));

        CDi = CL_values(iCorpo, i)^2 / (pi * config.AspectRatio(iCorpo) * e);
        CD_values(iCorpo, i) = CDi;
    end
end

% Plot CL vs CD for wing and tail
figure;
hold on;
for iCorpo = 1:config.NCorpi
    plot(CD_values(iCorpo, :), CL_values(iCorpo, :), 'o-', 'DisplayName', ['Body ' num2str(iCorpo)]);
end

xlabel('Drag Coefficient (C_D)');
ylabel('Lift Coefficient (C_L)');
title('Lift Coefficient vs Drag Coefficient for Wing and Tail');
legend('Profilo alare','Piani di coda');
grid on;
hold off;

% Define alpha range
alpha_range = linspace(-10, 10, 21);
CL_values = zeros(size(alpha_range));
CD_values = zeros(size(alpha_range));

for i = 1:length(alpha_range)
    alpha = alpha_range(i);
    U_Inf = U_Inf_Mag * [cosd(alpha)*cosd(beta), sind(beta), sind(alpha)*cosd(beta)];

    rowIndex = 0;
    for iCorpo = 1:config.NCorpi
        for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
            for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
                rowIndex = rowIndex + 1;
                NormalHere = Normals{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
                TermineNoto(rowIndex) = -dot(U_Inf, NormalHere);
            end
        end
    end

    Solution = linsolve(matriceA, TermineNoto);

    Lift = 0;
    for iCorpo = 1:config.NCorpi
        dSpan = config.SemiSpan(iCorpo) / config.SemiSpanwiseDiscr(iCorpo);
        GammaHere = reshape(Solution, config.ChordwiseDiscr(iCorpo), []);
        Lift = Lift + sum(sum(GammaHere)) * rho * U_Inf_Mag * dSpan;
    end

    CL_values(i) = Lift / (0.5 * rho * U_Inf_Mag^2 * sum(config.Surface));

    CDi = CL_values(i)^2 / (pi * config.AspectRatio(1) * e);
    CD_values(i) = CDi;
end

% Plot CL vs CD
figure;
plot(CD_values, CL_values, 'b-o');
xlabel('Drag Coefficient (C_D)');
ylabel('Lift Coefficient (C_L)');
title('Lift Coefficient vs Drag Coefficient');
grid on;


%% Compute and plot the lift distribution over the entire span
LiftDistribution = cell(config.NCorpi, 1);

for iCorpo = 1:config.NCorpi
    LiftDistribution{iCorpo} = sum(Gamma{iCorpo}, 1);
end

% Plot the circulation distribution over the span
figure(5);
hold on;
for iCorpo = 1:config.NCorpi
    y = linspace(-config.SemiSpan(iCorpo), config.SemiSpan(iCorpo), 2 * config.SemiSpanwiseDiscr(iCorpo));
    plot(y, LiftDistribution{iCorpo}, 'LineWidth', 1.5, 'DisplayName', ['Body ' num2str(iCorpo)]);
end

% Compute and plot the elliptic circulation distribution
for iCorpo = 1:config.NCorpi
    y = linspace(-config.SemiSpan(iCorpo), config.SemiSpan(iCorpo), 2 * config.SemiSpanwiseDiscr(iCorpo));
    Gamma_elliptic = max(LiftDistribution{iCorpo}) * sqrt(1 - (y / config.SemiSpan(iCorpo)).^2);
    plot(y, Gamma_elliptic, '--', 'LineWidth', 1.5, 'DisplayName', ['Elliptic Body ' num2str(iCorpo)]);
end

title('Circulation Distribution over the Span');
xlabel('Spanwise Position [m]');
ylabel('Circulation (\Gamma) [m^2/s]');
legend('Profilo alare','Piani di coda','Profilo alare ellittico','Piani di coda ellittici');
grid on;
hold off;

