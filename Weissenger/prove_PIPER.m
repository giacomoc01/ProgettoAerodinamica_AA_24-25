close all
clear all %#ok<*CLALL>
clc
%% Test Case 1
U_Inf_Mag = 1;
beta = 0;
alpha = 5;
U_Inf = [cosd(alpha)*cosd(beta) sind(beta) sind(alpha)*cosd(beta)] .* U_Inf_Mag;
rho = 1.225;
config.NCorpi = 1;
config.RootChord = 1.5962; %[m]
config.DihedralAngle = 0; % NON SO SE VA BENE!!
config.SweepAngle = 0; % [°]
config.TaperRatio = 1; 
config.AspectRatio = 6.944; 
config.Span = 10.73;
config.LEPosition_X = 0;
config.LEPosition_Y = 0;
config.LEPosition_Z = 0;
config.RotationAngle_X =[0 0];
config.RotationAngle_Y = [0;0];
config.RotationAngle_Z = [0 0];
% Discretization options
config.SemiSpanwiseDiscr = 20;
config.ChordwiseDiscr = 20;

%% Preliminary computations
% Computing the span
config.SemiSpan = config.Span./2;
% Computing the surface
config.Surface = 2 * (config.SemiSpan .* config.RootChord .* ( 1+config.TaperRatio)./2);
config.SurfaceProjected = config.Surface .* cosd(config.DihedralAngle);
% Computing the Tip chord
config.TipChord = config.RootChord .* config.TaperRatio;
% Compute MAC
config.MAC = (2/3) .* config.RootChord .* ( (1 + config.TaperRatio + config.TaperRatio.^2)./(1 + config.TaperRatio));

%% Create the geometry structure
ControlPoints = cell(config.NCorpi, 1);
InducedPoints = cell(config.NCorpi, 1);
Normals = cell(config.NCorpi, 1);
InfiniteVortices = cell(config.NCorpi, 1);
Vortices = cell(config.NCorpi, 1);
internalMesh = cell(config.NCorpi, 1);
WingExtremes = cell(config.NCorpi, 1);
for iCorpo = 1:config.NCorpi
    [ControlPoints{iCorpo}, InducedPoints{iCorpo}, Normals{iCorpo}, InfiniteVortices{iCorpo}, Vortices{iCorpo}, internalMesh{iCorpo}, WingExtremes{iCorpo}] = createStructure(config, iCorpo);
end
    
%% Matrices initialization
NPanelsTot = 2* config.SemiSpanwiseDiscr * config.ChordwiseDiscr';
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

%% Costruzione del termine noto
rowIndex = 0;
for iCorpo = 1:config.NCorpi
    
    % Cycle on all of its chordwise panels
    for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
        % Cycle on all of its spanwise panels
        for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
            
            % Update row index
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
Lift3D = 0;
for iCorpo = 1:config.NCorpi
    dSpan = config.SemiSpan(iCorpo) / config.SemiSpanwiseDiscr(iCorpo);
    for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
        for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
            GammaHere = Gamma{iCorpo}(ChordPanel_i, SpanPanel_i);
            Lift3D = Lift3D + rho * U_Inf_Mag * GammaHere * dSpan;
        end
    end
end
disp(['3D Lift: ' num2str(Lift3D) ' N']);


%% Compute 2D Lift distribution
Lift2D = cell(config.NCorpi, 1);
for iCorpo = 1:config.NCorpi
    Lift2D{iCorpo} = rho * U_Inf_Mag * sum(Gamma{iCorpo}, 1) * config.MAC(iCorpo);
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
CL_3D = Lift3D / (q_inf * config.Surface);

%% Calcolo della resistenza indotta
% Fattore di Oswald
e = 0.85;

% Calcolo del coefficiente di resistenza indotta
CDi = CL_3D^2 / (pi * config.AspectRatio * e);

% Calcolo della resistenza indotta
Di = CDi * q_inf * config.Surface;

% Stampa dei risultati
fprintf('Coefficiente di portanza: CL = %.4f\n', CL_3D);
fprintf('Coefficiente di resistenza indotta: CDi = %.4f\n', CDi);
fprintf('Resistenza indotta: Di = %.4f N\n', Di);

% Calcolo dell'efficienza aerodinamica
E = CL_3D/CDi;
fprintf('Efficienza aerodinamica: E = %.4f\n', E);

% Define angle range in degrees
alpha_range = -10:1:15;
CL_values = zeros(size(alpha_range));
CD_values = zeros(size(alpha_range));

% Calculate CL for each angle
for i = 1:length(alpha_range)
    % Convert angle to radians for calculations
    alpha = deg2rad(alpha_range(i));
    
    % Recompute U_inf direction for each angle
    U_Inf = [cos(alpha) 0 sin(alpha)] * U_Inf_Mag;
    
    % Recalculate termine noto for this angle
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
    
    % Solve system for this angle
    Solution = linsolve(matriceA, TermineNoto);

    % Calculate lift for this angle
    Lift = 0;
    for iCorpo = 1:config.NCorpi
        dSpan = config.SemiSpan(iCorpo) / config.SemiSpanwiseDiscr(iCorpo);
        GammaHere = reshape(Solution, config.ChordwiseDiscr(iCorpo), []);
        Lift = Lift + sum(sum(GammaHere)) * rho * U_Inf_Mag * dSpan;
    end
    
    % Calculate CL
    % Apply sweep angle correction
    sweep_correction_factor = cos(deg2rad(config.SweepAngle));
    CL_values(i) = (Lift / (0.5 * rho * U_Inf_Mag^2 * config.Surface)) * sweep_correction_factor;
    
    % Calculate induced drag coefficient for this angle
    CDi = CL_values(i)^2 / (pi * config.AspectRatio * e);
    CD_values(i) = CDi;
end

% Calculate cl_alpha using linear regression
alpha_rad = deg2rad(alpha_range);
cl_alpha = (CL_values(end) - CL_values(1)) / (alpha_rad(end) - alpha_rad(1));

% Plot CL vs alpha to verify linearity
figure;
plot(alpha_range, CL_values, 'b-o');
xlabel('Angle of Attack (degrees)');
ylabel('C_L');
title('Lift Coefficient vs Angle of Attack');
grid on;

% Calcolo dei coefficienti di portanza
q_inf = 0.5 * rho * U_Inf_Mag^2;

% CL 3D
CL_3D = Lift3D / (q_inf * config.Surface);


% CL 2D per ogni sezione
Lift2D = cell2mat(Lift2D);
for iCorpo = 1:config.NCorpi
    CL_2D = Lift2D(iCorpo) / (q_inf * config.MAC(iCorpo));
    fprintf('Coefficiente di portanza 2D sezione %d: CL = %.4f\n', iCorpo, CL_2D);
end

% Calculate CD for each angle
for i = 1:length(alpha_range)
    % Calculate induced drag coefficient for this angle
    CDi = CL_values(i)^2 / (pi * config.AspectRatio * e);
    CD_values(i) = CDi;
end

%% Plot CL vs CD for different sweep angles
% Define sweep angles to analyze
sweep_angles = [0, 15, 30];
colors = ['b', 'g', 'r'];
legend_entries = cell(1, length(sweep_angles));

figure;
hold on;

for k = 1:length(sweep_angles)
    config.SweepAngle = sweep_angles(k);
    CL_values = zeros(size(alpha_range));
    CD_values = zeros(size(alpha_range));
    
    for i = 1:length(alpha_range)
        alpha = deg2rad(alpha_range(i));
        U_Inf = [cos(alpha) 0 sin(alpha)] * U_Inf_Mag;
        
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
        
        sweep_correction_factor = cos(deg2rad(config.SweepAngle));
        CL_values(i) = (Lift / (0.5 * rho * U_Inf_Mag^2 * config.Surface)) * sweep_correction_factor;
        
        CDi = CL_values(i)^2 / (pi * config.AspectRatio * e);
        CD_values(i) = CDi;
    end
    
    plot(CD_values, CL_values, [colors(k) '-o'], 'DisplayName', ['Sweep Angle = ' num2str(config.SweepAngle) '°']);
    legend_entries{k} = ['Sweep Angle = ' num2str(config.SweepAngle) '°'];
end

xlabel('Drag Coefficient (C_D)');
ylabel('Lift Coefficient (C_L)');
title('Lift Coefficient vs Drag Coefficient for Different Sweep Angles');
legend(legend_entries, 'Location', 'Best');
grid on;
hold off;


%% Plot CL vs alpha for different sweep angles
% Define sweep angles to analyze
sweep_angles = [0, 15, 30];
colors = ['b', 'g', 'r'];
legend_entries = cell(1, length(sweep_angles));

figure;
hold on;

for k = 1:length(sweep_angles)
    config.SweepAngle = sweep_angles(k);
    CL_values = zeros(size(alpha_range));
    
    for i = 1:length(alpha_range)
        alpha = deg2rad(alpha_range(i));
        U_Inf = [cos(alpha) 0 sin(alpha)] * U_Inf_Mag;
        
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
        
        sweep_correction_factor = cos(deg2rad(config.SweepAngle));
        CL_values(i) = (Lift / (0.5 * rho * U_Inf_Mag^2 * config.Surface)) * sweep_correction_factor;
    end
    
    plot(alpha_range, CL_values, [colors(k) '-o'], 'DisplayName', ['Sweep Angle = ' num2str(config.SweepAngle) '°']);
    legend_entries{k} = ['Sweep Angle = ' num2str(config.SweepAngle) '°'];
end

xlabel('Angle of Attack (degrees)');
ylabel('Lift Coefficient (C_L)');
title('2D Lift Coefficient vs Angle of Attack for Different Sweep Angles');
legend(legend_entries, 'Location', 'Best');
grid on;
hold off;


%% Plot CL vs CD for different taper ratios
% Define taper ratios to analyze
taper_ratios = [0, 2, 4];
colors = ['b', 'g', 'r'];
legend_entries = cell(1, length(taper_ratios));

figure;
hold on;

for k = 1:length(taper_ratios)
    config.TaperRatio = taper_ratios(k);
    config.TipChord = config.RootChord * config.TaperRatio;
    config.Surface = 2 * (config.SemiSpan * config.RootChord * (1 + config.TaperRatio) / 2);
    config.SurfaceProjected = config.Surface * cosd(config.DihedralAngle);
    config.MAC = (2/3) * config.RootChord * ((1 + config.TaperRatio + config.TaperRatio^2) / (1 + config.TaperRatio));
    
    CL_values = zeros(size(alpha_range));
    CD_values = zeros(size(alpha_range));
    
    for i = 1:length(alpha_range)
        alpha = deg2rad(alpha_range(i));
        U_Inf = [cos(alpha) 0 sin(alpha)] * U_Inf_Mag;
        
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
        
        CL_values(i) = Lift / (0.5 * rho * U_Inf_Mag^2 * config.Surface);
        
        CDi = CL_values(i)^2 / (pi * config.AspectRatio * e);
        CD_values(i) = CDi;
    end
    
    plot(CD_values, CL_values, [colors(k) '-o'], 'DisplayName', ['Taper Ratio = ' num2str(config.TaperRatio)]);
    legend_entries{k} = ['Taper Ratio = ' num2str(config.TaperRatio)];
end

xlabel('Drag Coefficient (C_D)');
ylabel('Lift Coefficient (C_L)');
title('Lift Coefficient vs Drag Coefficient for Different Taper Ratios');
legend(legend_entries, 'Location', 'Best');
grid on;
hold off;

%% Plot CL vs alpha for different taper ratios
figure;
hold on;

for k = 1:length(taper_ratios)
    config.TaperRatio = taper_ratios(k);
    config.TipChord = config.RootChord * config.TaperRatio;
    config.Surface = 2 * (config.SemiSpan * config.RootChord * (1 + config.TaperRatio) / 2);
    config.SurfaceProjected = config.Surface * cosd(config.DihedralAngle);
    config.MAC = (2/3) * config.RootChord * ((1 + config.TaperRatio + config.TaperRatio^2) / (1 + config.TaperRatio));
    
    CL_values = zeros(size(alpha_range));
    
    for i = 1:length(alpha_range)
        alpha = deg2rad(alpha_range(i));
        U_Inf = [cos(alpha) 0 sin(alpha)] * U_Inf_Mag;
        
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
        
        CL_values(i) = Lift / (0.5 * rho * U_Inf_Mag^2 * config.Surface);
    end
    
    plot(alpha_range, CL_values, [colors(k) '-o'], 'DisplayName', ['Taper Ratio = ' num2str(config.TaperRatio)]);
    legend_entries{k} = ['Taper Ratio = ' num2str(config.TaperRatio)];
end

xlabel('Angle of Attack (degrees)');
ylabel('Lift Coefficient (C_L)');
title('2D Lift Coefficient vs Angle of Attack for Different Taper Ratios');
legend(legend_entries, 'Location', 'Best');
grid on;
hold off;

%% Compute and plot the lift distribution over the entire span
LiftDistribution = zeros(1, 2 * config.SemiSpanwiseDiscr);

for iCorpo = 1:config.NCorpi
    for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
        LiftDistribution(SpanPanel_i) = sum(Gamma{iCorpo}(:, SpanPanel_i));
    end
end

% Plot the circulation distribution over the span
figure;
y = linspace(-config.SemiSpan, config.SemiSpan, 2 * config.SemiSpanwiseDiscr);
plot(y, LiftDistribution, 'b-', 'LineWidth', 1.5);
hold on;

% Compute and plot the elliptic circulation distribution
Gamma_elliptic = max(LiftDistribution) * sqrt(1 - (y / config.SemiSpan).^2);
plot(y, Gamma_elliptic, 'r--', 'LineWidth', 1.5);

title('Circulation Distribution over the Span');
xlabel('Spanwise Position [m]');
ylabel('Circulation (\Gamma) [m^2/s]');
legend('Computed Circulation', 'Elliptic Circulation');
grid on;
hold off;

