% Read the file
filename = 'USA_a2_r200000_n9';
fileID = fopen(filename, 'r');

% Check if the file was opened successfully
if fileID == -1
    error('Error opening the file. Please check the file name and path.');
end

% Read the data, ignoring lines starting with '#'
data = textscan(fileID, '%f %f', 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'CommentStyle', '#');
fclose(fileID);

% Convert cell array to matrix
data = cell2mat(data);

% Split data into two sets
sep = round(length(data(:,1)) / 2) ;
data1 = data(1:213, :);
data2 = data(213+1:end, :);

% Remove points where x > 1.02
data1 = data1(data1(:,1) <= 110, :);
data2 = data2(data2(:,1) <= 110, :);

% Create the plot
figure;

% Plot the first set of data
plot(data1(:,1), data1(:,2), '-', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'b', 'DisplayName', 'Set 1: Dorso');

% Hold the plot to overlay the second set
hold on;

% Plot the second set of data
plot(data2(:,1), data2(:,2), '-', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'DisplayName', 'Set 2: Ventre');

% Add grid and axis labels
grid on;
title('USA25  Re=200000   N=9   Alpha=2', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('x/c (Posizione relativa lungo la corda)', 'FontSize', 12);
ylabel('C_f (Coefficiente di Attrito)', 'FontSize', 12);

% Add a legend
legend('Location', 'Best', 'FontSize', 11);

% Set axis properties for better visualization
set(gca, 'FontSize', 12, 'LineWidth', 1);
axis tight;

% Save the plot as a high-quality image
saveas(gcf, 'NACAAAAA.jpg');
