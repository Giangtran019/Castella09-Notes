% Load the dataset
data = readtable('/MATLAB Drive/entrepreneurship_decision_dataset.csv');

% Extract relevant columns
riskTolerance = data.Risk_Tolerance;
investmentAvailability = data.Investment_Availability;
growthRate = data.Growth_Rate;

% Combine into matrix for correlation
X = [riskTolerance, investmentAvailability, growthRate];

% Calculate Pearson correlation coefficients
corrMatrix = corr(X, 'Rows', 'complete');

% Display correlation matrix
disp('Correlation Matrix:');
disp(array2table(corrMatrix, ...
    'VariableNames', {'Risk_Tolerance', 'Investment_Availability', 'Growth_Rate'}, ...
    'RowNames', {'Risk_Tolerance', 'Investment_Availability', 'Growth_Rate'}));

% Visualize the correlation matrix as a heatmap
figure;
heatmap({'Risk_Tolerance', 'Investment_Availability', 'Growth_Rate'}, ...
        {'Risk_Tolerance', 'Investment_Availability', 'Growth_Rate'}, ...
        corrMatrix, ...
        'Colormap', parula, ...
        'ColorbarVisible', 'on');
title('Correlation Matrix Heatmap');