%Generate the normal-distributed data matrix
one_Vector = ones(50, 1);
X_cols = [normrnd(0, 1, 50, 1) normrnd(0, 1, 50, 1) normrnd(0, 1, 50, 1)];
X = [one_Vector X_cols];

true_Beta = [1 1 1 1]';
no_Errors_Y = X*true_Beta;
r_squared_x_values = linspace(0,1, 1000);
f_x_values = linspace(0,100, 1000);
sigma_squared = [1 2 4];
R_squared_values = []
F_values = []

%Outer loop goes through the three sigma-squared values, while inner loop
%simulates the error terms and R-squared values
for j = 1:3
    R_squared_empirical = []
for i = 1:1e5
    epsilon = normrnd(0, sqrt(sigma_squared(j)), 50, 1);
    Y = no_Errors_Y + epsilon;
    beta_hat = inv(X'*X)*X'*Y;
    Y_bar = mean(Y);
    Y_hat = X*beta_hat;
    ESS = sum((Y_hat - Y_bar).^2);
    TSS = sum((Y - Y_bar).^2);
    R_squared_empirical(i) = ESS/TSS;
    i + 1e5*(j-1)
end

R_squared_values(:, j) = R_squared_empirical;
F_values(:,j) = ((50-4)/(4-1))*R_squared_empirical./(1-R_squared_empirical);

end
figure;
 ksdensity(R_squared_values(:, 1), r_squared_x_values)
hold on;
 ksdensity(R_squared_values(:,2), r_squared_x_values)
 ksdensity(R_squared_values(:,3), r_squared_x_values)
legend('Sigma squared 1 pdf', 'Sigma squared 2 pdf', 'Sigma squared 4 pdf');

figure;
ksdensity(F_values(:,1), f_x_values)
hold on;
ksdensity(F_values(:,2), f_x_values)
ksdensity(F_values(:,3), f_x_values)
legend('Sigma squared 1 pdf', 'Sigma squared 2 pdf', 'Sigma squared 4 pdf');