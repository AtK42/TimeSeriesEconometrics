one_Vector = ones(100, 1);
sigma_squared = 2;
R_squared_values = []
for k = 2:10
    true_beta = zeros(k, 1);
    R_squared_empirical = []
    for i = 1:1e4
        X_cols = normrnd(0, 1, 100, k-1);
        X = [one_Vector X_cols];
        epsilon = normrnd(0, sqrt(sigma_squared), 100, 1);
        Y = X*true_beta + epsilon;
        beta_hat = inv(X'*X)*X'*Y;
        Y_bar = mean(Y);
        Y_hat = X*beta_hat;
        ESS = sum((Y_hat - Y_bar).^2);
        TSS = sum((Y - Y_bar).^2);
        R_squared_empirical(i) = ESS/TSS;
        i + (10000)*(k-2)
    end
    R_squared_values(:, k-1) = R_squared_empirical;
end
figure;
boxplot(R_squared_values)