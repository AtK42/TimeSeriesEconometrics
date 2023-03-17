%Specify samplesize and a
samplesize = 50;
a = 0.99;
true_beta = [1 2 3 4 5]';
%Run everything twice, once for T=20 and once for T=50
    one_Vector = ones(samplesize + 1, 1);
    X_cols = normrnd(0, 1, samplesize + 1, 4);
    X = [one_Vector X_cols];
    a_vector = [];
    beta_matrix =[];
    %Run the 1e4 iterations
    for i = 1:1e4
        epsilonzero = normrnd(0,1/(1-a^2));
        epsilon = [];
        epsilon(1) = epsilonzero;
        %Simulate the errors
        for j = 2:samplesize + 1
            epsilon(j) = a*epsilon(j-1) + randn;
        end
        Y = X*true_beta + epsilon';
        i
        %Iterate through the algorithm from page 225
        sigma = [];
        beta_estimate = X\Y;
        epsilon_estimate = Y - X*beta_estimate;  
        a_estimate = epsilon_estimate(1:samplesize)\epsilon_estimate(2:samplesize + 1);
        olda = inf;
        beta_estimate_old = inf;
        while max([abs(a_estimate - olda) abs(beta_estimate' -beta_estimate_old')]) > 1e-5 
            %Now we need to create the sigma matrix
            for row = 1:(samplesize + 1)
                for col = 1:(samplesize + 1)
                    sigma(row, col) = (a_estimate^(abs(row-col)))/(1 - a_estimate.^2);
                end
            end
            beta_estimate_old = beta_estimate;
            beta_estimate = inv(X'*inv(sigma)*X)*X'*inv(sigma)*Y;
            epsilon_estimate = Y - X*beta_estimate;
            olda = a_estimate;
            a_estimate = epsilon_estimate(1:samplesize)\epsilon_estimate(2:samplesize+1);
        end
        a_vector(i, 1) = a_estimate;
        beta_matrix(i, :) = beta_estimate;
    end
    figure;
    boxplot(a_vector - a)
     yline(0)
     str1 = ['Errors of a-estimates with true a = '  num2str(a) ' and samplesize ' num2str(samplesize)]
    title(str1)
    figure;
    boxplot(beta_matrix - true_beta');
    hold on;
    yline(0)
    str2 = ['Errors of the beta-estimate with true a = ' num2str(a) ' and samplesize ' num2str(samplesize)]
    title(str2)
    mean(a_vector-a)