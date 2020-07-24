%We will call the value that is proportional with E(R0), in (21), E*.
%E* = Sigma(pi(X | beta_i, gamma_i)*pi(beta_i, gamma_i)).
function E_R0s = R0_i(file, m)
    [I, R] = readFile(file);
    n = size(I,2) - 3;
    disp(n)
    weeks = ceil(n / 7);
    E_R0s = zeros(weeks,1);
    E_R0_log = zeros(weeks, 1);
  
    %Estimate E* in each weeek, start from March 6th in 'Singapore.csv'.
    for i = 1 : weeks
        E_R0s(i) = R_i(I, R, i, m);
        E_R0_log(i) = 1/abs(log(E_R0s(i)));
    end
    
    %%%Plot%%%
    figure(3)
    X1 = linspace(0, n, weeks + 1);
    Y1 = zeros(1, weeks + 1);
    for i = 2 : weeks + 1
        s = 7*(i - 2) + 1;
        e = 7*(i - 1) + 1;
        Y1(i) = mean((I(s + 1 : e) - I(s : e - 1) + 0.01) / (R(s + 1 : e) - R(s : e - 1) + 0.01));
        %Y1(i) = (I(e) - I(s))/(R(e) - R(s));
    end
    Y1 = Y1/max(Y1);
    
    X2 = linspace(0, n, weeks + 1);
    Y2 = zeros(1, weeks + 1);
    Y2(2:weeks + 1) = E_R0_log/max(E_R0_log);
    
    plot(X1,Y1);       
    hold on
    plot(X2, Y2, 'r');
    hold off
    legend('dI/dR', 'E(R0)');
end

%Estimate E* given the week.
function E_R0 = R_i(I, R, week, m)
    population = 5850343; %Singapore population
    n = 7;
    %Calculate Sigma(pi(X|b_i, c_i)*b/c)
    numerator = 0;
    %Calculate Sigma(pi(X))
    %pi_X = 0;
    
    %Calculate start (s) date and last date (e) of the week
    s = n*(week - 1) + 1;
    e = n*week;
    if e > size(I,2)
        e = size(I,2);
    end
    
    %Estimate beta, gamma in each day
    beta_mat = zeros(s - e);
    gamma_mat = zeros(s - e);
    for i = s + 1 : e
        beta_mat(i - s) = (I(i) + R(i) - I(i - 1) - R(i - 1))*population...
                        /(I(i - 1) * (population - I(i - 1) - R(i - 1)));
        gamma_mat(i - s) = (R(i) - R(i - 1))/I(i - 1); 
    end
    
    fprintf('I: %d -> %d\n', I(s), I(e));
    fprintf('R: %d -> %d\n', R(s), R(e));
    fprintf('Prior [b g]: %f %f with R0 = %f\n', mean(beta_mat), mean(gamma_mat), mean(beta_mat)/mean(gamma_mat));
    
    %Sampling based on mean of beta_mat and gamma_mat.
    %n = a for beta, = b for gamma.
    a = mean(beta_mat)*0.01;
    b = mean(gamma_mat)*0.01;
    sample = metropolis_hastings(mean(beta_mat)*a, mean(gamma_mat)*b, a, b, m);
    
    for i = 1 : size(sample, 1)
        %Calculate pi(X|beta_i, gamma_i), X(t(i)) = I(i).
        pi = prod(gampdf(I(s:e), sample(i,1), 1/sample(i,2)));
        numerator = numerator + pi*(sample(i,1)/sample(i,2));
        %pi_X = pi_X + pi;
    end
    %E_R0 = numerator/pi_X;
    E_R0 = numerator;
    fprintf('e^%f\n\n', log(E_R0));
    
end

function [I,R] = readFile(file)
    A = csvread(file, 1, 45);
    I = A(1,:);
    R = A(2,:);
end

