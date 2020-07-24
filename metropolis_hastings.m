function ret = metropolis_hastings(lambda_beta, lambda_gamma, v_beta, v_gamma,  N)    
%TODO: Create a sample of size N (beta, gamma)
%{
PARAM
    init_beta: initialized beta value
    init_gamma: initialized beta value
    N: Size of sample
OUTPUT
    Nx3 matrix:
        Column 1: beta
        Column 2: gamma
        Column 3: pi(beta, gamma)
%}

%Distribution used for sampling is Normal Distribution.

    beta_param = [lambda_beta, v_beta]; %%[shape, rate]
    gamma_param = [lambda_gamma v_gamma]; %%[shape, rate]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Define some lambda functions:%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %pi(beta, gamma): density function
    pi = @(X) ...
        gampdf(X(1), beta_param(1), 1/beta_param(2))...
        *gampdf(X(2), gamma_param(1), 1/gamma_param(2));
    %p(a, b | a*, b*)
    %p = @(q, ev) normpdf(q(1), ev(1), 1)*normpdf(q(2), ev(2), 1);
    p = @(q, ev, var1, var2) mvnpdf(q, ev, [var1 0; 0 var2]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %%%%%%%%%%%%%%%%%%%%%%
    %Metropolis-Hastrings%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ret = zeros(N,3);
    %a = abs(mvnrnd([0 0], [1/10000 0; 0 1/10000]));
    %ret(1,:) = [a(1), a(2), pi([a(1), a(2)])];
    
    ret(1,:) = [lambda_beta/v_beta, lambda_gamma/v_gamma, pi([lambda_beta/v_beta, lambda_gamma/v_gamma])];
    accept_times = 1;
    
    for i = 1 : N-1
        %Create random beta*, gamma* using normal distribution which mean value is
        %beta, gamma, respectively.
        var = accept_times/(N - accept_times);
        current = [ret(i,1), ret(i,2)];
        proposal = mvnrnd([current(1) current(2)], [var 0; 0 var]);
        %Calculate acceptance ratio: Normal Distribution is symmetric.
        ratio = pi(proposal)/pi(current);
        %Accept new candidate
        if(rand < min(1, ratio))
            accept_times = accept_times + 1;
            ret(i + 1,:) = [proposal(1), proposal(2), pi(proposal)];
        %Reject new candidate
        else
            ret(i + 1,:) = ret(i,:);
        end
    end    
    fprintf('Acceptance rate: %f percents\n', 100*accept_times/N);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%
    %Plot%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(1)
    subplot(2,3,1)
    x = ret(:,1);
    y = ret(:,2);
    scatter(x,y)
    xlabel('beta');
    ylabel('gamma');
    title('Scatter plot of (beta, gamma)')
    
    %figure(2)
    subplot(2,3,2)
    x = linspace(1,N,N);
    y = ret(:,1);
    plot(x,y)
    xlabel('iteration');
    ylabel('beta');
    title('Value of beta through each iteration');

    %figure(3)
    subplot(2,3,3)
    x = linspace(1,N,N);
    y = ret(:,2);
    plot(x,y)
    xlabel('iteration');
    ylabel('gamma');
    title('Value of gamma through each iteration');
    
    %figure(4)
    subplot(2,3,4)
    b = containers.Map('KeyType','int32','ValueType','int32');
    for i = 1 : N
        if(isKey(b, ret(i,1)))
            b(ret(i,1)) = b(ret(i,1)) + 1;
        else
            b(ret(i,1)) = 1;
        end
    end
    x = cell2mat(keys(b));
    y = cell2mat(values(b));
    plot(x,y)
    xlabel('beta');
    ylabel('frequency');
    title('Frequency plot of beta');

    %figure(5)
    subplot(2,3,5)
    g = containers.Map('KeyType','int32','ValueType','int32');
    for i = 1 : N
        if(isKey(g, ret(i,2)))
            g(ret(i,2)) = g(ret(i,2)) + 1;
        else
            g(ret(i,2)) = 1;
        end
    end
    x = cell2mat(keys(g));
    y = cell2mat(values(g));
    plot(x,y)
    xlabel('gamma');
    ylabel('frequency');
    title('Frequency plot of gamma');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    