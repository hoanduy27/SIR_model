%% Description %%
% TODO: Using Euler's algorithm to find solutions for SIR systems
% INPUT:
    % - t: time varaiable t
    % - beta: contact coefficient
    % - gamma: recovery coefficient
    % - S0: S(t0)
    % - I0: I(t0)
    % - R0: R(t0)
% OUTPUT: 
    % - ret: vector [S(t) I(t) R(t)]
function ret=SIR_Euler(t,beta,gamma,S0,I0,R0)
    %% Initial %%
    syms diff_S(I,S) diff_I(I,S) diff_R(I)
    N = S0+I0+R0; 
    diff_S(I,S) = -beta/N*I*S; % S'(t)
    diff_I(I,S) = beta/N*I*S - gamma*I; % I'(t)
    diff_R(I) = gamma*I; % R'(t)
    S_n = []; S_n(1) = S0; % initial for S_n series 
    I_n = []; I_n(1) = I0; % initial for I_n series
    R_n = []; R_n(1) = R0; % initial for R_n series
    
    %% Euler's algorithm %%
    step_size = 0.1 ; % We chose step size (delta t) = 0.1
    t_n = [0:step_size:t]; % t_n series
    for i=2:length(t_n) 
        % S_n+1 = S_n + delta_t*S'(t_n)
        S_n(i) = S_n(i-1) + step_size*diff_S(I_n(i-1),S_n(i-1)); 
        % I_n+1 = I_n + delta_t*I'(t_n)
        I_n(i) = I_n(i-1) + step_size*diff_I(I_n(i-1),S_n(i-1)); 
        % R_n+1 = R_n + delta_t*R'(t_n)
        R_n(i) = R_n(i-1) + step_size*diff_R(I_n(i-1)); 
    end 
    
    %% Return result %%
    ret=[];
    ret(1) = S_n(length(S_n)); % ret(1)=S(t)
    ret(2) = I_n(length(I_n)); % ret(2)=I(t)
    ret(3) = R_n(length(R_n)); % ret(3)=R(t)

    %% Draw Graph %%
    plot(t_n,S_n,'-b',t_n,I_n,'-r',t_n,R_n,'-g')
    legend('S(t)', 'I(t)', 'R(t)')
    xlabel('Time')
    ylabel('Number of people')
    set(gca,'XTick',0:1:t)
end     
