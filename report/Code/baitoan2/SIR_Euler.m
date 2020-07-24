%TODO: Using Euler's algorithm to find solutions for SIR systems
%INPUT:
    % - t: vector contain time varaiable t
    % - beta: infectious coefficient
    % - gama: recovery coefficient
    % - S0: value of S at time t0
    % - I0: value of I at time t0
    % - R0: value of R at time t0
%OUTPUT: 
    % - S_t: vector contain value of S at correspond time in vector t
    % - I_t: vector contain value of I at correspond time in vector t
    % - R_t: vector contain value of R at correspond time in vector t
function [S_t,I_t,R_t]=SIR_Euler(t,beta,gama,S0,I0,R0)
    syms diff_S(I,S) diff_I(I,S) diff_R(I)
    
    diff_S(I,S) = -beta*I*S; % S'(t)
    diff_I(I,S) = beta*I*S - gama*I; % I'(t)
    diff_R(I) = gama*I; % R'(t)
  
    S_t = []; S_t(1) = S0; % initial for S_t
    I_t = []; I_t(1) = I0; % initial for I_t
    R_t = []; R_t(1) = R0; % initial for R_t
    
    % Euler's algorithm
    for i=2:length(t) 
        % S_n+1 = S_n + delta_t*S'(t_n)
        S_t(i) = S_t(i-1) + (t(i)-t(i-1))*diff_S(I_t(i-1),S_t(i-1)); 
        
        % I_n+1 = I_n + delta_t*I'(t_n)
        I_t(i) = I_t(i-1) + (t(i)-t(i-1))*diff_I(I_t(i-1),S_t(i-1)); 
        
        % R_n+1 = R_n + delta_t*R'(t_n)
        R_t(i) = R_t(i-1) + (t(i)-t(i-1))*diff_R(I_t(i-1)); 
    end 
end     