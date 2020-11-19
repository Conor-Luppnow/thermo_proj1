function Proj1B
% Conor Luppnow
% This program:
%   - Calculates gE_exp vs x1 using VLE data.
%   - Determines the Wilson equation paramaters: lambda_12 and lambda_21
%   - Calculates gE_calc using the Wilson equation.
%   - plots gE_exp and gE_calc vs x1. 
%   - Calculates error (gE_exp - gE_calc)
%   - Allows for conclusion that the Wilson equation describes the
%     experimental data well.
%
    clc; clear; close;
    
    % set up variables and data.
    x1 = 0:0.05:1;
    x2 = 1-x1;
    T = 313.15;
    R_const = 8.314;
    v1 = 41.489;
    v2 = 18.156;
    
    % Set the model "mdlfunc" for use in nlinfit
    % b(1) = lambda_12
    % b(2) = lambda_21
    mdlfunc = @(b,x1) -R_const*T*(x1.* log(x1 + ((v2/v1)*exp(-b(1)/(R_const*T))).* x2) + x2.* log(x2 + ((v1/v2)*exp(-b(2)/(R_const*T))).* x1));

    % From 1A: A = 1286.46130
    A = 1286.46130;
    
    % calculate gE_exp using two-suffix Margules model.
    for i = 1:numel(x1)
        gE_exp(i) = A*x1(i)*x2(i);

    end

    % Use nlinfit to fit the data to the model mdlfunc.
    b0 = rand(1,2);

    [b,R,J,CovB,MSE] = nlinfit(x1,gE_exp,mdlfunc,b0);

    % calculate gE_calc using the Wilson equation.
    for j = 1:numel(x1)
        gE_calc(j) = -R_const*T*(x1(j)* log(x1(j) + ((v2/v1)*exp(-b(1)/(R_const*T)))* x2(j)) + x2(j)* log(x2(j) + ((v1/v2)*exp(-b(2)/(R_const*T)))* x1(j)));

    end

    % calculate the error (gE_exp - gE_calc)
    for k = 1:numel(x1)
        error(k) = gE_exp(k) - gE_calc(k);
    end

    % plot data, and print error values.
    plot(x1,gE_calc,'-k'); hold on
    plot(x1,gE_exp,'db'); grid on;
    xlabel('x1');
    ylabel('gE');
    legend('gE_c_a_l_c','gE_e_x_p')

    fprintf('Error: (gE_exp - gE_calc), for each x1\n');
    for q = 1:numel(x1)
        fprintf('For x1: %5.3f, error:%5.3f \n',x1(q),error(q));
    end
    fprintf('\n');
end