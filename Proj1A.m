function Proj1A
% Conor Luppnow
% This program:
%   - Determines the paramater A in the two-suffix Margules equation.
%   - Evaluates gamma_1 and gaamma_2 for each x1, and plots this data.
%   - Evalulates and plots ln(gamma_1/gamma_2) vs. x1 for use in the
%     thermodynamic consistency test, as well as the Wisniak test.
%   - Allows for the conclusion that the VLE data satisfies Wisniak’s test
%     criterion for the thermodynamic consistency test.
%   
    clc; clear; close;
    % set up variables and data
    x1 = 0:0.05:1;
    y1 = [0 0.275 0.436 0.543 0.618 0.675 0.720 0.756 0.786 0.811 0.833 0.853 0.871 0.888 0.903 0.918 0.933 0.949 0.964 0.981 1];
    P_exp = [7295 9562 11695 13682 15536 17256 18883 20390 21817 23150 24404 25604 26751 27845 28898 29925 30938 31939 32952 33979 35032];
    T = 313.15;
    R = 8.314;

    %function handle for minimizing P_exp vs P_calc
    fun = @(X) P(x1,X,P_exp);

    % Use fminsearch to obtain A paramater in two-suffix Margules model.
    options = optimset('Display', 'iter', 'MaxFunEvals', 1E10, 'MaxIter', 1E10);
    x0 = rand(1,2);
    [x, f, exitflag] = fminsearch(fun, x0, options);

    if exitflag > 0
        A = x(1);

    end
    fprintf('A = %5.5f and f = %5.5f \n',A,f);

    % calculate gamma_1,gamma_2, and ln(gamma_1/gamma_2)
    for j = 1:numel(x1)
        gamma_1(j) = exp((A/(R*T))*(1-x1(j))^2);
        gamma_2(j) = exp((A/(R*T))*(x1(j))^2);
        ln_gamma_ratio(j) = log(gamma_1(j) / gamma_2(j));
    end

    % Finally plot the data
    plot(x1,gamma_1,'or'); hold on
    plot(x1,gamma_2,'ok');  grid on

    plot(x1,ln_gamma_ratio,'ob');

    xlabel('x1');
    ylabel('\gamma_i');
    legend('\gamma_1','\gamma_2','ln(\gamma_1 / \gamma_2)','Location','North');

    % calculatea area to determine accuracy of VLE data.
    area = trapz(x1,ln_gamma_ratio);

    area_abs = trapz(x1,abs(ln_gamma_ratio));

    fprintf('area = %5.5f\n',area);

    % Evaluate D within the Wisniak test.

    D = 100 * (area/area_abs);

    fprintf('D = %5.5f\n',D);

end

% function to calculate A by minimizing the residuals of P_exp and P_calc.
function sser = P(x1,X,P_exp)
    P1_sat = 35430;
    P2_sat = 7.3837E3;

    T = 313.15;
    R = 8.314;
    A = X(1);
    sser = 0;

    for i = 1:numel(x1)
        P_calc(i) = x1(i)*exp((A/(R*T))*(1-x1(i))^2)*P1_sat + (1-x1(i)) * exp((A/(R*T))*(x1(i))^2) * P2_sat;

        sser = sser + (P_exp(i) - P_calc(i))^2;
    end
end





