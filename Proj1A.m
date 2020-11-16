function Proj1A
clc;
x1 = 0:0.05:1;
y1 = [0 0.275 0.436 0.543 0.618 0.675 0.720 0.756 0.786 0.811 0.833 0.853 0.871 0.888 0.903 0.918 0.933 0.949 0.964 0.981 1];
P_exp = [7295 9562 11695 13682 15536 17256 18883 20390 21817 23150 24404 25604 26751 27845 28898 29925 30938 31939 32952 33979 35032];


fun = @(X) P(x1,X,P_exp);

options = optimset('Display', 'iter', 'MaxFunEvals', 1E10, 'MaxIter', 1E10);
x0 = rand(1,2);
[x, f, exitflag] = fminsearch(fun, x0, options);

if exitflag > 0
    A = x(1);
    
end
fprintf('A = %5.5f and f = %5.5f \n',A,f);


end


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





