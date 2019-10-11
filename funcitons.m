function [f]=funcitons(x)
 
%f = 0.01*x^5 -2*x^4 + 500*(x-2)^2;
%f = x^3 + 5*x^2-3;                     % min in (-1,1) x0 = 0.5
f = -((2*x-5)^4 - (x^2-1)^3);          % max in (-10,0)
%f = -(8 + x^3 - 2*x - 2*exp(x));       % max in (-2,1)
%f = 3*x^2 + 12/x^3 - 5;                % min in (0.5,5)
%f = -(4*x*sin(x));                     % max in (0.5,3.14)  x0 = 1
%f = 2*(x-3)^2 + exp(0.5*x^2);       % min in (-2,3)   **
%f = 2*exp(x) - x^3 - 10*x;             % min in (0,4)
%f = x^2 - 10*exp(0.1*x);               % min in (-6,6)
%f = -(20*sin(x) - 15*x^2);             % max in (-4,4)
%f = abs(exp(x) - x^3);                 % find one root , x0 = 5

