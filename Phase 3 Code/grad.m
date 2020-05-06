
function [gradf]=grad(fun,x,R)
%{
This function returns the gradient vector of the function @fcn using a
central difference approximation. Inputs are the function @fcn,
x=[x1,..xn].
%}

dx=0.001;
for i=1:length(x)
    
    xtemp1=x;
    xtemp2=x;
    xtemp1(i)=xtemp1(i)+dx;
    xtemp2(i)=xtemp2(i)-dx;
    pfpx(i)=(fun(xtemp1,R)-fun(xtemp2,R))./(2*dx);
    
end
gradf=pfpx;
end