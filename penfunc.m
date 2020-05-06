%Epsilons
e1 = 1e-3;
e2 = 1e-3;
e3 = 1e-3;
e4 = 1e-3;

% Variable bounds
k = 1;
if k == 1
    N = 2;
    a = [13 100];
    b = [0 100];
    bnd = [a;b];
    x0 = [rand*(100-13)+13 rand*100]';
    
elseif k == 2
    N = 2;
    a = [0 10];
    b = [0 10];
    bnd = [a;b];
    x0 = [rand*5+5  rand*5+5]';
    
elseif k == 3
    N = 8;
    a = [100 10000];
    b = [1000 10000];
    c = [1000 10000];
    d = [10 1000];
    e = [10 1000];
    f = [10 1000];
    g = [10 1000];
    h = [10 1000];
    bnd = [a;b;c;d;e;f;g;h];
    x0 = [580 1400 5000 200 300 200 300 400]';
    
end

% Penalty Parameter
R_0 = 0.1;
k = 1; 

f_x0 = fx(x0);      

% Penalty function value at initial point
P_x0_R0 = pf(x0,R_0);    

% constraint 
cv = [g1(x0) g2(x0)]; 

%% Conjugate Gradient to calculate next point 
x1 = CG(x0,R_0,bnd,N,e1,e2,e3);         
f_x1 = fx(x1);      
P_x1_R0 = pf(x1,R_0);   
cv = [g1(x1) g2(x1)]; 

%print_seq(k,R_0,x1,P_x1_R0,CV);

k = k + 1;
R1 = 10*R_0;
P_x1_R1 = pf(x1,R1);   
cv = [g1(x1) g2(x1)];

x2 = DFP(x1,R1,bnd,N,e1,e2,e3);         %CG to calculate next point using x1
P_x2_R1 = pf(x2,R1);    
cv = [g1(x2) g2(x2)];

while abs(P_x2_R1 - P_x1_R0) > e4         %% termination condition
R_0 = R1;
P_x1_R0 = P_x2_R1;    
x1 = x2; 


R1 = 10*R_0;

k = k +1;
P_x1_R1 = pf(x1,R1);   
cv = [g1(x1) g2(x1)]; 

x2 = DFP(x1,R1,bnd,N,e1,e2,e3);         
cv = [g1(x2) g2(x2)]; 
P_x2_R1 = pf(x2,R1);    
end
fprintf("optima: (%.3f,%.3f)\n",x2(1),x2(2));
fprintf("f(optima) = %.4f\n",fx(x2));