
%% Penalty Function Method

clc;clear;
tol1 = 10e-5;       % tolerance or accuracy
tol2 = 10e-5;
R0 = 0.1;           % initial Penalty parameter
 %x0 =[rand(0,2), rand(0,2)];
 x0 = [20 5];
c = input('enter 0 for interior penalty, 1 for exterior penalty = ');
 if c==0
     c = 0.5;              % c to update penalty parameter
 else
     c = 10;
 end
 t = 1;
 x(t,:) = x0;
 R = R0;
 x1 = x(1);
 x2 = x(2);

 % constraints for bracket operator penalty
 
 s1 = max(0,-((x1-5).^2 +(x2-5).^2 -82.81));
 s2 = max(0,13-x1);
 s3 = max(0,-x2);
 s4 = max(0, x1-100);
 s5 = max(0, x2-100);
 k = s1^2 + s2^2 + s3^2 + s4^2 + s5^2;
 
 for t = 1
        x1 = x(t,1);
        x2 = x(t,2);
        
        F = penfun(x(t,:),0);
        P(1,t) = penfun(x(t,:),R);
        G(t,:) = mainsm(x(t,:),R);     % to perform unconstrained search for given R
        x(t+1,:) = G(t,:);
        t = t+1;
        P(1,t) = penfun(G(t-1,:),R);
        err = abs(P(1,t));
        R = c*R;
        x1 = x(t,1);
        x2 = x(t,2);
        disp(x);
       
 end

 
 a(:,1) = G(t-1,:)';

 while err>tol2
     G(t,:) = mainsm(x(t,:),R);
     x(t+1,:) = G(t,:);
     t = t+1;
     R = c*R;
     P(1,t) = penfun(G(t-1,:),R);
     err = abs(P(1,t)-P(1,t-1));
       x1 = x(t,1);
       x2 = x(t,2);
       s1 = ((x1-5).^2 +(x2-5).^2 -82.81);
       
     a(:,1) = G(t-1,:)';
     if abs(s1) <tol1
         break;
     elseif (abs(x(t,1)-x(t-1,1))<10e-5 && abs(x(t,2)-x(t-1,2))<10e-5)
         break;
     end
     
 end
 fval = penfun(x(t,:),0);
fprintf('The minimum point is (%4f,%4f) with a function value of %.4f  \n',a,fval)