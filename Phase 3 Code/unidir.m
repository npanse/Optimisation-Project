function [beta]= unidir(x_,s,R)


%% Function in terms of alpha

function [f_] = func_alpha(alpha,R)
    
    x_0 = x_ + alpha.*s;
    f_ = penfun(x_0,R);
end

%% BOUNDINGPHASE METHOD

function [a,b]=BoundingPhase(x0)

    del = 0.01;
    
    fx1 = f(x0-del,R);
    fx0 = f(x0,R);
    fx3 = f(x0+del,R);
  
    if(fx1 > fx3)
        del=del;
    else
        del=-del;
    end
    n=2;
    k=2;
    x(1)=x0;
    x(2)=x(1) + 2*del;
    fx1=f(x(2),R);
    
    while(fx1<fx0)
            fx0=fx1;
            x(k+1)=x(k)+2^k*del;
            %disp (x(k))
            n=k;
            k=k+1;
            fx1=f(x(k),R);
    end
    %disp("-----") ;
    a=x(n-1);
    b=x(k);
end

%% SECANT METHOD


function [z]= secant(a,b)

    a=0;
    b=5;
    fa=f(a,R);
    fb=f(b,R);
    nmax = 500;
    e = 10^(-4);

    %%DERIVATIVE

    f1a = grad(f,a,R);
    f1b = grad(f,b,R);
    
    % compute derivative of the function at the end points
    % according to the sign of the function proceed or finish the program

    if f1a==0
        z = a;
        z1 = fa;

    elseif f1b==0
        z = b;
        z1 = fb;

    elseif f1a<0 && f1b>0
        x1=a;
        x2=b;
        z = (x1+x2)/2; 
        %z = x2 - (grad(f,x2)./(grad(f,x2)-grad(f,x1))./(x2-x1));
    else
        %fprintf('Please try other values of a and b to locate the minima\n')
        z = a+b/2;
        return
    end

    % Iterate until the termination condition is reached.
    i=0;
    while(grad(f,z,R) > e)
        z = (x1+x2)/2;
        %z = x2 - (grad(f,x2))./((grad(f,x2)-grad(f,x1))./(x2-x1));
        %disp(grad(f,z));
        i = i+1;
        if grad(f,z,R)>0
            x1=a;
            x2=z;
        elseif grad(f,z,R)<0
            x1=z;
            x2=b;
        end
    end
    %disp("-----");
    z1 = f(z,R);
end



%% MAIN
    
    d = length(x_);
    alpha = 5*rand();
    %alpha = randi([-10 10],1);
    f = @func_alpha;
        
    % boundingphase method with initial guess x0 and del = 0.01
    [a,b]=BoundingPhase(alpha);   
    
    % secant method with output of boundingphase method as input
    if a<b
        [beta]= secant(a,b);
    else
        [beta]= secant(b,a);
    end
    % z = minima and z1 = minimum value of function
end    
