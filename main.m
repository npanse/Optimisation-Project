function [a b z z1]= main(x0)

%% BOUNDINGPHASE METHOD

function [a b k]=BoundingPhase(x0,del,funcName)

    fx1=feval(funcName,(x0-del));
    fx0=feval(funcName,x0);
    fx3=feval(funcName,(x0+del));
  

    if(fx1 > fx3)
        del=del;
    else
        del=-del;
    end
    n=2;
    k=2;
    x(1)=x0;
    x(2)=x(1) + 2*del;
    fx1=feval(funcName,x(2));
    
    while(fx1<fx0)
            fx0=fx1;
            x(k+1)=x(k)+2^k*del;
            %disp (x(k+1))
            n=k;
            k=k+1;
            fx1=feval(funcName, x(k));
    end

    a=x(n-1);
    b=x(k);
end

%% SECANT METHOD


function [z z1]= secant(a,b,funcName)

    fa=feval(funcName,a);
    fb=feval(funcName,b);

    e = 10^(-4);

    %%DERIVATIVE

    function [f]=deriv(x)

    del = 0.001;
    fx1=feval(funcName,(x+del));
    fx2=feval(funcName,(x-del));
    f=(fx1-fx2)/(2*del);

    end

    f1a = deriv(a);
    f1b = deriv(b);
    
    % compute derivative of the function at the end points
    % according to the sign of the function proceed or finish the program

    if f1a==0
        z = a;
        z1 = f1a;

    elseif f1b==0
        z = b;
        z1 = f1b;

    elseif f1a<0 && f1b>0
        x1=a;
        x2=b;
        z=x2-(deriv(x2))/((deriv(x2)-deriv(x1))/(x2-x1));
    else
        fprintf('Please try other values of a and b to locate the minima\n')
        return
    end

    % Iterate until the termination condition is reached.
     
    while(deriv(z) > e)
        z=x2-(deriv(x2))/((deriv(x2)-deriv(x1))/(x2-x1));
       % disp(z);
        if deriv(z)>0
            x1=a;
            x2=z;
        elseif deriv(z)<0
            x1=z;
            x2=b;
        end
    end
    z1=feval(funcName,z);
end


%% MAIN

    % Taking a function
    [f] = "funcitons";
    
    % boundingphase method with initial guess x0 and del = 0.01
    [a b k]=BoundingPhase(x0,0.001,f);   
    
    % secant method with output of boundingphase method as input
    if a<b
        [z z1]= secant(a,b,f);
    else
        [z z1]= secant(b,a,f);
    end
    % z = minima and z1 = minimum value of function
end    
