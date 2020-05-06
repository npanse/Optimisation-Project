function [x ,y] = mains2(n,fcn)
 
%% CONJUGATE GRADIENT METHOD  

    function [x]= cg(x0)
        tic

        % initialize
        x=x0;
        eps=1e-3;
        nmax=1000;

        % Initial Gradient 
        gradf = grad(fcn,x);
        disp(gradf);
        s = -gradf;
        d = dot(gradf,gradf);     % Dot product

        % Loop
        i = 0;    % iteration counter
        j = 0;    % Counter for inner while loop
        
        
        F(1)=fcn(x);
        formatSpec='i: %d    j=%d    f(x): %.6e \n';
        fprintf(formatSpec,i,j,F(1));

        %Iteration loop;

        for i = 1:nmax  
            %s = -grad(fcn,x);
            alpha = unidir(x,s,fcn);
            
            x = x + alpha.*s;
            if norm(gradf)>eps
                for j = 1:5  %termination condition

                    gradf = grad(fcn,x);
                    c = dot(gradf,gradf);
                    beta = c/d;
                    temp = s;
                    s = -gradf + beta.*s;
                    d = c;
                    %linear independency termination condition
                    angle = acosd(dot(s,temp)/(norm(s).*norm(temp)));

                    if angle<1            
                        break                   %breaks while loop
                    else
                        alpha=unidir(x,s,fcn);
                        x=x+alpha.*s;

                        %% Display
                        F(j+1)=fcn(x);
                        formatSpec='i: %d    j=%d    f(x): %.6e \n';
                        fprintf(formatSpec,i,j,F(j+1));

                    end
                    %plot(x,F,'.'); hold on
                    j=j+1;
                end
            end
            s = -grad(fcn,x);
            
            F(j+1)=fcn(x);
            
            if norm(gradf)<eps   %Breaks for loop

                f=fcn(x);      
                formatSpec='Function converged\nf(x): %.6e\n\n';
                fprintf(formatSpec,f);
                if F==fcn(x0)
                disp('***Initial point is a local minimum.***')
                end
                break
            end
            if i==nmax
                disp('Number of function evaluations exceeded.')
            end
        end
        toc
    end


%% INPUT
    a = -2;
    b = +2;
    
    %x0 = a + (b-a).*rand(1,n);
    x0 = [0,1.5];
    disp(x0);
    
    x = cg(x0);
    y = fcn(x);
    
end