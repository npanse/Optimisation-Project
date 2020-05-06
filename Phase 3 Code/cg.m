function [x] = mainsm(x0,R)
 
%% CONJUGATE GRADIENT METHOD  

    function [x]= cg(x0)
        tic

        % initialize
        x=x0;
        eps=1e-3;
        nmax=100;

        % Initial Gradient 
        gradf = grad(@penfun,x,R);
        %
        disp(gradf);
        
        d = dot(gradf,gradf);     % Dot product

        % Loop
        i = 0;    % iteration counter
        j = 0;    % Counter for inner loop
        u = 1;    % no. of unidirectional searches
        
        F(u)=penfun(x,R);
        formatSpec='i: %d    j=%d    f(x): %.6e \n';
        fprintf(formatSpec,i,j,F(u));

        %Iteration loop;

        for i = 1:nmax  
            s = -grad(@penfun,x,R);
            alpha = unidir(x,s,R);
            x0 = x;
            x = x + alpha.*s;
            
            if i==nmax
                disp('Number of function evaluations exceeded.')
            end
             
            for j = 1:5  %termination condition
                                    
                    gradf = grad(@penfun,x,R);
                    c = dot(gradf,gradf);
                    beta = c/d;
                    temp = s;
                    s = -gradf + beta.*s;
                    d = c;
                   
                    %linear dependency termination condition
                    angle = acosd(dot(s,temp)/(norm(s).*norm(temp)));

                    if angle<1            
                        break                   %breaks inner for loop
                    else
                        alpha=unidir(x,s,R);
                        x0 = x;
                        x = x + alpha.*s;
                        u = u + 1;
                        gradf = grad(@penfun,x,R);
                        if norm(gradf)<eps   
                            break           %Breaks for loop
                        end
                        if norm(x-x0)/norm(x0) < eps
                            break           %Breaks for loop 
                        end
                        %% Display
                        F(u)=penfun(x,R);
                        formatSpec='i: %d    j=%d    f(x): %.6e \n';
                        fprintf(formatSpec,i,j,F(u));

                    end
                    %plot(i,F,'.'); hold on
                    
                    j = j+1;   
            end
            
            if norm(grad(@penfun,x,R))<eps
                f=penfun(x,R);      
                formatSpec='Function converged\nf(x): %.6e\n\n';
                fprintf(formatSpec,f);
                break           %Breaks for loop
            end
            
            if norm(x-x0)/norm(x0) < eps
                f=penfun(x0,R);      
                formatSpec='Function converged\nf(x): %.6e\n\n';
                fprintf(formatSpec,f);
                break           %Breaks for loop 
            end
            
        end
        l = length(F);
        %plot(1:l,F);
        disp('Number of unidirectional searches = ');
        disp(u);
        
        toc
    end


%% INPUT
    a = -2;
    b = +2;
    
    %x0 = a + (b-a).*rand(1,n);
    %x0 = [1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5];
    %x_ = [0,0,0,0,0];
    disp(x0);
    x = cg(x0);
    y = penfun(x,R);
    
end