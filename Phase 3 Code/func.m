%% function used in penalty function method
% enter your objective and constraint function here
function P = penfun(x,R)
     x1 = x(1);
     x2 = x(2);
     s = ((x1-5).^2 +(x2-5).^2 -82.81);
     f = R.*s;       % R * constraint function
     P = (x1-10).^3 + (x2-20).^3 +f;  % Objective function
end