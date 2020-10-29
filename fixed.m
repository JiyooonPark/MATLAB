function [root,ea,iter]=fixed(g, x, es, maxit, varargin)
iter = 0; ea=100;
while (1)
    xold = x;
    x = g(x, varargin{:})
    iter = iter + 1;
    if x ~= 0, ea = abs((x - xold)/x) * 100, end
    if ea <= es || iter >= maxit, break, end
end
root = x;
end
% 여기는 x = g(x)형식임 
 