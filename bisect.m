function [root,fx,ea,iter]=bisect(func,xl,xu,es,maxit,varargin)

test = func(xl,varargin{:})*func(xu,varargin{:});
if test>0,error('no sign change'),end
iter = 0; xr = xl; ea = 100;
while (1)
    xrold = xr;
    xr = (xl + xu)/2
    iter = iter + 1;
    if xr ~= 0, ea = abs((xr - xrold)/xr) * 100, end
    test = func(xl,varargin{:})*func(xr,varargin{:});
    if test < 0
        xu = xr;
    elseif test > 0
        xl = xr;
    else
        ea = 0;
    end
    if ea <= es || iter >= maxit,break,end
end
root = xr; 
fx = func(xr, varargin{:});