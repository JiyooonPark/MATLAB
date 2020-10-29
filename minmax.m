function [n, m] = minmax(varargin)
if length(varargin) == 2
    n = 2;
    if varargin{1}>varargin{2}
        m = varargin{1};
    else
        m = varargin{2};
    end
else
    n = 3;
    m = varargin{3};
end
