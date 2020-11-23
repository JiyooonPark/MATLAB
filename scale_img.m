
function y = scale_img( x )

    a = 0;
    b = 255;
    xmin = min(x(:));
    xmax = max(x(:));

    y = (b - a)/(xmax - xmin)*(x - xmin) + a;

end


