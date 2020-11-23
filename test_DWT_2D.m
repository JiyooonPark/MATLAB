
% test_DWT_2D.m

% close all;
% clear;
% clc;

x = double( imread('boat.png') );

figure;
imshow( x, [0, 255] );
title('original image');

% Daubechies wavelet filters
lo = [1+sqrt(3), 3+sqrt(3), 3-sqrt(3), 1-sqrt(3)]/(4*sqrt(2));
hi = lo( end:-1:1 );
hi = hi .* [1, -1, 1, -1];

% DWT 2D
[a, h, v, d] = DWT_2D( x, lo, hi );

% IDWT 2D
x1 = IDWT_2D( a, h, v, d, lo, hi, size(x) );

% scale to 0-255
a = scale_img(a);
h = scale_img(h);
v = scale_img(v);
d = scale_img(d);

M = [a, h; v, d];
figure;
imshow( M, [0,255] );
title('coefficients');

figure;
imshow( x1, [0,255] );
title('reconstructed image');

max(abs(x(:) - x1(:)))



%%
function [a, h, v, d] = DWT_2D( x, lo, hi )

    [m, n] = size( x );
    len_f = length( lo );
    
    n1 = floor( (n + len_f - 1)/2 );
    rc = zeros( m, n1 );
    rd = zeros( m, n1 );
    
    % row DWT 
    for i = 1:m
        [rc(i,:), rd(i,:)] = DWT( x(i,:), lo, hi );
    end
    
    % col DWT
    m1 = floor( (m + len_f - 1)/2 );
    a = zeros( m1, n1 );
    h = zeros( m1, n1 );
    v = zeros( m1, n1 );
    d = zeros( m1, n1 );
    
    for j = 1:n1
        [tmp1, tmp2] = DWT( rc(:,j)', lo, hi );
        a(:,j) = tmp1';
        h(:,j) = tmp2';

        [tmp1, tmp2] = DWT( rd(:,j)', lo, hi );
        v(:,j) = tmp1';
        d(:,j) = tmp2';

    end
    
end

function x1 = IDWT_2D( a, h, v, d, lo, hi, size_x )

    n1 = size( a, 2 );
    rc = zeros( size_x(1), n1 );
    rd = zeros( size_x(1), n1 );
    
    for j = 1:n1
        tmp = IDWT( a(:,j)', h(:,j)', lo, hi, size_x(1) );
        rc(:,j) = tmp';
        
        tmp = IDWT( v(:,j)', d(:,j)', lo, hi, size_x(1) );
        rd(:,j) = tmp';
        
    end

    x1 = zeros( size_x );
    for i = 1:size_x(1)
        x1(i,:) = IDWT( rc(i,:), rd(i,:), lo, hi, size_x(2) );
    end
    
end



function [c, d] = DWT(x, lo, hi)

    len_lo = length(lo);
    len_ext = len_lo - 1;
    len_x = length(x);
    x_ext = [ x(len_ext:-1:1), x, x(len_x:-1:(len_x-len_ext+1)) ];
    
    lo_re = lo( end:-1:1 );
    hi_re = hi( end:-1:1 );
    c = conv(x_ext, lo_re, 'valid');
    d = conv(x_ext, hi_re, 'valid');

    c = c( 2:2:end );
    d = d( 2:2:end );

end

function x = IDWT(c, d, lo, hi, len_x)
% function [x, xc, xd] = IDWT(c, d, lo, hi, len_x)


    len_c = length(c);
    
    if mod( len_x, 2 ) == 0
        extra = 1;
    else
        extra = 0;
    end
    
    c_up = zeros(1, 2*len_c+extra);
    c_up(2:2:end) = c;
    
    d_up = zeros(1, 2*len_c+extra);
    d_up(2:2:end) = d;

    xc = conv(c_up, lo, 'valid');
    xd = conv(d_up, hi, 'valid');

    x = xc + xd;    
    
end


