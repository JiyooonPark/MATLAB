
% test_DWT_2D_comp.m

close all;
clear;
clc;

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

% Histogram
% figure;
% histogram( abs(a(:)) ); title('histogram of |a|');
% figure;
% histogram( abs(h(:)) ); title('histogram of |h|');
% figure;
% histogram( abs(v(:)) ); title('histogram of |v|');
% figure;
% histogram( abs(d(:)) ); title('histogram of |d|');

% Mean
fprintf( 'mean of |h| = %f\n', mean( abs(h(:)) ) );
fprintf( 'mean of |v| = %f\n', mean( abs(v(:)) ) );
fprintf( 'mean of |d| = %f\n', mean( abs(d(:)) ) );

% Small coefficients
tau = 20; % 91퍼센트를 감소시켜도 우리 눈에는 티가 잘 안난다
idx_h = abs(h) < tau; % 0 or 1 goes into idx_h
    % sum(idx_h(:));
    % sum(idx_h(:))/257*257*100 = 23.7203의 공간을 세이브 할 수 있다
idx_v = abs(v) < tau;
idx_d = abs(d) < tau;

fprintf( '# of |h| < %f = %d\n', tau, sum( idx_h(:) ) );
fprintf( '# of |v| < %f = %d\n', tau, sum( idx_v(:) ) );
fprintf( '# of |d| < %f = %d\n', tau, sum( idx_d(:) ) );

% Thresholding
% 해당 위치에 0을 넣는다 
h( idx_h ) = 0;
v( idx_v ) = 0;
d( idx_d ) = 0;

% IDWT 2D
x1 = IDWT_2D( a, h, v, d, lo, hi, size(x) );

a_s = scale_img(a);
h_s = scale_img(h);
v_s = scale_img(v);
d_s = scale_img(d);

M = [a_s, h_s; v_s, d_s];
figure;
imshow( M, [0,255] );
title('coefficients');

figure;
imshow( x1, [0,255] );
title('reconstructed image');

fprintf('max error = %f\n', max(abs(x(:) - x1(:))) );




%%
function [a, h, v, d] = DWT_2D( x, lo, hi )

    [m, n] = size( x );
    len_f = length( lo );
    
    n1 = floor( (n + len_f - 1)/2 );
    rc = zeros( m, n1 );
    rd = zeros( m, n1 );
    
    for i = 1:m
        [rc(i,:), rd(i,:)] = DWT( x(i,:), lo, hi );
    end

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

