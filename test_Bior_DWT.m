
% test_Bior_DWT.m
% 
% close all;
% clear;
% clc;

% original signal
K = 500;
%등간격으로 100개 그 범위 내에 넣기 
t = linspace(0, 2*pi, K);
% create signals ( different on each run)
x = sin(t) + 0.1*randn(1, K);
% figure;
% plot(t,x);

%% modify 
% lo_d = h tilde

% biothogonal wavelet filters
[lo_D, hi_D, lo_R, hi_R] = wfilters('bior4.4');
lo_D = lo_D( end :-1:1)
hi_D = hi_D( end :-1:1)
%%


% DWT, decomposition
% finding out c, d is our main goal

[c, d] = DWT(x, lo_D, hi_D);
%%
% IDWT, reconstruction
x1 = IDWT(c, d, lo_R, hi_R, length(x));

% reconstruction error
fprintf('max error = %10.4e\n', max( abs( x - x1 ) ));

% plot
figure;
subplot(4, 1, 1); % 4 rows 1 col 1st
plot(x);
title('original signal');

subplot(4, 1, 2); % 4 rows 1 col 1st
plot(d);
title('detail coefficients');

subplot(4, 1, 3);
plot(c);
title('core coefficients');

subplot(4, 1, 4);
plot(x1);
title('reconstructed signal');

%%

function [c, d] = DWT(x, lo, hi)

    len_lo = length(lo);
    len_ext = len_lo - 1;
    len_x = length(x);
    % try it yourself 
    x_ext = [ x(len_ext:-1:1), x, x(len_x:-1:(len_x-len_ext+1)) ];
    
    lo_re = lo( end:-1:1 ); % same as lo(-)
    hi_re = hi( end:-1:1 );
    c = conv(x_ext, lo_re, 'valid');
    d = conv(x_ext, hi_re, 'valid');
        
    %downsampling
    c = c( 2:2:end );
    d = d( 2:2:end );

end

function x = IDWT(c, d, lo, hi, len_x)

    len_c = length(c);
    
    % need 0 at the end?
    
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









