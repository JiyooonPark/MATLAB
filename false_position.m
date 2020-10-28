function [ r ] = false_position( f, a, b, N, eps_step, eps_abs )
    % Check that that neither end-point is a root
    % and if f(a) and f(b) have the same sign, throw an exception.

    if ( f(a) == 0 )
	r = a;
	return;
    elseif ( f(b) == 0 )
	r = b;
	return;
    elseif ( f(a) * f(b) > 0 )
        error( 'f(a) and f(b) do not have opposite signs' );
    end

    % We will iterate N times and if a root was not
    % found after N iterations, an exception will be thrown.

    c_old = Inf;

    for k = 1:N
        % Find the false position
        c = (a*f(b) + b*f(a))/(f(b) - f(a));

        % Check if we found a root or whether or not
        % we should continue with:
        %          [a, c] if f(a) and f(c) have opposite signs, or
        %          [c, b] if f(c) and f(b) have opposite signs.

        if ( f(c) == 0 )
            r = c;
            return;
        elseif ( f(c)*f(a) < 0 )
            b = c;
        else
            a = c;
        end

        % If |b - a| < eps_step, check whether or not
        %       |f(a)| < |f(b)| and |f(a)| < eps_abs and return 'a', or
        %       |f(b)| < eps_abs and return 'b'.

        if ( abs( c - c_old ) < eps_step )
            if ( abs( f(a) ) < abs( f(b) ) && abs( f(a) ) < eps_abs )
                r = a;
                return;
            elseif ( abs( f(b) ) < eps_abs )
                r = b;
                return;
            end
        end

	c_old = c;
    end

    error( 'the method did not converge' );
end