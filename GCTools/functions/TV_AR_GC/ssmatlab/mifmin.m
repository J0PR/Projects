function [y,fval,exitflag] = mifmin(r, den, x1, x2)

options = optimset('TolX',1e-12,'Display', 'off'); % Turn off Display
% options = optimset('TolX',1e-12,'Display', 'iter');  
[y,fval,exitflag] = fminbnd(@poly, x1, x2, options);

    function y = poly(x) % Compute the rational expression.
     y = polyval(r,x)/polyval(den,x);
    end
end
