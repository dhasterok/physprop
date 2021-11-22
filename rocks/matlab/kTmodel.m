function [dm,F] = kTmodel(m,x)

dm = forward(m,x);
F = frechet(m,x);

return

function dm = forward(m,x)
    dm = m(1)*(298./x).^m(2);
return

function F = frechet(m,x)
    F = [(298./x).^m(2) ...
         m(1)*(298./x).^m(2).*log(298./x)];
return