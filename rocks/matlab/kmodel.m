function [dm,F] = kmodel(k298,x)

n = [1.48;  % quartz
    -0.27;  % orthoclase
    -0.21;  % plagioclase
     1.54;  % muscovite
     1.54;  % biotite
     0.5;   % amphibole
     0.54;  % clinopyroxene
     0.3;   % orthopyroxene 0.49;  % olivine 1.24;  % spinel
     0.37]; % garnet

kmax = [0.8107;  % quartz
      0.1949;  % orthoclase
      0.0000;  % plagioclase
      0.1166;  % muscovite
      0.1166;  % biotite
      0.1725;  % amphibole
      0.1725;  % clinopyroxene
      0.1725;  % orthopyroxene 0.1725;  % olivine 0.0000;  % spinel
      0.1725]; % garnet

Tr = [558;  % quartz
      864;  % orthoclase
        0;  % plagioclase
      520;  % muscovite
      520;  % biotite
      762;  % amphibole
      762;  % clinopyroxene
      762;  % orthopyroxene 762;  % olivine 0;  % spinel
      762]; % garnet

a = [6.553e-3;  % quartz
     3.782e-3;  % orthoclase 
     0.000e-3;  % plagioclase 
     2.749e-3;  % muscovite 
     2.749e-3;  % biotite 
     3.900e-3;  % amphibole 
     3.900e-3;  % clinopyroxene 
     3.900e-3;  % orthopyroxene 3.900e-3;  % olivine 0.000e-3;  % spinel 
     3.900e-3]; % garnet

N = length(x(1,:)) - 1;
T = x(:,1);
w = x(:,2:N+1);

for i = 1:N
    klat(:,i) = k298(i)*(298./T).^n(i);
    krad(:,i) = kmax(i)*(1 + erf(a(i)*(T - Tr(i))));
end

keff = klat + krad;

dm = forward(T,w,n,keff);
F = frechet(T,w,n,keff);

return


% forward model
function dm = forward(T,w,n,keff)

dm = sum(w.*log(keff),2)./sum(w,2);

return


% Frechet matrix
function F = frechet(T,w,n,keff)

[nrow ncol] = size(keff);

for i = 1:ncol
    F(:,i) = w(:,i).*(298./T).^n(i)./keff(:,i);
end

return
