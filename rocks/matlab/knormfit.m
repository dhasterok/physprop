close all;
clear all;

data = load('ray_etal2006.dat');

ind = find(isnan(data(:,1)));
ind = [0;ind];

ns = length(ind)-1;

figure;
k0 = [];
tau = [];
kappa = [];
for i = 1:ns
    subplot(221); hold on;
    T = data(ind(i)+1:ind(i+1)-1,1);
    k = data(ind(i)+1:ind(i+1)-1,2);
    tau = [tau; T];
    kappa = [kappa; k];
    k0 = [k0; k(1)*ones(size(k))];
    
    plot(T,k);
    axis square;
    xlabel('T (K)');
    ylabel('k (W m^{-1} K^{-1})');
    ylim([0.5 2.5]);
    
    %m0 = [1.7; 0.3];
    %if i == 1
    %    plot(T,kTmodel(m0,T),'k-');
    %end
    
    %[m,rms] = invnewton('kTmodel',T,k,m0);
    
    A = [ones(size(T)) 1./T T.^2];
    
    m = inv(A'*A)*A'*k;
    
    M(:,i) = m;
    
    subplot(222); hold on;
    %plot(T,kTmodel(m,T));
    plot(T,A*m);
    axis square;
    xlabel('T (K)');
    ylabel('k (W m^{-1} K^{-1})');
    ylim([0.5 2.5]);
    
    subplot(223); hold on;
    plot(T,(k - A*m)./k(1));
    axis square;
    xlabel('T (K)');
    ylabel('k misfit (W m^{-1} K^{-1})');
end

n = mean(M,2);

n(1)
n(2)
n(3)

A = [ones(size(tau)) 1./tau tau.^2];
whos
ma = inv(A'*A)*A'*(kappa./k0);

for i = 1:ns
    subplot(221); hold on;
    T = data(ind(i)+1:ind(i+1)-1,1);
    k = data(ind(i)+1:ind(i+1)-1,2);

    subplot(224); hold on;
    plot(tau,kappa - k0.*(A*ma),'+');
    xlabel('T (K)');
    ylabel('k (W m^{-1} K^{-1})');
    %ylim([0.5 2.5]);
end

