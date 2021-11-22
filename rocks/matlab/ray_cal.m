close all;
clear all;

%data = load('ray_etal2006.dat');
%model = load('ray_etal2006.rpt');
data = load('cal.dat');
model = load('cal.rpt');

dind = [0; find(isnan(data(:,1)) == 1)];
mind = [0; find(isnan(model(:,1)) == 1)];

Tax = [200 900];
dax = [0 2];
kax = [1 5];
for i = 2:length(dind)
    dlo = dind(i-1)+1;
    dhi = dind(i)-1;

    mlo = mind(i-1)+1;
    mhi = mind(i)-1;

    T = model(mlo:mhi,2);
    rho = model(mlo:mhi,4);
    cp = model(mlo:mhi,5);
    k = model(mlo:mhi,6);
    kappa = model(mlo:mhi,7);

    Td = data(dlo:dhi,1);
    kappad = data(dlo:dhi,2);

    % diffusivity figure
    figure(1);
    subplot(221);
    plot(Td,kappad,'kx');
    hold on;
    plot(T,kappa*1e6,'r-');
    title(num2str(i-1));

    xlabel('T [K]');
    ylabel('\kappa [mm^2/s]');
    axis([Tax dax]);
    set(gca,'Box','on');
    hold off;

    subplot(222);
    plot(Td,1e-6*kappad.*rho.*cp,'kx');
    hold on;
    plot(T,kappa.*rho.*cp,'r-');

    xlabel('T [K]');
    ylabel('k [W/m/K]');
    axis([Tax kax]);
    set(gca,'Box','on');
    hold off;

    subplot(223);
    plot(Td,rho,'r-');
    xlabel('T [K]');
    ylabel('\rho [kg/m^3]');
    set(gca,'Box','on');

    subplot(224);
    plot(Td,cp,'r-');
    xlabel('T [K]');
    ylabel('C_{P} [J/kg/K]');
    set(gca,'Box','on');

    pause;
end
