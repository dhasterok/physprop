close all;
clear all;

%k298 = [8.79; 1.79; 2.2; 1.71; 1.39; 2.65; 4.25; 3.37; 4.75; 11.94; 4.97]; 
k298 = [8.79;  % quartz
        1.79;  % orthoclase
        2.20;  % plagioclase
        1.71;  % muscovite
        1.39;  % biotite
        3.80;  % amphibole
        4.25;  % clinopyroxene
        3.37;  % orthopyroxene
        4.97]; % garnet

mineral =  {'quartz'; 'orthoclase'; 'plagioclase'; 'muscovite';
    'biotite'; 'amphibole'; 'clinopyroxene'; 'orthopyroxene';
    'olivine'; 'spinel'; 'garnet'};

data = load('cal.csv');
mid = [ 0 NaN;  % quartz
        2 NaN;  % orthoclase
        4   5;  % plagioclase
        7 NaN;  % muscovite
        8 NaN;  % biotite
       13 NaN;  % amphibole
       17  18;  % clinopyroxene
       21  22;  % orthopyroxene % 23  24; 25 NaN;
       29 NaN]; % garnet

N = length(data(1,:)) - 1;
for i = 1:length(mid(:,1))
    if isnan(mid(i,2))
        ind = find(mid(i,1) == data(:,1));
        w(i,:) = data(ind,2:N+1);
    else
        ind = find(mid(i,1) == data(:,1) | mid(i,2) == data(:,1));
        w(i,:) = sum(data(ind,2:N+1),1);
    end
end

for i = 1:N
    w(:,i) = w(:,i)/sum(w(:,i));
end

w = w';

rpt = load('cal.rpt');

nans = find(isnan(rpt(:,1)));
ind = find(~isnan(rpt(:,1)));

cal = load('cal.dat');
cal(:,2) = 1e-6*cal(:,2).*rpt(:,4).*rpt(:,5);

T = cal(ind,1);
kobs = cal(ind,2);

for i = 1:length(nans)
    if i == 1
        top = 1;
        Top = 1;
    else
        top = nans(i-1)-i+2;
        Top = nans(i-1)+1;
    end
    bot = nans(i) - i;
    Bot = nans(i) - 1;

    [top bot Top Bot]
    for j = top:bot
        W(j,:) = w(i,:);
    end
    v(i,:) = [top bot];
end
ind = find(T <= 800);
[m,rms] = invnewton('kmodel',[T(ind) W(ind,:)],log(kobs(ind)),k298);
km = exp(kmodel(m,[T W]));


km298 = exp(m);

k0 = exp(kmodel(k298,[T W]));
for i = 1:length(v(:,1))
    top = v(i,1);
    bot = v(i,2);

    plot(T(top:bot),kobs(top:bot),'k.');
    hold on;
    plot(T(top:bot),k0(top:bot),'b-');
    plot(T(top:bot),km(top:bot),'r-');
    pause;
    hold off;
end
