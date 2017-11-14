% open file
nc=netcdf('/tmp/output.nc');
dims=3;

% Fourier analysis of the sun-planet distance for all the planets
figure
try1=[0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1];
i=1;
% indices=1:length(nc{'time'}(:));
% indices=1:1550;

% to reference by time:
times=nc{'time'}(:);
sim_year=times./365.25./86400;
indices=find(sim_year>10 & sim_year<=150);
indices=1:length(sim_year);

len1=length(indices);
NFFT=len1;

while NFFT>=len1
    L=fix(try1(i).*length(sim_year(indices)));
    NFFT=2^nextpow2(L);
    i=i+1;
end
Fs=1./abs((times(2)-times(1))./365.25./86400);
f = Fs/2*linspace(0,1,NFFT/2+1); 

% bodies={'Sun','Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune','Pluto'};
bodies={'Sun','Jupiter','Saturn','Uranus','Neptune','Pluto','2010 GB-174',...
    '2004 VN-112', '2000 CR-105', '2005 RH-52', '2003 HB-57', '2007 TG-422', ...
     '2002 GB-32', '2007 VJ-305', '2010 VZ-98', '2001 FP-185', '2012 VP-113','Sedna'};

h = waitbar(0,'Please wait...');
for i=2:length(bodies)
    %subplot(5,4,i-1);
    figure
    x=nc{'pos'}(:,i,1);
    y=nc{'pos'}(:,i,2);
    z=nc{'pos'}(:,i,3);
    dist=sqrt(x.^2+y.^2+z.^2);
    Y=fft(dist  ,NFFT)/L;

    plot(1./f,(2*abs(Y(1:NFFT/2+1))),'r')
    set(gca,'yscale','log','xscale','log')
    xlabel('Period (years)');ylabel('Power');
    text(0.1,0.9,bodies{i},'fontsize',15,'units','normalized');
    waitbar((i-1)/(length(bodies)-1),h)
end
close(h);

% close file
close(nc);
