% open file
nc=netcdf('/tmp/output.nc');

% Fourier analysis of the sun-planet distance for all the planets
figure
len1=length(nc{'time'}(:));
NFFT=len1;
try1=[0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1];
i=1;
while NFFT>=len1
    L=fix(try1(i).*length(nc{'time'}(:)));
    NFFT=2^nextpow2(L);
    i=i+1;
end
Fs=1./abs((nc{'time'}(2)-nc{'time'}(1))./365.25./86400);
f = Fs/2*linspace(0,1,NFFT/2+1); 

bodies={'Sun','Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune','Pluto'};
for i=2:length(bodies)
    subplot(3,3,i-1);
    
    Y=fft(sqrt(sum(squeeze(nc{'pos'}(:,i,:)).^2,2))  ,NFFT)/L;

    plot(1./f,(2*abs(Y(1:NFFT/2+1))),'r')
    set(gca,'yscale','log','xscale','log')
    xlabel('Period (years)');ylabel('Power');
    text(0.1,0.9,bodies{i},'fontsize',15,'units','normalized');
end

% close file
close(nc);
