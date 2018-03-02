% open file
nc=netcdf('/media/modeller/data/output.nc');


% plot the planet-sun distances, assuming the sun is the first element:
[r,c,p]=size(nc{'pos'}(:,:,:));
plot(nc{'time'}(:)./(365.25.*86400),sqrt(sum(squeeze(nc{'pos'}(:,:,:).^2),3)));
set(gca,'yscale','log');
xlabel('time (earth years)');
ylabel('planet-sun distance (m)');



% close file
close(nc);
