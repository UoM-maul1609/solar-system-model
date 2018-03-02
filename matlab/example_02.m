% open file
nc=netcdf('/media/modeller/data/output.nc');


% plot the orbits:
[r,c,p]=size(nc{'pos'}(:,:,:));
plot3(nc{'pos'}(:,:,1),nc{'pos'}(:,:,2),nc{'pos'}(:,:,3))

xlabel('x (m)');ylabel('y (m)');zlabel('z (m)')
grid on;

% close file
close(nc);
