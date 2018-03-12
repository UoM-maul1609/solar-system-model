object=13;
ms=1; %object-1;
plot_type=1; % 1 or 2
filenames={'/tmp/reference model.nc',...
    '/tmp/model12output.nc',...
    '/tmp/model13output.nc'};

for i=1:length(filenames)
    disp(['iteration: ',num2str(i)]);
    [dat(i)]=fourier_transform_basic(filenames{i},object);
end

t=1./dat(1).f(2:end);

switch plot_type
    case 1
    ycoord=1:length(filenames);
    dat2=[];
    for i=1:length(filenames)
        dat2=[dat2;dat(i).fft(ms,2:end)];
    end
    pcolor(t,ycoord,real(log10(dat2)));shading flat;set(gca,'xscale','log');
    
    case 2
    % or find the peak and plot

    arr=zeros(length(filenames),1)
    for i=1:length(filenames)
        ind=find(dat(i).fft(2:end)==max(dat(i).fft(2:end)));
        arr(i)=t(ind(1));
    end
    plot(arr);
    
    otherwise
        disp('error plot_type');
end