% calculate argument of perihelion
create_plot=false;

calculate_orbital_props=true;

if calculate_orbital_props
    object=7:18;
    %jj=1;
    %filenames={'/Volumes/My Passport/Holly Grey/Reference model/reference model.nc',...
    %    '/Volumes/My Passport/Holly Grey/Model 12/output.nc','/Volumes/My Passport/Holly Grey/Model 13/output.nc'};

    %filenames={'/Volumes/My Passport/Holly Grey/Model 25/output.nc'};
    col={'r','g','b'};

    nc=netcdf(filenames{jj});
    times=nc{'time'}(:);


    for i=object
        disp(['iteration ',num2str(i-min(object)+1)])
        x=nc{'pos'}(:,i,1);
        y=nc{'pos'}(:,i,2);
        z=nc{'pos'}(:,i,3);
        dist=sqrt(x.^2+y.^2+z.^2);

        vx=nc{'vel'}(:,i,1);
        vy=nc{'vel'}(:,i,2);
        vz=nc{'vel'}(:,i,3);

        % rxv
        omg=cross([x,y,z],[vx,vy,vz]);



        % perigee:
        [perigee,perigee_locs] = findpeaks(-dist,...
            'MinPeakProminence',0.1.*max(dist));
        perigee=-perigee;
        % apegee:
        [apegee,apegee_locs] = findpeaks(dist,...
            'MinPeakProminence',0.1.*max(dist));

        % ascending node
        [ascending_node,ascending_node_locs] = findpeaks([diff(z);0],...
            'MinPeakProminence',0.1.*max(diff(z)));

        len1=min([length(perigee_locs),length(ascending_node_locs),length(apegee_locs)]);

        perigee=perigee(1:len1);
        apegee=apegee(1:len1);

        sma=(perigee+apegee)./2;
        ascending_node=ascending_node(1:len1);
        perigee_locs=perigee_locs(1:len1);
        ascending_node_locs=ascending_node_locs(1:len1);
        apegee_locs=apegee_locs(1:len1);


        x_peri=x(perigee_locs);
        y_peri=y(perigee_locs);
        z_peri=z(perigee_locs);

        x_an=x(ascending_node_locs);
        y_an=y(ascending_node_locs);
        z_an=z(ascending_node_locs);

        % normal of the plane
        normal=cross([x_peri,y_peri,z_peri],[x_an,y_an,z_an]);
        %triple product
        one1=cross([x_an,y_an,z_an],normal);
        dp=x_peri.*one1(:,1)+y_peri.*one1(:,2)+z_peri.*one1(:,3);


        % calculate the longitude of ascending node:
    %     angle_peri=atan2(y_peri,x_peri);
    %     angle_an=atan2(y_an,x_an);
    %     long_an=angle_an-angle_peri;
    %     ind=find(long_an<0);
    %     long_an(ind)=2.*pi+long_an(ind);
        % https://en.wikipedia.org/wiki/Longitude_of_the_ascending_node
        k=zeros(len1,3);
        k(:,3)=1;
        h=omg(ascending_node_locs,:);
        n=cross(k,h);
        ny=n(:,2);
        long_an=acos(n(:,1)./sqrt(1e-15+sum(n.^2,2)));
        ind=find(ny<0);
        long_an(ind)=2.*pi-acos(n(ind,1)./sqrt(1e-15+sum(n(ind,:).^2,2)));

        orbit(i-min(object)+1).long_an=long_an.*180./pi;

        % calculate the argument of perihelion:
        %https://en.wikipedia.org/wiki/Argument_of_periapsis
        dot=x_peri.*x_an+y_peri.*y_an+z_peri.*z_an;
        mag_peri=sqrt(x_peri.^2+y_peri.^2+z_peri.^2);
        mag_an=sqrt(x_an.^2+y_an.^2+z_an.^2);
        arg_peri=acos(dot./(mag_peri.*mag_an+1e-15));
        orbit(i-min(object)+1).arg_peri=arg_peri;
        orbit(i-min(object)+1).arg_peri_deg=arg_peri.*180./pi;




    %     omg_peri=omg(perigee_locs,3);
    %     ind=find(omg_peri<0);
    %     orbit(i-min(object)+1).arg_peri_deg(ind)=360-arg_peri(ind).*180./pi;

        ind=find(normal(:,3)<0);
        orbit(i-min(object)+1).arg_peri_deg(ind)=360- ...
            orbit(i-min(object)+1).arg_peri_deg(ind);

        %     omg_peri=omg(perigee_locs,3);
    %     ind=find(omg_peri<0);
    %     orbit(i-min(object)+1).arg_peri_deg(ind)=360-...
    %         orbit(i-min(object)+1).arg_peri_deg(ind);

        %     orbit(i-min(object)+1).arg_peri_deg(ind)=-arg_peri(ind).*180./pi;
    %     orbit(i-min(object)+1).arg_peri_deg=orbit(i-min(object)+1).arg_peri_deg.* ...
    %         -sign(normal(:,3));


    orbit(i-min(object)+1).time_peri=times(perigee_locs); 

    %     orbit(i-min(object)+1).omg_peri=omg(perigee_locs);
        orbit(i-min(object)+1).x_peri=x_peri;
        orbit(i-min(object)+1).y_peri=y_peri;
        orbit(i-min(object)+1).z_peri=z_peri;
        orbit(i-min(object)+1).dp=dp;
        orbit(i-min(object)+1).sma=sma;
        orbit(i-min(object)+1).perigee=perigee(1:len1);
        orbit(i-min(object)+1).apegee=apegee(1:len1);


    end
    close(nc);

end



% plot(orbit(i).sma,orbit(i).arg_peri_deg,'.')





% % do a 2-d histogram of all argument of periapsi
% tbin=linspace(0,max(times),20);
% tbin2=[tbin,tbin(end)-tbin(end-1)];
% abin=linspace(0,360,15);
% abin2=[abin,abin(end)-abin(end-1)];
% bins=zeros(20,15);
% 
% meanbin=zeros(20,1);
% prctilebin=zeros(20,2);
% stdbin=zeros(20,1);
% for j=1:length(tbin)
%     dat1=[];
%     for k=1:length(abin)
%         
%         for i=1:length(object)
%             ind=find(orbit(i).time_peri>=tbin2(j) & ...
%                 orbit(i).time_peri<tbin2(j+1) &...
%                 orbit(i).arg_peri_deg>=abin2(k) & ...
%                 orbit(i).arg_peri_deg<abin2(k+1));
%             bins(j,k)=bins(j,k)+length(ind);
%             
%             dat1=[dat1;orbit(i).arg_peri_deg(ind)];
%         end
%         
%     end
%     meanbin(j)=nanmean(dat1);
%     prctilebin(j,1)=prctile(dat1,25);
%     prctilebin(j,2)=prctile(dat1,75);
%     stdbin(j)=nanstd(dat1);
% end

% pcolor(tbin,abin,bins');shading flat




% do a 2-d histogram of all argument of periapsi
au=1.496e+11;
smas=cat(1,orbit.sma)./au;
xbin=150:100:1150;
xbin2=[xbin,xbin(end)-xbin(end-1)];
abin=linspace(0,360,15);
abin2=[abin,abin(end)-abin(end-1)];
bins=zeros(length(xbin),15);

meanbin=zeros(length(xbin),1);
prctilebin=zeros(length(xbin),2);
stdbin=zeros(length(xbin),1);
n_sample=zeros(length(xbin),1);
for j=1:length(xbin)
    dat1=[];
    for k=1:length(abin)
        
        for i=1:length(object)
            ind=find(orbit(i).sma./au>=xbin2(j) & ...
                orbit(i).sma./au<xbin2(j+1) &...
                orbit(i).arg_peri_deg>=abin2(k) & ...
                orbit(i).arg_peri_deg<abin2(k+1));
            bins(j,k)=bins(j,k)+length(ind);
            
            dat1=[dat1;orbit(i).arg_peri_deg(ind)];
        end
        
    end
    meanbin(j)=nanmean(dat1);
    prctilebin(j,1)=prctile(dat1,25);
    prctilebin(j,2)=prctile(dat1,75);
    stdbin(j)=nanstd(dat1);
    n_sample(j)=length(dat1);
end
if create_plot
    errorbar(xbin,meanbin,stdbin,col{jj})
    hold on;
    ylim([0 360]);
end


S_control=stdbin(1);
S_compare=stdbin(2:end);
n_control=n_sample(1);
n_compare=n_sample(2:end);

F=S_control.^2./S_compare.^2;

ftheory=finv(0.95,n_control,n_compare);
ftheory1=finv(0.95,length(object).*ones(length(object),1),length(object).*ones(length(object),1));

fraction_lt_500=F(1:3)./ftheory1(1:3);
fraction_of_clustering=length(find(fraction_lt_500>1))./length(fraction_lt_500);


%
% subplot(211)
% plot(xbin,meanbin,col{j})
% hold on;
% subplot(212)
% plot(xbin,stdbin,col{j})
% hold on;
% pcolor(xbin,abin,bins');shading flat