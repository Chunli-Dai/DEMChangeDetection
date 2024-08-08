function  [Co] = diffdems_sub(dX4Sg, idregion, XYbg, f, fdir, rang0)
%Modified from ../codetest/diffdems.m
%study area boundary [xmin xmax ymin ymax]
%rang0=[331000  351000  3399000 3416000]; %pick your study region ;himalayas
dir='./';
Co=[];

%dX4Sg, translational offsets to apply to each DEM;
%idregion, index of DEMs within the target region;
%f,	list of filenames;
%fdir, directory of the filenames;

resr=2;%400; %Pick your resolution
tx=rang0(1):resr:rang0(2);ty=rang0(4):-resr:rang0(3);
xout=tx;yout=ty;

%get profiles
flagprof=0; %1 plot profiles; 0 not plot profiles;
flagdem=1; % 1 write coregistered dems; 2 write sequential DEM difference;
if flagprof==1 %
BBfile='BB1.dat';
out=load(BBfile); %out=[loni(:),lati(:),xi(:),yi(:)];
xi=out(:,3);yi=out(:,4);

%densify the line or get distances
clx=xi;cly=yi;dc=2.;
[clx2,cly2,S2]=interpcl(clx,cly,dc);
dsi=S2;xi=clx2;yi=cly2;
end

t=zeros(length(idregion),1); 

%sort idregion based on time
idblock=1:length(idregion);
id=idregion;
    nid=length(id);
    t=zeros(nid,1); 
    for j=1:nid %0%length(id)
        filename=f{id(j)};
%       ymd=filename2ymd(filename);
        infile= [fdir{id(j)},'/',f{id(j)}];
        [ymd,~]=strip2date(infile,3);
        %t(j)=datenum(ymd,'yyyymmdd');
        t(j)=datenum(ymd,'yyyymmddHHMMSS');
    end
    [~,idsort]=sort(t);idblock=idblock(idsort);

%loading DEM and apply offsets.
Legendstr=cell(length(idregion),1); idd=[];
%figure
%hold all
flagdem=1; % 1 write coregistered dems; 2 write sequential DEM difference;
if flagdem==1
   lenn=length(idregion);
elseif flagdem==2;
   lenn=length(idregion)-1;
end
for j=1:lenn  %length(idregion)-1 %[3,5,24,26] % 2012 Kamchatka Volcano 
        i=idregion(idblock(j));

        demdir=fdir{i};
        infile= [demdir,'/',f{i}];
        fprintf(['\n Loading DEM files in order (',num2str(j),'/',num2str(length(idregion)),'): ',infile])
        [ymd,~]=strip2date(infile,3);
        t(j)=datenum(ymd,'yyyymmddHHMMSS');
	XYbi=XYbg{i};

        if strcmp(infile(end-7:end),'meta.txt')
          str1='meta.txt';
        elseif strcmp(infile(end-7:end),'_mdf.txt')
          str1='mdf.txt';
        end
        infile1= strrep(infile,str1,'dem.tif');

        %data=readGeotiff(infile1,'map_subset',rang0);
	[data,~]=readdem(XYbi,infile,rang0,2);

        %Interpolation requires at least two sample points in each dimension.
        [m1,n1]=size(data.z);
        if min([m1,n1])<=1;idd=[idd,j];continue;end

        %refers to ~/arcticdemapp/river/rivergithub2/rivermainserial/riverprofsub.m
        p=dX4Sg(idblock(j),:); %pzxy

        fprintf(['\n Loading its translational offsets pzxy (',num2str(j),'/',num2str(length(idregion)),'): ',num2str(p)])
        data.z(data.z == -9999) = NaN; % %convert nodata values to nan for interpolation
        z2inp=interp2(data.x-p(2),data.y-p(3),double(data.z)-p(1),xout,yout','*linear',nan);
	clear dem1;dem1.x=xout;dem1.y=yout;dem1.z=z2inp;
	z2inp(isnan(z2inp))=-9999; %return to -9999; 

	%This should be what you need; the coregistered DEMs.
	demo.x=xout;demo.y=yout;demo.z=z2inp;

	%get profiles
	if flagprof==1 %
           tz=interp2(data.x-p(2),data.y-p(3),double(data.z)-p(1),xi,yi,'*linear',nan);
           if sum(~isnan(tz))<=10;idd=[idd,j];continue;end
	   hold all
           plot(dsi,tz,'.-','Linewidth',1)
	   saveas(gcf,'BBprofs','fig')
	end
	Legendstr{j}=ymd;

	if flagdem==1
          OutName=[ymd,'_codem.tif'];
	  projstr='polar stereo north';
	  writeGeotiff(OutName,xout,yout,double(demo.z),5,0,projstr)
	  continue 
	elseif flagdem==2
	%second DEM
	k=j+1;
	%k=3;
        i=idregion(idblock(k));

        demdir=fdir{i};
        infile= [demdir,'/',f{i}];
        fprintf(['\n Loading DEM files in order (',num2str(k),'/',num2str(length(idregion)),'): ',infile])
        [ymd2,~]=strip2date(infile,3);
       % t(k)=datenum(ymd,'yyyymmddHHMMSS');
	XYbi=XYbg{i};

        if strcmp(infile(end-7:end),'meta.txt')
          str1='meta.txt';
        elseif strcmp(infile(end-7:end),'_mdf.txt')
          str1='mdf.txt';
        end
        infile1= strrep(infile,str1,'dem.tif');

        %data=readGeotiff(infile1,'map_subset',rang0);
	[data,~]=readdem(XYbi,infile,rang0,2);

        %Interpolation requires at least two sample points in each dimension.
        [m1,n1]=size(data.z);
        if min([m1,n1])<=1;idd=[idd,k];continue;fprintf(['\n missing one pair:',num2str([j,k])]);end

        %refers to ~/arcticdemapp/river/rivergithub2/rivermainserial/riverprofsub.m
        p=dX4Sg(idblock(k),:); %pzxy

        fprintf(['\n Loading its translational offsets pzxy (',num2str(k),'/',num2str(length(idregion)),'): ',num2str(p)])
        data.z(data.z == -9999) = NaN; % %convert nodata values to nan for interpolation
        z2inp=interp2(data.x-p(2),data.y-p(3),double(data.z)-p(1),xout,yout','*linear',nan);
	%z2inp(isnan(z2inp))=-9999; %return to -9999; 
	clear dem2;dem2.x=xout;dem2.y=yout;dem2.z=z2inp;

	clear dem2m1; dem2m1=dem2;dem2m1.z=dem2.z-dem1.z;
	OutName=[ymd2,'m',ymd,'.tif'];
	projstr='polar stereo north';
	writeGeotiff(OutName,xout,yout,double(dem2m1.z),5,0,projstr)
	end %if flagdem

end
if flagprof==1 %
Legendstr(idd)=[] 
legend(Legendstr)
saveas(gcf,'BBprofs','fig')
end

%[status, cmdout]=system('rm -fr [e-Z]*');
%exit

return
end
