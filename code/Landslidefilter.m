function [co]=Landslidefilter(xout,yout,jump,jumpstd,timec,timestd,oflagc,odir)
%[co]=Landslidefilter(xout,yout,jump,jumpstd,timec,timestd,oflagc,odir);
% filter the jump
% load sv1.mat

constant %smlarea
Mstrip.x=xout;Mstrip.y=yout;
resx=mean(Mstrip.x(2:end)-Mstrip.x(1:end-1));resy=mean(Mstrip.y(2:end)-Mstrip.y(1:end-1));resr=mean([abs(resx),abs(resy)]);

Mf1=(timec~=0)&(abs(jump)>=3*jumpstd); % jump detected and significant;
% Mf1=oflagc==0|oflagc==1|oflagc==3|oflagc==4; % no data or no jump
% event time not homogeneous;
% using convolution to compute std, faster, 6 sec (vs 8 sec); Reference:/data/chunli/scripts/watermask.m 
kernel_size = 3; % dx=11*resr;
if 0 %old version
timec2=timec;timec2(timec==0)=min(timec(timec~=0)); 
uxmn=conv2(timec2,ones(kernel_size)/(kernel_size.^2),'same'); % obtain mean
uxstd=sqrt( abs(conv2(timec2.^2,ones(kernel_size)/(kernel_size.^2),'same') - uxmn.^2 )); % obtain standard deviation
else %try nan for no jump
timec2=timec;timec2(timec==0)=nan;
k=ones(kernel_size)/(kernel_size.^2);
uxmn = nanconv(timec2, k, 'noedge');
uxstd=sqrt( abs(nanconv(timec2.^2,k,'noedge') - uxmn.^2 )); % obtain standard deviation
end

Mf2=uxstd<1;% event time  homogeneous;
Mf3=(jumpstd) < 2; %10; % std reasonably small.
Mf4=jump<=-2; %abs(jump)>=2; %jump significant;Only negative

oflagbit=zeros(6,1);oflagbit(5)=1;
oflag0=bit2value(oflagbit);
%Mf5=(oflagc==5|oflagc==6);%oflagc==2|oflagc==6; % jump is at end or start.
Mf5=(oflagc==oflag0);

Mfc=Mf1&Mf2&Mf3&Mf4&Mf5;

Modj=Mfc;
Mfcb= bwareaopen(Modj, round(smlarea/resr/resr)); %remove small clusters
Mf=Mfcb;

% step for calculating local slopes....
if 0
data=readGeotiff(tilefile);
resgrad=data.info.map_info.dx;
 % get the slope
[px,py] = gradient(double(data.z),resgrad); % second parameter is the image resolution 
py=-py; % minus sign in front of py, since data.z has the y axis positive downwards.
[~,rho] = cart2pol(px,py);
slope=atan(double(rho))*180./pi;
sloper =interp2(data.x,data.y,slope,xout,yout','*linear');%toc;%4s
Mf6=sloper<10;%unlikely to have landslide.
end

% Mf=Mf4|Mf5;
% the high possibility based on DEM slope

%% Reflectance change detection
% transform DN to reflectance for panchromatic band only.

for jj=1:0 %[6,13:15] % 6 vs 15
%ymd=f{id(jj)}(6:13);
filename=f{id(jj)};ymd=filename2ymd(filename);
i=id(jj);  
infile= strrep([demdir,'/',f{i}],'meta.txt','ortho.tif');
data=readGeotiff(infile,'map_subset',rang0st);

metafile=strrep(infile,'ortho.tif','meta.txt');
if 0
%from scene filename get the strip meta.txt file
	mfile=strrep(metafile,'tif_results','strips');
	[demdir,filename,~]=fileparts(mfile);
	mfile=([demdir,'/',filename(1:48),char(42),'meta.txt']);%find the strip that contain the scene.
    [status , cmdout ]=system(['ls ',mfile]);
    if status ==0
        mfile=deblank(cmdout);
        c=textread((mfile),'%s','delimiter','\n');
        r=find(~cellfun(@isempty,strfind(c,filename(1:end-4))));
        if ~isempty(r)
        %found the strip metafile mfile
        end
    end% if status==0
    
end

[~,filename,~]=fileparts(metafile);
%assume in each strip, all scenes and each image in pairs have the same
%effectivebandwith and abscalfactor.
c=textread(metafile,'%s','delimiter','\n');
strg={['_abscalfact='],['_effbw='],['_Mean_sun_elevation=']};
for j=1:3
str=strg(j);
r=find(~cellfun(@isempty,strfind(c,str)));
str=c{r(1)};r1=strfind(str,'=');str(1:r1)='';
% Xbs=deblank(strrep(c1{r(1)},str,''));
t1(j)= sscanf(str, '%g', 1);
end
abscalfactor=t1(1);effectivebandwith=t1(2);meanSunEl=t1(3);
Theta=(90.-meanSunEl)*pi/180.;
% sun-earth distance polynomial function coefficients
% doy = day of the year: 
year=sscanf(filename(6:9), '%g', 1); month=sscanf(filename(10:11), '%g', 1); day=sscanf(filename(12:13), '%g', 1);
doy=juliandate(year,month,day)-juliandate(year,1,1)+1;
C = [1.8739e-26,-3.4455e-23,2.7359e-20,-1.2296e-17,3.0855e-15,-2.2412e-13,-5.8744e-11,6.9972e-10,2.5475e-06,-1.6415e-05,0.9833];
dES = polyval(C,doy); % earth-sun distance

satname=filename(1:4);
[GAIN,OFFSET,Esun]=readgainoffset(satname); %read Gain OFFset data, GainOffset.txt.

DN=double(data.z(:,:));
L=GAIN(end)*DN*(abscalfactor/effectivebandwith)+OFFSET(end);
rho=L*dES^2*pi/(Esun(end)*cos(Theta));


figure;set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.25 2.5 6 4]);
% surf(data.x*1e-3,data.y*1e-3,data.z);colormap gray; colorbar; shading interp; axis tight; view(0,90)
imagesc(data.x*1e-3,data.y*1e-3,rho);colormap gray; colorbar; shading interp; axis tight; view(0,-90)
hold all
plot(x0st*1e-3, y0st*1e-3,'r-','linewidth',6)
plot(xeq*1e-3, yeq*1e-3 ,'r*','Markersize',12)
axis equal;title([ymd,';j=',num2str(jj)]);
end


if 0
figure;set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]); 
hold all
xt=xout*1e-3;yt=yout*1e-3;
% imagesc(xt,yt,uxstd,'alphadata',~Mf2); colormap jet;colorbar; shading interp; axis equal; view(0,90)  
imagesc(xt,yt,uxstd); colormap jet;colorbar; shading interp; axis equal; view(0,90)  
set(gca,'FontSize', 18);
plot(xeq*1e-3, yeq*1e-3 ,'r*','Markersize',12)
xlabel('{\itx} (km)')
ylabel('y (km)')
% caxis([2008 2018 ])
% caxis([-1 1])
axis(rang0*1e-3)
end

figure;set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]); 
hold all
xt=xout*1e-3;yt=yout*1e-3;
imagesc(xt,yt,jump,'alphadata',Mf); colormap jet;colorbar; shading interp; axis equal; view(0,90)  
set(gca,'FontSize', 18);
% plot(xeq*1e-3, yeq*1e-3 ,'r*','Markersize',12)
xlabel('{\itx} (km)')
ylabel('y (km)')
caxis([-2 2])
caxis([-50 50])
% axis(rang0*1e-3)
title('Elevation change; flag=5, filtered ')
print('-dpdf','-r300','jumpfiltered')   
saveas(gcf,'jumpfiltered','fig')

[X,Y]=meshgrid(xout,yout);
[LAT,LON]=polarstereo_inv(X,Y,[], [],70,-45);
% id=X>=rang0b(1)&X<=rang0b(2)&Y>=rang0b(3)&Y<=rang0b(4);
out=[LAT(Mf) LON(Mf) jump(Mf) jumpstd(Mf) timec(Mf)];
% save -ascii latlonnov_barnesstrip8m.dat out
fid = fopen('jumpfiltered.dat','a'); %statistics
fprintf(fid,' %f %f %f %f %f \n',out');
fclose(fid);

%% write output
projstr='polar stereo north';
OutName=[odir,'_bitmask.tif']; % %1 is good (landslide); 0 is bad
writeGeotiff(OutName,xout,yout,int32(Mf),3,0,projstr)

co=[];

return

end
