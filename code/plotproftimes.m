function [Co]=plotproftimes(BB,xout,yout,jump3,jumpstd3,timec,timestd,data,dataaf,xq1,yq1)
% Given a line, plot the height profiles along this line
%from /home/dai.56/chunliwork/volcano/okmok/plotproftimes.m
%Input:BB, x,y coordinates of the given line (in meter);
%data: pre-event reference DEM;
%dataaf: after-event DEM
% xq1,yq1; points for check time series.
%output: figures of profiles

%addpath(genpath(['/Users/chunlidai/surge/home/dai.56/arcticdemapp/river/rivergithub2/']))

% Peak    5.5766804e+01   1.6031428e+02   1.0173529e+02   5.7031351e+00
Co=[];

if 0 %test
BB=getline; %manually draw lines.
BB=BB*1e3;

load sv1.mat %get xout, yout, jump
%a DEM file
infile='stripdata/SETSM_WV01_20140813_1020010034A1ED00_10200100319F9400_seg1_2m_v3.0_dem.tif';
% data=readGeotiff(infile,'map_subset',rang0);
rang0=[ -711300     -711000     -821100     -820900];
data=readGeotiff(infile,'map_subset',rang0);

rang0=[-711.4e3-500 -710.4e3+500 -821.9e3-500 -820.9e3+500]; %1km by 1km
% infile='../site1Eureka/stripdata/SETSM_WV01_20090810_1020010009A0A000_1020010008227C00_seg1_2m_v3.0_dem.tif';% %bad:large translational offset wrt the reference
infile='../site1Eureka/stripdata/SETSM_GE01_20090804_10504100049D6800_1050410002F81C00_seg1_2m_v3.0_dem.tif';
data=readGeotiff(infile,'map_subset',rang0);

infile='../site1Eureka/stripdata/SETSM_WV01_20170415_102001005F334C00_10200100601E3900_seg1_2m_v3.0_dem.tif';
dataaf=readGeotiff(infile,'map_subset',rang0);

end

pt=BB;%cooridnates in m

%densify the line
clx=pt(:,1);cly=pt(:,2);dc=2.;
[clx2,cly2,S2]=interpcl(clx,cly,dc);

dsi=S2;xi=clx2;yi=cly2;
% id=find(dsi<1e3);
% xi(id)=[];yi(id)=[];tz(id)=[];
hold on;plot(xi*1e-3,yi*1e-3,'b-','Markersize',12,'Linewidth',4)
[lati,loni]=polarstereo_inv(xi,yi,[], [],70,-45);
out=[loni(:),lati(:),xi(:),yi(:)];
save -ascii BB.dat out

%elevation change
tz =interp2(xout,yout,double(jump3),xi,yi,'*linear');
tzstd =interp2(xout,yout,double(jumpstd3),xi,yi,'*linear');

figure;plot(dsi,tz)
hold on;shadedErrorBar(dsi,tz,tzstd,'b',1)

%before and after event topography
th =interp2(data.x,data.y,double(data.z),xi,yi,'*linear');
thaf =interp2(dataaf.x,dataaf.y,double(dataaf.z),xi,yi,'*linear');%note this DEM might not be coregistered
hold on;plot(dsi,th,'-')

% event time
timep =interp2(xout,yout,double(timec),xi,yi,'*linear');
timepstd =interp2(xout,yout,double(timestd),xi,yi,'*linear');

%Plot elevation profile;
%note th+tz may not represent the actual after event topography
figure;
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0 0 3 1.5]); 
set(gcf, 'PaperSize', [ 3 1.5]); 
basevalue = 1540;
hold on;h=area(dsi,th);
h.FaceColor=[220,220,220]/255; %[0.5 0.5 0.5];
plot(dsi,th,'k-','Linewidth',1)
plot(dsi,th+tz,'r-','Linewidth',1) 
hold on;shadedErrorBar(dsi,th+tz,tzstd,'r',1)
plot(dsi,thaf,'g-','Linewidth',1)
box on
yyaxis left
xlabel('Distance (m)')
ylabel('Height (m)')
text(50,1590,'$B$','Interpreter','latex')
text(1350,1690,'$B^{\prime}$','Interpreter','latex')
axis([0 1500 1540 1750])
print('-dpng','-r400','topoprof') 
saveas(gcf,'topoprof','fig')

%Plot topography on left, event time profile on right
fig=figure;
left_color = [0 0 0]; %black
right_color = [0 0 1]; %b
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
hold on;h=area(dsi,th);
h.FaceColor=[220,220,220]/255; %[0.5 0.5 0.5];
plot(dsi,th,'k-','Linewidth',1)
hold all;
yyaxis left
xlabel('Distance (m)')
ylabel('Height (m)')
yyaxis right
plot(dsi,timep,'b.-')
hold on;shadedErrorBar(dsi,timep,timepstd,'b',1)
box on
yyaxis right
xlabel('Distance (m)')
ylabel('Event time (year)')
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0 0 3 1.5]); 
set(gcf, 'PaperSize', [ 3 1.5]); 
ofile='timeprof';
print('-dpng','-r400',ofile) 
saveas(gcf,ofile,'fig')

%plot the location of xq1, yq1 on the profile
for j=[2,8] %[4,9] %[1:length(xq1)]
    xeq=xq1(j);yeq=yq1(j);
    dist=sqrt( (yeq-cly2(:)).^2 + (xeq - clx2(:)).^2);
    [distmin,k]=min(dist);
    xobs2=S2(k);
%     hold on;plot(([xobs2,xobs2]),[20 65],'k-')
    hold on;plot(([xobs2,xobs2]),[10 56],'k-')  %[-12 4]
%     hold on;plot(([xobs2,xobs2]),[-12 4],'k-')
end
id=[1:length(xq1)];
[latq1,lonq1]=polarstereo_inv(xq1,yq1,[], [],70,-45); %[60.1768700940536 -141.184899979958]
out=[lonq1(:),latq1(:), xq1(:),yq1(:),id(:)];
save -ascii lonlatq1.dat out

%plot topography on the left and elevation change on the right yaxis
fig = figure;
left_color = [0 0 0]; %black
right_color = [1 0 0]; %r
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
hold on;h=area(dsi,th);
h.FaceColor=[220,220,220]/255; %[0.5 0.5 0.5];
plot(dsi,th,'k-','Linewidth',1)
yyaxis left
xlabel('Distance (m)')
ylabel('Height (m)')

yyaxis right
plot(dsi,tz,'r-','Linewidth',1)
hold on;shadedErrorBar(dsi,tz,tzstd,'r',1)
plot(dsi,zeros(size(dsi)),'r--','Linewidth',1)
box on
yyaxis right
xlabel('Distance (m)')
ylabel('Elevation change (m)')
% axis([0 1500 1540 1750])
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0 0 3 1.5]); 
set(gcf, 'PaperSize', [ 3 1.5]); 
print('-dpng','-r400','heightelevchange') 
saveas(gcf,'heightelevchange','fig')

%plot elevation change on the left yaxis, and time on the right yaxis
fig = figure;
left_color = [1 0 0]; %r [1 0 0]
right_color = [0 0 1]; %b
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(dsi,tz,'r-','Linewidth',1)
hold on;shadedErrorBar(dsi,tz,tzstd,'r',1)
plot(dsi,zeros(size(dsi)),'r--','Linewidth',1)
box on
yyaxis left
xlabel('Distance (m)')
ylabel('Elevation change (m)')

yyaxis right
plot(dsi,timep,'b.-')
hold on;shadedErrorBar(dsi,timep,timepstd,'b',1)
box on
yyaxis right
xlabel('Distance (m)')
ylabel('Event time (year)')
% axis([0 1500 1540 1750])
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0 0 3 1.5]); 
set(gcf, 'PaperSize', [ 3 1.5]); 
print('-dpng','-r400','Elevchangevstime') 
saveas(gcf,'Elevchangevstime','fig')


end

