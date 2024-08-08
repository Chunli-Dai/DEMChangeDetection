% if 0
%     %output in ../../plotproftimes.m
% load('../../outsv/version2/jump3b.mat')
% jump=jump3b;
% load('../../outsv/version2/lavamask3.mat'); %mp, a rough mask of main
% % deposit zone;
% %M the exact mask of deposit
% 
% rangeg(:,1)=[ -3436e3     -3434e3    2240e3 2244e3];%zone1; xmin xmax ymin ymax
% rangeg(:,2)=[ -3428e3 -3425e3 2240e3 2244e3];%zone2
% rangeg(:,3)=[ -3428e3 -3425e3   2233e3 2238e3];%zone3
% rangeg(:,4)=[ -3436e3     -3434e3   2233e3 2238e3];%zone4
% 
% end

function [volstd]= corrv2(xout,yout,jump,jumpstd,M,rangeg)
% version 2: do the correlation at four subzones
%Periodogram for correlation of the error of jump
%Function: Given elevation change map, mask of deposit, to estimate the
%           total volume and its uncertainty.
%
% Modified from /Users/chunlidai/Box/ESI2017/TolbachikVolcano/code/corrv2.m
% input: xout, x coordinate in meter, size of 1 by n;
%        yout, y coordinate in meter, size of 1 by n;
%        jump, elevation change in meter, size of n by n;
%        jumpstd, std of jump in meter, size of n by n;
%        M,   the final mask of deposit (1 inside, 0 outside), size of n by n;
%        rangeg, range of four zones that are outside of the main deposit zone, manually given;

macdir='/Users/chunlidai/surge/';
addpath(genpath([macdir,'/home/dai.56/arcticdemapp/volcano/code/']));


% res=2;
resx=mean(xout(2:end)-xout(1:end-1));
resy=mean(yout(2:end)-yout(1:end-1));
res=mean([abs(resx),abs(resy)]);

% reduce size for ploting figures
jumpr=jump(1:5:end,1:5:end);jumpstdr=jumpstd(1:5:end,1:5:end);
[X,Y]=meshgrid(xout(1:5:end),yout(1:5:end));
[LAT,LON]=polarstereo_inv(X,Y,[], [],70,-45);

if 0
out=[LAT(:) LON(:) jumpr(:) jumpstdr(:)];
save -ascii jumpr.dat out
end

figure;surf(LON,LAT,jumpr);shading interp;view(0,90)

%f
% jump(M)=0; %set zeros 
% mp=roipoly; %manually select zones
% [idy,idx]=find(mp==1);
% ranget=[yout(min(idy)),yout(max(idy)),xout(min(idx)),xout(max(idx))]
% 
% figure;imagesc(xout*1e-3,yout*1e-3,jump);colorbar;colormap jet; caxis([-100 100]);
% view(180,-90)

% rangeg(:,1)=[ 3451078     3446000    -1639780    -1633546];%zone1
% rangeg(:,2)=[ 3451722     3446000    -1629076    -1620298];%zone2
% rangeg(:,3)=[ 3462000     3456804    -1639312    -1631624];%zone3
% rangeg(:,4)=[ 3462000     3459490    -1628246    -1623728];%zone4

% k=4;
for k=1:4
ranget=rangeg(:,k);
x0=[ranget(1) ranget(2) ranget(2) ranget(1) ranget(1) ];y0=[ranget(4) ranget(4) ranget(3) ranget(3) ranget(4) ];
figure(5);hold on;plot(x0,y0,'k>-','linewidth',3);
outname=['corzone',num2str(k),'.dat'];
if 0
save -ascii corzone1.dat xylat
save -ascii corzone2.dat xylat
save -ascii corzone3.dat xylat
save -ascii corzone4.dat xylat
save(outname,'xylat','-ascii')
end

idy=find(yout>=ranget(3)&yout<=ranget(4));idx=find(xout>=ranget(1)&xout<=ranget(2));

% dg=jump;
dg=jump(min(idy):max(idy),min(idx):max(idx));
nanmean(dg(:)) %zone 4 -0.7562
nanstd(dg(:))  %zone 4 1.8
rms(dg(:))  %zone 4 1.9
figure;imagesc(xout(min(idx):max(idx))*1e3,yout(min(idy):max(idy))*1e3,dg);colorbar; colormap jet;caxis([-10 10])
view(180,-90)

DG=fft2(dg);dx=res;
[N1,N2]=size(dg);

Rggp=ifft2(conj(DG).*DG/N1/N2);
var1=Rggp(1,1);
var2=var(dg(:));
k1=(0:(N1-1));k2=(0:(N2-1));  %
[K2,K1]=meshgrid(k2,k1);
if 0
figure 
surf(K2,K1,(Rggp))
colormap jet,shading interp
axis equal
view(0,90)
box on
% hold on,[C,h]=contour(X,Y,dg);
xlabel('m2','fontsize',12);ylabel('m1','fontsize',12);
title('2-D auto-correlation of elevation change, unit: m^2','fontsize',14);  %??????
caxis([min(min(Rggp)),max(max(Rggp))]);
% caxis([min(min(log10(Rgg))),max(max(log10(Rgg)))]);
colorbar('horiz');
end

clear dgp
for i=1:round(N1/10)
    dgp(i:i)=Rggp(i,i);
end

cor=[[0:round(N1/10-1)]',dgp'/var1];%pixel distance, correlaton coefficient
% save cor.mat cor

% figure (3)
% hold all
% plot(cor(:,1),cor(:,2),'.-')


figure (17)
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0 0 6 3]); 
set(gcf, 'PaperSize', [ 6 3]); 
hold all;
box on
hold on
plot(cor(:,1)*res,cor(:,2),'-','Linewidth',2)
xlabel('Distance between pixels (meter)','fontsize',12);
ylabel('Correlation coefficient','fontsize',12);


leng(k)=length(cor(:,1));
corg{k}=cor;
end %k

nlen=min(leng);
for k=1:4
corall(:,k)=[corg{k}(1:nlen,2)];
end
corave=mean(corall,2);
cor=[corg{1}(1:nlen,1),corave];
save corave.mat cor

figure (17)
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0 0 6 3]); 
set(gcf, 'PaperSize', [ 6 3]); 
hold all;
box on
hold on
plot(cor(:,1)*res,cor(:,2),'-','Linewidth',2)
xlabel('Distance between pixels (meter)','fontsize',12);
ylabel('Correlation coefficient','fontsize',12);

legend('Area 1','Area 2','Area 3','Area 4','Average')
ofile=['corave'];
saveas(gcf,ofile,'fig')
print('-dpng','-r400',ofile) 
print('-dpdf','-r300',ofile) 


%% Step 2 adjust the std
if 0
load cor.mat;
% load jumpstd3.mat
load jump3.mat
load('lavamask2.mat')
end

% resr=2.;
resr=res;
%assume no correlation between pixels -> underestimate the uncertainty;
volstd=resr*resr*sqrt(sum(jumpstd(M).^2)); %m^3

[idy,idx]=find(M==1);
sigma=jumpstd(M);
npt=length(idy);

save test1.mat -v7.3

%Refers to ~/dai.56/chunli/scripts//Tol40m/corr.m
%/home/dai.56/dai.56/saturationTest/Tolbachik/work2m/corr.m

% cov=zeros(npt,npt);%too big %Mac maximum matrix size: 2e8;
if 0 %old method 1, super slow
    sumcov=0;
    for i=1:npt
        tic
        for j=1:npt
           dist=sqrt((idy(j)-idy(i))^2+(idx(j)-idx(i))^2);
           cori=interp1(cor(:,1),cor(:,2),dist,'linear',0);
           cov=sigma(i)*sigma(j)*cori;
           sumcov=sumcov+ cov;
    %        figure (10); hold on; plot(dist,cori,'b.')

        end
        fprintf([num2str(i),' processed']) 
        toc
    end

elseif 0 %approximate method 2
    %/home/dai.56/dai.56/saturationTest/Tolbachik/work2m/corr.m

    maxlag=200;%pixels
    sumcov=zeros(npt,1);
    di=400; %npt/400=30000; should be good approximation; tested in scripts/Tol40m/corr.m
    di=100;
    di=10;
    tic
    for i=1:di:npt
        flag=abs(idy-idy(i))<maxlag&abs(idx-idx(i))<maxlag;
        dist=sqrt((idy(flag)-idy(i)).^2+(idx(flag)-idx(i)).^2);
        sigmat=sigma(i)*sigma(flag);
        cori=interp1(cor(:,1),cor(:,2),dist,'linear',0);
        sumcov(i)=sigmat'*cori;
    end
    toc %3 minutes di=400; 11 minutes di =100;
    sumcov=sum(sumcov)*di;


elseif 0 % if npt*npt< 2e8/10 new method 3, matrix calculation, SUPER fast and accurate. but matrix too big -> crashes!
    tic
    pt=[idx,idy];
    D = pdist(pt);dist = squareform(D);    
    cori=interp1(cor(:,1),cor(:,2),dist,'linear',0);
    sigmaij=sigma(1:npt)*sigma(1:npt)';
    cov=sigmaij.*cori;
    sumcov=sum(cov(:));
    toc
% elseif 0
%     %a mix; matrix calculation for each i, fast. Not fast enough; 12 days;
%     sumcov=0;
%     for i=1:npt
% %         tic
%         J=1:npt;
%         dist=sqrt((idy(J)-idy(i)).^2+(idx(J)-idx(i)).^2);
%         cori=interp1(cor(:,1),cor(:,2),dist,'linear',0);
%         cov=sigma(i).*sigma(J).*cori;
%         sumcov=sumcov+ sum(cov(:));
% %         fprintf([num2str(i),' processed']) 
% %         toc
%     end
    
elseif 1 % mix method 4, matrix calculation, accurate and fast. No crashes; 
        % 10 days;
        %approximate with maxlag -> 2 days
    tic
    matrixmax=2e8;
    nmax=round((matrixmax/10)/npt);
    maxlag=npt+1; %300;

    JJ=1:npt;
    x2=idx(JJ);y2=idy(JJ);
    
    dn=20;
    idII=1:dn:npt; nblock=length(idII);
    IIs=zeros(size(idII));IIe=IIs;
    for ipt=1:length(idII)
    IIs(ipt)=idII(ipt);
    IIe(ipt)=min([npt,ipt*dn]);
    end
    
    sumcov=0;
    for ipt=1:length(idII)
%         tic
        II=IIs(ipt):IIe(ipt);
        x1=idx(II);y1=idy(II);
        
        [dist] = pdist2([x1, y1], [x2, y2]);
        flag=abs(dist)<maxlag;

%         cori=interp1(cor(:,1),cor(:,2),dist,'linear',0);
        cori=interp1(cor(:,1),cor(:,2),dist(flag),'linear',0);

        sigmaij=sigma(II(:))*sigma(JJ(:))';
%         cov=sigmaij.*cori;
        cov=sigmaij(flag).*cori;
        sumcov=sumcov+sum(cov(:));
%         toc
    end
    toc
end

volstd2=resr*resr*sqrt(sumcov)
volstd=volstd2;
%approximate method: 2.6694e+06 m^3 (di=400); 2.6764e6 (di=100)
%accurate method:
%landslide volume std 3.6181e5 m^3;
%subsidence volume std 9.6248e+05 (di=100; time: 1 minute); 9.5929e+05 (di=10; time 10 minutes);
			 (method 4, accurate; time: ) m^3;

if 0
% plot figure
figure;
plot(cor(:,1),cor(:,2),'.-')

figure;imagesc(xout*1e-3,yout*1e-3,jump);colorbar;colormap jet; caxis([-100 100]);
view(180,-90)
Md1 = imdilate(M, ones(3));
Me=logical(Md1-M);
[X,Y]=meshgrid(xout,yout);
hold on;plot(X(Me)*1e-3,Y(Me)*1e-3,'k.')

load('outline3.dat');
figure;plot(outline3(:,2),outline3(:,1),'k.')
axis equal
end

end %end of function
