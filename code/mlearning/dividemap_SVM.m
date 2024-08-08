
constant
codedir='/home/dai.56/arcticdemapp/landslide/code1/';  %Directory of codes.
flagorthotif=0; % save mosacked ortho image in Geotiff

regionstr='peel';

if flagorthotif==1
        %save orthoimage in geotiff 
        odirt1=['slumps_ortho/'];
        if ~exist(odirt1,'')
            mkdir(odirt1)
        end
end
%mkdir slumpsub_3bandsuint8jpg

addpath(genpath([codedir]));

load mat0.mat

% Divide target image to ~1 km size;
imagefiles_val={'/Users/chunlidai/Downloads/peel_pgcnew_matchfilter10m_allseasons/41_16_1_1_merge.tif'};
imagefiles_val={'41_16_1_1_merge.tif'};
%imagefiles_val={'32_33_2_2_jump.tif'};
changefile=imagefiles_val{1};
data0=readGeotiff(changefile);
resr=mean(diff(data0.x));

%only work on test zone in peel
Mx=data0.x <-2495300;
data0.x(~Mx)=[];data0.z(:,~Mx)=[];

rangi=[min(data0.x),max(data0.x), min(data0.y),max(data0.y)];
tic
[datao]=box2mosaic(rangi,3,f,fdir); % 100 m buffer zone.
toc

     %get orthoimage and DEM
     dem.x=datao.x;dem.y=datao.y;ortho=dem;delev=dem; %initialize
     dem.z=datao.z(:,:,2);ortho.z=datao.z(:,:,3);delev.z=datao.z(:,:,1);

        %save orthoimage in geotiff 
        if flagorthotif==1
        filenamet1i=[odirt1,'/','slump',num2str(ixy),'.tif'];
        projstr='polar stereo north';
        writeGeotiff(filenamet1i,datao.x,datao.y,uint16(datao.z(:,:,3)),2,0,projstr)
        end


     %make sure to set void data to NaNs;
     M1=delev.z==0;delev.z(M1)=nan;
     M1=dem.z==-9999; dem.z(M1)=nan;
     M1=ortho.z==0; ortho.z(M1)=nan;
     %data.z=delev.z;

save test_d_svm_part1_v2.mat -v7.3

%step: 500 m ; buffer 200 m on each size. for Peel Plateau
resr2m=2;
buff=0; 
%dx=500;
nw=9; 
dx=nw*resr2m;
dy=-dx;
% try 800 pixels (2m ): 1600 m ; For Eureka
%buff=200; dx=1200;dy=-dx;
%buff=100; dx=300;dy=-dx; %-> Chandi uses 250 pixels by 250 pixels;

rangx=data0.x(1):dx:data0.x(end);nx=length(rangx)-1;
rangy=data0.y(1):dy:data0.y(end);ny=length(rangy)-1;
ns=nx*ny;
ids=[];
featgt_par=cell(ny,nx);
countf=0;
poolsize=20;
poolobj=parpool(20);
parfor ix=1:nx
    for iy=1:ny
        ixy=iy+(ix-1)*ny;
	fprintf(['\n Working on ',num2str(ixy),'/',num2str(ns),' tiles.\n'])
        rangi=[rangx(ix)-buff rangx(ix+1)+buff rangy(iy+1)-buff rangy(iy)+buff];

	Mx=datao.x>=rangi(1)&datao.x<rangi(2);My=datao.y>=rangi(3)&datao.y<rangi(4);

            dem_win=[];delev_win=[];ortho_win=[];
            dem_win.x=dem.x(Mx);dem_win.y=dem.y(My);dem_win.z=dem.z(My,Mx);
            delev_win=dem_win; ortho_win=dem_win;
            delev_win.z=delev.z(My,Mx); ortho_win.z=ortho.z(My,Mx);


	M1=~(isnan(delev_win.z)|delev_win.z==0);
	n_good=sum(M1(:));
	nw_total=nw*nw;

	if n_good/nw_total<=0.8 
	   continue
	end

        [feati,featlabels]=extractfeature(dem_win,ortho_win,delev_win,1);        
	if 0 % not parallel
        %countf=countf+1;
        %featgt(:,countf)=feati;
	else %parallel
        featgt_par{iy,ix}=feati;
	end

    end
end
delete(poolobj)

if 0 % not needed
featgt=[];
countf=0;
for ix=1:nx
    for iy=1:ny
	if ~isempty(featgt_par{iy,ix})
	countf=countf+1;
	featgt(:,countf)=featgt_par{iy,ix};
	end
    end
end
featg=featgt';
end %if 0

save test_d_svm_v2.mat -v7.3

%save ixy.mat data0 ids buff -v7.3
save ixy_peel_svm.mat data0 ids buff dx rangx rangy ns featgt_par featg -v7.3
%save ixy_peel_svm.mat data0 ids buff dx rangx rangy ns featgt_par featg featlabels -v7.3


