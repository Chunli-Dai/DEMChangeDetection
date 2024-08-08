
%Aug 2023: Filter the shapefiles before ai (e.g., region18/change2polyg_18_m.shp) using three criteria. 
%Output: change2polyg_18_m_cri1to3.shp

constant

codedir='/home/chunlidai/blue/apps/landslide/code1/';  %Directory of codes.
%addpath(genpath([codedir]));
flagcopy=1; %1, copy slumpsub_region18_m/slumppeelsubi1.jpg to slumpsub_region18_m_cri1to3/
	    %0, do not copy files.

%Prepare: download all DEMs from PGC website using Downloaddems.m

%refers to Dividemap_candpar_m.m

[demdir,name,ext] =fileparts([strtrim(regiondir)]);
regid=name(11:12);
regionstr=['region',regid]; %'region09';

classtype_g = {'_vvl', '_vl', '_l', '_m', '_s', '_vs'};

for ic=3:4; %3:4; %1:6

%classtype=classtype_g{4};
   classtype=classtype_g{ic};

if strcmp(classtype,'_vvl')
   classtypein='_landslide';
else
   classtypein=classtype;
end

matfile=['ixy_',regionstr,'_cand',classtype,'.mat']; %'ixy_region18_cand_m.mat'; %['ixy_',regionstr,'_cand',classtype,'.mat'] ;
load(matfile);

shapefile=['change2polyg_',regid,classtypein,'.shp';]
oshapefile=['change2polyg_',regid,classtype,'_cri1to3.shp']; %after cri1to3

S1=shaperead(shapefile);
shpc=struct2cell(S1)';
areaj=[shpc{:,6}]';
M1=areaj>=2000;
S1s=S1(M1);S1=S1s;
if length(S1)~=ns
   fprintf(['\n ',classtype, ': loaded shapefile does not match the ns in the mat file: ',matfile]);
   exit
end

indirt1=['slumptif_',regionstr,classtype,'/'];
indir1=['slumpsub_',regionstr,classtype,'/'];
if flagcopy==1
  odir=['slumpsub_',regionstr,classtype,'_cri1to3/'];
  if ~exist(odir,'')
            mkdir(odir)
  end
end

%ids={};
ids=cell(ns,1);Mpcg2=cell(ns,1);

sz = getenv('SLURM_NTASKS');
sz=str2num(sz);
Mflag=zeros(ns,1); %1 overlap; 0 no overlap;
fprintf(['\n ',num2str(sz),' worker(s) allocated in job.slurm.\n'])
poolobj=parpool(sz);
parfor ix=1:ns

        fprintf(['\n Working on ',num2str(ix),'/',num2str(ns),' candidate cluster.\n'])
        
	ifile=[indirt1,'/','slump',num2str(ix),'.tif'];
	
	if exist(ifile,'file')

        [Mpc2]=getcri1to3sub(ifile);
	Mpc0=Mpcg{ix};
	Mpc2z=interp2(Mpc2.x,Mpc2.y,double(Mpc2.z),Mpc0.x,Mpc0.y','nearest',0);
	Mpc2.x=Mpc0.x;Mpc2.x=Mpc0.x;Mpc2.z=Mpc2z;

        Mpcg2{ix}=Mpc2;

	[n1,m1]=size(Mpcg{ix}.z);
	[n2,m2]=size(Mpcg2{ix}.z);
	if n1~=n2 ||m1~=m2
	    fprintf(['\n Getcri1to3.m ', ifile,' matrix size is not consistent, [n1 m1 n2 m2]=',num2str([n1 m1 n2 m2]),'.\n']) 
	else
% Combined the mask from three criteria with the input mask
	    Mpc2.z=Mpcg2{ix}.z&Mpcg{ix}.z;
   	    %Mcom{ix}=Mpc2; 
	    if sum(Mpc2.z(:))>1; 
		Mflag(ix)=1;

		if flagcopy==1
		%copy slumpsub_region18_m/slumppeelsubi1.jpg to slumpsub_region18_m_cri1to3/
		jfile=[indir1,'/slumppeelsubi',num2str(ix),'.jpg'];
		ofile=strrep(jfile,indir1,odir);

		str1=['cp  ',jfile,' ',ofile] 
		[status, cmdout]=system(str1);
		end

	    else; 
		Mflag(ix)=0;
	    end
	end %if size matches
	
	else

	    fprintf(['\n Getcri1to3.m ', ifile,' do not exsit. \n '])
	end % if exist
end
delete(poolobj)

% Save output shapefiles

Co=S1(Mflag==1);

fprintf(['Combined with criteria 1to3, there are a total of ',num2str(length(Co)), '  ',classtype,' clusters.\n'])

if ~isempty(Co)
 shapewrite(Co, oshapefile); % May 2023
end

ofile=['ixy_',regionstr,'_cand',classtype,'_cri1to3.mat'];
save(ofile,'Mflag', 'ns','-v7.3')

end %ic

% End of main function

