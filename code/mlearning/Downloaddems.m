
%Download DEMs for a region.

codedir='/home/chunlidai/blue/apps/landslide/code1/';  %Directory of codes.

addpath(genpath([codedir]));

constant
%tiledir='/blue/chunlidai/data/ArcticDEMmosaics/';

%rang0=[];
%[data0]=getmosaic(rang0,0);
%filename='../arcticdem_18_russia_cherskly/tilelistb_all';
%filename='../../arcticdem_21_russia_yakutiya_east/tilelistb';
filename='/orange/chunlidai/results/landslide/arcticdem_23_russia_yakutiya_west/tilelistb';

fid = fopen(filename);
n = linecount(fid);
fid = fopen(filename);
for i=1:n
    tilefileg{i}=[fgetl(fid)];
end

sz = getenv('SLURM_NTASKS');
sz=str2num(sz);
fprintf(['\n ',num2str(sz),' worker(s) allocated in job.slurm.\n'])
poolobj=parpool(sz);
parfor i=1:n
    tilefile=tilefileg{i};

    [dir,ifile,ext] =fileparts(tilefile);
    r=1;
    xid= sscanf(ifile(r+3:(r+4)), '%g', 1);
    yid= sscanf(ifile(r:(r+1)), '%g', 1);
    xids= sscanf(ifile(r+6), '%g', 1);
    yids= sscanf(ifile(r+8), '%g', 1);

    %07_42_2_1_2m_v4.1/07_42_2_1_2m_v4.1_dem.tif

    %tilefile=sprintf('%02d_%02d_%01d_%01d_2m_v3.0_reg_dem.tif',yid,xid,xids,yids);
    tilefile=sprintf('%02d_%02d_%01d_%01d_2m_v4.1_dem.tif',yid,xid,xids,yids);

        % find the dem tile file or download the data
    [status , cmdout ]=system(['find ',tiledir,' -name ',tilefile]); %status always 0, cmdout could be empty.
    if ~isempty(cmdout) && status ==0 % 
        tilefile=deblank(cmdout);
    else
        warning(['Tile file ',tilefile,' not found! Download it from website'])
        [dir,name,ext] =fileparts(tilefile);
    %     http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v2.0/43_59/43_59_2_1_5m_v2.0.tar
        tarfile=[name(1:17),'.tar'];
%       webtile=deblank(['http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v2.0/',name(1:5),'/',tarfile,'.gz']);
        %webtile=deblank(['http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v3.0/2m/',name(1:5),'/',tarfile,'.gz   --no-check-certificate']);
        webtile=deblank(['http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v4.1/2m/',name(1:5),'/',tarfile,'.gz   --no-check-certificate']) 
        %webtile=deblank(['http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v4.1/10m/',name(1:5),'/',tarfile,'.gz   --no-check-certificate']) 
        system(['wget  ',webtile,' -o downloadlog'])
        system(['gunzip -f ',tarfile,'.gz']);
        system(['tar -xvf ',tarfile]);
        system(['rm  ',tarfile]);
        system(['rm  ',name(1:17),'_mad.tif']);
        %collect all downloaded mosaic dems to tiledir/ % to do
        %system(['mv *dem_meta.txt  *reg.txt *.tar *_reg_dem.tif *_reg_matchtag.tif ',tiledir]);
        system(['mv ',name(1:17),'* ', tiledir]);
        [status , cmdout ]=system(['find ',tiledir,' -name ',tilefile]);
        tilefile=deblank(cmdout);
        if ~exist(tilefile,'file')
          % continue
        end
    end
end %for i

delete(poolobj)


