function	  [f, fdir,idregion,dX4Sg,oflag]=getdX4Sg(listfile);
% Get dX4Sg from listused.txt file.
% oflag: 0 bad, 1 successful

%  listfile=[dir1,tilenamell,'_listused.txt']; %34_23_2_2/34_23_2_2_25_25/34_23_2_2_25_25_listused.txt

%# px     py   pz   DEMFileName 
%# Strip DEM - [px py pz] =Coordinates aligned with the reference DEM (reference in bundle adjustment). 
%        0.09        -1.06         0.05  /scratch/sciteam/GS_bazu/mosaic_data/setsm_domain/ArcticDEM/region/arcticdem_07_canada_ellesmere/strips_v4/2m/WV01_20120414_102001001BA50B00_1020010017263600_2m_lsf_v030000//WV01_20120414_102001001BA50B00_1020010017263600_2m_lsf_seg1_meta.txt 
%        0.85        -0.07        -0.06  /scratch/sciteam/GS_bazu/mosaic_data/setsm_domain/ArcticDEM/region/arcticdem_07_canada_ellesmere/strips_v4/2m/WV01_20120414_102001001BA50B00_1020010017263600_2m_lsf_v040306//WV01_20120414_102001001BA50B00_1020010017263600_2m_lsf_seg1_meta.txt 

oflag=1;

try 
fid = fopen(listfile, 'rt');
listcell = textscan(fid, '%f %f %f %s ', 'headerlines', 2);
fclose(fid);

t1=cell2mat(listcell(1:3));  %px py pz
dX4Sg=[t1(:,3), t1(:,1:2)]; % z x y

n=length(listcell{1,4})
idregion=1:n;

fdir=cell(n,1); f=cell(n,1);

for i=1:n
    ifile=listcell{1,4}{i};	
   [demdir,name,ext] =fileparts([strtrim(ifile)]);
   f{i}=[name,ext];
   fdir{i}=[demdir,'/'];
end

catch e
   fprintf('getdX4Sg.m There was an error! The message was:\n%s',e.message);
   oflag=0;
   f=[];fdir=[]; idregion=[];dX4Sg=[];

end

return
end
