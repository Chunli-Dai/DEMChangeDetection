function [XYbg,projgdalg,f,fdir]=creatsublist(rang0,XYbg,projgdalg,f,fdir) 
%create file 'imagelist.dat' which include lists of images only within the selected tile.
constant

fprintf('\n Select images within the given tile, which takes:');
tic
x0=[rang0(1) rang0(2) rang0(2) rang0(1) rang0(1) ];y0=[rang0(4) rang0(4) rang0(3) rang0(3) rang0(4) ];
XYb0=[x0(:),y0(:)];
overlap_ratio1=zeros(length(fdir),1);
for j=1:length(fdir)
	projgdalj=projgdalg{j}; % e.g.,'epsg:32637'
	XYbj=XYbg{j};
	%
	if ~strcmp(projgdalj,projgdal);
        %  fprintf(['\n Image has different projection than that in the constant.m: ',ifile,', ' , projgdalj,' \n']);
        %  fprintf(['\n Converting the projection for the boundary. \n']);
	   xj=XYbj(:,1);yj=XYbj(:,2);
	%  xy to lat lon
	   [latj,lonj]=xy2latlon(xj,yj,projgdalj);

	%  lat lon to xy
           [xj,yj]=latlon2xy(latj,lonj,projgdal); 
	   XYbj=[xj(:),yj(:)];
        end

	%
	[overlap_ratio1(j),~,~]=getoverlap(XYb0,XYbj);
end

thres=20; %20; 0 include all, 20 

if exist('coverthres','var')
    thres=coverthres;
end

M1=overlap_ratio1<=thres;

XYbg(M1)=[]; f(M1)=[];fdir(M1)=[];projgdalg(M1)=[];

%in case of too many data, use only the summer and latest (in time) 100 measurements.
%Purpose: avoid snow issues.
%And keep only one image for each month. - not applied
for j=1:length(fdir)
            str{j}=['XXXX_',filename2date(f{(j)})]; %'WV02_20150617151345'
end
[idd,idkp]=selectdate(str,100);
XYbg(idd)=[]; f(idd)=[];fdir(idd)=[];projgdalg(idd)=[];

filename=[currentdir,'/imagelist.dat'];
fid = fopen(filename,'w');
for j=1:length(fdir)
   ifile=[fdir{j},'/',f{j}];
   fprintf(fid,'%s \n',ifile);
end
fclose(fid);

toc

return
end
