function [co]=mergejson(indir,odir)
%Merge two sets of training images and json files
% indir: a list of folder names to be merged
% odir: output folder name
%Input example:
%indir={'./slumpsIngmar/','../banks/pgcnew/slumpsbanks_band1val/','../ashleytuk/slumps/'};
%odir='./slumpsmerged/';
co=[];

odir1=[odir,'/train/'];odir2=[odir,'/val/'];
%if ~exist(odir,'dir')
    mkdir(odir)
    mkdir(odir1)
    mkdir(odir2)
%end

nk=length(indir);
strg_tra='';
strg_val='';
for k=1:nk
    groupid=['g',num2str(k),'.jpg'];
    idir=indir{k};

    for itv=1:2 %1 train; 2 val
	if itv ==1 %train
	   str_tv='train';
	   odiri=odir1;
	elseif itv==2 %val
	   str_tv='val';
	   odiri=odir2;
	end
    %cp slump1.jpg slump1g1.jpg
%    str1=['cp ',idir,'/train/*jpg ',odir1];
%    [status, cmdout]=system(str1);
    filename='listi';
    [status, cmdout]=system(['rm -f ',filename]);
    [status, cmdout]=system(['ls ',idir,'/',str_tv,'/*jpg >', filename ]);

    fid = fopen(filename);
    n = linecount(fid);
    fid = fopen(filename);
    for i=1:n
       ifile=[fgetl(fid)];
       [demdir,name,ext] =fileparts([strtrim(ifile)]);
       ofile=[odiri,name,ext];
       ofile=strrep(ofile,'.jpg',groupid);
       str1=['cp ',ifile,' ',ofile];
       [status, cmdout]=system(str1);
    end
    fclose(fid); 

    jsonfile=[idir,'/',str_tv,'/via_region_data.json'];
    if exist(jsonfile,'file')
    fid = fopen(jsonfile); 
    raw = fread(fid,inf); 
    str = char(raw'); 
    str =strrep(str,'.jpg',groupid);
    fclose(fid); 
    else; str=[];
    end

    % connect the start and end.
    if ~isempty(str)
    % all jason files before the last one.
    if k ~= nk
	str(end)=',';
    end
    % all jason files after the first one.
    if k ~=1
	str(1)='';
    end
    end

    if itv ==1 %train
    strg_tra=[strg_tra,str];
    elseif itv ==2 % val
    strg_val=[strg_val,str];
    end

    %val = jsondecode(str); %string to structure; some characters are modified, e.g. 0, dot.

    % Do the same for validation data.
    %str1=['cp ',idir,'/val/*jpg ',odir2];
	end %    for itv=1:2 %1 train; 2 val

end %k

    %write merged json files
    ofile=[odir1,'/via_region_data.json'];
    fid=fopen(ofile,'w');
    fprintf(fid, strg_tra);
    fclose(fid);

    ofile=[odir2,'/via_region_data.json'];
    fid=fopen(ofile,'w');
    fprintf(fid, strg_val);
    fclose(fid);

return
end
