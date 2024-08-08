function creatjson(filenames,filesize,XYb,ofile,category_idgc)
%creat structure and save it to a json file
%output similar to file:///Users/chunlidai/ArchivedBoxSync/ESI2019/machinelearning/balloon/val/via_region_data.json
% Input: filenames, filesize, polygon x, polygon y, ofilename

if 0 %test
str1={'a24631331976_defa3bb61f_k.jpg668058','teststr2','teststr3'};
str1{4}='teststr4';

%replace forbidden characters in the structure.
str1=strrep(str1,'.','dot');
% clear s1
% append new fields
for i=1:length(str1)
    field=str1{i};
    zero='zero';
    s1.(field).fileref='';
    s1.(field).size=668058;
    s1.(field).filename='24631331976_defa3bb61f_k.jpg';
    s1.(field).base64_img_data='';
    s1.(field).file_attributes=struct();
    
    s1.(field).regions.(zero).shape_attributes.name='polygon';
    s1.(field).regions.(zero).shape_attributes.all_points_x=[916,913,905,889,868,836,809,792,789,784,777,769,767,777,786,791,769,739,714,678,645,615,595,583,580,584,595,614,645,676,716,769,815,849,875,900,916,916];
    s1.(field).regions.(zero).shape_attributes.all_points_y=[515,583,616,656,696,737,753,767,777,785,785,778,768,766,760,755,755,743,728,702,670,629,588,539,500,458,425,394,360,342,329,331,347,371,398,442,504,515];
    
    s1.(field).regions.(zero).region_attributes=struct();
end
j1=jsonencode(s1);

%replace forbidden characters in the structure.
j1=strrep(j1,'dot','.');
j1=strrep(j1,'"zero"','"0"');
j1=strrep(j1,'a2','2');

fid=fopen('jason1.json','w');
fprintf(fid, j1);
fclose(fid);
end %if 0 %test

%%
n=length(filenames);
str1=cell(n,1);
for i=1:n
    str1{i}=[filenames{i},num2str(filesize(i))];
end 

%replace forbidden characters in the structure.
str1=strrep(str1,'.','dot');
% clear s1
% append new fields
for i=1:n
    
    field=str1{i};
    zero='zero';
    s1.(field).fileref='';
    s1.(field).size=filesize(i);
    s1.(field).filename=filenames{i};
    s1.(field).base64_img_data='';
    s1.(field).file_attributes=struct();
    
    if 0 %single polygon
        s1.(field).regions.(zero).shape_attributes.name='polygon';
        s1.(field).regions.(zero).shape_attributes.all_points_x=XYb{i}(:,1);
        s1.(field).regions.(zero).shape_attributes.all_points_y=XYb{i}(:,2);
        s1.(field).regions.(zero).region_attributes=struct();
    else %multiple polygons
        XYbi=XYb{i}; %cell 
	if 0 % exist('category_idgc','var')  %Warning: when using flagaug, size mismatch due to new XYbi.
          category_idg=category_idgc{i}.*ones(size(XYbi));
        else % assume only one class.
          category_idg=zeros(length(XYbi),1);
	end
%         regions=struct('shape_attributes',[],'region_attributes',[]);
        for k=1:length(XYbi)
            XYbi_k=XYbi{k};
	    category_id=category_idg(k);
%             s1.(field).regions{k}.shape_attributes.name='polygon';
%             regions(k).shape_attributes.name='polygon';
            s1.(field).regions.(['number',num2str(k-1)]).shape_attributes.name='polygon';
            s1.(field).regions.(['number',num2str(k-1)]).shape_attributes.all_points_x=XYbi_k(:,1);
            s1.(field).regions.(['number',num2str(k-1)]).shape_attributes.all_points_y=XYbi_k(:,2);
            s1.(field).regions.(['number',num2str(k-1)]).shape_attributes.category_id=category_id;
            s1.(field).regions.(['number',num2str(k-1)]).region_attributes=struct();
        end
%         s1.(field).regions=regions;
    end
    
end
j1=jsonencode(s1);

%replace forbidden characters in the structure.
j1=strrep(j1,'dot','.');
j1=strrep(j1,'"zero"','"0"');
j1=strrep(j1,'a2','2');

j1=strrep(j1,'number','');

fid=fopen(ofile,'w');
fprintf(fid, j1);
fclose(fid);

save t1.mat -v7.3

return
end
