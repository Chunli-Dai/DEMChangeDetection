function [co]= createmeta(infile)
%% Construct meta.txt file from DEM file
%see ../codetest/testotherdems.m 

%infile='20100913124829_SPI_11004_Iceland_V2_dem.tif';
co=[];

[XYbi,rangei]=dembd(infile);

[demdir,name,ext] =fileparts([strtrim(infile)]);
filename=[name,ext];

ndx=20;%20*8;
% figure;hold on;plot(XYbi(1:ndx:end,1),XYbi(1:ndx:end,2),'k.')
metafile=strrep(infile,'dem.tif','meta.txt');
fid1=fopen(metafile,'w');
fprintf(fid1,'%s\n', ['X:',num2str(round(XYbi(1:ndx:end,1)'))]);
fprintf(fid1,'%s\n', ['Y:',num2str(round(XYbi(1:ndx:end,2)'))]);
fprintf(fid1,'%s\n', ['scene 1 name=',infile]);
fprintf(fid1,'%s\n', ' ');
fprintf(fid1,'%s\n', ' ');
fprintf(fid1,'%s\n', ['Image 1=/WV02_',filename]);
fprintf(fid1,'%s\n', ['Image 2=/WV02_',filename]);
fclose(fid1);

return
end
