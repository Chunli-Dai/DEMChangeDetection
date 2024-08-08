addpath(genpath(['/home/dai.56/arcticdemapp/landslide/code1/']));

load sv1.mat
projstr='polar stereo north';
OutName='jump.tif';
writeGeotiff(OutName,xout,yout,double(jump),5,0,projstr)

OutName='jumpstd.tif';
writeGeotiff(OutName,xout,yout,double(jumpstd),5,0,projstr)

OutName='eventtime.tif';
writeGeotiff(OutName,xout,yout,double(timec),5,0,projstr)

OutName='eventtimestd.tif';
writeGeotiff(OutName,xout,yout,double(timestd),5,0,projstr)
