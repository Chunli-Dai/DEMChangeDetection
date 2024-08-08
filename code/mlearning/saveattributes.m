function [filenames,filesize,XYb,absfilenames]=saveattributes(filenames,filesize,XYb,absfilenames,countaug,dirfilename_aug,XYbi_aug)

        [~,fname,ext]=fileparts(dirfilename_aug);
        filenamei=[fname,ext];
        
        % get file attributes
        str1=imfinfo(dirfilename_aug);
        filesizei=str1.FileSize;

        filenames{countaug}=filenamei;
        filesize(countaug)=filesizei;
        XYb{countaug}=XYbi_aug;
        absfilenames{countaug}=dirfilename_aug;
end