function [co]=modifyjson(jsonfile,badfiles)
%Remove a list of filenames and attributes from a given json file
% jsonfile: given json file
% badfiles: a list of file information to be deleted.

co=[];

    if exist(jsonfile,'file')
    fid = fopen(jsonfile); 
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid); 
    else;
	fprintf(['\n json file does not exist: ',jsonfile,' .\n']) 
	str=[];
	return
    end

    % connect the start and end.
    if ~isempty(str)
%	str(1)='';

	ids =strfind(str,'"slump');
	ids=ids(1:2:end); %start id of each file
	ide =strfind(str,'}}}}');
	ide(end)=[];     %end id of each file
	ide=ide+3; 
	if length(ids)~=length(ide)
	   fprintf(['\n The index of start and end ids of each file is wrong! \n'])
	end
	nfile=length(ids); %number of all files

	idg_del=[];
	for i=1:length(badfiles)
	    fprintf(['\n Working on ',num2str(i),' out of ',num2str(length(badfiles)),' bad files. \n'])
	    filei=badfiles{i};
	    id=strfind(str,filei);
	    id=id-1; %index of "slump2g1.jpg", including the quote symbol.

        if ~isempty(id)
	    idf=find(ids==id(1));  %order number of the target file
	    id_del=[];
	    if idf~=nfile
     	        id_del=ids(idf):(ide(idf)+1); %plus 1 to include the comma ,
	    else  %the last file
	        id_del=(ids(idf)-1):(ide(idf)); %minus 1 to include the comma , before
	    end
	    idg_del=[idg_del(:);id_del(:)];
        else % 
            fprintf(['\n File information of ',filei,' is not found in json file: ',jsonfile,'.\n'])
        end % empty id
	end %for i

	str(idg_del)=[];

	% fix bug 1: in case of removing last two files, it yields extra comma: e.g. }}}},}
	if strcmp(str(end-5:end),'}}}},}')
		str(end-1)=''; %remove this extra comma
	end

	strg_m=str;
    else
	fprintf(['\n json file is empty: ',jsonfile,' .\n'])
        str=[];
        return
    end

    %val = jsondecode(str); %string to structure; some characters are modified, e.g. 0, dot.
    % modify val
    %strg_m=jsonencode(val);  

    jsonfile_bp=strrep(jsonfile,'.json','_bp1.json');  
    [status, cmdout]=system(['mv ',jsonfile, ' ',jsonfile_bp]);

    %write modified json files
    ofile=jsonfile;
    fid=fopen(ofile,'w');
    fprintf(fid, strg_m);
    fclose(fid);

return
end
