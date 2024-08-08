function [category_idsub]=getcateg(S1sub);
%get category id for input shapefile

    %gather catagory ids. 1-4 means stable, slumps, noise, gullies.
    idd=[];category_idsub=zeros(length(S1sub),1);
    for i1=1:length(S1sub)
      if isfield(S1sub(i1), 'zone')
        if contains(S1sub(i1).zone,'scar')
            category_idsub(i1)=2;
        elseif contains(S1sub(i1).zone,'nois') %two classes
            category_idsub(i1)=3;
        elseif contains(S1sub(i1).zone,'gull') %two classes
            category_idsub(i1)=4;
        else %no go / debri, 99
            category_idsub(i1)=99;
%    idd=[idd(:);i1];
        end
      else % no zone; % all scar
        category_idsub(i1)=2;
      end
    end
    S1sub(idd)=[];category_idsub(idd)=[];

return
end

