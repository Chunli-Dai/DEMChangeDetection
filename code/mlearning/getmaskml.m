function [maskml,Mpcg_ml2]=getmaskml(data0,Mpcg,ns,ids,maskdir2,flagold)
% reconstruct mask from machine learning 2.
% flagold, 1, old name convention; 0 new name convention.

resr2m=2;
mask=false(size(data0.z)); %false
maskx=data0.x;masky=data0.y;
Mpcg_ml2=cell(ns,1);
for ixy=1:ns %length(ids) %13%
%             ids{ix}=[idxs idxe idys idye];
ix=ixy;
id1=ids{ixy};
Mx=false(size(data0.x));My=false(size(data0.y));
Mx(id1(1):id1(2))=1; My(id1(3):id1(4))=1;

maskiold=mask(My,Mx); %Mx should be the same as that in Mpcg{ixy}.

% ofile=['slumppeelsubmask/slumppeelsubmaski',num2str(ixy),'.mat'];
if flagold==1
ofile=[maskdir2,'/slumppeelsubmaski',num2str(ixy),'.mat']; %old
else
ofile=[maskdir2,'/slumppeelsubi',num2str(ixy),'.mat'];  %new
end
if ~exist(ofile,'file')
%   fprintf(['\n ',ofile,'does not exist.']) ; % no slumps detected
    continue
end
A1=load(ofile);

A1=A1.arr;
% A1 should have the same size as Mpcg{ixy}.z;

%size matching, see box2mosaic.m: resr is 2m for maxlen<1km, 10m otherwise.
datai.x=data0.x(Mx); datai.y=data0.y(My);%datai.z=data0.z(My,Mx);
rangi=[min(datai.x) max(datai.x) min(datai.y) max(datai.y)]; %
maxlen=max([abs(rangi(2)-rangi(1))*1e-3, abs(rangi(4)-rangi(3))*1e-3]);% maximum length in km
if maxlen <= 1 % maximum side of the box is <=4 km; small tiles
%A1 2m to 10 m
tx=rangi(1):resr2m:rangi(2);ty=rangi(4):-resr2m:rangi(3);
A1 =interp2(tx,ty,double(A1),datai.x,datai.y','*linear',nan);%toc;%4s
end

[n1,m1]=size(maskiold);
[n2,m2]=size(A1);
if n1~=n2 ||m1~=m2
    fprintf([ofile,'matrix size is not consistent.\n']) 
else
    %rough match of mask, may introduce shifts;
%     A1=imresize(A1,size(maskiold));

    mask(My,Mx)=maskiold|A1;
end

Mpc=Mpcg{ix};
Mpc.z=A1; %or the accumulated mask: maskiold|A1.
Mpcg_ml2{ix}=Mpc; %update the mask from ML to Mpcg.
end %ixy

maskml.x=maskx;maskml.y=masky;maskml.z=mask;

return
end
