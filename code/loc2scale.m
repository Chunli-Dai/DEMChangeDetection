%given location (latitude), scale length, get the coordinates of start and
%ent points of the scale bar

%lon lat of the label of the scale bar
lat=(79+59/60+30/3600);lon=-(85+54/60);
len=200 ;%length of scale bar in meter

% lat=79.9967;lon=-(85+54/60+15/3600);
% len=50;

%constant parameters
R=6371e3; %Earth radius in meter

dlon=(len/(R*cos(lat*pi/180)))*180/pi; %degrees; the degree difference corresponding to the length.
lon1=lon-dlon/2;
lon2=lon+dlon/2;
dlat=1/4*len/R*180/pi; % scale bar would be 1/4*len below the label.
lats=lat-dlat; %latitude of scale bar; 

fprintf(['\n Location of the label: ',num2str([lon lat])])
fprintf(['\n Location of the scale bar:',num2str([lon1 lats lon2 lats])])

