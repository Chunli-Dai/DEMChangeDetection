
%stripdir='/fs/project/howat.4/EarthDEM/region*/strips_unf/2m/'; %Directory of strip DEM files.
stripdir='/home/chunlidai/blue/data/ArcticDEMdata/arcticdem_08_canada_baffin/';
codedir='/home/chunlidai/blue/apps/landslide/code1clean/';  %Directory of codes.

% 1 assume event time known; 0 assume time unknown
eqepoch=0;timefix=0; %Recommend 0 0

flagvolc=0; %Recommend 0. 1 volcano event, a duration of time; 0 landlside, a sudden change.
jumpflagc=1; %Recommend 1. 1: chose to estimte change no matter it's significant or not.
	      %otherwise: the algorithm decide not to estimate the change if it's not significant.
smlarea=200; %Recommend 200. To remove clusters smaller than this; square meters 
%In proposal, it sayes "disturbance areas commonly exceed 20 ha (200 000 m^2)"
% Ashley's polygons: minimum area 1e3 m^2, median area 11e3 m^2, mean area 24e3.

mons=7;mone=8; %mon>=5&mon<=10; %mons, start month of snow-free seasons. %mone, end month of snow free months.
		%suggest to use all seasons; For latitude 80N, use July and August for summer;
		%For latitude 70N, use June to September for summer.
mons=1;mone=12; %all season.

algorithmin=2;%Recommend 2. Fit model for time series. 
	%e.g.  1 linear (ice melting); 2 constant (landslides) ; 3 constant + linear (Okmok volcano)
poolsize=1; % Recommend 1.
resr=10;  %Recommend 2. Output spatial resolution in meters. 2; %4;%40.;

flagrock=0; %Recommend 0. 1 must use rock as control points (for ice melting); 0 it's ok if there are no rock control points (change detection).

flagplot=0; %Recommend 0.

flagfilter=1; %Recommend 1. 1 use the bitmask.tif to map out bad DEMs before change detection;
flagfiltercoreg=0;  %Recommend 0 (Not apply the filtering mask). 
		%1 \use the bitmask.tif to map out bad DEMs for setsm coregistration (coregflat=8);
		% Due to the improvement of coregistration, the cloudy image is not filtered out in the offset re-adjustment step.
	         % Check ~/chunliwork/landslide/alaska3c/site1/5kmoffsetErrMax_filter/ 

coregflag=3;%Recommend 8. 1 parameter (vertical), 3 parameter (Ian's); 7 parameter (MJ's); 8, MJ's setsm coregistration with 3 parameters.
%coregflag=3;% use 3 if we are using scheduler.F90 for bundling jobs. (coregflag=8) setsm is not compatible with scheduler.F90 since setsm use mpi. -> successful tested in arcticdem_09_canada_victoria with coregflag=8!

maxdt=1e9; % Recommend 1e9. in days;  filter out cross-track collected more than maxdt; %Recommend 1e9 (no limit) 
year_start=0;year_end=9999; % min max year of measurement; year >= years and year <= yeare %Recommend no limit [0, 9999]
flagoutput=0; %Recommend 0. if 1 output lots of figures and data for plotting, 0 save space;

maxpxpy=15; maxpz=20; maxsigma=15;% Recommend: maxpxpy=15; maxpz=20;maxsigma=15;ArcticDEM with coregflag=8;
%maxpxpy=15; maxpz=20; maxsigma=10;% Recommend: maxpxpy=15; maxpz=20;maxsigma=10;ArcticDEM;
% maxpxpy=50;maxpz=50; maxsigma=30; % other dems in coreg3.m adjustOffsets.m coregisterdems.m

flaggreenland=0; %Recommend 0. 1 use GIMP masks for greenland; 0 otherwise;
flagfilteradj=1;%Recommend 1. 1 apply filter in adjustOffsets.m, 0 do not apply filter (use all DEMs).
flagmatchtag=1; %Recommend 1. Whether or not to use matchtag as a filter to DEM; 1 apply filter; 0 not apply.

flag_diffdems=0; %Recommend 0. 1 output sequential DEM difference; 0 DO NOT output these.

demext='dem_10m.tif'; %'dem_10m.tif';%'dem.tif';
