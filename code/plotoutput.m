
	addpath(genpath(['/home/dai.56/arcticdemapp/landslide/code1/']));

	if 0
	%step1: prepare DEM mosaics
	rang0=[500000 700000 -1400000 -1200000];
	rang0=[500000-200e3 700000 -1400000-200e3 -1200000+200e3];
	rang0=[500000-100e3 700000 -1400000-100e3 -1200000+100e3];
	[data0]=getmosaic(rang0);

%	[status, cmdout]=system('psbasemap');

	%step2: convert cooridnates to longitude latitude. Attention: check no data values. -dstnodata nan
 	%gdalwarp $line  20000211000000_SRTM_utm_$line -t_srs epsg:32645 -tr 30 30
%	[status, cmdout]=system('gdalwarp givensite_rate.tif givensite_rate_lat.tif -t_srs epsg:4326');
	%convert no data value to nan
	%system('gdal_calc.py -A givensite_rate_lat.tif --outfile=t1.tif --NoDataValue=NaN --calc="(A)"');
	end %if 0

	%step3: plot GMT figures

	%gdal tools see ~/arcticdemapp/landslide/codetest/gdaltools.txt
	if 0
	str=['./plotelevd.gmt demmosaic_lat.tif givensite_rate_lat.tif '];
	str='./plotelevd.gmt demmosaic_lat.tif givensite_rate_lat.tif -R-25.3466944/-14.7424460/75.6233940/78.0411836 1.2';
	str='./plotelevd.gmt demmosaic_lat.tif merge_rate_lat.tif -R-25.6/-19/75.6/78. 1.2';

	str='./plotelevd.gmt gimpdem_300m_lat.tif greenland_rateNaN_lat.tif -R-25.6/-19/75.6/78. 1.2';
	fprintf(['\n ',str, '\n'])
%	[status, cmdout]=system([str])
	end

	%entire greenland
	system('rm boxfile.txt')
	str='./plotelevd_greenland.gmt gimpdem_300m_lat.tif greenland_rateNaN_filter1_filledIan_lat.tif -R-54./59/8/81r 2';
         [status, cmdout]=system([str]) 
	system('rm boxfile.txt')
	str='./plotelevd_greenland.gmt gimpdem_300m_lat.tif corrections_lat.tif -R-54./59/8/81r 0.4';
         [status, cmdout]=system([str]) 
exit

	str='./plotelevd.gmt gimpdem_300m_lat.tif allnov_lat.tif -R-90/8/58/84 100';
	str='./plotelevd_greenland.gmt gimpdem_300m_lat.tif allnov_lat.tif -R-90/58/8/84r 100'; %bp2 bad
	str='./plotelevd_greenland.gmt gimpdem_300m_lat.tif allnov_lat.tif -R-90/8/58/84 100';
	str='./plotelevd_greenland.gmt gimpdem_300m_lat.tif allnov_lat.tif -R-55.7/58/8/81r 100';
	str='./plotelevd_greenland.gmt gimpdem_300m_lat.tif allnov_lat.tif -R-54./59/8/81r 100'; % ^_^
	str='./plotelevd_greenland.gmt gimpdem_300m_lat.tif greenland_rateNaN_lat.tif -R-54./59/8/81r 2';
	system('ln -fs boxfile_ne.txt boxfile.txt')
        [status, cmdout]=system([str]) 
	str='./plotelevd_greenland.gmt gimpdem_300m_lat.tif greenland_rateNaN_lat_filter1.tif -R-54./59/8/81r 2';
%	str='./plotelevd_greenland.gmt gimpdem_300m_lat.tif greenland_rateNaN_lat_filter1_fill.tif -R-54./59/8/81r 2';
	%	str='./plotelevd.gmt gimpdem_300m_lat.tif allnov_lat.tif ';
	fprintf(['\n ',str, '\n'])
        [status, cmdout]=system([str]) 

	str='./plotelevd_greenland.gmt gimpdem_300m_lat.tif greenland_ratestdNaN_lat.tif -R-54./59/8/81r 0.4';
        [status, cmdout]=system([str]) 

	str='./plotelevd_greenland.gmt gimpdem_300m_lat.tif greenland_Tave_lat.tif -R-54./59/8/81r 20210101';
	system('rm boxfile.txt')
         [status, cmdout]=system([str]) 
	str='./plotelevd_greenland.gmt gimpdem_300m_lat.tif CountStripwRockNaNnov_lat.tif -R-54./59/8/81r 40'; % ^_^
	system('rm boxfile.txt')
         [status, cmdout]=system([str]) 

