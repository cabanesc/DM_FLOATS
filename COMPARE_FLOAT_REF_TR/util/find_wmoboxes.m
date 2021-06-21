
% function [ pa_wmo_numbers ] = find_wmoboxes( region, pa_wmo_boxes);
%
% This function finds  WMO boxes corresponding to the zone  region
% The WMO box numbers, between 90N and 90S, are stored in wmo_boxes.mat.
% The 1st column has the box numbers, the 2nd column denotes CTD data,
% the 3rd column denotes bottle data, the 4th column denotes Argo data.
% No data is denoted by 0. Otherwise 1.
%
% A. Wong, 16 August 2004
%
% C. Cabanes Nov. 2014 : extend la_x, so interp2 does not think longitudes in the range [5W 5E] are out-of-bound with matlab version >= R2012b

function [ pa_wmo_numbers ] = find_wmoboxes( region, pa_wmo_boxes);

pa_wmo_numbers = [ NaN.*ones( 25, 1 ), zeros( 25, 1 ), zeros( 25, 1 ), zeros( 25, 1 ) ] ;

la_lookup_x = [ ] ;
la_lookup_y = [ ] ;
vector_x = [] ;
vector_y = [] ;
%keyboard
la_x = [ -5:10:365 ] ; % 38 elements
for i=1:18
  la_lookup_x = [ la_lookup_x; la_x ] ;
end

la_y = [ 85:-10:-85 ] ; % 18 elements
for i=1:38
  la_lookup_y = [ la_lookup_y, la_y' ];
  vector_y = [ vector_y; la_y' ];
  vector_x = [ vector_x; la_x(i).*ones(18,1) ];
end

la_lookup_no = reshape( [ 1:648 ], 18, 36 ) ;
la_lookup_no=[la_lookup_no(:,end),la_lookup_no,la_lookup_no(:,1)];

ln_x=[region.lonmin:1:region.lonmax];
ln_x(ln_x<0)=ln_x(ln_x<0)+360;
ln_y=[region.latmin:1:region.latmax];
ik=0;
clear pa_wmo_numbers
for ix=1:length(ln_x)
    for iy=1:length(ln_y)
        ln_i = interp2( la_lookup_x, la_lookup_y, la_lookup_no, ln_x(ix), ln_y(iy), 'nearest' ) ;
        if( isnan(ln_i)==0 )
        ik=ik+1;
        pa_wmo_numbers(ik)=pa_wmo_boxes(ln_i,1);
        end
    end
end
pa_wmo_numbers=unique(pa_wmo_numbers);


