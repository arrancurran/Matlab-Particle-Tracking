% out=cntrd(im,est_pks,excl_dia,interactive)
% 
% PURPOSE:  calculates the centroid of bright spots to sub-pixel accuracy.
%  Inspired by Grier & Crocker's feature for IDL, but greatly simplified and optimized
%  for matlab
% 
% INPUT:
% im: image to process, particle should be bright spots on dark background with little noise
%   ofen an bandpass filtered brightfield image or a nice fluorescent image
%
% est_pks: locations of local maxima to pixel-level accuracy from pkfnd.m
%
% excl_dia: diamter of the window over which to average to calculate the centroid.  
%     should be big enough
%     to capture the whole particle but not so big that it captures others.  
%     if initial guess of center (from pkfnd) is far from the centroid, the
%     window will need to be larger than the particle size.  RECCOMMENDED
%     size is the long lengthscale used in bpass plus 2.
%     
%
% interactive:  OPTIONAL INPUT set this variable to one and it will show you the image used to calculate  
%    each centroid, the pixel-level peak and the centroid
%
% NOTE:
%  - if pkfnd, and cntrd return more then one location per particle then
%  you should try to filter your input more carefully.  If you still get
%  more than one peak for particle, use the optional excl_dia parameter in pkfnd
%  - If you want sub-pixel accuracy, you need to have a lot of pixels in your window (excl_dia>>1). 
%    To check for pixel bias, plot a histogram of the fractional parts of the resulting locations
%  - It is HIGHLY recommended to run in interactive mode to adjust the parameters before you
%    analyze a bunch of images.
%
% OUTPUT:  a N x 4 array containing, x, y and brightness for each feature
%           out(:,1) is the x-coordinates
%           out(:,2) is the y-coordinates
%           out(:,3) is the brightnesses
%           out(:,4) is the sqare of the radius of gyration

%{

CHANGELOG:

Feb 4 2005
Written by Eric R. Dufresne, Yale University.

May 2005
Inputs diamter instead of radius

Jun 2005
Added code from imdist/dist to make this stand alone. DB

Increased frame of reject locations around edge to 1.5*excl_dia ERD

By popular demand, 
1. altered input to be formatted in x,y space instead of row, column space. ERD
2. added forth column of output, rg^2. ERD

Aug 2005
Outputs had been shifted by [0.5,0.5] pixels.  No more! ERD

Aug 24 2005
Woops!  That last one was a red herring.  The real problem is the "ringing" from the output of bpass.
I fixed bpass (see note), and no longer need this kludge. Also, made it quite nice if est_pks=[]; ERD

Jun 2006
Added size and brightness output ot interactive mode. Also fixed bug in calculation of rg^2. ERD

Jun 2007
Small corrections to documentation. JWM

Jan 2023
Reformated to meet commenting and nomenclecture standards. AC
Removed interactive option. AC
Changed excl_rad such that we consider n pixels from the given peak where 2n + 1 is the input excl_dia.
Removed filtering image edges of peaks and this is already happening in pkfnd(). AC
Changed radius of giration calc. I dont trust rg = ( sum( tmp .* dst2, 'all' ) / norm ) ; since the applied mask influences
the estimated radius.
Changed the mask such that it is a pixellated circle of diameter = excl_dia rather than excl_dia - 1.




    % Create mask - window around trial location over which to calculate the centroid

    % msk_range = ( - excl_rad : excl_rad ) .^2 ;
    
    % cent_px = excl_rad + 1 ;

    % circ_msk_inv = zeros( excl_dia ) ;
    
    % for i = 1 : excl_dia
    %     circ_msk_inv( i, : ) = sqrt( ( i - cent_px ) ^2 + msk_range ) ;
    % end

    % ind = find( circ_msk_inv <= excl_rad ) ;

    % circ_msk_binary = zeros( excl_dia ) ;

    % circ_msk_binary( ind ) = 1.0 ;

    % dst2 = circ_msk_binary .* ( circ_msk_inv .^2 ) ;

%}

function particles = cntrd( img, est_pks, excl_dia, apply_mask )

    if rem( excl_dia, 2 ) == false 
        warning('Exclusion diameter (excl_dia) must be an odd integer.') ;
        out = est_pks ;
        return ;
    end

    if isempty( est_pks )
        warning('There were no estimated peaks (est_pks) provided. Maybe the threshold in pkfnd() is too high.')
        out = [ ] ;
        return;
    end

    % Number of pixels from the central pixel.
    excl_rad = floor( excl_dia / 2 ) ;

    if apply_mask == true

        cent_px = excl_rad + 1 ;
        % Create mask - window around trial location over which to calculate the centroid
        circ_msk_binary = zeros( excl_dia ) ;  
        circ_msk_binary( cent_px, cent_px ) = 1 ;    
        circ_msk_binary = bwdist( circ_msk_binary );
        circ_msk_binary = circ_msk_binary <= excl_rad;
    
    else

        circ_msk_binary = 1 ;
    end

    msk_ind_x = zeros( excl_dia ) ;

    for n = 1 : excl_dia, msk_ind_x( n, : ) = ( 0 : excl_dia - 1 ) ; end
    
    msk_ind_y = msk_ind_x' ;

    [ est_pks_num, ~ ] = size( est_pks ) ;

    particles = zeros( est_pks_num , 4 ) ;

    % Loop through all of the candidate positions
    for n = 1 : est_pks_num
        
        tmp = circ_msk_binary .* ... 
                img( ( est_pks( n, 2 ) - excl_rad  : est_pks( n, 2 ) + excl_rad ), ...
                     ( est_pks( n, 1 ) - excl_rad  : est_pks( n, 1 ) + excl_rad ) ) ;

        tot_br = sum( tmp, 'all' ) ;

        wght_ave_x = sum( tmp .* msk_ind_x, 'all' ) / tot_br - excl_rad ;

        wght_ave_y = sum( tmp .* msk_ind_y, 'all' ) / tot_br - excl_rad ;

        % Calculate an estimate of the particles diameter based on radius of gyration
        rad_gyr = 2 * sqrt( sum( (tmp) .^2 , 'all' ) / numel( tmp ) / 255 ) ;

        pk_val = max( tmp, [],'all' ) ;
        
        particles( n, : ) = [ est_pks( n, 1 ) + wght_ave_x, est_pks( n, 2 ) + wght_ave_y, pk_val, rad_gyr ] ;
                
    end
