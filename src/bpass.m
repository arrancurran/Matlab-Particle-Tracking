% 
% Three step filtering starting with a boxcar filter to remove long scale variation. 
% Boxcar filtered image, 'img_box' is then filtered with a gausian filter, to remove pixel noise.
% Finally, 'img_gaus', has any pixel values below 'baseline' set to zero.
% Any step can be skipped with a 'false' argument.
%
% img_out = bpass( img, box_filter, gaus_filter, baseline )
% img_out = bpass( img, box_filter, gaus_filter, baseline, display )
%
% img:          'array' 2D array of image pixel values.
%
% box_filter:   'bool' Set to 'true' for highpass filtering.
%               Set to 'false' to skip boxcar filtering of the input image.
%
% gaus_filter: 'int|bool' Characteristic length scale of noise in pixels.
%               Set to 'true' to apply an appropriate gausian filter based on the input image.
%               Or provide any positive integer for manual control of the gausian kernel. 
%               Set to 'false' to skip gausian filtering of the input image.
%
% baseline:     'int|bool' Reset any pixel values below 'baseline' to 0.
%               Set to 'false' to skip. 
%               An input of '0' will output the same result as 'false', but the code will scan the image for any values below 0.
%
% returns:      'array' 2D array of filtered image pixel values.
%
% Notes:        Performs a bandpass by convolving with an appropriate kernel. You can think of this as 
%               a two part process. First, a lowpassed image is produced by convolving the original 
%               with a gausian. Next, a second lowpassed image is produced by convolving the original 
%               with a boxcar function. By subtracting the boxcar version from the gausian version, we
%               are using the boxcar version to perform a highpass.

%{

Notes on convolution:

JWM: Do a 2D convolution with the kernels in two steps each. It is
possible to do the convolution in only one step per kernel with 

  gaus_conv = conv2(gaus_kernel',gaus_kernel,image,'same');
  boxcar_conv = conv2(box_kernel', box_kernel,image,'same');

but for some reason, this is slow. The whole operation could be reduced
to a single step using the associative and distributive properties of
convolution:

  filtered = conv2(image,...
    gaus_kernel'*gaus_kernel - box_kernel'*box_kernel,...
    'same');

But this is also comparatively slow (though inexplicably faster than the
above). It turns out that convolving with a column vector is faster than
convolving with a row vector, so instead of transposing the kernel, the
image is transposed twice.

CHANGELOG:

Feb 1993
Written by David G. Grier, The University of Chicago.

May 1995
Greatly revised version DGG.

Dec 1995
Added /field keyword JCC.

Aug 1999
Memory optimizations and fixed normalization, DGG.

Apr 2004-ish
Converted to Matlab by D.Blair.

June 2005
Fixed some bugs with conv2 to make sure the edges are removed D.B.
Removed inadvertent image shift ERD.

Aug 24 2005
Added threshold to output. Now sets all pixels with negative values equal to zero. Gets rid of ringing 
which was destroying sub-pixel accuracy, unless window size in cntrd was picked perfectly. Now centrd 
gets sub-pixel accuracy much more robustly ERD.

Jun 2007
Refactored for clarity and converted all convolutions to use column vector kernels for speed. Running
on my  macbook, the old version took ~1.3 seconds to do bpass(image,1,19) on a 1024 x 1024 image;
this version takes roughly half that. JWM

Jan 2023
Reformated to meet commenting and nomenclecture standards. AC
Added image scalling on input and ouptut so we process with full bandwidth (0-255) and we return full bandwidth. AC
Added check to see if image in has already been converted to 'double'. AC
Removed edge zero-ing as this happens in pkfnd() by removing peaks within edge rather than blacking out and potentially 
picking parts of particles.

%}

function img_out = bpass( img, box_filter, gaus_filter, baseline, display )

    if nargin < 4
        warning('No image filtering performed. Not enough arguments provided in bpass( img, box_filter, gaus_filter, baseline, display )')
        img_out = img ;
        return
    end

    if ~exist( 'display', 'var' ), display = false ; end

    if isa( img, 'double' ) ~= 1, img = double( img ) ; end

    normalize   = @( x ) x / sum( x ) ;
    scale2init8 = @( x ) ( x - min( x, [], 'all' ) ) ./ max( ( x - min( x, [], 'all' ) ), [], 'all' ) * 255 ;
    img         = scale2init8( img ) ;
    img_out     = img ;

    if box_filter
        box_kernel  = - ones( 3 ) / 9 ; box_kernel( 2, 2 ) = 8 / 9 ;
        img_box     = conv2( img_out, box_kernel, 'same' ) ;
        img_box     = scale2init8( img_box ) ;
        img_out     = img_box ;
    end

    if gaus_filter ~= false

        if gaus_filter == true
            % Fast Noise Variance Estimation, see https://doi.org/10.1006/cviu.1996.0060
            [ img_rows, img_cols ] = size( img ) ;

            gaus_sigma     = sum( abs( conv2( img, [ 1 -2 1 ; -2 4 -2 ; 1 -2 1 ] ) ), 'all'  ) ;
            gaus_filter    = round( gaus_sigma * sqrt( .5 * pi ) / ( 6 * ( img_rows - 2 ) * ( img_cols - 2 ) ) ) ;
        end

        if gaus_filter ~= 1 % Dont waste my time with good images!
            gaus_x      = - gaus_filter : gaus_filter ;
            gaus_kernel = normalize( exp( -( gaus_x / ( 2 * gaus_filter ) ) .^2 ) ) ;
            img_gaus    = conv2( img_out, gaus_kernel, 'same' ) ;
            img_gaus    = conv2( img_gaus, gaus_kernel', 'same' ) ;
            img_gaus    = scale2init8( img_gaus ) ;
            img_out     = img_gaus ;
        end
        
    end

    if baseline
        img_base = img_out ;
        img_base( img_base < baseline ) = 0 ; 
        img_out = img_base ;
    end

    if display == true

        fov = 512 ;
        figure_img = figure ; colormap( figure_img, 'gray') ; figure_hists = figure ;

        img_hist = @( x )  hist( x, min( x, [], 'all' ) : max( x, [], 'all' ) ) ;

        [ hist_raw, x_hist ] = img_hist( img ) ;
        figure_hists ; semilogy1_raw = semilogy( x_hist, sum( hist_raw, 2 ), 'ko' ) ;
        set(semilogy1_raw, 'DisplayName', 'Raw' ) ;
        hold on
        display_raw = subplot( 2, 2, 1, 'Parent', figure_img ) ; image( img( 1 : fov, 1 : fov ), 'Parent', display_raw) ;
        title( display_raw, 'Raw Image' ) ; set( display_raw, 'YTickLabel', [ ] ) ; set( display_raw, 'XTickLabel', [ ] ) ;

        if box_filter
            [ hist_box, x_hist ] = img_hist( img_box ) ;
            figure_hists ; semilogy1_box = semilogy( x_hist, sum( hist_box, 2 ), 'bo' ) ;
            set(semilogy1_box, 'DisplayName', 'Boxcar', 'MarkerFaceColor', 'b' ) ;
            display_box = subplot( 2, 2, 2, 'Parent', figure_img ) ; imagesc( img_box( 1 : fov, 1 : fov ), 'Parent', display_box) ;
            title( display_box, 'Boxcar Filtered Image' ) ;set( display_box, 'YTickLabel', [ ] ) ; set( display_box, 'XTickLabel', [ ] ) ;
        end

        if exist( 'img_gaus', 'var' )
            [ hist_g, x_hist ] = img_hist( img_gaus ) ;
            figure_hists ; semilogy1_g = semilogy( x_hist, sum( hist_g, 2 ), 'ro' ) ;
            set( semilogy1_g, 'DisplayName', 'Gausian', 'MarkerFaceColor', 'r' ) ;
            display_gaus = subplot( 2, 2, 3, 'Parent', figure_img ) ; image( img_gaus( 1 : fov, 1 : fov ), 'Parent', display_gaus) ;
            title( display_gaus, 'Gaussian Filtered Image' ) ;set( display_gaus, 'YTickLabel', [ ] ) ; set( display_gaus, 'XTickLabel', [ ] ) ;

        end

        [ hist_f, x_hist ] = img_hist( img_out ) ;
        figure_hists ; semilogy1_f = semilogy( x_hist, sum( hist_f, 2 ), 'go' ) ;
        set(semilogy1_f, 'DisplayName', 'Output' ) ;
        display_out = subplot( 2, 2, 4, 'Parent', figure_img ) ; image( img_out( 1 : fov, 1 : fov ), 'Parent', display_out) ;
        title( display_out, 'Output Image' ) ;set( display_out, 'YTickLabel', [ ] ) ; set( display_out, 'XTickLabel', [ ] ) ;

    end