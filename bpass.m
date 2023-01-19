% 
% Bandpass filter that suppresses pixel noise and long-wavelength variations.
%
% filtered_img = bpass( img, lnoise )
% filtered_img = bpass( img, lnoise, excl_dia )
% filtered_img = bpass( img, lnoise, excl_dia, threshold )
%
% img:          2D array of image pixel values.
%      
% lnoise:       Characteristic lengthscale of noise in pixels.
%
%               Additive noise averaged over this numnber of pixels should vanish. 
%               May be set to 0 or false, in which case only the highpass "background subtraction" 
%               operation is performed.
%
% excl_dia:     Integer length in pixels somewhat larger than a typical object.
%           
%               Can also be set to 0 or false, in which case only the lowpass "blurring"
%               operation defined by lnoise is done, without the background subtraction
%               defined by excl_dia.  Defaults to false.
%
% threshold:    After the convolution of the image, any pixel values below threshold are
%               reset to 0. Defualts to '0' so any negative values are reset.
%
%               Alternatively, set to -Inf so that no thresholding is performed.
%
% returns:      2D array of filtered image pixel values.
%
% Notes:        Performs a bandpass by convolving with an appropriate kernel. You can think of this as 
%               a two part process. First, a lowpassed image is produced by convolving the original 
%               with a gaussian. Next, a second lowpassed image is produced by convolving the original 
%               with a boxcar function. By subtracting the boxcar version from the gaussian version, we
%               are using the boxcar version to perform a highpass.
%
%                   
%

%{

Notes on convolution:

JWM: Do a 2D convolution with the kernels in two steps each. It is
possible to do the convolution in only one step per kernel with 

  gaussian_conv = conv2(gaussian_kernel',gaussian_kernel,image,'same');
  boxcar_conv = conv2(boxcar_kernel', boxcar_kernel,image,'same');

but for some reason, this is slow. The whole operation could be reduced
to a single step using the associative and distributive properties of
convolution:

  filtered = conv2(image,...
    gaussian_kernel'*gaussian_kernel - boxcar_kernel'*boxcar_kernel,...
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

function filtered_image = bpass( img, lnoise, excl_dia, threshold )

    if nargin < 2
        warning('No image filtering performed. Not enough arguments provided in bpass( image, lnoise )')
        filtered_image = img ;
        return
    end

    if nargin < 3, excl_dia = false ; end

    if excl_dia == true && rem( excl_dia, 2 ) == false 
        warning('No image filtering performed. Exclusion diameter (excl_dia) must be an odd integer.') ;
        filtered_image = img ;
        return ;
    end

    normalize = @( x ) x / sum( x ) ;

    if isa( img, 'double' ) ~= 1, img = double( img ) ; end

    img = img - min( img, [], 'all' ) ;
    img = img ./ max( img, [], 'all' ) * 255 ;
    
    if lnoise == 0

        gaussian_kernel = 1 ;

    else

        gaussian_x = - ceil( 5 * lnoise ) : ceil( 5 * lnoise ) ;   
        
        gaussian_kernel = normalize(...
            exp( -( gaussian_x / ( 2 * lnoise ) ) .^2 ) ) ;
        
    end

    gaussian_conv = conv2( img', gaussian_kernel', 'same' );
    gaussian_conv = conv2( gaussian_conv', gaussian_kernel', 'same' ) ;

    if excl_dia

        boxcar_kernel = normalize(...
            ones( 1, length( -round( excl_dia ) : round( excl_dia ) ) ) ) ;

        boxcar_conv = conv2( img', boxcar_kernel', 'same' ) ;
        boxcar_conv = conv2( boxcar_conv', boxcar_kernel', 'same' ) ;

        filtered_image = gaussian_conv - boxcar_conv ;

    else
        
        filtered_image = gaussian_conv ;

    end

    filtered_image = filtered_image - min( filtered_image, [], 'all' ) ;
    filtered_image = filtered_image ./ max( filtered_image, [], 'all' ) * 255 ;

    if nargin == 4
        filtered_image( filtered_image < threshold ) = 0 ;
    end

    