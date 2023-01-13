%Origanally used to plot centroids of two subsequent images used in
%generating the flow fields for Optica Supp. Mat.

thresh = 10 ;    % Threshold for bpass and pkfnd
lobject = 8 ;   % estimated size of colloid in pixels

fileID = fopen( ...
    ['/Users/arrancurran/Desktop/Gr test/'...
    'Gr test 2048 x 2048 x 1000.raw' ] ) ;
    
dat_raw = fread( fileID ) ;

dat_raw = reshape( dat_raw , 2048 , 2048 , 100 ) ;

img_0 = dat_raw( : , : , 1 ) ; img_0 = reshape( img_0 , 2048 , 2048 ) ;
img_1 = dat_raw( : , : , 2 ) ; img_1 = reshape( img_1 , 2048 , 2048 ) ;

img_bpass_0 = bpass( img_0 , 1 , lobject , thresh ) ; img_bpass_1 = bpass( img_1 , 1 , lobject , thresh ) ;
%{
    res = bpass( image_array, lnoise, lobject )
    INPUTS:
    lnoise:     Characteristic lengthscale of noise in pixels. Additive noise averaged over this length should
                vanish. May assume any positive floating value. May be set to 0 or false, in which case only the
                highpass "background subtraction" operation is performed.
    lobject:    (optional) Integer length in pixels somewhat larger than a typical object. Can also be set to
                0 or false, in which case only the lowpass "blurring" operation defined by lnoise is done,
                without the background subtraction defined by lobject.  Defaults to false.
    threshold:  (optional) By default, after the convolution, any negative pixels are reset to 0.  Threshold
                changes the threshhold for setting pixels to 0.  Positive values may be useful for removing
                stray noise or small particles.  Alternatively, can be set to -Inf so that no threshholding is
                performed at all.
%}

pk_0 = pkfnd( img_bpass_0 , thresh , lobject ) ; pk_1 = pkfnd( img_bpass_1 , thresh , lobject ) ;
%{
    res = pkfnd( image_array, thresh, lobject )
    INPUTS:
    threshold:  the minimum brightness of a pixel that might be local maxima. (NOTE: Make it big and the code runs
                faster but you might miss some particles.  Make it small and you'll get everything and it'll be slow.)
    lobject:    if your data's noisy, (e.g. a single particle has multiple local maxima), then set this optional
                keyword to a value slightly larger than the diameter of your blob.  if multiple peaks are found
                withing a radius of lobject/2 then the code will keep only the brightest.  Also gets rid of all
                peaks within lobject of boundary
    OUTPUT:     N x 2 array containing, [row,column] coordinates of local maxima
                out(:,1) are the x-coordinates of the maxima
                out(:,2) are the y-coordinates of the maxima
%}
pk_cntrd_0 = cntrd( img_bpass_0 , pk_0 , lobject + 2 ) ; pk_cntrd_1 = cntrd( img_bpass_1 , pk_1 , lobject + 2 ) ;

%plot( pk_cntrd_0( : , 3 ) , pk_cntrd_0( : , 4 ) , '.b' , pk_cntrd_1( : , 3 ) , pk_cntrd_1( : , 4 ) , '.r' )
%drawnow
% X Y frame
%position_lst = [ [ pk_cntrd_0( : , [ 1 2 ] ) pair * rot90( ones( 1 , length( pk_cntrd_0 ) ) ) ] ; [ pk_cntrd_1( : , [ 1 2 ] ) ( pair + 1 ) * rot90( ones( 1 , length( pk_cntrd_1 ) ) ) ] ] ;

%tr_lst = track( position_lst , 6 ) ;
% X Y frame ID

%for v = 1 : length( tr_lst ) - 1 ;
%    if tr_lst( v , 4 ) == tr_lst( v + 1 , 4 )
%        vel( ind , : ) = [ [ tr_lst( v , 1 : 2 ) tr_lst( v + 1 , 1 : 2 ) - tr_lst( v , 1 : 2 ) tr_lst( v , 3 : 4 ) ] ] ;
%        ind = ind + 1 ;
%    end
%end

subplot( 1 , 2 , 1 )
colormap( 'gray' ) , imagesc( img_bpass_0 ) ;
subplot( 1 , 2 , 2 )
contourf( img_bpass_0 , 'DisplayName' , 'img_bpass_0' , 'LineStyle' , 'none' ) ; colormap( 'gray' )
hold on ; plot( pk_cntrd_0( : , 1 ) , pk_cntrd_0( : , 2 ) , 'og' ) ; hold off
%drawnow
