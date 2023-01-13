set( 0 , 'DefaultFigureWindowStyle' , 'docked' ) ; clear all ; clc ; %close all ;

pairs = 999 ;
ind = 1 ; pair = 0 ;
position_lst = zeros( 1 , 3 ) ;
vel = zeros( 1e3 , 6 );
% X | Y | Vx | Vy | frm pair | prt no. must be larger than expected final vel - guess ...!?

thresh = 10 ;    % Threshold for bpass and pkfnd
lobject = 4 ;   % estimated sixe of colloid in pixels

for n = 1 : pairs ;
    clc ;
    pair = pair + 1
    
    img_0 = ( ( imread( [ '/Users/arrancurran/Documents/Oxford/Flow/22-10-13-3-27PM/ChanA_0001_0001_0001_' num2str( n , '%04i' ) ] , 'tif') ) ) ;
    img_1 = ( ( imread( [ '/Users/arrancurran/Documents/Oxford/Flow/22-10-13-3-27PM/ChanA_0001_0001_0001_' num2str( n + 1 , '%04i' ) ] , 'tif') ) ) ;
    
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
    position_lst = [ [ pk_cntrd_0( : , [ 1 2 ] ) pair * rot90( ones( 1 , length( pk_cntrd_0 ) ) ) ] ; [ pk_cntrd_1( : , [ 1 2 ] ) ( pair + 1 ) * rot90( ones( 1 , length( pk_cntrd_1 ) ) ) ] ] ;
    
    tr_lst = track( position_lst , 6 ) ;
    % X Y frame ID
    
    for v = 1 : length( tr_lst ) - 1 ;
        if tr_lst( v , 4 ) == tr_lst( v + 1 , 4 )
            vel( ind , : ) = [ [ tr_lst( v , 1 : 2 ) tr_lst( v + 1 , 1 : 2 ) - tr_lst( v , 1 : 2 ) tr_lst( v , 3 : 4 ) ] ] ;
            ind = ind + 1 ;
        end
    end
    %{
    subplot( 1 , 2 , 1 )
    colormap( 'gray' ) , imagesc( img_bpass_0 ) ;
    subplot( 1 , 2 , 2 )
    contourf( img_bpass_0 , 'DisplayName' , 'img_bpass_0' , 'LineStyle' , 'none' ) ; colormap( 'gray' )
    hold on ; plot( pk_cntrd_0( : , 1 ) , pk_cntrd_0( : , 2 ) , 'og' ) ; hold off
    quiver( vel( : , 1 ) , vel( : , 2 ) , vel( : , 3 ) , vel( : , 4 ) ) ;
    %drawnow
    %}
end

ind = find( vel( : , 1) == 0 & vel( : , 2 ) == 0 ) ;
vel( ind , : ) = [] ;
subImgSz = 12 ; % Spacing of vectors in px
block = size( img_0 ) / subImgSz ;
vel_block_Vx = zeros( block ) ;
vel_block_Vy = zeros( block ) ;
vel_temp = vel ;
for n = 1 : block( 1 ) ;
    for p = 1 : block( 2 ) ;
        ind = find( ...
            vel_temp( : , 1 ) >= ( n - 1 ) * subImgSz + 1  & ...
            vel_temp( : , 1 ) < n * subImgSz & ...
            vel_temp( : , 2 ) >= ( p - 1 ) * subImgSz + 1 & ...
            vel_temp( : , 2 ) < p * subImgSz ) ;
        vel_block_Vx( n , p ) = sum( vel_temp( ind , 3 ) ) ;
        vel_block_Vy( n , p ) = sum( vel_temp( ind , 4 ) ) ;
        vel_temp( ind , : ) =  [ ] ;
    end
end
[ X , Y ] = meshgrid( 1 : subImgSz : block( 1 ) * subImgSz ) ;
quiver( X , Y , -vel_block_Vx , -vel_block_Vy ) ;