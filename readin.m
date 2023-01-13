clear all ; clc ; %close all ;

thresh = 0 ;    % Threshold for bpass and pkfnd
lobject = 5 ;   % estimated size of colloid in pixels
tic

img_stack = TIFFStack( 'hellma_nofloat_16x16_CROPPED_CONTRAST_1000fps.tif' ) ;

[ i , j , frames ] = size( img_stack ) ;

position_lst = zeros( frames , 3 ) ;

img = zeros( i , j ) ;

x = 1 : 10 ;

tic
for n = 1 : frames
        
    img = double( img_stack( : , : , n ) ) ;
    
    cntrd_x = sum( sum( x .* img , 2) ) / sum( sum( img , 2 ) ) ;
    
    cntrd_y = sum( sum( rot90( flip( x ) ) .* img ) ) / sum( sum( img ) ) ;
    
    position_lst( n , : ) = [ cntrd_x cntrd_y n ] ; 
    
%     pk = pkfnd( img , thresh , 5 ) ;
% 
%     pk_cntrd = cntrd( img, pk , 3 ) ;
%     
%     position_lst( n , : ) = [ pk_cntrd( : , [ 1 2 ] ) n ] ; 
    
end

toc
