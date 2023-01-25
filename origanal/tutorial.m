


close all ; clear all ;


image_array = double( imread( '../test.tif' ) );

est_noise = floor( mean( mean( diff( image_array ), 2 ) ) ) ;

excl_dia = 15 ;
threshold = 41 ;

loop = 1 ;

time_bpass = zeros( loop, 1 ) ;
time_pkfnd = zeros( loop, 1 ) ;
time_cntrd = zeros( loop, 1 ) ;

for n = 1:loop
    tic
    filtered_image = bpass( image_array, est_noise, excl_dia ) ;
    time_bpass( n ) = toc;
    
    % colormap('gray'), imagesc( image_array ) ;
    tic
    est_pks = pkfnd( image_array, threshold, excl_dia ) ;
    time_pkfnd( n ) = toc ;
    
    tic
    particles = cntrd( image_array, est_pks, excl_dia );
    time_cntrd( n ) = toc ;
    
    % viscircles( [particles(:,1), particles(:,2)], particles(:,4) / 2 * .6, 'Color','g')  ;

    % figure ; plot(particles(:,4),particles(:,3),'.') 

end

time_bpass = mean( time_bpass )
time_pkfnd = mean( time_pkfnd )
time_cntrd = mean( time_cntrd )
total = time_bpass + time_pkfnd + time_cntrd
size( particles )