%
% Particle Tracking Tutorial
%
% Origanal: https://site.physics.georgetown.edu/matlab/tutorial.html
%
%

close all ; clear all ;

addpath( '../src/', '../img/' )

image_array = double( imread( '../img/tutorial_2048.tif' ) ) ;

% image_array = image_array( 1 : 512, 1 : 512 ) ;

% % img = img(1:512,1:512);
% [p,~] = size( image_array ) ;
% for n = 1 : p

%     limit = ( (p / 2 ) - n ) ;
%     pin( :, n ) = linspace( -limit, limit, p ) ;

% end
% pin = abs( pin ) ;
% pin = pin / max(pin, [], 'all' ) ;

% image_array = image_array - (120 * pin) ;

% image_array = image_array - min( image_array, [], 'all' ) ;

% image_array = image_array / max( image_array, [], 'all' ) * 255 ;

excl_dia = 15 ;
excl_rad = floor( excl_dia / 2 ) ;
baseline = 40 ;

loop = 1 ;

time_bpass = zeros( loop, 1 ) ;
time_pkfnd = zeros( loop, 1 ) ;
time_cntrd = zeros( loop, 1 ) ;

for n = 1 : loop

    tic
    filtered_image = bpass( image_array, true, false, baseline ) ;
    time_bpass( n ) = toc;
    
    image_array_fig = figure ; colormap('gray'), imagesc( image_array ) ; axis square ;
    
    tic
    est_pks = pkfnd( filtered_image, 40, excl_dia ) ;
    time_pkfnd = toc ;
    
    tic
    particles = cntrd( image_array, est_pks, excl_dia, true );
    time_cntrd( n ) = toc ;
    
end

time_bpass = mean( time_bpass )
time_pkfnd = mean( time_pkfnd )
time_cntrd = mean( time_cntrd )
total = time_bpass + time_pkfnd + time_cntrd
size( particles )


% for p = 1 : length( particles( :, 1 ) )
%     x = particles( p, 1 ) ;
%     y = particles( p, 2 ) ;
%     line( [ x - excl_rad, x + excl_rad ], [ y, y ], 'Color','green' )
%     line( [ x, x ], [ y - excl_rad, y + excl_rad ], 'Color','green' )
%     line( [ x - excl_rad, x - excl_rad ], [ y - excl_rad, y + excl_rad ], 'Color','green' )
%     line( [ x + excl_rad, x + excl_rad ], [ y - excl_rad, y + excl_rad ], 'Color','green' )
%     line( [ x - excl_rad, x + excl_rad ], [ y - excl_rad, y - excl_rad ], 'Color','green' )
%     line( [ x + excl_rad, x - excl_rad ], [ y + excl_rad, y + excl_rad ], 'Color','green' )
% end

% crc = viscircles( [ particles( :, 1 ), particles( :, 2 ) ], particles( :, 4 ) / 2, 'Color', 'g', 'EnhanceVisibility', false, 'LineWidth', 1 ) ;

% figure ; plot( particles(:,4), particles(:,3), 'o') 