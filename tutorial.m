%
% Particle Tracking Tutorial
%
% Origanal: https://site.physics.georgetown.edu/matlab/tutorial.html
%
%
%
%
%
%
%
%
close all ; clear all ;


image_array = imread( 'test.tif' ) ;

est_noise = floor( mean( mean( diff( image_array ), 2 ) ) ) ;

excl_dia = 9 ;
threshold = 80 ;
time_1=zeros(100,1);
time_2=zeros(100,1);
time_3=zeros(100,1);
for n = 1:1;
    tic
filtered_image = bpass( image_array, est_noise + 1, excl_dia ) ;
time_1(n) = toc;
tic
colormap('gray'), imagesc( image_array ) ;
est_pks = pkfnd( filtered_image, threshold, excl_dia ) ;
time_2(n) = toc;
tic
particles = cntrd( filtered_image, est_pks, excl_dia, true);

viscircles( [particles(:,1), particles(:,2)], particles(:,4), 'Color','g')  ;

figure ; plot(particles(:,4),particles(:,3),'.') 

time_3(n) = toc ;

end