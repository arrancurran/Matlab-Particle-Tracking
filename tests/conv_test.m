% Testing conv methods and their speeds
img = rand(2000,2000);
n = 2000 ;

normalize       = @( x ) x / sum( x ) ;

gaussian_x = - ceil( 2 ) : ceil( 2 ) ;

gaussian_kernel = normalize(...
            exp( -( gaussian_x / ( 2 * 2 ) ) .^2 ) ...
        ) ;

for k = 1:n
    tic
    img_gauss = conv2( img , gaussian_kernel, 'same' ) ;
    img_gauss = conv2( img_gauss , gaussian_kernel', 'same' ) ;
    two_step_man(k) = toc ;

end

for k = 1:n
    tic
    img_gauss = conv2( gaussian_kernel', gaussian_kernel, img, 'same' ) ;
    one_step_man(k) = toc ;
end

plot(two_step_man, 'k') ; hold on
plot(one_step_man, 'b') ; hold off

mean(two_step_man)
mean(one_step_man)
