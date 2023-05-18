clear; close all


N = 128;

I = zeros(N,N);
I(:,1:32) = 0.02;
I(:,33:64) = 1/3;
I(:,65:96) = .98;
I(:,97:end) = 2/3;
[x,y] = size(I);


sig = 0.6;
f = I + sig*randn(size(I)); 

pm.mu = 0.8;
pm.lambda = 1;

% Denoise
[uTV,error7] = denoiseTV(f,pm);

pm.maxDCA = 10;
[uL1L22_3,error6] =denoiseETV(f, 1.2, pm); % in the paper, we stop at two outer iterations

axes('position',[0.01,0.52,0.23,.4]),
h0 = imshow(I);
title(['Groundtruth, SSIM=1.0000'],'interpreter','latex','fontsize',22)
h0_axis = gca; 

axes('position',[0.25,0.52,0.23,.4]),
h1 = imshow(f);
title(['Noisy image, SSIM=',num2str(ssim(f,I),'%.4f')],'interpreter','latex','fontsize',22)
h1_axis = gca; 

axes('position',[0.49,0.52,0.23,.4]),
h2 = imshow(uTV); h2_axis = gca;  
title(['TV, SSIM=', num2str(ssim(uTV,I),'%.4f')],'interpreter','latex','fontsize',22)
 
axes('position',[0.73,0.52,0.23,.4]),
h8 = imshow(uL1L22_3); h8_axis = gca; 
 title(['Enhanced TV, SSIM=',num2str(ssim(uL1L22_3,I),'%.4f')],'interpreter','latex','fontsize',22)



axes('position',[0.04,0.05,0.2,.4]),
[a,b,c] = improfile(I,[1,128],[64,64]);
plot(a,c,'k') 
set(gca,'xlabel',[]),xlim([1,128]),set(gca,'FontSize',15)

axes('position',[0.28,0.05,0.2,.4]),
[a,b,c] = improfile(f,[1,128],[64,64]);
plot(a,c,'k') 
set(gca,'xlabel',[]),xlim([1,128]),set(gca,'FontSize',15)

axes('position',[0.52,0.05,0.2,.4]), 
[a,b,c] = improfile(uTV,[1,128],[64,64]);
plot(a,c,'k') 
set(gca,'xlabel',[]),xlim([1,128]),set(gca,'FontSize',15)

axes('position',[0.76,0.05,0.2,.4]),
[a,b,c] = improfile(uL1L22_3,[1,128],[64,64]);
plot(a,c,'k') 
set(gca,'xlabel',[]),xlim([1,128]),set(gca,'FontSize',15)


