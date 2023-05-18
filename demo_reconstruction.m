clear; close all;

N = 256;            % test image size
F = phantom(N);     % test image

% Sampling mask
L = 7;
MRIMask = MRImask(N, L);
Mask = fftshift(double(MRImask(N, L)));

data = Mask.*fft2(F)/N;

% Algorithm parameters
pm.mu = 1e3;
pm.lambda = 10;
pmTV = pm;

% Reconstruction
uTV = MRITV(Mask, data, pmTV);  fprintf('TV computed.\n')

pm.alpha = 1;
uL12 = MRIL12(Mask, data, pm); fprintf('anisotropic TV - isotropic TV computed.\n')

pm.alpha = 0.6;
uisoETV = MRIETVisotropic(Mask, data, pm); fprintf('isotropic enhanced TV computed.\n')

pm.alpha = 0.8;
uaniETV = MRIETV(Mask, data, pm); fprintf('anisotropic enhanced TV computed.\n')

% Display results
subplot(2,3,1), imshow(F,[]);  colormap('gray');
title('Phantom','interpreter','latex','fontsize',22), 

subplot(2,3,4), imshow(MRIMask); colormap('gray');
title('Sampling Mask','interpreter','latex','fontsize',22)

subplot(2,3,2), imshow(abs(uTV),[]); colormap('gray'); 
title(['TV, error = ',num2str(norm(uTV-F,'fro')/norm(F,'fro'),'%.4g')],'interpreter','latex','fontsize',22);  

subplot(2,3,3), imshow(abs(uL12),[]); colormap('gray'); 
title(['weighted aniTV - isoTV, error = ',num2str(norm(uL12-F,'fro')/norm(F,'fro'),'%.4g')],'interpreter','latex','fontsize',22);  

subplot(2,3,5), imshow(abs(uisoETV),[]); colormap('gray'); 
title(['iso-ETV, error = ',num2str(norm(uisoETV-F,'fro')/norm(F,'fro'),'%.4g')],'interpreter','latex','fontsize',22); 

subplot(2,3,6), imshow(abs(uaniETV),[]); colormap('gray');
title(['aniso-ETV, error = ',num2str(norm(uaniETV-F,'fro')/norm(F,'fro'),'%.4g')],'interpreter','latex','fontsize',22); 
