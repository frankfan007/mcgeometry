function [ img ] = i2m_open( img, arg )
% Close using iso2mesh functions, fill holes as middle step
%
% this function is written by:
% % Perdue KL, Diamond SG; 
% T1 magnetic resonance imaging head segmentation for diffuse optical tomography and electroencephalography. 
% J. Biomed. Opt. 19(2):026011.  doi:10.1117/1.JBO.19.2.026011.

%%
 img=thinbinvol(img, arg);
 img=fillholes3d(img, 0);
 img=thickenbinvol(img, arg);

end

