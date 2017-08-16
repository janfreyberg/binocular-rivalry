function [f, f_inv]=make_circular_checkerboard_pattern(sigma,spokes, resolution)
% usage 
% [f f_inv]=make_circular_checkerboard_pattern(2,4);
% $author Shripad Kondra, NBRC, India (28-09-2010)
if (nargin==0)
sigma=4;
spokes=6;
end
SUP=resolution; % This parameter controls the resolution
hsup=(SUP-1)/2;
[x,y]=meshgrid([-hsup:hsup]);
[THETA,r] = cart2pol(x,y);
r=(r./(SUP/2))*pi;
% r(r<0.04)=0; % uncomment to put a dot at the centre.
% r(r>(pi+0.01))=inf; % uncomment if you want to get exact circle
%%
f=sin(r*sigma);         % 1st concentric filter
f1=sin(THETA*spokes);   % 1st radial filter
f1=f1>=0;               % binarize
f11=f.*f1;              % point-wise multiply
f=sin(r*sigma+pi);      % 2nd concentric filter shifted by pi
f1=sin(THETA*spokes+pi);% 2nd radial filter shifted by pi
f1=f1>=0;               % binarize
f12 = f.*f1;            % point-wise multiply
f=(f11+f12)>=0;         % add the two filters and threshold
f_inv=(f11+f12)<=0;     % add the two filters and threshold
return

