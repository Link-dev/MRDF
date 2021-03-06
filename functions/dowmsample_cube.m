%This file is downloaded from https://github.com/qunmingwang/Code-for-S2-fusion
%Copyright belongs to Dr. Qunming Wang

function S=dowmsample_cube(cube,s,w,PSF);
[a0,b0,c0]=size(cube);
for k=1:c0
    S(:,:,k)=dowmsample_plane(cube(:,:,k),s,w,PSF);
end