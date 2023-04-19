function [Ztres,r,c,m,n,minval,maxval]=ImagePreProcess(flagPicture)

if flagPicture==1
    load trees
    Z=ind2gray(X(1:200,1:200),gray);
    figure;
    imshow(Z);
elseif flagPicture==2
    [X,map] = imread('Elephant.jpg');
    Z=rgb2gray(X);

    figure;
    imshow(X, map);
    figure;
    imshow(Z);
end

%% take DCT in 8x8 blocks and quantize it into 8-bit numbers
fun = @dct2;
temp = blkproc(Z,[8 8],fun);
% scale DCT to [0,1]
minval = min(min(temp));
maxval = max(max(temp));
xformed = (temp-minval)/(maxval-minval);
% quantize DCT coefficients to 256 levels
Zt=im2uint8(xformed);

%%%%%%%%%% CONVERSION TO A BIT STREAM %%%%%%%%%%%%%%%%

%% reshape DCT into 8x8 blocks for ease of transmission
[m, n] = size(Zt);
r=floor(m/8); c=floor(n/8);
Ztres = reshape(permute(reshape(Zt,8,r,8,c), [1 3 2 4]), 8,8,r*c);

end