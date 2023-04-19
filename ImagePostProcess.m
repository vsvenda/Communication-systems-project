function []=ImagePostProcess(newZtres,r,c,m,n,minval,maxval)

%% invert the reshaping operation
newZt = reshape(permute(reshape(newZtres,8,8,r,c), [1 3 2 4]), m,n);

%%
%%%%%%%%%% IMAGE POST-PROCESSING %%%%%%%%%%%%%%%%
temp=im2double(newZt)*(maxval-minval)+minval;
fun=@idct2;
newZ=blkproc(temp,[8 8],fun);

if n~=200
    newZ = uint8(newZ);
end

figure;
imshow(newZ);

end