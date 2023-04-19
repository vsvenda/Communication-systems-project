function [bits] = Blocks2bits(blocks, Ns, Ne)

blk3d = blocks(:,:,Ns:Ne);
blk1d = reshape(blk3d, 8*8*(Ne-Ns+1), 1);
bits2d = de2bi(blk1d, 8);
bits2dt = bits2d';
bits = bits2dt(:);

end
