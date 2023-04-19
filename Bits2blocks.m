function [newZtres] = Bits2blocks(newZtres, recbits, Ns, Ne)

len = length(recbits);
bits2d = reshape(recbits, 8, len/8);
bits2d = bits2d';
blk1d = bi2de(bits2d, 2);
blk3d = reshape(blk1d, 8, 8, Ne-Ns+1);
newZtres(:,:,Ns:Ne) = blk3d;

end
