function [recbits] = SamplingDetection(equalSig, flag)

% sampling
len = length(equalSig);

if flag == 1
    N = (len-31)/32;
    offset = 31;
elseif flag == 2
    N = (len-383)/32;
    offset = 223;
end


for i = 1:1:N
    flg = equalSig(offset+32*(i-1));
    if  flg > 0
        recbits(i) = uint8(1);
    elseif flg < 0
        recbits(i) = uint8(-1);
    end
end

recbits = recbits';

end