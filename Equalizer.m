function [equalSig] = Equalizer ( mstream, b, a, delta, flagEqualizer, flagPulseShape)

if flagEqualizer == 1   
    h = a;
    lenh = length(h);
    lenm = length(mstream);
    hfft = fft(h, 2*(lenm+lenh));
    
    Qmmse = conj(hfft)./((abs(hfft)).^2);
    if flagPulseShape==1
        Qmmse = Qmmse(:);
    end
    
    streamfft = fft(mstream, 2*(lenm+lenh));
    equalSig = ifft(streamfft.*Qmmse);
    
    equalSig = equalSig(1:length(mstream));
    
elseif flagEqualizer == 2   
    h = a;
    lenh = length(h);
    lenm = length(mstream);
    hfft = fft(h, 2*(lenm+lenh));
    
    Qmmse = conj(hfft)./((abs(hfft)).^2+2*delta^2);
    if flagPulseShape==1
        Qmmse = Qmmse(:);
    end
    
    streamfft = fft(mstream, 2*(lenm+lenh));
    equalSig = ifft(streamfft.*Qmmse);
    
    equalSig = equalSig(1:length(mstream));
end

end