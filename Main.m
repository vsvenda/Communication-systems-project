% MAIN SCRIPT

close all
clear all

N = 3;
flagPicture = 2;        % 1 - Trees; 2 - Elephant 
flagPulseShape = 1;     % 1 - half-sine; 2 - SRRC
delta = 0.4;            % amount of noise
flagEqualizer = 2;      % 1 - ZF; 2 - MMSE
alpha = 0.5;            % roll-off factor 
ts = 1/32;              % sample duration
T = 1;                  % period
K = 6 ;                 % truncation length
A = 1/1.414;            % normalization factor
flagChannel = 1;        % 1 - basic; 2 - outdoor; 3 - indoor

disp('****** Image pre-processing ******');
[Ztres,r,c,m,n,minval,maxval] = ImagePreProcess(flagPicture);

newZtres = uint8(zeros(size(Ztres)));

[b a] = FilterParameters(flagChannel);

for ntr = 1:N:r*c
    [bits] = Blocks2bits(Ztres, ntr, min(ntr+N-1, r*c));
    [symbols] = Modulation(bits, flagPulseShape);
    [chansymb] = Channel(symbols, a, b);
    [recsignal] = Noise(chansymb, delta);
    [mstream, matchfilt, SNR] = MatchedFilter(flagPulseShape, recsignal, T, K, A, alpha, ts, delta);
    [equalSig] = Equalizer(mstream, b, a, delta, flagEqualizer, flagPulseShape);
    [recbits] = SamplingDetection(equalSig, flagPulseShape);
    [newZtres] = Bits2blocks(newZtres, recbits, ntr, min(ntr+N-1, r*c));
    
    if abs(ntr - r*c) < N; Plots; 
    else; disp(ntr); end;
end

disp('****** Image post-procesing ******');
ImagePostProcess(newZtres,r,c,m,n,minval,maxval)
