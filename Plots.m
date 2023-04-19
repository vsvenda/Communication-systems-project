% SCRIPT FOR RESULT PLOTTING

% Modulation---------------------------------------------------------------
eyeDiagramOffset = 32;
if flagPulseShape==1
    figure, hold on, grid on, plot(symbols(1:30/ts));
    title('Symbols after modulation (half-sine)');

    eyediagram(symbols, 32, T, 15);
    title('Eye diagram after modulation (half-sine)');
elseif flagPulseShape==2
    figure, hold on, grid on, plot(symbols(1:30/ts));
    title('Symbols after modulation (SRRC)');

    eyediagram(symbols(2*K*eyeDiagramOffset+1:end-2), 32, T, 31);
    title('Eye diagram after modulation (SRRC)');
end

if flagPulseShape==1
    t=ts:ts:1;
    figure, plot(t,matchfilt), grid on, title('Half-sine pulse - impulse response');
    figure, freqz(matchfilt), grid on, title('Half-sine pulse - frequency response');
    Power = sum(matchfilt.^2);
elseif flagPulseShape==2
    alpha = 0.5;
    
    figure, hold on, grid on, title('SRRC pulse - impulse response');

    for K = 3:1:6
        t = -K*T+ts:ts:K*T;
        xt = (sin(pi*t/T*(1-alpha))+4*alpha*t/T.*cos(pi*t/T*(1+alpha)))./(pi*t/T.*(1-(4*alpha*t/T).^2));
        xt(find(t==0))=1-alpha+4*alpha/pi;
        xt(find(abs(t)==T/(4*alpha)))=alpha/sqrt(2)*((1+2/pi)*sin(pi/(4*alpha))+(1-2/pi)*cos(pi/(4*alpha)));
        plot(t, xt); 
        Power = sum(xt.^2);
    end
    
    figure, hold on, grid on, title('SRRC pulse - frequency response');
    for K = 3:1:6
        f = -(1+alpha)/(2*T):ts:(1+alpha)/(2*T);
        fedge1 = (1-alpha)/(2*T):ts:(1+alpha)/(2*T);
        fedge2 = -(1+alpha)/(2*T):ts:-(1-alpha)/(2*T);
        xf(find(abs(f)<=(1-alpha)/(2*T)))=sqrt(T);
        xf(find(and(f<=(1+alpha)/(2*T)+eps, f>=(1-alpha)/(2*T)-eps)))=sqrt(T)*cos(pi*T/(2*alpha)*(abs(fedge1)-(1-alpha)/(2*T)));
        xf(find(and(f<=-(1-alpha)/(2*T)+eps, f>=-(1+alpha)/(2*T)-eps)))=sqrt(T)*cos(pi*T/(2*alpha)*(abs(fedge2)-(1-alpha)/(2*T)));
        plot(f, xf); legend('alpha=0.5');
    end   
    
    figure, hold on, grid on, title('SRRC pulse - impulse response');

    K = 6;
    for alpha = 0.1:0.2:0.9
        t = -K*T+ts:ts:K*T;
        xt = (sin(pi*t/T*(1-alpha))+4*alpha*t/T.*cos(pi*t/T*(1+alpha)))./(pi*t/T.*(1-(4*alpha*t/T).^2));
        xt(find(t==0))=1-alpha+4*alpha/pi;
        xt(find(abs(t)==T/(4*alpha)))=alpha/sqrt(2)*((1+2/pi)*sin(pi/(4*alpha))+(1-2/pi)*cos(pi/(4*alpha)));
        plot(t, xt);
    end
    legend('alpha=0.1','alpha=0.3','alpha=0.5','alpha=0.7','alpha=0.9');
    
    figure, hold on, grid on, title('SRRC pulse - frequency response');
    for alpha = 0.1:0.2:0.9
        f = -(1+alpha)/(2*T):ts:(1+alpha)/(2*T);
        fedge1 = (1-alpha)/(2*T):ts:(1+alpha)/(2*T);
        fedge2 = -(1+alpha)/(2*T):ts:-(1-alpha)/(2*T);
        xf = zeros(size(f));
        xf(find(abs(f)<=(1-alpha)/(2*T)))=sqrt(T);
        xf(find(and(f<=(1+alpha)/(2*T)+eps, f>=(1-alpha)/(2*T)-eps)))=sqrt(T)*cos(pi*T/(2*alpha)*(abs(fedge1)-(1-alpha)/(2*T)));
        xf(find(and(f<=-(1-alpha)/(2*T)+eps, f>=-(1+alpha)/(2*T)-eps)))=sqrt(T)*cos(pi*T/(2*alpha)*(abs(fedge2)-(1-alpha)/(2*T)));
        plot(f, xf);
    end
    ylim([0 1.2]);
    legend('alpha=0.1','alpha=0.3','alpha=0.5','alpha=0.7','alpha=0.9');
end

% Channel------------------------------------------------------------------
if flagChannel==1
    hn=-3:0.01:5;
    for k=1:length(hn)
        if hn(k)==0
            ImpulseResponse(k)=1;
        elseif hn(k)==1
            ImpulseResponse(k)=1/2;
        elseif hn(k)==2
            ImpulseResponse(k)=3/4;
        elseif hn(k)==3
            ImpulseResponse(k)=-2/7;
        else
            ImpulseResponse(k)=0;
        end
    end
    b1 = [1 1/2 3/4 -2/7];
    a1 = 1;
elseif flagChannel==2
    hn=-3:0.01:30;
    for k=1:length(hn)
        if hn(k)==0
            ImpulseResponse(k)=0.5;
        elseif hn(k)==1
            ImpulseResponse(k)=1;
        elseif hn(k)==3
            ImpulseResponse(k)=0.63;
        elseif hn(k)==8
            ImpulseResponse(k)=0.25;
        elseif hn(k)==12
            ImpulseResponse(k)=0.16;
        elseif hn(k)==25
            ImpulseResponse(k)=0.1;
        else
            ImpulseResponse(k)=0;
        end
    end
    b1 = [0.5 1 0 0.63 zeros(1,4) 0.25 zeros(1,3) 0.16 zeros(1,13) 0.1];
    a1 = 1;
elseif flagChannel==3
    hn=-3:0.01:10;
    for k=1:length(hn)
        if hn(k)==0
            ImpulseResponse(k)=1;
        elseif hn(k)==1
            ImpulseResponse(k)=0.4365;
        elseif hn(k)==2
            ImpulseResponse(k)=0.1905;
        elseif hn(k)==3
            ImpulseResponse(k)=0.0832;
        elseif hn(k)==5
            ImpulseResponse(k)=0.0158;
        elseif hn(k)==7
            ImpulseResponse(k)=0.003;
        else
            ImpulseResponse(k)=0;
        end
    end
    b1 = [1 0.4365 0.1905 0.0832 0 0.0158 0 0.003];
    a1 = 1;
end
figure, plot(hn,ImpulseResponse), grid on, title('Channel - Impulse Response');

[h,w] = freqz(b1,a1,1000);
figure, subplot(2,1,1), plot(w/pi,20*log10(abs(h))), grid on;
title('Channel - Frequency response');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Magnitude (dB)');
subplot(2,1,2), plot(w/pi,angle(h)), grid on;
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Phase (rad)');

figure, hold on, plot(symbols(1:30/ts));
plot(chansymb(1:30/ts)), grid on, title('Modulated bit stream after the channel');

if flagPulseShape==1
    eyediagram(chansymb,32,T,15), title('Eye diagram - Channel output');
elseif flagPulseShape==2
    eyediagram(chansymb(2*K*eyeDiagramOffset+1:end-2),32,T,31), title('Eye diagram - Channel output');
end

% Noise--------------------------------------------------------------------
if flagPulseShape==1
    eyediagram(recsignal, 32, T, 15), title('Eye diagram with noise');
elseif flagPulseShape==2
    eyediagram(recsignal(2*K*eyeDiagramOffset+1:end-2),32,T,31), title('Eye diagram with noise');
end

% Matched filter-----------------------------------------------------------
if flagPulseShape==1
    t=ts:ts:1;
    figure, plot(t,matchfilt), grid on, title('Matched filter (half-sine) - impulse response');
    figure, freqz(matchfilt), grid on, title('Matched filter (half-sine) - frequency response');
    eyediagram(mstream(2*K*eyeDiagramOffset+1:end-32), 32, 2*T, 30), title('Eye diagram after the Matched Filter (half-sine);');
elseif flagPulseShape==2
    t=ts-K*T+T:ts:K*T+T;
    figure, plot(t,matchfilt), grid on, title('Matched filter (SRRC) - impulse response');
    figure, freqz(matchfilt), grid on, title('Matched filter (SRRC) - frequency response');
    eyediagram(mstream(2*K*eyeDiagramOffset+1:end-7*32),32,2*T,30), title('Eye diagram after the Matched Filter (SRRC);');  
end

% Equalizer----------------------------------------------------------------
if flagEqualizer==1
    if flagChannel==1
        b1=[1];
        a1=[1 1/2 3/4 -2/7];
        a2=[1 zeros(1,31) 1/2  zeros(1,31)  3/4  zeros(1,31) -2/7 zeros(1,31)];
    elseif flagChannel==2
        b1=[1];
        a1=[0.5 1 0 0.63 zeros(1,4) 0.25 zeros(1,3) 0.16 zeros(1,13) 0.1];
        a2=[0.5 zeros(1,31) 1 zeros(1,63) 0.63 zeros(1,127) 0.25 zeros(1,159)...
            0.16 zeros(1,415) 0.1];
    elseif flagChannel==3
        b1=[1];
        a1=[1 0.4365 0.1905 0.0832 0 0.0158 0 0.003];
        a2=[1 zeros(1,31) 0.4365 zeros(1,31) 0.1905 zeros(1,31) 0.0832...
            zeros(1,63) 0.0158 zeros(1,63) 0.003];
    end        
    [h,w]=freqz(b1,a1,1000);
    figure, subplot(2,1,1), plot(w/pi,20*log10(abs(h))), grid on;
    title('Equalizer (ZF) - Frequency response');
    xlabel('Normalized Frequency (\times\pi rad/sample)');
    ylabel('Magnitude (dB)');
    subplot(2,1,2), plot(w/pi,angle(h)), grid on;
    xlabel('Normalized Frequency (\times\pi rad/sample)');
    ylabel('Phase (rad)');
   
    [h,w]=freqz(b1,a2,1000);
    %figure, plot(abs(ifft(h))), title('Equalizer (ZF) - Impulse Response'); 
     if flagPulseShape==1
        eyediagram(equalSig(2*K*eyeDiagramOffset+1:end-32), 32, 2*T, 30), title('Equalizer (Half-sine, ZF) - Eye diagram'); 
    elseif flagPulseShape==2
        eyediagram(equalSig(2*K*eyeDiagramOffset+1:end-7*32),32,2*T,31), title('Equalizer (SRRC, ZF) - Eye diagram'); 
    end
    
elseif flagEqualizer==2
    b1 = [1 1/2 3/4 -2/7];
    a1 = 1;
    [h,w] = freqz(b1,a1,1000);
    Qmmse=conj(h)./(abs(h).^2+2*delta^2);
    figure, subplot(2,1,1), plot(w/pi,20*log10(abs(Qmmse))), grid on;
    title('Equalizer (MMSE) - Frequency response');
    xlabel('Normalized Frequency (\times\pi rad/sample)');
    ylabel('Magnitude (dB)');
    subplot(2,1,2), plot(w/pi,angle(Qmmse)), grid on;
    xlabel('Normalized Frequency (\times\pi rad/sample)');
    ylabel('Phase (rad)');
    b1 = [1 zeros(1,31) 1/2 zeros(1,31) 3/4 zeros(1,31) -2/7 zeros(1,31)];
    a1 = 1;
    [h,w]=freqz(b1,a1,1000);
    %figure, plot(abs(ifft(h))), title('Equalizer (MMSE) - Impulse Response');
    if flagPulseShape==1
        eyediagram(equalSig(2*K*eyeDiagramOffset+1:end-9*32), 32, 2*T, 30), title('Equalizer (Half-sine, MMSE) - Eye diagram'); 
    elseif flagPulseShape==2
        eyediagram(equalSig(2*K*eyeDiagramOffset+1:end-9*32),32,2*T,30), title('Equalizer (SRRC, MMSE) - Eye diagram'); 
    end
end 



figure
hold on
grid on

len = length(equalSig);

if flagPulseShape == 1
    offset = 31;
    nend = 30;
    tx = offset:32:len;
    val = 16*sign(equalSig(tx)).*ones(size(tx));
    h = stem(tx(1:nend), val(1:nend), ':kx'); hold on;
    hbase = h.BaseLine;
    hbase.LineStyle = ':';
    
    
elseif flagPulseShape == 2
    offset = 223;
    nend = 45;
    tx = offset:32:len;
    val = 32*sign(equalSig(tx)).*ones(size(tx));
    h = stem(tx(1:nend), val(1:nend), ':kx'); hold on;
    hbase = h.BaseLine;
    hbase.LineStyle = ':';
end

plot(mstream(1:nend*32));
plot(equalSig(1:nend*32));

legend('Sampled points', 'MF output', 'Equalizer output');
title('Equalizer output');
