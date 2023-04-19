function [ bitStreamFilter, MatchFilt, SNR] = MatchedFilter ( check, bitStreamChannel, T, K, A, alpha, ts, delta )

A = 1/1.414;
%HALF-SINE ****************************************************************
if check==1
    t=ts:ts:1; %We start at 1/fs so that we have 32 time-stamps
    for k=1:length(t)
        MatchFilt(k)=sin(pi/T*(T-t(k))); %Creating the Matched filter "signal"
    end
    bitStreamFilter=conv(bitStreamChannel,MatchFilt);
end
%RAISED-COSINE*************************************************************
if check==2

    for K = 2:1:6
        t = -K*T+ts:ts:K*T;
        xt = (sin(pi*t/T*(1-alpha))+4*alpha*t/T.*cos(pi*t/T*(1+alpha)))./(pi*t/T.*(1-(4*alpha*t/T).^2));
        xt(find(t==0))=1-alpha+4*alpha/pi;
        xt(find(abs(t)==T/(4*alpha)))=alpha/sqrt(2)*((1+2/pi)*sin(pi/(4*alpha))+(1-2/pi)*cos(pi/(4*alpha)));
        f = -(1+alpha)/(2*T):ts:(1+alpha)/(2*T);
        fedge1 = (1-alpha)/(2*T):ts:(1+alpha)/(2*T);
        fedge2 = -(1+alpha)/(2*T):ts:-(1-alpha)/(2*T);
        xf(find(abs(f)<=(1-alpha)/(2*T)))=sqrt(T);
        xf(find(and(f<=(1+alpha)/(2*T), f>=(1-alpha)/(2*T))))=sqrt(T)*cos(pi*T/(2*alpha)*(abs(fedge1)-(1-alpha)/(2*T)));
        xf(find(and(f<=-(1-alpha)/(2*T), f>=-(1+alpha)/(2*T))))=sqrt(T)*cos(pi*T/(2*alpha)*(abs(fedge2)-(1-alpha)/(2*T)));
    end
    %Matched filter output-------------------------------------------------
    MatchFilt = A*xt;
    bitStreamFilter=conv(bitStreamChannel,MatchFilt);
end

if delta~=0
    Power = bitStreamFilter.^2;
    totalP = sum(Power);
    SNR = 10*log(totalP/length(Power)/delta^2);
else
    SNR=0;
end


end

