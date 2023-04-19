function [symbols] = Modulation(bits, flag)

T = 1;
ts = T/32;
t = ts:ts:T;
A = 1/1.414;

lent = length(t);

ind0 = find(bits==0);
ind1 = find(bits==1);

if flag == 1
    % half-sine pulse shaping function
    for i = 1:1:length(ind0)
        ind = ind0(i)-1;
        symbols(ind*lent+1:(ind+1)*lent,1) = -sin(pi*t/T);
    end

    for i = 1:1:length(ind1)
        ind = ind1(i)-1;
        symbols(ind*lent+1:(ind+1)*lent,1) = sin(pi*t/T);
    end

elseif flag == 2
    alpha = 0.5;

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
    
    xt = A*xt;

    t = T-K*T+ts:ts:T+K*T;
    xrt = xt;

    xtRC = conv(xt, xrt)*2*K*T/length(t);

    abits(ind1) = 1;
    abits(ind0) = -1;
    
    % half-sine pulse shaping function
    for i = 1:1:length(abits)
        symb = zeros(1, 32);
        for j = -K+1:1:K
            if and(i-j > 0, i-j < length(abits)+1)
                kkk = abits(i-j);
                ppp = xt(((j+K-1)*T)/ts+1:((j+K)*T)/ts);
                symb = symb + kkk*ppp;
            end
        end
        symbols((i-1)*32+1:i*32) = symb;
    end
end

end