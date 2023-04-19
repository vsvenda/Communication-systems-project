function [b a] = FilterParameters(flagChannel)

if flagChannel == 1
        b = [1];
        a = [1 zeros(1,31) 1/2 zeros(1,31) 3/4 zeros(1,31) -2/7 zeros(1,31)];
elseif flagChannel == 2 
        b = [1];
        a = [0.5 zeros(1,31) 1 zeros(1,63) 0.63 zeros(1,159) 0.25...
           zeros(1,127) 0.16 zeros(1,415) 0.1 zeros(1,31)];
elseif flagChannel == 3
        b = [1];
        a = [1 zeros(1,31) 0.4365 zeros(1,31) 0.1905 zeros(1,31) 0.0832...
            zeros(1,63) 0.0158 zeros(1,63) 0.003 zeros(1,31)];       
end  

end