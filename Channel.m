function [ bitStreamChannel ] = Channel ( symbols, b, a )

bitStreamChannel = filter(b,a,symbols);

end

