function [recsignal] = Noise(chansymb, delta)

sz = size(chansymb);
noise = delta*randn(sz);
recsignal = chansymb + noise;

end