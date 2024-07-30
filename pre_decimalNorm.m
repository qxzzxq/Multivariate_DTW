function y = pre_decimalNorm(x)

d = floor(log10(max(x)))+1;
y = x/(10^d);
