function y = pre_tanhNorm(x)

y = 0.5*(tanh(0.01*(x-mean(x))/std(x,1))+1);