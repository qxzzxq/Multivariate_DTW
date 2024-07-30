function y = pre_sigmoidNorm(x)

y = 1./(1+exp((-(x-mean(x)))/std(x,1)));