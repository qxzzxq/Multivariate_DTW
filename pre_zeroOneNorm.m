function y = pre_zeroOneNorm(x)

% if std(x,1) == 0
%     y = zeros(length(x),1);
%     return y;
% end
y = (x-min(x))/(max(x)-min(x));