function y = pre_absOneNorm(x)

% if std(x,1) == 0
%     y = zeros(length(x),1);
%     return y;
% end
x_max = max(abs(x));
y = x/x_max;