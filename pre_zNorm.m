function y = pre_zNorm(x)

% if std(x,1) == 0
%     y = zeros(length(x),1);
%     return y;
% end
y = (x-mean(x))/std(x,1);