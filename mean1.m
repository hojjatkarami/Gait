function y = mean1(x,k,dim)
    % average across row dim = 1
if dim==2
    x=x';
end
m = size(x,1); % number of rows
n = size(x,2);   % number of columns
% y=zeros(n*k, m);

    for i = 0:k:n-k
        i;
        y(1:m,(i+k)/k) = mean(x(1:m,i+1:i+k),2);
        y;

    end
  if dim==2
      y=y';
  end

% y(1:m,(n-k+1)/k) = mean(x(1:m,n-k+1:n),dim);

end