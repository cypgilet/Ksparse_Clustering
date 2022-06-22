
% Convert class vector into indicator matrix
function Y = class2indicator(y,k)
  n = length(y); 
  Y = zeros(n,k);
  for i = 1:k 
     Y(:,i) = (y == i);  
  end
  Y = double(Y);  
end

