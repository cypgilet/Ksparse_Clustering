%Convert indicator matrix into class vector 
function y = indicator2class(Y)
    k = size(Y,2);
    n = size(Y,1); 
    y = zeros(n,1);
    for i = 1:k 
        y((Y(:,i) == 1)) = i ;
    end

end

