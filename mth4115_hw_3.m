clear 

matrix = [[1,2,2];[4,3,3];[3,2,4]];
naive_ge(matrix)

function arr = naive_ge(m)
%   for each column up to the n-1 ** must retain final nth pivot
    for j = 1 : length(m)-1
%       if zero pivot encountered end the script 
        if abs(m(j,j))<1e-6; error("zero pivot encountered"); end
%       for each row below the jth row being active to end of matrix 
        for i = j+1 : length(m)
%       multiplication/scalar factor for performing G.E. 
            mult = m(i,j)/m(j,j);
            wipe = m(j,:).*mult;
            m(i,:) = m(i,:) - wipe;
        end
    end
    
    arr = m;
    
%     for i = length(m) : -1 : 1
%         for j = i+1 : lenght(m)
%             b(i) = b(i) - a(i,j)*x(j);
%         end
%         x(i) = b(i)/a(i,i);
%     end
end