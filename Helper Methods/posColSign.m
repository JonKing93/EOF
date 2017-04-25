function[V] = posColSign(V)
%% Reverses the sign of columns whose elements are majority negative.

if ~ismatrix(V)
    error('posColSign is for 2D matrices');
end

[m,n] = size(V);

for k = 1:n
    if sum( V(:,k)<0 ) > (m/2)    
        V(:,k) = -1 * V(:,k);
    end
end
