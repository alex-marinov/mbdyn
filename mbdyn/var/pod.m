function [S,Aout,data,A,BB,mn] = pod(file, stp, ns, novel)

data = dlmread(file, ' ');
[N, m] = size(data)
if ( stp > N)
    error('too many steps required');
end
if ( ns > m)
    error('too many sv required');
end
if (novel == 1) 
    A = [];
    for i = 1:posgdl
        A =[A, data(1:stp,[1+(i-1)*12, 2+(i-1)*12 , 3+(i-1)*12,   4+(i-1)*12, 5+(i-1)*12 , 6+(i-1)*12])];
    end
    A = [A, data(1:stp, 12*posgdl+1:end)];
else
    A = data(1:stp, :);
end

[r,c] = size(A);
for i=1:c
    scl(i) = max(abs(A(:,i)));
    if (scl(i) ~= 0)
        A(:,i) = A(:,i)/scl(i);
    end
    mn(i) = mean(A(:,i));
    A(:,i) = A(:,i)-mn(i);
end

[U,E] = eigs(A*A',ns);
B =U'*A;
for i = 1:ns
    S(i) = norm(B(i,:));
    BB(i,:) = B(i,:)/S(i)^2;
end
Aout =BB*A';

% eigenvectors:
% plot(([1 1;sqrt(-1) -sqrt(-1)]*vv(:,[6 7])'*BB(1:12,3:12:108))')
