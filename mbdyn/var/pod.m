function [V,S,M,Aout,data,A,BB] = pod(file, stp, posgdl, ns, novel)

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
 %        A =[A, data(1:stp, [2+(i-1)*12 , 6+(i-1)*12])];
    end
    A = [A, data(1:stp, 12*posgdl+1:end)];
    jp = 6;
    %jp = 2;
else
    A = data(1:stp, :);
    jp = 12;
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
    Bmovx(i,1:posgdl) = B(i,1:jp:posgdl*jp)/S(i);
    Bmovy(i,1:posgdl) = B(i,2:jp:posgdl*jp)/S(i);
    Bmovz(i,1:posgdl) = B(i,3:jp:posgdl*jp)/S(i);
    Bmovp(i,1:posgdl) = B(i,4:jp:posgdl*jp)/S(i);
    Bmovq(i,1:posgdl) = B(i,5:jp:posgdl*jp)/S(i);
    Bmovr(i,1:posgdl) = B(i,6:jp:posgdl*jp)/S(i);
    BB(i,:) = B(i,:)/S(i)^2;
end
Aout =BB*A';
mx = mn(1:jp:posgdl*jp).*scl(1:jp:posgdl*jp);
my = mn(2:jp:posgdl*jp).*scl(2:jp:posgdl*jp);
mz = mn(3:jp:posgdl*jp).*scl(3:jp:posgdl*jp);
mp = mn(4:jp:posgdl*jp).*scl(4:jp:posgdl*jp);
mq = mn(5:jp:posgdl*jp).*scl(5:jp:posgdl*jp);
mr = mn(6:jp:posgdl*jp).*scl(6:jp:posgdl*jp);

for i=1:ns
 %   s = max(abs([Bmovx(i,:), Bmovy(i,:), Bmovz(i,:)]))
 %if (s > 1.e-15) 
        V((i-1)*6+1,:) = Bmovx(i,:).*scl(1:jp:posgdl*jp);
        V((i-1)*6+2,:) = Bmovy(i,:).*scl(2:jp:posgdl*jp);
        V((i-1)*6+3,:) = Bmovz(i,:).*scl(3:jp:posgdl*jp);
        V((i-1)*6+4,:) = Bmovp(i,:).*scl(4:jp:posgdl*jp);
        V((i-1)*6+5,:) = Bmovq(i,:).*scl(5:jp:posgdl*jp);
        V((i-1)*6+6,:) = Bmovr(i,:).*scl(6:jp:posgdl*jp);
        %  end
end
M = [mx;my;mz];
%M = mx;