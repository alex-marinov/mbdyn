
x = [-10:2:10]';
xq = [-10:0.1:10]';
y = x.^5    + 2 * (x - 200).^(-1) ;
y_p = 5 * x.^4  - 2 * (x - 200).^-2  ; 
yq = xq.^5  + 2 * (xq - 200).^-1  ;
yq_p = 5 * xq .^4 - 2 * (xq - 200).^-2 ;


tree = kd_buildtree(x,0);

[I,I_i] = f_regression_derivatives_kdtree(xq,x,tree,3,5,-2);
y_reg = I*y;
y_der = I_i{1}*y;
%[a,b,p] = knearn(-1,p,[],1);

figure
plot(xq,yq,'-','linewidth',2) ;
hold on;
plot(xq,y_reg,'-*r','linewidth',1.5) ;
xlabel('x'); ylabel('y')
legend('exact','MLS approximation') ;
grid on;

figure
plot(xq,yq_p,'-','linewidth',2) ;
hold on;
plot(xq,y_der,'-*r','linewidth',1.5) ;
xlabel('x'); ylabel('dy/dx')
legend('exact','MLS approximation') ;
grid on;

%% ERRORS:
figure
plot(xq,(yq-y_reg),'-*r','linewidth',1.5) ;
xlabel('x'); ylabel('relative error')
legend('function error') ;
grid on;

figure
plot(xq,(yq_p-y_der),'-*r','linewidth',1.5) ;
xlabel('x'); ylabel('dy/dx relative error')
legend('derivative error') ;
grid on;
