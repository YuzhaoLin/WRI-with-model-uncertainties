


clc,clear
%%
r = 40;
s = 10;
a = 1e-2; 
b = 1e-3;
dx = 20;
niter = 5;

%% 
for iter = 1:niter
    n = [iter*dx,iter*dx];
    % operators
    A  = rand(n(1)*n(2),n(1)*n(2));
    P  = rand(n(1)*n(2),r); 
    Q  = rand(n(1)*n(2),s);
    D  = rand(r,s);
    G  = @(u)getG(5,u,n);
    % FWI
    t1 = clock;
    U  =  A  \ Q;
    H  =  P*(D - P'*U);
    V  =  A' \ H;
    g  = 0;
    for k = 1:size(U,2)
        g  =  g + G(U(:,k))' * V(:,k);
    end
    t2 =  clock;
    tfwi(iter) = etime(t2,t1);

    % WRI
    t1 = clock;
    U  = (a * A' * A + b * P * P') \ ( a * A' * Q  + b * P * D );
    V  = (A * U - Q); 
    g  = 0;
    for k = 1:size(U,2)
        g  =  g + G(U(:,k))' * V(:,k);
    end
    t2 =  clock;
    twri(iter) = etime(t2,t1);

    % New formula \sigma = I
    t1 = clock;
    U  =  A  \ Q;
    Hw =  P' * ( A \ (A' \ P) ) ;
    H  =  P  * Hw * (D - P'*U);
    V  =  A' \ H;
    W  =  A  \ V;
    g  = 0;
    for k = 1:size(U,2)
        g  =  g -2*G(U(:,k))' * V(:,k) + 2*G(V(:,k))' * W(:,k);
    end
    t2 =  clock;
    tnewi(iter) = etime(t2,t1);

    % New formula \sigma = qq^*
    t1 = clock;
    U  = A  \ Q;
    D0 = P'*U;
    PM = eye(r,r);
    Hw  = zeros(r,r);
    for is = 1:size(U,2)
        Hw = Hw + PM - (PM * D0(:,is) * D0(:,is)' * PM)/(1 + D0(:,is)' * PM * D0(:,is)) ;
    end
    H  = P  * Hw * (D - P'*U);
    V  =  A' \ H;
    W  =  U * Q' * V;
    g  = 0;
    for k = 1:size(U,2)
        g  =  g -2*G(U(:,k))' * V(:,k) + 2*G(V(:,k))' * W(:,k);
    end
    t2 =  clock;
    tnewq(iter) = etime(t2,t1);
end

%% plot
figure;plot([1:niter]*dx,tfwi,'LineWidth',2); hold on
       plot([1:niter]*dx,twri,'r','LineWidth',2); hold on
       plot([1:niter]*dx,tnewi,'g','LineWidth',2); hold on
       plot([1:niter]*dx,tnewq,'y','LineWidth',2); 
       legend('FWI','WRI','New method with Identity matrix','New method with qq^* matrix','Location','NorthWest');
       xlabel('Grid','fontsize',18); 
       ylabel('Time [s]','fontsize',18); set(gca,'fontsize',15); 
       
print(1,'-depsc','-r300','time');
