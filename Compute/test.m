%%
n = [100,100];
r = 40;
s = 10;
a = 1e-2; 
b = 1e-3;

%% operators
A  = rand(n(1)*n(2),n(1)*n(2));
P  = rand(n(1)*n(2),r); 
Q  = rand(n(1)*n(2),s);
D  = rand(r,s);

%% FWI
tic;
U  =  A  \ Q;
V  =  A' \ (P*(D - P'*U));
toc;

%% WRI
tic;
U  = (a * A' * A + b * P * P') \ ( a * A' * Q  + b * P * D );
V  = (A * U - Q); 
toc;

%% New formula \sigma = I
tic;
U  =  A  \ Q;
H  =  P' * inv( A' * A ) * P;
V  =  A' \ (P*(D - P'*U));
W  =  A  \ V;
toc;

%% New formula \sigma = I
tic;
U  =  A  \ Q;
V  =  A' \ (P*(D - P'*U));
toc;

