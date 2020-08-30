function f = ricker(fs,ts,params)
%% ricker wavlet

t = 0:params.dt:params.T;
f = (1 - 2*(pi*fs*(t(:) - ts)).^2).*exp(-(pi*fs*(t(:) - ts)).^2);

end
