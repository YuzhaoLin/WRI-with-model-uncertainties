function out = Fmult(in,c,mode,params)
%% forward operator
%
%
%
%

%% basic parameters
t = 0:params.dt:params.T;
r = params.r;

nt = length(t);
nr = length(r);

f = -0.5/params.dt:.5/params.T:0.5/params.dt;
nf = length(f);

%% Green function 
[ff,rr] = ndgrid(ifftshift(f),r);
gh = exp(-1i*(2*pi*ff/c).*rr)./(1 + rr);

%% forward modling for wavelet or backward propagation for wavefield
if mode==1 % forward modling for wavelet 
    qh = fft(padarray(in(:), nf - nt, 'post'));
    uh = bsxfun(@times, gh, qh);
    out = ifft(uh,[],1);
    out = out(1:nt,:);
else % backward propagation for wavefield
    uh = fft(padarray(reshape(in,nt,nr), [nf - nt,0], 'post' ),[],1);
    qh = sum(conj(gh).*uh,2);
    out = ifft(qh,[],1);
    out = out(1:nt);
end
out = out(:);

end
