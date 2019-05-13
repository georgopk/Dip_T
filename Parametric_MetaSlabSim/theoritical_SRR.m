%An Accurate Approach of Mathematical Modeling of SRR and SR for Metamaterials
u=1e-3;
eps0=8.854*10^-12;
mu0=4*pi*10^-7;
h=1.52;
epsr=4.2;
len=10*u:5*u:70*u; %strip length
s=0.25*u:0.05*u:2*u; %gap between rings
w=1*u; %strip width
epsr_sub = 1+2*atan(h./(2*pi*(w+s)))*(epsr-1)/pi;
k = s/(s+w);
for i=1:length(len)
    Csrr = (2*len(i)-3*(w+s))*eps0.*epsr_sub*ellipticK(sqrt(1-k^2))/(2*ellipticK(k));
    rho = (s+w)/(len(i)-s-w);
    Lsrr = mu0*(len(i)-(w+s))*4.86*(log(0.98/rho)+1.84*rho);
    fr = 1./(2*pi*sqrt(Lsrr.*Csrr));
    plot(s, fr);
    hold on
end