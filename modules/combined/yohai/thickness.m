E1 = 4e4;
E0 = 2e5;
nu1 = 0.3;
nu0 = 0.3;

psic = 0.2;
sigmac = sqrt(2*E1*psic);
Gc = 0.04;

E1_ = E1/(1-nu1^2);
E0_ = E0/(1-nu0^2);
alpha = (E1_-E0_)/(E1_+E0_);
g = (1.258-0.4*alpha-0.26*alpha^3-0.3*alpha^4)/(1-alpha);
c = 2/pi/g;

Hc = c*E1_*Gc/sigmac^2