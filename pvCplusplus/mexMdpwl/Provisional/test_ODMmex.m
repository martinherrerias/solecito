
para(1).Iph = 11.842363;
para(1).Io = 2.394890881115644e-09;
para(1).Rs =0.315 ;
para(1).nVth = 2.000836411141469;
para(1).Rsh = 4.028255106927194e+02;
para(1).di2mutau = 0.0;
para(1).Vbi = 64.8;
para(1).Brev = 2.102623456790123e-06;

%{
para(2).Iph = 12.842363;
para(2).Io = 2.394890881115644e-09;
para(2).Rs =0.315 ;
para(2).nVth = 2.000836411141469;
para(2).Rsh = 4.028255106927194e+02;
para(2).di2mutau = 0.0;
para(2).Vbi = 64.8;
para(2).Brev = 2.102623456790123e-06;
%}

tic()
%[out1, out2, size] = mexVecODMpwlapprox(para, 0, [-20.2, 18.2]);
[out1, out2, size] = mexODMpwlapprox(para(1));
%[out1, out2, size] = mexODMpwlapprox(para(1), 0);
%mexVecODMpwlapprox(para, 0, [-20.2, 18.2]);
toc()
