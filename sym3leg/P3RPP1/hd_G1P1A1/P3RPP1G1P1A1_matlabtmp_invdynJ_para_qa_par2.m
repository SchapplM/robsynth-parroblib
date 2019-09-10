
%% Function calls and calculation
tauX = P3RPP1G1P1A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges);
Jinv = P3RPP1G1P1A1_Jinv(xP, qJ, pkin, koppelP, legFrame);
tauA  = (Jinv') \ tauX;
