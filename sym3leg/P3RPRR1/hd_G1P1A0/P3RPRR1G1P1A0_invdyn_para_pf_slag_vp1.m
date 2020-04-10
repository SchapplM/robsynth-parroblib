% Calculate vector of inverse dynamics forces for parallel robot
% P3RPRR1G1P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:23
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRR1G1P1A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:23:10
% EndTime: 2020-03-09 21:23:11
% DurationCPUTime: 1.23s
% Computational Cost: add. (6079->245), mult. (5092->400), div. (930->7), fcn. (3132->56), ass. (0->165)
t153 = m(2) + m(3);
t192 = m(1) * rSges(1,2);
t191 = m(2) * rSges(2,2);
t190 = m(3) * rSges(3,2);
t189 = pkin(1) * sin(pkin(7));
t138 = cos(pkin(7));
t188 = pkin(1) * t138;
t187 = pkin(2) * t138;
t151 = xDP(2);
t157 = 0.1e1 / pkin(3);
t175 = t151 * t157;
t139 = legFrame(3,3);
t115 = t139 + qJ(1,3);
t134 = qJ(1,3) + pkin(7);
t95 = t139 + t134;
t89 = qJ(3,3) + t95;
t79 = sin(t89);
t55 = -pkin(1) * sin(t115) - pkin(2) * sin(t95) - pkin(3) * t79;
t131 = pkin(7) + qJ(3,3);
t109 = sin(t131);
t145 = sin(qJ(3,3));
t75 = 0.1e1 / (pkin(1) * t109 + t145 * pkin(2));
t40 = t55 * t75 * t175;
t152 = xDP(1);
t174 = t152 * t157;
t82 = cos(t89);
t58 = -pkin(1) * cos(t115) - pkin(2) * cos(t95) - pkin(3) * t82;
t43 = t58 * t75 * t174;
t31 = t43 + t40;
t37 = (t151 * t79 + t152 * t82) * t75;
t22 = t37 + t31;
t186 = t22 * t31;
t140 = legFrame(2,3);
t116 = t140 + qJ(1,2);
t132 = pkin(7) + qJ(3,2);
t118 = qJ(1,2) + t132;
t90 = t140 + t118;
t80 = sin(t90);
t135 = qJ(1,2) + pkin(7);
t96 = t140 + t135;
t56 = -pkin(1) * sin(t116) - pkin(2) * sin(t96) - pkin(3) * t80;
t110 = sin(t132);
t146 = sin(qJ(3,2));
t76 = 0.1e1 / (pkin(1) * t110 + t146 * pkin(2));
t41 = t56 * t76 * t175;
t83 = cos(t90);
t59 = -pkin(1) * cos(t116) - pkin(2) * cos(t96) - pkin(3) * t83;
t44 = t59 * t76 * t174;
t32 = t44 + t41;
t38 = (t151 * t80 + t152 * t83) * t76;
t23 = t38 + t32;
t185 = t23 * t32;
t141 = legFrame(1,3);
t117 = t141 + qJ(1,1);
t136 = qJ(1,1) + pkin(7);
t97 = t141 + t136;
t91 = qJ(3,1) + t97;
t81 = sin(t91);
t57 = -pkin(1) * sin(t117) - pkin(2) * sin(t97) - pkin(3) * t81;
t133 = pkin(7) + qJ(3,1);
t111 = sin(t133);
t147 = sin(qJ(3,1));
t77 = 0.1e1 / (pkin(1) * t111 + t147 * pkin(2));
t42 = t57 * t77 * t175;
t84 = cos(t91);
t60 = -pkin(1) * cos(t117) - pkin(2) * cos(t97) - pkin(3) * t84;
t45 = t60 * t77 * t174;
t33 = t45 + t42;
t39 = (t151 * t81 + t152 * t84) * t77;
t24 = t39 + t33;
t184 = t24 * t33;
t182 = t157 * t55;
t181 = t157 * t56;
t180 = t157 * t57;
t179 = t157 * t58;
t178 = t157 * t59;
t177 = t157 * t60;
t127 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t92 = t127 * m(3) + Icges(3,3);
t176 = t157 * t92;
t173 = 0.2e1 * pkin(2) * pkin(3);
t159 = pkin(1) ^ 2;
t130 = pkin(2) ^ 2 + t159;
t172 = 0.2e1 * pkin(1);
t171 = pkin(1) * t190;
t170 = pkin(2) * t190;
t104 = pkin(2) + t188;
t112 = cos(t131);
t148 = cos(qJ(3,3));
t169 = pkin(1) * t112 + pkin(2) * t148;
t113 = cos(t132);
t149 = cos(qJ(3,2));
t168 = pkin(1) * t113 + pkin(2) * t149;
t114 = cos(t133);
t150 = cos(qJ(3,1));
t167 = pkin(1) * t114 + pkin(2) * t150;
t122 = sin(t140);
t125 = cos(t140);
t70 = -t122 * g(1) + t125 * g(2);
t73 = t125 * g(1) + t122 * g(2);
t166 = -cos(t118) * (rSges(3,1) * t70 - rSges(3,2) * t73) + sin(t118) * (rSges(3,1) * t73 + rSges(3,2) * t70);
t119 = qJ(1,3) + t131;
t121 = sin(t139);
t124 = cos(t139);
t69 = -t121 * g(1) + t124 * g(2);
t72 = t124 * g(1) + t121 * g(2);
t165 = -cos(t119) * (rSges(3,1) * t69 - rSges(3,2) * t72) + sin(t119) * (rSges(3,1) * t72 + rSges(3,2) * t69);
t120 = qJ(1,1) + t133;
t123 = sin(t141);
t126 = cos(t141);
t71 = -t123 * g(1) + t126 * g(2);
t74 = t126 * g(1) + t123 * g(2);
t164 = sin(t120) * (rSges(3,1) * t74 + rSges(3,2) * t71) - cos(t120) * (rSges(3,1) * t71 - rSges(3,2) * t74);
t163 = t169 * rSges(3,1);
t162 = t168 * rSges(3,1);
t161 = t167 * rSges(3,1);
t108 = m(2) * rSges(2,1) + m(3) * pkin(2);
t160 = Icges(1,3) + Icges(2,3) + Icges(3,3) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + (rSges(2,1) ^ 2 + t159 + (-0.2e1 * t189 + rSges(2,2)) * rSges(2,2)) * m(2) + 0.2e1 * t108 * t188;
t156 = pkin(3) ^ 2;
t144 = xDDP(1);
t143 = xDDP(2);
t107 = t147 * t170;
t106 = t146 * t170;
t105 = t145 * t170;
t93 = t130 + t127;
t88 = m(1) * rSges(1,1) + t153 * pkin(1);
t87 = t111 * t171;
t86 = t110 * t171;
t85 = t109 * t171;
t68 = t104 * rSges(3,1) - rSges(3,2) * t189;
t67 = rSges(3,1) * t189 + t104 * rSges(3,2);
t48 = Icges(3,3) - t107 - t87 + (t127 + t161) * m(3);
t47 = Icges(3,3) - t106 - t86 + (t127 + t162) * m(3);
t46 = Icges(3,3) - t105 - t85 + (t127 + t163) * m(3);
t36 = -0.2e1 * t107 - 0.2e1 * t87 + (0.2e1 * t161 + t93) * m(3) + t160;
t35 = -0.2e1 * t106 - 0.2e1 * t86 + (0.2e1 * t162 + t93) * m(3) + t160;
t34 = -0.2e1 * t105 - 0.2e1 * t85 + (0.2e1 * t163 + t93) * m(3) + t160;
t30 = (t60 * t176 + t48 * t84) * t77;
t29 = (t59 * t176 + t47 * t83) * t76;
t28 = (t58 * t176 + t46 * t82) * t75;
t27 = (t57 * t176 + t48 * t81) * t77;
t26 = (t56 * t176 + t47 * t80) * t76;
t25 = (t55 * t176 + t46 * t79) * t75;
t21 = t45 / 0.2e1 + t42 / 0.2e1 + t39;
t20 = t44 / 0.2e1 + t41 / 0.2e1 + t38;
t19 = t43 / 0.2e1 + t40 / 0.2e1 + t37;
t18 = (t48 * t177 + t36 * t84) * t77;
t17 = (t47 * t178 + t35 * t83) * t76;
t16 = (t46 * t179 + t34 * t82) * t75;
t15 = (t48 * t180 + t36 * t81) * t77;
t14 = (t47 * t181 + t35 * t80) * t76;
t13 = (t46 * t182 + t34 * t79) * t75;
t12 = (-pkin(3) * t184 + (-pkin(3) * t24 - t167 * t39) * t39) * t77;
t11 = (-pkin(3) * t185 + (-pkin(3) * t23 - t168 * t38) * t38) * t76;
t10 = (-pkin(3) * t186 + (-pkin(3) * t22 - t169 * t37) * t37) * t75;
t9 = (t21 * t150 * t173 + t39 * t130 + t24 * t156 + (pkin(3) * t114 * t21 + t39 * t187) * t172) * t157 * t77 * t39 + (t104 * t150 - t147 * t189 + pkin(3)) / (t104 * t147 + t150 * t189) * t184;
t8 = (t20 * t149 * t173 + t38 * t130 + t23 * t156 + (pkin(3) * t113 * t20 + t38 * t187) * t172) * t157 * t76 * t38 + (t104 * t149 - t146 * t189 + pkin(3)) / (t104 * t146 + t149 * t189) * t185;
t7 = (t19 * t148 * t173 + t37 * t130 + t22 * t156 + (pkin(3) * t112 * t19 + t37 * t187) * t172) * t157 * t75 * t37 + (t104 * t148 - t145 * t189 + pkin(3)) / (t104 * t145 + t148 * t189) * t186;
t6 = -t48 * t12 - t92 * t9 + ((pkin(2) * (t147 * rSges(3,1) + rSges(3,2) * t150) + (rSges(3,1) * t111 + rSges(3,2) * t114) * pkin(1)) * t39 ^ 2 + t164) * m(3);
t5 = -t47 * t11 - t92 * t8 + ((pkin(2) * (t146 * rSges(3,1) + rSges(3,2) * t149) + (rSges(3,1) * t110 + rSges(3,2) * t113) * pkin(1)) * t38 ^ 2 + t166) * m(3);
t4 = -t46 * t10 - t92 * t7 + ((pkin(2) * (t145 * rSges(3,1) + rSges(3,2) * t148) + (rSges(3,1) * t109 + rSges(3,2) * t112) * pkin(1)) * t37 ^ 2 + t165) * m(3);
t3 = -t36 * t12 - t48 * t9 + (-t71 * t108 + t74 * t191) * cos(t136) + (t108 * t74 + t71 * t191) * sin(t136) + (t74 * t192 - t71 * t88) * cos(qJ(1,1)) + sin(qJ(1,1)) * (t71 * t192 + t88 * t74) + (-0.2e1 * (t147 * t68 + t67 * t150) * t21 * t33 + t164) * m(3);
t2 = -t35 * t11 - t47 * t8 + (-t70 * t108 + t73 * t191) * cos(t135) + (t108 * t73 + t70 * t191) * sin(t135) + (t73 * t192 - t70 * t88) * cos(qJ(1,2)) + sin(qJ(1,2)) * (t70 * t192 + t88 * t73) + (-0.2e1 * (t146 * t68 + t67 * t149) * t20 * t32 + t166) * m(3);
t1 = -t34 * t10 - t46 * t7 + (-t69 * t108 + t72 * t191) * cos(t134) + (t108 * t72 + t69 * t191) * sin(t134) + (t72 * t192 - t69 * t88) * cos(qJ(1,3)) + sin(qJ(1,3)) * (t69 * t192 + t88 * t72) + (-0.2e1 * (t145 * t68 + t67 * t148) * t19 * t31 + t165) * m(3);
t49 = [(-g(1) + t144) * m(4) + ((t30 * t177 + t18 * t84) * t144 + (t18 * t81 + t30 * t180) * t143 + t84 * t3 + t6 * t177) * t77 + ((t17 * t83 + t29 * t178) * t144 + (t17 * t80 + t29 * t181) * t143 + t83 * t2 + t5 * t178) * t76 + ((t16 * t82 + t28 * t179) * t144 + (t16 * t79 + t28 * t182) * t143 + t82 * t1 + t4 * t179) * t75; (-g(2) + t143) * m(4) + ((t15 * t84 + t27 * t177) * t144 + (t15 * t81 + t27 * t180) * t143 + t81 * t3 + t6 * t180) * t77 + ((t14 * t83 + t26 * t178) * t144 + (t14 * t80 + t26 * t181) * t143 + t80 * t2 + t5 * t181) * t76 + ((t13 * t82 + t25 * t179) * t144 + (t13 * t79 + t25 * t182) * t143 + t79 * t1 + t4 * t182) * t75; (-g(3) + xDDP(3)) * (m(4) + 0.3e1 * t153);];
tauX  = t49;
