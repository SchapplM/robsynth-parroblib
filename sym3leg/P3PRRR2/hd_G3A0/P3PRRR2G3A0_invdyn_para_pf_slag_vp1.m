% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRR2G3A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
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
% Datum: 2020-03-09 21:20
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRR2G3A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G3A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR2G3A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRR2G3A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G3A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR2G3A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G3A0_invdyn_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR2G3A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRR2G3A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRR2G3A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G3A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G3A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:20:00
% EndTime: 2020-03-09 21:20:01
% DurationCPUTime: 0.93s
% Computational Cost: add. (2475->180), mult. (2338->323), div. (1374->5), fcn. (1686->36), ass. (0->143)
t110 = qJ(2,1) + qJ(3,1);
t113 = legFrame(1,2);
t101 = cos(t113);
t98 = sin(t113);
t66 = t98 * g(1) + t101 * g(2);
t69 = t101 * g(1) - t98 * g(2);
t169 = (rSges(3,1) * t69 + rSges(3,2) * t66) * cos(t110) + (rSges(3,1) * t66 - rSges(3,2) * t69) * sin(t110);
t109 = qJ(2,2) + qJ(3,2);
t112 = legFrame(2,2);
t100 = cos(t112);
t97 = sin(t112);
t65 = t97 * g(1) + t100 * g(2);
t68 = t100 * g(1) - t97 * g(2);
t168 = (rSges(3,1) * t68 + rSges(3,2) * t65) * cos(t109) + (rSges(3,1) * t65 - rSges(3,2) * t68) * sin(t109);
t108 = qJ(2,3) + qJ(3,3);
t111 = legFrame(3,2);
t96 = sin(t111);
t99 = cos(t111);
t64 = t96 * g(1) + t99 * g(2);
t67 = t99 * g(1) - t96 * g(2);
t167 = (rSges(3,1) * t67 + rSges(3,2) * t64) * cos(t108) + (rSges(3,1) * t64 - rSges(3,2) * t67) * sin(t108);
t166 = (pkin(1) * m(3));
t120 = cos(qJ(3,3));
t165 = t120 * pkin(1);
t121 = cos(qJ(3,2));
t164 = t121 * pkin(1);
t122 = cos(qJ(3,1));
t163 = t122 * pkin(1);
t117 = sin(qJ(3,3));
t105 = 0.1e1 / t117;
t124 = xDP(1);
t128 = 0.1e1 / pkin(2);
t130 = 0.1e1 / pkin(1);
t140 = t128 * t130;
t132 = t124 * t140;
t87 = -t111 + qJ(2,3);
t83 = qJ(3,3) + t87;
t74 = sin(t83);
t55 = pkin(2) * t74 + pkin(1) * sin(t87);
t46 = t55 * t105 * t132;
t123 = xDP(2);
t133 = t123 * t140;
t77 = cos(t83);
t58 = -pkin(2) * t77 - pkin(1) * cos(t87);
t49 = t58 * t105 * t133;
t31 = t46 + t49;
t141 = t124 * t130;
t142 = t123 * t130;
t34 = (-t141 * t74 + t142 * t77) * t105;
t22 = t34 + t31;
t162 = t22 * t31;
t118 = sin(qJ(3,2));
t106 = 0.1e1 / t118;
t88 = -t112 + qJ(2,2);
t84 = qJ(3,2) + t88;
t75 = sin(t84);
t56 = pkin(2) * t75 + pkin(1) * sin(t88);
t47 = t56 * t106 * t132;
t78 = cos(t84);
t59 = -pkin(2) * t78 - pkin(1) * cos(t88);
t50 = t59 * t106 * t133;
t32 = t47 + t50;
t35 = (-t141 * t75 + t142 * t78) * t106;
t23 = t35 + t32;
t161 = t23 * t32;
t119 = sin(qJ(3,1));
t107 = 0.1e1 / t119;
t89 = -t113 + qJ(2,1);
t85 = qJ(3,1) + t89;
t76 = sin(t85);
t57 = pkin(2) * t76 + pkin(1) * sin(t89);
t48 = t57 * t107 * t132;
t79 = cos(t85);
t60 = -pkin(2) * t79 - pkin(1) * cos(t89);
t51 = t60 * t107 * t133;
t33 = t48 + t51;
t36 = (-t141 * t76 + t142 * t79) * t107;
t24 = t36 + t33;
t160 = t24 * t33;
t152 = t128 * t55;
t151 = t128 * t56;
t150 = t128 * t57;
t149 = t128 * t58;
t148 = t128 * t59;
t147 = t128 * t60;
t102 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t73 = t102 * m(3) + Icges(3,3);
t146 = t128 * t73;
t145 = t105 * t130;
t144 = t106 * t130;
t143 = t107 * t130;
t139 = -2 * t166;
t138 = 0.2e1 * pkin(1) * pkin(2);
t137 = rSges(3,2) * t166;
t136 = rSges(3,1) * t165;
t135 = rSges(3,1) * t164;
t134 = rSges(3,1) * t163;
t131 = (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(2,3) + Icges(3,3);
t129 = pkin(1) ^ 2;
t127 = pkin(2) ^ 2;
t116 = xDDP(1);
t115 = xDDP(2);
t86 = t129 + t102;
t82 = t119 * t137;
t81 = t118 * t137;
t80 = t117 * t137;
t72 = t119 * rSges(3,1) + t122 * rSges(3,2);
t71 = t118 * rSges(3,1) + t121 * rSges(3,2);
t70 = t117 * rSges(3,1) + t120 * rSges(3,2);
t54 = Icges(3,3) - t82 + (t102 + t134) * m(3);
t53 = Icges(3,3) - t81 + (t102 + t135) * m(3);
t52 = Icges(3,3) - t80 + (t102 + t136) * m(3);
t45 = -0.2e1 * t82 + (t86 + 0.2e1 * t134) * m(3) + t131;
t44 = -0.2e1 * t81 + (t86 + 0.2e1 * t135) * m(3) + t131;
t43 = -0.2e1 * t80 + (t86 + 0.2e1 * t136) * m(3) + t131;
t30 = (t60 * t146 + t54 * t79) * t143;
t29 = (t59 * t146 + t53 * t78) * t144;
t28 = (t58 * t146 + t52 * t77) * t145;
t27 = (t57 * t146 - t54 * t76) * t143;
t26 = (t56 * t146 - t53 * t75) * t144;
t25 = (t55 * t146 - t52 * t74) * t145;
t21 = t48 / 0.2e1 + t51 / 0.2e1 + t36;
t20 = t47 / 0.2e1 + t50 / 0.2e1 + t35;
t19 = t46 / 0.2e1 + t49 / 0.2e1 + t34;
t18 = (t54 * t147 + t45 * t79) * t143;
t17 = (t53 * t148 + t44 * t78) * t144;
t16 = (t52 * t149 + t43 * t77) * t145;
t15 = (t54 * t150 - t45 * t76) * t143;
t14 = (t53 * t151 - t44 * t75) * t144;
t13 = (t52 * t152 - t43 * t74) * t145;
t12 = ((-pkin(2) * t24 - t36 * t163) * t36 - pkin(2) * t160) * t143;
t11 = ((-pkin(2) * t23 - t35 * t164) * t35 - pkin(2) * t161) * t144;
t10 = ((-pkin(2) * t22 - t34 * t165) * t34 - pkin(2) * t162) * t145;
t9 = ((t21 * t122 * t138 + t24 * t127 + t129 * t36) * t128 * t36 + (pkin(2) + t163) * t160) * t143;
t8 = ((t20 * t121 * t138 + t23 * t127 + t129 * t35) * t128 * t35 + (pkin(2) + t164) * t161) * t144;
t7 = ((t19 * t120 * t138 + t22 * t127 + t129 * t34) * t128 * t34 + (pkin(2) + t165) * t162) * t145;
t6 = -t54 * t12 - t73 * t9 + (pkin(1) * t36 ^ 2 * t72 + t169) * m(3);
t5 = -t53 * t11 - t73 * t8 + (pkin(1) * t35 ^ 2 * t71 + t168) * m(3);
t4 = -t52 * t10 - t73 * t7 + (pkin(1) * t34 ^ 2 * t70 + t167) * m(3);
t3 = -t45 * t12 - t54 * t9 + t72 * t33 * t21 * t139 + (t69 * t166 + m(2) * (rSges(2,1) * t69 + rSges(2,2) * t66)) * cos(qJ(2,1)) - (-t66 * t166 + m(2) * (-rSges(2,1) * t66 + rSges(2,2) * t69)) * sin(qJ(2,1)) + t169 * m(3);
t2 = -t44 * t11 - t53 * t8 + t71 * t32 * t20 * t139 + (t68 * t166 + m(2) * (rSges(2,1) * t68 + rSges(2,2) * t65)) * cos(qJ(2,2)) - (-t65 * t166 + m(2) * (-rSges(2,1) * t65 + rSges(2,2) * t68)) * sin(qJ(2,2)) + t168 * m(3);
t1 = -t43 * t10 - t52 * t7 + t70 * t31 * t19 * t139 + (t67 * t166 + m(2) * (rSges(2,1) * t67 + rSges(2,2) * t64)) * cos(qJ(2,3)) - (-t64 * t166 + m(2) * (-rSges(2,1) * t64 + rSges(2,2) * t67)) * sin(qJ(2,3)) + t167 * m(3);
t37 = [(-g(1) + t116) * m(4) + (((-t15 * t76 + t27 * t150) * t116 + (t27 * t147 + t15 * t79) * t115 - t76 * t3 + t6 * t150) * t107 + ((-t14 * t75 + t26 * t151) * t116 + (t14 * t78 + t26 * t148) * t115 - t75 * t2 + t5 * t151) * t106 + ((-t13 * t74 + t25 * t152) * t116 + (t13 * t77 + t25 * t149) * t115 - t74 * t1 + t4 * t152) * t105) * t130; (-g(2) + t115) * m(4) + (((t30 * t150 - t18 * t76) * t116 + (t30 * t147 + t18 * t79) * t115 + t79 * t3 + t6 * t147) * t107 + ((t29 * t151 - t17 * t75) * t116 + (t29 * t148 + t17 * t78) * t115 + t78 * t2 + t5 * t148) * t106 + ((t28 * t152 - t16 * t74) * t116 + (t28 * t149 + t16 * t77) * t115 + t77 * t1 + t4 * t149) * t105) * t130; (-g(3) + xDDP(3)) * ((3 * m(1)) + 0.3e1 * m(2) + 0.3e1 * m(3) + m(4));];
tauX  = t37;
