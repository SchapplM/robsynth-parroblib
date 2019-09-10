% Calculate vector of inverse dynamics forces for parallel robot
% P3RPR1G1P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [2x3]
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% m [3x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [3x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:58
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPR1G1P1A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(2+1,1),zeros(2+1,3),zeros(2+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_slag_vp1: qJ has to be [2x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_slag_vp1: Icges has to be [3x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:58:03
% EndTime: 2019-05-03 14:58:04
% DurationCPUTime: 1.02s
% Computational Cost: add. (3584->221), mult. (5291->378), div. (642->3), fcn. (4010->14), ass. (0->150)
t161 = 2 * rSges(2,3);
t132 = 0.1e1 / qJ(2,3);
t117 = sin(qJ(1,3));
t120 = cos(qJ(1,3));
t127 = pkin(1) + pkin(2);
t90 = -qJ(2,3) * t120 + t117 * t127;
t93 = qJ(2,3) * t117 + t120 * t127;
t111 = legFrame(3,3);
t96 = sin(t111);
t99 = cos(t111);
t55 = t90 * t99 + t93 * t96;
t58 = -t90 * t96 + t93 * t99;
t124 = xDP(3);
t125 = xDP(2);
t128 = xP(3);
t103 = sin(t128);
t104 = cos(t128);
t137 = koppelP(3,2);
t140 = koppelP(3,1);
t81 = -t103 * t137 + t104 * t140;
t64 = t124 * t81 + t125;
t126 = xDP(1);
t78 = t103 * t140 + t104 * t137;
t67 = -t124 * t78 + t126;
t19 = (t55 * t64 + t58 * t67) * t132;
t160 = 0.2e1 * t19;
t134 = 0.1e1 / qJ(2,2);
t112 = legFrame(2,3);
t100 = cos(t112);
t118 = sin(qJ(1,2));
t121 = cos(qJ(1,2));
t91 = -qJ(2,2) * t121 + t118 * t127;
t94 = qJ(2,2) * t118 + t121 * t127;
t97 = sin(t112);
t56 = t100 * t91 + t94 * t97;
t59 = t100 * t94 - t91 * t97;
t138 = koppelP(2,2);
t141 = koppelP(2,1);
t82 = -t103 * t138 + t104 * t141;
t65 = t124 * t82 + t125;
t79 = t103 * t141 + t104 * t138;
t68 = -t124 * t79 + t126;
t20 = (t56 * t65 + t59 * t68) * t134;
t159 = 0.2e1 * t20;
t136 = 0.1e1 / qJ(2,1);
t113 = legFrame(1,3);
t101 = cos(t113);
t119 = sin(qJ(1,1));
t122 = cos(qJ(1,1));
t92 = -qJ(2,1) * t122 + t119 * t127;
t95 = qJ(2,1) * t119 + t122 * t127;
t98 = sin(t113);
t57 = t101 * t92 + t95 * t98;
t60 = t101 * t95 - t92 * t98;
t139 = koppelP(1,2);
t142 = koppelP(1,1);
t83 = -t103 * t139 + t104 * t142;
t66 = t124 * t83 + t125;
t80 = t103 * t142 + t104 * t139;
t69 = -t124 * t80 + t126;
t21 = (t57 * t66 + t60 * t69) * t136;
t158 = 0.2e1 * t21;
t157 = 0.2e1 * t127;
t123 = pkin(1) + rSges(2,1);
t156 = m(2) * t123;
t155 = m(2) * t132;
t154 = m(2) * t134;
t153 = m(2) * t136;
t70 = t117 * t99 + t120 * t96;
t71 = -t117 * t96 + t120 * t99;
t31 = (t64 * t70 + t67 * t71) * t132;
t152 = t132 * t31;
t72 = t100 * t118 + t121 * t97;
t73 = t100 * t121 - t118 * t97;
t32 = (t65 * t72 + t68 * t73) * t134;
t151 = t134 * t32;
t74 = t101 * t119 + t122 * t98;
t75 = t101 * t122 - t119 * t98;
t33 = (t66 * t74 + t69 * t75) * t136;
t150 = t136 * t33;
t146 = (pkin(1) ^ 2);
t149 = -t146 + (-2 * pkin(1) - pkin(2)) * pkin(2);
t148 = rSges(2,3) ^ 2 + t146 + (2 * pkin(1) + rSges(2,1)) * rSges(2,1);
t147 = (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + Icges(2,2) + Icges(1,3);
t135 = qJ(2,1) ^ 2;
t133 = qJ(2,2) ^ 2;
t131 = qJ(2,3) ^ 2;
t130 = rSges(3,1);
t129 = rSges(3,2);
t116 = xDDP(1);
t115 = xDDP(2);
t114 = xDDP(3);
t110 = rSges(2,3) + qJ(2,1);
t109 = rSges(2,3) + qJ(2,2);
t108 = rSges(2,3) + qJ(2,3);
t106 = t124 ^ 2;
t89 = g(1) * t101 + g(2) * t98;
t88 = g(1) * t100 + g(2) * t97;
t87 = g(1) * t99 + g(2) * t96;
t86 = -g(1) * t98 + g(2) * t101;
t85 = -g(1) * t97 + g(2) * t100;
t84 = -g(1) * t96 + g(2) * t99;
t77 = -t103 * t129 + t104 * t130;
t76 = t103 * t130 + t104 * t129;
t63 = (qJ(2,1) * t161 + t135 + t148) * m(2) + t147;
t62 = (qJ(2,2) * t161 + t133 + t148) * m(2) + t147;
t61 = (qJ(2,3) * t161 + t131 + t148) * m(2) + t147;
t54 = -t106 * t83 - t114 * t80 + t116;
t53 = -t106 * t82 - t114 * t79 + t116;
t52 = -t106 * t81 - t114 * t78 + t116;
t51 = -t106 * t80 + t114 * t83 + t115;
t50 = -t106 * t79 + t114 * t82 + t115;
t49 = -t106 * t78 + t114 * t81 + t115;
t42 = (-t123 * t75 + t60) * t153;
t41 = (-t123 * t74 + t57) * t153;
t40 = (-t123 * t73 + t59) * t154;
t39 = (-t123 * t72 + t56) * t154;
t38 = (-t123 * t71 + t58) * t155;
t37 = (-t123 * t70 + t55) * t155;
t36 = (t74 * t83 - t75 * t80) * t136;
t35 = (t72 * t82 - t73 * t79) * t134;
t34 = (t70 * t81 - t71 * t78) * t132;
t30 = (-t156 * t60 + t63 * t75) * t136;
t29 = (-t156 * t57 + t63 * t74) * t136;
t28 = (-t156 * t59 + t62 * t73) * t134;
t27 = (-t156 * t56 + t62 * t72) * t134;
t26 = (-t156 * t58 + t61 * t71) * t132;
t25 = (-t156 * t55 + t61 * t70) * t132;
t24 = (t57 * t83 - t60 * t80) * t136;
t23 = (t56 * t82 - t59 * t79) * t134;
t22 = (t55 * t81 - t58 * t78) * t132;
t18 = (-t123 * t36 + t24) * m(2);
t17 = (-t123 * t35 + t23) * m(2);
t16 = (-t123 * t34 + t22) * m(2);
t15 = -t156 * t24 + t36 * t63;
t14 = -t156 * t23 + t35 * t62;
t13 = -t156 * t22 + t34 * t61;
t12 = (-t127 * t33 + t158) * t150;
t11 = (-t127 * t32 + t159) * t151;
t10 = (-t127 * t31 + t160) * t152;
t9 = (t21 * t157 + (-t135 + t149) * t33) * t150;
t8 = (t20 * t157 + (-t133 + t149) * t32) * t151;
t7 = (t19 * t157 + (-t131 + t149) * t31) * t152;
t6 = (-t110 * t33 ^ 2 - t119 * t89 + t12 * t123 + t122 * t86 - t9) * m(2);
t5 = (-t109 * t32 ^ 2 + t11 * t123 - t118 * t88 + t121 * t85 - t8) * m(2);
t4 = (-t108 * t31 ^ 2 + t10 * t123 - t117 * t87 + t120 * t84 - t7) * m(2);
t3 = -t63 * t12 + ((-rSges(1,1) * t86 + rSges(1,2) * t89) * t122 + t119 * (rSges(1,1) * t89 + rSges(1,2) * t86)) * m(1) + (t123 * t9 + t33 * t110 * t158 + (-t110 * t89 - t123 * t86) * t122 + t119 * (-t110 * t86 + t123 * t89)) * m(2);
t2 = -t62 * t11 + ((-rSges(1,1) * t85 + rSges(1,2) * t88) * t121 + t118 * (rSges(1,1) * t88 + rSges(1,2) * t85)) * m(1) + (t123 * t8 + t32 * t109 * t159 + (-t109 * t88 - t123 * t85) * t121 + t118 * (-t109 * t85 + t123 * t88)) * m(2);
t1 = -t61 * t10 + ((-rSges(1,1) * t84 + rSges(1,2) * t87) * t120 + t117 * (rSges(1,1) * t87 + rSges(1,2) * t84)) * m(1) + (t123 * t7 + t31 * t108 * t160 + (-t108 * t87 - t123 * t84) * t120 + t117 * (-t108 * t84 + t123 * t87)) * m(2);
t43 = [(-t106 * t77 - t76 * t114 - g(1) + t116) * m(3) + ((t30 * t75 + t42 * t60) * t54 + (t30 * t74 + t42 * t57) * t51 + t75 * t3 + t60 * t6) * t136 + ((t28 * t73 + t40 * t59) * t53 + (t28 * t72 + t40 * t56) * t50 + t73 * t2 + t59 * t5) * t134 + ((t26 * t71 + t38 * t58) * t52 + (t26 * t70 + t38 * t55) * t49 + t71 * t1 + t58 * t4) * t132; (-t106 * t76 + t114 * t77 - g(2) + t115) * m(3) + ((t29 * t75 + t41 * t60) * t54 + (t29 * t74 + t41 * t57) * t51 + t74 * t3 + t57 * t6) * t136 + ((t27 * t73 + t39 * t59) * t53 + (t27 * t72 + t39 * t56) * t50 + t72 * t2 + t56 * t5) * t134 + ((t25 * t71 + t37 * t58) * t52 + (t25 * t70 + t37 * t55) * t49 + t70 * t1 + t55 * t4) * t132; Icges(3,3) * t114 + t34 * t1 + t35 * t2 + t22 * t4 + t23 * t5 + t24 * t6 + t36 * t3 + (-t76 * t116 + t77 * t115 + (t129 ^ 2 + t130 ^ 2) * t114 + (g(1) * t130 + g(2) * t129) * t103 + (g(1) * t129 - g(2) * t130) * t104) * m(3) + ((t15 * t75 + t18 * t60) * t54 + (t15 * t74 + t18 * t57) * t51) * t136 + ((t14 * t73 + t17 * t59) * t53 + (t14 * t72 + t17 * t56) * t50) * t134 + ((t13 * t71 + t16 * t58) * t52 + (t13 * t70 + t16 * t55) * t49) * t132;];
tauX  = t43;
