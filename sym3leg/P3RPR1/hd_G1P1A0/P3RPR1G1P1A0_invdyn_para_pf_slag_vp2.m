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
% mrSges [3x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [3x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
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

function tauX = P3RPR1G1P1A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(2+1,1),zeros(2+1,3),zeros(2+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_slag_vp2: qJ has to be [2x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_slag_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_slag_vp2: Ifges has to be [3x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:58:06
% EndTime: 2019-05-03 14:58:07
% DurationCPUTime: 0.92s
% Computational Cost: add. (3438->207), mult. (5038->352), div. (642->3), fcn. (3938->14), ass. (0->149)
t157 = 2 * mrSges(2,3);
t135 = 0.1e1 / qJ(2,3);
t115 = legFrame(3,3);
t103 = sin(t115);
t106 = cos(t115);
t121 = sin(qJ(1,3));
t124 = cos(qJ(1,3));
t130 = pkin(1) + pkin(2);
t90 = -qJ(2,3) * t124 + t121 * t130;
t93 = qJ(2,3) * t121 + t124 * t130;
t55 = t103 * t93 + t106 * t90;
t58 = -t103 * t90 + t106 * t93;
t127 = xDP(3);
t128 = xDP(2);
t131 = xP(3);
t110 = sin(t131);
t111 = cos(t131);
t140 = koppelP(3,2);
t143 = koppelP(3,1);
t81 = -t110 * t140 + t111 * t143;
t61 = t127 * t81 + t128;
t129 = xDP(1);
t78 = t110 * t143 + t111 * t140;
t64 = -t127 * t78 + t129;
t19 = (t55 * t61 + t58 * t64) * t135;
t156 = 0.2e1 * t19;
t137 = 0.1e1 / qJ(2,2);
t116 = legFrame(2,3);
t104 = sin(t116);
t107 = cos(t116);
t122 = sin(qJ(1,2));
t125 = cos(qJ(1,2));
t91 = -qJ(2,2) * t125 + t122 * t130;
t94 = qJ(2,2) * t122 + t125 * t130;
t56 = t104 * t94 + t107 * t91;
t59 = -t104 * t91 + t107 * t94;
t141 = koppelP(2,2);
t144 = koppelP(2,1);
t82 = -t110 * t141 + t111 * t144;
t62 = t127 * t82 + t128;
t79 = t110 * t144 + t111 * t141;
t65 = -t127 * t79 + t129;
t20 = (t56 * t62 + t59 * t65) * t137;
t155 = 0.2e1 * t20;
t139 = 0.1e1 / qJ(2,1);
t117 = legFrame(1,3);
t105 = sin(t117);
t108 = cos(t117);
t123 = sin(qJ(1,1));
t126 = cos(qJ(1,1));
t92 = -qJ(2,1) * t126 + t123 * t130;
t95 = qJ(2,1) * t123 + t126 * t130;
t57 = t105 * t95 + t108 * t92;
t60 = -t105 * t92 + t108 * t95;
t142 = koppelP(1,2);
t145 = koppelP(1,1);
t83 = -t110 * t142 + t111 * t145;
t63 = t127 * t83 + t128;
t80 = t110 * t145 + t111 * t142;
t66 = -t127 * t80 + t129;
t21 = (t57 * t63 + t60 * t66) * t139;
t154 = 0.2e1 * t21;
t153 = 0.2e1 * t130;
t67 = t103 * t124 + t106 * t121;
t68 = -t103 * t121 + t106 * t124;
t25 = (t61 * t67 + t64 * t68) * t135;
t152 = t135 * t25;
t69 = t104 * t125 + t107 * t122;
t70 = -t104 * t122 + t107 * t125;
t26 = (t62 * t69 + t65 * t70) * t137;
t151 = t137 * t26;
t71 = t105 * t126 + t108 * t123;
t72 = -t105 * t123 + t108 * t126;
t27 = (t63 * t71 + t66 * t72) * t139;
t150 = t139 * t27;
t109 = m(2) * pkin(1) + mrSges(2,1);
t147 = (pkin(1) ^ 2);
t149 = -t147 + (-2 * pkin(1) - pkin(2)) * pkin(2);
t102 = m(2) * qJ(2,1) + mrSges(2,3);
t101 = m(2) * qJ(2,2) + mrSges(2,3);
t100 = m(2) * qJ(2,3) + mrSges(2,3);
t148 = 2 * pkin(1) * mrSges(2,1) + Ifges(2,2) + Ifges(1,3);
t138 = qJ(2,1) ^ 2;
t136 = qJ(2,2) ^ 2;
t134 = qJ(2,3) ^ 2;
t133 = mrSges(3,1);
t132 = mrSges(3,2);
t120 = xDDP(1);
t119 = xDDP(2);
t118 = xDDP(3);
t113 = t127 ^ 2;
t99 = mrSges(1,1) + t109;
t98 = -mrSges(1,2) + t102;
t97 = -mrSges(1,2) + t101;
t96 = -mrSges(1,2) + t100;
t89 = g(1) * t108 + g(2) * t105;
t88 = g(1) * t107 + g(2) * t104;
t87 = g(1) * t106 + g(2) * t103;
t86 = -g(1) * t105 + g(2) * t108;
t85 = -g(1) * t104 + g(2) * t107;
t84 = -g(1) * t103 + g(2) * t106;
t77 = m(2) * (t138 + t147) + qJ(2,1) * t157 + t148;
t76 = m(2) * (t136 + t147) + qJ(2,2) * t157 + t148;
t75 = m(2) * (t134 + t147) + qJ(2,3) * t157 + t148;
t74 = -t110 * t132 + t111 * t133;
t73 = t110 * t133 + t111 * t132;
t54 = -t113 * t83 - t118 * t80 + t120;
t53 = -t113 * t82 - t118 * t79 + t120;
t52 = -t113 * t81 - t118 * t78 + t120;
t51 = -t113 * t80 + t118 * t83 + t119;
t50 = -t113 * t79 + t118 * t82 + t119;
t49 = -t113 * t78 + t118 * t81 + t119;
t42 = (m(2) * t60 - t109 * t72) * t139;
t41 = (m(2) * t57 - t109 * t71) * t139;
t40 = (m(2) * t59 - t109 * t70) * t137;
t39 = (m(2) * t56 - t109 * t69) * t137;
t38 = (m(2) * t58 - t109 * t68) * t135;
t37 = (m(2) * t55 - t109 * t67) * t135;
t36 = (t71 * t83 - t72 * t80) * t139;
t35 = (t69 * t82 - t70 * t79) * t137;
t34 = (t67 * t81 - t68 * t78) * t135;
t33 = (-t109 * t60 + t72 * t77) * t139;
t32 = (-t109 * t57 + t71 * t77) * t139;
t31 = (-t109 * t59 + t70 * t76) * t137;
t30 = (-t109 * t56 + t69 * t76) * t137;
t29 = (-t109 * t58 + t68 * t75) * t135;
t28 = (-t109 * t55 + t67 * t75) * t135;
t24 = (t57 * t83 - t60 * t80) * t139;
t23 = (t56 * t82 - t59 * t79) * t137;
t22 = (t55 * t81 - t58 * t78) * t135;
t18 = m(2) * t24 - t109 * t36;
t17 = m(2) * t23 - t109 * t35;
t16 = m(2) * t22 - t109 * t34;
t15 = -t109 * t24 + t36 * t77;
t14 = -t109 * t23 + t35 * t76;
t13 = -t109 * t22 + t34 * t75;
t12 = (-t130 * t27 + t154) * t150;
t11 = (-t130 * t26 + t155) * t151;
t10 = (-t130 * t25 + t156) * t152;
t9 = (t21 * t153 + (-t138 + t149) * t27) * t150;
t8 = (t20 * t153 + (-t136 + t149) * t26) * t151;
t7 = (t19 * t153 + (-t134 + t149) * t25) * t152;
t6 = -t27 ^ 2 * t102 + t109 * t12 + (-t123 * t89 + t126 * t86 - t9) * m(2);
t5 = -t26 ^ 2 * t101 + t109 * t11 + (-t122 * t88 + t125 * t85 - t8) * m(2);
t4 = -t25 ^ 2 * t100 + t109 * t10 + (-t121 * t87 + t124 * t84 - t7) * m(2);
t3 = -t77 * t12 + t109 * t9 + t27 * t102 * t154 + (-t86 * t99 - t89 * t98) * t126 - (t86 * t98 - t89 * t99) * t123;
t2 = -t76 * t11 + t109 * t8 + t26 * t101 * t155 + (-t85 * t99 - t88 * t97) * t125 - (t85 * t97 - t88 * t99) * t122;
t1 = -t75 * t10 + t109 * t7 + t25 * t100 * t156 + (-t84 * t99 - t87 * t96) * t124 - (t84 * t96 - t87 * t99) * t121;
t43 = [-t113 * t74 - t73 * t118 + (t120 - g(1)) * m(3) + ((t33 * t72 + t42 * t60) * t54 + (t33 * t71 + t42 * t57) * t51 + t72 * t3 + t60 * t6) * t139 + ((t31 * t70 + t40 * t59) * t53 + (t31 * t69 + t40 * t56) * t50 + t70 * t2 + t59 * t5) * t137 + ((t29 * t68 + t38 * t58) * t52 + (t29 * t67 + t38 * t55) * t49 + t68 * t1 + t58 * t4) * t135; -t113 * t73 + t74 * t118 + (t119 - g(2)) * m(3) + ((t32 * t72 + t41 * t60) * t54 + (t32 * t71 + t41 * t57) * t51 + t71 * t3 + t57 * t6) * t139 + ((t30 * t70 + t39 * t59) * t53 + (t30 * t69 + t39 * t56) * t50 + t69 * t2 + t56 * t5) * t137 + ((t28 * t68 + t37 * t58) * t52 + (t28 * t67 + t37 * t55) * t49 + t67 * t1 + t55 * t4) * t135; t36 * t3 + t24 * t6 + t35 * t2 + t23 * t5 + t34 * t1 + t22 * t4 - t73 * t120 + t74 * t119 + Ifges(3,3) * t118 - (-g(1) * t133 - g(2) * t132) * t110 + t111 * (g(1) * t132 - g(2) * t133) + ((t15 * t72 + t18 * t60) * t54 + (t15 * t71 + t18 * t57) * t51) * t139 + ((t14 * t70 + t17 * t59) * t53 + (t14 * t69 + t17 * t56) * t50) * t137 + ((t13 * t68 + t16 * t58) * t52 + (t13 * t67 + t16 * t55) * t49) * t135;];
tauX  = t43;
