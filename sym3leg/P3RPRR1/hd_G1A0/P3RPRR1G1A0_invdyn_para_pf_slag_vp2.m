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
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
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

function tauX = P3RPRR1G1P1A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:23:13
% EndTime: 2020-03-09 21:23:14
% DurationCPUTime: 1.10s
% Computational Cost: add. (5863->231), mult. (4576->365), div. (930->7), fcn. (3252->56), ass. (0->149)
t150 = 2 * pkin(2);
t137 = m(2) + m(3);
t121 = sin(pkin(7));
t166 = pkin(1) * t121;
t122 = cos(pkin(7));
t165 = pkin(2) * t122;
t135 = xDP(2);
t140 = 0.1e1 / pkin(3);
t154 = t135 * t140;
t123 = legFrame(3,3);
t101 = t123 + qJ(1,3);
t118 = qJ(1,3) + pkin(7);
t91 = t123 + t118;
t87 = qJ(3,3) + t91;
t81 = sin(t87);
t58 = -pkin(1) * sin(t101) - pkin(2) * sin(t91) - pkin(3) * t81;
t115 = pkin(7) + qJ(3,3);
t129 = sin(qJ(3,3));
t78 = 0.1e1 / (pkin(1) * sin(t115) + t129 * pkin(2));
t40 = t58 * t78 * t154;
t136 = xDP(1);
t153 = t136 * t140;
t84 = cos(t87);
t61 = -pkin(1) * cos(t101) - pkin(2) * cos(t91) - pkin(3) * t84;
t43 = t61 * t78 * t153;
t31 = t43 + t40;
t37 = (t135 * t81 + t136 * t84) * t78;
t22 = t37 + t31;
t164 = t22 * t31;
t124 = legFrame(2,3);
t102 = t124 + qJ(1,2);
t116 = pkin(7) + qJ(3,2);
t105 = qJ(1,2) + t116;
t88 = t124 + t105;
t82 = sin(t88);
t119 = qJ(1,2) + pkin(7);
t92 = t124 + t119;
t59 = -pkin(1) * sin(t102) - pkin(2) * sin(t92) - pkin(3) * t82;
t130 = sin(qJ(3,2));
t79 = 0.1e1 / (pkin(1) * sin(t116) + t130 * pkin(2));
t41 = t59 * t79 * t154;
t85 = cos(t88);
t62 = -pkin(1) * cos(t102) - pkin(2) * cos(t92) - pkin(3) * t85;
t44 = t62 * t79 * t153;
t32 = t44 + t41;
t38 = (t135 * t82 + t136 * t85) * t79;
t23 = t38 + t32;
t163 = t23 * t32;
t125 = legFrame(1,3);
t103 = t125 + qJ(1,1);
t120 = qJ(1,1) + pkin(7);
t93 = t125 + t120;
t89 = qJ(3,1) + t93;
t83 = sin(t89);
t60 = -pkin(1) * sin(t103) - pkin(2) * sin(t93) - pkin(3) * t83;
t117 = pkin(7) + qJ(3,1);
t131 = sin(qJ(3,1));
t80 = 0.1e1 / (pkin(1) * sin(t117) + t131 * pkin(2));
t42 = t60 * t80 * t154;
t86 = cos(t89);
t63 = -pkin(1) * cos(t103) - pkin(2) * cos(t93) - pkin(3) * t86;
t45 = t63 * t80 * t153;
t33 = t45 + t42;
t39 = (t135 * t83 + t136 * t86) * t80;
t24 = t39 + t33;
t162 = t24 * t33;
t160 = t140 * t58;
t159 = t140 * t59;
t158 = t140 * t60;
t157 = t140 * t61;
t156 = t140 * t62;
t155 = t140 * t63;
t113 = m(3) * pkin(2) + mrSges(2,1);
t152 = pkin(3) * t150;
t151 = 0.2e1 * pkin(1);
t97 = t122 * pkin(1) + pkin(2);
t104 = qJ(1,3) + t115;
t107 = sin(t123);
t110 = cos(t123);
t72 = -t107 * g(1) + t110 * g(2);
t75 = t110 * g(1) + t107 * g(2);
t149 = -cos(t104) * (mrSges(3,1) * t72 - mrSges(3,2) * t75) + sin(t104) * (mrSges(3,1) * t75 + mrSges(3,2) * t72);
t108 = sin(t124);
t111 = cos(t124);
t73 = -t108 * g(1) + t111 * g(2);
t76 = t111 * g(1) + t108 * g(2);
t148 = -cos(t105) * (mrSges(3,1) * t73 - mrSges(3,2) * t76) + sin(t105) * (mrSges(3,1) * t76 + mrSges(3,2) * t73);
t106 = qJ(1,1) + t117;
t109 = sin(t125);
t112 = cos(t125);
t74 = -t109 * g(1) + t112 * g(2);
t77 = t112 * g(1) + t109 * g(2);
t147 = -cos(t106) * (mrSges(3,1) * t74 - mrSges(3,2) * t77) + sin(t106) * (mrSges(3,1) * t77 + mrSges(3,2) * t74);
t132 = cos(qJ(3,3));
t146 = t132 * mrSges(3,1) - mrSges(3,2) * t129;
t133 = cos(qJ(3,2));
t145 = t133 * mrSges(3,1) - mrSges(3,2) * t130;
t134 = cos(qJ(3,1));
t144 = t134 * mrSges(3,1) - mrSges(3,2) * t131;
t141 = pkin(2) ^ 2;
t142 = pkin(1) ^ 2;
t143 = (m(3) * t141) + t137 * t142 + Ifges(1,3) + Ifges(2,3) + Ifges(3,3);
t139 = pkin(3) ^ 2;
t128 = xDDP(1);
t127 = xDDP(2);
t114 = t141 + t142;
t100 = cos(t117);
t99 = cos(t116);
t98 = cos(t115);
t90 = t137 * pkin(1) + mrSges(1,1);
t71 = t97 * mrSges(3,1) - mrSges(3,2) * t166;
t70 = mrSges(3,1) * t166 + t97 * mrSges(3,2);
t51 = t131 * t71 + t70 * t134;
t50 = t130 * t71 + t70 * t133;
t49 = t129 * t71 + t70 * t132;
t48 = -t70 * t131 + t71 * t134 + Ifges(3,3);
t47 = -t70 * t130 + t71 * t133 + Ifges(3,3);
t46 = -t70 * t129 + t71 * t132 + Ifges(3,3);
t36 = t144 * t150 + ((t144 + t113) * t122 - (mrSges(3,1) * t131 + t134 * mrSges(3,2) + mrSges(2,2)) * t121) * t151 + t143;
t35 = t145 * t150 + ((t145 + t113) * t122 - (mrSges(3,1) * t130 + t133 * mrSges(3,2) + mrSges(2,2)) * t121) * t151 + t143;
t34 = t146 * t150 + ((t146 + t113) * t122 - (mrSges(3,1) * t129 + t132 * mrSges(3,2) + mrSges(2,2)) * t121) * t151 + t143;
t30 = (Ifges(3,3) * t155 + t48 * t86) * t80;
t29 = (Ifges(3,3) * t156 + t47 * t85) * t79;
t28 = (Ifges(3,3) * t157 + t46 * t84) * t78;
t27 = (Ifges(3,3) * t158 + t48 * t83) * t80;
t26 = (Ifges(3,3) * t159 + t47 * t82) * t79;
t25 = (Ifges(3,3) * t160 + t46 * t81) * t78;
t21 = t45 / 0.2e1 + t42 / 0.2e1 + t39;
t20 = t44 / 0.2e1 + t41 / 0.2e1 + t38;
t19 = t43 / 0.2e1 + t40 / 0.2e1 + t37;
t18 = (t48 * t155 + t36 * t86) * t80;
t17 = (t47 * t156 + t35 * t85) * t79;
t16 = (t46 * t157 + t34 * t84) * t78;
t15 = (t48 * t158 + t36 * t83) * t80;
t14 = (t47 * t159 + t35 * t82) * t79;
t13 = (t46 * t160 + t34 * t81) * t78;
t12 = (-pkin(3) * t162 + (-t24 * pkin(3) + (-pkin(1) * t100 - pkin(2) * t134) * t39) * t39) * t80;
t11 = (-pkin(3) * t163 + (-pkin(3) * t23 + (-pkin(1) * t99 - pkin(2) * t133) * t38) * t38) * t79;
t10 = (-pkin(3) * t164 + (-pkin(3) * t22 + (-pkin(1) * t98 - pkin(2) * t132) * t37) * t37) * t78;
t9 = (t21 * t134 * t152 + t39 * t114 + t24 * t139 + (pkin(3) * t100 * t21 + t39 * t165) * t151) * t140 * t80 * t39 + (-t131 * t166 + t134 * t97 + pkin(3)) / (t131 * t97 + t134 * t166) * t162;
t8 = (t20 * t133 * t152 + t38 * t114 + t23 * t139 + (pkin(3) * t20 * t99 + t38 * t165) * t151) * t140 * t79 * t38 + (-t130 * t166 + t97 * t133 + pkin(3)) / (t97 * t130 + t133 * t166) * t163;
t7 = (t19 * t132 * t152 + t37 * t114 + t22 * t139 + (pkin(3) * t19 * t98 + t37 * t165) * t151) * t140 * t78 * t37 + (-t129 * t166 + t97 * t132 + pkin(3)) / (t97 * t129 + t132 * t166) * t164;
t6 = t51 * t39 ^ 2 - Ifges(3,3) * t9 - t48 * t12 + t147;
t5 = t50 * t38 ^ 2 - Ifges(3,3) * t8 - t47 * t11 + t148;
t4 = t49 * t37 ^ 2 - Ifges(3,3) * t7 - t46 * t10 + t149;
t3 = -t36 * t12 - t48 * t9 - 0.2e1 * t51 * t21 * t33 + (t77 * mrSges(2,2) - t74 * t113) * cos(t120) + (t74 * mrSges(2,2) + t113 * t77) * sin(t120) + (mrSges(1,2) * t77 - t74 * t90) * cos(qJ(1,1)) + (t74 * mrSges(1,2) + t90 * t77) * sin(qJ(1,1)) + t147;
t2 = -t35 * t11 - t47 * t8 - 0.2e1 * t50 * t20 * t32 + (t76 * mrSges(2,2) - t73 * t113) * cos(t119) + (t73 * mrSges(2,2) + t113 * t76) * sin(t119) + (mrSges(1,2) * t76 - t73 * t90) * cos(qJ(1,2)) + (t73 * mrSges(1,2) + t90 * t76) * sin(qJ(1,2)) + t148;
t1 = -t34 * t10 - t46 * t7 - 0.2e1 * t49 * t19 * t31 + (t75 * mrSges(2,2) - t72 * t113) * cos(t118) + (t72 * mrSges(2,2) + t113 * t75) * sin(t118) + (mrSges(1,2) * t75 - t72 * t90) * cos(qJ(1,3)) + (t72 * mrSges(1,2) + t90 * t75) * sin(qJ(1,3)) + t149;
t52 = [(-g(1) + t128) * m(4) + ((t30 * t155 + t18 * t86) * t128 + (t30 * t158 + t18 * t83) * t127 + t86 * t3 + t6 * t155) * t80 + ((t29 * t156 + t17 * t85) * t128 + (t29 * t159 + t17 * t82) * t127 + t85 * t2 + t5 * t156) * t79 + ((t28 * t157 + t16 * t84) * t128 + (t16 * t81 + t28 * t160) * t127 + t84 * t1 + t4 * t157) * t78; (-g(2) + t127) * m(4) + ((t15 * t86 + t27 * t155) * t128 + (t15 * t83 + t27 * t158) * t127 + t83 * t3 + t6 * t158) * t80 + ((t14 * t85 + t26 * t156) * t128 + (t14 * t82 + t26 * t159) * t127 + t82 * t2 + t5 * t159) * t79 + ((t13 * t84 + t25 * t157) * t128 + (t13 * t81 + t25 * t160) * t127 + t81 * t1 + t4 * t160) * t78; (-g(3) + xDDP(3)) * (m(4) + 0.3e1 * t137);];
tauX  = t52;
