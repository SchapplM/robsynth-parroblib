% Calculate vector of inverse dynamics forces for parallel robot
% P3RRPRR1V1G2A0
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4]';
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
% Datum: 2022-11-04 17:08
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRPRR1V1G2A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-04 17:08:01
% EndTime: 2022-11-04 17:08:04
% DurationCPUTime: 2.16s
% Computational Cost: add. (7134->330), mult. (7797->607), div. (2637->7), fcn. (7074->18), ass. (0->206)
t222 = (pkin(1) * m(3));
t221 = (mrSges(2,3) - mrSges(1,2));
t127 = legFrame(3,2);
t104 = sin(t127);
t107 = cos(t127);
t146 = xDP(2);
t147 = xDP(1);
t139 = cos(qJ(2,3));
t116 = 0.1e1 / t139;
t148 = pkin(2) + pkin(1);
t122 = 0.1e1 / t148;
t189 = t116 * t122;
t64 = (t104 * t147 + t107 * t146) * t189;
t61 = t64 ^ 2;
t128 = legFrame(2,2);
t105 = sin(t128);
t108 = cos(t128);
t141 = cos(qJ(2,2));
t118 = 0.1e1 / t141;
t188 = t118 * t122;
t65 = (t105 * t147 + t108 * t146) * t188;
t62 = t65 ^ 2;
t129 = legFrame(1,2);
t106 = sin(t129);
t109 = cos(t129);
t143 = cos(qJ(2,1));
t120 = 0.1e1 / t143;
t187 = t120 * t122;
t66 = (t106 * t147 + t109 * t146) * t187;
t63 = t66 ^ 2;
t220 = 2 * mrSges(3,3);
t123 = pkin(3) + qJ(3,3);
t111 = 0.1e1 / t123;
t140 = cos(qJ(1,3));
t145 = xDP(3);
t133 = sin(qJ(2,3));
t134 = sin(qJ(1,3));
t184 = t134 * t139;
t73 = -t104 * t184 + t107 * t133;
t76 = t104 * t133 + t107 * t184;
t40 = (t140 * t145 + (t146 * t73 + t147 * t76) * t116) * t111;
t219 = 0.2e1 * t40;
t124 = pkin(3) + qJ(3,2);
t112 = 0.1e1 / t124;
t142 = cos(qJ(1,2));
t135 = sin(qJ(2,2));
t136 = sin(qJ(1,2));
t181 = t136 * t141;
t74 = -t105 * t181 + t108 * t135;
t77 = t105 * t135 + t108 * t181;
t41 = (t142 * t145 + (t146 * t74 + t147 * t77) * t118) * t112;
t218 = 0.2e1 * t41;
t125 = pkin(3) + qJ(3,1);
t113 = 0.1e1 / t125;
t144 = cos(qJ(1,1));
t137 = sin(qJ(2,1));
t138 = sin(qJ(1,1));
t178 = t138 * t143;
t75 = -t106 * t178 + t109 * t137;
t78 = t106 * t137 + t109 * t178;
t42 = (t144 * t145 + (t146 * t75 + t147 * t78) * t120) * t113;
t217 = 0.2e1 * t42;
t162 = (2 * mrSges(3,1) + t222) * pkin(1);
t211 = (Ifges(2,1) + Ifges(3,1));
t91 = Ifges(2,2) + Ifges(3,2) + t162 - t211;
t216 = 2 * t91;
t99 = mrSges(3,2) * pkin(1) - Ifges(2,4) - Ifges(3,4);
t215 = -2 * t99;
t121 = t148 ^ 2;
t214 = (m(3) * qJ(3,1));
t213 = (m(3) * qJ(3,2));
t212 = (m(3) * qJ(3,3));
t210 = (-Ifges(2,5) - Ifges(3,5));
t209 = (-Ifges(2,6) - Ifges(3,6));
t110 = mrSges(3,1) + t222;
t163 = mrSges(3,2) * qJ(3,3) + t209;
t166 = -pkin(1) * mrSges(3,3) - t210;
t67 = (-qJ(3,3) * t110 + t166) * t133 - t163 * t139;
t208 = t111 * t67;
t164 = mrSges(3,2) * qJ(3,2) + t209;
t68 = (-qJ(3,2) * t110 + t166) * t135 - t164 * t141;
t207 = t112 * t68;
t165 = mrSges(3,2) * qJ(3,1) + t209;
t69 = (-qJ(3,1) * t110 + t166) * t137 - t165 * t143;
t206 = t113 * t69;
t115 = t139 ^ 2;
t205 = t115 * t99;
t204 = t116 * t73;
t203 = t116 * t76;
t82 = mrSges(3,2) * t133 - t110 * t139;
t202 = t116 * t82;
t117 = t141 ^ 2;
t201 = t117 * t99;
t200 = t118 * t74;
t199 = t118 * t77;
t83 = mrSges(3,2) * t135 - t110 * t141;
t198 = t118 * t83;
t119 = t143 ^ 2;
t197 = t119 * t99;
t196 = t120 * t75;
t195 = t120 * t78;
t84 = mrSges(3,2) * t137 - t110 * t143;
t194 = t120 * t84;
t95 = Ifges(2,3) + Ifges(3,3) + t162;
t193 = t122 * t95;
t192 = t133 * t61;
t191 = t135 * t62;
t190 = t137 * t63;
t186 = t133 * t139;
t185 = t133 * t148;
t183 = t135 * t141;
t182 = t135 * t148;
t180 = t137 * t143;
t179 = t137 * t148;
t177 = t139 * t148;
t176 = t141 * t148;
t175 = t143 * t148;
t174 = 0.2e1 * t148;
t173 = Ifges(1,3) + t211;
t152 = -t123 * t140 + t134 * t177;
t55 = -t104 * t152 + t107 * t185;
t58 = t104 * t185 + t107 * t152;
t79 = t123 * t134 + t140 * t177;
t22 = (t145 * t79 + t146 * t55 + t147 * t58) * t111;
t151 = -t124 * t142 + t136 * t176;
t56 = -t105 * t151 + t108 * t182;
t59 = t105 * t182 + t108 * t151;
t80 = t124 * t136 + t142 * t176;
t23 = (t145 * t80 + t146 * t56 + t147 * t59) * t112;
t150 = -t125 * t144 + t138 * t175;
t57 = -t106 * t150 + t109 * t179;
t60 = t106 * t179 + t109 * t150;
t81 = t125 * t138 + t144 * t175;
t24 = (t145 * t81 + t146 * t57 + t147 * t60) * t113;
t172 = t67 * t189;
t171 = t68 * t188;
t170 = t69 * t187;
t169 = t122 * t140 * t67;
t168 = t122 * t142 * t68;
t167 = t122 * t144 * t69;
t103 = mrSges(3,3) + t214;
t102 = mrSges(3,3) + t213;
t101 = mrSges(3,3) + t212;
t88 = g(1) * t107 - g(2) * t104;
t161 = g(3) * t140 + t134 * t88;
t89 = g(1) * t108 - g(2) * t105;
t160 = g(3) * t142 + t136 * t89;
t90 = g(1) * t109 - g(2) * t106;
t159 = g(3) * t144 + t138 * t90;
t158 = mrSges(3,2) * t139 + t110 * t133;
t157 = mrSges(3,2) * t141 + t110 * t135;
t156 = mrSges(3,2) * t143 + t110 * t137;
t100 = mrSges(2,1) + t110;
t126 = mrSges(2,2) + mrSges(3,2);
t155 = -t100 * t139 + t126 * t133 - mrSges(1,1);
t154 = -t100 * t141 + t126 * t135 - mrSges(1,1);
t153 = -t100 * t143 + t126 * t137 - mrSges(1,1);
t132 = xDDP(1);
t131 = xDDP(2);
t130 = xDDP(3);
t98 = t103 + t221;
t97 = t102 + t221;
t96 = t101 + t221;
t87 = g(1) * t106 + g(2) * t109;
t86 = g(1) * t105 + g(2) * t108;
t85 = g(1) * t104 + g(2) * t107;
t54 = t180 * t215 + t119 * t91 + ((t220 + t214) * qJ(3,1)) + t173;
t53 = t183 * t215 + t117 * t91 + ((t220 + t213) * qJ(3,2)) + t173;
t52 = t186 * t215 + t115 * t91 + ((t220 + t212) * qJ(3,3)) + t173;
t39 = (t106 * t193 + t206 * t78) * t120;
t38 = (t105 * t193 + t207 * t77) * t118;
t37 = (t104 * t193 + t208 * t76) * t116;
t36 = (t109 * t193 + t206 * t75) * t120;
t35 = (t108 * t193 + t207 * t74) * t118;
t34 = (t107 * t193 + t208 * t73) * t116;
t33 = (t144 * t54 + t81 * t84) * t113;
t32 = (t142 * t53 + t80 * t83) * t112;
t31 = (t140 * t52 + t79 * t82) * t111;
t30 = (m(3) * t60 + t194 * t78) * t113;
t29 = (m(3) * t59 + t198 * t77) * t112;
t28 = (m(3) * t58 + t202 * t76) * t111;
t27 = (m(3) * t57 + t194 * t75) * t113;
t26 = (m(3) * t56 + t198 * t74) * t112;
t25 = (m(3) * t55 + t202 * t73) * t111;
t21 = t106 * t170 + (t195 * t54 + t60 * t84) * t113;
t20 = t105 * t171 + (t199 * t53 + t59 * t83) * t112;
t19 = t104 * t172 + (t203 * t52 + t58 * t82) * t111;
t18 = t109 * t170 + (t196 * t54 + t57 * t84) * t113;
t17 = t108 * t171 + (t200 * t53 + t56 * t83) * t112;
t16 = t107 * t172 + (t204 * t52 + t55 * t82) * t111;
t15 = (-t148 * t63 * t120 + (-t175 * t42 + 0.2e1 * t24) * t42) * t113;
t14 = (-t148 * t62 * t118 + (-t176 * t41 + 0.2e1 * t23) * t41) * t112;
t13 = (-t148 * t61 * t116 + (-t177 * t40 + 0.2e1 * t22) * t40) * t111;
t12 = (-t121 * t63 + ((-t119 * t121 - t125 ^ 2) * t42 + (t125 * t137 * t66 + t143 * t24) * t174) * t42) * t113;
t11 = (-t121 * t62 + ((-t117 * t121 - t124 ^ 2) * t41 + (t124 * t135 * t65 + t141 * t23) * t174) * t41) * t112;
t10 = (-t121 * t61 + ((-t115 * t121 - t123 ^ 2) * t40 + (t123 * t133 * t64 + t139 * t22) * t174) * t40) * t111;
t9 = -t69 * t15 + t95 * t120 * t190 + (-0.2e1 * t156 * t24 + (t180 * t91 + 0.2e1 * t197 - t99) * t42) * t42 + (t100 * t159 + t87 * t126) * t137 - (t87 * t100 - t126 * t159) * t143;
t8 = -t68 * t14 + t95 * t118 * t191 + (-0.2e1 * t157 * t23 + (t183 * t91 + 0.2e1 * t201 - t99) * t41) * t41 + (t100 * t160 + t86 * t126) * t135 - (t86 * t100 - t126 * t160) * t141;
t7 = -t67 * t13 + t95 * t116 * t192 + (-0.2e1 * t158 * t22 + (t186 * t91 + 0.2e1 * t205 - t99) * t40) * t40 + (t100 * t161 + t85 * t126) * t133 - (t85 * t100 - t126 * t161) * t139;
t6 = -t84 * t15 + (-t42 * t103 / 0.2e1 + t156 * t66) * t217 + (-g(3) * t138 + t144 * t90 - t12) * m(3);
t5 = -t83 * t14 + (-t41 * t102 / 0.2e1 + t157 * t65) * t218 + (-g(3) * t136 + t142 * t89 - t11) * m(3);
t4 = -t82 * t13 + (-t40 * t101 / 0.2e1 + t158 * t64) * t219 + (-g(3) * t134 + t140 * t88 - t10) * m(3);
t3 = -t54 * t15 - t84 * t12 - 0.4e1 * t42 * t66 * t197 - (t42 * t137 * t216 + t66 * (qJ(3,1) * mrSges(3,1) + pkin(1) * t103 + t210)) * t66 * t143 + (t103 * t24 + t66 * t99) * t217 + (t120 * t69 + t165) * t190 + (-t138 * t98 + t144 * t153) * t90 + (-t138 * t153 - t98 * t144) * g(3);
t2 = -t53 * t14 - t83 * t11 - 0.4e1 * t41 * t65 * t201 - (t41 * t135 * t216 + t65 * (qJ(3,2) * mrSges(3,1) + pkin(1) * t102 + t210)) * t65 * t141 + (t102 * t23 + t65 * t99) * t218 + (t118 * t68 + t164) * t191 + (-t136 * t97 + t142 * t154) * t89 + (-t136 * t154 - t97 * t142) * g(3);
t1 = -t52 * t13 - t82 * t10 - 0.4e1 * t40 * t64 * t205 - (t40 * t133 * t216 + t64 * (qJ(3,3) * mrSges(3,1) + pkin(1) * t101 + t210)) * t64 * t139 + (t101 * t22 + t64 * t99) * t219 + (t116 * t67 + t163) * t192 + (-t134 * t96 + t140 * t155) * t88 + (-t134 * t155 - t96 * t140) * g(3);
t43 = [(-g(1) + t132) * m(4) + ((t195 * t21 + t30 * t60) * t132 + (t196 * t21 + t30 * t57) * t131 + (t144 * t21 + t30 * t81) * t130 + t3 * t195 + t60 * t6) * t113 + ((t109 * t131 * t39 + (t132 * t39 + t9) * t106) * t120 + (t108 * t131 * t38 + (t132 * t38 + t8) * t105) * t118 + (t107 * t131 * t37 + (t132 * t37 + t7) * t104) * t116) * t122 + ((t199 * t20 + t29 * t59) * t132 + (t20 * t200 + t29 * t56) * t131 + (t142 * t20 + t29 * t80) * t130 + t2 * t199 + t59 * t5) * t112 + ((t19 * t203 + t28 * t58) * t132 + (t19 * t204 + t28 * t55) * t131 + (t140 * t19 + t28 * t79) * t130 + t1 * t203 + t58 * t4) * t111; (-g(2) + t131) * m(4) + ((t18 * t195 + t27 * t60) * t132 + (t18 * t196 + t27 * t57) * t131 + (t144 * t18 + t27 * t81) * t130 + t3 * t196 + t57 * t6) * t113 + ((t106 * t132 * t36 + (t131 * t36 + t9) * t109) * t120 + (t105 * t132 * t35 + (t131 * t35 + t8) * t108) * t118 + (t104 * t132 * t34 + (t131 * t34 + t7) * t107) * t116) * t122 + ((t17 * t199 + t26 * t59) * t132 + (t17 * t200 + t26 * t56) * t131 + (t142 * t17 + t26 * t80) * t130 + t2 * t200 + t56 * t5) * t112 + ((t16 * t203 + t25 * t58) * t132 + (t16 * t204 + t25 * t55) * t131 + (t140 * t16 + t25 * t79) * t130 + t1 * t204 + t55 * t4) * t111; (-g(3) + t130) * m(4) + (t81 * t6 + (t130 * t33 + t3) * t144 + (t81 * t130 + t57 * t131 + t60 * t132) * (m(3) * t81 + t144 * t84) * t113 + ((t106 * t167 + t33 * t78) * t132 + (t109 * t167 + t33 * t75) * t131) * t120) * t113 + (t80 * t5 + (t130 * t32 + t2) * t142 + (t80 * t130 + t56 * t131 + t59 * t132) * (m(3) * t80 + t142 * t83) * t112 + ((t105 * t168 + t32 * t77) * t132 + (t108 * t168 + t32 * t74) * t131) * t118) * t112 + (t79 * t4 + (t130 * t31 + t1) * t140 + (t79 * t130 + t55 * t131 + t58 * t132) * (m(3) * t79 + t140 * t82) * t111 + ((t104 * t169 + t31 * t76) * t132 + (t107 * t169 + t31 * t73) * t131) * t116) * t111;];
tauX  = t43;
