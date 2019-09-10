% Calculate vector of inverse dynamics forces for parallel robot
% P3RPP1G1P1A0
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
%   pkin=[a2,a3,d1]';
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
% Datum: 2019-05-03 14:53
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPP1G1P1A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPP1G1P1A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPP1G1P1A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPP1G1P1A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPP1G1P1A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPP1G1P1A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RPP1G1P1A0_invdyn_para_pf_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPP1G1P1A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPP1G1P1A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPP1G1P1A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPP1G1P1A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPP1G1P1A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:52:07
% EndTime: 2019-05-03 14:52:11
% DurationCPUTime: 3.83s
% Computational Cost: add. (27102->372), mult. (36802->558), div. (2031->3), fcn. (14708->14), ass. (0->241)
t174 = -qJ(3,3) - pkin(1);
t192 = xP(3);
t165 = sin(t192);
t166 = cos(t192);
t201 = koppelP(3,2);
t204 = koppelP(3,1);
t117 = t165 * t204 + t166 * t201;
t187 = xDP(3);
t189 = xDP(1);
t100 = -t117 * t187 + t189;
t195 = (qJ(3,3) ^ 2);
t208 = (pkin(1) ^ 2);
t250 = 1 + t208;
t253 = 2 * qJ(3,3);
t144 = pkin(1) * t253 + t195 + t250;
t196 = (qJ(2,3) ^ 2);
t263 = t144 + t196;
t135 = 1 / t263;
t181 = sin(qJ(1,3));
t184 = cos(qJ(1,3));
t238 = t174 * t184;
t229 = qJ(2,3) * t238;
t103 = t144 * t181 + t229;
t239 = t174 * t181;
t147 = qJ(2,3) * t239;
t106 = t144 * t184 - t147;
t171 = legFrame(3,3);
t156 = sin(t171);
t159 = cos(t171);
t70 = t103 * t159 + t106 * t156;
t73 = -t103 * t156 + t106 * t159;
t120 = -t165 * t201 + t166 * t204;
t188 = xDP(2);
t97 = t120 * t187 + t188;
t34 = (t100 * t73 + t70 * t97) * t135;
t272 = t174 * t34;
t175 = -qJ(3,2) - pkin(1);
t202 = koppelP(2,2);
t205 = koppelP(2,1);
t118 = t165 * t205 + t166 * t202;
t101 = -t118 * t187 + t189;
t197 = qJ(3,2) ^ 2;
t254 = 2 * qJ(3,2);
t145 = pkin(1) * t254 + t197 + t250;
t198 = (qJ(2,2) ^ 2);
t262 = t145 + t198;
t136 = 1 / t262;
t182 = sin(qJ(1,2));
t185 = cos(qJ(1,2));
t236 = t175 * t185;
t230 = qJ(2,2) * t236;
t104 = t145 * t182 + t230;
t237 = t175 * t182;
t148 = qJ(2,2) * t237;
t107 = t145 * t185 - t148;
t172 = legFrame(2,3);
t157 = sin(t172);
t160 = cos(t172);
t71 = t104 * t160 + t107 * t157;
t74 = -t104 * t157 + t107 * t160;
t121 = -t165 * t202 + t166 * t205;
t98 = t121 * t187 + t188;
t35 = (t101 * t74 + t71 * t98) * t136;
t271 = t175 * t35;
t176 = -qJ(3,1) - pkin(1);
t203 = koppelP(1,2);
t206 = koppelP(1,1);
t119 = t165 * t206 + t166 * t203;
t102 = -t119 * t187 + t189;
t199 = qJ(3,1) ^ 2;
t255 = 2 * qJ(3,1);
t146 = pkin(1) * t255 + t199 + t250;
t200 = (qJ(2,1) ^ 2);
t261 = t146 + t200;
t137 = 1 / t261;
t183 = sin(qJ(1,1));
t186 = cos(qJ(1,1));
t234 = t176 * t186;
t231 = qJ(2,1) * t234;
t105 = t146 * t183 + t231;
t235 = t176 * t183;
t149 = qJ(2,1) * t235;
t108 = t146 * t186 - t149;
t173 = legFrame(1,3);
t158 = sin(t173);
t161 = cos(t173);
t72 = t105 * t161 + t108 * t158;
t75 = -t105 * t158 + t108 * t161;
t122 = -t165 * t203 + t166 * t206;
t99 = t122 * t187 + t188;
t36 = (t102 * t75 + t72 * t99) * t137;
t270 = t176 * t36;
t167 = 1 + t196;
t168 = 1 + t198;
t169 = 1 + t200;
t131 = qJ(2,1) * t186 + t235;
t134 = qJ(2,1) * t183 - t234;
t93 = t131 * t161 - t134 * t158;
t96 = t131 * t158 + t134 * t161;
t266 = (t102 * t93 + t96 * t99) * t137;
t130 = qJ(2,2) * t185 + t237;
t133 = qJ(2,2) * t182 - t236;
t92 = t130 * t160 - t133 * t157;
t95 = t130 * t157 + t133 * t160;
t265 = (t101 * t92 + t95 * t98) * t136;
t129 = qJ(2,3) * t184 + t239;
t132 = qJ(2,3) * t181 - t238;
t91 = t129 * t159 - t132 * t156;
t94 = t129 * t156 + t132 * t159;
t264 = (t100 * t91 + t94 * t97) * t135;
t177 = mrSges(3,2) + mrSges(2,3);
t190 = m(2) + m(3);
t221 = qJ(2,1) * t190 + t177;
t220 = qJ(2,2) * t190 + t177;
t219 = qJ(2,3) * t190 + t177;
t260 = 2 * m(3);
t259 = 2 * pkin(1);
t258 = 0.2e1 * t264;
t257 = 0.2e1 * t265;
t256 = 0.2e1 * t266;
t252 = 0.2e1 * t177;
t251 = -3 * t208;
t162 = t174 * m(3);
t163 = t175 * m(3);
t164 = t176 * m(3);
t249 = (-mrSges(2,2) + mrSges(3,3));
t109 = t167 * t181 - t229;
t112 = -t167 * t184 - t147;
t76 = t109 * t159 - t112 * t156;
t79 = t109 * t156 + t112 * t159;
t40 = (t100 * t76 + t79 * t97) * t135;
t110 = t168 * t182 - t230;
t113 = -t168 * t185 - t148;
t77 = t110 * t160 - t113 * t157;
t80 = t110 * t157 + t113 * t160;
t41 = (t101 * t77 + t80 * t98) * t136;
t111 = t169 * t183 - t231;
t114 = -t169 * t186 - t149;
t78 = t111 * t161 - t114 * t158;
t81 = t111 * t158 + t114 * t161;
t42 = (t102 * t78 + t81 * t99) * t137;
t247 = qJ(2,1) * t42;
t245 = qJ(2,2) * t41;
t243 = qJ(2,3) * t40;
t242 = t135 * t264;
t241 = t136 * t265;
t240 = t137 * t266;
t233 = Ifges(2,1) + Ifges(3,1) + Ifges(1,3);
t232 = m(2) * pkin(1) + t249;
t228 = mrSges(1,1) + t232;
t138 = -t232 + t162;
t139 = -t232 + t163;
t140 = -t232 + t164;
t207 = pkin(1) * t208;
t194 = mrSges(4,1);
t193 = mrSges(4,2);
t180 = xDDP(1);
t179 = xDDP(2);
t178 = xDDP(3);
t170 = t187 ^ 2;
t155 = m(3) * qJ(2,1) + mrSges(3,2);
t154 = m(3) * qJ(2,2) + mrSges(3,2);
t153 = m(3) * qJ(2,3) + mrSges(3,2);
t152 = -t164 + mrSges(3,3);
t151 = -t163 + mrSges(3,3);
t150 = -t162 + mrSges(3,3);
t143 = -mrSges(1,2) + t221;
t142 = -mrSges(1,2) + t220;
t141 = -mrSges(1,2) + t219;
t128 = g(1) * t161 + g(2) * t158;
t127 = g(1) * t160 + g(2) * t157;
t126 = g(1) * t159 + g(2) * t156;
t125 = -g(1) * t158 + g(2) * t161;
t124 = -g(1) * t157 + g(2) * t160;
t123 = -g(1) * t156 + g(2) * t159;
t116 = -t165 * t193 + t166 * t194;
t115 = t165 * t194 + t166 * t193;
t90 = (m(3) * t199) + (mrSges(3,3) * t255) + (t208 + t200) * t190 + ((m(3) * qJ(3,1) + t249) * t259) + qJ(2,1) * t252 + t233;
t89 = (m(3) * t197) + (mrSges(3,3) * t254) + (t208 + t198) * t190 + ((m(3) * qJ(3,2) + t249) * t259) + qJ(2,2) * t252 + t233;
t88 = (m(3) * t195) + (mrSges(3,3) * t253) + (t208 + t196) * t190 + ((m(3) * qJ(3,3) + t249) * t259) + qJ(2,3) * t252 + t233;
t87 = -t119 * t178 - t122 * t170 + t180;
t86 = -t118 * t178 - t121 * t170 + t180;
t85 = -t117 * t178 - t120 * t170 + t180;
t84 = -t119 * t170 + t122 * t178 + t179;
t83 = -t118 * t170 + t121 * t178 + t179;
t82 = -t117 * t170 + t120 * t178 + t179;
t63 = (t140 * t96 + t190 * t81) * t137;
t62 = (t139 * t95 + t190 * t80) * t136;
t61 = (t138 * t94 + t190 * t79) * t135;
t60 = (t140 * t93 + t190 * t78) * t137;
t59 = (t139 * t92 + t190 * t77) * t136;
t58 = (t138 * t91 + t190 * t76) * t135;
t57 = (m(3) * t72 + t155 * t96) * t137;
t56 = (m(3) * t71 + t154 * t95) * t136;
t55 = (m(3) * t70 + t153 * t94) * t135;
t54 = (m(3) * t75 + t155 * t93) * t137;
t53 = (m(3) * t74 + t154 * t92) * t136;
t52 = (m(3) * t73 + t153 * t91) * t135;
t51 = (-t119 * t93 + t122 * t96) * t137;
t50 = (-t118 * t92 + t121 * t95) * t136;
t49 = (-t117 * t91 + t120 * t94) * t135;
t45 = (-t119 * t78 + t122 * t81) * t137;
t44 = (-t118 * t77 + t121 * t80) * t136;
t43 = (-t117 * t76 + t120 * t79) * t135;
t39 = (-t119 * t75 + t122 * t72) * t137;
t38 = (-t118 * t74 + t121 * t71) * t136;
t37 = (-t117 * t73 + t120 * t70) * t135;
t33 = (t140 * t81 + t155 * t72 + t90 * t96) * t137;
t32 = (t139 * t80 + t154 * t71 + t89 * t95) * t136;
t31 = (t138 * t79 + t153 * t70 + t88 * t94) * t135;
t30 = (t140 * t78 + t155 * t75 + t90 * t93) * t137;
t29 = (t139 * t77 + t154 * t74 + t89 * t92) * t136;
t28 = (t138 * t76 + t153 * t73 + t88 * t91) * t135;
t27 = t140 * t51 + t190 * t45;
t26 = t139 * t50 + t190 * t44;
t25 = t138 * t49 + t190 * t43;
t24 = m(3) * t39 + t155 * t51;
t23 = m(3) * t38 + t154 * t50;
t22 = m(3) * t37 + t153 * t49;
t21 = t140 * t45 + t155 * t39 + t51 * t90;
t20 = t139 * t44 + t154 * t38 + t50 * t89;
t19 = t138 * t43 + t153 * t37 + t49 * t88;
t18 = 0.2e1 * (t247 - t270) * t240;
t17 = 0.2e1 * (t245 - t271) * t241;
t16 = 0.2e1 * (t243 - t272) * t242;
t15 = (-0.2e1 * t176 * t247 - t36 + (-t200 - t169) * t36 - t261 * t266 * qJ(2,1)) * t240;
t14 = (-0.2e1 * t175 * t245 - t35 + (-t198 - t168) * t35 - t262 * t265 * qJ(2,2)) * t241;
t13 = (-0.2e1 * t174 * t243 - t34 + (-t196 - t167) * t34 - t263 * t264 * qJ(2,3)) * t242;
t12 = (0.2e1 * t146 * t42 + 0.2e1 * qJ(2,1) * t270 + (-t207 + (-t199 - t200 + t251) * qJ(3,1) - qJ(3,1) + (-3 * t199 - t169) * pkin(1)) * t266) * t240;
t11 = (0.2e1 * t145 * t41 + 0.2e1 * qJ(2,2) * t271 + (-t207 + (-t197 - t198 + t251) * qJ(3,2) - qJ(3,2) + (-3 * t197 - t168) * pkin(1)) * t265) * t241;
t10 = (0.2e1 * t144 * t40 + 0.2e1 * qJ(2,3) * t272 + (-t207 + (-t195 - t196 + t251) * qJ(3,3) - qJ(3,3) + (-3 * t195 - t167) * pkin(1)) * t264) * t242;
t9 = -t140 * t18 - (t221 * t266 + t36 * t260) * t266 + (t125 * t186 - t128 * t183 - t15) * t190;
t8 = -t139 * t17 - (t220 * t265 + t35 * t260) * t265 + (t124 * t185 - t127 * t182 - t14) * t190;
t7 = -t138 * t16 - (t219 * t264 + t34 * t260) * t264 + (t123 * t184 - t126 * t181 - t13) * t190;
t6 = -t152 * t266 ^ 2 - t155 * t18 + (-t125 * t183 - t128 * t186 + t256 * t42 - t12) * m(3);
t5 = -t151 * t265 ^ 2 - t154 * t17 + (-t124 * t182 - t127 * t185 + t257 * t41 - t11) * m(3);
t4 = -t150 * t264 ^ 2 - t153 * t16 + (-t123 * t181 - t126 * t184 + t258 * t40 - t10) * m(3);
t3 = -t90 * t18 - t140 * t15 - t155 * t12 + (t36 * t152 + t221 * t42) * t256 + (-t143 * t128 - (-t164 + t228) * t125) * t186 - t183 * ((-mrSges(1,1) + t140) * t128 + t125 * t143);
t2 = -t89 * t17 - t139 * t14 - t154 * t11 + (t35 * t151 + t220 * t41) * t257 + (-t142 * t127 - (-t163 + t228) * t124) * t185 - t182 * ((-mrSges(1,1) + t139) * t127 + t124 * t142);
t1 = -t88 * t16 - t138 * t13 - t153 * t10 + (t34 * t150 + t219 * t40) * t258 + (-t141 * t126 - (-t162 + t228) * t123) * t184 - t181 * ((-mrSges(1,1) + t138) * t126 + t123 * t141);
t46 = [-t115 * t178 - t170 * t116 + (t180 - g(1)) * m(4) + ((t30 * t93 + t54 * t75 + t60 * t78) * t87 + (t30 * t96 + t54 * t72 + t60 * t81) * t84 + t93 * t3 + t78 * t9 + t75 * t6) * t137 + ((t29 * t92 + t53 * t74 + t59 * t77) * t86 + (t29 * t95 + t53 * t71 + t59 * t80) * t83 + t92 * t2 + t77 * t8 + t74 * t5) * t136 + ((t28 * t91 + t52 * t73 + t58 * t76) * t85 + (t28 * t94 + t52 * t70 + t58 * t79) * t82 + t91 * t1 + t76 * t7 + t73 * t4) * t135; -t170 * t115 + t116 * t178 + (t179 - g(2)) * m(4) + ((t33 * t93 + t57 * t75 + t63 * t78) * t87 + (t33 * t96 + t57 * t72 + t63 * t81) * t84 + t96 * t3 + t81 * t9 + t72 * t6) * t137 + ((t32 * t92 + t56 * t74 + t62 * t77) * t86 + (t32 * t95 + t56 * t71 + t62 * t80) * t83 + t95 * t2 + t80 * t8 + t71 * t5) * t136 + ((t31 * t91 + t55 * t73 + t61 * t76) * t85 + (t31 * t94 + t55 * t70 + t61 * t79) * t82 + t94 * t1 + t79 * t7 + t70 * t4) * t135; t51 * t3 + t45 * t9 + t39 * t6 + t50 * t2 + t44 * t8 + t38 * t5 + t49 * t1 + t43 * t7 + t37 * t4 - t115 * t180 + t116 * t179 + Ifges(4,3) * t178 - (-g(1) * t194 - g(2) * t193) * t165 + t166 * (g(1) * t193 - g(2) * t194) + ((t21 * t93 + t24 * t75 + t27 * t78) * t87 + (t21 * t96 + t24 * t72 + t27 * t81) * t84) * t137 + ((t20 * t92 + t23 * t74 + t26 * t77) * t86 + (t20 * t95 + t23 * t71 + t26 * t80) * t83) * t136 + ((t19 * t91 + t22 * t73 + t25 * t76) * t85 + (t19 * t94 + t22 * t70 + t25 * t79) * t82) * t135;];
tauX  = t46;
