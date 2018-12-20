% Calculate vector of inverse dynamics forces for parallel robot
% P3PRP2A0
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
%   pkin=[a2,a3,d2]';
% m [4x1]
%   mass of all robot links (including platform)
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
% Datum: 2018-12-20 17:39
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauX = P3PRP2A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRP2A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRP2A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRP2A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRP2A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRP2A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PRP2A0_invdyn_para_pf_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRP2A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRP2A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRP2A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRP2A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRP2A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:38:57
% EndTime: 2018-12-20 17:39:00
% DurationCPUTime: 3.48s
% Computational Cost: add. (26085->404), mult. (49288->613), div. (2214->3), fcn. (23453->14), ass. (0->248)
t191 = (qJ(3,1) ^ 2);
t199 = (pkin(2) ^ 2);
t257 = -t199 - 1;
t151 = -t191 - t257;
t182 = cos(qJ(2,1));
t168 = t182 ^ 2;
t179 = sin(qJ(2,1));
t227 = t179 * t182;
t214 = pkin(2) * t227;
t221 = t191 + t199;
t261 = 2 * qJ(3,1);
t105 = 0.1e1 / (t151 * t168 + t214 * t261 - t221 - 0.1e1);
t186 = xP(3);
t162 = sin(t186);
t163 = cos(t186);
t194 = koppelP(1,2);
t197 = koppelP(1,1);
t131 = -t162 * t194 + t163 * t197;
t183 = xDP(3);
t184 = xDP(2);
t108 = t131 * t183 + t184;
t173 = legFrame(1,3);
t160 = cos(t173);
t157 = sin(t173);
t240 = qJ(3,1) * t160;
t219 = pkin(2) * t240;
t206 = t157 * t191 - t219;
t241 = qJ(3,1) * t157;
t141 = pkin(2) * t241;
t224 = t191 * t160 + t141;
t90 = t206 * t182 - t179 * (t160 + t224);
t245 = t90 * t108;
t128 = t162 * t197 + t163 * t194;
t185 = xDP(1);
t111 = -t128 * t183 + t185;
t87 = t224 * t182 - t179 * (-t157 - t206);
t248 = t87 * t111;
t271 = (t245 + t248) * t105;
t190 = (qJ(3,2) ^ 2);
t150 = -t190 - t257;
t181 = cos(qJ(2,2));
t167 = t181 ^ 2;
t178 = sin(qJ(2,2));
t228 = t178 * t181;
t215 = pkin(2) * t228;
t222 = t190 + t199;
t260 = 2 * qJ(3,2);
t104 = 0.1e1 / (t150 * t167 + t215 * t260 - t222 - 0.1e1);
t193 = koppelP(2,2);
t196 = koppelP(2,1);
t130 = -t162 * t193 + t163 * t196;
t107 = t130 * t183 + t184;
t172 = legFrame(2,3);
t159 = cos(t172);
t156 = sin(t172);
t236 = qJ(3,2) * t159;
t218 = pkin(2) * t236;
t205 = t156 * t190 - t218;
t237 = qJ(3,2) * t156;
t140 = pkin(2) * t237;
t225 = t190 * t159 + t140;
t89 = t205 * t181 - t178 * (t159 + t225);
t246 = t89 * t107;
t127 = t162 * t196 + t163 * t193;
t110 = -t127 * t183 + t185;
t86 = t225 * t181 - t178 * (-t156 - t205);
t249 = t86 * t110;
t270 = (t246 + t249) * t104;
t189 = (qJ(3,3) ^ 2);
t149 = -t189 - t257;
t180 = cos(qJ(2,3));
t166 = t180 ^ 2;
t177 = sin(qJ(2,3));
t229 = t177 * t180;
t216 = pkin(2) * t229;
t223 = t189 + t199;
t259 = 2 * qJ(3,3);
t103 = 0.1e1 / (t149 * t166 + t216 * t259 - t223 - 0.1e1);
t192 = koppelP(3,2);
t195 = koppelP(3,1);
t129 = -t162 * t192 + t163 * t195;
t106 = t129 * t183 + t184;
t171 = legFrame(3,3);
t158 = cos(t171);
t155 = sin(t171);
t232 = qJ(3,3) * t158;
t217 = pkin(2) * t232;
t204 = t155 * t189 - t217;
t233 = qJ(3,3) * t155;
t139 = pkin(2) * t233;
t226 = t189 * t158 + t139;
t88 = t204 * t180 - t177 * (t158 + t226);
t247 = t88 * t106;
t126 = t162 * t195 + t163 * t192;
t109 = -t126 * t183 + t185;
t85 = t226 * t180 - t177 * (-t155 - t204);
t250 = t85 * t109;
t269 = (t247 + t250) * t103;
t268 = 0.2e1 * pkin(2);
t230 = qJ(3,3) * t180;
t97 = 0.2e1 * t155 * t230 - t177 * (pkin(2) * t155 + t232);
t98 = -0.2e1 * t158 * t230 + t177 * (pkin(2) * t158 - t233);
t67 = (t106 * t98 + t109 * t97) * t103;
t267 = t67 ^ 2;
t234 = qJ(3,2) * t181;
t100 = -0.2e1 * t159 * t234 + t178 * (pkin(2) * t159 - t237);
t99 = 0.2e1 * t156 * t234 - t178 * (pkin(2) * t156 + t236);
t68 = (t100 * t107 + t110 * t99) * t104;
t266 = t68 ^ 2;
t238 = qJ(3,1) * t182;
t101 = 0.2e1 * t157 * t238 - t179 * (pkin(2) * t157 + t240);
t102 = -0.2e1 * t160 * t238 + t179 * (pkin(2) * t160 - t241);
t69 = (t101 * t111 + t102 * t108) * t105;
t265 = t69 ^ 2;
t264 = 0.2e1 * t67;
t263 = 0.2e1 * t68;
t262 = 0.2e1 * t69;
t258 = m(3) * pkin(2);
t256 = m(3) * t180;
t255 = m(3) * t181;
t254 = m(3) * t182;
t253 = t103 * t67;
t252 = t104 * t68;
t251 = t105 * t69;
t198 = pkin(2) * t199;
t58 = t189 * t67;
t212 = -t199 * t67 + t268 * t269 - t58;
t220 = -t199 + t257;
t10 = (t212 * t230 + ((t198 + (1 + t189) * pkin(2)) * t67 - t269 + t220 * t269) * t177) * t253;
t135 = g(1) * t158 + g(2) * t155;
t244 = t10 + t135;
t59 = t190 * t68;
t211 = -t199 * t68 + t268 * t270 - t59;
t11 = (t211 * t234 + ((t198 + (1 + t190) * pkin(2)) * t68 - t270 + t220 * t270) * t178) * t252;
t136 = g(1) * t159 + g(2) * t156;
t243 = t11 + t136;
t60 = t191 * t69;
t210 = -t199 * t69 + t268 * t271 - t60;
t12 = (t210 * t238 + ((t198 + (1 + t191) * pkin(2)) * t69 - t271 + t220 * t271) * t179) * t251;
t137 = g(1) * t160 + g(2) * t157;
t242 = t12 + t137;
t239 = qJ(3,1) * t168;
t235 = qJ(3,2) * t167;
t231 = qJ(3,3) * t166;
t161 = mrSges(3,1) + t258;
t154 = m(3) * qJ(3,1) + mrSges(3,3);
t153 = m(3) * qJ(3,2) + mrSges(3,3);
t152 = m(3) * qJ(3,3) + mrSges(3,3);
t213 = mrSges(3,1) * t268 + Ifges(3,2) + Ifges(2,3);
t188 = mrSges(4,1);
t187 = mrSges(4,2);
t176 = xDDP(1);
t175 = xDDP(2);
t174 = xDDP(3);
t169 = t183 ^ 2;
t164 = m(1) + m(2) + m(3);
t148 = mrSges(2,1) + t161;
t147 = -mrSges(2,2) + t154;
t146 = -mrSges(2,2) + t153;
t145 = -mrSges(2,2) + t152;
t138 = -t258 / 0.2e1 - mrSges(2,1) / 0.2e1 - mrSges(3,1) / 0.2e1;
t134 = -g(1) * t157 + g(2) * t160;
t133 = -g(1) * t156 + g(2) * t159;
t132 = -g(1) * t155 + g(2) * t158;
t125 = m(3) * t221 + mrSges(3,3) * t261 + t213;
t124 = m(3) * t222 + mrSges(3,3) * t260 + t213;
t123 = m(3) * t223 + mrSges(3,3) * t259 + t213;
t122 = -t187 * t162 + t163 * t188;
t121 = t188 * t162 + t163 * t187;
t120 = t151 * t157 + 0.2e1 * t219;
t119 = t150 * t156 + 0.2e1 * t218;
t118 = t149 * t155 + 0.2e1 * t217;
t117 = t151 * t160 - 0.2e1 * t141;
t116 = t150 * t159 - 0.2e1 * t140;
t115 = t149 * t158 - 0.2e1 * t139;
t114 = t147 * t179 + t148 * t182;
t113 = t146 * t178 + t148 * t181;
t112 = t145 * t177 + t148 * t180;
t96 = -t128 * t174 - t131 * t169 + t176;
t95 = -t127 * t174 - t130 * t169 + t176;
t94 = -t126 * t174 - t129 * t169 + t176;
t93 = -t128 * t169 + t131 * t174 + t175;
t92 = -t127 * t169 + t130 * t174 + t175;
t91 = -t126 * t169 + t129 * t174 + t175;
t84 = -t117 * t227 + t120 * t168 - t157 * t199 - t157 - t219;
t83 = -t116 * t228 + t119 * t167 - t156 * t199 - t156 - t218;
t82 = -t115 * t229 + t118 * t166 - t155 * t199 - t155 - t217;
t81 = t117 * t168 + t120 * t227 - t160 * t199 + t141 - t160;
t80 = t116 * t167 + t119 * t228 - t159 * t199 + t140 - t159;
t79 = t115 * t166 + t118 * t229 - t158 * t199 + t139 - t158;
t72 = (-t101 * t128 + t102 * t131) * t105;
t71 = (t100 * t130 - t127 * t99) * t104;
t70 = (-t126 * t97 + t129 * t98) * t103;
t66 = pkin(2) * t68;
t65 = pkin(2) * t67;
t61 = t69 * pkin(2);
t57 = (-t128 * t87 + t131 * t90) * t105;
t56 = (-t127 * t86 + t130 * t89) * t104;
t55 = (-t126 * t85 + t129 * t88) * t103;
t51 = (-t128 * t81 + t131 * t84) * t105;
t50 = (-t127 * t80 + t130 * t83) * t104;
t49 = (-t126 * t79 + t129 * t82) * t103;
t48 = (-t102 * t161 + (-t182 * t84 + t90) * m(3)) * t105;
t47 = (-t100 * t161 + (-t181 * t83 + t89) * m(3)) * t104;
t46 = (-t161 * t98 + (-t180 * t82 + t88) * m(3)) * t103;
t45 = (-t101 * t161 + (-t182 * t81 + t87) * m(3)) * t105;
t44 = (-t161 * t99 + (-t181 * t80 + t86) * m(3)) * t104;
t43 = (-t161 * t97 + (-t180 * t79 + t85) * m(3)) * t103;
t42 = (t102 * t114 + t164 * t84 - t254 * t90) * t105;
t41 = (t100 * t113 + t164 * t83 - t255 * t89) * t104;
t40 = (t112 * t98 + t164 * t82 - t256 * t88) * t103;
t39 = (t101 * t114 + t164 * t81 - t254 * t87) * t105;
t38 = (t113 * t99 + t164 * t80 - t255 * t86) * t104;
t37 = (t112 * t97 + t164 * t79 - t256 * t85) * t103;
t36 = (t102 * t125 + t114 * t84 - t161 * t90) * t105;
t35 = (t100 * t124 + t113 * t83 - t161 * t89) * t104;
t34 = (t112 * t82 + t123 * t98 - t161 * t88) * t103;
t33 = (t101 * t125 + t114 * t81 - t161 * t87) * t105;
t32 = (t113 * t80 + t124 * t99 - t161 * t86) * t104;
t31 = (t112 * t79 + t123 * t97 - t161 * t85) * t103;
t30 = t66 - t270;
t29 = t65 - t269;
t28 = t61 - t271;
t27 = -t161 * t72 + (-t182 * t51 + t57) * m(3);
t26 = -t161 * t71 + (-t181 * t50 + t56) * m(3);
t25 = -t161 * t70 + (-t180 * t49 + t55) * m(3);
t24 = t114 * t72 + t164 * t51 - t254 * t57;
t23 = t113 * t71 + t164 * t50 - t255 * t56;
t22 = t112 * t70 + t164 * t49 - t256 * t55;
t21 = t114 * t51 + t125 * t72 - t161 * t57;
t20 = t113 * t50 + t124 * t71 - t161 * t56;
t19 = t112 * t49 + t123 * t70 - t161 * t55;
t18 = ((0.2e1 * (t66 + (-t249 / 0.2e1 - t246 / 0.2e1) * t104) * t235 - (pkin(2) * t30 - t59) * t228 - qJ(3,2) * t270) * t104 + (-qJ(3,2) + t215 - t235) * t104 * t270) * t68;
t17 = ((0.2e1 * (t65 + (-t250 / 0.2e1 - t247 / 0.2e1) * t103) * t231 - (pkin(2) * t29 - t58) * t229 - qJ(3,3) * t269) * t103 + (-qJ(3,3) + t216 - t231) * t103 * t269) * t67;
t16 = ((0.2e1 * (t61 + (-t248 / 0.2e1 - t245 / 0.2e1) * t105) * t239 - (pkin(2) * t28 - t60) * t227 - qJ(3,1) * t271) * t105 + (-qJ(3,1) + t214 - t239) * t105 * t271) * t69;
t15 = ((t28 - t271) * t227 + (-t168 * t69 - t210 + t69) * qJ(3,1)) * t251;
t14 = ((t30 - t270) * t228 + (-t167 * t68 - t211 + t68) * qJ(3,2)) * t252;
t13 = ((t29 - t269) * t229 + (-t166 * t67 - t212 + t67) * qJ(3,3)) * t253;
t9 = -t265 * t154 + t161 * t16 + (-t179 * t134 + t182 * t242 - t15) * m(3);
t8 = -t266 * t153 + t161 * t18 + (-t178 * t133 + t181 * t243 - t14) * m(3);
t7 = -t267 * t152 + t161 * t17 + (-t177 * t132 + t180 * t244 - t13) * m(3);
t6 = -t114 * t12 - t125 * t16 + t161 * t15 + t271 * t154 * t262 + (-t134 * t147 - t137 * t148) * t182 + (t134 * t148 - t137 * t147) * t179;
t5 = -t113 * t11 - t124 * t18 + t161 * t14 + t270 * t153 * t263 + (-t133 * t146 - t136 * t148) * t181 + (t133 * t148 - t136 * t146) * t178;
t4 = -t112 * t10 - t123 * t17 + t161 * t13 + t269 * t152 * t264 + (-t132 * t145 - t135 * t148) * t180 + (t132 * t148 - t135 * t145) * t177;
t3 = -t114 * t16 + (m(3) * t271 + t138 * t69) * t179 * t262 + (m(3) * t15 + t147 * t265) * t182 - t242 * t164;
t2 = -t113 * t18 + (m(3) * t270 + t138 * t68) * t178 * t263 + (m(3) * t14 + t146 * t266) * t181 - t243 * t164;
t1 = -t112 * t17 + (m(3) * t269 + t138 * t67) * t177 * t264 + (m(3) * t13 + t145 * t267) * t180 - t244 * t164;
t52 = [-t121 * t174 - t169 * t122 + (t176 - g(1)) * m(4) + ((t101 * t33 + t39 * t81 + t45 * t87) * t96 + (t102 * t33 + t39 * t84 + t45 * t90) * t93 + t81 * t3 + t101 * t6 + t87 * t9) * t105 + ((t32 * t99 + t38 * t80 + t44 * t86) * t95 + (t100 * t32 + t38 * t83 + t44 * t89) * t92 + t80 * t2 + t99 * t5 + t86 * t8) * t104 + ((t31 * t97 + t37 * t79 + t43 * t85) * t94 + (t31 * t98 + t37 * t82 + t43 * t88) * t91 + t79 * t1 + t97 * t4 + t85 * t7) * t103; -t169 * t121 + t122 * t174 + (t175 - g(2)) * m(4) + ((t101 * t36 + t42 * t81 + t48 * t87) * t96 + (t102 * t36 + t42 * t84 + t48 * t90) * t93 + t84 * t3 + t102 * t6 + t90 * t9) * t105 + ((t35 * t99 + t41 * t80 + t47 * t86) * t95 + (t100 * t35 + t41 * t83 + t47 * t89) * t92 + t83 * t2 + t100 * t5 + t89 * t8) * t104 + ((t34 * t97 + t40 * t79 + t46 * t85) * t94 + (t34 * t98 + t40 * t82 + t46 * t88) * t91 + t82 * t1 + t98 * t4 + t88 * t7) * t103; t51 * t3 + t72 * t6 + t57 * t9 + t50 * t2 + t71 * t5 + t56 * t8 + t49 * t1 + t70 * t4 + t55 * t7 - t121 * t176 + t122 * t175 + Ifges(4,3) * t174 - (-g(1) * t188 - g(2) * t187) * t162 + t163 * (g(1) * t187 - g(2) * t188) + ((t101 * t21 + t24 * t81 + t27 * t87) * t96 + (t102 * t21 + t24 * t84 + t27 * t90) * t93) * t105 + ((t20 * t99 + t23 * t80 + t26 * t86) * t95 + (t100 * t20 + t23 * t83 + t26 * t89) * t92) * t104 + ((t19 * t97 + t22 * t79 + t25 * t85) * t94 + (t19 * t98 + t22 * t82 + t25 * t88) * t91) * t103;];
tauX  = t52;
