% Calculate vector of inverse dynamics forces for parallel robot
% P3PRP1A0
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
% Datum: 2019-05-03 14:42
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRP1A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRP1A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRP1A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRP1A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRP1A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRP1A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PRP1A0_invdyn_para_pf_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRP1A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRP1A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRP1A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRP1A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRP1A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:41:40
% EndTime: 2019-05-03 14:41:43
% DurationCPUTime: 3.63s
% Computational Cost: add. (26085->404), mult. (49288->610), div. (2214->3), fcn. (23453->14), ass. (0->251)
t190 = (pkin(2) ^ 2);
t161 = 1 + t190;
t182 = (qJ(3,1) ^ 2);
t142 = -t182 + t161;
t173 = cos(qJ(2,1));
t159 = t173 ^ 2;
t170 = sin(qJ(2,1));
t227 = t170 * t173;
t217 = pkin(2) * t227;
t224 = t182 + t190;
t260 = 2 * qJ(3,1);
t105 = 0.1e1 / (t142 * t159 + t217 * t260 - t224 - 0.1e1);
t177 = xP(3);
t153 = sin(t177);
t154 = cos(t177);
t185 = koppelP(1,2);
t188 = koppelP(1,1);
t128 = -t153 * t185 + t154 * t188;
t174 = xDP(3);
t175 = xDP(2);
t108 = t128 * t174 + t175;
t164 = legFrame(1,3);
t148 = sin(t164);
t151 = cos(t164);
t241 = qJ(3,1) * t148;
t216 = pkin(2) * t241;
t201 = t151 * t182 - t216;
t240 = qJ(3,1) * t151;
t211 = pkin(2) * t240;
t202 = -t148 * t182 - t211;
t90 = t201 * t173 - t170 * (t148 - t202);
t245 = t90 * t108;
t125 = t153 * t188 + t154 * t185;
t176 = xDP(1);
t111 = -t125 * t174 + t176;
t87 = t202 * t173 + t170 * (-t151 - t201);
t248 = t87 * t111;
t270 = (t245 + t248) * t105;
t181 = (qJ(3,2) ^ 2);
t141 = -t181 + t161;
t172 = cos(qJ(2,2));
t158 = t172 ^ 2;
t169 = sin(qJ(2,2));
t228 = t169 * t172;
t218 = pkin(2) * t228;
t225 = t181 + t190;
t259 = 2 * qJ(3,2);
t104 = 0.1e1 / (t141 * t158 + t218 * t259 - t225 - 0.1e1);
t184 = koppelP(2,2);
t187 = koppelP(2,1);
t127 = -t153 * t184 + t154 * t187;
t107 = t127 * t174 + t175;
t163 = legFrame(2,3);
t147 = sin(t163);
t150 = cos(t163);
t237 = qJ(3,2) * t147;
t215 = pkin(2) * t237;
t198 = t150 * t181 - t215;
t236 = qJ(3,2) * t150;
t212 = pkin(2) * t236;
t199 = -t147 * t181 - t212;
t89 = t198 * t172 - t169 * (t147 - t199);
t246 = t89 * t107;
t124 = t153 * t187 + t154 * t184;
t110 = -t124 * t174 + t176;
t86 = t199 * t172 + t169 * (-t150 - t198);
t249 = t86 * t110;
t269 = (t246 + t249) * t104;
t180 = (qJ(3,3) ^ 2);
t140 = -t180 + t161;
t171 = cos(qJ(2,3));
t157 = t171 ^ 2;
t168 = sin(qJ(2,3));
t229 = t168 * t171;
t219 = pkin(2) * t229;
t226 = t180 + t190;
t258 = 2 * qJ(3,3);
t103 = 0.1e1 / (t140 * t157 + t219 * t258 - t226 - 0.1e1);
t183 = koppelP(3,2);
t186 = koppelP(3,1);
t126 = -t153 * t183 + t154 * t186;
t106 = t126 * t174 + t175;
t162 = legFrame(3,3);
t146 = sin(t162);
t149 = cos(t162);
t233 = qJ(3,3) * t146;
t214 = pkin(2) * t233;
t195 = t149 * t180 - t214;
t232 = qJ(3,3) * t149;
t213 = pkin(2) * t232;
t196 = -t146 * t180 - t213;
t88 = t195 * t171 - t168 * (t146 - t196);
t247 = t88 * t106;
t123 = t153 * t186 + t154 * t183;
t109 = -t123 * t174 + t176;
t85 = t196 * t171 + t168 * (-t149 - t195);
t250 = t85 * t109;
t268 = (t247 + t250) * t103;
t267 = 0.2e1 * pkin(2);
t230 = qJ(3,3) * t171;
t220 = -0.2e1 * t230;
t97 = t146 * t220 + t168 * (pkin(2) * t146 - t232);
t98 = t149 * t220 + t168 * (pkin(2) * t149 + t233);
t67 = (t106 * t97 + t109 * t98) * t103;
t266 = t67 ^ 2;
t234 = qJ(3,2) * t172;
t221 = -0.2e1 * t234;
t100 = t150 * t221 + t169 * (pkin(2) * t150 + t237);
t99 = t147 * t221 + t169 * (pkin(2) * t147 - t236);
t68 = (t100 * t110 + t107 * t99) * t104;
t265 = t68 ^ 2;
t238 = qJ(3,1) * t173;
t222 = -0.2e1 * t238;
t101 = t148 * t222 + t170 * (pkin(2) * t148 - t240);
t102 = t151 * t222 + t170 * (pkin(2) * t151 + t241);
t69 = (t101 * t108 + t102 * t111) * t105;
t264 = t69 ^ 2;
t263 = 0.2e1 * t67;
t262 = 0.2e1 * t68;
t261 = 0.2e1 * t69;
t257 = m(3) * pkin(2);
t256 = m(3) * t171;
t255 = m(3) * t172;
t254 = m(3) * t173;
t253 = t103 * t67;
t252 = t104 * t68;
t251 = t105 * t69;
t189 = pkin(2) * t190;
t58 = t180 * t67;
t209 = -t190 * t67 + t267 * t268 - t58;
t223 = -t190 - t161;
t10 = (t209 * t230 + ((t189 + (1 + t180) * pkin(2)) * t67 - t268 + t223 * t268) * t168) * t253;
t129 = -g(1) * t146 + g(2) * t149;
t244 = t10 + t129;
t59 = t181 * t68;
t208 = -t190 * t68 + t267 * t269 - t59;
t11 = (t208 * t234 + ((t189 + (1 + t181) * pkin(2)) * t68 - t269 + t223 * t269) * t169) * t252;
t130 = -g(1) * t147 + g(2) * t150;
t243 = t11 + t130;
t60 = t182 * t69;
t207 = -t190 * t69 + t267 * t270 - t60;
t12 = (t207 * t238 + ((t189 + (1 + t182) * pkin(2)) * t69 - t270 + t223 * t270) * t170) * t251;
t131 = -g(1) * t148 + g(2) * t151;
t242 = t12 + t131;
t239 = qJ(3,1) * t159;
t235 = qJ(3,2) * t158;
t231 = qJ(3,3) * t157;
t152 = mrSges(3,1) + t257;
t145 = m(3) * qJ(3,1) + mrSges(3,3);
t144 = m(3) * qJ(3,2) + mrSges(3,3);
t143 = m(3) * qJ(3,3) + mrSges(3,3);
t210 = mrSges(3,1) * t267 + Ifges(3,2) + Ifges(2,3);
t203 = -t142 * t148 + 0.2e1 * t211;
t200 = -t147 * t141 + 0.2e1 * t212;
t197 = -t140 * t146 + 0.2e1 * t213;
t179 = mrSges(4,1);
t178 = mrSges(4,2);
t167 = xDDP(1);
t166 = xDDP(2);
t165 = xDDP(3);
t160 = t174 ^ 2;
t155 = m(1) + m(2) + m(3);
t139 = mrSges(2,1) + t152;
t138 = -mrSges(2,2) + t145;
t137 = -mrSges(2,2) + t144;
t136 = -mrSges(2,2) + t143;
t135 = -t257 / 0.2e1 - mrSges(2,1) / 0.2e1 - mrSges(3,1) / 0.2e1;
t134 = g(1) * t151 + g(2) * t148;
t133 = g(1) * t150 + g(2) * t147;
t132 = g(1) * t149 + g(2) * t146;
t122 = t224 * m(3) + mrSges(3,3) * t260 + t210;
t121 = t225 * m(3) + mrSges(3,3) * t259 + t210;
t120 = t226 * m(3) + mrSges(3,3) * t258 + t210;
t119 = -t178 * t153 + t154 * t179;
t118 = t179 * t153 + t154 * t178;
t117 = t142 * t151 + 0.2e1 * t216;
t116 = t141 * t150 + 0.2e1 * t215;
t115 = t140 * t149 + 0.2e1 * t214;
t114 = t138 * t170 + t139 * t173;
t113 = t137 * t169 + t139 * t172;
t112 = t136 * t168 + t139 * t171;
t96 = -t125 * t165 - t128 * t160 + t167;
t95 = -t124 * t165 - t127 * t160 + t167;
t94 = -t123 * t165 - t126 * t160 + t167;
t93 = -t125 * t160 + t128 * t165 + t166;
t92 = -t124 * t160 + t127 * t165 + t166;
t91 = -t123 * t160 + t126 * t165 + t166;
t84 = -t117 * t227 + t190 * t148 + t203 * t159 + t148 - t211;
t83 = t117 * t159 - t151 * t190 + t203 * t227 - t151 - t216;
t82 = -t116 * t228 + t190 * t147 + t200 * t158 + t147 - t212;
t81 = t116 * t158 - t150 * t190 + t200 * t228 - t150 - t215;
t80 = -t115 * t229 + t190 * t146 + t197 * t157 + t146 - t213;
t79 = t115 * t157 - t149 * t190 + t197 * t229 - t149 - t214;
t72 = (t101 * t128 - t102 * t125) * t105;
t71 = (-t100 * t124 + t127 * t99) * t104;
t70 = (-t123 * t98 + t126 * t97) * t103;
t66 = pkin(2) * t68;
t65 = pkin(2) * t67;
t61 = t69 * pkin(2);
t57 = (-t125 * t87 + t128 * t90) * t105;
t56 = (-t124 * t86 + t127 * t89) * t104;
t55 = (-t123 * t85 + t126 * t88) * t103;
t51 = (-t125 * t84 + t128 * t83) * t105;
t50 = (-t124 * t82 + t127 * t81) * t104;
t49 = (-t123 * t80 + t126 * t79) * t103;
t48 = (-t102 * t152 + (-t173 * t84 + t87) * m(3)) * t105;
t47 = (-t101 * t152 + (-t173 * t83 + t90) * m(3)) * t105;
t46 = (-t100 * t152 + (-t172 * t82 + t86) * m(3)) * t104;
t45 = (-t152 * t99 + (-t172 * t81 + t89) * m(3)) * t104;
t44 = (-t152 * t98 + (-t171 * t80 + t85) * m(3)) * t103;
t43 = (-t152 * t97 + (-t171 * t79 + t88) * m(3)) * t103;
t42 = (t102 * t114 + t155 * t84 - t87 * t254) * t105;
t41 = (t101 * t114 + t155 * t83 - t90 * t254) * t105;
t40 = (t100 * t113 + t155 * t82 - t86 * t255) * t104;
t39 = (t113 * t99 + t155 * t81 - t89 * t255) * t104;
t38 = (t112 * t98 + t155 * t80 - t85 * t256) * t103;
t37 = (t112 * t97 + t155 * t79 - t88 * t256) * t103;
t36 = (t102 * t122 + t114 * t84 - t152 * t87) * t105;
t35 = (t101 * t122 + t114 * t83 - t152 * t90) * t105;
t34 = (t100 * t121 + t113 * t82 - t152 * t86) * t104;
t33 = (t113 * t81 + t121 * t99 - t152 * t89) * t104;
t32 = (t112 * t80 + t120 * t98 - t152 * t85) * t103;
t31 = (t112 * t79 + t120 * t97 - t152 * t88) * t103;
t30 = t66 - t269;
t29 = t65 - t268;
t28 = t61 - t270;
t27 = -t152 * t72 + (-t173 * t51 + t57) * m(3);
t26 = -t152 * t71 + (-t172 * t50 + t56) * m(3);
t25 = -t152 * t70 + (-t171 * t49 + t55) * m(3);
t24 = t114 * t72 + t155 * t51 - t57 * t254;
t23 = t113 * t71 + t155 * t50 - t56 * t255;
t22 = t112 * t70 + t155 * t49 - t55 * t256;
t21 = t114 * t51 + t122 * t72 - t152 * t57;
t20 = t113 * t50 + t121 * t71 - t152 * t56;
t19 = t112 * t49 + t120 * t70 - t152 * t55;
t18 = ((0.2e1 * (t66 + (-t249 / 0.2e1 - t246 / 0.2e1) * t104) * t235 - (pkin(2) * t30 - t59) * t228 - qJ(3,2) * t269) * t104 + (-qJ(3,2) + t218 - t235) * t104 * t269) * t68;
t17 = ((0.2e1 * (t65 + (-t250 / 0.2e1 - t247 / 0.2e1) * t103) * t231 - (pkin(2) * t29 - t58) * t229 - qJ(3,3) * t268) * t103 + (-qJ(3,3) + t219 - t231) * t103 * t268) * t67;
t16 = ((0.2e1 * (t61 + (-t248 / 0.2e1 - t245 / 0.2e1) * t105) * t239 - (pkin(2) * t28 - t60) * t227 - qJ(3,1) * t270) * t105 + (-qJ(3,1) + t217 - t239) * t105 * t270) * t69;
t15 = ((t28 - t270) * t227 + (-t159 * t69 - t207 + t69) * qJ(3,1)) * t251;
t14 = ((t30 - t269) * t228 + (-t158 * t68 - t208 + t68) * qJ(3,2)) * t252;
t13 = ((t29 - t268) * t229 + (-t157 * t67 - t209 + t67) * qJ(3,3)) * t253;
t9 = -t264 * t145 + t152 * t16 + (-t134 * t170 + t242 * t173 - t15) * m(3);
t8 = -t265 * t144 + t152 * t18 + (-t133 * t169 + t243 * t172 - t14) * m(3);
t7 = -t266 * t143 + t152 * t17 + (-t132 * t168 + t244 * t171 - t13) * m(3);
t6 = -t114 * t12 - t122 * t16 + t152 * t15 + t270 * t145 * t261 + (-t131 * t139 - t134 * t138) * t173 - (t131 * t138 - t134 * t139) * t170;
t5 = -t113 * t11 - t121 * t18 + t152 * t14 + t269 * t144 * t262 + (-t130 * t139 - t133 * t137) * t172 - (t130 * t137 - t133 * t139) * t169;
t4 = -t112 * t10 - t120 * t17 + t152 * t13 + t268 * t143 * t263 + (-t129 * t139 - t132 * t136) * t171 - (t129 * t136 - t132 * t139) * t168;
t3 = -t114 * t16 + (m(3) * t270 + t135 * t69) * t170 * t261 + (m(3) * t15 + t264 * t138) * t173 - t242 * t155;
t2 = -t113 * t18 + (m(3) * t269 + t135 * t68) * t169 * t262 + (m(3) * t14 + t265 * t137) * t172 - t243 * t155;
t1 = -t112 * t17 + (m(3) * t268 + t135 * t67) * t168 * t263 + (m(3) * t13 + t266 * t136) * t171 - t244 * t155;
t52 = [-t118 * t165 - t160 * t119 + (t167 - g(1)) * m(4) + ((t102 * t36 + t42 * t84 + t48 * t87) * t96 + (t101 * t36 + t42 * t83 + t48 * t90) * t93 + t84 * t3 + t102 * t6 + t87 * t9) * t105 + ((t100 * t34 + t40 * t82 + t46 * t86) * t95 + (t34 * t99 + t40 * t81 + t46 * t89) * t92 + t82 * t2 + t100 * t5 + t86 * t8) * t104 + ((t32 * t98 + t38 * t80 + t44 * t85) * t94 + (t32 * t97 + t38 * t79 + t44 * t88) * t91 + t80 * t1 + t98 * t4 + t85 * t7) * t103; -t160 * t118 + t119 * t165 + (t166 - g(2)) * m(4) + ((t102 * t35 + t41 * t84 + t47 * t87) * t96 + (t101 * t35 + t41 * t83 + t47 * t90) * t93 + t83 * t3 + t101 * t6 + t90 * t9) * t105 + ((t100 * t33 + t39 * t82 + t45 * t86) * t95 + (t33 * t99 + t39 * t81 + t45 * t89) * t92 + t81 * t2 + t99 * t5 + t89 * t8) * t104 + ((t31 * t98 + t37 * t80 + t43 * t85) * t94 + (t31 * t97 + t37 * t79 + t43 * t88) * t91 + t79 * t1 + t97 * t4 + t88 * t7) * t103; t51 * t3 + t72 * t6 + t57 * t9 + t50 * t2 + t71 * t5 + t56 * t8 + t49 * t1 + t70 * t4 + t55 * t7 - t118 * t167 + t119 * t166 + Ifges(4,3) * t165 - (-g(1) * t179 - g(2) * t178) * t153 + t154 * (g(1) * t178 - g(2) * t179) + ((t102 * t21 + t24 * t84 + t27 * t87) * t96 + (t101 * t21 + t24 * t83 + t27 * t90) * t93) * t105 + ((t100 * t20 + t23 * t82 + t26 * t86) * t95 + (t20 * t99 + t23 * t81 + t26 * t89) * t92) * t104 + ((t19 * t98 + t22 * t80 + t25 * t85) * t94 + (t19 * t97 + t22 * t79 + t25 * t88) * t91) * t103;];
tauX  = t52;
