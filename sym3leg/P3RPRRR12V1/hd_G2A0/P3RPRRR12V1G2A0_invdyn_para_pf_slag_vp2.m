% Calculate vector of inverse dynamics forces for parallel robot
% P3RPRRR12V1G2A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2020-08-06 18:25
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRRR12V1G2A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G2A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:24:32
% EndTime: 2020-08-06 18:24:35
% DurationCPUTime: 3.03s
% Computational Cost: add. (8730->392), mult. (12342->675), div. (3243->7), fcn. (10101->18), ass. (0->235)
t142 = sin(qJ(3,3));
t106 = t142 * pkin(3) + qJ(2,3);
t103 = 0.1e1 / t106;
t154 = xDP(3);
t155 = xDP(2);
t156 = xDP(1);
t128 = 0.1e1 / t142;
t213 = t103 * t128;
t136 = legFrame(3,2);
t120 = sin(t136);
t123 = cos(t136);
t148 = cos(qJ(3,3));
t191 = t148 * qJ(2,3);
t149 = cos(qJ(1,3));
t245 = pkin(3) * t149;
t246 = pkin(3) * t148;
t250 = (-pkin(5) - pkin(6));
t126 = pkin(1) - t250;
t143 = sin(qJ(1,3));
t91 = qJ(2,3) * t149 - t126 * t143;
t52 = (t120 * t246 - t91 * t123) * t142 + (t148 - 0.1e1) * (t148 + 0.1e1) * t123 * t245 + t120 * t191;
t131 = t148 ^ 2;
t58 = (t91 * t120 + t123 * t246) * t142 + (-t131 + 0.1e1) * t120 * t245 + t123 * t191;
t82 = t143 * t106 + t126 * t149;
t262 = 0.2e1 * t82 * t103 * t154 + 0.2e1 * (t155 * t58 + t156 * t52) * t213;
t144 = sin(qJ(3,2));
t107 = t144 * pkin(3) + qJ(2,2);
t104 = 0.1e1 / t107;
t129 = 0.1e1 / t144;
t212 = t104 * t129;
t137 = legFrame(2,2);
t124 = cos(t137);
t150 = cos(qJ(3,2));
t121 = sin(t137);
t208 = t121 * t150;
t151 = cos(qJ(1,2));
t244 = pkin(3) * t151;
t145 = sin(qJ(1,2));
t92 = qJ(2,2) * t151 - t126 * t145;
t53 = (pkin(3) * t208 - t92 * t124) * t144 + (t150 - 0.1e1) * (t150 + 0.1e1) * t124 * t244 + qJ(2,2) * t208;
t132 = t150 ^ 2;
t203 = t124 * t150;
t59 = (pkin(3) * t203 + t92 * t121) * t144 + (-t132 + 0.1e1) * t121 * t244 + qJ(2,2) * t203;
t83 = t145 * t107 + t126 * t151;
t261 = 0.2e1 * t83 * t104 * t154 + 0.2e1 * (t155 * t59 + t156 * t53) * t212;
t146 = sin(qJ(3,1));
t108 = t146 * pkin(3) + qJ(2,1);
t105 = 0.1e1 / t108;
t130 = 0.1e1 / t146;
t211 = t105 * t130;
t138 = legFrame(1,2);
t125 = cos(t138);
t152 = cos(qJ(3,1));
t122 = sin(t138);
t206 = t122 * t152;
t153 = cos(qJ(1,1));
t243 = pkin(3) * t153;
t147 = sin(qJ(1,1));
t93 = qJ(2,1) * t153 - t126 * t147;
t54 = (pkin(3) * t206 - t93 * t125) * t146 + (t152 - 0.1e1) * (t152 + 0.1e1) * t125 * t243 + qJ(2,1) * t206;
t133 = t152 ^ 2;
t201 = t125 * t152;
t60 = (pkin(3) * t201 + t93 * t122) * t146 + (-t133 + 0.1e1) * t122 * t243 + qJ(2,1) * t201;
t84 = t147 * t108 + t126 * t153;
t260 = 0.2e1 * t84 * t105 * t154 + 0.2e1 * (t155 * t60 + t156 * t54) * t211;
t259 = 0.2e1 * t148;
t258 = 0.2e1 * t150;
t257 = 0.2e1 * t152;
t158 = m(2) + m(3);
t256 = t158 * qJ(2,1);
t255 = t158 * qJ(2,2);
t254 = t158 * qJ(2,3);
t253 = 2 * pkin(1);
t252 = -0.2e1 * pkin(3);
t251 = -2 * Ifges(3,4);
t166 = 0.1e1 / pkin(3);
t199 = t128 * t166;
t76 = (-t120 * t156 - t123 * t155) * t199;
t249 = pkin(3) * t76;
t197 = t129 * t166;
t77 = (-t121 * t156 - t124 * t155) * t197;
t248 = pkin(3) * t77;
t195 = t130 * t166;
t78 = (-t122 * t156 - t125 * t155) * t195;
t247 = pkin(3) * t78;
t242 = -mrSges(1,2) + mrSges(2,3);
t241 = (mrSges(3,3) - mrSges(2,2));
t67 = (t149 * t154 + (-t120 * t155 + t123 * t156) * t143) * t103;
t240 = t126 * t67;
t68 = (t151 * t154 + (-t121 * t155 + t124 * t156) * t145) * t104;
t239 = t126 * t68;
t69 = (t153 * t154 + (-t122 * t155 + t125 * t156) * t147) * t105;
t238 = t126 * t69;
t157 = pkin(1) + pkin(5);
t112 = t157 * mrSges(3,2) - Ifges(3,6);
t113 = t157 * mrSges(3,1) - Ifges(3,5);
t79 = t112 * t142 - t113 * t148;
t96 = t148 * mrSges(3,1) - t142 * mrSges(3,2);
t237 = (t149 * t79 + t82 * t96) * t213;
t236 = t128 * t52;
t235 = t128 * t58;
t99 = -m(2) * pkin(1) - t157 * m(3) - t241;
t64 = (t149 * t99 + t158 * t82) * t103;
t234 = t128 * t64;
t233 = t128 * t96;
t232 = t128 * t99;
t80 = t112 * t144 - t113 * t150;
t97 = t150 * mrSges(3,1) - t144 * mrSges(3,2);
t231 = (t151 * t80 + t83 * t97) * t212;
t230 = t129 * t53;
t229 = t129 * t59;
t65 = (t151 * t99 + t158 * t83) * t104;
t228 = t129 * t65;
t227 = t129 * t97;
t226 = t129 * t99;
t81 = t112 * t146 - t113 * t152;
t98 = t152 * mrSges(3,1) - t146 * mrSges(3,2);
t225 = (t153 * t81 + t84 * t98) * t211;
t224 = t130 * t54;
t223 = t130 * t60;
t66 = (t153 * t99 + t158 * t84) * t105;
t222 = t130 * t66;
t221 = t130 * t98;
t220 = t130 * t99;
t219 = t142 * mrSges(3,1);
t218 = t144 * mrSges(3,1);
t217 = t146 * mrSges(3,1);
t216 = qJ(2,1) * t146;
t215 = qJ(2,2) * t144;
t214 = qJ(2,3) * t142;
t210 = t120 * t143;
t209 = t121 * t145;
t207 = t122 * t147;
t205 = t123 * t143;
t204 = t124 * t145;
t202 = t125 * t147;
t200 = t128 * t158;
t198 = t129 * t158;
t196 = t130 * t158;
t135 = Ifges(3,1) - Ifges(3,2);
t194 = t135 * t142;
t193 = t135 * t144;
t192 = t135 * t146;
t190 = t76 * t246;
t189 = t150 * t248;
t188 = t152 * t247;
t187 = mrSges(3,2) * t214;
t186 = mrSges(3,2) * t215;
t185 = mrSges(3,2) * t216;
t184 = t120 * t199;
t183 = t121 * t197;
t182 = t122 * t195;
t181 = t123 * t199;
t180 = t124 * t197;
t179 = t125 * t195;
t88 = t123 * g(1) - t120 * g(2);
t178 = g(3) * t143 - t149 * t88;
t89 = t124 * g(1) - t121 * g(2);
t177 = g(3) * t145 - t151 * t89;
t90 = t125 * g(1) - t122 * g(2);
t176 = g(3) * t147 - t153 * t90;
t175 = mrSges(3,2) * t148 + t219;
t174 = mrSges(3,2) * t150 + t218;
t173 = mrSges(3,2) * t152 + t217;
t172 = t175 + t242 + t254;
t171 = t174 + t242 + t255;
t170 = t173 + t242 + t256;
t95 = mrSges(1,1) - t99;
t164 = pkin(5) ^ 2;
t165 = pkin(3) ^ 2;
t167 = pkin(1) ^ 2;
t169 = (t250 * t253) - t164 - t165 - t167 + ((-2 * pkin(5) - pkin(6)) * pkin(6));
t168 = m(3) * t164 + 2 * mrSges(3,3) * pkin(5) + (m(3) * pkin(5) + t241) * t253 + Ifges(2,1) + Ifges(3,2) + Ifges(1,3);
t162 = qJ(2,1) ^ 2;
t161 = qJ(2,2) ^ 2;
t160 = qJ(2,3) ^ 2;
t141 = xDDP(1);
t140 = xDDP(2);
t139 = xDDP(3);
t111 = mrSges(2,3) + t256;
t110 = mrSges(2,3) + t255;
t109 = mrSges(2,3) + t254;
t94 = t95 * g(3);
t87 = t122 * g(1) + t125 * g(2);
t86 = t121 * g(1) + t124 * g(2);
t85 = t120 * g(1) + t123 * g(2);
t75 = t78 ^ 2;
t74 = t77 ^ 2;
t73 = t76 ^ 2;
t63 = t69 ^ 2;
t62 = t68 ^ 2;
t61 = t67 ^ 2;
t57 = t135 * t133 + (t162 + t167) * t158 + 0.2e1 * (mrSges(2,3) + t217) * qJ(2,1) + (mrSges(3,2) * qJ(2,1) - Ifges(3,4) * t146) * t257 + t168;
t56 = t135 * t132 + (t161 + t167) * t158 + 0.2e1 * (mrSges(2,3) + t218) * qJ(2,2) + (mrSges(3,2) * qJ(2,2) - Ifges(3,4) * t144) * t258 + t168;
t55 = t135 * t131 + (t160 + t167) * t158 + 0.2e1 * (mrSges(2,3) + t219) * qJ(2,3) + (mrSges(3,2) * qJ(2,3) - Ifges(3,4) * t142) * t259 + t168;
t42 = t63 + t75;
t41 = t62 + t74;
t40 = t61 + t73;
t39 = (t153 * t57 + t84 * t99) * t105;
t38 = (t151 * t56 + t83 * t99) * t104;
t37 = (t149 * t55 + t82 * t99) * t103;
t36 = -t98 * t179 + (t60 * t196 - t99 * t207) * t105;
t35 = -t97 * t180 + (t59 * t198 - t99 * t209) * t104;
t34 = -t96 * t181 + (t58 * t200 - t99 * t210) * t103;
t33 = -t98 * t182 + (t54 * t196 + t99 * t202) * t105;
t32 = -t97 * t183 + (t53 * t198 + t99 * t204) * t104;
t31 = -t96 * t184 + (t52 * t200 + t99 * t205) * t103;
t30 = -Ifges(3,3) * t179 + (-t81 * t207 + t60 * t221) * t105;
t29 = -Ifges(3,3) * t180 + (-t80 * t209 + t59 * t227) * t104;
t28 = -Ifges(3,3) * t181 + (-t79 * t210 + t58 * t233) * t103;
t27 = -Ifges(3,3) * t182 + (t81 * t202 + t54 * t221) * t105;
t26 = -Ifges(3,3) * t183 + (t80 * t204 + t53 * t227) * t104;
t25 = -Ifges(3,3) * t184 + (t79 * t205 + t52 * t233) * t103;
t21 = -t81 * t179 + (-t57 * t207 + t60 * t220) * t105;
t20 = -t80 * t180 + (-t56 * t209 + t59 * t226) * t104;
t19 = -t79 * t181 + (-t55 * t210 + t58 * t232) * t103;
t18 = -t81 * t182 + (t57 * t202 + t54 * t220) * t105;
t17 = -t80 * t183 + (t56 * t204 + t53 * t226) * t104;
t16 = -t79 * t184 + (t55 * t205 + t52 * t232) * t103;
t15 = (0.2e1 * t188 + t260 - t238) * t69 * t105;
t14 = (0.2e1 * t189 + t261 - t239) * t68 * t104;
t13 = (0.2e1 * t190 + t262 - t240) * t67 * t103;
t12 = (((t152 * t238 - t247) * t146 - qJ(2,1) * t78) * t130 * t247 + t238 * t260 + (t126 * t188 + (t133 * t165 + t216 * t252 - t162 + t169) * t69) * t69) * t105;
t11 = (((t150 * t239 - t248) * t144 - qJ(2,2) * t77) * t129 * t248 + t239 * t261 + (t126 * t189 + (t132 * t165 + t215 * t252 - t161 + t169) * t68) * t68) * t104;
t10 = (((t148 * t240 - t249) * t142 - qJ(2,3) * t76) * t128 * t249 + t240 * t262 + (t126 * t190 + (t131 * t165 + t214 * t252 - t160 + t169) * t67) * t67) * t103;
t9 = -t81 * t15 - t98 * t12 - (t133 * t251 + Ifges(3,4) - t185) * t63 + t146 * (mrSges(3,1) * t87 + t176 * mrSges(3,2)) + (-Ifges(3,3) * t75 * t130 + t63 * t192 + mrSges(3,2) * t87 + (-qJ(2,1) * t63 - t176) * mrSges(3,1)) * t152;
t8 = -t80 * t14 - t97 * t11 - (t132 * t251 + Ifges(3,4) - t186) * t62 + t144 * (mrSges(3,1) * t86 + t177 * mrSges(3,2)) + (-Ifges(3,3) * t74 * t129 + t62 * t193 + mrSges(3,2) * t86 + (-qJ(2,2) * t62 - t177) * mrSges(3,1)) * t150;
t7 = -t79 * t13 - t96 * t10 - (t131 * t251 + Ifges(3,4) - t187) * t61 + t142 * (mrSges(3,1) * t85 + t178 * mrSges(3,2)) + (-Ifges(3,3) * t73 * t128 + t61 * t194 + mrSges(3,2) * t85 + (-qJ(2,3) * t61 - t178) * mrSges(3,1)) * t148;
t6 = -t42 * t217 - t111 * t63 - t99 * t15 + (-t12 - t176) * t158 + (-mrSges(3,2) * t42 - t75 * t221) * t152;
t5 = -t41 * t218 - t110 * t62 - t99 * t14 + (-t11 - t177) * t158 + (-mrSges(3,2) * t41 - t74 * t227) * t150;
t4 = -t40 * t219 - t109 * t61 - t99 * t13 + (-t10 - t178) * t158 + (-mrSges(3,2) * t40 - t73 * t233) * t148;
t3 = -t99 * t12 + t94 * t147 - t57 * t15 + (-t170 * t147 - t95 * t153) * t90 - t170 * t153 * g(3) + (t113 * t146 + (-t130 * t81 + t112) * t152) * t75 + ((t111 + t173) * t260 + (-0.2e1 * t185 + (-0.4e1 * t133 + 0.2e1) * Ifges(3,4) + (mrSges(3,1) * qJ(2,1) - t192) * t257) * t78) * t69;
t2 = -t99 * t11 - t56 * t14 + t94 * t145 + (-t171 * t145 - t95 * t151) * t89 - t171 * t151 * g(3) + (t113 * t144 + (-t129 * t80 + t112) * t150) * t74 + ((t110 + t174) * t261 + (-0.2e1 * t186 + (-0.4e1 * t132 + 0.2e1) * Ifges(3,4) + (mrSges(3,1) * qJ(2,2) - t193) * t258) * t77) * t68;
t1 = -t99 * t10 - t55 * t13 + t94 * t143 + (-t172 * t143 - t95 * t149) * t88 - t172 * t149 * g(3) + (t113 * t142 + (-t128 * t79 + t112) * t148) * t73 + ((t109 + t175) * t262 + (-0.2e1 * t187 + (-0.4e1 * t131 + 0.2e1) * Ifges(3,4) + (mrSges(3,1) * qJ(2,3) - t194) * t259) * t76) * t67;
t22 = [(-g(1) + t141) * m(4) + ((t18 * t202 + t33 * t224) * t141 + (-t18 * t207 + t33 * t223) * t140 + (t153 * t18 + t33 * t84) * t139 + t3 * t202 + t6 * t224) * t105 + ((-t125 * t140 * t27 + (-t141 * t27 - t9) * t122) * t130 + (-t124 * t140 * t26 + (-t141 * t26 - t8) * t121) * t129 + (-t123 * t140 * t25 + (-t141 * t25 - t7) * t120) * t128) * t166 + ((t17 * t204 + t32 * t230) * t141 + (-t17 * t209 + t32 * t229) * t140 + (t151 * t17 + t32 * t83) * t139 + t2 * t204 + t5 * t230) * t104 + ((t16 * t205 + t31 * t236) * t141 + (-t16 * t210 + t31 * t235) * t140 + (t149 * t16 + t31 * t82) * t139 + t1 * t205 + t4 * t236) * t103; (-g(2) + t140) * m(4) + ((t21 * t202 + t36 * t224) * t141 + (-t21 * t207 + t36 * t223) * t140 + (t153 * t21 + t36 * t84) * t139 - t3 * t207 + t6 * t223) * t105 + ((-t122 * t141 * t30 + (-t140 * t30 - t9) * t125) * t130 + (-t121 * t141 * t29 + (-t140 * t29 - t8) * t124) * t129 + (-t120 * t141 * t28 + (-t140 * t28 - t7) * t123) * t128) * t166 + ((t20 * t204 + t35 * t230) * t141 + (-t20 * t209 + t35 * t229) * t140 + (t151 * t20 + t35 * t83) * t139 - t2 * t209 + t5 * t229) * t104 + ((t19 * t205 + t34 * t236) * t141 + (-t19 * t210 + t34 * t235) * t140 + (t149 * t19 + t34 * t82) * t139 - t1 * t210 + t4 * t235) * t103; (-g(3) + t139) * m(4) + ((-t120 * t237 - t121 * t231 - t122 * t225) * t141 + (-t123 * t237 - t124 * t231 - t125 * t225) * t140) * t166 + ((t39 * t202 + t54 * t222) * t141 + (-t39 * t207 + t60 * t222) * t140 + (t153 * t39 + t66 * t84) * t139 + t153 * t3 + t84 * t6) * t105 + ((t38 * t204 + t53 * t228) * t141 + (-t38 * t209 + t59 * t228) * t140 + (t151 * t38 + t65 * t83) * t139 + t151 * t2 + t83 * t5) * t104 + ((t37 * t205 + t52 * t234) * t141 + (-t37 * t210 + t58 * t234) * t140 + (t149 * t37 + t64 * t82) * t139 + t149 * t1 + t82 * t4) * t103;];
tauX  = t22;
