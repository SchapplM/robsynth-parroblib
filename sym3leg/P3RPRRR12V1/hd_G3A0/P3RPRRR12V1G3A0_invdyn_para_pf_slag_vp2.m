% Calculate vector of inverse dynamics forces for parallel robot
% P3RPRRR12V1G3A0
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
% Datum: 2020-08-06 18:28
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRRR12V1G3A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:28:04
% EndTime: 2020-08-06 18:28:07
% DurationCPUTime: 3.31s
% Computational Cost: add. (8730->392), mult. (12501->678), div. (3243->7), fcn. (10260->18), ass. (0->232)
t139 = sin(qJ(3,3));
t103 = pkin(3) * t139 + qJ(2,3);
t100 = 0.1e1 / t103;
t151 = xDP(3);
t152 = xDP(2);
t153 = xDP(1);
t125 = 0.1e1 / t139;
t215 = t100 * t125;
t133 = legFrame(3,2);
t120 = cos(t133);
t145 = cos(qJ(3,3));
t128 = t145 ^ 2;
t140 = sin(qJ(1,3));
t248 = (-pkin(5) - pkin(6));
t123 = pkin(1) - t248;
t146 = cos(qJ(1,3));
t171 = qJ(2,3) * t140 + t123 * t146;
t117 = sin(t133);
t212 = t117 * t145;
t52 = t171 * t120 * t139 + qJ(2,3) * t212 + (t139 * t212 + (-t128 + 0.1e1) * t120 * t140) * pkin(3);
t206 = t120 * t145;
t55 = (pkin(3) * t206 - t171 * t117) * t139 + t140 * pkin(3) * (t145 - 0.1e1) * (t145 + 0.1e1) * t117 + qJ(2,3) * t206;
t82 = t103 * t146 - t123 * t140;
t260 = 0.2e1 * t82 * t100 * t151 + 0.2e1 * (t152 * t55 + t153 * t52) * t215;
t141 = sin(qJ(3,2));
t104 = pkin(3) * t141 + qJ(2,2);
t101 = 0.1e1 / t104;
t126 = 0.1e1 / t141;
t214 = t101 * t126;
t134 = legFrame(2,2);
t121 = cos(t134);
t147 = cos(qJ(3,2));
t129 = t147 ^ 2;
t142 = sin(qJ(1,2));
t148 = cos(qJ(1,2));
t172 = qJ(2,2) * t142 + t123 * t148;
t118 = sin(t134);
t210 = t118 * t147;
t53 = t172 * t121 * t141 + qJ(2,2) * t210 + (t141 * t210 + (-t129 + 0.1e1) * t121 * t142) * pkin(3);
t204 = t121 * t147;
t56 = (pkin(3) * t204 - t172 * t118) * t141 + t142 * pkin(3) * (t147 - 0.1e1) * (t147 + 0.1e1) * t118 + qJ(2,2) * t204;
t83 = t104 * t148 - t123 * t142;
t259 = 0.2e1 * t83 * t101 * t151 + 0.2e1 * (t152 * t56 + t153 * t53) * t214;
t143 = sin(qJ(3,1));
t105 = pkin(3) * t143 + qJ(2,1);
t102 = 0.1e1 / t105;
t127 = 0.1e1 / t143;
t213 = t102 * t127;
t135 = legFrame(1,2);
t122 = cos(t135);
t149 = cos(qJ(3,1));
t130 = t149 ^ 2;
t144 = sin(qJ(1,1));
t150 = cos(qJ(1,1));
t173 = qJ(2,1) * t144 + t123 * t150;
t119 = sin(t135);
t208 = t119 * t149;
t54 = t173 * t122 * t143 + qJ(2,1) * t208 + (t143 * t208 + (-t130 + 0.1e1) * t122 * t144) * pkin(3);
t202 = t122 * t149;
t57 = (pkin(3) * t202 - t173 * t119) * t143 + t144 * pkin(3) * (t149 - 0.1e1) * (t149 + 0.1e1) * t119 + qJ(2,1) * t202;
t84 = t105 * t150 - t123 * t144;
t258 = 0.2e1 * t84 * t102 * t151 + 0.2e1 * (t152 * t57 + t153 * t54) * t213;
t257 = 0.2e1 * t145;
t256 = 0.2e1 * t147;
t255 = 0.2e1 * t149;
t155 = m(2) + m(3);
t254 = qJ(2,1) * t155;
t253 = qJ(2,2) * t155;
t252 = qJ(2,3) * t155;
t251 = 2 * pkin(1);
t250 = -0.2e1 * pkin(3);
t249 = -2 * Ifges(3,4);
t163 = 0.1e1 / pkin(3);
t199 = t125 * t163;
t76 = (-t117 * t153 - t120 * t152) * t199;
t247 = pkin(3) * t76;
t197 = t126 * t163;
t77 = (-t118 * t153 - t121 * t152) * t197;
t246 = pkin(3) * t77;
t195 = t127 * t163;
t78 = (-t119 * t153 - t122 * t152) * t195;
t245 = pkin(3) * t78;
t244 = -mrSges(1,2) + mrSges(2,3);
t243 = (mrSges(3,3) - mrSges(2,2));
t242 = mrSges(3,1) * t139;
t241 = mrSges(3,1) * t141;
t240 = mrSges(3,1) * t143;
t67 = (-t140 * t151 + (-t117 * t152 + t120 * t153) * t146) * t100;
t239 = t123 * t67;
t68 = (-t142 * t151 + (-t118 * t152 + t121 * t153) * t148) * t101;
t238 = t123 * t68;
t69 = (-t144 * t151 + (-t119 * t152 + t122 * t153) * t150) * t102;
t237 = t123 * t69;
t154 = pkin(1) + pkin(5);
t109 = mrSges(3,2) * t154 - Ifges(3,6);
t110 = mrSges(3,1) * t154 - Ifges(3,5);
t79 = t109 * t139 - t110 * t145;
t93 = mrSges(3,1) * t145 - t139 * mrSges(3,2);
t236 = (-t140 * t79 + t82 * t93) * t215;
t235 = t125 * t52;
t234 = t125 * t55;
t96 = -m(2) * pkin(1) - m(3) * t154 - t243;
t64 = (-t140 * t96 + t155 * t82) * t100;
t233 = t125 * t64;
t232 = t125 * t93;
t231 = t125 * t96;
t80 = t109 * t141 - t110 * t147;
t94 = mrSges(3,1) * t147 - t141 * mrSges(3,2);
t230 = (-t142 * t80 + t83 * t94) * t214;
t229 = t126 * t53;
t228 = t126 * t56;
t65 = (-t142 * t96 + t155 * t83) * t101;
t227 = t126 * t65;
t226 = t126 * t94;
t225 = t126 * t96;
t81 = t109 * t143 - t110 * t149;
t95 = mrSges(3,1) * t149 - t143 * mrSges(3,2);
t224 = (-t144 * t81 + t84 * t95) * t213;
t223 = t127 * t54;
t222 = t127 * t57;
t66 = (-t144 * t96 + t155 * t84) * t102;
t221 = t127 * t66;
t220 = t127 * t95;
t219 = t127 * t96;
t218 = qJ(2,1) * t143;
t217 = qJ(2,2) * t141;
t216 = qJ(2,3) * t139;
t211 = t117 * t146;
t209 = t118 * t148;
t207 = t119 * t150;
t205 = t120 * t146;
t203 = t121 * t148;
t201 = t122 * t150;
t200 = t125 * t155;
t198 = t126 * t155;
t196 = t127 * t155;
t132 = Ifges(3,1) - Ifges(3,2);
t194 = t132 * t139;
t193 = t132 * t141;
t192 = t132 * t143;
t191 = t145 * t247;
t190 = t147 * t246;
t189 = t149 * t245;
t188 = mrSges(3,2) * t216;
t187 = mrSges(3,2) * t217;
t186 = mrSges(3,2) * t218;
t185 = t117 * t199;
t184 = t118 * t197;
t183 = t119 * t195;
t182 = t120 * t199;
t181 = t121 * t197;
t180 = t122 * t195;
t88 = g(1) * t120 - g(2) * t117;
t179 = g(3) * t146 + t140 * t88;
t89 = g(1) * t121 - g(2) * t118;
t178 = g(3) * t148 + t142 * t89;
t90 = g(1) * t122 - g(2) * t119;
t177 = g(3) * t150 + t144 * t90;
t176 = mrSges(3,2) * t145 + t242;
t175 = mrSges(3,2) * t147 + t241;
t174 = mrSges(3,2) * t149 + t240;
t170 = -t176 - t244 - t252;
t169 = -t175 - t244 - t253;
t168 = -t174 - t244 - t254;
t167 = mrSges(1,1) - t96;
t161 = pkin(5) ^ 2;
t162 = pkin(3) ^ 2;
t164 = pkin(1) ^ 2;
t166 = (t248 * t251) - t161 - t162 - t164 + ((-2 * pkin(5) - pkin(6)) * pkin(6));
t165 = m(3) * t161 + 2 * mrSges(3,3) * pkin(5) + (m(3) * pkin(5) + t243) * t251 + Ifges(2,1) + Ifges(3,2) + Ifges(1,3);
t159 = qJ(2,1) ^ 2;
t158 = qJ(2,2) ^ 2;
t157 = qJ(2,3) ^ 2;
t138 = xDDP(1);
t137 = xDDP(2);
t136 = xDDP(3);
t108 = mrSges(2,3) + t254;
t107 = mrSges(2,3) + t253;
t106 = mrSges(2,3) + t252;
t91 = t167 * g(3);
t87 = g(1) * t119 + g(2) * t122;
t86 = g(1) * t118 + g(2) * t121;
t85 = g(1) * t117 + g(2) * t120;
t75 = t78 ^ 2;
t74 = t77 ^ 2;
t73 = t76 ^ 2;
t63 = t69 ^ 2;
t62 = t68 ^ 2;
t61 = t67 ^ 2;
t60 = t132 * t130 + (t159 + t164) * t155 + 0.2e1 * (mrSges(2,3) + t240) * qJ(2,1) + (mrSges(3,2) * qJ(2,1) - Ifges(3,4) * t143) * t255 + t165;
t59 = t132 * t129 + (t158 + t164) * t155 + 0.2e1 * (mrSges(2,3) + t241) * qJ(2,2) + (mrSges(3,2) * qJ(2,2) - Ifges(3,4) * t141) * t256 + t165;
t58 = t132 * t128 + (t157 + t164) * t155 + 0.2e1 * (mrSges(2,3) + t242) * qJ(2,3) + (mrSges(3,2) * qJ(2,3) - Ifges(3,4) * t139) * t257 + t165;
t42 = t63 + t75;
t41 = t62 + t74;
t40 = t61 + t73;
t39 = (-t144 * t60 + t84 * t96) * t102;
t38 = (-t142 * t59 + t83 * t96) * t101;
t37 = (-t140 * t58 + t82 * t96) * t100;
t36 = -t95 * t180 + (t57 * t196 - t96 * t207) * t102;
t35 = -t94 * t181 + (t56 * t198 - t96 * t209) * t101;
t34 = -t93 * t182 + (t55 * t200 - t96 * t211) * t100;
t33 = -t95 * t183 + (t54 * t196 + t96 * t201) * t102;
t32 = -t94 * t184 + (t53 * t198 + t96 * t203) * t101;
t31 = -t93 * t185 + (t52 * t200 + t96 * t205) * t100;
t30 = -Ifges(3,3) * t180 + (-t81 * t207 + t57 * t220) * t102;
t29 = -Ifges(3,3) * t181 + (-t80 * t209 + t56 * t226) * t101;
t28 = -Ifges(3,3) * t182 + (-t79 * t211 + t55 * t232) * t100;
t27 = -Ifges(3,3) * t183 + (t81 * t201 + t54 * t220) * t102;
t26 = -Ifges(3,3) * t184 + (t80 * t203 + t53 * t226) * t101;
t25 = -Ifges(3,3) * t185 + (t79 * t205 + t52 * t232) * t100;
t21 = -t81 * t180 + (-t60 * t207 + t57 * t219) * t102;
t20 = -t80 * t181 + (-t59 * t209 + t56 * t225) * t101;
t19 = -t79 * t182 + (-t58 * t211 + t55 * t231) * t100;
t18 = -t81 * t183 + (t60 * t201 + t54 * t219) * t102;
t17 = -t80 * t184 + (t59 * t203 + t53 * t225) * t101;
t16 = -t79 * t185 + (t58 * t205 + t52 * t231) * t100;
t15 = (0.2e1 * t189 + t258 - t237) * t69 * t102;
t14 = (0.2e1 * t190 + t259 - t238) * t68 * t101;
t13 = (0.2e1 * t191 + t260 - t239) * t67 * t100;
t12 = (((t149 * t237 - t245) * t143 - qJ(2,1) * t78) * t127 * t245 + t237 * t258 + (t123 * t189 + (t130 * t162 + t218 * t250 - t159 + t166) * t69) * t69) * t102;
t11 = (((t147 * t238 - t246) * t141 - qJ(2,2) * t77) * t126 * t246 + t238 * t259 + (t123 * t190 + (t129 * t162 + t217 * t250 - t158 + t166) * t68) * t68) * t101;
t10 = (((t145 * t239 - t247) * t139 - qJ(2,3) * t76) * t125 * t247 + t239 * t260 + (t123 * t191 + (t128 * t162 + t216 * t250 - t157 + t166) * t67) * t67) * t100;
t9 = -t81 * t15 - t95 * t12 - (t130 * t249 + Ifges(3,4) - t186) * t63 + t143 * (mrSges(3,1) * t87 + t177 * mrSges(3,2)) + (-Ifges(3,3) * t75 * t127 + t63 * t192 + mrSges(3,2) * t87 + (-qJ(2,1) * t63 - t177) * mrSges(3,1)) * t149;
t8 = -t80 * t14 - t94 * t11 - (t129 * t249 + Ifges(3,4) - t187) * t62 + t141 * (mrSges(3,1) * t86 + t178 * mrSges(3,2)) + (-Ifges(3,3) * t74 * t126 + t62 * t193 + mrSges(3,2) * t86 + (-qJ(2,2) * t62 - t178) * mrSges(3,1)) * t147;
t7 = -t79 * t13 - t93 * t10 - (t128 * t249 + Ifges(3,4) - t188) * t61 + t139 * (mrSges(3,1) * t85 + t179 * mrSges(3,2)) + (-Ifges(3,3) * t73 * t125 + t61 * t194 + mrSges(3,2) * t85 + (-qJ(2,3) * t61 - t179) * mrSges(3,1)) * t145;
t6 = -t42 * t240 - t108 * t63 - t96 * t15 + (-t12 - t177) * t155 + (-mrSges(3,2) * t42 - t75 * t220) * t149;
t5 = -t41 * t241 - t107 * t62 - t96 * t14 + (-t11 - t178) * t155 + (-mrSges(3,2) * t41 - t74 * t226) * t147;
t4 = -t40 * t242 - t106 * t61 - t96 * t13 + (-t10 - t179) * t155 + (-mrSges(3,2) * t40 - t73 * t232) * t145;
t3 = -t96 * t12 - t60 * t15 + t91 * t150 + (t144 * t167 + t168 * t150) * t90 - t144 * t168 * g(3) + (t110 * t143 + (-t127 * t81 + t109) * t149) * t75 + ((t108 + t174) * t258 + (-0.2e1 * t186 + (-0.4e1 * t130 + 0.2e1) * Ifges(3,4) + (mrSges(3,1) * qJ(2,1) - t192) * t255) * t78) * t69;
t2 = -t96 * t11 - t59 * t14 + t91 * t148 + (t142 * t167 + t169 * t148) * t89 - t142 * t169 * g(3) + (t110 * t141 + (-t126 * t80 + t109) * t147) * t74 + ((t107 + t175) * t259 + (-0.2e1 * t187 + (-0.4e1 * t129 + 0.2e1) * Ifges(3,4) + (mrSges(3,1) * qJ(2,2) - t193) * t256) * t77) * t68;
t1 = -t96 * t10 - t58 * t13 + t91 * t146 + (t140 * t167 + t170 * t146) * t88 - t140 * t170 * g(3) + (t110 * t139 + (-t125 * t79 + t109) * t145) * t73 + ((t106 + t176) * t260 + (-0.2e1 * t188 + (-0.4e1 * t128 + 0.2e1) * Ifges(3,4) + (mrSges(3,1) * qJ(2,3) - t194) * t257) * t76) * t67;
t22 = [(-g(1) + t138) * m(4) + ((t18 * t201 + t33 * t223) * t138 + (-t18 * t207 + t33 * t222) * t137 + (-t144 * t18 + t33 * t84) * t136 + t3 * t201 + t6 * t223) * t102 + ((-t122 * t137 * t27 + (-t138 * t27 - t9) * t119) * t127 + (-t121 * t137 * t26 + (-t138 * t26 - t8) * t118) * t126 + (-t120 * t137 * t25 + (-t138 * t25 - t7) * t117) * t125) * t163 + ((t17 * t203 + t32 * t229) * t138 + (-t17 * t209 + t32 * t228) * t137 + (-t142 * t17 + t32 * t83) * t136 + t2 * t203 + t5 * t229) * t101 + ((t16 * t205 + t31 * t235) * t138 + (-t16 * t211 + t31 * t234) * t137 + (-t140 * t16 + t31 * t82) * t136 + t1 * t205 + t4 * t235) * t100; (-g(2) + t137) * m(4) + ((t21 * t201 + t36 * t223) * t138 + (-t21 * t207 + t36 * t222) * t137 + (-t144 * t21 + t36 * t84) * t136 - t3 * t207 + t6 * t222) * t102 + ((-t119 * t138 * t30 + (-t137 * t30 - t9) * t122) * t127 + (-t118 * t138 * t29 + (-t137 * t29 - t8) * t121) * t126 + (-t117 * t138 * t28 + (-t137 * t28 - t7) * t120) * t125) * t163 + ((t20 * t203 + t35 * t229) * t138 + (-t20 * t209 + t35 * t228) * t137 + (-t142 * t20 + t35 * t83) * t136 - t2 * t209 + t5 * t228) * t101 + ((t19 * t205 + t34 * t235) * t138 + (-t19 * t211 + t34 * t234) * t137 + (-t140 * t19 + t34 * t82) * t136 - t1 * t211 + t4 * t234) * t100; (-g(3) + t136) * m(4) + ((-t117 * t236 - t118 * t230 - t119 * t224) * t138 + (-t120 * t236 - t121 * t230 - t122 * t224) * t137) * t163 + ((t39 * t201 + t54 * t221) * t138 + (-t39 * t207 + t57 * t221) * t137 + (-t144 * t39 + t66 * t84) * t136 - t144 * t3 + t84 * t6) * t102 + ((t38 * t203 + t53 * t227) * t138 + (-t38 * t209 + t56 * t227) * t137 + (-t142 * t38 + t65 * t83) * t136 - t142 * t2 + t83 * t5) * t101 + ((t37 * t205 + t52 * t233) * t138 + (-t37 * t211 + t55 * t233) * t137 + (-t140 * t37 + t64 * t82) * t136 - t140 * t1 + t82 * t4) * t100;];
tauX  = t22;
