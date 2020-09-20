% Calculate vector of inverse dynamics forces for parallel robot
% P3RPRR1G3P3A0
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
% Datum: 2020-03-09 21:27
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRR1G3P3A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:26:49
% EndTime: 2020-03-09 21:26:51
% DurationCPUTime: 2.22s
% Computational Cost: add. (17346->339), mult. (12111->494), div. (2034->7), fcn. (8358->74), ass. (0->217)
t203 = 2 * pkin(2);
t156 = pkin(7) + qJ(3,3);
t170 = sin(qJ(3,3));
t93 = 0.1e1 / (pkin(1) * sin(t156) + t170 * pkin(2));
t227 = t93 / 0.2e1;
t157 = pkin(7) + qJ(3,2);
t172 = sin(qJ(3,2));
t94 = 0.1e1 / (pkin(1) * sin(t157) + t172 * pkin(2));
t226 = t94 / 0.2e1;
t158 = pkin(7) + qJ(3,1);
t174 = sin(qJ(3,1));
t95 = 0.1e1 / (pkin(1) * sin(t158) + t174 * pkin(2));
t225 = t95 / 0.2e1;
t183 = m(2) + m(3);
t248 = mrSges(3,1) * g(3);
t247 = mrSges(3,2) * g(3);
t246 = g(3) * mrSges(1,2);
t245 = g(3) * mrSges(2,2);
t162 = sin(pkin(7));
t244 = pkin(1) * t162;
t163 = cos(pkin(7));
t243 = pkin(2) * t163;
t181 = xDP(2);
t186 = 0.1e1 / pkin(3);
t207 = t186 / 0.2e1;
t201 = t181 * t207;
t159 = qJ(1,3) + pkin(7);
t164 = legFrame(3,2);
t119 = t164 + t159;
t120 = -t164 + t159;
t142 = qJ(1,3) + t164;
t143 = qJ(1,3) - t164;
t112 = qJ(3,3) + t119;
t113 = qJ(3,3) + t120;
t76 = -sin(t112) + sin(t113);
t46 = t76 * pkin(3) + (-sin(t119) + sin(t120)) * pkin(2) + (-sin(t142) + sin(t143)) * pkin(1);
t196 = t46 * t201;
t139 = qJ(1,3) + t156;
t125 = sin(t139);
t180 = xDP(3);
t212 = t125 * t180;
t182 = xDP(1);
t208 = t182 / 0.2e1;
t209 = t181 / 0.2e1;
t79 = cos(t113) + cos(t112);
t233 = (t208 * t79 + t209 * t76) * t93;
t200 = t182 * t207;
t49 = -t79 * pkin(3) + (-cos(t120) - cos(t119)) * pkin(2) + (-cos(t143) - cos(t142)) * pkin(1);
t43 = t49 * t93 * t200;
t206 = t180 * t186;
t133 = sin(t159);
t171 = sin(qJ(1,3));
t73 = pkin(1) * t171 + pkin(2) * t133 + pkin(3) * t125;
t239 = t73 * t93;
t67 = t206 * t239;
t236 = t43 + t67;
t16 = (-t196 - t212) * t93 + t233 + t236;
t25 = -t196 * t93 + t236;
t242 = t16 * t25;
t160 = qJ(1,2) + pkin(7);
t165 = legFrame(2,2);
t121 = t165 + t160;
t122 = -t165 + t160;
t144 = qJ(1,2) + t165;
t145 = qJ(1,2) - t165;
t114 = qJ(3,2) + t121;
t115 = qJ(3,2) + t122;
t77 = -sin(t114) + sin(t115);
t47 = t77 * pkin(3) + (-sin(t121) + sin(t122)) * pkin(2) + (-sin(t144) + sin(t145)) * pkin(1);
t195 = t47 * t201;
t140 = qJ(1,2) + t157;
t126 = sin(t140);
t211 = t126 * t180;
t80 = cos(t115) + cos(t114);
t232 = (t208 * t80 + t209 * t77) * t94;
t50 = -t80 * pkin(3) + (-cos(t122) - cos(t121)) * pkin(2) + (-cos(t145) - cos(t144)) * pkin(1);
t44 = t50 * t94 * t200;
t134 = sin(t160);
t173 = sin(qJ(1,2));
t74 = pkin(1) * t173 + pkin(2) * t134 + pkin(3) * t126;
t238 = t74 * t94;
t68 = t206 * t238;
t235 = t44 + t68;
t17 = (-t195 - t211) * t94 + t232 + t235;
t26 = -t195 * t94 + t235;
t241 = t17 * t26;
t161 = qJ(1,1) + pkin(7);
t166 = legFrame(1,2);
t123 = t166 + t161;
t124 = -t166 + t161;
t146 = qJ(1,1) + t166;
t147 = qJ(1,1) - t166;
t116 = qJ(3,1) + t123;
t117 = qJ(3,1) + t124;
t78 = -sin(t116) + sin(t117);
t48 = t78 * pkin(3) + (-sin(t123) + sin(t124)) * pkin(2) + (-sin(t146) + sin(t147)) * pkin(1);
t194 = t48 * t201;
t141 = qJ(1,1) + t158;
t127 = sin(t141);
t210 = t127 * t180;
t81 = cos(t117) + cos(t116);
t231 = (t208 * t81 + t209 * t78) * t95;
t51 = -t81 * pkin(3) + (-cos(t124) - cos(t123)) * pkin(2) + (-cos(t147) - cos(t146)) * pkin(1);
t45 = t51 * t95 * t200;
t135 = sin(t161);
t175 = sin(qJ(1,1));
t75 = pkin(1) * t175 + pkin(2) * t135 + pkin(3) * t127;
t237 = t75 * t95;
t69 = t206 * t237;
t234 = t45 + t69;
t18 = (-t194 - t210) * t95 + t231 + t234;
t27 = -t194 * t95 + t234;
t240 = t18 * t27;
t230 = t125 * t93;
t229 = t126 * t94;
t228 = t127 * t95;
t224 = t186 * t46;
t223 = t186 * t47;
t222 = t186 * t48;
t221 = t186 * t49;
t220 = t186 * t50;
t219 = t186 * t51;
t176 = cos(qJ(3,3));
t131 = pkin(1) * t163 + pkin(2);
t85 = mrSges(3,1) * t244 + mrSges(3,2) * t131;
t86 = mrSges(3,1) * t131 - mrSges(3,2) * t244;
t55 = -t170 * t85 + t176 * t86 + Ifges(3,3);
t218 = t186 * t55;
t177 = cos(qJ(3,2));
t56 = -t172 * t85 + t177 * t86 + Ifges(3,3);
t217 = t186 * t56;
t178 = cos(qJ(3,1));
t57 = -t174 * t85 + t178 * t86 + Ifges(3,3);
t216 = t186 * t57;
t215 = t186 * t73;
t214 = t186 * t74;
t213 = t186 * t75;
t154 = m(3) * pkin(2) + mrSges(2,1);
t205 = pkin(3) * t203;
t204 = 0.2e1 * pkin(1);
t148 = sin(t164);
t149 = sin(t165);
t150 = sin(t166);
t151 = cos(t164);
t152 = cos(t165);
t153 = cos(t166);
t202 = (t148 * t151 + t149 * t152 + t150 * t153) * t183;
t90 = g(1) * t151 - g(2) * t148;
t199 = (mrSges(3,2) * t90 + t248) * cos(t139) + t125 * (mrSges(3,1) * t90 - t247);
t91 = g(1) * t152 - g(2) * t149;
t198 = (mrSges(3,2) * t91 + t248) * cos(t140) + t126 * (mrSges(3,1) * t91 - t247);
t92 = g(1) * t153 - g(2) * t150;
t197 = (mrSges(3,2) * t92 + t248) * cos(t141) + t127 * (mrSges(3,1) * t92 - t247);
t193 = -t201 / 0.2e1;
t192 = mrSges(3,1) * t176 - mrSges(3,2) * t170;
t191 = mrSges(3,1) * t177 - mrSges(3,2) * t172;
t190 = mrSges(3,1) * t178 - mrSges(3,2) * t174;
t187 = pkin(2) ^ 2;
t188 = pkin(1) ^ 2;
t189 = (m(3) * t187) + t183 * t188 + Ifges(1,3) + Ifges(2,3) + Ifges(3,3);
t185 = pkin(3) ^ 2;
t169 = xDDP(1);
t168 = xDDP(2);
t167 = xDDP(3);
t155 = t187 + t188;
t138 = cos(t158);
t137 = cos(t157);
t136 = cos(t156);
t132 = g(3) * t154;
t118 = pkin(1) * t183 + mrSges(1,1);
t111 = g(3) * t118;
t89 = g(1) * t150 + g(2) * t153;
t88 = g(1) * t149 + g(2) * t152;
t87 = g(1) * t148 + g(2) * t151;
t60 = t174 * t86 + t178 * t85;
t59 = t172 * t86 + t177 * t85;
t58 = t170 * t86 + t176 * t85;
t54 = t190 * t203 + ((t190 + t154) * t163 - (mrSges(3,1) * t174 + mrSges(3,2) * t178 + mrSges(2,2)) * t162) * t204 + t189;
t53 = t191 * t203 + ((t191 + t154) * t163 - (mrSges(3,1) * t172 + mrSges(3,2) * t177 + mrSges(2,2)) * t162) * t204 + t189;
t52 = t192 * t203 + ((t192 + t154) * t163 - (mrSges(3,1) * t170 + mrSges(3,2) * t176 + mrSges(2,2)) * t162) * t204 + t189;
t42 = (Ifges(3,3) * t213 - t127 * t57) * t95;
t41 = (Ifges(3,3) * t214 - t126 * t56) * t94;
t40 = (Ifges(3,3) * t215 - t125 * t55) * t93;
t39 = -t210 * t95 + t231;
t38 = -t211 * t94 + t232;
t37 = -t212 * t93 + t233;
t36 = (-t127 * t54 + t213 * t57) * t95;
t35 = (-t126 * t53 + t214 * t56) * t94;
t34 = (-t125 * t52 + t215 * t55) * t93;
t33 = (Ifges(3,3) * t219 + t57 * t81) * t225;
t32 = (Ifges(3,3) * t220 + t56 * t80) * t226;
t31 = (Ifges(3,3) * t221 + t55 * t79) * t227;
t30 = (-Ifges(3,3) * t222 + t57 * t78) * t225;
t29 = (-Ifges(3,3) * t223 + t56 * t77) * t226;
t28 = (-Ifges(3,3) * t224 + t55 * t76) * t227;
t24 = (t216 * t51 + t54 * t81) * t225;
t23 = (t217 * t50 + t53 * t80) * t226;
t22 = (t218 * t49 + t52 * t79) * t227;
t21 = (-t216 * t48 + t54 * t78) * t225;
t20 = (-t217 * t47 + t53 * t77) * t226;
t19 = (-t218 * t46 + t52 * t76) * t227;
t15 = t45 / 0.2e1 + t69 / 0.2e1 + (t193 * t48 - t210) * t95 + t231;
t14 = t44 / 0.2e1 + t68 / 0.2e1 + (t193 * t47 - t211) * t94 + t232;
t13 = t43 / 0.2e1 + t67 / 0.2e1 + (t193 * t46 - t212) * t93 + t233;
t12 = (-pkin(3) * t240 + (-t18 * pkin(3) + (-pkin(1) * t138 - pkin(2) * t178) * t39) * t39) * t95;
t11 = (-pkin(3) * t241 + (-pkin(3) * t17 + (-pkin(1) * t137 - pkin(2) * t177) * t38) * t38) * t94;
t10 = (-pkin(3) * t242 + (-pkin(3) * t16 + (-pkin(1) * t136 - pkin(2) * t176) * t37) * t37) * t93;
t9 = (t15 * t178 * t205 + t39 * t155 + t18 * t185 + (pkin(3) * t138 * t15 + t243 * t39) * t204) * t186 * t95 * t39 + (t131 * t178 - t174 * t244 + pkin(3)) / (t131 * t174 + t178 * t244) * t240;
t8 = (t14 * t177 * t205 + t38 * t155 + t17 * t185 + (pkin(3) * t137 * t14 + t243 * t38) * t204) * t186 * t94 * t38 + (t131 * t177 - t172 * t244 + pkin(3)) / (t131 * t172 + t177 * t244) * t241;
t7 = (t13 * t176 * t205 + t37 * t155 + t16 * t185 + (pkin(3) * t13 * t136 + t243 * t37) * t204) * t186 * t93 * t37 + (t131 * t176 - t170 * t244 + pkin(3)) / (t131 * t170 + t176 * t244) * t242;
t6 = t39 ^ 2 * t60 - Ifges(3,3) * t9 - t12 * t57 + t197;
t5 = t38 ^ 2 * t59 - Ifges(3,3) * t8 - t11 * t56 + t198;
t4 = t37 ^ 2 * t58 - Ifges(3,3) * t7 - t55 * t10 + t199;
t3 = -t54 * t12 - t57 * t9 - 0.2e1 * t60 * t15 * t27 + (mrSges(2,2) * t92 + t132) * cos(t161) + (t154 * t92 - t245) * t135 + (mrSges(1,2) * t92 + t111) * cos(qJ(1,1)) + t175 * (t118 * t92 - t246) + t197;
t2 = -t53 * t11 - t56 * t8 - 0.2e1 * t59 * t14 * t26 + (mrSges(2,2) * t91 + t132) * cos(t160) + (t154 * t91 - t245) * t134 + (mrSges(1,2) * t91 + t111) * cos(qJ(1,2)) + t173 * (t118 * t91 - t246) + t198;
t1 = -t52 * t10 - t55 * t7 - 0.2e1 * t58 * t13 * t25 + (mrSges(2,2) * t90 + t132) * cos(t159) + (t154 * t90 - t245) * t133 + (mrSges(1,2) * t90 + t111) * cos(qJ(1,3)) + t171 * (t118 * t90 - t246) + t199;
t61 = [(-g(1) + t169) * m(4) + t202 * t168 + (-t22 * t230 - t23 * t229 - t24 * t228 + (t237 * t33 + t238 * t32 + t239 * t31) * t186) * t167 + (-t148 * t87 - t149 * t88 - t150 * t89 + (t148 ^ 2 + t149 ^ 2 + t150 ^ 2) * t169) * t183 + ((t219 * t33 + t24 * t81) * t169 + (-t222 * t33 + t24 * t78) * t168 + t81 * t3 + t6 * t219) * t225 + ((t220 * t32 + t23 * t80) * t169 + (-t223 * t32 + t23 * t77) * t168 + t80 * t2 + t5 * t220) * t226 + ((t22 * t79 + t221 * t31) * t169 + (t22 * t76 - t224 * t31) * t168 + t79 * t1 + t4 * t221) * t227; (-g(2) + t168) * m(4) + t202 * t169 + (-t19 * t230 - t20 * t229 - t21 * t228 + (t237 * t30 + t238 * t29 + t239 * t28) * t186) * t167 + (-t151 * t87 - t152 * t88 - t153 * t89 + (t151 ^ 2 + t152 ^ 2 + t153 ^ 2) * t168) * t183 + ((t21 * t81 + t219 * t30) * t169 + (t21 * t78 - t222 * t30) * t168 + t78 * t3 - t6 * t222) * t225 + ((t20 * t80 + t220 * t29) * t169 + (t20 * t77 - t223 * t29) * t168 + t77 * t2 - t5 * t223) * t226 + ((t19 * t79 + t221 * t28) * t169 + (t19 * t76 - t224 * t28) * t168 + t76 * t1 - t4 * t224) * t227; -t1 * t230 - t2 * t229 - t3 * t228 - g(3) * m(4) + (t6 * t237 + t5 * t238 + t4 * t239) * t186 + (-t34 * t230 - t35 * t229 - t36 * t228 + m(4) + (t237 * t42 + t238 * t41 + t239 * t40) * t186) * t167 + ((t219 * t42 + t36 * t81) * t169 + (-t222 * t42 + t36 * t78) * t168) * t225 + ((t220 * t41 + t35 * t80) * t169 + (-t223 * t41 + t35 * t77) * t168) * t226 + ((t221 * t40 + t34 * t79) * t169 + (-t224 * t40 + t34 * t76) * t168) * t227;];
tauX  = t61;
