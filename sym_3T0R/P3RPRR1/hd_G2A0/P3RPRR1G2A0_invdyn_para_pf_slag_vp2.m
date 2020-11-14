% Calculate vector of inverse dynamics forces for parallel robot
% P3RPRR1G2A0
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
% Datum: 2020-03-09 21:25
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRR1G2A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G2A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRR1G2A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRR1G2A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G2A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRR1G2A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G2A0_invdyn_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRR1G2A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRR1G2A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRR1G2A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G2A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G2A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:24:53
% EndTime: 2020-03-09 21:24:55
% DurationCPUTime: 2.19s
% Computational Cost: add. (17346->327), mult. (12111->481), div. (2034->7), fcn. (8358->74), ass. (0->207)
t206 = 2 * pkin(2);
t156 = pkin(7) + qJ(3,3);
t170 = sin(qJ(3,3));
t96 = 0.1e1 / (pkin(1) * sin(t156) + t170 * pkin(2));
t227 = t96 / 0.2e1;
t157 = pkin(7) + qJ(3,2);
t171 = sin(qJ(3,2));
t97 = 0.1e1 / (pkin(1) * sin(t157) + t171 * pkin(2));
t226 = t97 / 0.2e1;
t158 = pkin(7) + qJ(3,1);
t172 = sin(qJ(3,1));
t98 = 0.1e1 / (pkin(1) * sin(t158) + t172 * pkin(2));
t225 = t98 / 0.2e1;
t187 = m(2) + m(3);
t162 = sin(pkin(7));
t238 = pkin(1) * t162;
t163 = cos(pkin(7));
t237 = pkin(2) * t163;
t186 = xDP(1);
t190 = 0.1e1 / pkin(3);
t210 = t190 / 0.2e1;
t203 = t186 * t210;
t159 = qJ(1,3) + pkin(7);
t164 = legFrame(3,2);
t122 = t164 + t159;
t123 = -t164 + t159;
t142 = qJ(1,3) + t164;
t143 = qJ(1,3) - t164;
t115 = qJ(3,3) + t122;
t116 = qJ(3,3) + t123;
t79 = sin(t115) + sin(t116);
t46 = t79 * pkin(3) + (sin(t122) + sin(t123)) * pkin(2) + (sin(t142) + sin(t143)) * pkin(1);
t199 = t46 * t96 * t203;
t185 = xDP(2);
t204 = t185 * t210;
t82 = -cos(t116) + cos(t115);
t49 = -t82 * pkin(3) + (cos(t123) - cos(t122)) * pkin(2) + (cos(t143) - cos(t142)) * pkin(1);
t43 = t49 * t96 * t204;
t184 = xDP(3);
t209 = t184 * t190;
t139 = qJ(1,3) + t156;
t128 = cos(t139);
t136 = cos(t159);
t174 = cos(qJ(1,3));
t76 = -t174 * pkin(1) - pkin(2) * t136 - pkin(3) * t128;
t233 = t76 * t96;
t67 = t209 * t233;
t25 = t43 + t67 - t199;
t211 = t186 / 0.2e1;
t212 = t185 / 0.2e1;
t230 = t128 * t96;
t37 = t184 * t230 + (t211 * t79 + t212 * t82) * t96;
t16 = t25 + t37;
t236 = t16 * t25;
t160 = qJ(1,2) + pkin(7);
t165 = legFrame(2,2);
t124 = t165 + t160;
t125 = -t165 + t160;
t144 = qJ(1,2) + t165;
t145 = qJ(1,2) - t165;
t117 = qJ(3,2) + t124;
t118 = qJ(3,2) + t125;
t80 = sin(t117) + sin(t118);
t47 = t80 * pkin(3) + (sin(t124) + sin(t125)) * pkin(2) + (sin(t144) + sin(t145)) * pkin(1);
t198 = t47 * t97 * t203;
t83 = -cos(t118) + cos(t117);
t50 = -t83 * pkin(3) + (cos(t125) - cos(t124)) * pkin(2) + (cos(t145) - cos(t144)) * pkin(1);
t44 = t50 * t97 * t204;
t140 = qJ(1,2) + t157;
t129 = cos(t140);
t137 = cos(t160);
t176 = cos(qJ(1,2));
t77 = -t176 * pkin(1) - pkin(2) * t137 - pkin(3) * t129;
t232 = t77 * t97;
t68 = t209 * t232;
t26 = t44 + t68 - t198;
t229 = t129 * t97;
t38 = t184 * t229 + (t211 * t80 + t212 * t83) * t97;
t17 = t26 + t38;
t235 = t17 * t26;
t161 = qJ(1,1) + pkin(7);
t166 = legFrame(1,2);
t126 = t166 + t161;
t127 = -t166 + t161;
t146 = qJ(1,1) + t166;
t147 = qJ(1,1) - t166;
t119 = qJ(3,1) + t126;
t120 = qJ(3,1) + t127;
t81 = sin(t119) + sin(t120);
t48 = t81 * pkin(3) + (sin(t126) + sin(t127)) * pkin(2) + (sin(t146) + sin(t147)) * pkin(1);
t197 = t48 * t98 * t203;
t84 = -cos(t120) + cos(t119);
t51 = -t84 * pkin(3) + (cos(t127) - cos(t126)) * pkin(2) + (cos(t147) - cos(t146)) * pkin(1);
t45 = t51 * t98 * t204;
t141 = qJ(1,1) + t158;
t130 = cos(t141);
t138 = cos(t161);
t178 = cos(qJ(1,1));
t78 = -t178 * pkin(1) - pkin(2) * t138 - pkin(3) * t130;
t231 = t78 * t98;
t69 = t209 * t231;
t27 = t45 + t69 - t197;
t228 = t130 * t98;
t39 = t184 * t228 + (t211 * t81 + t212 * t84) * t98;
t18 = t27 + t39;
t234 = t18 * t27;
t224 = t190 * t46;
t223 = t190 * t47;
t222 = t190 * t48;
t221 = t190 * t49;
t220 = t190 * t50;
t219 = t190 * t51;
t173 = cos(qJ(3,3));
t131 = t163 * pkin(1) + pkin(2);
t88 = mrSges(3,1) * t238 + t131 * mrSges(3,2);
t89 = t131 * mrSges(3,1) - mrSges(3,2) * t238;
t55 = -t88 * t170 + t89 * t173 + Ifges(3,3);
t218 = t190 * t55;
t175 = cos(qJ(3,2));
t56 = -t88 * t171 + t89 * t175 + Ifges(3,3);
t217 = t190 * t56;
t177 = cos(qJ(3,1));
t57 = -t88 * t172 + t89 * t177 + Ifges(3,3);
t216 = t190 * t57;
t215 = t190 * t76;
t214 = t190 * t77;
t213 = t190 * t78;
t154 = m(3) * pkin(2) + mrSges(2,1);
t208 = pkin(3) * t206;
t207 = 0.2e1 * pkin(1);
t148 = sin(t164);
t149 = sin(t165);
t150 = sin(t166);
t151 = cos(t164);
t152 = cos(t165);
t153 = cos(t166);
t205 = (t148 * t151 + t149 * t152 + t150 * t153) * t187;
t180 = mrSges(3,2) * g(3);
t182 = mrSges(3,1) * g(3);
t93 = t151 * g(1) - t148 * g(2);
t202 = -t128 * (mrSges(3,1) * t93 - t180) + sin(t139) * (mrSges(3,2) * t93 + t182);
t94 = t152 * g(1) - t149 * g(2);
t201 = -t129 * (mrSges(3,1) * t94 - t180) + sin(t140) * (mrSges(3,2) * t94 + t182);
t95 = t153 * g(1) - t150 * g(2);
t200 = -t130 * (mrSges(3,1) * t95 - t180) + sin(t141) * (mrSges(3,2) * t95 + t182);
t196 = t173 * mrSges(3,1) - mrSges(3,2) * t170;
t195 = t175 * mrSges(3,1) - mrSges(3,2) * t171;
t194 = t177 * mrSges(3,1) - mrSges(3,2) * t172;
t191 = pkin(2) ^ 2;
t192 = pkin(1) ^ 2;
t193 = (m(3) * t191) + t187 * t192 + Ifges(1,3) + Ifges(2,3) + Ifges(3,3);
t189 = pkin(3) ^ 2;
t181 = mrSges(1,2) * g(3);
t179 = g(3) * mrSges(2,2);
t169 = xDDP(1);
t168 = xDDP(2);
t167 = xDDP(3);
t155 = t191 + t192;
t135 = cos(t158);
t134 = cos(t157);
t133 = cos(t156);
t132 = t154 * g(3);
t121 = t187 * pkin(1) + mrSges(1,1);
t114 = t121 * g(3);
t92 = t150 * g(1) + t153 * g(2);
t91 = t149 * g(1) + t152 * g(2);
t90 = t148 * g(1) + t151 * g(2);
t60 = t172 * t89 + t88 * t177;
t59 = t171 * t89 + t88 * t175;
t58 = t170 * t89 + t88 * t173;
t54 = t194 * t206 + ((t194 + t154) * t163 - (mrSges(3,1) * t172 + t177 * mrSges(3,2) + mrSges(2,2)) * t162) * t207 + t193;
t53 = t195 * t206 + ((t195 + t154) * t163 - (mrSges(3,1) * t171 + t175 * mrSges(3,2) + mrSges(2,2)) * t162) * t207 + t193;
t52 = t196 * t206 + ((t196 + t154) * t163 - (mrSges(3,1) * t170 + t173 * mrSges(3,2) + mrSges(2,2)) * t162) * t207 + t193;
t42 = (Ifges(3,3) * t213 + t130 * t57) * t98;
t41 = (Ifges(3,3) * t214 + t129 * t56) * t97;
t40 = (Ifges(3,3) * t215 + t128 * t55) * t96;
t36 = (t130 * t54 + t57 * t213) * t98;
t35 = (t129 * t53 + t56 * t214) * t97;
t34 = (t128 * t52 + t55 * t215) * t96;
t33 = (Ifges(3,3) * t219 + t57 * t84) * t225;
t32 = (Ifges(3,3) * t220 + t56 * t83) * t226;
t31 = (Ifges(3,3) * t221 + t55 * t82) * t227;
t30 = (-Ifges(3,3) * t222 + t57 * t81) * t225;
t29 = (-Ifges(3,3) * t223 + t56 * t80) * t226;
t28 = (-Ifges(3,3) * t224 + t55 * t79) * t227;
t24 = (t51 * t216 + t54 * t84) * t225;
t23 = (t50 * t217 + t53 * t83) * t226;
t22 = (t49 * t218 + t52 * t82) * t227;
t21 = (-t48 * t216 + t54 * t81) * t225;
t20 = (-t47 * t217 + t53 * t80) * t226;
t19 = (-t46 * t218 + t52 * t79) * t227;
t15 = -t197 / 0.2e1 + t45 / 0.2e1 + t69 / 0.2e1 + t39;
t14 = -t198 / 0.2e1 + t44 / 0.2e1 + t68 / 0.2e1 + t38;
t13 = -t199 / 0.2e1 + t43 / 0.2e1 + t67 / 0.2e1 + t37;
t12 = (-pkin(3) * t234 + (-t18 * pkin(3) + (-pkin(1) * t135 - pkin(2) * t177) * t39) * t39) * t98;
t11 = (-pkin(3) * t235 + (-pkin(3) * t17 + (-pkin(1) * t134 - pkin(2) * t175) * t38) * t38) * t97;
t10 = (-pkin(3) * t236 + (-pkin(3) * t16 + (-pkin(1) * t133 - pkin(2) * t173) * t37) * t37) * t96;
t9 = (t15 * t177 * t208 + t39 * t155 + t18 * t189 + (pkin(3) * t135 * t15 + t39 * t237) * t207) * t190 * t98 * t39 + (t177 * t131 - t172 * t238 + pkin(3)) / (t172 * t131 + t177 * t238) * t234;
t8 = (t14 * t175 * t208 + t38 * t155 + t17 * t189 + (pkin(3) * t134 * t14 + t38 * t237) * t207) * t190 * t97 * t38 + (t131 * t175 - t171 * t238 + pkin(3)) / (t131 * t171 + t175 * t238) * t235;
t7 = (t13 * t173 * t208 + t37 * t155 + t16 * t189 + (pkin(3) * t13 * t133 + t37 * t237) * t207) * t190 * t96 * t37 + (t131 * t173 - t170 * t238 + pkin(3)) / (t131 * t170 + t173 * t238) * t236;
t6 = t39 ^ 2 * t60 - Ifges(3,3) * t9 - t57 * t12 + t200;
t5 = t38 ^ 2 * t59 - Ifges(3,3) * t8 - t56 * t11 + t201;
t4 = t37 ^ 2 * t58 - Ifges(3,3) * t7 - t55 * t10 + t202;
t3 = -t54 * t12 - t57 * t9 - 0.2e1 * t60 * t15 * t27 + (-t95 * t154 + t179) * t138 + (t95 * mrSges(2,2) + t132) * sin(t161) + (-t95 * t121 + t181) * t178 + sin(qJ(1,1)) * (t95 * mrSges(1,2) + t114) + t200;
t2 = -t53 * t11 - t56 * t8 - 0.2e1 * t59 * t14 * t26 + (-t94 * t154 + t179) * t137 + (t94 * mrSges(2,2) + t132) * sin(t160) + (-t94 * t121 + t181) * t176 + sin(qJ(1,2)) * (t94 * mrSges(1,2) + t114) + t201;
t1 = -t52 * t10 - t55 * t7 - 0.2e1 * t58 * t13 * t25 + (-t93 * t154 + t179) * t136 + (t93 * mrSges(2,2) + t132) * sin(t159) + (-t93 * t121 + t181) * t174 + sin(qJ(1,3)) * (t93 * mrSges(1,2) + t114) + t202;
t61 = [(-g(1) + t169) * m(4) + t205 * t168 + (t19 * t230 + t20 * t229 + t21 * t228 + (t30 * t231 + t29 * t232 + t28 * t233) * t190) * t167 + (-t148 * t90 - t149 * t91 - t150 * t92 + (t148 ^ 2 + t149 ^ 2 + t150 ^ 2) * t169) * t187 + ((t21 * t81 - t30 * t222) * t169 + (t21 * t84 + t30 * t219) * t168 + t81 * t3 - t6 * t222) * t225 + ((t20 * t80 - t29 * t223) * t169 + (t20 * t83 + t29 * t220) * t168 + t80 * t2 - t5 * t223) * t226 + ((t19 * t79 - t28 * t224) * t169 + (t19 * t82 + t28 * t221) * t168 + t79 * t1 - t4 * t224) * t227; (-g(2) + t168) * m(4) + t205 * t169 + (t22 * t230 + t23 * t229 + t24 * t228 + (t33 * t231 + t32 * t232 + t31 * t233) * t190) * t167 + (-t151 * t90 - t152 * t91 - t153 * t92 + (t151 ^ 2 + t152 ^ 2 + t153 ^ 2) * t168) * t187 + ((-t33 * t222 + t24 * t81) * t169 + (t33 * t219 + t24 * t84) * t168 + t84 * t3 + t6 * t219) * t225 + ((-t32 * t223 + t23 * t80) * t169 + (t32 * t220 + t23 * t83) * t168 + t83 * t2 + t5 * t220) * t226 + ((t22 * t79 - t31 * t224) * t169 + (t22 * t82 + t31 * t221) * t168 + t82 * t1 + t4 * t221) * t227; t1 * t230 + t2 * t229 + t3 * t228 - g(3) * m(4) + (t6 * t231 + t5 * t232 + t4 * t233) * t190 + (t34 * t230 + t35 * t229 + t36 * t228 + m(4) + (t42 * t231 + t41 * t232 + t40 * t233) * t190) * t167 + ((-t42 * t222 + t36 * t81) * t169 + (t42 * t219 + t36 * t84) * t168) * t225 + ((-t41 * t223 + t35 * t80) * t169 + (t41 * t220 + t35 * t83) * t168) * t226 + ((-t40 * t224 + t34 * t79) * t169 + (t40 * t221 + t34 * t82) * t168) * t227;];
tauX  = t61;
