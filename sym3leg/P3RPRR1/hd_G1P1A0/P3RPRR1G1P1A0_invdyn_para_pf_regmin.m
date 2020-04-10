% Calculate minimal parameter regressor of inverse dynamics forces for
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tauX_reg [3x8]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:23
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3RPRR1G1P1A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G1P1A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:23:17
% EndTime: 2020-03-09 21:23:18
% DurationCPUTime: 1.49s
% Computational Cost: add. (17900->225), mult. (12290->384), div. (1908->11), fcn. (7980->44), ass. (0->190)
t173 = 1 / pkin(3);
t229 = 2 * t173;
t148 = pkin(7) + qJ(3,3);
t131 = sin(t148);
t228 = pkin(1) * t131;
t149 = pkin(7) + qJ(3,2);
t132 = sin(t149);
t227 = pkin(1) * t132;
t150 = pkin(7) + qJ(3,1);
t133 = sin(t150);
t226 = pkin(1) * t133;
t134 = cos(t148);
t225 = pkin(1) * t134;
t135 = cos(t149);
t224 = pkin(1) * t135;
t136 = cos(t150);
t223 = pkin(1) * t136;
t222 = pkin(1) * sin(pkin(7));
t152 = cos(pkin(7));
t221 = pkin(2) * t152;
t170 = xDP(2);
t171 = xDP(1);
t154 = legFrame(2,3);
t138 = t154 + qJ(1,2);
t122 = pkin(7) + t138;
t119 = qJ(3,2) + t122;
t107 = sin(t119);
t59 = -pkin(1) * sin(t138) - pkin(2) * sin(t122) - pkin(3) * t107;
t110 = cos(t119);
t62 = -pkin(1) * cos(t138) - pkin(2) * cos(t122) - pkin(3) * t110;
t40 = t59 * t170 + t62 * t171;
t160 = sin(qJ(3,2));
t95 = t160 * pkin(2) + t227;
t90 = 0.1e1 / t95;
t220 = t40 * t90;
t153 = legFrame(3,3);
t137 = t153 + qJ(1,3);
t121 = pkin(7) + t137;
t118 = qJ(3,3) + t121;
t106 = sin(t118);
t58 = -pkin(1) * sin(t137) - pkin(2) * sin(t121) - pkin(3) * t106;
t109 = cos(t118);
t61 = -pkin(1) * cos(t137) - pkin(2) * cos(t121) - pkin(3) * t109;
t41 = t58 * t170 + t61 * t171;
t158 = sin(qJ(3,3));
t94 = t158 * pkin(2) + t228;
t88 = 0.1e1 / t94;
t219 = t41 * t88;
t155 = legFrame(1,3);
t139 = t155 + qJ(1,1);
t123 = pkin(7) + t139;
t120 = qJ(3,1) + t123;
t108 = sin(t120);
t60 = -pkin(1) * sin(t139) - pkin(2) * sin(t123) - pkin(3) * t108;
t111 = cos(t120);
t63 = -pkin(1) * cos(t139) - pkin(2) * cos(t123) - pkin(3) * t111;
t42 = t60 * t170 + t63 * t171;
t162 = sin(qJ(3,1));
t96 = t162 * pkin(2) + t226;
t92 = 0.1e1 / t96;
t218 = t42 * t92;
t217 = t58 * t88;
t216 = t59 * t90;
t215 = t60 * t92;
t214 = t61 * t88;
t213 = t62 * t90;
t212 = t63 * t92;
t67 = t106 * t170 + t109 * t171;
t64 = t67 ^ 2;
t89 = 0.1e1 / t94 ^ 2;
t211 = t64 * t89;
t68 = t107 * t170 + t110 * t171;
t65 = t68 ^ 2;
t91 = 0.1e1 / t95 ^ 2;
t210 = t65 * t91;
t69 = t108 * t170 + t111 * t171;
t66 = t69 ^ 2;
t93 = 0.1e1 / t96 ^ 2;
t209 = t66 * t93;
t49 = t67 * t88;
t208 = t67 * t89;
t50 = t68 * t90;
t207 = t68 * t91;
t51 = t69 * t92;
t206 = t69 * t93;
t205 = t106 * t88;
t204 = t107 * t90;
t203 = t108 * t92;
t202 = t109 * t88;
t201 = t110 * t90;
t200 = t111 * t92;
t164 = cos(qJ(3,3));
t199 = t164 * t89;
t166 = cos(qJ(3,2));
t198 = t166 * t91;
t168 = cos(qJ(3,1));
t197 = t168 * t93;
t196 = g(1) * t109 + g(2) * t106;
t195 = g(1) * t110 + g(2) * t107;
t194 = g(1) * t111 + g(2) * t108;
t156 = xDDP(2);
t193 = t156 * t173;
t157 = xDDP(1);
t192 = t157 * t173;
t191 = 0.2e1 * pkin(2) * pkin(3);
t190 = 0.2e1 * pkin(1);
t185 = t173 * t219;
t38 = t49 + t185;
t189 = (-pkin(3) * t38 + (-pkin(2) * t164 - t225) * t49) * t208;
t184 = t173 * t218;
t39 = t51 + t184;
t188 = (-t39 * pkin(3) + (-pkin(2) * t168 - t223) * t51) * t206;
t186 = t173 * t220;
t37 = t50 + t186;
t187 = (-pkin(3) * t37 + (-pkin(2) * t166 - t224) * t50) * t207;
t183 = g(1) * t106 - g(2) * t109;
t182 = g(1) * t107 - g(2) * t110;
t181 = g(1) * t108 - g(2) * t111;
t34 = t49 + t185 / 0.2e1;
t180 = -0.2e1 * t34 * t185;
t35 = t50 + t186 / 0.2e1;
t179 = -0.2e1 * t35 * t186;
t36 = t51 + t184 / 0.2e1;
t178 = -0.2e1 * t36 * t184;
t32 = t38 * t89 * t41;
t70 = t156 * t205;
t73 = t157 * t202;
t22 = t32 + t70 + t73 - t189;
t33 = t39 * t93 * t42;
t72 = t156 * t203;
t75 = t157 * t200;
t24 = t33 + t72 + t75 - t188;
t31 = t37 * t91 * t40;
t71 = t156 * t204;
t74 = t157 * t201;
t23 = t31 + t71 + t74 - t187;
t130 = t152 * pkin(1) + pkin(2);
t147 = pkin(1) ^ 2 + pkin(2) ^ 2;
t172 = pkin(3) ^ 2;
t177 = t193 * t217 + t192 * t214 + (-t38 * (t164 * t130 - t158 * t222 + pkin(3)) / (t158 * t130 + t164 * t222) * t219 - (t34 * t164 * t191 + t147 * t49 + t38 * t172 + (pkin(3) * t134 * t34 + t49 * t221) * t190) * t208) * t173;
t176 = t193 * t216 + t192 * t213 + (-t37 * (t166 * t130 - t160 * t222 + pkin(3)) / (t160 * t130 + t166 * t222) * t220 - (t35 * t166 * t191 + t147 * t50 + t37 * t172 + (pkin(3) * t135 * t35 + t50 * t221) * t190) * t207) * t173;
t175 = t193 * t215 + t192 * t212 + (-t39 * (t168 * t130 - t162 * t222 + pkin(3)) / (t162 * t130 + t168 * t222) * t218 - (t36 * t168 * t191 + t147 * t51 + t39 * t172 + (pkin(3) * t136 * t36 + t51 * t221) * t190) * t206) * t173;
t174 = 0.1e1 / pkin(3) ^ 2;
t169 = cos(qJ(1,1));
t167 = cos(qJ(1,2));
t165 = cos(qJ(1,3));
t163 = sin(qJ(1,1));
t161 = sin(qJ(1,2));
t159 = sin(qJ(1,3));
t146 = (xDDP(3) - g(3));
t145 = cos(t155);
t144 = cos(t154);
t143 = cos(t153);
t142 = sin(t155);
t141 = sin(t154);
t140 = sin(t153);
t87 = t145 * g(1) + t142 * g(2);
t86 = t144 * g(1) + t141 * g(2);
t85 = t143 * g(1) + t140 * g(2);
t84 = t142 * g(1) - t145 * g(2);
t83 = t141 * g(1) - t144 * g(2);
t82 = t140 * g(1) - t143 * g(2);
t57 = -t84 * t163 + t87 * t169;
t56 = -t83 * t161 + t86 * t167;
t55 = -t82 * t159 + t85 * t165;
t54 = t87 * t163 + t84 * t169;
t53 = t86 * t161 + t83 * t167;
t52 = t85 * t159 + t82 * t165;
t21 = (g(1) * t163 - g(2) * t169) * t145 + (g(1) * t169 + g(2) * t163) * t142 + pkin(1) * t24;
t20 = (g(1) * t161 - g(2) * t167) * t144 + (g(1) * t167 + g(2) * t161) * t141 + pkin(1) * t23;
t19 = (g(1) * t159 - g(2) * t165) * t143 + (t165 * g(1) + g(2) * t159) * t140 + pkin(1) * t22;
t18 = (-t162 * t24 + t66 * t197) * pkin(2) + (-t133 * t24 + t136 * t209) * pkin(1) + t194;
t17 = (-t160 * t23 + t65 * t198) * pkin(2) + (-t132 * t23 + t135 * t210) * pkin(1) + t195;
t16 = (t158 * t211 + t164 * t22) * pkin(2) + (t131 * t211 + t134 * t22) * pkin(1) + t183;
t15 = (t160 * t210 + t166 * t23) * pkin(2) + (t132 * t210 + t135 * t23) * pkin(1) + t182;
t14 = (-t158 * t22 + t64 * t199) * pkin(2) + (-t131 * t22 + t134 * t211) * pkin(1) + t196;
t13 = (t162 * t209 + t168 * t24) * pkin(2) + (t133 * t209 + t136 * t24) * pkin(1) + t181;
t12 = t175 + t24;
t11 = t176 + t23;
t10 = t177 + t22;
t9 = t175 + 0.2e1 * t33 + 0.2e1 * t72 + 0.2e1 * t75 - 0.2e1 * t188;
t8 = t176 + 0.2e1 * t31 + 0.2e1 * t71 + 0.2e1 * t74 - 0.2e1 * t187;
t7 = t177 + 0.2e1 * t32 + 0.2e1 * t70 + 0.2e1 * t73 - 0.2e1 * t189;
t6 = t178 * t223 - t9 * t226 - pkin(2) * (t162 * t9 + (t174 * t42 + t69 * t229) * t42 * t197) + t194;
t5 = t179 * t224 - t8 * t227 - pkin(2) * (t160 * t8 + (t174 * t40 + t68 * t229) * t40 * t198) + t195;
t4 = t180 * t225 - t7 * t228 - pkin(2) * (t158 * t7 + (t174 * t41 + t67 * t229) * t41 * t199) + t196;
t3 = pkin(2) * (t162 * t178 + t9 * t168) + (t133 * t178 + t9 * t136) * pkin(1) + t181;
t2 = pkin(2) * (t160 * t179 + t8 * t166) + (t132 * t179 + t8 * t135) * pkin(1) + t182;
t1 = pkin(2) * (t158 * t180 + t7 * t164) + (t131 * t180 + t7 * t134) * pkin(1) + t183;
t25 = [t24 * t200 + t23 * t201 + t22 * t202, t54 * t200 + t53 * t201 + t52 * t202, t57 * t200 + t56 * t201 + t55 * t202, (t19 * t202 + t20 * t201 + t21 * t200) * pkin(1), t10 * t202 + t11 * t201 + t12 * t200 + (t10 * t214 + t11 * t213 + t12 * t212) * t173, t1 * t202 + t2 * t201 + t3 * t200 + (t13 * t212 + t15 * t213 + t16 * t214) * t173, t4 * t202 + t5 * t201 + t6 * t200 + (t14 * t214 + t17 * t213 + t18 * t212) * t173, t157 - g(1); t24 * t203 + t23 * t204 + t22 * t205, t54 * t203 + t53 * t204 + t52 * t205, t57 * t203 + t56 * t204 + t55 * t205, (t19 * t205 + t20 * t204 + t21 * t203) * pkin(1), t10 * t205 + t11 * t204 + t12 * t203 + (t10 * t217 + t11 * t216 + t12 * t215) * t173, t1 * t205 + t2 * t204 + t3 * t203 + (t13 * t215 + t15 * t216 + t16 * t217) * t173, t4 * t205 + t5 * t204 + t6 * t203 + (t14 * t217 + t17 * t216 + t18 * t215) * t173, t156 - g(2); 0, 0, 0, 3 * t146, 0, 0, 0, t146;];
tauX_reg  = t25;
