% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PRRRR1G2P3A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tauX_reg [3x12]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:16
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3PRRRR1G2P3A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G2P3A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR1G2P3A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR1G2P3A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G2P3A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR1G2P3A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G2P3A0_invdyn_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G2P3A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G2P3A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:16:31
% EndTime: 2020-03-09 21:16:33
% DurationCPUTime: 2.36s
% Computational Cost: add. (2539->238), mult. (6963->520), div. (3429->17), fcn. (7395->18), ass. (0->231)
t127 = cos(qJ(3,3));
t103 = 0.1e1 / t127;
t122 = sin(qJ(2,3));
t95 = 0.1e1 / t122;
t250 = t103 * t95;
t128 = cos(qJ(2,3));
t118 = xDDP(3);
t119 = xDDP(2);
t120 = xDDP(1);
t146 = t127 ^ 2;
t104 = 0.1e1 / t146;
t105 = t103 * t104;
t135 = 0.1e1 / pkin(2);
t217 = t105 * t135;
t133 = xDP(3);
t121 = sin(qJ(3,3));
t214 = t121 * t128;
t134 = xDP(1);
t256 = xDP(2);
t115 = legFrame(3,2);
t88 = sin(t115);
t91 = cos(t115);
t70 = t134 * t91 - t256 * t88;
t58 = -t127 * t133 + t70 * t214;
t55 = t58 ^ 2;
t96 = 0.1e1 / t122 ^ 2;
t262 = t55 * t96;
t213 = t122 * t127;
t222 = t91 * t121;
t61 = t88 * t213 - t222;
t225 = t88 * t121;
t64 = t91 * t213 + t225;
t67 = t70 ^ 2;
t37 = -t88 * g(1) - t91 * g(2) + (t128 * t118 + (t119 * t64 + t120 * t61) * t103 + (t67 + t262) * t217) * t95;
t233 = t37 * t128;
t28 = g(3) * t122 + t233;
t265 = t28 * t250;
t129 = cos(qJ(3,2));
t107 = 0.1e1 / t129;
t124 = sin(qJ(2,2));
t98 = 0.1e1 / t124;
t247 = t107 * t98;
t130 = cos(qJ(2,2));
t149 = t129 ^ 2;
t108 = 0.1e1 / t149;
t109 = t107 * t108;
t216 = t109 * t135;
t123 = sin(qJ(3,2));
t212 = t123 * t130;
t116 = legFrame(2,2);
t89 = sin(t116);
t92 = cos(t116);
t72 = t134 * t92 - t256 * t89;
t59 = -t129 * t133 + t72 * t212;
t56 = t59 ^ 2;
t99 = 0.1e1 / t124 ^ 2;
t261 = t56 * t99;
t211 = t124 * t129;
t221 = t92 * t123;
t62 = t89 * t211 - t221;
t224 = t89 * t123;
t65 = t92 * t211 + t224;
t68 = t72 ^ 2;
t38 = -t89 * g(1) - t92 * g(2) + (t130 * t118 + (t119 * t65 + t120 * t62) * t107 + (t68 + t261) * t216) * t98;
t232 = t38 * t130;
t29 = g(3) * t124 + t232;
t264 = t29 * t247;
t126 = sin(qJ(2,1));
t101 = 0.1e1 / t126;
t131 = cos(qJ(3,1));
t111 = 0.1e1 / t131;
t219 = t101 * t111;
t132 = cos(qJ(2,1));
t152 = t131 ^ 2;
t112 = 0.1e1 / t152;
t113 = t111 * t112;
t215 = t113 * t135;
t102 = 0.1e1 / t126 ^ 2;
t125 = sin(qJ(3,1));
t210 = t125 * t132;
t117 = legFrame(1,2);
t90 = sin(t117);
t93 = cos(t117);
t74 = t134 * t93 - t256 * t90;
t60 = -t131 * t133 + t74 * t210;
t57 = t60 ^ 2;
t260 = t102 * t57;
t209 = t126 * t131;
t220 = t93 * t125;
t63 = t90 * t209 - t220;
t223 = t90 * t125;
t66 = t93 * t209 + t223;
t69 = t74 ^ 2;
t39 = -t90 * g(1) - t93 * g(2) + (t132 * t118 + (t119 * t66 + t120 * t63) * t111 + (t69 + t260) * t215) * t101;
t30 = g(3) * t126 + t39 * t132;
t263 = t30 * t219;
t136 = 0.1e1 / pkin(2) ^ 2;
t259 = -0.2e1 * t104;
t258 = -0.2e1 * t108;
t257 = -0.2e1 * t112;
t255 = t70 * t95;
t254 = t72 * t98;
t253 = t101 * t74;
t252 = t103 * t88;
t251 = t103 * t91;
t249 = t107 * t89;
t248 = t107 * t92;
t246 = t111 * t90;
t245 = t111 * t93;
t172 = t119 * t88 - t120 * t91;
t194 = t95 * t214;
t181 = t104 * t194;
t243 = t128 * t95;
t19 = (-t172 * t181 + (-t95 * t118 + (-(t121 * t122 * t70 + t58 * t243) * t96 * t58 + (-t121 * t58 - t128 * t70) * t255) * t217) * t103) * t135;
t244 = t127 * t19;
t171 = t119 * t89 - t120 * t92;
t192 = t98 * t212;
t180 = t108 * t192;
t241 = t130 * t98;
t20 = (-t171 * t180 + (-t98 * t118 + (-(t123 * t124 * t72 + t59 * t241) * t99 * t59 + (-t123 * t59 - t130 * t72) * t254) * t216) * t107) * t135;
t242 = t129 * t20;
t170 = t119 * t90 - t120 * t93;
t188 = t101 * t210;
t179 = t112 * t188;
t218 = t101 * t132;
t21 = (-t170 * t179 + (-t101 * t118 + (-(t125 * t126 * t74 + t60 * t218) * t102 * t60 + (-t125 * t60 - t132 * t74) * t253) * t215) * t111) * t135;
t240 = t131 * t21;
t239 = t136 * t55;
t238 = t136 * t56;
t237 = t136 * t57;
t236 = t136 * t67;
t235 = t136 * t68;
t234 = t136 * t69;
t193 = t121 * t236;
t52 = t103 * t135 * t172 + t105 * t193;
t231 = t52 * t121;
t230 = t52 * t127;
t191 = t123 * t235;
t53 = t107 * t135 * t171 + t109 * t191;
t229 = t53 * t123;
t228 = t53 * t129;
t190 = t125 * t234;
t54 = t111 * t135 * t170 + t113 * t190;
t227 = t54 * t125;
t226 = t54 * t131;
t208 = t61 * t250;
t207 = t64 * t250;
t106 = 0.1e1 / t146 ^ 2;
t206 = t106 * t262;
t205 = t62 * t247;
t204 = t65 * t247;
t110 = 0.1e1 / t149 ^ 2;
t203 = t110 * t261;
t202 = t63 * t219;
t201 = t66 * t219;
t114 = 0.1e1 / t152 ^ 2;
t200 = t114 * t260;
t199 = t104 * t243;
t198 = t106 * t239;
t197 = t108 * t241;
t196 = t110 * t238;
t195 = t114 * t237;
t189 = t112 * t218;
t187 = t136 * t58 * t255;
t186 = t136 * t59 * t254;
t185 = t136 * t60 * t253;
t184 = t121 * t206;
t183 = t123 * t203;
t182 = t125 * t200;
t178 = t88 * t181;
t177 = t91 * t181;
t176 = t89 * t180;
t175 = t92 * t180;
t174 = t90 * t179;
t173 = t93 * t179;
t169 = t187 * t259;
t168 = t186 * t258;
t31 = g(3) * t128 - t37 * t122;
t76 = t91 * g(1) - t88 * g(2);
t167 = t31 * t121 + t76 * t127 - t194 * t28;
t33 = g(3) * t130 - t38 * t124;
t77 = t92 * g(1) - t89 * g(2);
t166 = t33 * t123 + t77 * t129 - t192 * t29;
t165 = t185 * t257;
t35 = g(3) * t132 - t39 * t126;
t78 = t93 * g(1) - t90 * g(2);
t164 = t35 * t125 + t78 * t131 - t188 * t30;
t40 = -t104 * t193 + t230;
t163 = t40 * t181 - t19;
t41 = -t108 * t191 + t228;
t162 = t41 * t180 - t20;
t42 = -t112 * t190 + t226;
t161 = t42 * t179 - t21;
t46 = t103 * t236 + t231;
t160 = -t103 * t19 + t46 * t199;
t47 = t107 * t235 + t229;
t159 = -t107 * t20 + t47 * t197;
t48 = t111 * t234 + t227;
t158 = -t111 * t21 + t48 * t189;
t94 = t121 ^ 2;
t157 = -t28 * t94 * t199 - t103 * (-t76 * t121 + t31 * t127);
t97 = t123 ^ 2;
t156 = -t29 * t97 * t197 - t107 * (-t77 * t123 + t33 * t129);
t100 = t125 ^ 2;
t155 = -t100 * t30 * t189 - t111 * (-t78 * t125 + t35 * t131);
t137 = t135 * t136;
t51 = (t112 * t69 + t200) * t136;
t50 = (t108 * t68 + t203) * t136;
t49 = (t104 * t67 + t206) * t136;
t45 = (t257 + t114) * t102 * t237;
t44 = (t258 + t110) * t99 * t238;
t43 = (t259 + t106) * t96 * t239;
t18 = -t101 * t195 + t132 * t21;
t17 = t130 * t20 - t98 * t196;
t16 = t128 * t19 - t95 * t198;
t15 = -t102 * t132 * t195 - t21 * t126;
t14 = -t99 * t130 * t196 - t20 * t124;
t13 = -t96 * t128 * t198 - t19 * t122;
t12 = t100 * t21 + t125 * t165;
t11 = t123 * t168 + t97 * t20;
t10 = t121 * t169 + t94 * t19;
t9 = t125 * t240 + (-0.2e1 * t111 + t113) * t185;
t8 = t123 * t242 + (-0.2e1 * t107 + t109) * t186;
t7 = t121 * t244 + (-0.2e1 * t103 + t105) * t187;
t6 = (t51 * t125 - t226) * t126 - t132 * (t125 * t21 + t165);
t5 = (t50 * t123 - t228) * t124 - t130 * (t123 * t20 + t168);
t4 = (t49 * t121 - t230) * t122 - t128 * (t121 * t19 + t169);
t3 = (-t51 * t131 - t227) * t126 + (0.2e1 * t125 * t113 * t185 + t240) * t132;
t2 = (-t50 * t129 - t229) * t124 + (0.2e1 * t123 * t109 * t186 + t242) * t130;
t1 = (-t49 * t127 - t231) * t122 + (0.2e1 * t121 * t105 * t187 + t244) * t128;
t22 = [t39 * t202 + t38 * t205 + t37 * t208, (t173 * t21 + t20 * t175 + t19 * t177) * t135, t18 * t202 + t16 * t208 + t17 * t205 + (t173 * t30 + t29 * t175 + t28 * t177) * t135, t15 * t202 + t13 * t208 + t14 * t205 + (t173 * t35 + t33 * t175 + t31 * t177) * t135, (t93 * t182 + t92 * t183 + t91 * t184) * t137 + (t10 * t177 + t11 * t175 + t12 * t173) * t135, (0.2e1 * t173 * t9 + 0.2e1 * t8 * t175 + 0.2e1 * t7 * t177 - t45 * t245 - t44 * t248 - t43 * t251) * t135, (t158 * t220 + t159 * t221 + t160 * t222) * t135, (t161 * t93 + t162 * t92 + t163 * t91) * t135, (-t54 * t245 - t53 * t248 - t52 * t251) * t135, t1 * t208 + t3 * t202 + t2 * t205 + (-t164 * t245 - t166 * t248 - t167 * t251) * t135, t6 * t202 + t4 * t208 + t5 * t205 + (t155 * t93 + t156 * t92 + t157 * t91) * t135, t120 - g(1); t39 * t201 + t38 * t204 + t37 * t207, (-t174 * t21 - t20 * t176 - t19 * t178) * t135, t18 * t201 + t16 * t207 + t17 * t204 + (-t174 * t30 - t29 * t176 - t28 * t178) * t135, t15 * t201 + t13 * t207 + t14 * t204 + (-t174 * t35 - t33 * t176 - t31 * t178) * t135, (-t90 * t182 - t89 * t183 - t88 * t184) * t137 + (-t10 * t178 - t11 * t176 - t12 * t174) * t135, (-0.2e1 * t9 * t174 - 0.2e1 * t8 * t176 - 0.2e1 * t7 * t178 + t45 * t246 + t44 * t249 + t43 * t252) * t135, (-t158 * t223 - t159 * t224 - t160 * t225) * t135, (-t161 * t90 - t162 * t89 - t163 * t88) * t135, (t54 * t246 + t53 * t249 + t52 * t252) * t135, t1 * t207 + t3 * t201 + t2 * t204 + (t164 * t246 + t166 * t249 + t167 * t252) * t135, t6 * t201 + t4 * t207 + t5 * t204 + (-t155 * t90 - t156 * t89 - t157 * t88) * t135, t119 - g(2); t39 * t218 + t98 * t232 + t95 * t233, (-t19 * t250 - t20 * t247 - t21 * t219) * t135, t18 * t218 + t16 * t243 + t17 * t241 + (-t263 - t264 - t265) * t135, t15 * t218 + t13 * t243 + t14 * t241 + (-t35 * t219 - t33 * t247 - t31 * t250) * t135, (-t10 * t250 - t11 * t247 - t12 * t219) * t135, 0.2e1 * (-t9 * t219 - t8 * t247 - t7 * t250) * t135, (-t48 * t219 - t47 * t247 - t46 * t250) * t135, (-t42 * t219 - t41 * t247 - t40 * t250) * t135, 0, t1 * t243 + t3 * t218 + t2 * t241 + (-t101 * t30 - t28 * t95 - t29 * t98) * t135, t6 * t218 + t4 * t243 + t5 * t241 + (t121 * t265 + t123 * t264 + t125 * t263) * t135, t118 - g(3);];
tauX_reg  = t22;
