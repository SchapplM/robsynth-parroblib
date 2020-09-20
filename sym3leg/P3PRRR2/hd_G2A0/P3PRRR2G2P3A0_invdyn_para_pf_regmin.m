% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PRRR2G2P3A0
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
%   pkin=[a3,a4]';
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
% Datum: 2020-03-09 21:21
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3PRRR2G2P3A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G2P3A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR2G2P3A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRR2G2P3A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G2P3A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR2G2P3A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G2P3A0_invdyn_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G2P3A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G2P3A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:21:37
% EndTime: 2020-03-09 21:21:39
% DurationCPUTime: 1.87s
% Computational Cost: add. (10699->247), mult. (11604->452), div. (2781->10), fcn. (9129->36), ass. (0->241)
t176 = 0.1e1 / pkin(2);
t166 = cos(qJ(3,3));
t268 = t166 * pkin(1);
t285 = (pkin(2) + t268) * t176;
t168 = cos(qJ(3,2));
t267 = t168 * pkin(1);
t284 = (pkin(2) + t267) * t176;
t170 = cos(qJ(3,1));
t266 = t170 * pkin(1);
t283 = (pkin(2) + t266) * t176;
t172 = xDP(3);
t282 = -2 * t172;
t281 = -g(1) / 0.2e1;
t280 = g(2) / 0.2e1;
t151 = qJ(2,3) + qJ(3,3);
t154 = legFrame(3,2);
t124 = t154 + t151;
t109 = sin(t124);
t279 = t109 / 0.2e1;
t152 = qJ(2,2) + qJ(3,2);
t155 = legFrame(2,2);
t126 = t155 + t152;
t111 = sin(t126);
t278 = t111 / 0.2e1;
t153 = qJ(2,1) + qJ(3,1);
t156 = legFrame(1,2);
t128 = t156 + t153;
t113 = sin(t128);
t277 = t113 / 0.2e1;
t125 = -t154 + t151;
t116 = cos(t125);
t276 = t116 / 0.2e1;
t127 = -t155 + t152;
t118 = cos(t127);
t275 = t118 / 0.2e1;
t129 = -t156 + t153;
t120 = cos(t129);
t274 = t120 / 0.2e1;
t178 = 0.1e1 / pkin(1);
t273 = t178 / 0.4e1;
t272 = pkin(2) * t166;
t271 = pkin(2) * t168;
t270 = pkin(2) * t170;
t269 = pkin(2) * t172;
t158 = xDDP(2);
t265 = t158 - g(2);
t159 = xDDP(1);
t264 = t159 - g(1);
t175 = pkin(2) ^ 2;
t160 = sin(qJ(3,3));
t146 = 0.1e1 / t160 ^ 2;
t179 = 0.1e1 / pkin(1) ^ 2;
t235 = t146 * t179;
t110 = sin(t125);
t115 = cos(t124);
t136 = sin(t151);
t173 = xDP(2);
t174 = xDP(1);
t49 = t136 * t282 + (t109 - t110) * t174 + (t115 + t116) * t173;
t191 = -t49 * t235 / 0.2e1;
t145 = 0.1e1 / t160;
t257 = t49 * t145;
t200 = t257 / 0.2e1;
t157 = xDDP(3);
t226 = t176 * t178;
t204 = t157 * t226;
t161 = sin(qJ(2,3));
t260 = t145 * (t161 * pkin(1) + pkin(2) * t136);
t133 = pkin(1) + t272;
t167 = cos(qJ(2,3));
t139 = sin(t154);
t142 = cos(t154);
t73 = t139 * t174 + t142 * t173;
t40 = (-t73 * t133 + t160 * t269) * t167 + t161 * (pkin(2) * t73 * t160 + t133 * t172);
t207 = t40 * t226;
t194 = t145 * t207;
t43 = t178 * t200;
t34 = t43 + t194 / 0.2e1;
t37 = t43 + t194;
t263 = (t37 * t175 + (0.2e1 * t34 * t272 + t200) * pkin(1)) * t176 * t191 + t204 * t260;
t162 = sin(qJ(3,2));
t148 = 0.1e1 / t162 ^ 2;
t233 = t148 * t179;
t112 = sin(t127);
t117 = cos(t126);
t137 = sin(t152);
t50 = t137 * t282 + (t111 - t112) * t174 + (t117 + t118) * t173;
t190 = -t50 * t233 / 0.2e1;
t147 = 0.1e1 / t162;
t256 = t50 * t147;
t199 = t256 / 0.2e1;
t163 = sin(qJ(2,2));
t259 = t147 * (t163 * pkin(1) + pkin(2) * t137);
t134 = pkin(1) + t271;
t169 = cos(qJ(2,2));
t140 = sin(t155);
t143 = cos(t155);
t74 = t140 * t174 + t143 * t173;
t41 = (-t74 * t134 + t162 * t269) * t169 + t163 * (pkin(2) * t74 * t162 + t134 * t172);
t206 = t41 * t226;
t193 = t147 * t206;
t44 = t178 * t199;
t35 = t44 + t193 / 0.2e1;
t38 = t44 + t193;
t262 = (t38 * t175 + (0.2e1 * t35 * t271 + t199) * pkin(1)) * t176 * t190 + t204 * t259;
t164 = sin(qJ(3,1));
t150 = 0.1e1 / t164 ^ 2;
t231 = t150 * t179;
t114 = sin(t129);
t119 = cos(t128);
t138 = sin(t153);
t51 = t138 * t282 + (t113 - t114) * t174 + (t119 + t120) * t173;
t189 = -t51 * t231 / 0.2e1;
t149 = 0.1e1 / t164;
t255 = t51 * t149;
t198 = t255 / 0.2e1;
t165 = sin(qJ(2,1));
t258 = t149 * (t165 * pkin(1) + pkin(2) * t138);
t135 = pkin(1) + t270;
t171 = cos(qJ(2,1));
t141 = sin(t156);
t144 = cos(t156);
t75 = t141 * t174 + t144 * t173;
t42 = (-t75 * t135 + t164 * t269) * t171 + t165 * (pkin(2) * t75 * t164 + t135 * t172);
t205 = t42 * t226;
t192 = t149 * t205;
t45 = t178 * t198;
t36 = t45 + t192 / 0.2e1;
t39 = t45 + t192;
t261 = (t39 * t175 + (0.2e1 * t36 * t270 + t198) * pkin(1)) * t176 * t189 + t204 * t258;
t254 = t136 * t145;
t253 = t136 * t157;
t252 = t137 * t147;
t251 = t137 * t157;
t250 = t138 * t149;
t249 = t138 * t157;
t248 = t139 * t145;
t247 = t139 * t159;
t246 = t140 * t147;
t245 = t140 * t159;
t244 = t141 * t149;
t243 = t141 * t159;
t242 = t142 * t145;
t241 = t142 * t158;
t240 = t143 * t147;
t239 = t143 * t158;
t238 = t144 * t149;
t237 = t144 * t158;
t236 = t145 * t178;
t234 = t147 * t178;
t232 = t149 * t178;
t230 = t157 * t178;
t229 = t161 * t160;
t228 = t163 * t162;
t227 = t165 * t164;
t70 = -pkin(2) * t229 + t133 * t167;
t225 = t70 * t248;
t76 = t167 * t166 - t229;
t224 = t76 * t248;
t71 = -pkin(2) * t228 + t134 * t169;
t223 = t71 * t246;
t77 = t169 * t168 - t228;
t222 = t77 * t246;
t72 = -pkin(2) * t227 + t135 * t171;
t221 = t72 * t244;
t78 = t171 * t170 - t227;
t220 = t78 * t244;
t219 = t70 * t242;
t218 = t76 * t242;
t217 = t71 * t240;
t216 = t77 * t240;
t215 = t72 * t238;
t214 = t78 * t238;
t213 = t76 * t236;
t212 = t40 * t235;
t211 = t77 * t234;
t210 = t41 * t233;
t209 = t78 * t232;
t208 = t42 * t231;
t203 = t49 ^ 2 * t273;
t202 = t50 ^ 2 * t273;
t201 = t51 ^ 2 * t273;
t28 = t37 * t212;
t31 = -t166 * t257 / 0.2e1 - pkin(2) * t37;
t55 = t213 * t247;
t58 = t213 * t241;
t197 = t31 * t191 + t28 + t55 + t58;
t29 = t38 * t210;
t32 = -t168 * t256 / 0.2e1 - pkin(2) * t38;
t56 = t211 * t245;
t59 = t211 * t239;
t196 = t32 * t190 + t29 + t56 + t59;
t30 = t39 * t208;
t33 = -t170 * t255 / 0.2e1 - t39 * pkin(2);
t57 = t209 * t243;
t60 = t209 * t237;
t195 = t33 * t189 + t30 + t57 + t60;
t188 = g(1) * t276 + g(2) * t279 + t110 * t280 + t115 * t281 + g(3) * cos(t151);
t187 = g(1) * t275 + g(2) * t278 + t112 * t280 + t117 * t281 + g(3) * cos(t152);
t186 = g(1) * t274 + g(2) * t277 + t114 * t280 + t119 * t281 + g(3) * cos(t153);
t185 = g(1) * t279 + g(2) * t276 - g(3) * t136 + t110 * t281 + t115 * t280;
t184 = g(1) * t278 + g(2) * t275 - g(3) * t137 + t112 * t281 + t117 * t280;
t183 = g(1) * t277 + g(2) * t274 - g(3) * t138 + t114 * t281 + t119 * t280;
t182 = (-t241 - t247) * t70 * t176;
t181 = (-t239 - t245) * t71 * t176;
t180 = (-t237 - t243) * t72 * t176;
t177 = 0.1e1 / pkin(2) ^ 2;
t81 = t141 * g(1) + t144 * g(2);
t80 = t140 * g(1) + t143 * g(2);
t79 = t139 * g(1) + t142 * g(2);
t66 = g(3) * t171 + t81 * t165;
t65 = g(3) * t169 + t80 * t163;
t64 = g(3) * t167 + t79 * t161;
t63 = -g(3) * t165 + t81 * t171;
t62 = -g(3) * t163 + t80 * t169;
t61 = -g(3) * t161 + t79 * t167;
t54 = t265 * t141 - t264 * t144;
t53 = t265 * t140 - t264 * t143;
t52 = t265 * t139 - t264 * t142;
t21 = -t230 * t250 + t195;
t20 = -t230 * t252 + t196;
t19 = -t230 * t254 + t197;
t18 = t149 * t201 + t21 * t266 + t186;
t17 = t147 * t202 + t20 * t267 + t187;
t16 = t145 * t203 + t19 * t268 + t188;
t15 = -t164 * t21 * pkin(1) + t170 * t150 * t201 + t183;
t14 = -t162 * t20 * pkin(1) + t168 * t148 * t202 + t184;
t13 = -t160 * t19 * pkin(1) + t166 * t146 * t203 + t185;
t12 = 0.2e1 * t30 + 0.2e1 * t57 + 0.2e1 * t60 + (-t39 * t42 * t283 - t33 * t51) * t231 + (t180 - 0.2e1 * t249) * t232 + t261;
t11 = 0.2e1 * t29 + 0.2e1 * t56 + 0.2e1 * t59 + (-t38 * t41 * t284 - t32 * t50) * t233 + (t181 - 0.2e1 * t251) * t234 + t262;
t10 = 0.2e1 * t28 + 0.2e1 * t55 + 0.2e1 * t58 + (-t37 * t40 * t285 - t31 * t49) * t235 + (t182 - 0.2e1 * t253) * t236 + t263;
t9 = -t30 * t283 + (t180 - t249) * t232 + t195 + t261;
t8 = -t29 * t284 + (t181 - t251) * t234 + t196 + t262;
t7 = -t28 * t285 + (t182 - t253) * t236 + t197 + t263;
t6 = pkin(1) * (t12 * t170 - 0.2e1 * t36 * t205) + t186;
t5 = pkin(1) * (t11 * t168 - 0.2e1 * t35 * t206) + t187;
t4 = pkin(1) * (t10 * t166 - 0.2e1 * t34 * t207) + t188;
t3 = -pkin(1) * (t164 * t12 + (t176 * t51 + t177 * t42) * t170 * t208) + t183;
t2 = -pkin(1) * (t162 * t11 + (t176 * t50 + t177 * t41) * t168 * t210) + t184;
t1 = -pkin(1) * (t160 * t10 + (t176 * t49 + t177 * t40) * t166 * t212) + t185;
t22 = [-t142 * t52 - t143 * t53 - t144 * t54, (t19 * t224 + t20 * t222 + t21 * t220) * t178, (t66 * t220 + t65 * t222 + t64 * t224) * t178, (t63 * t220 + t62 * t222 + t61 * t224) * t178, (t7 * t224 + t8 * t222 + t9 * t220 + (-t9 * t221 - t8 * t223 - t7 * t225) * t176) * t178, (t4 * t224 + t5 * t222 + t6 * t220 + (-t16 * t225 - t17 * t223 - t18 * t221) * t176) * t178, (t1 * t224 + t2 * t222 + t3 * t220 + (-t13 * t225 - t14 * t223 - t15 * t221) * t176) * t178, t264; t139 * t52 + t140 * t53 + t141 * t54, (t19 * t218 + t20 * t216 + t21 * t214) * t178, (t66 * t214 + t65 * t216 + t64 * t218) * t178, (t63 * t214 + t62 * t216 + t61 * t218) * t178, (t7 * t218 + t8 * t216 + t9 * t214 + (-t9 * t215 - t8 * t217 - t7 * t219) * t176) * t178, (t4 * t218 + t5 * t216 + t6 * t214 + (-t16 * t219 - t17 * t217 - t18 * t215) * t176) * t178, (t1 * t218 + t2 * t216 + t3 * t214 + (-t13 * t219 - t14 * t217 - t15 * t215) * t176) * t178, t265; 0, (-t19 * t254 - t20 * t252 - t21 * t250) * t178, (-t66 * t250 - t65 * t252 - t64 * t254) * t178, (-t63 * t250 - t62 * t252 - t61 * t254) * t178, (-t7 * t254 - t8 * t252 - t9 * t250 + (t9 * t258 + t8 * t259 + t7 * t260) * t176) * t178, (-t4 * t254 - t5 * t252 - t6 * t250 + (t16 * t260 + t17 * t259 + t18 * t258) * t176) * t178, (-t1 * t254 - t2 * t252 - t3 * t250 + (t13 * t260 + t14 * t259 + t15 * t258) * t176) * t178, t157 - g(3);];
tauX_reg  = t22;
