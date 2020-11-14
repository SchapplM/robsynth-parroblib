% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RRPRR1V1G2A0
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tauX_reg [3x13]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:34
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3RRPRR1V1G2A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_regmin: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR1V1G2A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:33:57
% EndTime: 2020-08-06 19:34:03
% DurationCPUTime: 5.16s
% Computational Cost: add. (8650->343), mult. (13197->667), div. (2346->14), fcn. (11253->18), ass. (0->288)
t136 = pkin(3) + qJ(3,3);
t116 = 0.1e1 / t136 ^ 2;
t151 = cos(qJ(2,3));
t124 = 0.1e1 / t151;
t254 = t116 * t124;
t139 = legFrame(3,2);
t106 = sin(t139);
t109 = cos(t139);
t146 = sin(qJ(1,3));
t152 = cos(qJ(1,3));
t157 = xDP(3);
t158 = xDP(2);
t159 = xDP(1);
t145 = sin(qJ(2,3));
t85 = t106 * t159 + t109 * t158;
t292 = t145 * t85;
t49 = (t152 * t157 + (-t106 * t158 + t109 * t159) * t146) * t151 + t292;
t220 = t49 * t254;
t160 = pkin(1) + pkin(2);
t237 = t146 * t151;
t241 = t136 * t152;
t232 = t151 * t160;
t91 = t146 * t136 + t152 * t232;
t31 = (-t159 * t241 + (t145 * t158 + t159 * t237) * t160) * t109 + (t158 * t241 - (-t145 * t159 + t158 * t237) * t160) * t106 + t91 * t157;
t309 = t109 * g(1);
t312 = t106 * g(2);
t97 = t309 - t312;
t73 = g(3) * t152 + t97 * t146;
t338 = -0.2e1 * t31 * t220 + t73;
t137 = pkin(3) + qJ(3,2);
t119 = 0.1e1 / t137 ^ 2;
t153 = cos(qJ(2,2));
t127 = 0.1e1 / t153;
t251 = t119 * t127;
t140 = legFrame(2,2);
t107 = sin(t140);
t110 = cos(t140);
t148 = sin(qJ(1,2));
t154 = cos(qJ(1,2));
t147 = sin(qJ(2,2));
t86 = t107 * t159 + t110 * t158;
t291 = t147 * t86;
t50 = (t154 * t157 + (-t107 * t158 + t110 * t159) * t148) * t153 + t291;
t215 = t50 * t251;
t235 = t148 * t153;
t240 = t137 * t154;
t231 = t153 * t160;
t92 = t148 * t137 + t154 * t231;
t32 = (-t159 * t240 + (t147 * t158 + t159 * t235) * t160) * t110 + (t158 * t240 - (-t147 * t159 + t158 * t235) * t160) * t107 + t92 * t157;
t308 = t110 * g(1);
t311 = t107 * g(2);
t98 = t308 - t311;
t74 = g(3) * t154 + t98 * t148;
t337 = -0.2e1 * t32 * t215 + t74;
t138 = pkin(3) + qJ(3,1);
t122 = 0.1e1 / t138 ^ 2;
t155 = cos(qJ(2,1));
t130 = 0.1e1 / t155;
t248 = t122 * t130;
t141 = legFrame(1,2);
t108 = sin(t141);
t111 = cos(t141);
t150 = sin(qJ(1,1));
t156 = cos(qJ(1,1));
t149 = sin(qJ(2,1));
t87 = t108 * t159 + t111 * t158;
t290 = t149 * t87;
t51 = (t156 * t157 + (-t108 * t158 + t111 * t159) * t150) * t155 + t290;
t210 = t51 * t248;
t233 = t150 * t155;
t239 = t138 * t156;
t230 = t155 * t160;
t93 = t150 * t138 + t156 * t230;
t33 = (-t159 * t239 + (t149 * t158 + t159 * t233) * t160) * t111 + (t158 * t239 - (-t149 * t159 + t158 * t233) * t160) * t108 + t93 * t157;
t307 = t111 * g(1);
t310 = t108 * g(2);
t99 = t307 - t310;
t75 = g(3) * t156 + t99 * t150;
t336 = -0.2e1 * t33 * t210 + t75;
t115 = 0.1e1 / t136;
t255 = t115 * t152;
t316 = g(3) * t146;
t70 = -t97 * t152 + t316;
t221 = t70 * t255;
t118 = 0.1e1 / t137;
t252 = t118 * t154;
t315 = g(3) * t148;
t71 = -t98 * t154 + t315;
t216 = t71 * t252;
t121 = 0.1e1 / t138;
t249 = t121 * t156;
t314 = g(3) * t150;
t72 = -t99 * t156 + t314;
t211 = t72 * t249;
t172 = t155 ^ 2;
t131 = 0.1e1 / t172;
t132 = t130 * t131;
t134 = 0.1e1 / t160;
t142 = xDDP(3);
t143 = xDDP(2);
t144 = xDDP(1);
t123 = t121 * t122;
t224 = t51 * t123 * t33;
t78 = -t108 * t233 + t149 * t111;
t81 = t149 * t108 + t111 * t233;
t84 = t87 ^ 2;
t21 = -t130 * t224 + (t134 * t84 * t132 + t156 * t142 + (t81 * t144 + t78 * t143 - (-t160 * t51 + t33) * t122 * t51) * t130) * t121;
t285 = t21 * qJ(3,1);
t335 = -t285 + t336;
t170 = t153 ^ 2;
t128 = 0.1e1 / t170;
t129 = t127 * t128;
t120 = t118 * t119;
t225 = t50 * t120 * t32;
t77 = -t107 * t235 + t147 * t110;
t80 = t147 * t107 + t110 * t235;
t83 = t86 ^ 2;
t20 = -t127 * t225 + (t134 * t83 * t129 + t154 * t142 + (t80 * t144 + t77 * t143 - (-t160 * t50 + t32) * t119 * t50) * t127) * t118;
t287 = t20 * qJ(3,2);
t334 = -t287 + t337;
t168 = t151 ^ 2;
t125 = 0.1e1 / t168;
t126 = t124 * t125;
t117 = t115 * t116;
t226 = t49 * t117 * t31;
t76 = -t106 * t237 + t145 * t109;
t79 = t145 * t106 + t109 * t237;
t82 = t85 ^ 2;
t19 = -t124 * t226 + (t134 * t82 * t126 + t152 * t142 + (t79 * t144 + t76 * t143 - (-t160 * t49 + t31) * t116 * t49) * t124) * t115;
t289 = t19 * qJ(3,3);
t333 = -t289 + t338;
t323 = t125 - 0.2e1;
t322 = t128 - 0.2e1;
t321 = t131 - 0.2e1;
t320 = pkin(1) * t134;
t319 = pkin(1) * t151;
t318 = pkin(1) * t153;
t317 = pkin(1) * t155;
t161 = pkin(1) ^ 2;
t246 = t125 * t145;
t298 = t134 * t49;
t184 = pkin(1) * t246 * t298;
t195 = t160 * t226;
t133 = t160 ^ 2;
t22 = (-t49 + t292) * t136 * t124 + (-t133 * t49 + t160 * t31) * t151 * t115;
t263 = qJ(3,3) * t125;
t271 = t91 * t142;
t178 = t146 * t232 - t241;
t238 = t145 * t160;
t61 = t106 * t238 + t178 * t109;
t274 = t61 * t144;
t58 = -t178 * t106 + t109 * t238;
t277 = t58 * t143;
t135 = 0.1e1 / t160 ^ 2;
t295 = t135 * t82;
t205 = t145 * t295;
t46 = t126 * t205 + (t106 * t144 + t109 * t143) * t134 * t124;
t280 = t46 * t145;
t283 = (-t145 * t49 + t85) * t124 ^ 2;
t1 = t168 * t161 * t19 - (pkin(1) * t280 + t333) * qJ(3,3) - 0.2e1 * ((t309 / 0.2e1 - t312 / 0.2e1) * t152 + t263 * t295 / 0.2e1 - t316 / 0.2e1 - t22 * t220 / 0.2e1 - t195 / 0.2e1 + (t274 / 0.2e1 + t277 / 0.2e1 + t271 / 0.2e1 + (t184 + t283 / 0.2e1) * t85) * t115) * t319;
t313 = t1 * t124;
t244 = t128 * t147;
t297 = t134 * t50;
t183 = pkin(1) * t244 * t297;
t193 = t160 * t225;
t23 = (-t50 + t291) * t137 * t127 + (-t133 * t50 + t160 * t32) * t153 * t118;
t264 = qJ(3,2) * t128;
t270 = t92 * t142;
t177 = t148 * t231 - t240;
t236 = t147 * t160;
t62 = t107 * t236 + t177 * t110;
t273 = t62 * t144;
t60 = -t177 * t107 + t110 * t236;
t275 = t60 * t143;
t294 = t135 * t83;
t204 = t147 * t294;
t47 = t129 * t204 + (t107 * t144 + t110 * t143) * t134 * t127;
t279 = t47 * t147;
t282 = (-t147 * t50 + t86) * t127 ^ 2;
t2 = t170 * t161 * t20 - (pkin(1) * t279 + t334) * qJ(3,2) - 0.2e1 * ((t308 / 0.2e1 - t311 / 0.2e1) * t154 + t264 * t294 / 0.2e1 - t315 / 0.2e1 - t23 * t215 / 0.2e1 - t193 / 0.2e1 + (t273 / 0.2e1 + t275 / 0.2e1 + t270 / 0.2e1 + (t183 + t282 / 0.2e1) * t86) * t118) * t318;
t306 = t127 * t2;
t242 = t131 * t149;
t296 = t134 * t51;
t182 = pkin(1) * t242 * t296;
t191 = t160 * t224;
t24 = (-t51 + t290) * t138 * t130 + (-t133 * t51 + t160 * t33) * t155 * t121;
t265 = qJ(3,1) * t131;
t269 = t93 * t142;
t176 = t150 * t230 - t239;
t234 = t149 * t160;
t63 = t108 * t234 + t176 * t111;
t272 = t63 * t144;
t59 = -t176 * t108 + t111 * t234;
t276 = t59 * t143;
t293 = t135 * t84;
t203 = t149 * t293;
t48 = t132 * t203 + (t108 * t144 + t111 * t143) * t134 * t130;
t278 = t48 * t149;
t281 = (-t149 * t51 + t87) * t130 ^ 2;
t3 = t172 * t161 * t21 - (pkin(1) * t278 + t335) * qJ(3,1) - 0.2e1 * ((t307 / 0.2e1 - t310 / 0.2e1) * t156 + t265 * t293 / 0.2e1 - t314 / 0.2e1 - t24 * t210 / 0.2e1 - t191 / 0.2e1 + (t272 / 0.2e1 + t276 / 0.2e1 + t269 / 0.2e1 + (t182 + t281 / 0.2e1) * t87) * t121) * t317;
t305 = t130 * t3;
t304 = t115 * t70;
t43 = t49 ^ 2;
t303 = t116 * t43;
t302 = t118 * t71;
t44 = t50 ^ 2;
t301 = t119 * t44;
t300 = t121 * t72;
t45 = t51 ^ 2;
t299 = t122 * t45;
t288 = t19 * t145;
t286 = t20 * t147;
t284 = t21 * t149;
t94 = t106 * g(1) + t109 * g(2);
t268 = t94 * t151;
t95 = t107 * g(1) + t110 * g(2);
t267 = t95 * t153;
t96 = t108 * g(1) + t111 * g(2);
t266 = t96 * t155;
t262 = t106 * t124;
t261 = t107 * t127;
t260 = t108 * t130;
t259 = t109 * t124;
t258 = t110 * t127;
t257 = t111 * t130;
t256 = t115 * t124;
t253 = t118 * t127;
t250 = t121 * t130;
t247 = t124 * t145;
t245 = t127 * t147;
t243 = t130 * t149;
t223 = t76 * t256;
t222 = t79 * t256;
t219 = t117 * t125 * t43;
t218 = t77 * t253;
t217 = t80 * t253;
t214 = t120 * t128 * t44;
t213 = t78 * t250;
t212 = t81 * t250;
t209 = t123 * t131 * t45;
t208 = t19 * t247;
t207 = t20 * t245;
t206 = t21 * t243;
t199 = t115 * t85 * t298;
t13 = t151 * t288 - t323 * t199;
t202 = 0.2e1 * t13 * t256;
t198 = t118 * t86 * t297;
t14 = t153 * t286 - t322 * t198;
t201 = 0.2e1 * t14 * t253;
t197 = t121 * t87 * t296;
t15 = t155 * t284 - t321 * t197;
t200 = 0.2e1 * t15 * t250;
t190 = t247 * t304;
t189 = t246 * t303;
t188 = t245 * t302;
t187 = t244 * t301;
t186 = t243 * t300;
t185 = t242 * t299;
t37 = t124 * t295 + t280;
t38 = t127 * t294 + t279;
t39 = t130 * t293 + t278;
t175 = t106 * t208 + t107 * t207 + t108 * t206;
t174 = t109 * t208 + t110 * t207 + t111 * t206;
t57 = t96 * t149 + t155 * t75;
t56 = t149 * t75 - t266;
t55 = t95 * t147 + t153 * t74;
t54 = t147 * t74 - t267;
t53 = t94 * t145 + t151 * t73;
t52 = t145 * t73 - t268;
t36 = -t131 * t203 + t48 * t155;
t35 = -t128 * t204 + t47 * t153;
t34 = -t125 * t205 + t46 * t151;
t30 = t321 * t299;
t29 = t322 * t301;
t28 = t323 * t303;
t18 = (0.2e1 * t130 * t197 + t284) * t149;
t17 = (0.2e1 * t127 * t198 + t286) * t147;
t16 = (0.2e1 * t124 * t199 + t288) * t145;
t12 = -t39 * pkin(1) + 0.2e1 * t285 - t336;
t11 = -t38 * pkin(1) + 0.2e1 * t287 - t337;
t10 = -t37 * pkin(1) + 0.2e1 * t289 - t338;
t9 = (pkin(1) * t45 * t248 + t335) * t149 - t266 + t48 * pkin(1);
t8 = (pkin(1) * t44 * t251 + t334) * t147 - t267 + t47 * pkin(1);
t7 = (pkin(1) * t43 * t254 + t333) * t145 - t268 + t46 * pkin(1);
t6 = -t191 - t21 * t317 + (-t130 * t24 * t51 - t45 * t265) * t122 + (t269 + t276 + t272 + (0.2e1 * t182 + t281) * t87) * t121 - t72;
t5 = -t193 - t20 * t318 + (-t127 * t23 * t50 - t44 * t264) * t119 + (t270 + t275 + t273 + (0.2e1 * t183 + t282) * t86) * t118 - t71;
t4 = -t195 - t19 * t319 + (-t124 * t22 * t49 - t43 * t263) * t116 + (t271 + t277 + t274 + (0.2e1 * t184 + t283) * t85) * t115 - t70;
t25 = [t19 * t222 + t20 * t217 + t21 * t212, t72 * t212 + t71 * t217 + t70 * t222, t75 * t212 + t74 * t217 + t73 * t222, t16 * t222 + t17 * t217 + t18 * t212 + (-t106 * t189 - t107 * t187 - t108 * t185) * t134, (t30 * t260 + t29 * t261 + t28 * t262) * t134 + t79 * t202 + t80 * t201 + t81 * t200, t134 * t175 + t212 * t39 + t217 * t38 + t222 * t37, t34 * t222 + t35 * t217 + t36 * t212 + (t106 * t19 + t107 * t20 + t108 * t21) * t134, (t48 * t260 + t47 * t261 + t46 * t262) * t134, t79 * t304 + t80 * t302 + t81 * t300 + (t260 * t56 + t261 * t54 + t262 * t52) * t134, -t79 * t190 - t80 * t188 - t81 * t186 + (t260 * t57 + t261 * t55 + t262 * t53) * t134, t10 * t222 + t11 * t217 + t12 * t212 - t175 * t320 - t209 * t63 - t214 * t62 - t219 * t61, (t81 * t305 + t6 * t63) * t121 + (t80 * t306 + t5 * t62) * t118 + (t79 * t313 + t4 * t61) * t115 + (t260 * t9 + t261 * t8 + t262 * t7) * t320, t144 - g(1); t19 * t223 + t20 * t218 + t21 * t213, t72 * t213 + t71 * t218 + t70 * t223, t75 * t213 + t74 * t218 + t73 * t223, t16 * t223 + t17 * t218 + t18 * t213 + (-t109 * t189 - t110 * t187 - t111 * t185) * t134, (t30 * t257 + t29 * t258 + t28 * t259) * t134 + t76 * t202 + t77 * t201 + t78 * t200, t134 * t174 + t213 * t39 + t218 * t38 + t223 * t37, t34 * t223 + t35 * t218 + t36 * t213 + (t109 * t19 + t110 * t20 + t111 * t21) * t134, (t48 * t257 + t47 * t258 + t46 * t259) * t134, t76 * t304 + t77 * t302 + t78 * t300 + (t257 * t56 + t258 * t54 + t259 * t52) * t134, -t76 * t190 - t77 * t188 - t78 * t186 + (t257 * t57 + t258 * t55 + t259 * t53) * t134, t10 * t223 + t11 * t218 + t12 * t213 - t174 * t320 - t209 * t59 - t214 * t60 - t219 * t58, (t78 * t305 + t59 * t6) * t121 + (t77 * t306 + t5 * t60) * t118 + (t76 * t313 + t4 * t58) * t115 + (t257 * t9 + t258 * t8 + t259 * t7) * t320, t143 - g(2); t19 * t255 + t20 * t252 + t21 * t249, t211 + t216 + t221, t75 * t249 + t74 * t252 + t73 * t255, t16 * t255 + t17 * t252 + t18 * t249, 0.2e1 * t13 * t255 + 0.2e1 * t14 * t252 + 0.2e1 * t15 * t249, t39 * t249 + t38 * t252 + t37 * t255, t36 * t249 + t35 * t252 + t34 * t255, 0, t151 * t221 + t153 * t216 + t155 * t211, -t145 * t221 - t147 * t216 - t149 * t211, t10 * t255 + t11 * t252 + t12 * t249 - t209 * t93 - t214 * t92 - t219 * t91, (t156 * t3 + t6 * t93) * t121 + (t154 * t2 + t5 * t92) * t118 + (t1 * t152 + t4 * t91) * t115, t142 - g(3);];
tauX_reg  = t25;
