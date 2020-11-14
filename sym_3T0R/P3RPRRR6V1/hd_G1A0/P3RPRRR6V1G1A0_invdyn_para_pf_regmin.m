% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RPRRR6V1G1A0
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
% tauX_reg [3x12]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:32
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3RPRRR6V1G1A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G1A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:32:11
% EndTime: 2020-08-06 18:32:18
% DurationCPUTime: 6.21s
% Computational Cost: add. (33598->461), mult. (29534->716), div. (2391->20), fcn. (14805->86), ass. (0->365)
t235 = sin(qJ(1,3));
t241 = cos(qJ(1,3));
t139 = -t235 * g(1) + t241 * g(2);
t140 = t241 * g(1) + t235 * g(2);
t228 = legFrame(3,3);
t200 = sin(t228);
t203 = cos(t228);
t226 = sin(pkin(7));
t227 = cos(pkin(7));
t256 = 0.1e1 / pkin(3) ^ 2;
t247 = xDP(2);
t248 = xDP(1);
t249 = (pkin(6) + pkin(5));
t110 = -2 * pkin(2) * t247 - 2 * t248 * t249;
t123 = -pkin(2) * t248 + t247 * t249;
t193 = qJ(1,3) + t228;
t165 = pkin(7) + t193;
t145 = sin(t165);
t168 = sin(t193);
t171 = cos(t193);
t152 = qJ(3,3) + t165;
t133 = cos(t152);
t153 = -qJ(3,3) + t165;
t134 = cos(t153);
t315 = t133 + t134;
t127 = sin(t152);
t128 = sin(t153);
t318 = t127 + t128;
t148 = cos(t165);
t394 = 0.2e1 * t148;
t397 = -2 * pkin(1);
t61 = t123 * t394 + t110 * t145 + (t168 * t247 + t171 * t248) * t397 + (-t247 * t318 - t248 * t315) * pkin(3);
t212 = pkin(7) + qJ(3,3);
t181 = sin(t212);
t215 = -pkin(7) + qJ(3,3);
t184 = sin(t215);
t250 = 0.2e1 * qJ(3,3);
t219 = sin(t250);
t196 = pkin(3) * t219;
t234 = sin(qJ(3,3));
t417 = 2 * pkin(2);
t97 = t234 * t417 + t196 + (t181 + t184) * pkin(1);
t370 = t61 ^ 2 / t97 ^ 2;
t293 = t256 * t370;
t321 = -t139 * t203 + t140 * t200;
t199 = t227 * pkin(1);
t240 = cos(qJ(3,3));
t352 = t240 * pkin(3) + pkin(2);
t124 = t199 + t352;
t117 = 0.1e1 / t124;
t232 = xDDP(2);
t233 = xDDP(1);
t118 = 0.1e1 / t124 ^ 2;
t88 = 0.1e1 / t97;
t290 = t118 * t234 * t88;
t94 = -t145 * t248 + t148 * t247;
t85 = t94 ^ 2;
t347 = t118 * t85;
t371 = t226 * pkin(1);
t401 = t249 + t371;
t40 = 0.2e1 * t94 * t61 * t290 + (-t145 * t233 + t148 * t232 - t347 * t401) * t117;
t386 = t40 * pkin(1);
t423 = (-0.2e1 * t386 - t321) * t227 + (pkin(1) * t293 - t139 * t200 - t140 * t203) * t226 + pkin(5) * t293 - t40 * t417;
t237 = sin(qJ(1,2));
t243 = cos(qJ(1,2));
t141 = -t237 * g(1) + t243 * g(2);
t142 = t243 * g(1) + t237 * g(2);
t229 = legFrame(2,3);
t201 = sin(t229);
t204 = cos(t229);
t194 = qJ(1,2) + t229;
t166 = pkin(7) + t194;
t146 = sin(t166);
t169 = sin(t194);
t172 = cos(t194);
t156 = qJ(3,2) + t166;
t135 = cos(t156);
t157 = -qJ(3,2) + t166;
t136 = cos(t157);
t314 = t135 + t136;
t129 = sin(t156);
t130 = sin(t157);
t317 = t129 + t130;
t149 = cos(t166);
t393 = 0.2e1 * t149;
t62 = t123 * t393 + t110 * t146 + (t169 * t247 + t172 * t248) * t397 + (-t247 * t317 - t248 * t314) * pkin(3);
t213 = pkin(7) + qJ(3,2);
t182 = sin(t213);
t216 = -pkin(7) + qJ(3,2);
t185 = sin(t216);
t251 = 0.2e1 * qJ(3,2);
t220 = sin(t251);
t197 = pkin(3) * t220;
t236 = sin(qJ(3,2));
t98 = t236 * t417 + t197 + (t182 + t185) * pkin(1);
t369 = t62 ^ 2 / t98 ^ 2;
t292 = t256 * t369;
t320 = -t141 * t204 + t142 * t201;
t242 = cos(qJ(3,2));
t351 = t242 * pkin(3) + pkin(2);
t125 = t199 + t351;
t119 = 0.1e1 / t125;
t120 = 0.1e1 / t125 ^ 2;
t90 = 0.1e1 / t98;
t289 = t120 * t236 * t90;
t95 = -t146 * t248 + t149 * t247;
t86 = t95 ^ 2;
t345 = t120 * t86;
t41 = 0.2e1 * t95 * t62 * t289 + (-t146 * t233 + t149 * t232 - t345 * t401) * t119;
t385 = t41 * pkin(1);
t422 = (-0.2e1 * t385 - t320) * t227 + (pkin(1) * t292 - t141 * t201 - t142 * t204) * t226 + pkin(5) * t292 - t41 * t417;
t239 = sin(qJ(1,1));
t245 = cos(qJ(1,1));
t143 = -g(1) * t239 + g(2) * t245;
t144 = g(1) * t245 + g(2) * t239;
t230 = legFrame(1,3);
t202 = sin(t230);
t205 = cos(t230);
t195 = qJ(1,1) + t230;
t167 = pkin(7) + t195;
t147 = sin(t167);
t170 = sin(t195);
t173 = cos(t195);
t160 = qJ(3,1) + t167;
t137 = cos(t160);
t161 = -qJ(3,1) + t167;
t138 = cos(t161);
t313 = t137 + t138;
t131 = sin(t160);
t132 = sin(t161);
t316 = t131 + t132;
t150 = cos(t167);
t392 = 0.2e1 * t150;
t63 = t123 * t392 + t110 * t147 + (t170 * t247 + t173 * t248) * t397 + (-t247 * t316 - t248 * t313) * pkin(3);
t214 = pkin(7) + qJ(3,1);
t183 = sin(t214);
t217 = -pkin(7) + qJ(3,1);
t186 = sin(t217);
t252 = 0.2e1 * qJ(3,1);
t221 = sin(t252);
t198 = pkin(3) * t221;
t238 = sin(qJ(3,1));
t99 = t238 * t417 + t198 + (t183 + t186) * pkin(1);
t368 = t63 ^ 2 / t99 ^ 2;
t291 = t256 * t368;
t319 = -t143 * t205 + t144 * t202;
t244 = cos(qJ(3,1));
t350 = t244 * pkin(3) + pkin(2);
t126 = t199 + t350;
t121 = 0.1e1 / t126;
t122 = 0.1e1 / t126 ^ 2;
t92 = 0.1e1 / t99;
t288 = t122 * t238 * t92;
t96 = -t147 * t248 + t150 * t247;
t87 = t96 ^ 2;
t343 = t122 * t87;
t42 = 0.2e1 * t96 * t63 * t288 + (-t147 * t233 + t150 * t232 - t343 * t401) * t121;
t384 = t42 * pkin(1);
t421 = (t319 + 0.2e1 * t384) * t227 - (pkin(1) * t291 - t143 * t202 - t144 * t205) * t226 - pkin(5) * t291 + t42 * t417;
t416 = 0.2e1 * t240;
t415 = 0.2e1 * t242;
t414 = 0.2e1 * t244;
t180 = t199 + pkin(2);
t323 = t238 * t180;
t103 = 0.1e1 / (t198 + 0.2e1 * t323);
t365 = t63 * t92;
t413 = 0.2e1 * t103 * t365;
t324 = t236 * t180;
t102 = 0.1e1 / (t197 + 0.2e1 * t324);
t366 = t62 * t90;
t412 = 0.2e1 * t102 * t366;
t325 = t234 * t180;
t101 = 0.1e1 / (t196 + 0.2e1 * t325);
t367 = t61 * t88;
t411 = 0.2e1 * t101 * t367;
t222 = t240 ^ 2;
t254 = pkin(3) ^ 2;
t255 = 0.1e1 / pkin(3);
t387 = 2 * t249;
t277 = (pkin(1) ^ 2) + (pkin(2) ^ 2) + (t249 ^ 2) + t371 * t387;
t330 = t401 * t234;
t282 = t330 * t367;
t296 = t255 * t367;
t283 = t124 * t296;
t348 = t117 * t94;
t287 = -t348 / 0.2e1;
t309 = 0.2e1 * t199;
t322 = pkin(3) * t417;
t388 = -2 * t249;
t396 = -2 * pkin(2);
t76 = -t315 * pkin(3) + t145 * t388 + t148 * t396 + t171 * t397;
t355 = t76 * t88;
t73 = -t318 * pkin(3) + t145 * t396 + t148 * t387 + t168 * t397;
t358 = t73 * t88;
t31 = (t232 * t358 + t233 * t355 + 0.2e1 * (-t282 + (t222 * t254 + t240 * t322 + t309 * t352 + t277) * t348) / (t196 / 0.2e1 + t325) * t287 + (-t240 * t283 + t330 * t348) * t411) * t255;
t342 = t31 * t234;
t19 = t240 * t293 + t342;
t410 = -t19 / 0.2e1;
t223 = t242 ^ 2;
t329 = t401 * t236;
t280 = t329 * t366;
t295 = t255 * t366;
t281 = t125 * t295;
t346 = t119 * t95;
t286 = -t346 / 0.2e1;
t77 = -t314 * pkin(3) + t146 * t388 + t149 * t396 + t172 * t397;
t354 = t77 * t90;
t74 = -t317 * pkin(3) + t146 * t396 + t149 * t387 + t169 * t397;
t357 = t74 * t90;
t32 = (t232 * t357 + t233 * t354 + 0.2e1 * (-t280 + (t223 * t254 + t242 * t322 + t309 * t351 + t277) * t346) / (t197 / 0.2e1 + t324) * t286 + (-t242 * t281 + t329 * t346) * t412) * t255;
t341 = t32 * t236;
t21 = t242 * t292 + t341;
t409 = -t21 / 0.2e1;
t224 = t244 ^ 2;
t328 = t401 * t238;
t278 = t328 * t365;
t294 = t255 * t365;
t279 = t126 * t294;
t344 = t121 * t96;
t285 = -t344 / 0.2e1;
t78 = -t313 * pkin(3) + t147 * t388 + t150 * t396 + t173 * t397;
t353 = t78 * t92;
t75 = -t316 * pkin(3) + t147 * t396 + t150 * t387 + t170 * t397;
t356 = t75 * t92;
t33 = (t232 * t356 + t233 * t353 + 0.2e1 * (-t278 + (t224 * t254 + t244 * t322 + t309 * t350 + t277) * t344) / (t198 / 0.2e1 + t323) * t285 + (-t244 * t279 + t328 * t344) * t413) * t255;
t340 = t33 * t238;
t23 = t244 * t291 + t340;
t408 = -t23 / 0.2e1;
t159 = t252 + t167;
t162 = -0.2e1 * qJ(3,1) + t167;
t178 = qJ(3,1) + t195;
t179 = -qJ(3,1) + t195;
t311 = 2 * pkin(1);
t359 = (t316 * t417 + (sin(t179) + sin(t178)) * t311 + t313 * t388 + (sin(t162) + sin(t159) + 0.2e1 * t147) * pkin(3)) * t92;
t407 = t359 / 0.2e1;
t155 = t251 + t166;
t158 = -0.2e1 * qJ(3,2) + t166;
t176 = qJ(3,2) + t194;
t177 = -qJ(3,2) + t194;
t360 = (t317 * t417 + (sin(t177) + sin(t176)) * t311 + t314 * t388 + (sin(t158) + sin(t155) + 0.2e1 * t146) * pkin(3)) * t90;
t406 = t360 / 0.2e1;
t151 = t250 + t165;
t154 = -0.2e1 * qJ(3,3) + t165;
t174 = qJ(3,3) + t193;
t175 = -qJ(3,3) + t193;
t361 = (t318 * t417 + (sin(t175) + sin(t174)) * t311 + t315 * t388 + (sin(t154) + sin(t151) + 0.2e1 * t145) * pkin(3)) * t88;
t405 = t361 / 0.2e1;
t362 = (t313 * t417 + (cos(t179) + cos(t178)) * t311 + t316 * t387 + (cos(t162) + cos(t159) + t392) * pkin(3)) * t92;
t404 = t362 / 0.2e1;
t363 = (t314 * t417 + (cos(t177) + cos(t176)) * t311 + t317 * t387 + (cos(t158) + cos(t155) + t393) * pkin(3)) * t90;
t403 = t363 / 0.2e1;
t364 = (t315 * t417 + (cos(t175) + cos(t174)) * t311 + t318 * t387 + (cos(t154) + cos(t151) + t394) * pkin(3)) * t88;
t402 = t364 / 0.2e1;
t231 = xDDP(3);
t349 = (t231 - g(3));
t391 = -0.2e1 * t222;
t390 = -0.2e1 * t223;
t389 = -0.2e1 * t224;
t383 = t127 / 0.2e1;
t382 = t129 / 0.2e1;
t381 = t131 / 0.2e1;
t380 = t134 / 0.2e1;
t379 = t136 / 0.2e1;
t378 = t138 / 0.2e1;
t377 = t184 / 0.2e1;
t376 = t185 / 0.2e1;
t375 = t186 / 0.2e1;
t374 = cos(t212) / 0.2e1;
t373 = cos(t213) / 0.2e1;
t372 = cos(t214) / 0.2e1;
t339 = t40 * t234;
t338 = t41 * t236;
t337 = t42 * t238;
t336 = t117 * t145;
t335 = t117 * t148;
t334 = t119 * t146;
t333 = t119 * t149;
t332 = t121 * t147;
t331 = t121 * t150;
t327 = t232 / 0.2e1;
t326 = t233 / 0.2e1;
t312 = 2 * t349;
t100 = pkin(2) * t309 + t254 / 0.2e1 + t277;
t305 = ((-t100 * t348 + t282) * t416 + (-cos(t250) + (t391 - 0.1e1) * t180 * t117) * t94 * pkin(3)) * t348;
t304 = ((-t100 * t346 + t280) * t415 + (-cos(t251) + (t390 - 0.1e1) * t180 * t119) * t95 * pkin(3)) * t346;
t303 = ((-t100 * t344 + t278) * t414 + (-cos(t252) + (t389 - 0.1e1) * t180 * t121) * t96 * pkin(3)) * t344;
t302 = t88 * t339;
t301 = t90 * t338;
t300 = t92 * t337;
t299 = t240 * t40 * t88;
t298 = t242 * t41 * t90;
t297 = t244 * t42 * t92;
t276 = t296 * t348;
t275 = t295 * t346;
t274 = t294 * t344;
t273 = t240 * t85 * t290;
t272 = t242 * t86 * t289;
t271 = t244 * t87 * t288;
t270 = 0.2e1 * t276;
t269 = 0.2e1 * t275;
t268 = 0.2e1 * t274;
t49 = t219 * t287 * t401 + t283;
t52 = t326 * t364;
t55 = t327 * t361;
t264 = -t101 * t305 + t49 * t411 + t52 + t55;
t50 = t220 * t286 * t401 + t281;
t53 = t326 * t363;
t56 = t327 * t360;
t263 = -t102 * t304 + t50 * t412 + t53 + t56;
t51 = t221 * t285 * t401 + t279;
t54 = t326 * t362;
t57 = t327 * t359;
t262 = -t103 * t303 + t51 * t413 + t54 + t57;
t192 = cos(t217);
t191 = cos(t216);
t190 = cos(t215);
t116 = t205 * g(1) + t202 * g(2);
t115 = t204 * g(1) + t201 * g(2);
t114 = t203 * g(1) + t200 * g(2);
t113 = t202 * g(1) - t205 * g(2);
t112 = t201 * g(1) - t204 * g(2);
t111 = t200 * g(1) - t203 * g(2);
t84 = -t113 * t239 + t116 * t245;
t83 = -t112 * t237 + t115 * t243;
t82 = -t111 * t235 + t114 * t241;
t81 = t113 * t245 + t116 * t239;
t80 = t112 * t243 + t115 * t237;
t79 = t111 * t241 + t114 * t235;
t72 = (t389 + 0.1e1) * t343;
t71 = (t390 + 0.1e1) * t345;
t70 = (t391 + 0.1e1) * t347;
t39 = t319 + t384;
t38 = t320 + t385;
t37 = t321 + t386;
t36 = (t244 * t268 + t337) * t238;
t35 = (t242 * t269 + t338) * t236;
t34 = (t240 * t270 + t339) * t234;
t30 = t33 * t244;
t29 = t32 * t242;
t28 = t31 * t240;
t27 = t337 * t414 + (0.4e1 * t224 - 0.2e1) * t274;
t26 = t338 * t415 + (0.4e1 * t223 - 0.2e1) * t275;
t25 = t339 * t416 + (0.4e1 * t222 - 0.2e1) * t276;
t24 = -t238 * t291 + t30;
t22 = -t236 * t292 + t29;
t20 = -t234 * t293 + t28;
t18 = t262 + t349;
t17 = t263 + t349;
t16 = t264 + t349;
t15 = pkin(2) * t268 + t33 * pkin(5) + (t226 * t33 + t227 * t268) * pkin(1);
t14 = pkin(2) * t269 + t32 * pkin(5) + (t226 * t32 + t227 * t269) * pkin(1);
t13 = pkin(2) * t270 + t31 * pkin(5) + (t226 * t31 + t227 * t270) * pkin(1);
t12 = -t238 * t18 + (t132 / 0.2e1 + t381) * g(2) + (t378 + t137 / 0.2e1) * g(1) + (-pkin(5) * t244 + (t375 - t183 / 0.2e1) * pkin(1)) * t42 + (pkin(2) * t244 + (t192 / 0.2e1 + t372) * pkin(1)) * t343;
t11 = -t236 * t17 + (t130 / 0.2e1 + t382) * g(2) + (t379 + t135 / 0.2e1) * g(1) + (-pkin(5) * t242 + (t376 - t182 / 0.2e1) * pkin(1)) * t41 + (pkin(2) * t242 + (t191 / 0.2e1 + t373) * pkin(1)) * t345;
t10 = -t234 * t16 + (t128 / 0.2e1 + t383) * g(2) + (t380 + t133 / 0.2e1) * g(1) + (-pkin(5) * t240 + (t377 - t181 / 0.2e1) * pkin(1)) * t40 + (pkin(2) * t240 + (t190 / 0.2e1 + t374) * pkin(1)) * t347;
t9 = (0.2e1 * t54 + 0.2e1 * t57 + (0.4e1 * t365 * t51 - 0.2e1 * t303) * t103 + t312) * t244 / 0.2e1 - t238 * (-pkin(2) * t343 + t42 * pkin(5)) + (t378 - t137 / 0.2e1) * g(2) + (-t132 / 0.2e1 + t381) * g(1) + ((-t192 / 0.2e1 + t372) * t42 + (t375 + t183 / 0.2e1) * t343) * pkin(1);
t8 = (0.2e1 * t53 + 0.2e1 * t56 + (0.4e1 * t366 * t50 - 0.2e1 * t304) * t102 + t312) * t242 / 0.2e1 - t236 * (-pkin(2) * t345 + t41 * pkin(5)) + (t379 - t135 / 0.2e1) * g(2) + (-t130 / 0.2e1 + t382) * g(1) + ((-t191 / 0.2e1 + t373) * t41 + (t376 + t182 / 0.2e1) * t345) * pkin(1);
t7 = (0.2e1 * t52 + 0.2e1 * t55 + (0.4e1 * t367 * t49 - 0.2e1 * t305) * t101 + t312) * t240 / 0.2e1 - t234 * (-pkin(2) * t347 + t40 * pkin(5)) + (t380 - t133 / 0.2e1) * g(2) + (-t128 / 0.2e1 + t383) * g(1) + ((-t190 / 0.2e1 + t374) * t40 + (t377 + t181 / 0.2e1) * t347) * pkin(1);
t6 = -t238 * t15 + t421 * t244;
t5 = -t244 * t15 - t421 * t238;
t4 = -t236 * t14 - t422 * t242;
t3 = -t242 * t14 + t422 * t236;
t2 = -t234 * t13 - t423 * t240;
t1 = -t240 * t13 + t423 * t234;
t43 = [-t332 * t42 - t334 * t41 - t336 * t40, -t332 * t81 - t334 * t80 - t336 * t79, -t332 * t84 - t334 * t83 - t336 * t82, t16 * t402 + t17 * t403 + t18 * t404 + (-t332 * t39 - t334 * t38 - t336 * t37) * pkin(1), -t34 * t336 - t35 * t334 - t36 * t332 + (-t271 * t78 - t272 * t77 - t273 * t76) * t255, -t25 * t336 - t26 * t334 - t27 * t332 + (t353 * t72 + t354 * t71 + t355 * t70) * t255, -t19 * t336 - t21 * t334 - t23 * t332 + (t300 * t78 + t301 * t77 + t302 * t76) * t255, -t20 * t336 - t22 * t334 - t24 * t332 + (t297 * t78 + t298 * t77 + t299 * t76) * t255, (t31 * t355 + t32 * t354 + t33 * t353) * t255, -t2 * t336 - t4 * t334 - t6 * t332 + (t353 * t9 + t354 * t8 + t355 * t7) * t255 + t20 * t402 + t22 * t403 + t24 * t404, -t1 * t336 - t3 * t334 - t5 * t332 + (t10 * t355 + t11 * t354 + t12 * t353) * t255 + t364 * t410 + t363 * t409 + t362 * t408, t233 - g(1); t331 * t42 + t333 * t41 + t335 * t40, t331 * t81 + t333 * t80 + t335 * t79, t331 * t84 + t333 * t83 + t335 * t82, t16 * t405 + t17 * t406 + t18 * t407 + (t331 * t39 + t333 * t38 + t335 * t37) * pkin(1), t34 * t335 + t35 * t333 + t36 * t331 + (-t271 * t75 - t272 * t74 - t273 * t73) * t255, t25 * t335 + t26 * t333 + t27 * t331 + (t356 * t72 + t357 * t71 + t358 * t70) * t255, t19 * t335 + t21 * t333 + t23 * t331 + (t300 * t75 + t301 * t74 + t302 * t73) * t255, t20 * t335 + t22 * t333 + t24 * t331 + (t297 * t75 + t298 * t74 + t299 * t73) * t255, (t31 * t358 + t32 * t357 + t33 * t356) * t255, t2 * t335 + t4 * t333 + t6 * t331 + (t356 * t9 + t357 * t8 + t358 * t7) * t255 + t20 * t405 + t22 * t406 + t24 * t407, t1 * t335 + t3 * t333 + t5 * t331 + (t10 * t358 + t11 * t357 + t12 * t356) * t255 + t361 * t410 + t360 * t409 + t359 * t408, t232 - g(2); 0, 0, 0, -(3 * g(3)) + (3 * t231) + t262 + t263 + t264, 0, 0, 0, 0, 0, t28 + t29 + t30 + (-t234 * t370 - t236 * t369 - t238 * t368) * t256, -t342 - t341 - t340 + (-t240 * t370 - t242 * t369 - t244 * t368) * t256, t349;];
tauX_reg  = t43;
