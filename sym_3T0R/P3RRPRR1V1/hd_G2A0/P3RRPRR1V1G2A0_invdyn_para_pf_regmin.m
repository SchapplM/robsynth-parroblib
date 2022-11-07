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
% tauX_reg [3x15]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2022-11-04 17:08
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
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
% StartTime: 2022-11-04 17:08:08
% EndTime: 2022-11-04 17:08:13
% DurationCPUTime: 6.06s
% Computational Cost: add. (14110->424), mult. (21426->787), div. (3702->17), fcn. (18183->18), ass. (0->347)
t179 = sin(qJ(1,3));
t185 = cos(qJ(1,3));
t172 = legFrame(3,2);
t139 = cos(t172);
t373 = t139 * g(1);
t136 = sin(t172);
t376 = t136 * g(2);
t235 = t373 - t376;
t103 = g(3) * t185 + t235 * t179;
t169 = pkin(3) + qJ(3,3);
t146 = 0.1e1 / t169 ^ 2;
t184 = cos(qJ(2,3));
t155 = 0.1e1 / t184;
t193 = pkin(1) + pkin(2);
t301 = t184 * t193;
t127 = t179 * t169 + t185 * t301;
t178 = sin(qJ(2,3));
t190 = xDP(3);
t191 = xDP(2);
t192 = xDP(1);
t308 = t179 * t184;
t313 = t169 * t185;
t52 = (-t192 * t313 + (t178 * t191 + t192 * t308) * t193) * t139 + (t191 * t313 - (-t178 * t192 + t191 * t308) * t193) * t136 + t127 * t190;
t115 = t136 * t192 + t139 * t191;
t310 = t178 * t115;
t70 = (t185 * t190 + (-t136 * t191 + t139 * t192) * t179) * t184 + t310;
t292 = t155 * t52 * t70;
t409 = -0.2e1 * t146 * t292 + t103;
t181 = sin(qJ(1,2));
t187 = cos(qJ(1,2));
t173 = legFrame(2,2);
t140 = cos(t173);
t372 = t140 * g(1);
t137 = sin(t173);
t375 = t137 * g(2);
t234 = t372 - t375;
t104 = g(3) * t187 + t234 * t181;
t170 = pkin(3) + qJ(3,2);
t149 = 0.1e1 / t170 ^ 2;
t186 = cos(qJ(2,2));
t159 = 0.1e1 / t186;
t300 = t186 * t193;
t128 = t181 * t170 + t187 * t300;
t180 = sin(qJ(2,2));
t305 = t181 * t186;
t312 = t170 * t187;
t53 = (-t192 * t312 + (t180 * t191 + t192 * t305) * t193) * t140 + (t191 * t312 - (-t180 * t192 + t191 * t305) * t193) * t137 + t128 * t190;
t116 = t137 * t192 + t140 * t191;
t307 = t180 * t116;
t71 = (t187 * t190 + (-t137 * t191 + t140 * t192) * t181) * t186 + t307;
t291 = t159 * t53 * t71;
t408 = -0.2e1 * t149 * t291 + t104;
t183 = sin(qJ(1,1));
t189 = cos(qJ(1,1));
t174 = legFrame(1,2);
t141 = cos(t174);
t371 = t141 * g(1);
t138 = sin(t174);
t374 = t138 * g(2);
t233 = t371 - t374;
t105 = g(3) * t189 + t233 * t183;
t171 = pkin(3) + qJ(3,1);
t152 = 0.1e1 / t171 ^ 2;
t188 = cos(qJ(2,1));
t163 = 0.1e1 / t188;
t299 = t188 * t193;
t129 = t183 * t171 + t189 * t299;
t182 = sin(qJ(2,1));
t302 = t183 * t188;
t311 = t171 * t189;
t54 = (-t192 * t311 + (t182 * t191 + t192 * t302) * t193) * t141 + (t191 * t311 - (-t182 * t192 + t191 * t302) * t193) * t138 + t129 * t190;
t117 = t138 * t192 + t141 * t191;
t304 = t182 * t117;
t72 = (t189 * t190 + (-t138 * t191 + t141 * t192) * t183) * t188 + t304;
t290 = t163 * t54 * t72;
t407 = -0.2e1 * t152 * t290 + t105;
t106 = -t136 * t308 + t178 * t139;
t109 = t136 * t178 + t139 * t308;
t112 = t115 ^ 2;
t145 = 0.1e1 / t169;
t147 = t145 * t146;
t154 = t184 ^ 2;
t157 = t155 / t154;
t167 = 0.1e1 / t193;
t175 = xDDP(3);
t176 = xDDP(2);
t177 = xDDP(1);
t37 = -t147 * t292 + (t167 * t112 * t157 + t185 * t175 + (t109 * t177 + t106 * t176 - (-t193 * t70 + t52) * t146 * t70) * t155) * t145;
t349 = t37 * qJ(3,3);
t209 = -t349 + t409;
t107 = -t137 * t305 + t180 * t140;
t110 = t137 * t180 + t140 * t305;
t113 = t116 ^ 2;
t148 = 0.1e1 / t170;
t150 = t148 * t149;
t158 = t186 ^ 2;
t161 = t159 / t158;
t38 = -t150 * t291 + (t167 * t113 * t161 + t187 * t175 + (t110 * t177 + t107 * t176 - (-t193 * t71 + t53) * t149 * t71) * t159) * t148;
t348 = t38 * qJ(3,2);
t210 = -t348 + t408;
t108 = -t138 * t302 + t182 * t141;
t111 = t138 * t182 + t141 * t302;
t114 = t117 ^ 2;
t151 = 0.1e1 / t171;
t153 = t151 * t152;
t162 = t188 ^ 2;
t165 = t163 / t162;
t39 = -t153 * t290 + (t167 * t114 * t165 + t189 * t175 + (t111 * t177 + t108 * t176 - (-t193 * t72 + t54) * t152 * t72) * t163) * t151;
t347 = t39 * qJ(3,1);
t211 = -t347 + t407;
t379 = g(3) * t179;
t397 = -t235 * t185 + t379;
t378 = g(3) * t181;
t396 = -t234 * t187 + t378;
t377 = g(3) * t183;
t395 = -t233 * t189 + t377;
t314 = t151 * t189;
t260 = t395 * t314;
t316 = t148 * t187;
t261 = t396 * t316;
t318 = t145 * t185;
t262 = t397 * t318;
t329 = t117 * t151;
t51 = (-t182 * t72 + t117) * t163 ^ 2 * t329;
t215 = t183 * t299 - t311;
t303 = t182 * t193;
t86 = -t215 * t138 + t141 * t303;
t80 = t86 * t151 * t176;
t90 = t138 * t303 + t215 * t141;
t84 = t90 * t151 * t177;
t99 = t129 * t151 * t175;
t406 = -t51 - t80 - t84 - t99 + t395;
t330 = t116 * t148;
t50 = (-t180 * t71 + t116) * t159 ^ 2 * t330;
t216 = t181 * t300 - t312;
t306 = t180 * t193;
t87 = -t216 * t137 + t140 * t306;
t81 = t87 * t148 * t176;
t89 = t137 * t306 + t216 * t140;
t83 = t89 * t148 * t177;
t98 = t128 * t148 * t175;
t405 = -t50 - t81 - t83 - t98 + t396;
t331 = t115 * t145;
t49 = (-t178 * t70 + t115) * t155 ^ 2 * t331;
t217 = t179 * t301 - t313;
t309 = t178 * t193;
t85 = -t217 * t136 + t139 * t309;
t79 = t85 * t145 * t176;
t88 = t136 * t309 + t217 * t139;
t82 = t88 * t145 * t177;
t97 = t127 * t145 * t175;
t404 = -t49 - t79 - t82 - t97 + t397;
t156 = 0.1e1 / t184 ^ 2;
t168 = 0.1e1 / t193 ^ 2;
t334 = t112 * t168;
t238 = t156 * qJ(3,3) * t334;
t400 = -t238 + t404;
t160 = 0.1e1 / t186 ^ 2;
t333 = t113 * t168;
t237 = t160 * qJ(3,2) * t333;
t399 = -t237 + t405;
t164 = 0.1e1 / t188 ^ 2;
t332 = t114 * t168;
t236 = t164 * qJ(3,1) * t332;
t398 = -t236 + t406;
t388 = -2 * pkin(1);
t387 = 2 * pkin(1);
t268 = t178 * t334;
t67 = t157 * t268 + (t136 * t177 + t139 * t176) * t167 * t155;
t386 = t67 * pkin(1);
t267 = t180 * t333;
t68 = t161 * t267 + (t137 * t177 + t140 * t176) * t167 * t159;
t385 = t68 * pkin(1);
t266 = t182 * t332;
t69 = t165 * t266 + (t138 * t177 + t141 * t176) * t167 * t163;
t384 = t69 * pkin(1);
t383 = t156 - 0.2e1;
t382 = t160 - 0.2e1;
t381 = t164 - 0.2e1;
t380 = pkin(1) * t167;
t370 = t145 * t397;
t64 = t70 ^ 2;
t369 = t146 * t64;
t368 = t148 * t396;
t65 = t71 ^ 2;
t367 = t149 * t65;
t366 = t151 * t395;
t66 = t72 ^ 2;
t365 = t152 * t66;
t364 = t154 * t37;
t363 = t156 * t64;
t362 = t158 * t38;
t361 = t160 * t65;
t360 = t162 * t39;
t359 = t164 * t66;
t358 = t178 * t37;
t357 = t180 * t38;
t356 = t182 * t39;
t355 = t184 * t37;
t354 = t186 * t38;
t353 = t188 * t39;
t352 = t193 * t52;
t351 = t193 * t53;
t350 = t193 * t54;
t346 = t67 * t178;
t345 = t67 * t184;
t344 = t68 * t180;
t343 = t68 * t186;
t342 = t69 * t182;
t341 = t69 * t188;
t340 = t106 * t155;
t339 = t107 * t159;
t338 = t108 * t163;
t337 = t109 * t155;
t336 = t110 * t159;
t335 = t111 * t163;
t130 = t136 * g(1) + t139 * g(2);
t328 = t130 * t178;
t118 = t130 * t184;
t131 = t137 * g(1) + t140 * g(2);
t327 = t131 * t180;
t121 = t131 * t186;
t132 = t138 * g(1) + t141 * g(2);
t326 = t132 * t182;
t124 = t132 * t188;
t325 = t136 * t155;
t324 = t137 * t159;
t323 = t138 * t163;
t322 = t139 * t155;
t321 = t140 * t159;
t320 = t141 * t163;
t319 = t145 * t155;
t317 = t148 * t159;
t315 = t151 * t163;
t298 = pkin(1) * t369;
t297 = pkin(1) * t367;
t296 = pkin(1) * t365;
t295 = pkin(1) * t355;
t294 = pkin(1) * t354;
t293 = pkin(1) * t353;
t289 = t146 * t363;
t288 = t147 * t363;
t287 = t147 * t352;
t286 = t149 * t361;
t285 = t150 * t361;
t284 = t150 * t351;
t283 = t152 * t359;
t282 = t153 * t359;
t281 = t153 * t350;
t280 = t155 * t358;
t279 = t159 * t357;
t278 = t163 * t356;
t166 = t193 ^ 2;
t277 = (t169 * (-t70 + t310) * t155 + (-t166 * t70 + t352) * t184 * t145) * t146 * t155;
t276 = (t170 * (-t71 + t307) * t159 + (-t166 * t71 + t351) * t186 * t148) * t149 * t159;
t275 = (t171 * (-t72 + t304) * t163 + (-t166 * t72 + t350) * t188 * t151) * t152 * t163;
t274 = t106 * t319;
t273 = t107 * t317;
t272 = t108 * t315;
t271 = t109 * t319;
t270 = t110 * t317;
t269 = t111 * t315;
t265 = t167 * t331;
t264 = t167 * t330;
t263 = t167 * t329;
t250 = t70 * t265;
t26 = t178 * t355 - t383 * t250;
t259 = 0.2e1 * t26 * t319;
t249 = t71 * t264;
t25 = t180 * t354 - t382 * t249;
t258 = 0.2e1 * t25 * t317;
t248 = t72 * t263;
t27 = t182 * t353 - t381 * t248;
t257 = 0.2e1 * t27 * t315;
t256 = t155 * t298;
t255 = t159 * t297;
t254 = t163 * t296;
t247 = t178 * t289;
t246 = t180 * t286;
t245 = t182 * t283;
t244 = t397 * t274;
t243 = t396 * t273;
t242 = t395 * t272;
t241 = t397 * t271;
t240 = t396 * t270;
t239 = t395 * t269;
t232 = t156 * t178 * t265;
t231 = t160 * t180 * t264;
t230 = t164 * t182 * t263;
t226 = pkin(1) * t232;
t225 = pkin(1) * t231;
t224 = pkin(1) * t230;
t58 = t155 * t334 + t346;
t59 = t159 * t333 + t344;
t60 = t163 * t332 + t342;
t214 = t277 + t287;
t213 = t276 + t284;
t212 = t275 + t281;
t35 = 0.2e1 * t155 * t250 + t358;
t34 = 0.2e1 * t159 * t249 + t357;
t36 = 0.2e1 * t163 * t248 + t356;
t205 = t136 * t280 + t137 * t279 + t138 * t278;
t204 = t139 * t280 + t140 * t279 + t141 * t278;
t194 = pkin(1) ^ 2;
t78 = t105 * t188 + t326;
t77 = t105 * t182 - t124;
t76 = t104 * t186 + t327;
t75 = t104 * t180 - t121;
t74 = t103 * t184 + t328;
t73 = t103 * t178 - t118;
t57 = -t164 * t266 + t341;
t56 = -t160 * t267 + t343;
t55 = -t156 * t268 + t345;
t48 = t381 * t365;
t47 = t382 * t367;
t46 = t383 * t369;
t33 = 0.2e1 * t71 * t231 - t354;
t32 = t36 * t182;
t31 = 0.2e1 * t72 * t230 - t353;
t30 = t34 * t180;
t29 = t35 * t178;
t28 = 0.2e1 * t70 * t232 - t355;
t24 = -t60 * pkin(1) + 0.2e1 * t347 - t407;
t23 = -t59 * pkin(1) + 0.2e1 * t348 - t408;
t22 = -t58 * pkin(1) + 0.2e1 * t349 - t409;
t21 = (t211 + t254) * t182 - t124 + t384;
t20 = (t211 + 0.2e1 * t254) * t182 - t124 + 0.2e1 * t384;
t19 = (t210 + t255) * t180 - t121 + t385;
t18 = (t210 + 0.2e1 * t255) * t180 - t121 + 0.2e1 * t385;
t17 = (t209 + t256) * t178 - t118 + t386;
t16 = (t209 + 0.2e1 * t256) * t178 - t118 + 0.2e1 * t386;
t15 = t211 * t188 - t381 * t296 + t326;
t14 = t210 * t186 - t382 * t297 + t327;
t13 = t209 * t184 - t383 * t298 + t328;
t12 = -qJ(3,1) * t283 - t293 + (-t212 + 0.2e1 * t224) * t72 - t406;
t11 = -qJ(3,2) * t286 - t294 + (-t213 + 0.2e1 * t225) * t71 - t405;
t10 = -qJ(3,3) * t289 - t295 + (-t214 + 0.2e1 * t226) * t70 - t404;
t9 = -qJ(3,1) * t342 + t360 * t387 + ((t212 - 0.4e1 * t224) * t72 + t398) * t188;
t8 = -qJ(3,2) * t344 + t362 * t387 + ((t213 - 0.4e1 * t225) * t71 + t399) * t186;
t7 = -qJ(3,3) * t346 + t364 * t387 + ((t214 - 0.4e1 * t226) * t70 + t400) * t184;
t6 = (-0.2e1 * t293 - t398) * t182 - qJ(3,1) * t341 + (-t212 * t182 + (0.2e1 * t164 - 0.4e1) * pkin(1) * t263) * t72;
t5 = (-0.2e1 * t294 - t399) * t180 - qJ(3,2) * t343 + (-t213 * t180 + (0.2e1 * t160 - 0.4e1) * pkin(1) * t264) * t71;
t4 = (-0.2e1 * t295 - t400) * t178 - qJ(3,3) * t345 + (-t214 * t178 + (0.2e1 * t156 - 0.4e1) * pkin(1) * t265) * t70;
t3 = t194 * t360 + ((t371 / 0.2e1 - t374 / 0.2e1) * t189 + t236 / 0.2e1 - t377 / 0.2e1 + t84 / 0.2e1 + t80 / 0.2e1 + t99 / 0.2e1 + t51 / 0.2e1 + (t224 - t275 / 0.2e1 - t281 / 0.2e1) * t72) * t188 * t388 - qJ(3,1) * (pkin(1) * t342 + t211);
t2 = t194 * t362 + ((t372 / 0.2e1 - t375 / 0.2e1) * t187 + t237 / 0.2e1 - t378 / 0.2e1 + t83 / 0.2e1 + t81 / 0.2e1 + t98 / 0.2e1 + t50 / 0.2e1 + (t225 - t276 / 0.2e1 - t284 / 0.2e1) * t71) * t186 * t388 - qJ(3,2) * (pkin(1) * t344 + t210);
t1 = t194 * t364 + ((t373 / 0.2e1 - t376 / 0.2e1) * t185 + t238 / 0.2e1 - t379 / 0.2e1 + t82 / 0.2e1 + t79 / 0.2e1 + t97 / 0.2e1 + t49 / 0.2e1 + (t226 - t277 / 0.2e1 - t287 / 0.2e1) * t70) * t184 * t388 - qJ(3,3) * (pkin(1) * t346 + t209);
t40 = [t39 * t269 + t38 * t270 + t37 * t271, t239 + t240 + t241, t103 * t271 + t104 * t270 + t105 * t269, t29 * t271 + t30 * t270 + t32 * t269 + (-t136 * t247 - t137 * t246 - t138 * t245) * t167, (t48 * t323 + t47 * t324 + t46 * t325) * t167 + t109 * t259 + t110 * t258 + t111 * t257, t205 * t167 + t60 * t269 + t59 * t270 + t58 * t271, t55 * t271 + t56 * t270 + t57 * t269 + (t136 * t37 + t137 * t38 + t138 * t39) * t167, (t69 * t323 + t68 * t324 + t67 * t325) * t167, t109 * t370 + t110 * t368 + t111 * t366 + (t77 * t323 + t75 * t324 + t73 * t325) * t167, -t178 * t241 - t180 * t240 - t182 * t239 + (t78 * t323 + t76 * t324 + t74 * t325) * t167, (t31 * t90 + t9 * t335) * t151 + (t33 * t89 + t8 * t336) * t148 + (t28 * t88 + t7 * t337) * t145 + (t16 * t325 + t18 * t324 + t20 * t323) * t167, (t6 * t335 + t36 * t90) * t151 + (t5 * t336 + t34 * t89) * t148 + (t4 * t337 + t35 * t88) * t145 + (t13 * t325 + t14 * t324 + t15 * t323) * t167, -t205 * t380 + t22 * t271 + t23 * t270 + t24 * t269 - t90 * t282 - t89 * t285 - t88 * t288, (t12 * t90 + t3 * t335) * t151 + (t11 * t89 + t2 * t336) * t148 + (t1 * t337 + t10 * t88) * t145 + (t17 * t325 + t19 * t324 + t21 * t323) * t380, t177 - g(1); t39 * t272 + t38 * t273 + t37 * t274, t242 + t243 + t244, t103 * t274 + t104 * t273 + t105 * t272, t29 * t274 + t30 * t273 + t32 * t272 + (-t139 * t247 - t140 * t246 - t141 * t245) * t167, (t48 * t320 + t47 * t321 + t46 * t322) * t167 + t106 * t259 + t107 * t258 + t108 * t257, t204 * t167 + t60 * t272 + t59 * t273 + t58 * t274, t55 * t274 + t56 * t273 + t57 * t272 + (t139 * t37 + t140 * t38 + t141 * t39) * t167, (t69 * t320 + t68 * t321 + t67 * t322) * t167, t106 * t370 + t107 * t368 + t108 * t366 + (t77 * t320 + t75 * t321 + t73 * t322) * t167, -t178 * t244 - t180 * t243 - t182 * t242 + (t78 * t320 + t76 * t321 + t74 * t322) * t167, (t31 * t86 + t9 * t338) * t151 + (t33 * t87 + t8 * t339) * t148 + (t28 * t85 + t7 * t340) * t145 + (t16 * t322 + t18 * t321 + t20 * t320) * t167, (t6 * t338 + t36 * t86) * t151 + (t5 * t339 + t34 * t87) * t148 + (t4 * t340 + t35 * t85) * t145 + (t13 * t322 + t14 * t321 + t15 * t320) * t167, -t204 * t380 + t22 * t274 + t23 * t273 + t24 * t272 - t86 * t282 - t87 * t285 - t85 * t288, (t12 * t86 + t3 * t338) * t151 + (t11 * t87 + t2 * t339) * t148 + (t1 * t340 + t10 * t85) * t145 + (t17 * t322 + t19 * t321 + t21 * t320) * t380, t176 - g(2); t39 * t314 + t38 * t316 + t37 * t318, t260 + t261 + t262, t103 * t318 + t104 * t316 + t105 * t314, t29 * t318 + t30 * t316 + t32 * t314, 0.2e1 * t25 * t316 + 0.2e1 * t26 * t318 + 0.2e1 * t27 * t314, t60 * t314 + t59 * t316 + t58 * t318, t57 * t314 + t56 * t316 + t55 * t318, 0, t184 * t262 + t186 * t261 + t188 * t260, -t178 * t262 - t180 * t261 - t182 * t260, (t129 * t31 + t189 * t9) * t151 + (t128 * t33 + t187 * t8) * t148 + (t127 * t28 + t185 * t7) * t145, (t129 * t36 + t189 * t6) * t151 + (t128 * t34 + t187 * t5) * t148 + (t127 * t35 + t185 * t4) * t145, -t127 * t288 - t128 * t285 - t129 * t282 + t22 * t318 + t23 * t316 + t24 * t314, (t12 * t129 + t189 * t3) * t151 + (t11 * t128 + t187 * t2) * t148 + (t1 * t185 + t10 * t127) * t145, t175 - g(3);];
tauX_reg  = t40;
