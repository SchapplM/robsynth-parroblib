% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PRRRR8V2G3A0
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2020-08-06 18:05
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3PRRRR8V2G3A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_regmin: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G3A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:04:49
% EndTime: 2020-08-06 18:05:02
% DurationCPUTime: 12.41s
% Computational Cost: add. (66759->478), mult. (135321->923), div. (5832->17), fcn. (125931->22), ass. (0->393)
t250 = cos(qJ(2,3));
t258 = pkin(7) + pkin(6);
t215 = t250 * t258;
t244 = sin(qJ(2,3));
t192 = pkin(2) * t244 - t215;
t234 = sin(pkin(4));
t236 = cos(pkin(4));
t243 = sin(qJ(3,3));
t369 = t236 * t243;
t169 = pkin(3) * t369 + t192 * t234;
t249 = cos(qJ(3,3));
t380 = t234 * t244;
t230 = t249 ^ 2;
t447 = pkin(3) * t230;
t136 = pkin(2) * t369 + t169 * t249 + t380 * t447;
t127 = 0.1e1 / t136;
t252 = cos(qJ(2,2));
t216 = t252 * t258;
t246 = sin(qJ(2,2));
t193 = pkin(2) * t246 - t216;
t245 = sin(qJ(3,2));
t367 = t236 * t245;
t170 = pkin(3) * t367 + t193 * t234;
t251 = cos(qJ(3,2));
t378 = t234 * t246;
t231 = t251 ^ 2;
t446 = pkin(3) * t231;
t137 = pkin(2) * t367 + t170 * t251 + t378 * t446;
t130 = 0.1e1 / t137;
t254 = cos(qJ(2,1));
t217 = t254 * t258;
t248 = sin(qJ(2,1));
t194 = pkin(2) * t248 - t217;
t247 = sin(qJ(3,1));
t365 = t236 * t247;
t171 = pkin(3) * t365 + t194 * t234;
t253 = cos(qJ(3,1));
t376 = t234 * t248;
t232 = t253 ^ 2;
t445 = pkin(3) * t232;
t138 = pkin(2) * t365 + t171 * t253 + t376 * t445;
t133 = 0.1e1 / t138;
t237 = legFrame(3,2);
t221 = sin(t237);
t224 = cos(t237);
t241 = xDDP(2);
t242 = xDDP(1);
t477 = t127 * (-t221 * t241 + t224 * t242);
t238 = legFrame(2,2);
t222 = sin(t238);
t225 = cos(t238);
t476 = t130 * (-t222 * t241 + t225 * t242);
t239 = legFrame(1,2);
t223 = sin(t239);
t226 = cos(t239);
t475 = t133 * (-t223 * t241 + t226 * t242);
t235 = cos(pkin(8));
t292 = g(1) * t226 - g(2) * t223;
t233 = sin(pkin(8));
t441 = t233 * g(3);
t168 = t292 * t235 - t441;
t382 = t233 * t236;
t187 = g(1) * t234 + g(2) * t382;
t188 = g(1) * t382 - t234 * g(2);
t371 = t235 * t236;
t201 = g(3) * t371;
t364 = t236 * t248;
t180 = t233 * t254 + t235 * t364;
t214 = t248 * t258;
t197 = pkin(2) * t254 + t214;
t377 = t234 * t247;
t274 = pkin(3) * t377 - t194 * t236;
t120 = -t180 * t445 - t197 * t233 * t253 + (pkin(2) * t377 + t274 * t253) * t235;
t260 = 0.1e1 / pkin(3);
t256 = xDP(2);
t257 = xDP(1);
t183 = t223 * t256 - t226 * t257;
t255 = xDP(3);
t208 = t255 * t233;
t159 = t183 * t235 + t208;
t370 = t235 * t255;
t162 = t183 * t233 - t370;
t442 = pkin(3) * t253;
t211 = pkin(2) + t442;
t186 = t211 * t248 - t217;
t117 = (t211 * t254 + t214) * t159 * t236 - t162 * t186;
t191 = t211 * t365;
t373 = t234 * t253;
t150 = t186 * t373 + t191;
t417 = t117 / t150;
t323 = t260 * t417;
t325 = t133 * t417;
t111 = t159 * t373 + (t159 * t364 + t162 * t254) * t247;
t134 = 0.1e1 / t138 ^ 2;
t218 = pkin(2) ^ 2 + t258 ^ 2;
t259 = pkin(3) ^ 2;
t324 = t247 * t417;
t459 = 0.2e1 * pkin(2);
t352 = pkin(3) * t459;
t420 = t111 * t133;
t346 = t111 * t134 * (-t258 * t324 + (t232 * t259 + t253 * t352 + t218) * t420);
t240 = xDDP(3);
t393 = t133 * t240;
t141 = t197 * t235 + t274 * t233;
t177 = t233 * t364 - t235 * t254;
t383 = t233 * t234;
t449 = pkin(2) * t247;
t102 = (t177 * t223 + t226 * t376) * t445 + (-t141 * t223 + t171 * t226) * t253 + (-t223 * t383 + t226 * t236) * t449;
t426 = t102 * t133;
t99 = -(t177 * t226 - t223 * t376) * t445 + (t141 * t226 + t171 * t223) * t253 + (t223 * t236 + t226 * t383) * t449;
t438 = t133 * t99;
t334 = t258 * t420;
t295 = t247 * t334;
t87 = t295 - t417;
t69 = t253 * t346 + t120 * t393 + (pkin(2) * t323 - t253 * t87) * t325 + t242 * t438 + t241 * t426;
t60 = t187 * t223 - t188 * t226 - t234 * t69 - t201;
t429 = t60 * t254;
t51 = t168 * t248 - t429;
t293 = g(1) * t225 - g(2) * t222;
t167 = t293 * t235 - t441;
t366 = t236 * t246;
t179 = t233 * t252 + t235 * t366;
t213 = t246 * t258;
t196 = pkin(2) * t252 + t213;
t379 = t234 * t245;
t275 = pkin(3) * t379 - t193 * t236;
t119 = -t179 * t446 - t196 * t233 * t251 + (pkin(2) * t379 + t275 * t251) * t235;
t182 = t222 * t256 - t225 * t257;
t158 = t182 * t235 + t208;
t161 = t182 * t233 - t370;
t443 = pkin(3) * t251;
t210 = pkin(2) + t443;
t185 = t210 * t246 - t216;
t116 = (t210 * t252 + t213) * t158 * t236 - t161 * t185;
t190 = t210 * t367;
t374 = t234 * t251;
t149 = t185 * t374 + t190;
t418 = t116 / t149;
t326 = t260 * t418;
t328 = t130 * t418;
t110 = t158 * t374 + (t158 * t366 + t161 * t252) * t245;
t131 = 0.1e1 / t137 ^ 2;
t327 = t245 * t418;
t421 = t110 * t130;
t347 = t110 * t131 * (-t258 * t327 + (t231 * t259 + t251 * t352 + t218) * t421);
t397 = t130 * t240;
t140 = t196 * t235 + t275 * t233;
t176 = t233 * t366 - t235 * t252;
t450 = pkin(2) * t245;
t101 = (t176 * t222 + t225 * t378) * t446 + (-t140 * t222 + t170 * t225) * t251 + (-t222 * t383 + t225 * t236) * t450;
t427 = t101 * t130;
t98 = -(t176 * t225 - t222 * t378) * t446 + (t140 * t225 + t170 * t222) * t251 + (t222 * t236 + t225 * t383) * t450;
t439 = t130 * t98;
t336 = t258 * t421;
t296 = t245 * t336;
t86 = t296 - t418;
t68 = t251 * t347 + t119 * t397 + (pkin(2) * t326 - t251 * t86) * t328 + t242 * t439 + t241 * t427;
t59 = t187 * t222 - t188 * t225 - t234 * t68 - t201;
t430 = t59 * t252;
t50 = t167 * t246 - t430;
t294 = g(1) * t224 - g(2) * t221;
t166 = t294 * t235 - t441;
t368 = t236 * t244;
t178 = t233 * t250 + t235 * t368;
t212 = t244 * t258;
t195 = pkin(2) * t250 + t212;
t381 = t234 * t243;
t276 = pkin(3) * t381 - t192 * t236;
t118 = -t178 * t447 - t195 * t233 * t249 + (pkin(2) * t381 + t276 * t249) * t235;
t181 = t221 * t256 - t224 * t257;
t157 = t181 * t235 + t208;
t160 = t181 * t233 - t370;
t444 = pkin(3) * t249;
t209 = pkin(2) + t444;
t184 = t209 * t244 - t215;
t115 = (t209 * t250 + t212) * t157 * t236 - t160 * t184;
t189 = t209 * t369;
t375 = t234 * t249;
t148 = t184 * t375 + t189;
t419 = t115 / t148;
t329 = t260 * t419;
t330 = t127 * t419;
t109 = t157 * t375 + (t157 * t368 + t160 * t250) * t243;
t128 = 0.1e1 / t136 ^ 2;
t307 = t243 * t419;
t422 = t109 * t127;
t348 = t109 * t128 * (-t258 * t307 + (t230 * t259 + t249 * t352 + t218) * t422);
t401 = t127 * t240;
t139 = t195 * t235 + t276 * t233;
t175 = t233 * t368 - t235 * t250;
t451 = pkin(2) * t243;
t100 = (t175 * t221 + t224 * t380) * t447 + (-t139 * t221 + t169 * t224) * t249 + (-t221 * t383 + t224 * t236) * t451;
t428 = t100 * t127;
t97 = -(t175 * t224 - t221 * t380) * t447 + (t139 * t224 + t169 * t221) * t249 + (t221 * t236 + t224 * t383) * t451;
t440 = t127 * t97;
t338 = t258 * t422;
t297 = t243 * t338;
t85 = t297 - t419;
t67 = t249 * t348 + t118 * t401 + (pkin(2) * t329 - t249 * t85) * t330 + t242 * t440 + t241 * t428;
t58 = t187 * t221 - t188 * t224 - t234 * t67 - t201;
t431 = t58 * t250;
t49 = t166 * t244 - t431;
t285 = t329 * t422;
t152 = t175 * t243 + t233 * t375;
t154 = t178 * t243 + t235 * t375;
t339 = t250 * t422;
t357 = t244 * t249;
t34 = t152 * t401 - ((t234 * t339 + t236 * t329) * t447 + ((-t307 + t338) * t244 + pkin(2) * t339) * t375 + t236 * t85) * t127 * t422 - (t234 * t250 * t329 + (t236 * t230 - t357 * t381 - t236) * t422) * t330 - t154 * t477;
t358 = t243 * t249;
t454 = 0.2e1 * t230 - 0.1e1;
t22 = t454 * t285 + t34 * t358;
t471 = 0.2e1 * t22;
t284 = t326 * t421;
t153 = t176 * t245 + t233 * t374;
t155 = t179 * t245 + t235 * t374;
t337 = t252 * t421;
t355 = t246 * t251;
t35 = t153 * t397 - ((t234 * t337 + t236 * t326) * t446 + ((-t327 + t336) * t246 + pkin(2) * t337) * t374 + t236 * t86) * t130 * t421 - (t234 * t252 * t326 + (t236 * t231 - t355 * t379 - t236) * t421) * t328 - t155 * t476;
t356 = t245 * t251;
t453 = 0.2e1 * t231 - 0.1e1;
t23 = t453 * t284 + t35 * t356;
t470 = 0.2e1 * t23;
t283 = t323 * t420;
t354 = t247 * t253;
t151 = t177 * t247 + t233 * t373;
t156 = t180 * t247 + t235 * t373;
t335 = t254 * t420;
t353 = t248 * t253;
t36 = t151 * t393 - ((t234 * t335 + t236 * t323) * t445 + ((-t324 + t334) * t248 + pkin(2) * t335) * t373 + t236 * t87) * t133 * t420 - (t234 * t254 * t323 + (t236 * t232 - t353 * t377 - t236) * t420) * t325 - t156 * t475;
t452 = 0.2e1 * t232 - 0.1e1;
t24 = t452 * t283 + t36 * t354;
t469 = 0.2e1 * t24;
t425 = t109 ^ 2 * t128;
t424 = t110 ^ 2 * t131;
t423 = t111 ^ 2 * t134;
t455 = pkin(6) / 0.2e1;
t172 = pkin(3) * t357 + t192;
t359 = t240 * t260;
t360 = t236 * t260;
t372 = t234 * t260;
t363 = t236 * t250;
t124 = (t233 * t363 + t235 * t244) * t444 + t195 * t382 + t192 * t235;
t410 = t124 * t127;
t121 = (t233 * t244 - t235 * t363) * t444 - t195 * t371 + t192 * t233;
t413 = t121 * t260;
t448 = pkin(2) * t260;
t70 = t359 * t410 - t348 * t360 - (-t236 * t297 + (-t243 * t172 * t372 + t236 * (t249 * t448 + t230)) * t419) / (t172 * t375 + t189) * t329 + t413 * t477;
t458 = -0.2e1 * pkin(2) * t285 - 0.2e1 * t70 * t455;
t173 = pkin(3) * t355 + t193;
t362 = t236 * t252;
t125 = (t233 * t362 + t235 * t246) * t443 + t196 * t382 + t193 * t235;
t408 = t125 * t130;
t122 = (t233 * t246 - t235 * t362) * t443 - t196 * t371 + t193 * t233;
t412 = t122 * t260;
t71 = t359 * t408 - t347 * t360 - (-t236 * t296 + (-t245 * t173 * t372 + t236 * (t251 * t448 + t231)) * t418) / (t173 * t374 + t190) * t326 + t412 * t476;
t457 = -0.2e1 * pkin(2) * t284 - 0.2e1 * t71 * t455;
t174 = pkin(3) * t353 + t194;
t361 = t236 * t254;
t126 = (t233 * t361 + t235 * t248) * t442 + t197 * t382 + t194 * t235;
t406 = t126 * t133;
t123 = (t233 * t248 - t235 * t361) * t442 - t197 * t371 + t194 * t233;
t411 = t123 * t260;
t72 = t359 * t406 - t346 * t360 - (-t236 * t295 + (-t247 * t174 * t372 + t236 * (t253 * t448 + t232)) * t417) / (t174 * t373 + t191) * t323 + t411 * t475;
t456 = -0.2e1 * pkin(2) * t283 - 0.2e1 * t72 * t455;
t437 = t250 * t34;
t436 = t252 * t35;
t435 = t254 * t36;
t434 = t34 * t243;
t433 = t35 * t245;
t432 = t36 * t247;
t416 = t118 * t127;
t415 = t119 * t130;
t414 = t120 * t133;
t409 = t124 * t260;
t407 = t125 * t260;
t405 = t126 * t260;
t404 = t127 * t152;
t403 = t127 * t221;
t402 = t127 * t224;
t400 = t130 * t153;
t399 = t130 * t222;
t398 = t130 * t225;
t396 = t133 * t151;
t395 = t133 * t223;
t394 = t133 * t226;
t392 = t166 * t250;
t391 = t167 * t252;
t390 = t168 * t254;
t345 = t127 * t434;
t344 = t127 * t249 * t34;
t343 = t130 * t433;
t342 = t130 * t251 * t35;
t341 = t133 * t432;
t340 = t133 * t253 * t36;
t261 = 0.1e1 / pkin(3) ^ 2;
t333 = t115 ^ 2 / t148 ^ 2 * t261;
t332 = t116 ^ 2 / t149 ^ 2 * t261;
t331 = t117 ^ 2 / t150 ^ 2 * t261;
t322 = t121 * t403;
t321 = t121 * t402;
t319 = t122 * t399;
t318 = t122 * t398;
t316 = t123 * t395;
t315 = t123 * t394;
t313 = t154 * t403;
t312 = t154 * t402;
t311 = t155 * t399;
t310 = t155 * t398;
t309 = t156 * t395;
t308 = t156 * t394;
t306 = t121 * t345;
t305 = t121 * t344;
t304 = t122 * t343;
t303 = t122 * t342;
t302 = t123 * t341;
t301 = t123 * t340;
t300 = t127 * t358 * t425;
t299 = t130 * t356 * t424;
t298 = t133 * t354 * t423;
t291 = t244 * (t333 + t425) - t437;
t290 = t246 * (t332 + t424) - t436;
t289 = t248 * (t331 + t423) - t435;
t288 = t121 * t300;
t287 = t122 * t299;
t286 = t123 * t298;
t282 = 0.2e1 * t285;
t281 = 0.2e1 * t284;
t280 = 0.2e1 * t283;
t220 = t235 * g(3);
t163 = t294 * t233 + t220;
t64 = -g(1) * t221 - g(2) * t224 + t67;
t19 = pkin(2) * t425 - pkin(6) * t34 + (-t163 * t236 - t234 * t64) * t244 + t392;
t61 = -t163 * t234 + t236 * t64;
t13 = t19 * t243 + t249 * t61;
t270 = -pkin(6) * t333 + t34 * t459;
t7 = (t270 - t431) * t249 + t243 * t458 + t166 * t357;
t279 = t13 * t413 - t154 * t7;
t164 = t293 * t233 + t220;
t65 = -g(1) * t222 - g(2) * t225 + t68;
t20 = pkin(2) * t424 - pkin(6) * t35 + (-t164 * t236 - t234 * t65) * t246 + t391;
t62 = -t164 * t234 + t236 * t65;
t15 = t20 * t245 + t251 * t62;
t269 = -pkin(6) * t332 + t35 * t459;
t8 = (t269 - t430) * t251 + t245 * t457 + t167 * t355;
t278 = t15 * t412 - t155 * t8;
t165 = t292 * t233 + t220;
t66 = -g(1) * t223 - g(2) * t226 + t69;
t21 = pkin(2) * t423 - pkin(6) * t36 + (-t165 * t236 - t234 * t66) * t248 + t390;
t63 = -t165 * t234 + t236 * t66;
t17 = t21 * t247 + t253 * t63;
t268 = -pkin(6) * t331 + t36 * t459;
t9 = (t268 - t429) * t253 + t247 * t456 + t168 * t353;
t277 = -t156 * t9 + t17 * t411;
t10 = t249 * t458 + (-t270 - t49) * t243;
t14 = t19 * t249 - t243 * t61;
t273 = -t10 * t154 + t14 * t413;
t11 = t251 * t457 + (-t269 - t50) * t245;
t16 = t20 * t251 - t245 * t62;
t272 = -t11 * t155 + t16 * t412;
t12 = t253 * t456 + (-t268 - t51) * t247;
t18 = t21 * t253 - t247 * t63;
t271 = -t12 * t156 + t18 * t411;
t84 = t452 * t423;
t83 = t453 * t424;
t82 = t454 * t425;
t54 = t248 * t60 + t390;
t53 = t246 * t59 + t391;
t52 = t244 * t58 + t392;
t48 = -t247 * t331 + t72 * t253;
t47 = t72 * t247 + t253 * t331;
t46 = -t245 * t332 + t71 * t251;
t45 = t71 * t245 + t251 * t332;
t44 = -t243 * t333 + t70 * t249;
t43 = t70 * t243 + t249 * t333;
t39 = t248 * t72 + t254 * t280;
t38 = t246 * t71 + t252 * t281;
t37 = t244 * t70 + t250 * t282;
t33 = t246 * t35 + t252 * t424;
t32 = t246 * t424 - t436;
t31 = t248 * t36 + t254 * t423;
t30 = t244 * t34 + t250 * t425;
t29 = t248 * t423 - t435;
t28 = t244 * t425 - t437;
t27 = (t253 * t280 + t432) * t247;
t26 = (t251 * t281 + t433) * t245;
t25 = (t249 * t282 + t434) * t243;
t6 = (-t247 * t39 - t289 * t253) * t234 + t236 * t48;
t5 = (t289 * t247 - t253 * t39) * t234 - t236 * t47;
t4 = (-t245 * t38 - t290 * t251) * t234 + t236 * t46;
t3 = (t290 * t245 - t251 * t38) * t234 - t236 * t45;
t2 = (-t243 * t37 - t291 * t249) * t234 + t236 * t44;
t1 = (t291 * t243 - t249 * t37) * t234 - t236 * t43;
t40 = [t66 * t438 + t65 * t439 + t64 * t440, -t36 * t308 - t35 * t310 - t34 * t312, -t49 * t312 - t50 * t310 - t51 * t308 + (-t28 * t440 - t29 * t438 - t32 * t439) * t234, -t52 * t312 - t53 * t310 - t54 * t308 + (-t30 * t440 - t31 * t438 - t33 * t439) * t234, -t25 * t312 - t26 * t310 - t27 * t308 + (-t224 * t288 - t225 * t287 - t226 * t286) * t260, (-t315 * t84 - t318 * t83 - t321 * t82) * t260 - 0.2e1 * t22 * t312 - 0.2e1 * t23 * t310 - 0.2e1 * t24 * t308, -t43 * t312 - t45 * t310 - t47 * t308 + (t224 * t306 + t225 * t304 + t226 * t302) * t260, -t44 * t312 - t46 * t310 - t48 * t308 + (t224 * t305 + t225 * t303 + t226 * t301) * t260, (t315 * t72 + t318 * t71 + t321 * t70) * t260, (t226 * t277 + t6 * t99) * t133 + (t225 * t278 + t4 * t98) * t130 + (t2 * t97 + t224 * t279) * t127, (t226 * t271 + t5 * t99) * t133 + (t225 * t272 + t3 * t98) * t130 + (t1 * t97 + t224 * t273) * t127, t242 - g(1); t66 * t426 + t65 * t427 + t64 * t428, t36 * t309 + t35 * t311 + t34 * t313, t49 * t313 + t50 * t311 + t51 * t309 + (-t28 * t428 - t29 * t426 - t32 * t427) * t234, t52 * t313 + t53 * t311 + t54 * t309 + (-t30 * t428 - t31 * t426 - t33 * t427) * t234, t25 * t313 + t26 * t311 + t27 * t309 + (t221 * t288 + t222 * t287 + t223 * t286) * t260, (t316 * t84 + t319 * t83 + t322 * t82) * t260 + t313 * t471 + t311 * t470 + t309 * t469, t43 * t313 + t45 * t311 + t47 * t309 + (-t221 * t306 - t222 * t304 - t223 * t302) * t260, t44 * t313 + t46 * t311 + t48 * t309 + (-t221 * t305 - t222 * t303 - t223 * t301) * t260, (-t316 * t72 - t319 * t71 - t322 * t70) * t260, (t102 * t6 - t223 * t277) * t133 + (t101 * t4 - t222 * t278) * t130 + (t100 * t2 - t221 * t279) * t127, (t102 * t5 - t223 * t271) * t133 + (t101 * t3 - t222 * t272) * t130 + (t1 * t100 - t221 * t273) * t127, t241 - g(2); t66 * t414 + t65 * t415 + t64 * t416, t34 * t404 + t35 * t400 + t36 * t396, t49 * t404 + t50 * t400 + t51 * t396 + (-t28 * t416 - t29 * t414 - t32 * t415) * t234, t52 * t404 + t53 * t400 + t54 * t396 + (-t30 * t416 - t31 * t414 - t33 * t415) * t234, t25 * t404 + t26 * t400 + t27 * t396 + (-t124 * t300 - t125 * t299 - t126 * t298) * t260, (-t406 * t84 - t408 * t83 - t410 * t82) * t260 + t404 * t471 + t400 * t470 + t396 * t469, t43 * t404 + t45 * t400 + t47 * t396 + (t124 * t345 + t125 * t343 + t126 * t341) * t260, t44 * t404 + t46 * t400 + t48 * t396 + (t124 * t344 + t125 * t342 + t126 * t340) * t260, (t406 * t72 + t408 * t71 + t410 * t70) * t260, (t120 * t6 + t151 * t9 + t17 * t405) * t133 + (t119 * t4 + t15 * t407 + t153 * t8) * t130 + (t118 * t2 + t13 * t409 + t152 * t7) * t127, (t12 * t151 + t120 * t5 + t18 * t405) * t133 + (t11 * t153 + t119 * t3 + t16 * t407) * t130 + (t1 * t118 + t10 * t152 + t14 * t409) * t127, t240 - g(3);];
tauX_reg  = t40;
