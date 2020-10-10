% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PRRRR8V2G2A0
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
% Datum: 2020-08-06 17:50
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3PRRRR8V2G2A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_regmin: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:49:43
% EndTime: 2020-08-06 17:49:56
% DurationCPUTime: 12.53s
% Computational Cost: add. (66759->478), mult. (135321->923), div. (5832->17), fcn. (125931->22), ass. (0->399)
t250 = cos(qJ(2,3));
t258 = pkin(7) + pkin(6);
t215 = t250 * t258;
t244 = sin(qJ(2,3));
t192 = pkin(2) * t244 - t215;
t234 = sin(pkin(4));
t236 = cos(pkin(4));
t243 = sin(qJ(3,3));
t378 = t236 * t243;
t169 = pkin(3) * t378 + t192 * t234;
t249 = cos(qJ(3,3));
t392 = t234 * t244;
t230 = t249 ^ 2;
t453 = pkin(3) * t230;
t136 = pkin(2) * t378 + t169 * t249 + t392 * t453;
t127 = 0.1e1 / t136;
t252 = cos(qJ(2,2));
t216 = t252 * t258;
t246 = sin(qJ(2,2));
t193 = pkin(2) * t246 - t216;
t245 = sin(qJ(3,2));
t376 = t236 * t245;
t170 = pkin(3) * t376 + t193 * t234;
t251 = cos(qJ(3,2));
t390 = t234 * t246;
t231 = t251 ^ 2;
t452 = pkin(3) * t231;
t137 = pkin(2) * t376 + t170 * t251 + t390 * t452;
t130 = 0.1e1 / t137;
t254 = cos(qJ(2,1));
t217 = t254 * t258;
t248 = sin(qJ(2,1));
t194 = pkin(2) * t248 - t217;
t247 = sin(qJ(3,1));
t374 = t236 * t247;
t171 = pkin(3) * t374 + t194 * t234;
t253 = cos(qJ(3,1));
t388 = t234 * t248;
t232 = t253 ^ 2;
t451 = pkin(3) * t232;
t138 = pkin(2) * t374 + t171 * t253 + t388 * t451;
t133 = 0.1e1 / t138;
t239 = legFrame(1,2);
t223 = sin(t239);
t226 = cos(t239);
t256 = xDP(2);
t257 = xDP(1);
t183 = t223 * t256 - t226 * t257;
t233 = sin(pkin(8));
t255 = xDP(3);
t208 = t255 * t233;
t235 = cos(pkin(8));
t159 = t183 * t235 + t208;
t379 = t235 * t255;
t162 = t183 * t233 - t379;
t448 = pkin(3) * t253;
t211 = pkin(2) + t448;
t186 = t211 * t248 - t217;
t214 = t248 * t258;
t117 = t162 * (t211 * t254 + t214) * t236 + t159 * t186;
t191 = t211 * t374;
t385 = t234 * t253;
t150 = t186 * t385 + t191;
t261 = 0.1e1 / pkin(3) ^ 2;
t340 = t117 ^ 2 * t261 / t150 ^ 2;
t373 = t236 * t248;
t177 = t233 * t373 - t235 * t254;
t153 = t177 * t247 + t233 * t385;
t180 = t233 * t254 + t235 * t373;
t380 = t235 * t253;
t156 = -t180 * t247 - t234 * t380;
t260 = 0.1e1 / pkin(3);
t426 = t117 / t150;
t331 = t260 * t426;
t332 = t247 * t426;
t333 = t133 * t426;
t111 = t162 * t385 - (t159 * t254 - t162 * t373) * t247;
t429 = t111 * t133;
t343 = t258 * t429;
t344 = t254 * t429;
t362 = t248 * t253;
t389 = t234 * t247;
t240 = xDDP(3);
t402 = t133 * t240;
t241 = xDDP(2);
t242 = xDDP(1);
t466 = (t223 * t241 - t226 * t242) * t133;
t304 = t247 * t343;
t87 = t304 - t426;
t36 = t156 * t402 - ((t234 * t344 + t236 * t331) * t451 + ((-t332 + t343) * t248 + pkin(2) * t344) * t385 + t236 * t87) * t133 * t429 - (t254 * t234 * t331 + (t236 * t232 - t362 * t389 - t236) * t429) * t333 + t153 * t466;
t465 = 0.2e1 * pkin(2);
t383 = t235 * t236;
t187 = -t234 * g(1) + g(2) * t383;
t188 = g(1) * t383 + g(2) * t234;
t395 = t233 * t236;
t201 = g(3) * t395;
t197 = pkin(2) * t254 + t214;
t283 = pkin(3) * t389 - t194 * t236;
t120 = -t177 * t451 + t197 * t380 + (pkin(2) * t389 + t283 * t253) * t233;
t134 = 0.1e1 / t138 ^ 2;
t218 = pkin(2) ^ 2 + t258 ^ 2;
t259 = pkin(3) ^ 2;
t361 = pkin(3) * t465;
t355 = t111 * t134 * (-t258 * t332 + (t232 * t259 + t253 * t361 + t218) * t429);
t141 = -t197 * t233 + t283 * t235;
t394 = t234 * t235;
t455 = pkin(2) * t247;
t102 = (t180 * t226 + t223 * t388) * t451 + (-t141 * t226 + t171 * t223) * t253 + (t223 * t236 - t226 * t394) * t455;
t435 = t102 * t133;
t101 = -(t180 * t223 - t226 * t388) * t451 + (t141 * t223 + t171 * t226) * t253 + (t223 * t394 + t226 * t236) * t455;
t436 = t101 * t133;
t69 = t253 * t355 + t120 * t402 + (pkin(2) * t331 - t87 * t253) * t333 + t241 * t436 + t242 * t435;
t271 = t187 * t223 - t188 * t226 + t234 * t69 + t201;
t57 = t271 * t254;
t489 = -pkin(6) * t340 + t36 * t465 + t57;
t238 = legFrame(2,2);
t222 = sin(t238);
t225 = cos(t238);
t182 = t222 * t256 - t225 * t257;
t158 = t182 * t235 + t208;
t161 = t182 * t233 - t379;
t449 = pkin(3) * t251;
t210 = pkin(2) + t449;
t185 = t210 * t246 - t216;
t213 = t246 * t258;
t116 = t161 * (t210 * t252 + t213) * t236 + t158 * t185;
t190 = t210 * t376;
t386 = t234 * t251;
t149 = t185 * t386 + t190;
t341 = t116 ^ 2 * t261 / t149 ^ 2;
t375 = t236 * t246;
t176 = t233 * t375 - t235 * t252;
t152 = t176 * t245 + t233 * t386;
t179 = t233 * t252 + t235 * t375;
t381 = t235 * t251;
t155 = -t179 * t245 - t234 * t381;
t427 = t116 / t149;
t334 = t260 * t427;
t335 = t245 * t427;
t336 = t130 * t427;
t110 = t161 * t386 - (t158 * t252 - t161 * t375) * t245;
t430 = t110 * t130;
t345 = t258 * t430;
t346 = t252 * t430;
t364 = t246 * t251;
t391 = t234 * t245;
t406 = t130 * t240;
t467 = t130 * (t222 * t241 - t225 * t242);
t305 = t245 * t345;
t86 = t305 - t427;
t35 = t155 * t406 - ((t234 * t346 + t236 * t334) * t452 + ((-t335 + t345) * t246 + pkin(2) * t346) * t386 + t236 * t86) * t130 * t430 - (t252 * t234 * t334 + (t236 * t231 - t364 * t391 - t236) * t430) * t336 + t152 * t467;
t196 = pkin(2) * t252 + t213;
t284 = pkin(3) * t391 - t193 * t236;
t119 = -t176 * t452 + t196 * t381 + (pkin(2) * t391 + t284 * t251) * t233;
t131 = 0.1e1 / t137 ^ 2;
t356 = t110 * t131 * (-t258 * t335 + (t231 * t259 + t251 * t361 + t218) * t430);
t140 = -t196 * t233 + t284 * t235;
t456 = pkin(2) * t245;
t100 = (t179 * t225 + t222 * t390) * t452 + (-t140 * t225 + t170 * t222) * t251 + (t222 * t236 - t225 * t394) * t456;
t437 = t100 * t130;
t99 = -(t179 * t222 - t225 * t390) * t452 + (t140 * t222 + t170 * t225) * t251 + (t222 * t394 + t225 * t236) * t456;
t444 = t130 * t99;
t68 = t251 * t356 + t119 * t406 + (pkin(2) * t334 - t86 * t251) * t336 + t241 * t444 + t242 * t437;
t272 = t187 * t222 - t188 * t225 + t234 * t68 + t201;
t56 = t272 * t252;
t488 = -pkin(6) * t341 + t35 * t465 + t56;
t377 = t236 * t244;
t175 = t233 * t377 - t235 * t250;
t387 = t234 * t249;
t151 = t175 * t243 + t233 * t387;
t178 = t233 * t250 + t235 * t377;
t382 = t235 * t249;
t154 = -t178 * t243 - t234 * t382;
t237 = legFrame(3,2);
t221 = sin(t237);
t224 = cos(t237);
t181 = t221 * t256 - t224 * t257;
t157 = t181 * t235 + t208;
t160 = t181 * t233 - t379;
t450 = pkin(3) * t249;
t209 = pkin(2) + t450;
t184 = t209 * t244 - t215;
t212 = t244 * t258;
t115 = t160 * (t209 * t250 + t212) * t236 + t157 * t184;
t189 = t209 * t378;
t148 = t184 * t387 + t189;
t428 = t115 / t148;
t337 = t260 * t428;
t338 = t243 * t428;
t339 = t127 * t428;
t109 = t160 * t387 - (t157 * t250 - t160 * t377) * t243;
t431 = t109 * t127;
t347 = t258 * t431;
t348 = t250 * t431;
t366 = t244 * t249;
t393 = t234 * t243;
t410 = t127 * t240;
t468 = t127 * (t221 * t241 - t224 * t242);
t306 = t243 * t347;
t85 = t306 - t428;
t34 = t154 * t410 - ((t234 * t348 + t236 * t337) * t453 + ((-t338 + t347) * t244 + pkin(2) * t348) * t387 + t236 * t85) * t127 * t431 - (t250 * t234 * t337 + (t236 * t230 - t366 * t393 - t236) * t431) * t339 + t151 * t468;
t342 = t115 ^ 2 * t261 / t148 ^ 2;
t195 = pkin(2) * t250 + t212;
t285 = pkin(3) * t393 - t192 * t236;
t118 = -t175 * t453 + t195 * t382 + (pkin(2) * t393 + t285 * t249) * t233;
t128 = 0.1e1 / t136 ^ 2;
t357 = t109 * t128 * (-t258 * t338 + (t230 * t259 + t249 * t361 + t218) * t431);
t139 = -t195 * t233 + t285 * t235;
t457 = pkin(2) * t243;
t98 = (t178 * t224 + t221 * t392) * t453 + (-t139 * t224 + t169 * t221) * t249 + (t221 * t236 - t224 * t394) * t457;
t445 = t127 * t98;
t97 = -(t178 * t221 - t224 * t392) * t453 + (t139 * t221 + t169 * t224) * t249 + (t221 * t394 + t224 * t236) * t457;
t446 = t127 * t97;
t67 = t249 * t357 + t118 * t410 + (pkin(2) * t337 - t85 * t249) * t339 + t241 * t446 + t242 * t445;
t273 = t187 * t221 - t188 * t224 + t234 * t67 + t201;
t55 = t273 * t250;
t487 = -pkin(6) * t342 + t34 * t465 + t55;
t220 = g(3) * t235;
t303 = g(1) * t224 - g(2) * t221;
t469 = t303 * t233 + t220;
t483 = t244 * t469;
t302 = g(1) * t225 - g(2) * t222;
t470 = t302 * t233 + t220;
t482 = t246 * t470;
t301 = g(1) * t226 - g(2) * t223;
t471 = t301 * t233 + t220;
t481 = t248 * t471;
t480 = t469 * t250;
t479 = t470 * t252;
t478 = t471 * t254;
t294 = t337 * t431;
t367 = t243 * t249;
t460 = 0.2e1 * t230 - 0.1e1;
t22 = t460 * t294 + t34 * t367;
t477 = 0.2e1 * t22;
t293 = t334 * t430;
t365 = t245 * t251;
t459 = 0.2e1 * t231 - 0.1e1;
t23 = t459 * t293 + t35 * t365;
t476 = 0.2e1 * t23;
t292 = t331 * t429;
t363 = t247 * t253;
t458 = 0.2e1 * t232 - 0.1e1;
t24 = t458 * t292 + t36 * t363;
t475 = 0.2e1 * t24;
t434 = t109 ^ 2 * t128;
t433 = t110 ^ 2 * t131;
t432 = t111 ^ 2 * t134;
t461 = pkin(6) / 0.2e1;
t172 = pkin(3) * t366 + t192;
t368 = t240 * t260;
t369 = t236 * t260;
t384 = t234 * t260;
t372 = t236 * t250;
t124 = (t233 * t372 + t235 * t244) * t450 + t195 * t395 + t192 * t235;
t416 = t124 * t260;
t121 = (t233 * t244 - t235 * t372) * t450 - t195 * t383 + t233 * t192;
t422 = t121 * t127;
t454 = pkin(2) * t260;
t70 = t368 * t422 - t357 * t369 - (-t236 * t306 + (-t172 * t243 * t384 + t236 * (t249 * t454 + t230)) * t428) / (t172 * t387 + t189) * t337 + t416 * t468;
t464 = -0.2e1 * pkin(2) * t294 - 0.2e1 * t70 * t461;
t173 = pkin(3) * t364 + t193;
t371 = t236 * t252;
t125 = (t233 * t371 + t235 * t246) * t449 + t196 * t395 + t193 * t235;
t415 = t125 * t260;
t122 = (t233 * t246 - t235 * t371) * t449 - t196 * t383 + t233 * t193;
t420 = t122 * t130;
t71 = t368 * t420 - t356 * t369 - (-t236 * t305 + (-t173 * t245 * t384 + t236 * (t251 * t454 + t231)) * t427) / (t173 * t386 + t190) * t334 + t415 * t467;
t463 = -0.2e1 * pkin(2) * t293 - 0.2e1 * t71 * t461;
t174 = pkin(3) * t362 + t194;
t370 = t236 * t254;
t126 = (t233 * t370 + t235 * t248) * t448 + t197 * t395 + t194 * t235;
t414 = t126 * t260;
t123 = (t233 * t248 - t235 * t370) * t448 - t197 * t383 + t233 * t194;
t418 = t123 * t133;
t72 = t368 * t418 - t355 * t369 - (-t236 * t304 + (-t174 * t247 * t384 + t236 * (t253 * t454 + t232)) * t426) / (t174 * t385 + t191) * t331 + t414 * t466;
t462 = -0.2e1 * pkin(2) * t292 - 0.2e1 * t72 * t461;
t447 = g(3) * t233;
t443 = t250 * t34;
t442 = t252 * t35;
t441 = t254 * t36;
t440 = t34 * t243;
t439 = t35 * t245;
t438 = t36 * t247;
t425 = t118 * t127;
t424 = t119 * t130;
t423 = t120 * t133;
t421 = t121 * t260;
t419 = t122 * t260;
t417 = t123 * t260;
t413 = t127 * t154;
t412 = t127 * t221;
t411 = t127 * t224;
t409 = t130 * t155;
t408 = t130 * t222;
t407 = t130 * t225;
t405 = t133 * t156;
t404 = t133 * t223;
t403 = t133 * t226;
t354 = t127 * t440;
t353 = t127 * t249 * t34;
t352 = t130 * t439;
t351 = t130 * t251 * t35;
t350 = t133 * t438;
t349 = t133 * t253 * t36;
t330 = t124 * t412;
t329 = t124 * t411;
t327 = t125 * t408;
t326 = t125 * t407;
t324 = t126 * t404;
t323 = t126 * t403;
t321 = t151 * t412;
t320 = t151 * t411;
t319 = t152 * t408;
t318 = t152 * t407;
t317 = t153 * t404;
t316 = t153 * t403;
t315 = t124 * t354;
t314 = t124 * t353;
t313 = t125 * t352;
t312 = t125 * t351;
t311 = t126 * t350;
t310 = t126 * t349;
t309 = t127 * t367 * t434;
t308 = t130 * t365 * t433;
t307 = t133 * t363 * t432;
t300 = t244 * (t342 + t434) - t443;
t299 = t246 * (t341 + t433) - t442;
t298 = t248 * (t340 + t432) - t441;
t297 = t124 * t309;
t296 = t125 * t308;
t295 = t126 * t307;
t291 = 0.2e1 * t294;
t290 = 0.2e1 * t293;
t289 = 0.2e1 * t292;
t163 = t303 * t235 - t447;
t64 = -g(1) * t221 - g(2) * t224 + t67;
t19 = pkin(2) * t434 - pkin(6) * t34 + (t163 * t236 - t234 * t64) * t244 + t480;
t61 = t163 * t234 + t236 * t64;
t14 = t19 * t249 - t243 * t61;
t7 = t249 * t464 + (-t483 - t487) * t243;
t288 = t14 * t416 + t151 * t7;
t164 = t302 * t235 - t447;
t65 = -g(1) * t222 - g(2) * t225 + t68;
t20 = pkin(2) * t433 - pkin(6) * t35 + (t164 * t236 - t234 * t65) * t246 + t479;
t62 = t164 * t234 + t236 * t65;
t16 = t20 * t251 - t245 * t62;
t8 = t251 * t463 + (-t482 - t488) * t245;
t287 = t152 * t8 + t16 * t415;
t165 = t301 * t235 - t447;
t66 = -g(1) * t223 - g(2) * t226 + t69;
t21 = pkin(2) * t432 - pkin(6) * t36 + (t165 * t236 - t234 * t66) * t248 + t478;
t63 = t165 * t234 + t236 * t66;
t18 = t21 * t253 - t247 * t63;
t9 = t253 * t462 + (-t481 - t489) * t247;
t286 = t153 * t9 + t18 * t414;
t10 = t243 * t464 + t249 * t487 + t469 * t366;
t13 = t19 * t243 + t249 * t61;
t282 = t10 * t151 + t13 * t416;
t11 = t245 * t463 + t488 * t251 + t470 * t364;
t15 = t20 * t245 + t251 * t62;
t281 = t11 * t152 + t15 * t415;
t12 = t247 * t462 + t489 * t253 + t471 * t362;
t17 = t21 * t247 + t253 * t63;
t280 = t12 * t153 + t17 * t414;
t84 = t458 * t432;
t83 = t459 * t433;
t82 = t460 * t434;
t54 = t57 + t481;
t53 = t56 + t482;
t52 = t55 + t483;
t51 = -t248 * t271 + t478;
t50 = -t246 * t272 + t479;
t49 = -t244 * t273 + t480;
t48 = -t247 * t340 + t72 * t253;
t47 = t72 * t247 + t253 * t340;
t46 = -t245 * t341 + t71 * t251;
t45 = t71 * t245 + t251 * t341;
t44 = -t243 * t342 + t70 * t249;
t43 = t70 * t243 + t249 * t342;
t39 = t248 * t72 + t254 * t289;
t38 = t246 * t71 + t252 * t290;
t37 = t244 * t70 + t250 * t291;
t33 = t246 * t35 + t252 * t433;
t32 = t246 * t433 - t442;
t31 = t248 * t36 + t254 * t432;
t30 = t244 * t34 + t250 * t434;
t29 = t248 * t432 - t441;
t28 = t244 * t434 - t443;
t27 = (t253 * t289 + t438) * t247;
t26 = (t251 * t290 + t439) * t245;
t25 = (t249 * t291 + t440) * t243;
t6 = (-t247 * t39 - t298 * t253) * t234 + t236 * t48;
t5 = (t298 * t247 - t253 * t39) * t234 - t236 * t47;
t4 = (-t245 * t38 - t299 * t251) * t234 + t236 * t46;
t3 = (t299 * t245 - t251 * t38) * t234 - t236 * t45;
t2 = (-t243 * t37 - t300 * t249) * t234 + t236 * t44;
t1 = (t300 * t243 - t249 * t37) * t234 - t236 * t43;
t40 = [t66 * t435 + t65 * t437 + t64 * t445, -t36 * t316 - t35 * t318 - t34 * t320, -t52 * t320 - t53 * t318 - t54 * t316 + (-t28 * t445 - t29 * t435 - t32 * t437) * t234, -t49 * t320 - t50 * t318 - t51 * t316 + (-t30 * t445 - t31 * t435 - t33 * t437) * t234, -t25 * t320 - t26 * t318 - t27 * t316 + (t224 * t297 + t225 * t296 + t226 * t295) * t260, (t323 * t84 + t326 * t83 + t329 * t82) * t260 - 0.2e1 * t22 * t320 - 0.2e1 * t23 * t318 - 0.2e1 * t24 * t316, -t43 * t320 - t45 * t318 - t47 * t316 + (-t224 * t315 - t225 * t313 - t226 * t311) * t260, -t44 * t320 - t46 * t318 - t48 * t316 + (-t224 * t314 - t225 * t312 - t226 * t310) * t260, (-t323 * t72 - t326 * t71 - t329 * t70) * t260, (t102 * t6 - t226 * t280) * t133 + (t100 * t4 - t225 * t281) * t130 + (t2 * t98 - t224 * t282) * t127, (t102 * t5 - t226 * t286) * t133 + (t100 * t3 - t225 * t287) * t130 + (t1 * t98 - t224 * t288) * t127, t242 - g(1); t66 * t436 + t65 * t444 + t64 * t446, t36 * t317 + t35 * t319 + t34 * t321, t52 * t321 + t53 * t319 + t54 * t317 + (-t28 * t446 - t29 * t436 - t32 * t444) * t234, t49 * t321 + t50 * t319 + t51 * t317 + (-t30 * t446 - t31 * t436 - t33 * t444) * t234, t25 * t321 + t26 * t319 + t27 * t317 + (-t221 * t297 - t222 * t296 - t223 * t295) * t260, (-t324 * t84 - t327 * t83 - t330 * t82) * t260 + t321 * t477 + t319 * t476 + t317 * t475, t43 * t321 + t45 * t319 + t47 * t317 + (t221 * t315 + t222 * t313 + t223 * t311) * t260, t44 * t321 + t46 * t319 + t48 * t317 + (t221 * t314 + t222 * t312 + t223 * t310) * t260, (t324 * t72 + t327 * t71 + t330 * t70) * t260, (t101 * t6 + t223 * t280) * t133 + (t222 * t281 + t4 * t99) * t130 + (t2 * t97 + t221 * t282) * t127, (t101 * t5 + t223 * t286) * t133 + (t222 * t287 + t3 * t99) * t130 + (t1 * t97 + t221 * t288) * t127, t241 - g(2); t66 * t423 + t65 * t424 + t64 * t425, t34 * t413 + t35 * t409 + t36 * t405, t52 * t413 + t53 * t409 + t54 * t405 + (-t28 * t425 - t29 * t423 - t32 * t424) * t234, t49 * t413 + t50 * t409 + t51 * t405 + (-t30 * t425 - t31 * t423 - t33 * t424) * t234, t25 * t413 + t26 * t409 + t27 * t405 + (-t121 * t309 - t122 * t308 - t123 * t307) * t260, (-t84 * t418 - t83 * t420 - t82 * t422) * t260 + t413 * t477 + t409 * t476 + t405 * t475, t43 * t413 + t45 * t409 + t47 * t405 + (t121 * t354 + t122 * t352 + t123 * t350) * t260, t44 * t413 + t46 * t409 + t48 * t405 + (t121 * t353 + t122 * t351 + t123 * t349) * t260, (t72 * t418 + t71 * t420 + t70 * t422) * t260, (t12 * t156 + t120 * t6 + t17 * t417) * t133 + (t11 * t155 + t119 * t4 + t15 * t419) * t130 + (t10 * t154 + t118 * t2 + t13 * t421) * t127, (t120 * t5 + t156 * t9 + t18 * t417) * t133 + (t119 * t3 + t155 * t8 + t16 * t419) * t130 + (t1 * t118 + t14 * t421 + t154 * t7) * t127, t240 - g(3);];
tauX_reg  = t40;
