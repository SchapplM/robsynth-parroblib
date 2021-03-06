% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RRRRR2G1A0
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
% tauX_reg [3x14]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:05
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3RRRRR2G1A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G1A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRRRR2G1A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRRRR2G1A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G1A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR2G1A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G1A0_invdyn_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G1A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G1A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:04:23
% EndTime: 2020-03-09 21:04:30
% DurationCPUTime: 7.30s
% Computational Cost: add. (28886->509), mult. (49659->940), div. (13599->25), fcn. (42342->60), ass. (0->456)
t298 = 1 / pkin(1);
t532 = 2 * t298;
t271 = xDDP(3);
t531 = -t271 / 0.2e1;
t284 = cos(qJ(2,3));
t275 = sin(qJ(2,3));
t238 = 0.1e1 / t275;
t293 = xDP(2);
t294 = xDP(1);
t133 = t275 * t293 + t284 * t294;
t134 = t275 * t294 - t284 * t293;
t268 = legFrame(3,3);
t214 = qJ(1,3) + t268;
t205 = sin(t214);
t208 = cos(t214);
t283 = cos(qJ(3,3));
t274 = sin(qJ(3,3));
t292 = xDP(3);
t428 = t274 * t292;
t91 = t428 + (t133 * t208 - t134 * t205) * t283;
t468 = t238 * t91;
t246 = t283 ^ 2;
t492 = pkin(2) * t246;
t248 = 0.1e1 / t283 ^ 2;
t296 = 0.1e1 / pkin(2);
t424 = t296 * t298;
t373 = t248 * t424;
t343 = t238 * t373;
t480 = t284 * pkin(1);
t481 = t283 * pkin(2);
t163 = t480 + t481;
t235 = pkin(1) * t293;
t236 = pkin(1) * t294;
t82 = -t163 * t428 + (-(t133 * t481 + t236) * t208 - (-t134 * t481 + t235) * t205) * t283;
t326 = t82 * t343;
t247 = 0.1e1 / t283;
t382 = t238 * t247 * t298;
t85 = t91 * t382;
t79 = t85 + t326;
t530 = (-t284 * t468 - t79 * t492) * t91;
t287 = cos(qJ(2,2));
t278 = sin(qJ(2,2));
t241 = 0.1e1 / t278;
t135 = t278 * t293 + t287 * t294;
t136 = t278 * t294 - t287 * t293;
t269 = legFrame(2,3);
t215 = qJ(1,2) + t269;
t206 = sin(t215);
t209 = cos(t215);
t286 = cos(qJ(3,2));
t277 = sin(qJ(3,2));
t427 = t277 * t292;
t92 = t427 + (t135 * t209 - t136 * t206) * t286;
t467 = t241 * t92;
t251 = t286 ^ 2;
t491 = pkin(2) * t251;
t253 = 0.1e1 / t286 ^ 2;
t371 = t253 * t424;
t341 = t241 * t371;
t476 = t287 * pkin(1);
t477 = t286 * pkin(2);
t166 = t476 + t477;
t83 = -t166 * t427 + (-(t135 * t477 + t236) * t209 - (-t136 * t477 + t235) * t206) * t286;
t325 = t83 * t341;
t252 = 0.1e1 / t286;
t379 = t241 * t252 * t298;
t86 = t92 * t379;
t80 = t86 + t325;
t529 = (-t287 * t467 - t80 * t491) * t92;
t290 = cos(qJ(2,1));
t281 = sin(qJ(2,1));
t244 = 0.1e1 / t281;
t137 = t281 * t293 + t290 * t294;
t138 = t281 * t294 - t290 * t293;
t270 = legFrame(1,3);
t216 = qJ(1,1) + t270;
t207 = sin(t216);
t210 = cos(t216);
t289 = cos(qJ(3,1));
t280 = sin(qJ(3,1));
t426 = t280 * t292;
t93 = t426 + (t137 * t210 - t138 * t207) * t289;
t466 = t244 * t93;
t256 = t289 ^ 2;
t490 = pkin(2) * t256;
t258 = 0.1e1 / t289 ^ 2;
t369 = t258 * t424;
t339 = t244 * t369;
t472 = t290 * pkin(1);
t473 = t289 * pkin(2);
t169 = t472 + t473;
t84 = -t169 * t426 + (-(t137 * t473 + t236) * t210 - (-t138 * t473 + t235) * t207) * t289;
t324 = t84 * t339;
t257 = 0.1e1 / t289;
t376 = t244 * t257 * t298;
t87 = t93 * t376;
t81 = t87 + t324;
t528 = (-t290 * t466 - t81 * t490) * t93;
t300 = t283 * t246;
t249 = 0.1e1 / t300;
t527 = t238 * t249;
t302 = t286 * t251;
t254 = 0.1e1 / t302;
t526 = t241 * t254;
t304 = t289 * t256;
t259 = 0.1e1 / t304;
t525 = t244 * t259;
t450 = t163 * t274;
t311 = t343 * t450;
t364 = t275 * t428;
t410 = t284 * t492;
t239 = 0.1e1 / t275 ^ 2;
t299 = 0.1e1 / pkin(1) ^ 2;
t441 = t239 * t299;
t465 = t249 * t91;
t295 = pkin(2) ^ 2;
t70 = t295 * t79 * t300;
t73 = t85 + t326 / 0.2e1;
t332 = (t70 + (0.2e1 * t73 * t410 + (-t364 + t468) * t283 * t247) * pkin(1)) * t441 * t465;
t333 = -t292 * t424 / 0.2e1;
t273 = xDDP(1);
t360 = t424 / 0.2e1;
t334 = t273 * t360;
t272 = xDDP(2);
t335 = t272 * t360;
t297 = 0.1e1 / pkin(2) ^ 2;
t390 = t297 / t246 ^ 2 * t82;
t356 = (t70 + (t79 * t410 - t364) * pkin(1)) * t390;
t211 = qJ(2,3) + t214;
t199 = qJ(3,3) + t211;
t187 = cos(t199);
t200 = -qJ(3,3) + t211;
t188 = cos(t200);
t521 = -0.2e1 * pkin(1);
t112 = t208 * t521 + (-t187 - t188) * pkin(2);
t262 = qJ(2,3) + qJ(3,3);
t217 = sin(t262);
t263 = qJ(2,3) - qJ(3,3);
t218 = sin(t263);
t414 = t217 + t218;
t157 = 0.1e1 / t414;
t456 = t112 * t157;
t181 = sin(t199);
t182 = sin(t200);
t109 = t205 * t521 + (-t181 - t182) * pkin(2);
t459 = t109 * t157;
t496 = -t296 / 0.2e1;
t425 = t292 * t296;
t374 = t247 * t425;
t338 = t284 * t374;
t495 = pkin(1) * t275;
t67 = (-t274 * t79 * t495 + t247 * t292) * t283 + pkin(1) * t338;
t320 = t332 * t496 - t356 * t441 / 0.2e1 + t311 * t531 + t67 * t333 * t527 + t335 * t459 + t334 * t456;
t430 = t272 * t298;
t193 = sin(t211);
t447 = t193 * t238;
t124 = t430 * t447;
t429 = t273 * t298;
t196 = cos(t211);
t444 = t196 * t238;
t127 = t429 * t444;
t344 = t274 * t382;
t142 = t271 * t344;
t261 = t292 ^ 2;
t365 = t261 * t424;
t160 = t365 * t527;
t398 = t82 * t441;
t55 = t79 * t247 * t398;
t323 = t124 + t127 + t142 + t160 + t55;
t381 = t248 * t441;
t329 = t381 * t530;
t37 = t323 - t329;
t19 = t37 + t320;
t229 = sin(t268);
t232 = cos(t268);
t432 = t261 * t297;
t361 = t432 / 0.2e1;
t285 = cos(qJ(1,3));
t478 = t285 * g(2);
t479 = t285 * g(1);
t276 = sin(qJ(1,3));
t483 = t276 * g(2);
t489 = g(1) * t276;
t520 = 0.2e1 * pkin(1);
t524 = t284 * (-(t478 - t489) * t232 + (t479 + t483) * t229 + t19 * t520) - 0.2e1 * ((-t479 / 0.2e1 - t483 / 0.2e1) * t232 + (-t478 / 0.2e1 + t489 / 0.2e1) * t229 + pkin(1) * (t248 * t361 + (t296 * t465 + t390 / 0.2e1) * t398)) * t275;
t449 = t166 * t277;
t310 = t341 * t449;
t363 = t278 * t427;
t409 = t287 * t491;
t242 = 0.1e1 / t278 ^ 2;
t440 = t242 * t299;
t464 = t254 * t92;
t71 = t295 * t80 * t302;
t74 = t86 + t325 / 0.2e1;
t331 = (t71 + (0.2e1 * t74 * t409 + (-t363 + t467) * t286 * t252) * pkin(1)) * t440 * t464;
t389 = t297 / t251 ^ 2 * t83;
t355 = (t71 + (t80 * t409 - t363) * pkin(1)) * t389;
t212 = qJ(2,2) + t215;
t201 = qJ(3,2) + t212;
t189 = cos(t201);
t202 = -qJ(3,2) + t212;
t190 = cos(t202);
t113 = t209 * t521 + (-t189 - t190) * pkin(2);
t264 = qJ(2,2) + qJ(3,2);
t219 = sin(t264);
t265 = qJ(2,2) - qJ(3,2);
t220 = sin(t265);
t413 = t219 + t220;
t158 = 0.1e1 / t413;
t455 = t113 * t158;
t183 = sin(t201);
t184 = sin(t202);
t110 = t206 * t521 + (-t183 - t184) * pkin(2);
t458 = t110 * t158;
t372 = t252 * t425;
t337 = t287 * t372;
t494 = pkin(1) * t278;
t68 = (-t277 * t80 * t494 + t252 * t292) * t286 + pkin(1) * t337;
t319 = t331 * t496 - t355 * t440 / 0.2e1 + t310 * t531 + t68 * t333 * t526 + t335 * t458 + t334 * t455;
t194 = sin(t212);
t446 = t194 * t241;
t125 = t430 * t446;
t197 = cos(t212);
t443 = t197 * t241;
t128 = t429 * t443;
t342 = t277 * t379;
t143 = t271 * t342;
t161 = t365 * t526;
t397 = t83 * t440;
t56 = t80 * t252 * t397;
t322 = t125 + t128 + t143 + t161 + t56;
t378 = t253 * t440;
t328 = t378 * t529;
t38 = t322 - t328;
t20 = t38 + t319;
t230 = sin(t269);
t233 = cos(t269);
t288 = cos(qJ(1,2));
t474 = t288 * g(2);
t475 = t288 * g(1);
t279 = sin(qJ(1,2));
t482 = t279 * g(2);
t488 = g(1) * t279;
t523 = t287 * (-(t474 - t488) * t233 + (t475 + t482) * t230 + t20 * t520) - 0.2e1 * ((-t475 / 0.2e1 - t482 / 0.2e1) * t233 + (-t474 / 0.2e1 + t488 / 0.2e1) * t230 + pkin(1) * (t253 * t361 + (t296 * t464 + t389 / 0.2e1) * t397)) * t278;
t448 = t169 * t280;
t309 = t339 * t448;
t362 = t281 * t426;
t408 = t290 * t490;
t245 = 0.1e1 / t281 ^ 2;
t439 = t245 * t299;
t463 = t259 * t93;
t72 = t295 * t81 * t304;
t75 = t87 + t324 / 0.2e1;
t330 = (t72 + (0.2e1 * t75 * t408 + (-t362 + t466) * t289 * t257) * pkin(1)) * t439 * t463;
t388 = t297 / t256 ^ 2 * t84;
t354 = (t72 + (t81 * t408 - t362) * pkin(1)) * t388;
t213 = qJ(2,1) + t216;
t203 = qJ(3,1) + t213;
t191 = cos(t203);
t204 = -qJ(3,1) + t213;
t192 = cos(t204);
t114 = t210 * t521 + (-t191 - t192) * pkin(2);
t266 = qJ(2,1) + qJ(3,1);
t221 = sin(t266);
t267 = qJ(2,1) - qJ(3,1);
t222 = sin(t267);
t412 = t221 + t222;
t159 = 0.1e1 / t412;
t454 = t114 * t159;
t185 = sin(t203);
t186 = sin(t204);
t111 = t207 * t521 + (-t185 - t186) * pkin(2);
t457 = t111 * t159;
t370 = t257 * t425;
t336 = t290 * t370;
t493 = pkin(1) * t281;
t69 = (-t280 * t81 * t493 + t257 * t292) * t289 + pkin(1) * t336;
t318 = t330 * t496 - t354 * t439 / 0.2e1 + t309 * t531 + t69 * t333 * t525 + t335 * t457 + t334 * t454;
t195 = sin(t213);
t445 = t195 * t244;
t126 = t430 * t445;
t198 = cos(t213);
t442 = t198 * t244;
t129 = t429 * t442;
t340 = t280 * t376;
t144 = t271 * t340;
t162 = t365 * t525;
t396 = t84 * t439;
t57 = t81 * t257 * t396;
t321 = t126 + t129 + t144 + t162 + t57;
t375 = t258 * t439;
t327 = t375 * t528;
t39 = t321 - t327;
t21 = t39 + t318;
t231 = sin(t270);
t234 = cos(t270);
t291 = cos(qJ(1,1));
t484 = g(2) * t291;
t282 = sin(qJ(1,1));
t485 = g(2) * t282;
t486 = g(1) * t291;
t487 = g(1) * t282;
t522 = t290 * (-(t484 - t487) * t234 + (t485 + t486) * t231 + t21 * t520) - 0.2e1 * ((-t486 / 0.2e1 - t485 / 0.2e1) * t234 + (-t484 / 0.2e1 + t487 / 0.2e1) * t231 + pkin(1) * (t258 * t361 + (t296 * t463 + t388 / 0.2e1) * t396)) * t281;
t518 = -pkin(1) / 0.2e1;
t517 = pkin(1) / 0.2e1;
t516 = g(1) / 0.2e1;
t515 = -g(2) / 0.2e1;
t514 = pkin(1) * t37;
t513 = pkin(1) * t38;
t512 = pkin(1) * t39;
t368 = t274 * t432;
t431 = t271 * t296;
t121 = t247 * t431 + t249 * t368;
t511 = pkin(1) * (t121 * t275 + 0.2e1 * t79 * t338);
t367 = t277 * t432;
t122 = t252 * t431 + t254 * t367;
t510 = pkin(1) * (t122 * t278 + 0.2e1 * t80 * t337);
t366 = t280 * t432;
t123 = t257 * t431 + t259 * t366;
t509 = pkin(1) * (t123 * t281 + 0.2e1 * t81 * t336);
t508 = t181 / 0.2e1;
t507 = t183 / 0.2e1;
t506 = t185 / 0.2e1;
t505 = -t188 / 0.2e1;
t504 = -t190 / 0.2e1;
t503 = -t192 / 0.2e1;
t502 = -t218 / 0.2e1;
t501 = -t220 / 0.2e1;
t500 = -t222 / 0.2e1;
t223 = cos(t262);
t499 = t223 / 0.2e1;
t225 = cos(t264);
t498 = t225 / 0.2e1;
t227 = cos(t266);
t497 = t227 / 0.2e1;
t392 = t79 * t425;
t16 = (t127 / 0.2e1 + t124 / 0.2e1 + t142 / 0.2e1 - t329 / 0.2e1 + t55 / 0.2e1 + t160 / 0.2e1 + t320) * t274 + t392;
t471 = t16 * t274;
t391 = t80 * t425;
t17 = (t128 / 0.2e1 + t125 / 0.2e1 + t143 / 0.2e1 - t328 / 0.2e1 + t56 / 0.2e1 + t161 / 0.2e1 + t319) * t277 + t391;
t470 = t17 * t277;
t387 = t81 * t425;
t18 = (t129 / 0.2e1 + t126 / 0.2e1 + t144 / 0.2e1 - t327 / 0.2e1 + t57 / 0.2e1 + t162 / 0.2e1 + t318) * t280 + t387;
t469 = t18 * t280;
t115 = t121 * t274 + t247 * t432;
t453 = t115 * t157;
t116 = t122 * t277 + t252 * t432;
t452 = t116 * t158;
t117 = t123 * t280 + t257 * t432;
t451 = t117 * t159;
t438 = t247 * t274;
t437 = t247 * t296;
t436 = t252 * t277;
t435 = t252 * t296;
t434 = t257 * t280;
t433 = t257 * t296;
t423 = g(2) * t508 + t187 * t516;
t422 = g(2) * t507 + t189 * t516;
t421 = g(2) * t506 + t191 * t516;
t420 = g(1) * t508 + t187 * t515;
t419 = g(1) * t507 + t189 * t515;
t418 = g(1) * t506 + t191 * t515;
t417 = g(1) * t196 + g(2) * t193;
t416 = g(1) * t197 + g(2) * t194;
t415 = g(1) * t198 + g(2) * t195;
t407 = t157 * t471;
t406 = t158 * t470;
t405 = t159 * t469;
t404 = t16 * t274 ^ 2 * t238;
t403 = t238 * t471;
t402 = t17 * t277 ^ 2 * t241;
t401 = t241 * t470;
t400 = t18 * t280 ^ 2 * t244;
t399 = t244 * t469;
t88 = t91 ^ 2;
t395 = t248 * t298 * t88;
t89 = t92 ^ 2;
t394 = t253 * t298 * t89;
t90 = t93 ^ 2;
t393 = t258 * t298 * t90;
t386 = t248 * t450;
t385 = t253 * t449;
t384 = t258 * t448;
t383 = t238 * t438;
t380 = t241 * t436;
t377 = t244 * t434;
t359 = g(1) * t193 - g(2) * t196;
t358 = g(1) * t194 - g(2) * t197;
t357 = g(1) * t195 - g(2) * t198;
t353 = t239 * t395;
t352 = t88 * t381;
t351 = t242 * t394;
t350 = t89 * t378;
t349 = t245 * t393;
t348 = t90 * t375;
t347 = t238 * t386;
t346 = t241 * t385;
t345 = t244 * t384;
t317 = g(2) * t505 + t182 * t516;
t316 = g(2) * t504 + t184 * t516;
t315 = g(2) * t503 + t186 * t516;
t314 = g(1) * t505 + t182 * t515;
t313 = g(1) * t504 + t184 * t515;
t312 = g(1) * t503 + t186 * t515;
t228 = cos(t267);
t226 = cos(t265);
t224 = cos(t263);
t141 = t234 * g(1) + t231 * g(2);
t140 = t233 * g(1) + t230 * g(2);
t139 = t232 * g(1) + t229 * g(2);
t132 = t231 * g(1) - t234 * g(2);
t131 = t230 * g(1) - t233 * g(2);
t130 = t229 * g(1) - t232 * g(2);
t108 = t123 * t289 - t258 * t366;
t107 = t122 * t286 - t253 * t367;
t106 = t121 * t283 - t248 * t368;
t105 = -t132 * t282 + t141 * t291;
t104 = -t131 * t279 + t140 * t288;
t103 = -t130 * t276 + t139 * t285;
t102 = t132 * t291 + t141 * t282;
t101 = t131 * t288 + t140 * t279;
t100 = t130 * t285 + t139 * t276;
t78 = t81 ^ 2;
t77 = t80 ^ 2;
t76 = t79 ^ 2;
t36 = t244 * t393 + t39 * t472 + t357;
t35 = t238 * t395 + t37 * t480 + t359;
t34 = t241 * t394 + t38 * t476 + t358;
t33 = t290 * t349 - t39 * t493 + t415;
t32 = t284 * t353 - t37 * t495 + t417;
t31 = t287 * t351 - t38 * t494 + t416;
t30 = ((t221 - t222) * t39 + (-t227 + t228) * t348) * t518 + t312 + t421;
t29 = ((t219 - t220) * t38 + (-t225 + t226) * t350) * t518 + t313 + t422;
t28 = ((t217 - t218) * t37 + (-t223 + t224) * t352) * t518 + t314 + t423;
t27 = ((t227 + t228) * t39 + t412 * t348) * t517 + t315 + t418;
t26 = ((t225 + t226) * t38 + t413 * t350) * t517 + t316 + t419;
t25 = ((t223 + t224) * t37 + t414 * t352) * t517 + t317 + t420;
t24 = (-t258 * t528 - t354) * t439 + (-t330 + ((-t259 * t292 * t69 - t271 * t384) * t244 + (t111 * t272 + t114 * t273) * t159) * t298) * t296 + t321;
t23 = (-t253 * t529 - t355) * t440 + (-t331 + ((-t254 * t292 * t68 - t271 * t385) * t241 + (t110 * t272 + t113 * t273) * t158) * t298) * t296 + t322;
t22 = (-t248 * t530 - t356) * t441 + (-t332 + ((-t249 * t292 * t67 - t271 * t386) * t238 + (t109 * t272 + t112 * t273) * t157) * t298) * t296 + t323;
t15 = (-t75 * t84 * t369 + t21 * t290) * t520 + t357;
t14 = (-t83 * t74 * t371 + t20 * t287) * t520 + t358;
t13 = (-t73 * t82 * t373 + t19 * t284) * t520 + t359;
t12 = (t75 * t290 * t324 + t281 * t21) * t521 + t415;
t11 = (t74 * t287 * t325 + t278 * t20) * t521 + t416;
t10 = (t73 * t284 * t326 + t275 * t19) * t521 + t417;
t9 = -t81 * t370 + (t24 * t280 + 0.2e1 * t387) * t289;
t8 = -t80 * t372 + (t23 * t277 + 0.2e1 * t391) * t286;
t7 = -t79 * t374 + (t22 * t274 + 0.2e1 * t392) * t283;
t6 = -t280 * t509 + t522 * t289;
t5 = -t522 * t280 - t289 * t509;
t4 = -t277 * t510 + t523 * t286;
t3 = -t523 * t277 - t286 * t510;
t2 = -t274 * t511 + t524 * t283;
t1 = -t524 * t274 - t283 * t511;
t40 = [(t37 * t444 + t38 * t443 + t39 * t442) * t298, (t100 * t444 + t101 * t443 + t102 * t442) * t298, (t103 * t444 + t104 * t443 + t105 * t442) * t298, (t22 * t444 + t23 * t443 + t24 * t442 + (t22 * t456 + t23 * t455 + t24 * t454) * t296) * t298, (t13 * t444 + t14 * t443 + t15 * t442 + (t34 * t455 + t35 * t456 + t36 * t454) * t296) * t298, (t10 * t444 + t11 * t443 + t12 * t442 + (t31 * t455 + t32 * t456 + t33 * t454) * t296) * t298, (t196 * t403 + t197 * t401 + t198 * t399 + (t112 * t407 + t113 * t406 + t114 * t405) * t296) * t532, (t7 * t444 + t8 * t443 + t9 * t442 + (t9 * t454 + t8 * t455 + t7 * t456) * t296) * t532, (t115 * t444 + t116 * t443 + t117 * t442 + (t112 * t453 + t113 * t452 + t114 * t451) * t296) * t298, (t106 * t444 + t107 * t443 + t108 * t442 + (t106 * t456 + t107 * t455 + t108 * t454) * t296) * t298, 0, (t2 * t444 + t4 * t443 + t6 * t442 + (t25 * t456 + t26 * t455 + t27 * t454) * t296) * t298, (t1 * t444 + t3 * t443 + t5 * t442 + (t28 * t456 + t29 * t455 + t30 * t454) * t296) * t298, t273 - g(1); (t37 * t447 + t38 * t446 + t39 * t445) * t298, (t100 * t447 + t101 * t446 + t102 * t445) * t298, (t103 * t447 + t104 * t446 + t105 * t445) * t298, (t22 * t447 + t23 * t446 + t24 * t445 + (t22 * t459 + t23 * t458 + t24 * t457) * t296) * t298, (t13 * t447 + t14 * t446 + t15 * t445 + (t34 * t458 + t35 * t459 + t36 * t457) * t296) * t298, (t10 * t447 + t11 * t446 + t12 * t445 + (t31 * t458 + t32 * t459 + t33 * t457) * t296) * t298, (t193 * t403 + t194 * t401 + t195 * t399 + (t109 * t407 + t110 * t406 + t111 * t405) * t296) * t532, (t7 * t447 + t8 * t446 + t9 * t445 + (t9 * t457 + t8 * t458 + t7 * t459) * t296) * t532, (t115 * t447 + t116 * t446 + t117 * t445 + (t109 * t453 + t110 * t452 + t111 * t451) * t296) * t298, (t106 * t447 + t107 * t446 + t108 * t445 + (t106 * t459 + t107 * t458 + t108 * t457) * t296) * t298, 0, (t2 * t447 + t4 * t446 + t6 * t445 + (t25 * t459 + t26 * t458 + t27 * t457) * t296) * t298, (t1 * t447 + t3 * t446 + t5 * t445 + (t28 * t459 + t29 * t458 + t30 * t457) * t296) * t298, t272 - g(2); (t37 * t383 + t377 * t39 + t38 * t380) * t298, (t100 * t383 + t101 * t380 + t102 * t377) * t298, (t103 * t383 + t104 * t380 + t105 * t377) * t298, (t22 * t383 + t23 * t380 + t24 * t377 + (-t22 * t347 - t23 * t346 - t24 * t345) * t296) * t298, (t13 * t383 + t14 * t380 + t15 * t377 + (-t34 * t346 - t345 * t36 - t347 * t35) * t296) * t298, (t10 * t383 + t11 * t380 + t12 * t377 + (-t31 * t346 - t32 * t347 - t33 * t345) * t296) * t298, (t247 * t404 + t252 * t402 + t257 * t400) * t532 + (-t274 * t76 - t277 * t77 - t280 * t78 + (-t163 * t248 * t404 - t166 * t253 * t402 - t169 * t258 * t400) * t532) * t296, (t257 * t78 * (-0.2e1 * t256 + 0.1e1) + t252 * t77 * (-0.2e1 * t251 + 0.1e1) + t247 * t76 * (-0.2e1 * t246 + 0.1e1)) * t296 + (t7 * t383 + t8 * t380 + t9 * t377 + (-t345 * t9 - t346 * t8 - t347 * t7) * t296) * t532, (t115 * t383 + t116 * t380 + t117 * t377) * t298 + (t22 * t438 + t23 * t436 + t24 * t434 + (-t115 * t347 - t116 * t346 - t117 * t345) * t298) * t296, (t106 * t383 + t107 * t380 + t108 * t377) * t298 + (t22 + t23 + t24 + (-t106 * t347 - t107 * t346 - t108 * t345) * t298) * t296, (t121 * t247 + t122 * t252 + t123 * t257) * t296, t6 * t340 - t27 * t309 + (-g(3) * t289 + (t221 / 0.2e1 + t500) * t349 + (-t228 / 0.2e1 + t497) * t512 - t315 + t418) * t433 + t4 * t342 - t26 * t310 + (-g(3) * t286 + (t219 / 0.2e1 + t501) * t351 + (-t226 / 0.2e1 + t498) * t513 - t316 + t419) * t435 + t2 * t344 - t25 * t311 + (-g(3) * t283 + (t217 / 0.2e1 + t502) * t353 + (-t224 / 0.2e1 + t499) * t514 - t317 + t420) * t437, t5 * t340 - t30 * t309 + (g(3) * t280 + (t228 / 0.2e1 + t497) * t349 + (-t221 / 0.2e1 + t500) * t512 - t312 + t421) * t433 + t3 * t342 - t29 * t310 + (g(3) * t277 + (t226 / 0.2e1 + t498) * t351 + (-t219 / 0.2e1 + t501) * t513 - t313 + t422) * t435 + t1 * t344 - t28 * t311 + (g(3) * t274 + (t224 / 0.2e1 + t499) * t353 + (-t217 / 0.2e1 + t502) * t514 - t314 + t423) * t437, t271 - g(3);];
tauX_reg  = t40;
