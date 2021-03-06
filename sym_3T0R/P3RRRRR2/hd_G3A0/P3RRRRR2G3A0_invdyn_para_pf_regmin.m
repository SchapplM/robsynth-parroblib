% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RRRRR2G3A0
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
% Datum: 2020-03-09 21:13
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3RRRRR2G3A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G3A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:12:04
% EndTime: 2020-03-09 21:12:15
% DurationCPUTime: 11.34s
% Computational Cost: add. (33316->523), mult. (83325->1052), div. (14535->22), fcn. (78822->42), ass. (0->462)
t288 = sin(qJ(1,3));
t297 = cos(qJ(1,3));
t280 = legFrame(3,2);
t250 = sin(t280);
t253 = cos(t280);
t341 = t253 * g(1) - g(2) * t250;
t163 = g(3) * t297 + t341 * t288;
t305 = xDP(2);
t306 = xDP(1);
t172 = t250 * t306 + t253 * t305;
t286 = sin(qJ(3,3));
t287 = sin(qJ(2,3));
t295 = cos(qJ(3,3));
t296 = cos(qJ(2,3));
t304 = xDP(3);
t326 = t250 * t305 - t253 * t306;
t465 = t288 * t304;
t109 = ((t326 * t288 - t297 * t304) * t287 - t296 * (t326 * t297 + t465)) * t295 + t286 * t172;
t266 = 0.1e1 / t295;
t260 = 0.1e1 / t287 ^ 2;
t311 = 0.1e1 / pkin(1) ^ 2;
t487 = t260 * t311;
t394 = t266 * t487;
t259 = 0.1e1 / t287;
t491 = t259 * t266;
t535 = pkin(2) * t295;
t310 = 0.1e1 / pkin(1);
t488 = t259 * t310;
t396 = t266 * t488;
t100 = t109 * t396;
t267 = 0.1e1 / t295 ^ 2;
t308 = 0.1e1 / pkin(2);
t457 = t308 * t310;
t385 = t267 * t457;
t350 = t259 * t385;
t181 = t287 * t288 - t297 * t296;
t182 = t297 * t287 + t296 * t288;
t546 = pkin(1) * t296;
t451 = t286 * t546;
t531 = pkin(2) * t306;
t532 = pkin(2) * t305;
t265 = t295 ^ 2;
t538 = pkin(2) * t265;
t539 = pkin(1) * t306;
t540 = pkin(1) * t305;
t97 = (-t326 * t181 + t182 * t304) * t538 + ((-t286 * t532 - t297 * t539) * t253 + (-t286 * t531 + t297 * t540) * t250 + pkin(1) * t465) * t295 - t172 * t451;
t332 = t97 * t350;
t88 = t100 + t332;
t70 = -t109 * t296 * t491 - t88 * t535;
t335 = t70 * t109 * t394;
t283 = xDDP(3);
t467 = t283 * t310;
t277 = qJ(1,3) + qJ(2,3);
t247 = sin(t277);
t509 = t247 * t259;
t344 = t467 * t509;
t496 = t253 * t286;
t512 = t181 * t295;
t151 = t250 * t512 + t496;
t284 = xDDP(2);
t139 = t151 * t284 * t396;
t502 = t250 * t286;
t152 = -t253 * t512 + t502;
t285 = xDDP(1);
t140 = t152 * t285 * t396;
t169 = t172 ^ 2;
t268 = t266 / t265;
t349 = t259 * t268 * t457;
t157 = t169 * t349;
t67 = t88 * t97 * t394;
t365 = t139 + t140 + t157 + t67;
t49 = -t335 - t344 + t365;
t46 = t49 * pkin(1);
t577 = t46 + t163;
t294 = sin(qJ(1,1));
t303 = cos(qJ(1,1));
t282 = legFrame(1,2);
t252 = sin(t282);
t255 = cos(t282);
t339 = t255 * g(1) - g(2) * t252;
t167 = g(3) * t303 + t339 * t294;
t174 = t252 * t306 + t255 * t305;
t292 = sin(qJ(3,1));
t293 = sin(qJ(2,1));
t301 = cos(qJ(3,1));
t302 = cos(qJ(2,1));
t324 = t252 * t305 - t255 * t306;
t461 = t294 * t304;
t111 = ((t324 * t294 - t303 * t304) * t293 - t302 * (t324 * t303 + t461)) * t301 + t292 * t174;
t274 = 0.1e1 / t301;
t264 = 0.1e1 / t293 ^ 2;
t477 = t264 * t311;
t386 = t274 * t477;
t263 = 0.1e1 / t293;
t481 = t263 * t274;
t533 = pkin(2) * t301;
t478 = t263 * t310;
t388 = t274 * t478;
t103 = t111 * t388;
t275 = 0.1e1 / t301 ^ 2;
t383 = t275 * t457;
t346 = t263 * t383;
t185 = t293 * t294 - t303 * t302;
t186 = t303 * t293 + t302 * t294;
t542 = pkin(1) * t302;
t449 = t292 * t542;
t273 = t301 ^ 2;
t536 = pkin(2) * t273;
t99 = (-t324 * t185 + t186 * t304) * t536 + ((-t292 * t532 - t303 * t539) * t255 + (-t292 * t531 + t303 * t540) * t252 + pkin(1) * t461) * t301 - t174 * t449;
t330 = t99 * t346;
t90 = t103 + t330;
t72 = -t111 * t302 * t481 - t90 * t533;
t333 = t72 * t111 * t386;
t279 = qJ(1,1) + qJ(2,1);
t249 = sin(t279);
t505 = t249 * t263;
t342 = t467 * t505;
t492 = t255 * t292;
t510 = t185 * t301;
t155 = t252 * t510 + t492;
t143 = t155 * t284 * t388;
t498 = t252 * t292;
t156 = -t255 * t510 + t498;
t144 = t156 * t285 * t388;
t171 = t174 ^ 2;
t276 = t274 / t273;
t345 = t263 * t276 * t457;
t159 = t171 * t345;
t69 = t90 * t99 * t386;
t363 = t143 + t144 + t159 + t69;
t51 = -t333 - t342 + t363;
t47 = t51 * pkin(1);
t576 = t47 + t167;
t291 = sin(qJ(1,2));
t300 = cos(qJ(1,2));
t281 = legFrame(2,2);
t251 = sin(t281);
t254 = cos(t281);
t340 = t254 * g(1) - g(2) * t251;
t165 = g(3) * t300 + t340 * t291;
t173 = t251 * t306 + t254 * t305;
t289 = sin(qJ(3,2));
t290 = sin(qJ(2,2));
t298 = cos(qJ(3,2));
t299 = cos(qJ(2,2));
t325 = t251 * t305 - t254 * t306;
t463 = t291 * t304;
t110 = ((t325 * t291 - t300 * t304) * t290 - t299 * (t325 * t300 + t463)) * t298 + t289 * t173;
t270 = 0.1e1 / t298;
t262 = 0.1e1 / t290 ^ 2;
t482 = t262 * t311;
t390 = t270 * t482;
t261 = 0.1e1 / t290;
t486 = t261 * t270;
t534 = pkin(2) * t298;
t483 = t261 * t310;
t392 = t270 * t483;
t102 = t110 * t392;
t271 = 0.1e1 / t298 ^ 2;
t384 = t271 * t457;
t348 = t261 * t384;
t183 = t290 * t291 - t300 * t299;
t184 = t300 * t290 + t299 * t291;
t544 = pkin(1) * t299;
t450 = t289 * t544;
t269 = t298 ^ 2;
t537 = pkin(2) * t269;
t98 = (-t325 * t183 + t184 * t304) * t537 + ((-t289 * t532 - t300 * t539) * t254 + (-t289 * t531 + t300 * t540) * t251 + pkin(1) * t463) * t298 - t173 * t450;
t331 = t98 * t348;
t89 = t102 + t331;
t71 = -t110 * t299 * t486 - t89 * t534;
t334 = t71 * t110 * t390;
t278 = qJ(1,2) + qJ(2,2);
t248 = sin(t278);
t507 = t248 * t261;
t343 = t467 * t507;
t494 = t254 * t289;
t511 = t183 * t298;
t153 = t251 * t511 + t494;
t141 = t153 * t284 * t392;
t500 = t251 * t289;
t154 = -t254 * t511 + t500;
t142 = t154 * t285 * t392;
t170 = t173 ^ 2;
t272 = t270 / t269;
t347 = t261 * t272 * t457;
t158 = t170 * t347;
t68 = t89 * t98 * t390;
t364 = t141 + t142 + t158 + t68;
t50 = -t334 - t343 + t364;
t48 = pkin(1) * t50;
t575 = t48 + t165;
t426 = t110 ^ 2 * t271 * t310;
t101 = t262 * t426;
t166 = -g(3) * t291 + t340 * t300;
t574 = t101 + t166;
t427 = t109 ^ 2 * t267 * t310;
t104 = t260 * t427;
t164 = -g(3) * t288 + t341 * t297;
t573 = t104 + t164;
t425 = t111 ^ 2 * t275 * t310;
t105 = t264 * t425;
t168 = -g(3) * t294 + t339 * t303;
t572 = t105 + t168;
t448 = t181 * t538;
t545 = pkin(1) * t297;
t130 = -t250 * t448 + (-pkin(2) * t496 + t250 * t545) * t295 - t253 * t451;
t124 = t130 * t284 * t350;
t133 = t253 * t448 + (-pkin(2) * t502 - t253 * t545) * t295 - t250 * t451;
t127 = t133 * t285 * t350;
t382 = t283 * t457;
t160 = pkin(1) * t288 + t182 * t535;
t412 = t160 * t491;
t148 = t382 * t412;
t474 = t268 * t308;
t434 = t97 * t474;
t338 = t434 * t487;
t378 = -t308 * t311 / 0.2e1;
t307 = pkin(2) ^ 2;
t466 = t286 * t287;
t400 = t172 * t466;
t460 = t295 * t296;
t524 = t265 * t88;
t556 = 0.2e1 * pkin(2);
t82 = t100 + t332 / 0.2e1;
t564 = t267 * (t307 * t524 + (t82 * t460 * t556 + (t109 * t259 - t400) * t266) * pkin(1));
t475 = t266 * t308;
t353 = t172 * t296 * t475;
t567 = t172 * ((-pkin(1) * t88 * t466 + t172 * t266) * t295 + pkin(1) * t353);
t61 = pkin(2) * t524 + (-t400 * t475 + t88 * t460) * pkin(1);
t323 = t260 * t109 * t378 * t564 - t61 * t338 / 0.2e1 + t124 / 0.2e1 + t127 / 0.2e1 + t148 / 0.2e1 - t349 * t567 / 0.2e1;
t309 = 0.1e1 / pkin(2) ^ 2;
t515 = t169 * t309;
t456 = t309 * t311;
t527 = 0.2e1 * t109 * t338 + t97 ^ 2 * t260 / t265 ^ 2 * t456;
t557 = 0.2e1 * pkin(1);
t571 = (-(t49 + t323) * t557 - t163) * t296 + t287 * (pkin(1) * (t267 * t515 + t527) - t164);
t447 = t183 * t537;
t543 = pkin(1) * t300;
t131 = -t251 * t447 + (-pkin(2) * t494 + t251 * t543) * t298 - t254 * t450;
t125 = t131 * t284 * t348;
t134 = t254 * t447 + (-pkin(2) * t500 - t254 * t543) * t298 - t251 * t450;
t128 = t134 * t285 * t348;
t161 = pkin(1) * t291 + t184 * t534;
t411 = t161 * t486;
t149 = t382 * t411;
t471 = t272 * t308;
t431 = t98 * t471;
t337 = t431 * t482;
t464 = t289 * t290;
t399 = t173 * t464;
t459 = t298 * t299;
t522 = t269 * t89;
t83 = t102 + t331 / 0.2e1;
t563 = t271 * (t307 * t522 + (t83 * t459 * t556 + (t110 * t261 - t399) * t270) * pkin(1));
t472 = t270 * t308;
t352 = t173 * t299 * t472;
t566 = t173 * ((-pkin(1) * t89 * t464 + t173 * t270) * t298 + pkin(1) * t352);
t62 = pkin(2) * t522 + (-t399 * t472 + t89 * t459) * pkin(1);
t322 = t262 * t110 * t378 * t563 - t62 * t337 / 0.2e1 + t125 / 0.2e1 + t128 / 0.2e1 + t149 / 0.2e1 - t347 * t566 / 0.2e1;
t514 = t170 * t309;
t526 = 0.2e1 * t110 * t337 + t98 ^ 2 * t262 / t269 ^ 2 * t456;
t570 = (-(t50 + t322) * t557 - t165) * t299 + t290 * (pkin(1) * (t271 * t514 + t526) - t166);
t446 = t185 * t536;
t541 = pkin(1) * t303;
t132 = -t252 * t446 + (-pkin(2) * t492 + t252 * t541) * t301 - t255 * t449;
t126 = t132 * t284 * t346;
t135 = t255 * t446 + (-pkin(2) * t498 - t255 * t541) * t301 - t252 * t449;
t129 = t135 * t285 * t346;
t162 = pkin(1) * t294 + t186 * t533;
t410 = t162 * t481;
t150 = t382 * t410;
t468 = t276 * t308;
t429 = t99 * t468;
t336 = t429 * t477;
t462 = t292 * t293;
t398 = t174 * t462;
t458 = t301 * t302;
t520 = t273 * t90;
t84 = t103 + t330 / 0.2e1;
t562 = t275 * (t307 * t520 + (t84 * t458 * t556 + (t111 * t263 - t398) * t274) * pkin(1));
t469 = t274 * t308;
t351 = t174 * t302 * t469;
t565 = t174 * ((-pkin(1) * t90 * t462 + t174 * t274) * t301 + pkin(1) * t351);
t63 = pkin(2) * t520 + (-t398 * t469 + t90 * t458) * pkin(1);
t321 = t264 * t111 * t378 * t562 - t63 * t336 / 0.2e1 + t126 / 0.2e1 + t129 / 0.2e1 + t150 / 0.2e1 - t345 * t565 / 0.2e1;
t513 = t171 * t309;
t525 = 0.2e1 * t111 * t336 + t99 ^ 2 * t264 / t273 ^ 2 * t456;
t569 = (-(t51 + t321) * t557 - t167) * t302 + t293 * (pkin(1) * (t275 * t513 + t525) - t168);
t568 = 0.2e1 * t310;
t558 = -0.2e1 * pkin(1);
t554 = g(1) / 0.2e1;
t553 = -g(2) / 0.2e1;
t233 = -t280 + t277;
t552 = sin(t233) / 0.2e1;
t235 = -t281 + t278;
t551 = sin(t235) / 0.2e1;
t237 = -t282 + t279;
t550 = sin(t237) / 0.2e1;
t232 = t280 + t277;
t549 = cos(t232) / 0.2e1;
t234 = t281 + t278;
t548 = cos(t234) / 0.2e1;
t236 = t282 + t279;
t547 = cos(t236) / 0.2e1;
t523 = t266 * t70;
t521 = t270 * t71;
t519 = t274 * t72;
t508 = t247 * t283;
t506 = t248 * t283;
t504 = t249 * t283;
t503 = t250 * t266;
t501 = t251 * t270;
t499 = t252 * t274;
t497 = t253 * t266;
t495 = t254 * t270;
t493 = t255 * t274;
t490 = t259 * t267;
t489 = t259 * t286;
t485 = t261 * t271;
t484 = t261 * t289;
t480 = t263 * t275;
t479 = t263 * t292;
t476 = t266 * t286;
t473 = t270 * t289;
t470 = t274 * t292;
t403 = t286 * t515;
t136 = t268 * t403 + (t250 * t285 + t253 * t284) * t475;
t455 = (t88 * t353 + t136 * t287 / 0.2e1) * t558;
t402 = t289 * t514;
t137 = t272 * t402 + (t251 * t285 + t254 * t284) * t472;
t454 = (t89 * t352 + t137 * t290 / 0.2e1) * t558;
t401 = t292 * t513;
t138 = t276 * t401 + (t252 * t285 + t255 * t284) * t469;
t453 = (t90 * t351 + t138 * t293 / 0.2e1) * t558;
t445 = t172 * t308 * t88;
t444 = t173 * t308 * t89;
t369 = t474 * t567;
t374 = t61 * t434;
t381 = t124 + t127 + t148;
t435 = t308 * t564;
t24 = (-t369 - t508) * t488 + (-t374 + (-t435 - t523) * t109) * t487 + t365 + t381;
t443 = t24 * t476;
t85 = t88 ^ 2;
t442 = (-0.2e1 * t265 + 0.1e1) * t266 * t85;
t86 = t89 ^ 2;
t441 = (-0.2e1 * t269 + 0.1e1) * t270 * t86;
t87 = t90 ^ 2;
t440 = (-0.2e1 * t273 + 0.1e1) * t274 * t87;
t368 = t471 * t566;
t373 = t62 * t431;
t380 = t125 + t128 + t149;
t432 = t308 * t563;
t25 = (-t368 - t506) * t483 + (-t373 + (-t432 - t521) * t110) * t482 + t364 + t380;
t439 = t25 * t473;
t40 = t573 * t287 + t577 * t296;
t438 = t40 * t491;
t41 = t574 * t290 + t575 * t299;
t437 = t41 * t486;
t42 = t572 * t293 + t576 * t302;
t436 = t42 * t481;
t367 = t468 * t565;
t372 = t63 * t429;
t379 = t126 + t129 + t150;
t430 = t308 * t562;
t27 = (-t367 - t504) * t478 + (-t372 + (-t430 - t519) * t111) * t477 + t363 + t379;
t433 = t27 * t470;
t428 = t90 * t174 * t308;
t424 = t130 * t490;
t423 = t131 * t485;
t422 = t132 * t480;
t421 = t133 * t490;
t420 = t134 * t485;
t419 = t135 * t480;
t418 = t151 * t491;
t417 = t152 * t491;
t416 = t153 * t486;
t415 = t154 * t486;
t414 = t155 * t481;
t413 = t156 * t481;
t409 = t163 * t491;
t408 = t164 * t491;
t407 = t165 * t486;
t406 = t166 * t486;
t405 = t167 * t481;
t404 = t168 * t481;
t397 = t259 * t476;
t395 = t267 * t489;
t393 = t261 * t473;
t391 = t271 * t484;
t389 = t263 * t470;
t387 = t275 * t479;
t16 = (t140 / 0.2e1 + t139 / 0.2e1 - t344 / 0.2e1 - t335 / 0.2e1 + t67 / 0.2e1 + t157 / 0.2e1 + t323) * t286 + t445;
t371 = t16 * t397;
t17 = (t142 / 0.2e1 + t141 / 0.2e1 - t343 / 0.2e1 - t334 / 0.2e1 + t68 / 0.2e1 + t158 / 0.2e1 + t322) * t289 + t444;
t370 = t17 * t393;
t18 = (t144 / 0.2e1 + t143 / 0.2e1 - t342 / 0.2e1 - t333 / 0.2e1 + t69 / 0.2e1 + t159 / 0.2e1 + t321) * t292 + t428;
t366 = t18 * t389;
t362 = t130 * t395;
t361 = t131 * t391;
t360 = t132 * t387;
t359 = t133 * t395;
t358 = t134 * t391;
t357 = t135 * t387;
t356 = t160 * t397;
t355 = t161 * t393;
t354 = t162 * t389;
t217 = sin(t232);
t224 = cos(t233);
t329 = g(1) * t552 + g(2) * t549 + t217 * t554 + t224 * t553 + g(3) * cos(t277);
t219 = sin(t234);
t226 = cos(t235);
t328 = g(1) * t551 + g(2) * t548 + t219 * t554 + t226 * t553 + g(3) * cos(t278);
t221 = sin(t236);
t228 = cos(t237);
t327 = g(1) * t550 + g(2) * t547 + t221 * t554 + t228 * t553 + g(3) * cos(t279);
t320 = g(1) * t549 + g(2) * t552 - g(3) * t247 + t217 * t553 + t224 * t554;
t319 = g(1) * t548 + g(2) * t551 - g(3) * t248 + t219 * t553 + t226 * t554;
t318 = g(1) * t547 + g(2) * t550 - g(3) * t249 + t221 * t553 + t228 * t554;
t192 = g(1) * t252 + g(2) * t255;
t190 = g(1) * t251 + g(2) * t254;
t188 = g(1) * t250 + g(2) * t253;
t117 = t138 * t292 + t274 * t513;
t116 = t137 * t289 + t270 * t514;
t115 = t136 * t286 + t266 * t515;
t114 = t138 * t301 - t275 * t401;
t113 = t137 * t298 - t271 * t402;
t112 = t136 * t295 - t267 * t403;
t45 = t263 * t425 + t51 * t542 + t327;
t44 = t261 * t426 + t50 * t544 + t328;
t43 = t259 * t427 + t49 * t546 + t329;
t39 = t302 * t105 - t293 * t47 + t318;
t38 = t299 * t101 - t290 * t48 + t319;
t37 = t296 * t104 - t287 * t46 + t320;
t36 = -t293 * t576 + t572 * t302;
t35 = -t287 * t577 + t573 * t296;
t34 = -t290 * t575 + t574 * t299;
t33 = t192 * t292 + t301 * t36;
t32 = -t192 * t301 + t292 * t36;
t31 = t188 * t286 + t295 * t35;
t30 = -t188 * t295 + t286 * t35;
t29 = t190 * t289 + t298 * t34;
t28 = -t190 * t298 + t289 * t34;
t26 = 0.2e1 * t143 + 0.2e1 * t144 + 0.2e1 * t159 + 0.2e1 * t69 + (-t367 - 0.2e1 * t504) * t478 + (-t372 + (-t430 - 0.2e1 * t519) * t111) * t477 + t379;
t23 = 0.2e1 * t139 + 0.2e1 * t140 + 0.2e1 * t157 + 0.2e1 * t67 + (-t369 - 0.2e1 * t508) * t488 + (-t374 + (-t435 - 0.2e1 * t523) * t109) * t487 + t381;
t22 = 0.2e1 * t141 + 0.2e1 * t142 + 0.2e1 * t158 + 0.2e1 * t68 + (-t368 - 0.2e1 * t506) * t483 + (-t373 + (-t432 - 0.2e1 * t521) * t110) * t482 + t380;
t15 = -pkin(1) * (t293 * t26 + t525 * t302) + t318;
t14 = -pkin(1) * (t290 * t22 + t526 * t299) + t319;
t13 = -pkin(1) * (t287 * t23 + t527 * t296) + t320;
t12 = (-0.2e1 * t84 * t99 * t383 + t26 * t302) * pkin(1) + t327;
t11 = (-0.2e1 * t83 * t98 * t384 + t22 * t299) * pkin(1) + t328;
t10 = (-0.2e1 * t82 * t97 * t385 + t23 * t296) * pkin(1) + t329;
t9 = t27 * t292 * t301 + (-t274 + 0.2e1 * t301) * t428;
t8 = t25 * t289 * t298 + (-t270 + 0.2e1 * t298) * t444;
t7 = t24 * t286 * t295 + (-t266 + 0.2e1 * t295) * t445;
t6 = t569 * t292 + t301 * t453;
t5 = t292 * t453 - t569 * t301;
t4 = t570 * t289 + t298 * t454;
t3 = t289 * t454 - t570 * t298;
t2 = t571 * t286 + t295 * t455;
t1 = t286 * t455 - t571 * t295;
t19 = [(t51 * t413 + t50 * t415 + t49 * t417) * t310, (t152 * t409 + t154 * t407 + t156 * t405) * t310, (t152 * t408 + t154 * t406 + t156 * t404) * t310, (t24 * t417 + t25 * t415 + t27 * t413 + (t24 * t421 + t25 * t420 + t27 * t419) * t308) * t310, (t10 * t417 + t11 * t415 + t12 * t413 + (t419 * t45 + t420 * t44 + t421 * t43) * t308) * t310, (t13 * t417 + t14 * t415 + t15 * t413 + (t37 * t421 + t38 * t420 + t39 * t419) * t308) * t310, (t152 * t371 + t154 * t370 + t156 * t366) * t568 + (-t85 * t502 - t86 * t500 - t87 * t498 + (t16 * t359 + t17 * t358 + t18 * t357) * t568) * t308, (t250 * t442 + t251 * t441 + t252 * t440) * t308 + (t7 * t417 + t8 * t415 + t9 * t413 + (t419 * t9 + t420 * t8 + t421 * t7) * t308) * t568, (t115 * t417 + t116 * t415 + t117 * t413) * t310 + (t250 * t443 + t251 * t439 + t252 * t433 + (t115 * t421 + t116 * t420 + t117 * t419) * t310) * t308, (t112 * t417 + t113 * t415 + t114 * t413) * t310 + (t24 * t250 + t25 * t251 + t252 * t27 + (t112 * t421 + t113 * t420 + t114 * t419) * t310) * t308, (t136 * t503 + t137 * t501 + t138 * t499) * t308, (t1 * t417 + t3 * t415 + t413 * t5) * t310 + (t30 * t503 + t28 * t501 + t32 * t499 + (t133 * t438 + t134 * t437 + t135 * t436) * t310) * t308, (t2 * t417 + t4 * t415 + t413 * t6) * t310 + (t31 * t503 + t29 * t501 + t33 * t499 + (-t357 * t42 - t358 * t41 - t359 * t40) * t310) * t308, t285 - g(1); (t51 * t414 + t50 * t416 + t49 * t418) * t310, (t151 * t409 + t153 * t407 + t155 * t405) * t310, (t151 * t408 + t153 * t406 + t155 * t404) * t310, (t24 * t418 + t25 * t416 + t27 * t414 + (t24 * t424 + t25 * t423 + t27 * t422) * t308) * t310, (t10 * t418 + t11 * t416 + t12 * t414 + (t422 * t45 + t423 * t44 + t424 * t43) * t308) * t310, (t13 * t418 + t14 * t416 + t15 * t414 + (t37 * t424 + t38 * t423 + t39 * t422) * t308) * t310, (t151 * t371 + t153 * t370 + t155 * t366) * t568 + (-t85 * t496 - t86 * t494 - t87 * t492 + (t16 * t362 + t17 * t361 + t18 * t360) * t568) * t308, (t253 * t442 + t254 * t441 + t255 * t440) * t308 + (t7 * t418 + t8 * t416 + t9 * t414 + (t422 * t9 + t423 * t8 + t424 * t7) * t308) * t568, (t115 * t418 + t116 * t416 + t117 * t414) * t310 + (t253 * t443 + t254 * t439 + t255 * t433 + (t115 * t424 + t116 * t423 + t117 * t422) * t310) * t308, (t112 * t418 + t113 * t416 + t114 * t414) * t310 + (t24 * t253 + t25 * t254 + t255 * t27 + (t112 * t424 + t113 * t423 + t114 * t422) * t310) * t308, (t136 * t497 + t137 * t495 + t138 * t493) * t308, (t1 * t418 + t3 * t416 + t414 * t5) * t310 + (t30 * t497 + t28 * t495 + t32 * t493 + (t130 * t438 + t131 * t437 + t132 * t436) * t310) * t308, (t2 * t418 + t4 * t416 + t414 * t6) * t310 + (t31 * t497 + t29 * t495 + t33 * t493 + (-t360 * t42 - t361 * t41 - t362 * t40) * t310) * t308, t284 - g(2); (-t49 * t509 - t50 * t507 - t51 * t505) * t310, (-t163 * t509 - t165 * t507 - t167 * t505) * t310, (-t164 * t509 - t166 * t507 - t168 * t505) * t310, (-t24 * t509 - t25 * t507 - t27 * t505 + (t24 * t412 + t25 * t411 + t27 * t410) * t308) * t310, (-t10 * t509 - t11 * t507 - t12 * t505 + (t410 * t45 + t411 * t44 + t412 * t43) * t308) * t310, (-t13 * t509 - t14 * t507 - t15 * t505 + (t37 * t412 + t38 * t411 + t39 * t410) * t308) * t310, (-t16 * t247 * t489 - t17 * t248 * t484 - t18 * t249 * t479 + (t16 * t356 + t17 * t355 + t18 * t354) * t308) * t568, (-t7 * t509 - t8 * t507 - t9 * t505 + (t410 * t9 + t411 * t8 + t412 * t7) * t308) * t568, (-t115 * t509 - t116 * t507 - t117 * t505 + (t115 * t412 + t116 * t411 + t117 * t410) * t308) * t310, (-t112 * t509 - t113 * t507 - t114 * t505 + (t112 * t412 + t113 * t411 + t114 * t410) * t308) * t310, 0, (-t1 * t509 - t3 * t507 - t5 * t505 + (t160 * t259 * t40 + t161 * t261 * t41 + t162 * t263 * t42) * t308) * t310, (-t2 * t509 - t4 * t507 - t6 * t505 + (-t354 * t42 - t355 * t41 - t356 * t40) * t308) * t310, t283 - g(3);];
tauX_reg  = t19;
