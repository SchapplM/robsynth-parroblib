% Calculate minimal parameter regressor of inverse dynamics forces for
% P4PRRRR8V1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% xDP [4x1]
%   Generalized platform velocities
% xDDP [4x1]
%   Generalized platform accelerations
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tauX_reg [4x15]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:06
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P4PRRRR8V1G2A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V1G2A0_invdyn_para_pf_regmin: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR8V1G2A0_invdyn_para_pf_regmin: xDP has to be [4x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [4 1]), ...
  'P4PRRRR8V1G2A0_invdyn_para_pf_regmin: xDDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V1G2A0_invdyn_para_pf_regmin: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR8V1G2A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P4PRRRR8V1G2A0_invdyn_para_pf_regmin: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V1G2A0_invdyn_para_pf_regmin: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V1G2A0_invdyn_para_pf_regmin: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:03:27
% EndTime: 2020-08-07 11:03:47
% DurationCPUTime: 20.97s
% Computational Cost: add. (119109->668), mult. (280584->1292), div. (11524->18), fcn. (219964->30), ass. (0->526)
t301 = cos(pkin(6));
t279 = t301 * g(3);
t299 = sin(pkin(6));
t307 = legFrame(4,2);
t280 = sin(t307);
t284 = cos(t307);
t404 = g(1) * t284 - g(2) * t280;
t221 = t404 * t299 + t279;
t304 = sin(qJ(2,4));
t306 = cos(qJ(2,4));
t300 = sin(pkin(3));
t302 = cos(pkin(3));
t536 = t301 * t302;
t251 = -t300 * g(1) + g(2) * t536;
t252 = g(1) * t536 + t300 * g(2);
t602 = t299 * g(3);
t264 = t302 * t602;
t305 = cos(qJ(3,4));
t609 = pkin(2) * t305;
t254 = pkin(5) * t304 + t306 * t609;
t522 = t304 * t305;
t253 = pkin(2) * t522 - t306 * pkin(5);
t303 = sin(qJ(3,4));
t610 = pkin(2) * t303;
t363 = -t253 * t302 + t300 * t610;
t143 = -t299 * t254 + t301 * t363;
t448 = t253 * t300 + t302 * t610;
t123 = t143 * t280 + t284 * t448;
t124 = -t143 * t284 + t280 * t448;
t144 = t301 * t254 + t299 * t363;
t331 = xP(4);
t288 = sin(t331);
t289 = cos(t331);
t332 = koppelP(4,2);
t336 = koppelP(4,1);
t243 = t288 * t336 + t289 * t332;
t247 = -t288 * t332 + t289 * t336;
t327 = xDP(4);
t298 = t327 ^ 2;
t311 = xDDP(4);
t313 = xDDP(2);
t171 = -t298 * t243 + t247 * t311 + t313;
t314 = xDDP(1);
t175 = -t243 * t311 - t298 * t247 + t314;
t290 = 0.1e1 / t305;
t312 = xDDP(3);
t329 = xDP(2);
t330 = xDP(1);
t141 = (-t243 * t327 + t330) * t284 - t280 * (t247 * t327 + t329);
t328 = xDP(3);
t269 = t328 * t301;
t137 = t141 * t299 + t269;
t527 = t302 * t328;
t535 = t302 * t304;
t540 = t300 * t305;
t116 = ((-t141 * t535 - t306 * t328) * t299 + (t141 * t306 - t304 * t527) * t301) * t303 - t137 * t540;
t611 = 0.1e1 / t448;
t578 = t116 * t611;
t497 = t611 * t578;
t541 = t299 * t328;
t133 = t141 * t301 - t541;
t534 = t302 * t306;
t104 = -(t304 * t133 + t137 * t534) * t609 + pkin(5) * (t133 * t306 - t137 * t535);
t191 = 0.1e1 / t448 ^ 2;
t586 = t104 * t191;
t340 = pkin(5) ^ 2;
t341 = pkin(2) ^ 2;
t554 = t611 * t303;
t469 = t290 * t554;
t87 = -pkin(5) * t104 * t469 + (t290 * t340 + t305 * t341) * t578;
t512 = pkin(5) * t578;
t444 = t303 * t512;
t587 = t104 * t611;
t91 = (t444 - t587) * t290;
t72 = (t123 * t171 + t124 * t175 + t144 * t312) * t611 + (t87 * t497 - t91 * t586) * t290;
t351 = t251 * t280 - t252 * t284 + t72 * t300 + t264;
t55 = t304 * t221 + t306 * t351;
t326 = cos(qJ(2,1));
t320 = sin(qJ(2,1));
t325 = cos(qJ(3,1));
t513 = t320 * t325;
t257 = pkin(2) * t513 - t326 * pkin(5);
t319 = sin(qJ(3,1));
t606 = pkin(2) * t319;
t445 = t257 * t300 + t302 * t606;
t201 = 0.1e1 / t445 ^ 2;
t614 = 0.1e1 / t445;
t324 = cos(qJ(2,2));
t318 = sin(qJ(2,2));
t323 = cos(qJ(3,2));
t516 = t318 * t323;
t256 = pkin(2) * t516 - t324 * pkin(5);
t317 = sin(qJ(3,2));
t607 = pkin(2) * t317;
t446 = t256 * t300 + t302 * t607;
t199 = 0.1e1 / t446 ^ 2;
t613 = 0.1e1 / t446;
t322 = cos(qJ(2,3));
t316 = sin(qJ(2,3));
t321 = cos(qJ(3,3));
t519 = t316 * t321;
t255 = pkin(2) * t519 - t322 * pkin(5);
t315 = sin(qJ(3,3));
t608 = pkin(2) * t315;
t447 = t255 * t300 + t302 * t608;
t197 = 0.1e1 / t447 ^ 2;
t612 = 0.1e1 / t447;
t347 = t325 ^ 2;
t297 = 0.1e1 / t347;
t335 = koppelP(1,2);
t339 = koppelP(1,1);
t246 = t288 * t339 + t289 * t335;
t250 = -t288 * t335 + t289 * t339;
t310 = legFrame(1,2);
t283 = sin(t310);
t287 = cos(t310);
t140 = (-t246 * t327 + t330) * t287 - t283 * (t250 * t327 + t329);
t136 = t140 * t299 + t269;
t531 = t302 * t320;
t537 = t300 * t325;
t122 = ((-t140 * t531 - t326 * t328) * t299 + (t140 * t326 - t320 * t527) * t301) * t319 - t136 * t537;
t575 = t122 ^ 2 * t201;
t494 = t297 * t575;
t545 = t614 * t319;
t372 = t494 * t545;
t238 = t299 * t531 - t301 * t326;
t528 = t302 * t326;
t603 = pkin(2) * t325;
t449 = t238 * pkin(5) + (t299 * t528 + t301 * t320) * t603;
t356 = t449 * t372;
t346 = t323 ^ 2;
t295 = 0.1e1 / t346;
t334 = koppelP(2,2);
t338 = koppelP(2,1);
t245 = t288 * t338 + t289 * t334;
t249 = -t288 * t334 + t289 * t338;
t309 = legFrame(2,2);
t282 = sin(t309);
t286 = cos(t309);
t139 = (-t245 * t327 + t330) * t286 - t282 * (t249 * t327 + t329);
t135 = t139 * t299 + t269;
t532 = t302 * t318;
t538 = t300 * t323;
t121 = ((-t139 * t532 - t324 * t328) * t299 + (t139 * t324 - t318 * t527) * t301) * t317 - t135 * t538;
t576 = t121 ^ 2 * t199;
t495 = t295 * t576;
t548 = t613 * t317;
t373 = t495 * t548;
t529 = t302 * t324;
t604 = pkin(2) * t323;
t450 = -pkin(5) * (-t299 * t532 + t301 * t324) + (t299 * t529 + t301 * t318) * t604;
t357 = t450 * t373;
t345 = t321 ^ 2;
t293 = 0.1e1 / t345;
t333 = koppelP(3,2);
t337 = koppelP(3,1);
t244 = t288 * t337 + t289 * t333;
t248 = -t288 * t333 + t289 * t337;
t308 = legFrame(3,2);
t281 = sin(t308);
t285 = cos(t308);
t142 = (-t244 * t327 + t330) * t285 - t281 * (t248 * t327 + t329);
t138 = t142 * t299 + t269;
t533 = t302 * t316;
t539 = t300 * t321;
t120 = ((-t142 * t533 - t322 * t328) * t299 + (t142 * t322 - t316 * t527) * t301) * t315 - t138 * t539;
t577 = t120 ^ 2 * t197;
t496 = t293 * t577;
t551 = t612 * t315;
t374 = t496 * t551;
t530 = t302 * t322;
t605 = pkin(2) * t321;
t451 = -pkin(5) * (-t299 * t533 + t301 * t322) + (t299 * t530 + t301 * t316) * t605;
t358 = t451 * t374;
t344 = t305 ^ 2;
t291 = 0.1e1 / t344;
t579 = t116 ^ 2 * t191;
t498 = t291 * t579;
t375 = t498 * t554;
t237 = -t299 * t535 + t301 * t306;
t452 = -pkin(5) * t237 + (t299 * t534 + t301 * t304) * t609;
t359 = t452 * t375;
t601 = t313 - g(2);
t600 = t314 - g(1);
t525 = t303 * t304;
t355 = t302 * t525 + t540;
t524 = t303 * t306;
t179 = t299 * t355 - t301 * t524;
t342 = 0.1e1 / pkin(2);
t553 = t611 * t342;
t468 = t290 * t553;
t382 = t104 * t300 * t468;
t555 = t611 * t290;
t470 = t284 * t555;
t412 = t179 * t470;
t543 = t280 * t290;
t471 = t611 * t543;
t180 = -t299 * t524 - t301 * t355;
t478 = t180 * t555;
t526 = t302 * t342;
t33 = t179 * t171 * t471 - t175 * t412 + t312 * t478 + (-(t302 * t91 + (pkin(2) * (t300 * t306 * t578 + t526 * t587) * t344 - t300 * (t104 * t554 - t512) * t522) * t290) * t497 - (t306 * t382 + (-t300 * t525 + (-t290 + t305) * t302) * t578) * t586) * t291;
t599 = t611 * t33;
t134 = t142 * t301 - t541;
t108 = -(t316 * t134 + t138 * t530) * t605 + pkin(5) * (t322 * t134 - t138 * t533);
t172 = -t298 * t244 + t248 * t311 + t313;
t176 = -t244 * t311 - t298 * t248 + t314;
t292 = 0.1e1 / t321;
t550 = t612 * t342;
t463 = t292 * t550;
t380 = t108 * t300 * t463;
t521 = t315 * t316;
t354 = t302 * t521 + t539;
t520 = t315 * t322;
t181 = t299 * t354 - t301 * t520;
t552 = t612 * t292;
t465 = t285 * t552;
t410 = t181 * t465;
t466 = t281 * t552;
t411 = t181 * t466;
t184 = -t299 * t520 - t301 * t354;
t477 = t184 * t552;
t574 = t120 * t612;
t493 = t612 * t574;
t511 = pkin(5) * t574;
t584 = t108 * t197;
t585 = t108 * t612;
t443 = t315 * t511;
t93 = (t443 - t585) * t292;
t34 = t172 * t411 - t176 * t410 + t312 * t477 + (-(t302 * t93 + (pkin(2) * (t300 * t322 * t574 + t526 * t585) * t345 - t300 * (t108 * t551 - t511) * t519) * t292) * t493 + (-t322 * t380 + (t300 * t521 + (t292 - t321) * t302) * t574) * t584) * t293;
t598 = t612 * t34;
t131 = t139 * t301 - t541;
t109 = -(t318 * t131 + t135 * t529) * t604 + pkin(5) * (t324 * t131 - t135 * t532);
t173 = -t298 * t245 + t249 * t311 + t313;
t177 = -t245 * t311 - t298 * t249 + t314;
t294 = 0.1e1 / t323;
t547 = t613 * t342;
t458 = t294 * t547;
t378 = t109 * t300 * t458;
t518 = t317 * t318;
t353 = t302 * t518 + t538;
t517 = t317 * t324;
t182 = t299 * t353 - t301 * t517;
t549 = t613 * t294;
t460 = t286 * t549;
t408 = t182 * t460;
t461 = t282 * t549;
t409 = t182 * t461;
t185 = -t299 * t517 - t301 * t353;
t476 = t185 * t549;
t573 = t121 * t613;
t492 = t613 * t573;
t510 = pkin(5) * t573;
t582 = t109 * t199;
t583 = t109 * t613;
t442 = t317 * t510;
t94 = (t442 - t583) * t294;
t35 = t173 * t409 - t177 * t408 + t312 * t476 + (-(t302 * t94 + (pkin(2) * (t300 * t324 * t573 + t526 * t583) * t346 - t300 * (t109 * t548 - t510) * t516) * t294) * t492 + (-t324 * t378 + (t300 * t518 + (t294 - t323) * t302) * t573) * t582) * t295;
t597 = t613 * t35;
t132 = t140 * t301 - t541;
t110 = -(t320 * t132 + t136 * t528) * t603 + pkin(5) * (t326 * t132 - t136 * t531);
t174 = -t298 * t246 + t250 * t311 + t313;
t178 = -t246 * t311 - t298 * t250 + t314;
t515 = t319 * t320;
t352 = t302 * t515 + t537;
t514 = t319 * t326;
t183 = t299 * t352 - t301 * t514;
t296 = 0.1e1 / t325;
t544 = t614 * t342;
t453 = t296 * t544;
t376 = t110 * t300 * t453;
t546 = t614 * t296;
t455 = t287 * t546;
t407 = t183 * t455;
t542 = t283 * t296;
t456 = t614 * t542;
t186 = -t299 * t514 - t301 * t352;
t475 = t186 * t546;
t572 = t122 * t614;
t491 = t614 * t572;
t509 = pkin(5) * t572;
t580 = t110 * t201;
t581 = t110 * t614;
t441 = t319 * t509;
t95 = (t441 - t581) * t296;
t36 = t183 * t174 * t456 - t178 * t407 + t312 * t475 + (-(t302 * t95 + (pkin(2) * (t300 * t326 * t572 + t526 * t581) * t347 - t300 * (t110 * t545 - t509) * t513) * t296) * t491 + (-t326 * t376 + (t300 * t515 + (t296 - t325) * t302) * t572) * t580) * t297;
t596 = t614 * t36;
t595 = t306 * t33;
t594 = t322 * t34;
t593 = t324 * t35;
t592 = t326 * t36;
t591 = t33 * t303;
t590 = t34 * t315;
t589 = t35 * t317;
t588 = t36 * t319;
t571 = t123 * t611;
t570 = t124 * t611;
t258 = pkin(5) * t316 + t322 * t605;
t362 = -t255 * t302 + t300 * t608;
t145 = -t299 * t258 + t301 * t362;
t125 = t145 * t281 + t285 * t447;
t569 = t125 * t612;
t126 = -t145 * t285 + t281 * t447;
t568 = t126 * t612;
t259 = pkin(5) * t318 + t324 * t604;
t361 = -t256 * t302 + t300 * t607;
t146 = -t299 * t259 + t301 * t361;
t127 = t146 * t282 + t286 * t446;
t567 = t127 * t613;
t128 = -t146 * t286 + t282 * t446;
t566 = t128 * t613;
t260 = pkin(5) * t320 + t326 * t603;
t360 = -t257 * t302 + t300 * t606;
t147 = -t299 * t260 + t301 * t360;
t129 = t147 * t283 + t287 * t445;
t565 = t129 * t614;
t130 = -t147 * t287 + t283 * t445;
t564 = t130 * t614;
t563 = t144 * t611;
t148 = t301 * t258 + t299 * t362;
t562 = t148 * t612;
t149 = t301 * t259 + t299 * t361;
t561 = t149 * t613;
t150 = t301 * t260 + t299 * t360;
t560 = t150 * t614;
t559 = t179 * t611;
t558 = t181 * t612;
t557 = t182 * t613;
t556 = t183 * t614;
t508 = t452 * t599;
t507 = t451 * t598;
t506 = t450 * t597;
t505 = t449 * t596;
t403 = g(1) * t285 - g(2) * t281;
t222 = t403 * t299 + t279;
t464 = t292 * t551;
t88 = -pkin(5) * t108 * t464 + (t292 * t340 + t321 * t341) * t574;
t76 = (t125 * t172 + t126 * t176 + t148 * t312) * t612 + (t88 * t493 - t93 * t584) * t292;
t350 = t251 * t281 - t252 * t285 + t76 * t300 + t264;
t60 = t316 * t222 + t322 * t350;
t504 = t60 * t558;
t402 = g(1) * t286 - g(2) * t282;
t223 = t402 * t299 + t279;
t459 = t294 * t548;
t89 = -pkin(5) * t109 * t459 + (t294 * t340 + t323 * t341) * t573;
t77 = (t127 * t173 + t128 * t177 + t149 * t312) * t613 + (t89 * t492 - t94 * t582) * t294;
t349 = t251 * t282 - t252 * t286 + t77 * t300 + t264;
t61 = t318 * t223 + t324 * t349;
t503 = t61 * t557;
t343 = 0.1e1 / pkin(2) ^ 2;
t502 = t104 ^ 2 * t191 * t343;
t501 = t108 ^ 2 * t197 * t343;
t500 = t109 ^ 2 * t199 * t343;
t499 = t110 ^ 2 * t201 * t343;
t371 = t243 * t284 + t280 * t247;
t490 = t371 * t559;
t368 = t246 * t287 + t283 * t250;
t489 = t368 * t556;
t370 = t244 * t285 + t281 * t248;
t488 = t370 * t558;
t369 = t245 * t286 + t282 * t249;
t487 = t369 * t557;
t155 = -(-t299 * t304 + t301 * t534) * t609 - pkin(5) * (t299 * t306 + t301 * t535);
t486 = t155 * t555;
t485 = t452 * t371 * t611;
t158 = -(-t299 * t316 + t301 * t530) * t605 - pkin(5) * (t299 * t322 + t301 * t533);
t484 = t158 * t552;
t159 = -(-t299 * t318 + t301 * t529) * t604 - pkin(5) * (t299 * t324 + t301 * t532);
t483 = t159 * t549;
t160 = -(-t299 * t320 + t301 * t528) * t603 - pkin(5) * (t299 * t326 + t301 * t531);
t482 = t160 * t546;
t481 = t451 * t370 * t612;
t480 = t450 * t369 * t613;
t479 = t449 * t368 * t614;
t187 = -t237 * t303 + t299 * t540;
t474 = t187 * t543;
t188 = t238 * t319 + t299 * t537;
t473 = t188 * t542;
t472 = t611 * t553;
t467 = t612 * t550;
t462 = t613 * t547;
t457 = t614 * t544;
t454 = t296 * t545;
t440 = t60 * t477;
t439 = t61 * t476;
t401 = g(1) * t287 - g(2) * t283;
t224 = t401 * t299 + t279;
t90 = -pkin(5) * t110 * t454 + (t296 * t340 + t325 * t341) * t572;
t78 = (t129 * t174 + t130 * t178 + t150 * t312) * t614 + (t90 * t491 - t95 * t580) * t296;
t348 = t251 * t283 - t252 * t287 + t78 * t300 + t264;
t62 = t320 * t224 + t326 * t348;
t438 = t62 * t475;
t437 = t33 * t469;
t436 = t55 * t469;
t435 = t34 * t464;
t434 = t35 * t459;
t433 = t36 * t454;
t432 = t116 * t472;
t431 = t120 * t467;
t430 = t121 * t462;
t429 = t122 * t457;
t428 = t290 * t490;
t427 = t296 * t489;
t426 = t292 * t488;
t425 = t294 * t487;
t424 = t452 * t471;
t423 = t452 * t470;
t422 = t290 * t485;
t421 = t451 * t466;
t420 = t451 * t465;
t419 = t292 * t481;
t418 = t450 * t461;
t417 = t450 * t460;
t416 = t294 * t480;
t415 = t449 * t456;
t414 = t449 * t455;
t413 = t296 * t479;
t406 = t611 * t474;
t405 = t614 * t473;
t400 = t452 * t437;
t399 = t451 * t435;
t398 = t450 * t434;
t397 = t449 * t433;
t396 = t179 * t436;
t395 = t60 * t426;
t394 = t61 * t425;
t393 = t62 * t427;
t392 = t60 * t411;
t391 = t61 * t409;
t390 = t60 * t410;
t389 = t61 * t408;
t388 = t62 * t407;
t387 = (t502 + t579) * t291 * t304 - t595;
t386 = (t501 + t577) * t293 * t316 - t594;
t385 = (t500 + t576) * t295 * t318 - t593;
t384 = (t499 + t575) * t297 * t320 - t592;
t383 = t104 * t432;
t381 = t108 * t431;
t379 = t109 * t430;
t377 = t110 * t429;
t367 = 0.2e1 * t383;
t366 = 0.2e1 * t381;
t365 = 0.2e1 * t379;
t364 = 0.2e1 * t377;
t242 = -t288 * t311 - t289 * t298;
t241 = -t288 * t298 + t289 * t311;
t220 = t401 * t301 - t602;
t219 = t402 * t301 - t602;
t218 = t403 * t301 - t602;
t217 = t404 * t301 - t602;
t194 = t326 * t224;
t193 = t324 * t223;
t192 = t322 * t222;
t189 = t306 * t221;
t114 = (t129 * t250 - t130 * t246) * t614;
t113 = (t127 * t249 - t128 * t245) * t613;
t112 = (t125 * t248 - t126 * t244) * t612;
t111 = (t123 * t247 - t124 * t243) * t611;
t102 = (t297 - 0.2e1) * t575;
t101 = (t295 - 0.2e1) * t576;
t100 = (t293 - 0.2e1) * t577;
t99 = (t291 - 0.2e1) * t579;
t75 = -t283 * g(1) - t287 * g(2) + t78;
t74 = -t282 * g(1) - t286 * g(2) + t77;
t73 = -t281 * g(1) - t285 * g(2) + t76;
t71 = -t280 * g(1) - t284 * g(2) + t72;
t70 = -t300 * t220 - t75 * t302;
t69 = -t300 * t219 - t74 * t302;
t68 = -t300 * t218 - t73 * t302;
t67 = -t300 * t217 - t71 * t302;
t65 = t194 + (t220 * t302 - t300 * t75) * t320;
t64 = t193 + (t219 * t302 - t300 * t74) * t318;
t63 = t192 + (t218 * t302 - t300 * t73) * t316;
t59 = -t320 * t348 + t194;
t58 = -t318 * t349 + t193;
t57 = -t316 * t350 + t192;
t56 = t189 + (t217 * t302 - t300 * t71) * t304;
t53 = -t304 * t351 + t189;
t52 = (-t302 * t90 * t429 - (-t319 * t257 * t376 + t302 * (-t296 * t441 + t325 * t581)) * t110 * t457) * t297 + (t160 * t312 + (t174 * t283 - t178 * t287) * t449) * t453;
t51 = (-t302 * t89 * t430 - (-t317 * t256 * t378 + t302 * (-t294 * t442 + t323 * t583)) * t109 * t462) * t295 + (t159 * t312 + (t173 * t282 - t177 * t286) * t450) * t458;
t50 = (-t302 * t88 * t431 - (-t315 * t255 * t380 + t302 * (-t292 * t443 + t321 * t585)) * t108 * t467) * t293 + (t158 * t312 + (t172 * t281 - t176 * t285) * t451) * t463;
t49 = (-t302 * t87 * t432 - (-t303 * t253 * t382 + t302 * (-t290 * t444 + t305 * t587)) * t104 * t472) * t291 + (t155 * t312 + (t171 * t280 - t175 * t284) * t452) * t468;
t48 = t296 * t499 + t52 * t319;
t47 = t294 * t500 + t51 * t317;
t46 = t292 * t501 + t50 * t315;
t45 = -t297 * t319 * t499 + t52 * t325;
t44 = -t295 * t317 * t500 + t51 * t323;
t43 = -t293 * t315 * t501 + t50 * t321;
t42 = t290 * t502 + t49 * t303;
t41 = -t291 * t303 * t502 + t49 * t305;
t40 = t297 * t326 * t364 + t320 * t52;
t39 = t295 * t324 * t365 + t318 * t51;
t38 = t293 * t322 * t366 + t316 * t50;
t37 = t291 * t306 * t367 + t304 * t49;
t32 = t318 * t35 + t324 * t495;
t31 = t318 * t495 - t593;
t30 = t320 * t36 + t326 * t494;
t29 = t316 * t34 + t322 * t496;
t28 = t320 * t494 - t592;
t27 = t316 * t496 - t594;
t26 = t304 * t33 + t306 * t498;
t25 = t304 * t498 - t595;
t24 = (t296 * t364 + t588) * t319;
t23 = (t294 * t365 + t589) * t317;
t22 = (t292 * t366 + t590) * t315;
t21 = (t290 * t367 + t591) * t303;
t20 = t70 * t319 + t65 * t325;
t19 = t65 * t319 - t70 * t325;
t18 = t69 * t317 + t64 * t323;
t17 = t64 * t317 - t69 * t323;
t16 = t68 * t315 + t63 * t321;
t15 = t63 * t315 - t68 * t321;
t14 = 0.2e1 * t325 * t588 + (-0.2e1 * t297 + 0.4e1) * t377;
t13 = 0.2e1 * t323 * t589 + (-0.2e1 * t295 + 0.4e1) * t379;
t12 = 0.2e1 * t321 * t590 + (-0.2e1 * t293 + 0.4e1) * t381;
t11 = 0.2e1 * t305 * t591 + (-0.2e1 * t291 + 0.4e1) * t383;
t10 = t67 * t303 + t56 * t305;
t9 = t56 * t303 - t67 * t305;
t8 = (t384 * t319 - t325 * t40) * t300 - t302 * t48;
t7 = (t385 * t317 - t323 * t39) * t300 - t302 * t47;
t6 = (t386 * t315 - t321 * t38) * t300 - t302 * t46;
t5 = (-t319 * t40 - t384 * t325) * t300 + t302 * t45;
t4 = (-t317 * t39 - t385 * t323) * t300 + t302 * t44;
t3 = (-t315 * t38 - t386 * t321) * t300 + t302 * t43;
t2 = (t387 * t303 - t305 * t37) * t300 - t302 * t42;
t1 = (-t303 * t37 - t387 * t305) * t300 + t302 * t41;
t54 = [t75 * t564 + t74 * t566 + t73 * t568 + t71 * t570, -t33 * t412 - t34 * t410 - t35 * t408 - t36 * t407, -t55 * t412 - t390 - t389 - t388 + (-t25 * t570 - t27 * t568 - t28 * t564 - t31 * t566) * t300, -t53 * t412 - t57 * t410 - t58 * t408 - t59 * t407 + (-t26 * t570 - t29 * t568 - t30 * t564 - t32 * t566) * t300, -t21 * t412 - t22 * t410 - t23 * t408 - t24 * t407 + (t284 * t359 + t285 * t358 + t286 * t357 + t287 * t356) * t342, -t11 * t412 - t12 * t410 - t13 * t408 - t14 * t407 + (-t100 * t420 - t101 * t417 - t102 * t414 - t423 * t99) * t342, -t42 * t412 - t46 * t410 - t47 * t408 - t48 * t407 + (-t284 * t400 - t285 * t399 - t286 * t398 - t287 * t397) * t342, -t41 * t412 - t43 * t410 - t44 * t408 - t45 * t407 + (-t284 * t508 - t285 * t507 - t286 * t506 - t287 * t505) * t342, (-t414 * t52 - t417 * t51 - t420 * t50 - t423 * t49) * t342, -t284 * t55 * t559 - t285 * t504 - t286 * t503 - t287 * t62 * t556 + t1 * t570 + t3 * t568 + t4 * t566 + t5 * t564 + (-t15 * t420 - t17 * t417 - t19 * t414 - t423 * t9) * t342, t8 * t564 + t319 * t388 + t7 * t566 + t317 * t389 + t6 * t568 + t315 * t390 + t2 * t570 + t284 * t396 + (-t10 * t423 - t16 * t420 - t18 * t417 - t20 * t414) * t342, 0, t242, -t241, t600; t75 * t565 + t74 * t567 + t73 * t569 + t71 * t571, t33 * t406 + t34 * t411 + t35 * t409 + t36 * t405, t392 + t391 + t55 * t406 + t62 * t405 + (-t25 * t571 - t27 * t569 - t28 * t565 - t31 * t567) * t300, t57 * t411 + t58 * t409 + t53 * t406 + t59 * t405 + (-t26 * t571 - t29 * t569 - t30 * t565 - t32 * t567) * t300, t22 * t411 + t23 * t409 + t21 * t406 + t24 * t405 + (-t280 * t359 - t281 * t358 - t282 * t357 - t283 * t356) * t342, t11 * t406 + t12 * t411 + t13 * t409 + t14 * t405 + (t100 * t421 + t101 * t418 + t102 * t415 + t424 * t99) * t342, t46 * t411 + t47 * t409 + t42 * t406 + t48 * t405 + (t280 * t400 + t281 * t399 + t282 * t398 + t283 * t397) * t342, t43 * t411 + t44 * t409 + t41 * t406 + t45 * t405 + (t280 * t508 + t281 * t507 + t282 * t506 + t283 * t505) * t342, (t415 * t52 + t418 * t51 + t421 * t50 + t424 * t49) * t342, t281 * t504 + t282 * t503 + t3 * t569 + t4 * t567 + (t283 * t188 * t62 + t129 * t5) * t614 + (t187 * t280 * t55 + t123 * t1) * t611 + (t15 * t421 + t17 * t418 + t19 * t415 + t424 * t9) * t342, -t315 * t392 - t317 * t391 + t6 * t569 + t7 * t567 + (-t319 * t473 * t62 + t129 * t8) * t614 + (-t303 * t474 * t55 + t123 * t2) * t611 + (t10 * t424 + t16 * t421 + t18 * t418 + t20 * t415) * t342, 0, t241, t242, t601; t75 * t560 + t74 * t561 + t73 * t562 + t71 * t563, t33 * t478 + t34 * t477 + t35 * t476 + t36 * t475, t55 * t478 + t440 + t439 + t438 + (-t25 * t563 - t27 * t562 - t28 * t560 - t31 * t561) * t300, t53 * t478 + t57 * t477 + t58 * t476 + t59 * t475 + (-t26 * t563 - t29 * t562 - t30 * t560 - t32 * t561) * t300, t21 * t478 + t22 * t477 + t23 * t476 + t24 * t475 + (-t155 * t375 - t158 * t374 - t159 * t373 - t160 * t372) * t342, t11 * t478 + t12 * t477 + t13 * t476 + t14 * t475 + (t100 * t484 + t101 * t483 + t102 * t482 + t486 * t99) * t342, t42 * t478 + t46 * t477 + t47 * t476 + t48 * t475 + (t155 * t437 + t158 * t435 + t159 * t434 + t160 * t433) * t342, t41 * t478 + t43 * t477 + t44 * t476 + t45 * t475 + (t155 * t599 + t158 * t598 + t159 * t597 + t160 * t596) * t342, (t482 * t52 + t483 * t51 + t484 * t50 + t486 * t49) * t342, t1 * t563 + t3 * t562 + t4 * t561 + t5 * t560 + t180 * t611 * t55 + t184 * t612 * t60 + t185 * t613 * t61 + t186 * t614 * t62 + (t15 * t484 + t17 * t483 + t19 * t482 + t486 * t9) * t342, t8 * t560 - t319 * t438 + t7 * t561 - t317 * t439 + t6 * t562 - t315 * t440 + t2 * t563 - t180 * t436 + (t10 * t486 + t16 * t484 + t18 * t483 + t20 * t482) * t342, 0, 0, 0, t312 - g(3); t111 * t71 + t112 * t73 + t113 * t74 + t114 * t75, t33 * t428 + t34 * t426 + t35 * t425 + t36 * t427, t55 * t428 + t393 + t395 + t394 + (-t111 * t25 - t112 * t27 - t113 * t31 - t114 * t28) * t300, t53 * t428 + t59 * t427 + t57 * t426 + t58 * t425 + (-t111 * t26 - t112 * t29 - t113 * t32 - t114 * t30) * t300, t21 * t428 + t24 * t427 + t22 * t426 + t23 * t425 + (-t368 * t356 - t369 * t357 - t370 * t358 - t371 * t359) * t342, t11 * t428 + t12 * t426 + t13 * t425 + t14 * t427 + (t100 * t419 + t101 * t416 + t102 * t413 + t422 * t99) * t342, t42 * t428 + t48 * t427 + t46 * t426 + t47 * t425 + (t413 * t588 + t416 * t589 + t419 * t590 + t422 * t591) * t342, t41 * t428 + t45 * t427 + t43 * t426 + t44 * t425 + (t33 * t485 + t34 * t481 + t35 * t480 + t36 * t479) * t342, (t413 * t52 + t416 * t51 + t419 * t50 + t422 * t49) * t342, t55 * t490 + t62 * t489 + t60 * t488 + t61 * t487 + t111 * t1 + t112 * t3 + t113 * t4 + t114 * t5 + (t15 * t419 + t17 * t416 + t19 * t413 + t422 * t9) * t342, t114 * t8 - t319 * t393 + t113 * t7 - t317 * t394 + t112 * t6 - t315 * t395 + t111 * t2 - t371 * t396 + (t10 * t422 + t16 * t419 + t18 * t416 + t20 * t413) * t342, t311, -t288 * t600 + t289 * t601, -t288 * t601 - t289 * t600, 0;];
tauX_reg  = t54;
