% Calculate vector of inverse dynamics forces for parallel robot
% P4PRRRR8V2G3A0
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
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% tauX [4x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:30
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P4PRRRR8V2G3A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_slag_vp1: xDP has to be [4x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [4 1]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_slag_vp1: xDDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:24:52
% EndTime: 2020-08-07 11:25:08
% DurationCPUTime: 16.76s
% Computational Cost: add. (79957->784), mult. (172971->1408), div. (9248->9), fcn. (162006->30), ass. (0->494)
t313 = sin(qJ(2,1));
t319 = cos(qJ(2,1));
t327 = pkin(7) + pkin(6);
t228 = pkin(2) * t313 - t319 * t327;
t293 = sin(pkin(4));
t295 = cos(pkin(4));
t312 = sin(qJ(3,1));
t433 = t295 * t312;
t172 = pkin(3) * t433 + t228 * t293;
t318 = cos(qJ(3,1));
t448 = t293 * t313;
t290 = t318 ^ 2;
t539 = pkin(3) * t290;
t140 = 0.1e1 / (pkin(2) * t433 + t172 * t318 + t448 * t539);
t311 = sin(qJ(2,2));
t317 = cos(qJ(2,2));
t227 = pkin(2) * t311 - t317 * t327;
t310 = sin(qJ(3,2));
t435 = t295 * t310;
t171 = pkin(3) * t435 + t227 * t293;
t316 = cos(qJ(3,2));
t450 = t293 * t311;
t289 = t316 ^ 2;
t540 = pkin(3) * t289;
t139 = 0.1e1 / (pkin(2) * t435 + t171 * t316 + t450 * t540);
t309 = sin(qJ(2,3));
t315 = cos(qJ(2,3));
t226 = pkin(2) * t309 - t315 * t327;
t308 = sin(qJ(3,3));
t437 = t295 * t308;
t170 = pkin(3) * t437 + t226 * t293;
t314 = cos(qJ(3,3));
t452 = t293 * t309;
t288 = t314 ^ 2;
t541 = pkin(3) * t288;
t138 = 0.1e1 / (pkin(2) * t437 + t170 * t314 + t452 * t541);
t297 = sin(qJ(2,4));
t299 = cos(qJ(2,4));
t224 = pkin(2) * t297 - t299 * t327;
t296 = sin(qJ(3,4));
t440 = t295 * t296;
t169 = pkin(3) * t440 + t224 * t293;
t298 = cos(qJ(3,4));
t456 = t293 * t297;
t286 = t298 ^ 2;
t542 = pkin(3) * t286;
t137 = 0.1e1 / (pkin(2) * t440 + t169 * t298 + t456 * t542);
t303 = legFrame(1,2);
t274 = sin(t303);
t278 = cos(t303);
t218 = g(1) * t274 + g(2) * t278;
t222 = g(1) * t278 - g(2) * t274;
t294 = cos(pkin(8));
t269 = g(3) * t294;
t292 = sin(pkin(8));
t382 = -t222 * t292 - t269;
t582 = t218 * t293 + t382 * t295;
t302 = legFrame(2,2);
t273 = sin(t302);
t277 = cos(t302);
t217 = g(1) * t273 + g(2) * t277;
t221 = g(1) * t277 - g(2) * t273;
t384 = -t221 * t292 - t269;
t581 = t217 * t293 + t384 * t295;
t301 = legFrame(3,2);
t272 = sin(t301);
t276 = cos(t301);
t216 = g(1) * t272 + g(2) * t276;
t220 = g(1) * t276 - g(2) * t272;
t386 = -t220 * t292 - t269;
t580 = t216 * t293 + t386 * t295;
t300 = legFrame(4,2);
t271 = sin(t300);
t275 = cos(t300);
t215 = g(1) * t271 + g(2) * t275;
t219 = g(1) * t275 - g(2) * t271;
t388 = -t219 * t292 - t269;
t579 = t215 * t293 + t388 * t295;
t570 = -rSges(3,1) * t318 + rSges(3,2) * t312;
t569 = -rSges(3,1) * t316 + rSges(3,2) * t310;
t568 = -rSges(3,1) * t314 + rSges(3,2) * t308;
t567 = -rSges(3,1) * t298 + rSges(3,2) * t296;
t561 = m(3) * rSges(3,1);
t416 = rSges(3,2) * t561;
t263 = -Icges(3,4) + t416;
t328 = pkin(2) * m(3);
t411 = t328 / 0.2e1;
t394 = rSges(3,1) * t411;
t566 = t263 * t286 + t296 * t394;
t565 = t263 * t288 + t308 * t394;
t564 = t263 * t289 + t310 * t394;
t563 = t263 * t290 + t312 * t394;
t562 = 0.2e1 * pkin(2);
t439 = t295 * t297;
t194 = t292 * t439 - t294 * t299;
t455 = t293 * t298;
t158 = t194 * t296 + t292 * t455;
t195 = t292 * t299 + t294 * t439;
t159 = t195 * t296 + t294 * t455;
t324 = xDP(3);
t329 = xP(4);
t284 = sin(t329);
t285 = cos(t329);
t332 = koppelP(4,2);
t336 = koppelP(4,1);
t207 = t284 * t336 + t285 * t332;
t211 = -t284 * t332 + t285 * t336;
t323 = xDP(4);
t325 = xDP(2);
t326 = xDP(1);
t376 = (t211 * t323 + t325) * t271 - (-t207 * t323 + t326) * t275;
t98 = (t158 * t324 + t376 * t159) * t137;
t97 = t98 ^ 2;
t436 = t295 * t309;
t199 = t292 * t436 - t294 * t315;
t447 = t293 * t314;
t160 = t199 * t308 + t292 * t447;
t202 = t292 * t315 + t294 * t436;
t163 = t202 * t308 + t294 * t447;
t333 = koppelP(3,2);
t337 = koppelP(3,1);
t208 = t284 * t337 + t285 * t333;
t212 = -t284 * t333 + t285 * t337;
t375 = (t212 * t323 + t325) * t272 - (-t208 * t323 + t326) * t276;
t102 = (t160 * t324 + t375 * t163) * t138;
t99 = t102 ^ 2;
t434 = t295 * t311;
t200 = t292 * t434 - t294 * t317;
t445 = t293 * t316;
t161 = t200 * t310 + t292 * t445;
t203 = t292 * t317 + t294 * t434;
t164 = t203 * t310 + t294 * t445;
t334 = koppelP(2,2);
t338 = koppelP(2,1);
t209 = t284 * t338 + t285 * t334;
t213 = -t284 * t334 + t285 * t338;
t374 = (t213 * t323 + t325) * t273 - (-t209 * t323 + t326) * t277;
t103 = (t161 * t324 + t374 * t164) * t139;
t100 = t103 ^ 2;
t432 = t295 * t313;
t201 = t292 * t432 - t294 * t319;
t443 = t293 * t318;
t162 = t201 * t312 + t292 * t443;
t204 = t292 * t319 + t294 * t432;
t165 = t204 * t312 + t294 * t443;
t335 = koppelP(1,2);
t339 = koppelP(1,1);
t210 = t284 * t339 + t285 * t335;
t214 = -t284 * t335 + t285 * t339;
t373 = (t214 * t323 + t325) * t274 - (-t210 * t323 + t326) * t278;
t104 = (t162 * t324 + t373 * t165) * t140;
t101 = t104 ^ 2;
t225 = pkin(2) * t299 + t297 * t327;
t438 = t295 * t299;
t441 = t294 * t295;
t538 = pkin(3) * t298;
t125 = (t292 * t297 - t294 * t438) * t538 - t225 * t441 + t224 * t292;
t458 = t292 * t295;
t126 = (t292 * t438 + t294 * t297) * t538 + t225 * t458 + t294 * t224;
t343 = 0.1e1 / pkin(3);
t497 = t137 * t343;
t90 = (-t376 * t125 + t126 * t324) * t497;
t560 = pkin(3) * t90;
t229 = pkin(2) * t315 + t309 * t327;
t431 = t295 * t315;
t537 = pkin(3) * t314;
t128 = (t292 * t309 - t294 * t431) * t537 - t229 * t441 + t226 * t292;
t131 = (t292 * t431 + t294 * t309) * t537 + t229 * t458 + t294 * t226;
t496 = t138 * t343;
t94 = (-t375 * t128 + t131 * t324) * t496;
t559 = pkin(3) * t94;
t230 = pkin(2) * t317 + t311 * t327;
t430 = t295 * t317;
t536 = pkin(3) * t316;
t129 = (t292 * t311 - t294 * t430) * t536 - t230 * t441 + t227 * t292;
t132 = (t292 * t430 + t294 * t311) * t536 + t230 * t458 + t294 * t227;
t495 = t139 * t343;
t95 = (-t374 * t129 + t132 * t324) * t495;
t558 = pkin(3) * t95;
t231 = pkin(2) * t319 + t313 * t327;
t429 = t295 * t319;
t535 = pkin(3) * t318;
t130 = (t292 * t313 - t294 * t429) * t535 - t231 * t441 + t228 * t292;
t133 = (t292 * t429 + t294 * t313) * t535 + t231 * t458 + t294 * t228;
t494 = t140 * t343;
t96 = (-t373 * t130 + t133 * t324) * t494;
t557 = pkin(3) * t96;
t340 = rSges(3,2) ^ 2;
t341 = rSges(3,1) ^ 2;
t244 = (-t340 + t341) * m(3) + Icges(3,2) - Icges(3,1);
t556 = t244 / 0.2e1;
t320 = pkin(6) + rSges(3,3);
t548 = m(3) * t320;
t250 = rSges(3,2) * t548 - Icges(3,6);
t555 = -t250 / 0.4e1;
t251 = rSges(3,1) * t548 - Icges(3,5);
t554 = t251 / 0.4e1;
t553 = -t263 / 0.2e1;
t527 = rSges(3,2) * t298;
t237 = rSges(3,1) * t296 + t527;
t157 = -t237 * t456 - t295 * t567;
t552 = m(3) * t157;
t523 = rSges(3,2) * t314;
t241 = rSges(3,1) * t308 + t523;
t166 = -t241 * t452 - t295 * t568;
t551 = m(3) * t166;
t522 = rSges(3,2) * t316;
t242 = rSges(3,1) * t310 + t522;
t167 = -t242 * t450 - t295 * t569;
t550 = m(3) * t167;
t521 = rSges(3,2) * t318;
t243 = rSges(3,1) * t312 + t521;
t168 = -t243 * t448 - t295 * t570;
t549 = m(3) * t168;
t547 = m(3) * t343;
t546 = pkin(2) * t296;
t545 = pkin(2) * t308;
t544 = pkin(2) * t310;
t543 = pkin(2) * t312;
t534 = g(3) * t292;
t344 = pkin(2) ^ 2;
t266 = t327 ^ 2 + t344;
t342 = pkin(3) ^ 2;
t517 = t296 * t90;
t415 = pkin(3) * t517;
t423 = pkin(3) * t562;
t533 = (-t327 * t415 + (t286 * t342 + t298 * t423 + t266) * t98) * t98;
t515 = t308 * t94;
t414 = pkin(3) * t515;
t520 = t102 * (-t327 * t414 + (t288 * t342 + t314 * t423 + t266) * t102);
t514 = t310 * t95;
t413 = pkin(3) * t514;
t519 = t103 * (-t327 * t413 + (t289 * t342 + t316 * t423 + t266) * t103);
t513 = t312 * t96;
t412 = pkin(3) * t513;
t518 = t104 * (-t327 * t412 + (t290 * t342 + t318 * t423 + t266) * t104);
t516 = t299 * t98;
t512 = t327 * t98;
t511 = t102 * t315;
t510 = t102 * t327;
t509 = t103 * t317;
t508 = t103 * t327;
t507 = t104 * t319;
t506 = t104 * t327;
t505 = t125 * t343;
t504 = t126 * t343;
t503 = t128 * t343;
t502 = t129 * t343;
t501 = t130 * t343;
t500 = t131 * t343;
t499 = t132 * t343;
t498 = t133 * t343;
t259 = m(2) * rSges(2,1) + t328;
t193 = -m(3) * t567 + t259;
t246 = m(2) * rSges(2,2) - t548;
t153 = t193 * t299 - t246 * t297;
t493 = t153 * t293;
t196 = -m(3) * t568 + t259;
t154 = t196 * t315 - t246 * t309;
t492 = t154 * t293;
t197 = -m(3) * t569 + t259;
t155 = t197 * t317 - t246 * t311;
t491 = t155 * t293;
t198 = -m(3) * t570 + t259;
t156 = t198 * t319 - t246 * t313;
t490 = t156 * t293;
t489 = t159 * t271;
t488 = t159 * t275;
t487 = t163 * t272;
t486 = t163 * t276;
t485 = t164 * t273;
t484 = t164 * t277;
t483 = t165 * t274;
t482 = t165 * t278;
t173 = -t250 * t298 - t296 * t251;
t481 = t173 * t343;
t174 = -t250 * t314 - t308 * t251;
t480 = t174 * t343;
t175 = -t250 * t316 - t310 * t251;
t479 = t175 * t343;
t176 = -t250 * t318 - t312 * t251;
t478 = t176 * t343;
t476 = t215 * t295;
t474 = t216 * t295;
t472 = t217 * t295;
t470 = t218 * t295;
t245 = (t340 + t341) * m(3) + Icges(3,3);
t465 = t245 * t343;
t460 = rSges(3,2) * t293 * t269;
t459 = t292 * t293;
t457 = t293 * t296;
t454 = t293 * t299;
t453 = t293 * t308;
t451 = t293 * t310;
t449 = t293 * t312;
t446 = t293 * t315;
t444 = t293 * t317;
t442 = t293 * t319;
t428 = t295 * t343;
t427 = t297 * t298;
t426 = t309 * t314;
t425 = t311 * t316;
t424 = t313 * t318;
t418 = -t416 / 0.2e1 + Icges(3,4) / 0.2e1;
t417 = -0.2e1 * rSges(3,2) * pkin(2);
t410 = t157 * t547;
t409 = t166 * t547;
t408 = t167 * t547;
t407 = t168 * t547;
t406 = t296 * t512;
t405 = t308 * t510;
t404 = t310 * t508;
t403 = t312 * t506;
t402 = t271 * t505;
t401 = t275 * t505;
t400 = t272 * t503;
t399 = t276 * t503;
t398 = t273 * t502;
t397 = t277 * t502;
t396 = t274 * t501;
t395 = t278 * t501;
t268 = rSges(3,2) * t411;
t393 = (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,1) + Icges(2,3);
t387 = -t219 * t294 + t534;
t385 = -t220 * t294 + t534;
t383 = -t221 * t294 + t534;
t381 = -t222 * t294 + t534;
t372 = t207 * t275 + t211 * t271;
t371 = t208 * t276 + t212 * t272;
t370 = t209 * t277 + t213 * t273;
t369 = t210 * t278 + t214 * t274;
t247 = t320 ^ 2 + t340 + t344;
t264 = 0.2e1 * t263;
t280 = t561 * t562;
t127 = t244 * t286 + (-t264 * t296 + t280) * t298 + (t296 * t417 + t247) * m(3) + t393;
t41 = t406 - t560;
t13 = (((t295 * t90 + t98 * t454) * t542 + ((-t415 + t512) * t297 + pkin(2) * t516) * t455 + t295 * t41) * t98 + (t90 * t454 + (t286 * t295 - t427 * t457 - t295) * t98) * t560) * t137;
t189 = pkin(3) * t427 + t224;
t17 = t137 * t428 * t533 + (-t295 * t406 + (-t189 * t457 + t295 * (pkin(2) * t298 + t542)) * t90) / (t189 * t455 + (pkin(2) + t538) * t440) * t90;
t21 = (-t298 * t533 - (pkin(2) * t90 - t298 * t41) * t560) * t137;
t5 = -t21 * t493 - t127 * t13 - t173 * t17 - 0.4e1 * t90 * ((t555 * t296 + t554 * t298) * t90 + ((t296 * t556 + t268) * t298 + t553 + t566) * t98) + (-t193 * t579 - t387 * t246) * t299 - t297 * (t387 * t193 - t246 * t579);
t223 = (t341 / 0.2e1 - t340 / 0.2e1) * m(3) - Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1;
t348 = t579 * t297 - t387 * t299;
t9 = -t21 * t552 - t173 * t13 - t245 * t17 + 0.2e1 * ((t223 * t296 + t268) * t298 + t418 + t566) * t97 + m(3) * (((t388 * t293 - t476) * rSges(3,1) + t348 * rSges(3,2)) * t298 + (t460 + (t219 * t459 + t476) * rSges(3,2) + t348 * rSges(3,1)) * t296);
t368 = t159 * t5 - t9 * t505;
t43 = t405 - t559;
t14 = (((t102 * t446 + t295 * t94) * t541 + ((-t414 + t510) * t309 + pkin(2) * t511) * t447 + t295 * t43) * t102 + (t94 * t446 + (t288 * t295 - t426 * t453 - t295) * t102) * t559) * t138;
t190 = pkin(3) * t426 + t226;
t18 = t138 * t428 * t520 + (-t295 * t405 + (-t190 * t453 + (pkin(2) * t314 + t541) * t295) * t94) / (t190 * t447 + (pkin(2) + t537) * t437) * t94;
t22 = (-t314 * t520 - (pkin(2) * t94 - t314 * t43) * t559) * t138;
t347 = t580 * t309 - t385 * t315;
t10 = -t22 * t551 - t174 * t14 - t245 * t18 + 0.2e1 * ((t223 * t308 + t268) * t314 + t418 + t565) * t99 + m(3) * (((t386 * t293 - t474) * rSges(3,1) + t347 * rSges(3,2)) * t314 + (t460 + (t220 * t459 + t474) * rSges(3,2) + t347 * rSges(3,1)) * t308);
t134 = t244 * t288 + (-t264 * t308 + t280) * t314 + (t308 * t417 + t247) * m(3) + t393;
t6 = -t22 * t492 - t134 * t14 - t174 * t18 - 0.4e1 * t94 * ((t555 * t308 + t554 * t314) * t94 + ((t308 * t556 + t268) * t314 + t553 + t565) * t102) + (-t196 * t580 - t385 * t246) * t315 - t309 * (t385 * t196 - t246 * t580);
t367 = -t10 * t503 + t163 * t6;
t44 = t404 - t558;
t15 = (((t103 * t444 + t295 * t95) * t540 + ((-t413 + t508) * t311 + pkin(2) * t509) * t445 + t295 * t44) * t103 + (t95 * t444 + (t289 * t295 - t425 * t451 - t295) * t103) * t558) * t139;
t191 = pkin(3) * t425 + t227;
t19 = t139 * t428 * t519 + (-t295 * t404 + (-t191 * t451 + (pkin(2) * t316 + t540) * t295) * t95) / (t191 * t445 + (pkin(2) + t536) * t435) * t95;
t23 = (-t316 * t519 - (pkin(2) * t95 - t316 * t44) * t558) * t139;
t346 = t581 * t311 - t383 * t317;
t11 = -t23 * t550 - t175 * t15 - t245 * t19 + 0.2e1 * ((t223 * t310 + t268) * t316 + t418 + t564) * t100 + m(3) * (((t384 * t293 - t472) * rSges(3,1) + t346 * rSges(3,2)) * t316 + (t460 + (t221 * t459 + t472) * rSges(3,2) + t346 * rSges(3,1)) * t310);
t135 = t244 * t289 + (-t264 * t310 + t280) * t316 + (t310 * t417 + t247) * m(3) + t393;
t7 = -t23 * t491 - t135 * t15 - t175 * t19 - 0.4e1 * t95 * ((t555 * t310 + t554 * t316) * t95 + ((t310 * t556 + t268) * t316 + t553 + t564) * t103) + (-t197 * t581 - t383 * t246) * t317 - t311 * (t383 * t197 - t246 * t581);
t366 = -t11 * t502 + t164 * t7;
t45 = t403 - t557;
t16 = (((t104 * t442 + t295 * t96) * t539 + ((-t412 + t506) * t313 + pkin(2) * t507) * t443 + t295 * t45) * t104 + (t96 * t442 + (t290 * t295 - t424 * t449 - t295) * t104) * t557) * t140;
t192 = pkin(3) * t424 + t228;
t20 = t140 * t428 * t518 + (-t295 * t403 + (-t192 * t449 + (pkin(2) * t318 + t539) * t295) * t96) / (t192 * t443 + (pkin(2) + t535) * t433) * t96;
t24 = (-t318 * t518 - (pkin(2) * t96 - t318 * t45) * t557) * t140;
t345 = t582 * t313 - t381 * t319;
t12 = -t24 * t549 - t176 * t16 - t245 * t20 + 0.2e1 * ((t223 * t312 + t268) * t318 + t418 + t563) * t101 + m(3) * (((t382 * t293 - t470) * rSges(3,1) + t345 * rSges(3,2)) * t318 + (t460 + (t222 * t459 + t470) * rSges(3,2) + t345 * rSges(3,1)) * t312);
t136 = t244 * t290 + (-t264 * t312 + t280) * t318 + (t312 * t417 + t247) * m(3) + t393;
t8 = -t24 * t490 - t136 * t16 - t176 * t20 - 0.4e1 * t96 * ((t555 * t312 + t554 * t318) * t96 + ((t312 * t556 + t268) * t318 + t553 + t563) * t104) + (-t198 * t582 - t381 * t246) * t319 - t313 * (t381 * t198 - t246 * t582);
t365 = -t12 * t501 + t165 * t8;
t364 = pkin(3) * t457 - t224 * t295;
t363 = pkin(3) * t453 - t226 * t295;
t362 = pkin(3) * t451 - t227 * t295;
t361 = pkin(3) * t449 - t228 * t295;
t360 = t125 * t481 - t127 * t159;
t359 = t125 * t465 - t159 * t173;
t358 = t128 * t480 - t134 * t163;
t357 = t128 * t465 - t163 * t174;
t356 = t129 * t479 - t135 * t164;
t355 = t129 * t465 - t164 * t175;
t354 = t130 * t478 - t136 * t165;
t353 = t130 * t465 - t165 * t176;
t352 = t125 * t410 - t159 * t493;
t351 = t128 * t409 - t163 * t492;
t350 = t129 * t408 - t164 * t491;
t349 = t130 * t407 - t165 * t490;
t331 = rSges(4,1);
t330 = rSges(4,2);
t307 = xDDP(1);
t306 = xDDP(2);
t305 = xDDP(3);
t304 = xDDP(4);
t291 = t323 ^ 2;
t287 = m(1) + m(2) + m(3);
t206 = -t284 * t330 + t285 * t331;
t205 = t284 * t331 + t285 * t330;
t152 = -t210 * t304 - t214 * t291 + t307;
t151 = -t209 * t304 - t213 * t291 + t307;
t150 = -t208 * t304 - t212 * t291 + t307;
t149 = -t207 * t304 - t211 * t291 + t307;
t148 = -t210 * t291 + t214 * t304 + t306;
t147 = -t209 * t291 + t213 * t304 + t306;
t146 = -t208 * t291 + t212 * t304 + t306;
t145 = -t207 * t291 + t211 * t304 + t306;
t144 = t231 * t294 + t361 * t292;
t143 = t230 * t294 + t362 * t292;
t142 = t229 * t294 + t363 * t292;
t141 = t225 * t294 + t364 * t292;
t124 = -t204 * t539 - t231 * t292 * t318 + (pkin(2) * t449 + t361 * t318) * t294;
t123 = -t203 * t540 - t230 * t292 * t316 + (pkin(2) * t451 + t362 * t316) * t294;
t122 = -t202 * t541 - t229 * t292 * t314 + (pkin(2) * t453 + t363 * t314) * t294;
t121 = -t195 * t542 - t225 * t292 * t298 + (pkin(2) * t457 + t364 * t298) * t294;
t120 = -(t201 * t278 - t274 * t448) * t539 + (t144 * t278 + t172 * t274) * t318 + (t274 * t295 + t278 * t459) * t543;
t119 = (t201 * t274 + t278 * t448) * t539 + (-t144 * t274 + t172 * t278) * t318 + (-t274 * t459 + t278 * t295) * t543;
t118 = -(t200 * t277 - t273 * t450) * t540 + (t143 * t277 + t171 * t273) * t316 + (t273 * t295 + t277 * t459) * t544;
t117 = (t200 * t273 + t277 * t450) * t540 + (-t143 * t273 + t171 * t277) * t316 + (-t273 * t459 + t277 * t295) * t544;
t116 = -(t199 * t276 - t272 * t452) * t541 + (t142 * t276 + t170 * t272) * t314 + (t272 * t295 + t276 * t459) * t545;
t115 = (t199 * t272 + t276 * t452) * t541 + (-t142 * t272 + t170 * t276) * t314 + (-t272 * t459 + t276 * t295) * t545;
t114 = -(t194 * t275 - t271 * t456) * t542 + (t141 * t275 + t169 * t271) * t298 + (t271 * t295 + t275 * t459) * t546;
t113 = (t194 * t271 + t275 * t456) * t542 + (-t141 * t271 + t169 * t275) * t298 + (-t271 * t459 + t275 * t295) * t546;
t112 = t369 * t165 * t140;
t111 = t370 * t164 * t139;
t110 = t371 * t163 * t138;
t109 = t372 * t159 * t137;
t108 = t369 * t130 * t494;
t107 = t370 * t129 * t495;
t106 = t371 * t128 * t496;
t105 = t372 * t125 * t497;
t93 = t96 ^ 2;
t92 = t95 ^ 2;
t91 = t94 ^ 2;
t89 = (t124 * t549 + t133 * t465 + t162 * t176) * t140;
t88 = (t123 * t550 + t132 * t465 + t161 * t175) * t139;
t87 = (t122 * t551 + t131 * t465 + t160 * t174) * t138;
t86 = t90 ^ 2;
t85 = (t124 * t287 + t133 * t407 + t162 * t490) * t140;
t84 = (t123 * t287 + t132 * t408 + t161 * t491) * t139;
t83 = (t122 * t287 + t131 * t409 + t160 * t492) * t138;
t82 = (t121 * t552 + t126 * t465 + t158 * t173) * t137;
t81 = (t121 * t287 + t126 * t410 + t158 * t493) * t137;
t80 = (t119 * t214 - t120 * t210) * t140;
t79 = (t117 * t213 - t118 * t209) * t139;
t78 = (t115 * t212 - t116 * t208) * t138;
t77 = (t124 * t490 + t133 * t478 + t136 * t162) * t140;
t76 = (t123 * t491 + t132 * t479 + t135 * t161) * t139;
t75 = (t122 * t492 + t131 * t480 + t134 * t160) * t138;
t74 = (t113 * t211 - t114 * t207) * t137;
t73 = (t121 * t493 + t126 * t481 + t127 * t158) * t137;
t72 = (t120 * t549 + t353 * t278) * t140;
t71 = (t119 * t549 - t353 * t274) * t140;
t70 = (t118 * t550 + t355 * t277) * t139;
t69 = (t117 * t550 - t355 * t273) * t139;
t68 = (t116 * t551 + t357 * t276) * t138;
t67 = (t115 * t551 - t357 * t272) * t138;
t66 = (t120 * t287 + t349 * t278) * t140;
t65 = (t119 * t287 - t349 * t274) * t140;
t64 = (t118 * t287 + t350 * t277) * t139;
t63 = (t117 * t287 - t350 * t273) * t139;
t62 = (t116 * t287 + t351 * t276) * t138;
t61 = (t115 * t287 - t351 * t272) * t138;
t60 = (t114 * t552 + t359 * t275) * t137;
t59 = (t113 * t552 - t359 * t271) * t137;
t58 = (t114 * t287 + t352 * t275) * t137;
t57 = (t113 * t287 - t352 * t271) * t137;
t56 = (t120 * t490 + t354 * t278) * t140;
t55 = (t119 * t490 - t354 * t274) * t140;
t54 = (t118 * t491 + t356 * t277) * t139;
t53 = (t117 * t491 - t356 * t273) * t139;
t52 = (t116 * t492 + t358 * t276) * t138;
t51 = (t115 * t492 - t358 * t272) * t138;
t50 = (t114 * t493 + t360 * t275) * t137;
t49 = (t113 * t493 - t360 * t271) * t137;
t40 = -t108 * t245 + t112 * t176 + t80 * t549;
t39 = -t107 * t245 + t111 * t175 + t79 * t550;
t38 = -t106 * t245 + t110 * t174 + t78 * t551;
t37 = -t108 * t549 + t112 * t490 + t287 * t80;
t36 = -t107 * t550 + t111 * t491 + t287 * t79;
t35 = -t106 * t551 + t110 * t492 + t287 * t78;
t34 = -t105 * t245 + t109 * t173 + t74 * t552;
t33 = -t105 * t552 + t109 * t493 + t287 * t74;
t32 = -t108 * t176 + t112 * t136 + t80 * t490;
t31 = -t107 * t175 + t111 * t135 + t79 * t491;
t30 = -t106 * t174 + t110 * t134 + t78 * t492;
t29 = -t105 * t173 + t109 * t127 + t74 * t493;
t4 = (-t156 * t16 + (-t246 * t319 - t259 * t313) * t101) * t293 + (-t24 - t218) * t287 + (-t168 * t20 + (-0.2e1 * (rSges(3,1) * t513 + t96 * t521) * t507 + t570 * t313 * (t101 + t93)) * t293 - t93 * t295 * t243) * m(3);
t3 = (-t15 * t155 + (-t246 * t317 - t259 * t311) * t100) * t293 + (-t23 - t217) * t287 + (-t167 * t19 + (-0.2e1 * (rSges(3,1) * t514 + t95 * t522) * t509 + t569 * t311 * (t100 + t92)) * t293 - t92 * t295 * t242) * m(3);
t2 = (-t14 * t154 + (-t246 * t315 - t259 * t309) * t99) * t293 + (-t22 - t216) * t287 + (-t166 * t18 + (-0.2e1 * (rSges(3,1) * t515 + t94 * t523) * t511 + t568 * t309 * (t99 + t91)) * t293 - t91 * t295 * t241) * m(3);
t1 = (-t13 * t153 + (-t246 * t299 - t259 * t297) * t97) * t293 + (-t21 - t215) * t287 + (-t157 * t17 + (-0.2e1 * (rSges(3,1) * t517 + t90 * t527) * t516 + t567 * t297 * (t97 + t86)) * t293 - t86 * t295 * t237) * m(3);
t25 = [(-t205 * t304 - t291 * t206 - g(1) + t307) * m(4) + ((t119 * t66 - t72 * t396 + t56 * t483) * t148 + (t124 * t66 + t162 * t56 + t72 * t498) * t305 + (t66 * t152 + t4) * t120 + ((-t165 * t56 + t72 * t501) * t152 - t365) * t278) * t140 + ((t117 * t64 - t70 * t398 + t54 * t485) * t147 + (t123 * t64 + t161 * t54 + t70 * t499) * t305 + (t64 * t151 + t3) * t118 + ((-t164 * t54 + t70 * t502) * t151 - t366) * t277) * t139 + ((t115 * t62 - t68 * t400 + t52 * t487) * t146 + (t122 * t62 + t160 * t52 + t68 * t500) * t305 + (t62 * t150 + t2) * t116 + ((-t163 * t52 + t68 * t503) * t150 - t367) * t276) * t138 + ((t113 * t58 - t60 * t402 + t50 * t489) * t145 + (t121 * t58 + t158 * t50 + t60 * t504) * t305 + (t58 * t149 + t1) * t114 + ((-t159 * t50 + t60 * t505) * t149 - t368) * t275) * t137; (-t291 * t205 + t206 * t304 - g(2) + t306) * m(4) + ((t124 * t65 + t162 * t55 + t71 * t498) * t305 + (t120 * t65 + t395 * t71 - t55 * t482) * t152 + (t65 * t148 + t4) * t119 + ((t165 * t55 - t71 * t501) * t148 + t365) * t274) * t140 + ((t123 * t63 + t161 * t53 + t69 * t499) * t305 + (t118 * t63 + t397 * t69 - t53 * t484) * t151 + (t63 * t147 + t3) * t117 + ((t164 * t53 - t69 * t502) * t147 + t366) * t273) * t139 + ((t122 * t61 + t160 * t51 + t67 * t500) * t305 + (t116 * t61 + t399 * t67 - t51 * t486) * t150 + (t61 * t146 + t2) * t115 + ((t163 * t51 - t67 * t503) * t146 + t367) * t272) * t138 + ((t114 * t57 + t401 * t59 - t49 * t488) * t149 + (t121 * t57 + t158 * t49 + t59 * t504) * t305 + (t57 * t145 + t1) * t113 + ((t159 * t49 - t59 * t505) * t145 + t368) * t271) * t137; (t305 - g(3)) * m(4) + ((t119 * t85 - t396 * t89 + t77 * t483) * t148 + (t124 * t85 + t162 * t77 + t89 * t498) * t305 + (t120 * t85 + t395 * t89 - t77 * t482) * t152 + t162 * t8 + t124 * t4 + t12 * t498) * t140 + ((t117 * t84 - t398 * t88 + t76 * t485) * t147 + (t123 * t84 + t161 * t76 + t88 * t499) * t305 + (t118 * t84 + t397 * t88 - t76 * t484) * t151 + t161 * t7 + t123 * t3 + t11 * t499) * t139 + ((t115 * t83 - t400 * t87 + t75 * t487) * t146 + (t122 * t83 + t160 * t75 + t87 * t500) * t305 + (t116 * t83 + t399 * t87 - t75 * t486) * t150 + t160 * t6 + t122 * t2 + t10 * t500) * t138 + ((t113 * t81 - t402 * t82 + t73 * t489) * t145 + (t114 * t81 + t401 * t82 - t73 * t488) * t149 + (t121 * t81 + t158 * t73 + t82 * t504) * t305 + t121 * t1 + t158 * t5 + t9 * t504) * t137; Icges(4,3) * t304 + t74 * t1 - t106 * t10 - t105 * t9 - t107 * t11 - t108 * t12 + t109 * t5 + t110 * t6 + t111 * t7 + t112 * t8 + t78 * t2 + t79 * t3 + t80 * t4 + ((g(1) * t331 + g(2) * t330) * t284 + (g(1) * t330 - g(2) * t331) * t285 + (t330 ^ 2 + t331 ^ 2) * t304 + t206 * t306 - t205 * t307) * m(4) + ((t119 * t37 + t32 * t483 - t396 * t40) * t148 + (t124 * t37 + t162 * t32 + t40 * t498) * t305 + (t120 * t37 - t32 * t482 + t395 * t40) * t152) * t140 + ((t117 * t36 + t31 * t485 - t39 * t398) * t147 + (t123 * t36 + t161 * t31 + t39 * t499) * t305 + (t118 * t36 - t31 * t484 + t39 * t397) * t151) * t139 + ((t115 * t35 + t30 * t487 - t38 * t400) * t146 + (t122 * t35 + t160 * t30 + t38 * t500) * t305 + (t116 * t35 - t30 * t486 + t38 * t399) * t150) * t138 + ((t113 * t33 + t29 * t489 - t34 * t402) * t145 + (t114 * t33 - t29 * t488 + t34 * t401) * t149 + (t121 * t33 + t158 * t29 + t34 * t504) * t305) * t137;];
tauX  = t25;
