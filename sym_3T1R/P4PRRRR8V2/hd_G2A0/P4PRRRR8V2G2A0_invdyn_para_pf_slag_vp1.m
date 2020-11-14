% Calculate vector of inverse dynamics forces for parallel robot
% P4PRRRR8V2G2A0
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
% Datum: 2020-08-07 11:22
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P4PRRRR8V2G2A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_slag_vp1: xDP has to be [4x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [4 1]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_slag_vp1: xDDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:15:55
% EndTime: 2020-08-07 11:16:11
% DurationCPUTime: 16.54s
% Computational Cost: add. (79957->788), mult. (172971->1414), div. (9248->9), fcn. (162006->30), ass. (0->496)
t314 = sin(qJ(2,1));
t320 = cos(qJ(2,1));
t328 = pkin(7) + pkin(6);
t229 = pkin(2) * t314 - t320 * t328;
t294 = sin(pkin(4));
t296 = cos(pkin(4));
t313 = sin(qJ(3,1));
t434 = t296 * t313;
t172 = pkin(3) * t434 + t229 * t294;
t319 = cos(qJ(3,1));
t453 = t294 * t314;
t291 = t319 ^ 2;
t541 = pkin(3) * t291;
t140 = 0.1e1 / (pkin(2) * t434 + t172 * t319 + t453 * t541);
t312 = sin(qJ(2,2));
t318 = cos(qJ(2,2));
t228 = pkin(2) * t312 - t318 * t328;
t311 = sin(qJ(3,2));
t436 = t296 * t311;
t171 = pkin(3) * t436 + t228 * t294;
t317 = cos(qJ(3,2));
t455 = t294 * t312;
t290 = t317 ^ 2;
t542 = pkin(3) * t290;
t139 = 0.1e1 / (pkin(2) * t436 + t171 * t317 + t455 * t542);
t310 = sin(qJ(2,3));
t316 = cos(qJ(2,3));
t227 = pkin(2) * t310 - t316 * t328;
t309 = sin(qJ(3,3));
t438 = t296 * t309;
t170 = pkin(3) * t438 + t227 * t294;
t315 = cos(qJ(3,3));
t457 = t294 * t310;
t289 = t315 ^ 2;
t543 = pkin(3) * t289;
t138 = 0.1e1 / (pkin(2) * t438 + t170 * t315 + t457 * t543);
t298 = sin(qJ(2,4));
t300 = cos(qJ(2,4));
t225 = pkin(2) * t298 - t300 * t328;
t297 = sin(qJ(3,4));
t441 = t296 * t297;
t169 = pkin(3) * t441 + t225 * t294;
t299 = cos(qJ(3,4));
t461 = t294 * t298;
t287 = t299 ^ 2;
t544 = pkin(3) * t287;
t137 = 0.1e1 / (pkin(2) * t441 + t169 * t299 + t461 * t544);
t304 = legFrame(1,2);
t275 = sin(t304);
t279 = cos(t304);
t218 = g(1) * t275 + g(2) * t279;
t222 = g(1) * t279 - g(2) * t275;
t293 = sin(pkin(8));
t270 = g(3) * t293;
t295 = cos(pkin(8));
t382 = t222 * t295 - t270;
t584 = t218 * t294 + t382 * t296;
t303 = legFrame(2,2);
t274 = sin(t303);
t278 = cos(t303);
t217 = g(1) * t274 + g(2) * t278;
t221 = g(1) * t278 - g(2) * t274;
t384 = t221 * t295 - t270;
t583 = t217 * t294 + t384 * t296;
t302 = legFrame(3,2);
t273 = sin(t302);
t277 = cos(t302);
t216 = g(1) * t273 + g(2) * t277;
t220 = g(1) * t277 - g(2) * t273;
t386 = t220 * t295 - t270;
t582 = t216 * t294 + t386 * t296;
t301 = legFrame(4,2);
t272 = sin(t301);
t276 = cos(t301);
t215 = g(1) * t272 + g(2) * t276;
t219 = g(1) * t276 - g(2) * t272;
t388 = t219 * t295 - t270;
t581 = t215 * t294 + t388 * t296;
t572 = -rSges(3,1) * t319 + rSges(3,2) * t313;
t571 = -rSges(3,1) * t317 + rSges(3,2) * t311;
t570 = -rSges(3,1) * t315 + rSges(3,2) * t309;
t569 = -rSges(3,1) * t299 + rSges(3,2) * t297;
t563 = m(3) * rSges(3,1);
t417 = rSges(3,2) * t563;
t264 = -Icges(3,4) + t417;
t329 = pkin(2) * m(3);
t412 = t329 / 0.2e1;
t395 = rSges(3,1) * t412;
t568 = t264 * t287 + t297 * t395;
t567 = t264 * t289 + t309 * t395;
t566 = t264 * t290 + t311 * t395;
t565 = t264 * t291 + t313 * t395;
t564 = 0.2e1 * pkin(2);
t440 = t296 * t298;
t194 = t293 * t440 - t295 * t300;
t460 = t294 * t299;
t158 = t194 * t297 + t293 * t460;
t195 = t293 * t300 + t295 * t440;
t445 = t295 * t299;
t159 = -t195 * t297 - t294 * t445;
t325 = xDP(3);
t330 = xP(4);
t285 = sin(t330);
t286 = cos(t330);
t333 = koppelP(4,2);
t337 = koppelP(4,1);
t207 = t285 * t337 + t286 * t333;
t211 = -t285 * t333 + t286 * t337;
t324 = xDP(4);
t326 = xDP(2);
t327 = xDP(1);
t377 = (t211 * t324 + t326) * t272 - (-t207 * t324 + t327) * t276;
t98 = (t377 * t158 + t159 * t325) * t137;
t97 = t98 ^ 2;
t437 = t296 * t310;
t199 = t293 * t437 - t295 * t316;
t452 = t294 * t315;
t160 = t199 * t309 + t293 * t452;
t202 = t293 * t316 + t295 * t437;
t444 = t295 * t315;
t163 = -t202 * t309 - t294 * t444;
t334 = koppelP(3,2);
t338 = koppelP(3,1);
t208 = t285 * t338 + t286 * t334;
t212 = -t285 * t334 + t286 * t338;
t376 = (t212 * t324 + t326) * t273 - (-t208 * t324 + t327) * t277;
t102 = (t376 * t160 + t163 * t325) * t138;
t99 = t102 ^ 2;
t435 = t296 * t312;
t200 = t293 * t435 - t295 * t318;
t450 = t294 * t317;
t161 = t200 * t311 + t293 * t450;
t203 = t293 * t318 + t295 * t435;
t443 = t295 * t317;
t164 = -t203 * t311 - t294 * t443;
t335 = koppelP(2,2);
t339 = koppelP(2,1);
t209 = t285 * t339 + t286 * t335;
t213 = -t285 * t335 + t286 * t339;
t375 = (t213 * t324 + t326) * t274 - (-t209 * t324 + t327) * t278;
t103 = (t375 * t161 + t164 * t325) * t139;
t100 = t103 ^ 2;
t433 = t296 * t314;
t201 = t293 * t433 - t295 * t320;
t448 = t294 * t319;
t162 = t201 * t313 + t293 * t448;
t204 = t293 * t320 + t295 * t433;
t442 = t295 * t319;
t165 = -t204 * t313 - t294 * t442;
t336 = koppelP(1,2);
t340 = koppelP(1,1);
t210 = t285 * t340 + t286 * t336;
t214 = -t285 * t336 + t286 * t340;
t374 = (t214 * t324 + t326) * t275 - (-t210 * t324 + t327) * t279;
t104 = (t374 * t162 + t165 * t325) * t140;
t101 = t104 ^ 2;
t226 = pkin(2) * t300 + t298 * t328;
t439 = t296 * t300;
t446 = t295 * t296;
t540 = pkin(3) * t299;
t125 = (t293 * t298 - t295 * t439) * t540 - t226 * t446 + t225 * t293;
t464 = t293 * t296;
t126 = (t293 * t439 + t295 * t298) * t540 + t226 * t464 + t225 * t295;
t344 = 0.1e1 / pkin(3);
t90 = (t125 * t325 + t377 * t126) * t344 * t137;
t562 = pkin(3) * t90;
t230 = pkin(2) * t316 + t310 * t328;
t432 = t296 * t316;
t539 = pkin(3) * t315;
t128 = (t293 * t310 - t295 * t432) * t539 - t230 * t446 + t227 * t293;
t131 = (t293 * t432 + t295 * t310) * t539 + t230 * t464 + t227 * t295;
t94 = (t128 * t325 + t376 * t131) * t344 * t138;
t561 = pkin(3) * t94;
t231 = pkin(2) * t318 + t312 * t328;
t431 = t296 * t318;
t538 = pkin(3) * t317;
t129 = (t293 * t312 - t295 * t431) * t538 - t231 * t446 + t228 * t293;
t132 = (t293 * t431 + t295 * t312) * t538 + t231 * t464 + t228 * t295;
t95 = (t129 * t325 + t375 * t132) * t344 * t139;
t560 = pkin(3) * t95;
t232 = pkin(2) * t320 + t314 * t328;
t430 = t296 * t320;
t537 = pkin(3) * t319;
t130 = (t293 * t314 - t295 * t430) * t537 - t232 * t446 + t229 * t293;
t133 = (t293 * t430 + t295 * t314) * t537 + t232 * t464 + t229 * t295;
t96 = (t130 * t325 + t374 * t133) * t344 * t140;
t559 = pkin(3) * t96;
t341 = rSges(3,2) ^ 2;
t342 = rSges(3,1) ^ 2;
t245 = (-t341 + t342) * m(3) + Icges(3,2) - Icges(3,1);
t558 = t245 / 0.2e1;
t321 = pkin(6) + rSges(3,3);
t550 = m(3) * t321;
t251 = rSges(3,2) * t550 - Icges(3,6);
t557 = -t251 / 0.4e1;
t252 = rSges(3,1) * t550 - Icges(3,5);
t556 = t252 / 0.4e1;
t555 = -t264 / 0.2e1;
t529 = rSges(3,2) * t299;
t238 = rSges(3,1) * t297 + t529;
t157 = -t238 * t461 - t296 * t569;
t554 = m(3) * t157;
t525 = rSges(3,2) * t315;
t242 = rSges(3,1) * t309 + t525;
t166 = -t242 * t457 - t296 * t570;
t553 = m(3) * t166;
t524 = rSges(3,2) * t317;
t243 = rSges(3,1) * t311 + t524;
t167 = -t243 * t455 - t296 * t571;
t552 = m(3) * t167;
t523 = rSges(3,2) * t319;
t244 = rSges(3,1) * t313 + t523;
t168 = -t244 * t453 - t296 * t572;
t551 = m(3) * t168;
t549 = m(3) * t344;
t548 = pkin(2) * t297;
t547 = pkin(2) * t309;
t546 = pkin(2) * t311;
t545 = pkin(2) * t313;
t536 = g(3) * t295;
t345 = pkin(2) ^ 2;
t267 = t328 ^ 2 + t345;
t343 = pkin(3) ^ 2;
t519 = t297 * t90;
t416 = pkin(3) * t519;
t424 = pkin(3) * t564;
t535 = (-t328 * t416 + (t287 * t343 + t299 * t424 + t267) * t98) * t98;
t517 = t309 * t94;
t415 = pkin(3) * t517;
t522 = t102 * (-t328 * t415 + (t289 * t343 + t315 * t424 + t267) * t102);
t516 = t311 * t95;
t414 = pkin(3) * t516;
t521 = t103 * (-t328 * t414 + (t290 * t343 + t317 * t424 + t267) * t103);
t515 = t313 * t96;
t413 = pkin(3) * t515;
t520 = t104 * (-t328 * t413 + (t291 * t343 + t319 * t424 + t267) * t104);
t518 = t300 * t98;
t514 = t328 * t98;
t513 = t102 * t316;
t512 = t102 * t328;
t511 = t103 * t318;
t510 = t103 * t328;
t509 = t104 * t320;
t508 = t104 * t328;
t507 = t125 * t344;
t506 = t126 * t344;
t505 = t128 * t344;
t504 = t129 * t344;
t503 = t130 * t344;
t502 = t131 * t344;
t501 = t132 * t344;
t500 = t133 * t344;
t260 = m(2) * rSges(2,1) + t329;
t193 = -m(3) * t569 + t260;
t247 = m(2) * rSges(2,2) - t550;
t153 = t193 * t300 - t247 * t298;
t499 = t153 * t294;
t196 = -m(3) * t570 + t260;
t154 = t196 * t316 - t247 * t310;
t498 = t154 * t294;
t197 = -m(3) * t571 + t260;
t155 = t197 * t318 - t247 * t312;
t497 = t155 * t294;
t198 = -m(3) * t572 + t260;
t156 = t198 * t320 - t247 * t314;
t496 = t156 * t294;
t495 = t158 * t272;
t494 = t158 * t276;
t493 = t160 * t273;
t492 = t160 * t277;
t491 = t161 * t274;
t490 = t161 * t278;
t489 = t162 * t275;
t488 = t162 * t279;
t173 = -t251 * t299 - t297 * t252;
t487 = t173 * t344;
t174 = -t251 * t315 - t309 * t252;
t486 = t174 * t344;
t175 = -t251 * t317 - t311 * t252;
t485 = t175 * t344;
t176 = -t251 * t319 - t313 * t252;
t484 = t176 * t344;
t482 = t215 * t296;
t480 = t216 * t296;
t478 = t217 * t296;
t476 = t218 * t296;
t246 = (t341 + t342) * m(3) + Icges(3,3);
t471 = t246 * t344;
t470 = t247 * t293;
t465 = rSges(3,2) * t270 * t294;
t463 = t294 * t295;
t462 = t294 * t297;
t459 = t294 * t300;
t458 = t294 * t309;
t456 = t294 * t311;
t454 = t294 * t313;
t451 = t294 * t316;
t449 = t294 * t318;
t447 = t294 * t320;
t429 = t296 * t344;
t428 = t298 * t299;
t427 = t310 * t315;
t426 = t312 * t317;
t425 = t314 * t319;
t419 = -t417 / 0.2e1 + Icges(3,4) / 0.2e1;
t418 = -0.2e1 * rSges(3,2) * pkin(2);
t411 = t157 * t549;
t410 = t166 * t549;
t409 = t167 * t549;
t408 = t168 * t549;
t407 = t297 * t514;
t406 = t309 * t512;
t405 = t311 * t510;
t404 = t313 * t508;
t403 = t272 * t506;
t402 = t276 * t506;
t401 = t273 * t502;
t400 = t277 * t502;
t399 = t274 * t501;
t398 = t278 * t501;
t397 = t275 * t500;
t396 = t279 * t500;
t269 = rSges(3,2) * t412;
t394 = (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,1) + Icges(2,3);
t389 = t219 * t293 + t536;
t387 = t220 * t293 + t536;
t385 = t221 * t293 + t536;
t383 = t222 * t293 + t536;
t248 = t321 ^ 2 + t341 + t345;
t265 = 0.2e1 * t264;
t281 = t563 * t564;
t127 = t245 * t287 + (-t265 * t297 + t281) * t299 + (t297 * t418 + t248) * m(3) + t394;
t41 = t407 - t562;
t13 = (((t296 * t90 + t98 * t459) * t544 + ((-t416 + t514) * t298 + pkin(2) * t518) * t460 + t41 * t296) * t98 + (t90 * t459 + (t287 * t296 - t428 * t462 - t296) * t98) * t562) * t137;
t189 = pkin(3) * t428 + t225;
t17 = t137 * t429 * t535 + (-t296 * t407 + (-t189 * t462 + t296 * (pkin(2) * t299 + t544)) * t90) / (t189 * t460 + (pkin(2) + t540) * t441) * t90;
t21 = (-t299 * t535 - (pkin(2) * t90 - t299 * t41) * t562) * t137;
t223 = t247 * t536;
t5 = -t21 * t499 - t127 * t13 - t173 * t17 - 0.4e1 * ((t557 * t297 + t556 * t299) * t90 + ((t297 * t558 + t269) * t299 + t555 + t568) * t98) * t90 + (-t193 * t581 + t219 * t470 + t223) * t300 - t298 * (-t389 * t193 - t247 * t581);
t224 = (t342 / 0.2e1 - t341 / 0.2e1) * m(3) - Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1;
t349 = t581 * t298 + t389 * t300;
t9 = -t21 * t554 - t173 * t13 - t246 * t17 + 0.2e1 * ((t224 * t297 + t269) * t299 + t419 + t568) * t97 + m(3) * (((t388 * t294 - t482) * rSges(3,1) + t349 * rSges(3,2)) * t299 + t297 * (t465 + (-t219 * t463 + t482) * rSges(3,2) + t349 * rSges(3,1)));
t373 = -t158 * t5 - t9 * t506;
t43 = t406 - t561;
t14 = (((t102 * t451 + t296 * t94) * t543 + ((-t415 + t512) * t310 + pkin(2) * t513) * t452 + t43 * t296) * t102 + (t94 * t451 + (t289 * t296 - t427 * t458 - t296) * t102) * t561) * t138;
t190 = pkin(3) * t427 + t227;
t18 = t138 * t429 * t522 + (-t296 * t406 + (-t190 * t458 + t296 * (pkin(2) * t315 + t543)) * t94) / (t190 * t452 + (pkin(2) + t539) * t438) * t94;
t22 = (-t315 * t522 - (pkin(2) * t94 - t315 * t43) * t561) * t138;
t348 = t582 * t310 + t387 * t316;
t10 = -t22 * t553 - t174 * t14 - t246 * t18 + 0.2e1 * ((t224 * t309 + t269) * t315 + t419 + t567) * t99 + m(3) * (((t386 * t294 - t480) * rSges(3,1) + t348 * rSges(3,2)) * t315 + t309 * (t465 + (-t220 * t463 + t480) * rSges(3,2) + t348 * rSges(3,1)));
t134 = t245 * t289 + (-t265 * t309 + t281) * t315 + (t309 * t418 + t248) * m(3) + t394;
t6 = -t22 * t498 - t134 * t14 - t174 * t18 - 0.4e1 * ((t557 * t309 + t556 * t315) * t94 + ((t309 * t558 + t269) * t315 + t555 + t567) * t102) * t94 + (-t196 * t582 + t220 * t470 + t223) * t316 - t310 * (-t387 * t196 - t247 * t582);
t372 = -t10 * t502 - t160 * t6;
t44 = t405 - t560;
t15 = (((t103 * t449 + t296 * t95) * t542 + ((-t414 + t510) * t312 + pkin(2) * t511) * t450 + t44 * t296) * t103 + (t95 * t449 + (t290 * t296 - t426 * t456 - t296) * t103) * t560) * t139;
t191 = pkin(3) * t426 + t228;
t19 = t139 * t429 * t521 + (-t296 * t405 + (-t191 * t456 + t296 * (pkin(2) * t317 + t542)) * t95) / (t191 * t450 + (pkin(2) + t538) * t436) * t95;
t23 = (-t317 * t521 - (pkin(2) * t95 - t317 * t44) * t560) * t139;
t347 = t583 * t312 + t385 * t318;
t11 = -t23 * t552 - t175 * t15 - t246 * t19 + 0.2e1 * ((t224 * t311 + t269) * t317 + t419 + t566) * t100 + m(3) * (((t384 * t294 - t478) * rSges(3,1) + t347 * rSges(3,2)) * t317 + t311 * (t465 + (-t221 * t463 + t478) * rSges(3,2) + t347 * rSges(3,1)));
t135 = t245 * t290 + (-t265 * t311 + t281) * t317 + (t311 * t418 + t248) * m(3) + t394;
t7 = -t23 * t497 - t135 * t15 - t175 * t19 - 0.4e1 * ((t557 * t311 + t556 * t317) * t95 + ((t311 * t558 + t269) * t317 + t555 + t566) * t103) * t95 + (-t197 * t583 + t221 * t470 + t223) * t318 - t312 * (-t385 * t197 - t247 * t583);
t371 = -t11 * t501 - t161 * t7;
t45 = t404 - t559;
t16 = (((t104 * t447 + t296 * t96) * t541 + ((-t413 + t508) * t314 + pkin(2) * t509) * t448 + t45 * t296) * t104 + (t96 * t447 + (t291 * t296 - t425 * t454 - t296) * t104) * t559) * t140;
t192 = pkin(3) * t425 + t229;
t20 = t140 * t429 * t520 + (-t296 * t404 + (-t192 * t454 + t296 * (pkin(2) * t319 + t541)) * t96) / (t192 * t448 + (pkin(2) + t537) * t434) * t96;
t24 = (-t319 * t520 - (pkin(2) * t96 - t319 * t45) * t559) * t140;
t346 = t584 * t314 + t383 * t320;
t12 = -t24 * t551 - t176 * t16 - t246 * t20 + 0.2e1 * ((t224 * t313 + t269) * t319 + t419 + t565) * t101 + m(3) * (((t382 * t294 - t476) * rSges(3,1) + t346 * rSges(3,2)) * t319 + t313 * (t465 + (-t222 * t463 + t476) * rSges(3,2) + t346 * rSges(3,1)));
t136 = t245 * t291 + (-t265 * t313 + t281) * t319 + (t313 * t418 + t248) * m(3) + t394;
t8 = -t24 * t496 - t136 * t16 - t176 * t20 - 0.4e1 * ((t557 * t313 + t556 * t319) * t96 + ((t313 * t558 + t269) * t319 + t555 + t565) * t104) * t96 + (-t198 * t584 + t222 * t470 + t223) * t320 - t314 * (-t383 * t198 - t247 * t584);
t370 = -t12 * t500 - t162 * t8;
t369 = pkin(3) * t462 - t225 * t296;
t368 = pkin(3) * t458 - t227 * t296;
t367 = pkin(3) * t456 - t228 * t296;
t366 = pkin(3) * t454 - t229 * t296;
t365 = t126 * t487 + t127 * t158;
t364 = t126 * t471 + t158 * t173;
t363 = t131 * t486 + t134 * t160;
t362 = t131 * t471 + t160 * t174;
t361 = t132 * t485 + t135 * t161;
t360 = t132 * t471 + t161 * t175;
t359 = t133 * t484 + t136 * t162;
t358 = t133 * t471 + t162 * t176;
t357 = t137 * (t207 * t276 + t211 * t272);
t356 = t138 * (t208 * t277 + t212 * t273);
t355 = t139 * (t209 * t278 + t213 * t274);
t354 = t140 * (t210 * t279 + t214 * t275);
t353 = t126 * t411 + t158 * t499;
t352 = t131 * t410 + t160 * t498;
t351 = t132 * t409 + t161 * t497;
t350 = t133 * t408 + t162 * t496;
t332 = rSges(4,1);
t331 = rSges(4,2);
t308 = xDDP(1);
t307 = xDDP(2);
t306 = xDDP(3);
t305 = xDDP(4);
t292 = t324 ^ 2;
t288 = m(1) + m(2) + m(3);
t206 = -t285 * t331 + t332 * t286;
t205 = t285 * t332 + t331 * t286;
t152 = -t210 * t305 - t214 * t292 + t308;
t151 = -t209 * t305 - t213 * t292 + t308;
t150 = -t208 * t305 - t212 * t292 + t308;
t149 = -t207 * t305 - t211 * t292 + t308;
t148 = -t210 * t292 + t214 * t305 + t307;
t147 = -t209 * t292 + t213 * t305 + t307;
t146 = -t208 * t292 + t212 * t305 + t307;
t145 = -t207 * t292 + t211 * t305 + t307;
t144 = -t232 * t293 + t366 * t295;
t143 = -t231 * t293 + t367 * t295;
t142 = -t230 * t293 + t368 * t295;
t141 = -t226 * t293 + t369 * t295;
t124 = -t201 * t541 + t232 * t442 + (pkin(2) * t454 + t366 * t319) * t293;
t123 = -t200 * t542 + t231 * t443 + (pkin(2) * t456 + t367 * t317) * t293;
t122 = -t199 * t543 + t230 * t444 + (pkin(2) * t458 + t368 * t315) * t293;
t121 = -t194 * t544 + t226 * t445 + (pkin(2) * t462 + t369 * t299) * t293;
t120 = -(t204 * t275 - t279 * t453) * t541 + (t144 * t275 + t172 * t279) * t319 + (t275 * t463 + t279 * t296) * t545;
t119 = -(t203 * t274 - t278 * t455) * t542 + (t143 * t274 + t171 * t278) * t317 + (t274 * t463 + t278 * t296) * t546;
t118 = -(t202 * t273 - t277 * t457) * t543 + (t142 * t273 + t170 * t277) * t315 + (t273 * t463 + t277 * t296) * t547;
t117 = (t204 * t279 + t275 * t453) * t541 + (-t144 * t279 + t172 * t275) * t319 + (t275 * t296 - t279 * t463) * t545;
t116 = (t203 * t278 + t274 * t455) * t542 + (-t143 * t278 + t171 * t274) * t317 + (t274 * t296 - t278 * t463) * t546;
t115 = (t202 * t277 + t273 * t457) * t543 + (-t142 * t277 + t170 * t273) * t315 + (t273 * t296 - t277 * t463) * t547;
t114 = -(t195 * t272 - t276 * t461) * t544 + (t141 * t272 + t169 * t276) * t299 + (t272 * t463 + t276 * t296) * t548;
t113 = (t195 * t276 + t272 * t461) * t544 + (-t141 * t276 + t169 * t272) * t299 + (t272 * t296 - t276 * t463) * t548;
t112 = t162 * t354;
t111 = t161 * t355;
t110 = t160 * t356;
t109 = t158 * t357;
t108 = t354 * t500;
t107 = t355 * t501;
t106 = t356 * t502;
t105 = t357 * t506;
t93 = t96 ^ 2;
t92 = t95 ^ 2;
t91 = t94 ^ 2;
t89 = (t124 * t551 + t130 * t471 + t165 * t176) * t140;
t88 = (t123 * t552 + t129 * t471 + t164 * t175) * t139;
t87 = (t122 * t553 + t128 * t471 + t163 * t174) * t138;
t86 = t90 ^ 2;
t85 = (t124 * t288 + t130 * t408 + t165 * t496) * t140;
t84 = (t123 * t288 + t129 * t409 + t164 * t497) * t139;
t83 = (t122 * t288 + t128 * t410 + t163 * t498) * t138;
t82 = (t121 * t554 + t125 * t471 + t159 * t173) * t137;
t81 = (t121 * t288 + t125 * t411 + t159 * t499) * t137;
t80 = (-t117 * t210 + t120 * t214) * t140;
t79 = (-t116 * t209 + t119 * t213) * t139;
t78 = (-t115 * t208 + t118 * t212) * t138;
t77 = (t124 * t496 + t130 * t484 + t136 * t165) * t140;
t76 = (t123 * t497 + t129 * t485 + t135 * t164) * t139;
t75 = (t122 * t498 + t128 * t486 + t134 * t163) * t138;
t74 = (-t113 * t207 + t114 * t211) * t137;
t73 = (t121 * t499 + t125 * t487 + t127 * t159) * t137;
t72 = (t120 * t551 + t358 * t275) * t140;
t71 = (t119 * t552 + t360 * t274) * t139;
t70 = (t118 * t553 + t362 * t273) * t138;
t69 = (t117 * t551 - t358 * t279) * t140;
t68 = (t116 * t552 - t360 * t278) * t139;
t67 = (t115 * t553 - t362 * t277) * t138;
t66 = (t120 * t288 + t350 * t275) * t140;
t65 = (t119 * t288 + t351 * t274) * t139;
t64 = (t118 * t288 + t352 * t273) * t138;
t63 = (t117 * t288 - t350 * t279) * t140;
t62 = (t116 * t288 - t351 * t278) * t139;
t61 = (t115 * t288 - t352 * t277) * t138;
t60 = (t114 * t554 + t364 * t272) * t137;
t59 = (t113 * t554 - t364 * t276) * t137;
t58 = (t114 * t288 + t353 * t272) * t137;
t57 = (t113 * t288 - t353 * t276) * t137;
t56 = (t120 * t496 + t359 * t275) * t140;
t55 = (t119 * t497 + t361 * t274) * t139;
t54 = (t118 * t498 + t363 * t273) * t138;
t53 = (t117 * t496 - t359 * t279) * t140;
t52 = (t116 * t497 - t361 * t278) * t139;
t51 = (t115 * t498 - t363 * t277) * t138;
t50 = (t114 * t499 + t365 * t272) * t137;
t49 = (t113 * t499 - t365 * t276) * t137;
t40 = t108 * t246 + t112 * t176 + t80 * t551;
t39 = t107 * t246 + t111 * t175 + t79 * t552;
t38 = t106 * t246 + t110 * t174 + t78 * t553;
t37 = t108 * t551 + t112 * t496 + t288 * t80;
t36 = t107 * t552 + t111 * t497 + t288 * t79;
t35 = t106 * t553 + t110 * t498 + t288 * t78;
t34 = t105 * t246 + t109 * t173 + t74 * t554;
t33 = t105 * t554 + t109 * t499 + t288 * t74;
t32 = t108 * t176 + t112 * t136 + t80 * t496;
t31 = t107 * t175 + t111 * t135 + t79 * t497;
t30 = t106 * t174 + t110 * t134 + t78 * t498;
t29 = t105 * t173 + t109 * t127 + t74 * t499;
t4 = (-t156 * t16 + (-t247 * t320 - t260 * t314) * t101) * t294 + (-t24 - t218) * t288 + (-t168 * t20 + (-0.2e1 * (rSges(3,1) * t515 + t96 * t523) * t509 + t572 * t314 * (t101 + t93)) * t294 - t93 * t296 * t244) * m(3);
t3 = (-t15 * t155 + (-t247 * t318 - t260 * t312) * t100) * t294 + (-t23 - t217) * t288 + (-t167 * t19 + (-0.2e1 * (rSges(3,1) * t516 + t95 * t524) * t511 + t571 * t312 * (t100 + t92)) * t294 - t92 * t296 * t243) * m(3);
t2 = (-t154 * t14 + (-t247 * t316 - t260 * t310) * t99) * t294 + (-t22 - t216) * t288 + (-t166 * t18 + (-0.2e1 * (rSges(3,1) * t517 + t94 * t525) * t513 + t570 * t310 * (t99 + t91)) * t294 - t91 * t296 * t242) * m(3);
t1 = (-t13 * t153 + (-t247 * t300 - t260 * t298) * t97) * t294 + (-t21 - t215) * t288 + (-t157 * t17 + (-0.2e1 * (rSges(3,1) * t519 + t90 * t529) * t518 + t569 * t298 * (t97 + t86)) * t294 - t86 * t296 * t238) * m(3);
t25 = [(-t205 * t305 - t292 * t206 - g(1) + t308) * m(4) + ((t120 * t63 + t397 * t69 + t53 * t489) * t148 + (t124 * t63 + t165 * t53 + t69 * t503) * t306 + (t63 * t152 + t4) * t117 + ((-t162 * t53 - t69 * t500) * t152 + t370) * t279) * t140 + ((t119 * t62 + t399 * t68 + t52 * t491) * t147 + (t123 * t62 + t164 * t52 + t68 * t504) * t306 + (t62 * t151 + t3) * t116 + ((-t161 * t52 - t68 * t501) * t151 + t371) * t278) * t139 + ((t122 * t61 + t163 * t51 + t67 * t505) * t306 + (t118 * t61 + t401 * t67 + t51 * t493) * t146 + (t61 * t150 + t2) * t115 + ((-t160 * t51 - t67 * t502) * t150 + t372) * t277) * t138 + ((t121 * t57 + t159 * t49 + t59 * t507) * t306 + (t114 * t57 + t403 * t59 + t49 * t495) * t145 + (t57 * t149 + t1) * t113 + ((-t158 * t49 - t59 * t506) * t149 + t373) * t276) * t137; (-t292 * t205 + t206 * t305 - g(2) + t307) * m(4) + ((t117 * t66 - t396 * t72 - t56 * t488) * t152 + (t124 * t66 + t165 * t56 + t72 * t503) * t306 + (t66 * t148 + t4) * t120 + ((t162 * t56 + t72 * t500) * t148 - t370) * t275) * t140 + ((t116 * t65 - t398 * t71 - t55 * t490) * t151 + (t123 * t65 + t164 * t55 + t71 * t504) * t306 + (t65 * t147 + t3) * t119 + ((t161 * t55 + t71 * t501) * t147 - t371) * t274) * t139 + ((t115 * t64 - t400 * t70 - t54 * t492) * t150 + (t122 * t64 + t163 * t54 + t70 * t505) * t306 + (t64 * t146 + t2) * t118 + ((t160 * t54 + t70 * t502) * t146 - t372) * t273) * t138 + ((t113 * t58 - t402 * t60 - t50 * t494) * t149 + (t121 * t58 + t159 * t50 + t60 * t507) * t306 + (t58 * t145 + t1) * t114 + ((t158 * t50 + t60 * t506) * t145 - t373) * t272) * t137; (t306 - g(3)) * m(4) + (t165 * t8 + t124 * t4 + t12 * t503 + (t120 * t85 + t397 * t89 + t77 * t489) * t148 + (t117 * t85 - t396 * t89 - t77 * t488) * t152 + (t124 * t85 + t165 * t77 + t89 * t503) * t306) * t140 + (t164 * t7 + t123 * t3 + t11 * t504 + (t116 * t84 - t398 * t88 - t76 * t490) * t151 + (t123 * t84 + t164 * t76 + t88 * t504) * t306 + (t119 * t84 + t399 * t88 + t76 * t491) * t147) * t139 + (t163 * t6 + t122 * t2 + t10 * t505 + (t115 * t83 - t400 * t87 - t75 * t492) * t150 + (t122 * t83 + t163 * t75 + t87 * t505) * t306 + (t118 * t83 + t401 * t87 + t75 * t493) * t146) * t138 + (t159 * t5 + t121 * t1 + t9 * t507 + (t113 * t81 - t402 * t82 - t73 * t494) * t149 + (t121 * t81 + t159 * t73 + t82 * t507) * t306 + (t114 * t81 + t403 * t82 + t73 * t495) * t145) * t137; Icges(4,3) * t305 + t74 * t1 + t106 * t10 + t105 * t9 + t107 * t11 + t108 * t12 + t109 * t5 + t110 * t6 + t111 * t7 + t112 * t8 + t78 * t2 + t79 * t3 + t80 * t4 + (t206 * t307 - t205 * t308 + (g(1) * t331 - g(2) * t332) * t286 + t285 * (g(1) * t332 + g(2) * t331) + (t331 ^ 2 + t332 ^ 2) * t305) * m(4) + ((t120 * t37 + t32 * t489 + t397 * t40) * t148 + (t117 * t37 - t32 * t488 - t396 * t40) * t152 + (t124 * t37 + t165 * t32 + t40 * t503) * t306) * t140 + ((t116 * t36 - t31 * t490 - t39 * t398) * t151 + (t123 * t36 + t164 * t31 + t39 * t504) * t306 + (t119 * t36 + t31 * t491 + t39 * t399) * t147) * t139 + ((t115 * t35 - t30 * t492 - t38 * t400) * t150 + (t122 * t35 + t163 * t30 + t38 * t505) * t306 + (t118 * t35 + t30 * t493 + t38 * t401) * t146) * t138 + ((t113 * t33 - t29 * t494 - t34 * t402) * t149 + (t121 * t33 + t159 * t29 + t34 * t507) * t306 + (t114 * t33 + t29 * t495 + t34 * t403) * t145) * t137;];
tauX  = t25;
