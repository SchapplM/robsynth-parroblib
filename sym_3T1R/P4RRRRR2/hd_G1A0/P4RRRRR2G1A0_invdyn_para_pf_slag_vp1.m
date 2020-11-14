% Calculate vector of inverse dynamics forces for parallel robot
% P4RRRRR2G1A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2020-08-07 17:26
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P4RRRRR2G1A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4RRRRR2G1A0_invdyn_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4RRRRR2G1A0_invdyn_para_pf_slag_vp1: xDP has to be [4x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [4 1]), ...
  'P4RRRRR2G1A0_invdyn_para_pf_slag_vp1: xDDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4RRRRR2G1A0_invdyn_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4RRRRR2G1A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4RRRRR2G1A0_invdyn_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4RRRRR2G1A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4RRRRR2G1A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P4RRRRR2G1A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4RRRRR2G1A0_invdyn_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4RRRRR2G1A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 17:24:02
% EndTime: 2020-08-07 17:24:15
% DurationCPUTime: 12.72s
% Computational Cost: add. (35429->727), mult. (38047->1269), div. (13140->23), fcn. (24966->114), ass. (0->554)
t657 = 2 * pkin(1);
t656 = -2 * rSges(3,1);
t655 = 2 * rSges(3,2);
t631 = m(3) * rSges(3,2);
t253 = rSges(3,1) * t631 - Icges(3,4);
t380 = 2 * qJ(3,1);
t312 = sin(t380);
t315 = cos(t380);
t639 = rSges(3,2) ^ 2;
t640 = rSges(3,1) ^ 2;
t223 = (-t639 + t640) * m(3) - Icges(3,1) + Icges(3,2);
t625 = t223 / 0.2e1;
t654 = -t253 * t312 + t315 * t625;
t379 = 2 * qJ(3,2);
t311 = sin(t379);
t314 = cos(t379);
t653 = -t253 * t311 + t314 * t625;
t378 = 2 * qJ(3,3);
t310 = sin(t378);
t313 = cos(t378);
t652 = -t253 * t310 + t313 * t625;
t377 = 2 * qJ(3,4);
t300 = sin(t377);
t301 = cos(t377);
t651 = -t253 * t300 + t301 * t625;
t630 = m(3) * rSges(3,3);
t249 = m(2) * rSges(2,2) - t630;
t359 = sin(qJ(2,1));
t365 = cos(qJ(2,1));
t369 = m(2) * rSges(2,1);
t358 = sin(qJ(3,1));
t364 = cos(qJ(3,1));
t646 = -rSges(3,1) * t364 + rSges(3,2) * t358;
t650 = -t249 * t359 - (m(3) * t646 - t369) * t365;
t357 = sin(qJ(2,2));
t363 = cos(qJ(2,2));
t356 = sin(qJ(3,2));
t362 = cos(qJ(3,2));
t645 = -rSges(3,1) * t362 + rSges(3,2) * t356;
t649 = -t249 * t357 - (m(3) * t645 - t369) * t363;
t355 = sin(qJ(2,3));
t361 = cos(qJ(2,3));
t354 = sin(qJ(3,3));
t360 = cos(qJ(3,3));
t644 = -rSges(3,1) * t360 + rSges(3,2) * t354;
t648 = -t249 * t355 - (m(3) * t644 - t369) * t361;
t347 = sin(qJ(2,4));
t349 = cos(qJ(2,4));
t346 = sin(qJ(3,4));
t348 = cos(qJ(3,4));
t643 = -rSges(3,1) * t348 + rSges(3,2) * t346;
t647 = -t249 * t347 - (m(3) * t643 - t369) * t349;
t642 = -2 * pkin(1);
t638 = -4 * t253;
t637 = 2 * t253;
t636 = -m(3) / 0.2e1;
t635 = m(3) / 0.2e1;
t634 = m(3) * pkin(1);
t633 = m(1) * rSges(1,2);
t632 = m(3) * rSges(3,1);
t374 = xP(4);
t298 = sin(t374);
t299 = cos(t374);
t381 = koppelP(4,2);
t385 = koppelP(4,1);
t183 = -t298 * t381 + t299 * t385;
t370 = xDP(4);
t372 = xDP(2);
t145 = t183 * t370 + t372;
t179 = t298 * t385 + t299 * t381;
t373 = xDP(1);
t149 = -t179 * t370 + t373;
t371 = xDP(3);
t307 = 0.1e1 / t348;
t305 = 0.1e1 / t347;
t393 = 1 / pkin(1);
t530 = t305 * t393;
t463 = t346 * t530;
t423 = t307 * t463;
t304 = qJ(1,4) + qJ(2,4);
t342 = legFrame(4,3);
t245 = t342 + t304;
t229 = cos(t245);
t470 = t229 * t530;
t228 = sin(t245);
t471 = t228 * t530;
t82 = t145 * t471 + t149 * t470 + t371 * t423;
t81 = t82 ^ 2;
t629 = pkin(1) * t81;
t382 = koppelP(3,2);
t386 = koppelP(3,1);
t184 = -t298 * t382 + t299 * t386;
t146 = t184 * t370 + t372;
t180 = t298 * t386 + t299 * t382;
t150 = -t180 * t370 + t373;
t320 = 0.1e1 / t360;
t316 = 0.1e1 / t355;
t526 = t316 * t393;
t462 = t354 * t526;
t422 = t320 * t462;
t343 = legFrame(3,3);
t263 = qJ(1,3) + t343;
t246 = qJ(2,3) + t263;
t236 = cos(t246);
t466 = t236 * t526;
t233 = sin(t246);
t469 = t233 * t526;
t87 = t146 * t469 + t150 * t466 + t371 * t422;
t84 = t87 ^ 2;
t628 = pkin(1) * t84;
t383 = koppelP(2,2);
t387 = koppelP(2,1);
t185 = -t298 * t383 + t299 * t387;
t147 = t185 * t370 + t372;
t181 = t298 * t387 + t299 * t383;
t151 = -t181 * t370 + t373;
t324 = 0.1e1 / t362;
t317 = 0.1e1 / t357;
t524 = t317 * t393;
t461 = t356 * t524;
t421 = t324 * t461;
t338 = qJ(2,2) + qJ(1,2);
t344 = legFrame(2,3);
t247 = t344 + t338;
t237 = cos(t247);
t465 = t237 * t524;
t234 = sin(t247);
t468 = t234 * t524;
t88 = t147 * t468 + t151 * t465 + t371 * t421;
t85 = t88 ^ 2;
t627 = pkin(1) * t85;
t384 = koppelP(1,2);
t388 = koppelP(1,1);
t186 = -t298 * t384 + t299 * t388;
t148 = t186 * t370 + t372;
t182 = t298 * t388 + t299 * t384;
t152 = -t182 * t370 + t373;
t328 = 0.1e1 / t364;
t318 = 0.1e1 / t359;
t522 = t318 * t393;
t460 = t358 * t522;
t420 = t328 * t460;
t345 = legFrame(1,3);
t265 = qJ(1,1) + t345;
t248 = qJ(2,1) + t265;
t238 = cos(t248);
t464 = t238 * t522;
t235 = sin(t248);
t467 = t235 * t522;
t89 = t148 * t467 + t152 * t464 + t371 * t420;
t86 = t89 ^ 2;
t626 = pkin(1) * t86;
t624 = pkin(1) * t347;
t623 = pkin(1) * t349;
t622 = pkin(1) * t355;
t621 = pkin(1) * t357;
t620 = pkin(1) * t359;
t619 = pkin(1) * t361;
t618 = pkin(1) * t363;
t617 = pkin(1) * t365;
t306 = t348 ^ 2;
t616 = pkin(2) * t306;
t319 = t360 ^ 2;
t615 = pkin(2) * t319;
t323 = t362 ^ 2;
t614 = pkin(2) * t323;
t327 = t364 ^ 2;
t613 = pkin(2) * t327;
t612 = pkin(2) * t348;
t611 = pkin(2) * t360;
t610 = pkin(2) * t362;
t609 = pkin(2) * t364;
t331 = t370 ^ 2;
t608 = m(4) * t331;
t350 = xDDP(4);
t607 = m(4) * t350;
t352 = xDDP(2);
t606 = m(4) * t352;
t353 = xDDP(1);
t605 = m(4) * t353;
t308 = 0.1e1 / t348 ^ 2;
t390 = 0.1e1 / pkin(2);
t528 = t308 * t390;
t475 = (t612 + t623) * t528;
t531 = t305 * t346;
t415 = t475 * t531;
t407 = t393 * t415;
t403 = t371 * t407;
t302 = qJ(2,4) + qJ(3,4);
t261 = qJ(1,4) + t302;
t230 = t342 + t261;
t303 = qJ(2,4) - qJ(3,4);
t262 = qJ(1,4) + t303;
t231 = t342 + t262;
t254 = qJ(1,4) + t342;
t137 = sin(t254) * t642 + (-sin(t231) - sin(t230)) * pkin(2);
t255 = sin(t302);
t256 = sin(t303);
t203 = 0.1e1 / (t255 + t256);
t507 = t390 * t393;
t479 = t203 * t507;
t431 = t137 * t479;
t93 = t145 * t431;
t138 = cos(t254) * t642 + (-cos(t230) - cos(t231)) * pkin(2);
t430 = t138 * t479;
t94 = t149 * t430;
t34 = t94 / 0.2e1 + t93 / 0.2e1 - t403 / 0.2e1 + t82;
t69 = t93 + t94 - t403;
t604 = t69 * t34;
t321 = 0.1e1 / t360 ^ 2;
t520 = t321 * t390;
t474 = (t611 + t619) * t520;
t527 = t316 * t354;
t414 = t474 * t527;
t406 = t393 * t414;
t402 = t371 * t406;
t239 = qJ(3,3) + t246;
t240 = -qJ(3,3) + t246;
t139 = sin(t263) * t642 + (-sin(t240) - sin(t239)) * pkin(2);
t333 = qJ(2,3) + qJ(3,3);
t266 = sin(t333);
t334 = qJ(2,3) - qJ(3,3);
t267 = sin(t334);
t205 = 0.1e1 / (t266 + t267);
t478 = t205 * t507;
t429 = t139 * t478;
t95 = t146 * t429;
t142 = cos(t263) * t642 + (-cos(t239) - cos(t240)) * pkin(2);
t426 = t142 * t478;
t98 = t150 * t426;
t42 = t98 / 0.2e1 + t95 / 0.2e1 - t402 / 0.2e1 + t87;
t70 = t95 + t98 - t402;
t603 = t70 * t42;
t325 = 0.1e1 / t362 ^ 2;
t518 = t325 * t390;
t473 = (t610 + t618) * t518;
t525 = t317 * t356;
t413 = t473 * t525;
t405 = t393 * t413;
t401 = t371 * t405;
t336 = qJ(2,2) + qJ(3,2);
t294 = qJ(1,2) + t336;
t241 = t344 + t294;
t337 = qJ(2,2) - qJ(3,2);
t295 = qJ(1,2) + t337;
t242 = t344 + t295;
t264 = qJ(1,2) + t344;
t140 = sin(t264) * t642 + (-sin(t242) - sin(t241)) * pkin(2);
t269 = sin(t336);
t270 = sin(t337);
t206 = 0.1e1 / (t269 + t270);
t477 = t206 * t507;
t428 = t140 * t477;
t96 = t147 * t428;
t143 = cos(t264) * t642 + (-cos(t241) - cos(t242)) * pkin(2);
t425 = t143 * t477;
t99 = t151 * t425;
t43 = t99 / 0.2e1 + t96 / 0.2e1 - t401 / 0.2e1 + t88;
t71 = t96 + t99 - t401;
t602 = t71 * t43;
t243 = qJ(3,1) + t248;
t244 = -qJ(3,1) + t248;
t144 = cos(t265) * t642 + (-cos(t243) - cos(t244)) * pkin(2);
t339 = qJ(2,1) + qJ(3,1);
t272 = sin(t339);
t340 = qJ(2,1) - qJ(3,1);
t273 = sin(t340);
t207 = 0.1e1 / (t272 + t273);
t476 = t207 * t507;
t424 = t144 * t476;
t100 = t152 * t424;
t329 = 0.1e1 / t364 ^ 2;
t516 = t329 * t390;
t472 = (t609 + t617) * t516;
t523 = t318 * t358;
t412 = t472 * t523;
t404 = t393 * t412;
t400 = t371 * t404;
t141 = sin(t265) * t642 + (-sin(t244) - sin(t243)) * pkin(2);
t427 = t141 * t476;
t97 = t148 * t427;
t44 = t100 / 0.2e1 + t97 / 0.2e1 - t400 / 0.2e1 + t89;
t72 = t100 + t97 - t400;
t601 = t72 * t44;
t284 = sin(t342);
t288 = cos(t342);
t191 = -g(1) * t284 + g(2) * t288;
t600 = rSges(3,1) * t191;
t285 = sin(t343);
t289 = cos(t343);
t192 = -g(1) * t285 + g(2) * t289;
t599 = rSges(3,1) * t192;
t286 = sin(t344);
t290 = cos(t344);
t193 = -g(1) * t286 + g(2) * t290;
t598 = rSges(3,1) * t193;
t287 = sin(t345);
t291 = cos(t345);
t194 = -g(1) * t287 + g(2) * t291;
t597 = rSges(3,1) * t194;
t195 = g(1) * t288 + g(2) * t284;
t596 = rSges(3,1) * t195;
t196 = g(1) * t289 + g(2) * t285;
t595 = rSges(3,1) * t196;
t197 = g(1) * t290 + g(2) * t286;
t594 = rSges(3,1) * t197;
t198 = g(1) * t291 + g(2) * t287;
t593 = rSges(3,1) * t198;
t584 = rSges(3,3) * t191;
t583 = rSges(3,3) * t192;
t582 = rSges(3,3) * t193;
t581 = rSges(3,3) * t194;
t580 = rSges(3,3) * t195;
t579 = rSges(3,3) * t196;
t578 = rSges(3,3) * t197;
t577 = rSges(3,3) * t198;
t576 = (-t179 * t331 + t183 * t350 + t352) * t393;
t575 = (-t180 * t331 + t184 * t350 + t352) * t393;
t574 = (-t181 * t331 + t185 * t350 + t352) * t393;
t573 = (-t182 * t331 + t186 * t350 + t352) * t393;
t572 = (-t179 * t350 - t183 * t331 + t353) * t393;
t571 = (-t180 * t350 - t184 * t331 + t353) * t393;
t570 = (-t181 * t350 - t185 * t331 + t353) * t393;
t569 = (-t182 * t350 - t186 * t331 + t353) * t393;
t564 = t203 * t390;
t563 = t205 * t390;
t562 = t206 * t390;
t561 = t207 * t390;
t560 = t223 * t300;
t559 = t223 * t310;
t558 = t223 * t311;
t557 = t223 * t312;
t556 = t228 * t305;
t555 = t229 * t305;
t553 = t233 * t316;
t552 = t234 * t317;
t551 = t235 * t318;
t550 = t236 * t316;
t549 = t237 * t317;
t548 = t238 * t318;
t251 = -rSges(3,2) * t630 + Icges(3,6);
t543 = t251 * t308;
t542 = t251 * t321;
t541 = t251 * t325;
t540 = t251 * t329;
t538 = t253 * t301;
t534 = t253 * t313;
t533 = t253 * t314;
t532 = t253 * t315;
t529 = t307 * t390;
t521 = t320 * t390;
t519 = t324 * t390;
t517 = t328 * t390;
t332 = t371 ^ 2;
t515 = t332 * t390;
t514 = t332 / pkin(2) ^ 2;
t513 = t346 * t347;
t351 = xDDP(3);
t512 = t351 * t393;
t511 = t354 * t355;
t510 = t356 * t357;
t509 = t358 * t359;
t508 = t371 * t390;
t505 = pkin(1) * t632;
t504 = pkin(1) * t631;
t503 = t634 / 0.2e1;
t502 = pkin(1) * t508;
t501 = t349 * t616;
t500 = t361 * t615;
t499 = t363 * t614;
t498 = t365 * t613;
t36 = t69 + t82;
t497 = t36 * t508;
t48 = t70 + t87;
t496 = t48 * t508;
t49 = t71 + t88;
t495 = t49 * t508;
t50 = t72 + t89;
t494 = t50 * t508;
t493 = -2 * t505;
t492 = t639 + t640;
t491 = t137 * t564;
t490 = t138 * t564;
t489 = t139 * t563;
t488 = t140 * t562;
t487 = t141 * t561;
t486 = t142 * t563;
t485 = t143 * t562;
t484 = t144 * t561;
t252 = rSges(3,1) * t630 - Icges(3,5);
t153 = t251 * t348 - t252 * t346;
t483 = t153 * t564;
t154 = t251 * t360 - t252 * t354;
t482 = t154 * t563;
t155 = t251 * t362 - t252 * t356;
t481 = t155 * t562;
t156 = t251 * t364 - t252 * t358;
t480 = t156 * t561;
t459 = t346 * t514;
t458 = t354 * t514;
t457 = t356 * t514;
t456 = t358 * t514;
t455 = t371 * t513;
t454 = t371 * t511;
t453 = t371 * t510;
t452 = t371 * t509;
t451 = t252 * t514;
t450 = m(1) * rSges(1,1) + m(2) * pkin(1) + t634;
t445 = t514 / 0.2e1;
t444 = rSges(3,3) + t624;
t443 = rSges(3,3) + t622;
t442 = rSges(3,3) + t621;
t441 = rSges(3,3) + t620;
t440 = m(3) * (t308 * t445 + t604) * t624;
t439 = m(3) * (t321 * t445 + t603) * t622;
t438 = m(3) * (t325 * t445 + t602) * t621;
t437 = m(3) * (t329 * t445 + t601) * t620;
t436 = -t505 / 0.2e1;
t411 = t249 * t349 + t347 * t369;
t410 = t249 * t361 + t355 * t369;
t409 = t249 * t363 + t357 * t369;
t408 = t249 * t365 + t359 * t369;
t399 = Icges(2,3) + ((rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2)) + Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1 + (2 * rSges(3,3) ^ 2 + t492) * t635;
t105 = t399 + t651;
t106 = t399 + t652;
t107 = t399 + t653;
t108 = t399 + t654;
t392 = pkin(1) ^ 2;
t398 = Icges(1,3) + ((m(3) + m(2)) * t392) + ((rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1)) + t399;
t397 = t364 * t327;
t396 = t362 * t323;
t395 = t360 * t319;
t394 = t348 * t306;
t389 = pkin(2) ^ 2;
t376 = rSges(4,1);
t375 = rSges(4,2);
t341 = qJ(1,1) + qJ(2,1);
t335 = qJ(1,3) + qJ(2,3);
t330 = 0.1e1 / t397;
t326 = 0.1e1 / t396;
t322 = 0.1e1 / t395;
t309 = 0.1e1 / t394;
t297 = qJ(1,1) + t340;
t296 = qJ(1,1) + t339;
t293 = qJ(1,3) + t334;
t292 = qJ(1,3) + t333;
t283 = cos(t341);
t282 = cos(t340);
t281 = cos(t339);
t280 = cos(t338);
t279 = cos(t337);
t278 = cos(t336);
t277 = cos(t335);
t276 = cos(t334);
t275 = cos(t333);
t274 = sin(t341);
t271 = sin(t338);
t268 = sin(t335);
t260 = cos(t304);
t259 = cos(t303);
t258 = cos(t302);
t257 = sin(t304);
t232 = t492 * m(3) + Icges(3,3);
t174 = rSges(3,2) * t198;
t173 = rSges(3,2) * t197;
t172 = rSges(3,2) * t196;
t171 = rSges(3,2) * t195;
t170 = rSges(3,2) * t194;
t169 = rSges(3,2) * t193;
t168 = rSges(3,2) * t192;
t167 = rSges(3,2) * t191;
t163 = -t298 * t375 + t299 * t376;
t162 = t298 * t376 + t299 * t375;
t136 = (-t441 * t631 + Icges(3,6)) * t364 - t358 * (t441 * t632 - Icges(3,5));
t135 = (-t442 * t631 + Icges(3,6)) * t362 - t356 * (t442 * t632 - Icges(3,5));
t134 = (-t443 * t631 + Icges(3,6)) * t360 - t354 * (t443 * t632 - Icges(3,5));
t133 = m(2) * (rSges(2,1) * t198 + rSges(2,2) * t194);
t132 = m(2) * (rSges(2,1) * t197 + rSges(2,2) * t193);
t131 = m(2) * (rSges(2,1) * t196 + rSges(2,2) * t192);
t130 = m(2) * (rSges(2,1) * t195 + rSges(2,2) * t191);
t129 = m(2) * (-rSges(2,1) * t194 + rSges(2,2) * t198);
t128 = m(2) * (-rSges(2,1) * t193 + rSges(2,2) * t197);
t127 = m(2) * (-rSges(2,1) * t192 + rSges(2,2) * t196);
t126 = m(2) * (-rSges(2,1) * t191 + rSges(2,2) * t195);
t125 = (-t444 * t631 + Icges(3,6)) * t348 - t346 * (t444 * t632 - Icges(3,5));
t104 = (-t182 * t238 + t186 * t235) * t522;
t103 = (-t181 * t237 + t185 * t234) * t524;
t102 = (-t180 * t236 + t184 * t233) * t526;
t101 = (-t179 * t229 + t183 * t228) * t530;
t92 = t650 * pkin(1) + t108;
t91 = t649 * pkin(1) + t107;
t90 = t648 * pkin(1) + t106;
t83 = t647 * pkin(1) + t105;
t80 = t650 * t657 + t398 + t654;
t79 = t649 * t657 + t398 + t653;
t78 = t648 * t657 + t398 + t652;
t77 = t647 * t657 + t398 + t651;
t76 = (t141 * t186 - t144 * t182) * t476;
t75 = (t140 * t185 - t143 * t181) * t477;
t74 = (t139 * t184 - t142 * t180) * t478;
t73 = (t137 * t183 - t138 * t179) * t479;
t68 = t156 * t517 + (-t108 * t472 + t328 * t92) * t460;
t67 = t155 * t519 + (-t107 * t473 + t324 * t91) * t461;
t66 = t154 * t521 + (-t106 * t474 + t320 * t90) * t462;
t65 = (t108 * t484 + t92 * t548) * t393;
t64 = (t107 * t485 + t91 * t549) * t393;
t63 = (t106 * t486 + t90 * t550) * t393;
t62 = (t108 * t487 + t92 * t551) * t393;
t61 = (t107 * t488 + t91 * t552) * t393;
t60 = (t106 * t489 + t90 * t553) * t393;
t59 = t153 * t529 + (-t105 * t475 + t307 * t83) * t463;
t58 = (t105 * t490 + t83 * t555) * t393;
t57 = (t105 * t491 + t83 * t556) * t393;
t56 = (t92 * t484 + t80 * t548) * t393;
t55 = (t91 * t485 + t79 * t549) * t393;
t54 = (t90 * t486 + t78 * t550) * t393;
t53 = (t92 * t487 + t80 * t551) * t393;
t52 = (t91 * t488 + t79 * t552) * t393;
t51 = (t90 * t489 + t78 * t553) * t393;
t41 = (t83 * t490 + t77 * t555) * t393;
t40 = (t83 * t491 + t77 * t556) * t393;
t39 = t389 * t50 * t397;
t38 = t389 * t49 * t396;
t37 = t389 * t48 * t395;
t33 = t136 * t517 + (t328 * t80 - t92 * t472) * t460;
t32 = t135 * t519 + (t324 * t79 - t91 * t473) * t461;
t31 = t134 * t521 + (t320 * t78 - t90 * t474) * t462;
t30 = t389 * t36 * t394;
t29 = t125 * t529 + (t307 * t77 - t83 * t475) * t463;
t28 = t104 * t92 + t108 * t76;
t27 = t103 * t91 + t107 * t75;
t26 = t102 * t90 + t106 * t74;
t25 = t101 * t83 + t105 * t73;
t24 = t104 * t80 + t76 * t92;
t23 = t103 * t79 + t75 * t91;
t22 = t102 * t78 + t74 * t90;
t21 = t101 * t77 + t73 * t83;
t16 = ((-t364 * t89 * t617 - t50 * t613) * t328 * t89 - t50 * t72 * t609 - t330 * t515) * t522;
t15 = ((-t362 * t88 * t618 - t49 * t614) * t324 * t88 - t49 * t71 * t610 - t326 * t515) * t524;
t14 = ((-t360 * t87 * t619 - t48 * t615) * t320 * t87 - t48 * t70 * t611 - t322 * t515) * t526;
t13 = ((-t348 * t82 * t623 - t36 * t616) * t307 * t82 - t36 * t69 * t612 - t309 * t515) * t530;
t12 = (((-pkin(1) * t50 * t509 + t328 * t371) * t364 + t365 * t328 * t502) * t330 * t508 + ((t39 + t44 * t498 * t657 + (-pkin(1) * t328 * t452 + t392 * t89) * t364) * t89 + (t39 + (t50 * t498 - t452) * pkin(1)) * t72) * t516) * t522;
t11 = (((-pkin(1) * t49 * t510 + t324 * t371) * t362 + t363 * t324 * t502) * t326 * t508 + ((t38 + t43 * t499 * t657 + (-pkin(1) * t324 * t453 + t392 * t88) * t362) * t88 + (t38 + (t49 * t499 - t453) * pkin(1)) * t71) * t518) * t524;
t10 = (((-pkin(1) * t48 * t511 + t320 * t371) * t360 + t361 * t320 * t502) * t322 * t508 + ((t37 + t42 * t500 * t657 + (-pkin(1) * t320 * t454 + t392 * t87) * t360) * t87 + (t37 + (t48 * t500 - t454) * pkin(1)) * t70) * t520) * t526;
t9 = (((-pkin(1) * t36 * t513 + t307 * t371) * t348 + t349 * t307 * t502) * t309 * t508 + ((t30 + t34 * t501 * t657 + (-pkin(1) * t307 * t455 + t392 * t82) * t348) * t82 + (t30 + (t36 * t501 - t455) * pkin(1)) * t69) * t528) * t530;
t8 = -t108 * t12 + t129 * t283 + t274 * t133 - t92 * t16 + (t156 * t330 - t540) * t456 + t408 * t626 + (-t451 + (-0.2e1 * t532 - t557) * t494) * t328 + ((t646 * t194 - t577) * t283 + t274 * (-t646 * t198 - t581) + ((-t282 / 0.2e1 + t281 / 0.2e1) * rSges(3,2) + (t273 / 0.2e1 + t272 / 0.2e1) * rSges(3,1)) * t626) * m(3);
t7 = -t107 * t11 + t128 * t280 + t271 * t132 - t91 * t15 + (t155 * t326 - t541) * t457 + t409 * t627 + (-t451 + (-0.2e1 * t533 - t558) * t495) * t324 + ((t645 * t193 - t578) * t280 + t271 * (-t645 * t197 - t582) + ((-t279 / 0.2e1 + t278 / 0.2e1) * rSges(3,2) + (t270 / 0.2e1 + t269 / 0.2e1) * rSges(3,1)) * t627) * m(3);
t6 = -t106 * t10 + t127 * t277 + t268 * t131 - t90 * t14 + (t154 * t322 - t542) * t458 + t410 * t628 + (-t451 + (-0.2e1 * t534 - t559) * t496) * t320 + ((t644 * t192 - t579) * t277 + t268 * (-t644 * t196 - t583) + ((-t276 / 0.2e1 + t275 / 0.2e1) * rSges(3,2) + (t267 / 0.2e1 + t266 / 0.2e1) * rSges(3,1)) * t628) * m(3);
t5 = -t105 * t9 + t126 * t260 - t83 * t13 + t257 * t130 + (t153 * t309 - t543) * t459 + t411 * t629 + (-t451 + (-0.2e1 * t538 - t560) * t497) * t307 + ((t643 * t191 - t580) * t260 + t257 * (-t643 * t195 - t584) + ((-t259 / 0.2e1 + t258 / 0.2e1) * rSges(3,2) + (t256 / 0.2e1 + t255 / 0.2e1) * rSges(3,1)) * t629) * m(3);
t4 = -t80 * t16 - t92 * t12 + (-t329 * t451 + t437 * t656) * t364 + sin(qJ(1,1)) * (t194 * t633 + t450 * t198) - t408 * t657 * t601 + (t437 * t655 + (t136 * t330 - t540) * t514) * t358 + (t364 * t638 + (0.2e1 * (-t223 * t358 - t365 * t504) * t364 + t365 * t358 * t493 + t637) * t328) * t494 + (-m(3) * t577 + t129) * t283 + t274 * (-m(3) * t581 + t133) + (-t450 * t194 + t198 * t633) * cos(qJ(1,1)) + ((t174 + t597) * cos(t297) + (t170 - t593) * sin(t297)) * t636 + ((t174 - t597) * cos(t296) + (t170 + t593) * sin(t296)) * t635;
t3 = -t79 * t15 - t91 * t11 + (-t325 * t451 + t438 * t656) * t362 + sin(qJ(1,2)) * (t193 * t633 + t450 * t197) - t409 * t657 * t602 + (t438 * t655 + (t135 * t326 - t541) * t514) * t356 + (t362 * t638 + (0.2e1 * (-t223 * t356 - t363 * t504) * t362 + t363 * t356 * t493 + t637) * t324) * t495 + (-m(3) * t578 + t128) * t280 + t271 * (-m(3) * t582 + t132) + (-t450 * t193 + t197 * t633) * cos(qJ(1,2)) + ((t173 + t598) * cos(t295) + (t169 - t594) * sin(t295)) * t636 + ((t173 - t598) * cos(t294) + (t169 + t594) * sin(t294)) * t635;
t2 = -t78 * t14 - t90 * t10 + (-t321 * t451 + t439 * t656) * t360 + sin(qJ(1,3)) * (t192 * t633 + t450 * t196) - t410 * t657 * t603 + (t439 * t655 + (t134 * t322 - t542) * t514) * t354 + (t360 * t638 + (0.2e1 * (-t223 * t354 - t361 * t504) * t360 + t361 * t354 * t493 + t637) * t320) * t496 + (-m(3) * t579 + t127) * t277 + t268 * (-m(3) * t583 + t131) + (-t450 * t192 + t196 * t633) * cos(qJ(1,3)) + ((t172 + t599) * cos(t293) + (t168 - t595) * sin(t293)) * t636 + ((t172 - t599) * cos(t292) + (t168 + t595) * sin(t292)) * t635;
t1 = -t77 * t13 - t83 * t9 + (-t308 * t451 + t440 * t656) * t348 + sin(qJ(1,4)) * (t191 * t633 + t450 * t195) + (t440 * t655 + (t125 * t309 - t543) * t514) * t346 - t411 * t657 * t604 + (t348 * t638 + (0.2e1 * (-t223 * t346 - t349 * t504) * t348 + t349 * t346 * t493 + t637) * t307) * t497 + (-m(3) * t580 + t126) * t260 + t257 * (-m(3) * t584 + t130) + (-t450 * t191 + t195 * t633) * cos(qJ(1,4)) + ((t171 + t600) * cos(t262) + (t167 - t596) * sin(t262)) * t636 + ((t171 - t600) * cos(t261) + (t167 + t596) * sin(t261)) * t635;
t17 = [(t65 * t484 + t56 * t548) * t569 + (t65 * t487 + t56 * t551) * t573 + t4 * t464 + t8 * t424 + (t64 * t485 + t55 * t549) * t570 + (t64 * t488 + t55 * t552) * t574 + t3 * t465 + t7 * t425 + (t63 * t486 + t54 * t550) * t571 + (t63 * t489 + t54 * t553) * t575 + t2 * t466 + t6 * t426 + (t41 * t555 + t58 * t490) * t572 + (t41 * t556 + t58 * t491) * t576 + t1 * t470 + t5 * t430 + t605 - t162 * t607 - t163 * t608 - m(4) * g(1) + (-t65 * t412 + (t56 * t523 + (t136 * t548 + t144 * t480) * t390) * t328 - t64 * t413 + (t55 * t525 + (t135 * t549 + t143 * t481) * t390) * t324 - t63 * t414 + (t54 * t527 + (t134 * t550 + t142 * t482) * t390) * t320 - t58 * t415 + (t41 * t531 + (t125 * t555 + t138 * t483) * t390) * t307) * t512; (t484 * t62 + t53 * t548) * t569 + (t487 * t62 + t53 * t551) * t573 + t4 * t467 + t8 * t427 + (t485 * t61 + t52 * t549) * t570 + (t488 * t61 + t52 * t552) * t574 + t3 * t468 + t7 * t428 + (t486 * t60 + t51 * t550) * t571 + (t489 * t60 + t51 * t553) * t575 + t2 * t469 + t6 * t429 + (t40 * t555 + t490 * t57) * t572 + (t40 * t556 + t491 * t57) * t576 + t1 * t471 + t5 * t431 + t606 + t163 * t607 - t162 * t608 - m(4) * g(2) + (-t62 * t412 + (t53 * t523 + (t136 * t551 + t141 * t480) * t390) * t328 - t61 * t413 + (t52 * t525 + (t135 * t552 + t140 * t481) * t390) * t324 - t60 * t414 + (t51 * t527 + (t134 * t553 + t139 * t482) * t390) * t320 - t57 * t415 + (t40 * t531 + (t125 * t556 + t137 * t483) * t390) * t307) * t512; t4 * t420 + t3 * t421 + t2 * t422 + t1 * t423 - t8 * t404 - t7 * t405 - t6 * t406 - t5 * t407 + (t33 * t548 + t484 * t68) * t569 + (t32 * t549 + t485 * t67) * t570 + (t31 * t550 + t486 * t66) * t571 + (t29 * t555 + t490 * t59) * t572 + (t33 * t551 + t487 * t68) * t573 + (t32 * t552 + t488 * t67) * t574 + (t31 * t553 + t489 * t66) * t575 + (-t136 * t16 - t156 * t12 + t232 * t330 * t456 + (g(3) * t646 + (t194 * t274 + t198 * t283) * (rSges(3,1) * t358 + rSges(3,2) * t364)) * m(3) + (t557 / 0.2e1 + t532) * t50 ^ 2 + (t273 * t436 + (rSges(3,1) * t272 + (t281 + t282) * rSges(3,2)) * t503) * t86) * t517 + (-t135 * t15 - t155 * t11 + t232 * t326 * t457 + (g(3) * t645 + (t193 * t271 + t197 * t280) * (rSges(3,1) * t356 + rSges(3,2) * t362)) * m(3) + (t558 / 0.2e1 + t533) * t49 ^ 2 + (t270 * t436 + (rSges(3,1) * t269 + (t278 + t279) * rSges(3,2)) * t503) * t85) * t519 + (-t134 * t14 - t154 * t10 + t232 * t322 * t458 + (g(3) * t644 + (t192 * t268 + t196 * t277) * (rSges(3,1) * t354 + rSges(3,2) * t360)) * m(3) + (t559 / 0.2e1 + t534) * t48 ^ 2 + (t267 * t436 + (rSges(3,1) * t266 + (t275 + t276) * rSges(3,2)) * t503) * t84) * t521 + (-t125 * t13 - t153 * t9 + t232 * t309 * t459 + (g(3) * t643 + (t191 * t257 + t195 * t260) * (rSges(3,1) * t346 + rSges(3,2) * t348)) * m(3) + (t560 / 0.2e1 + t538) * t36 ^ 2 + (t256 * t436 + (rSges(3,1) * t255 + (t258 + t259) * rSges(3,2)) * t503) * t81) * t529 + (t29 * t556 + t491 * t59) * t576 - m(4) * g(3) + (m(4) + (t33 * t328 - t68 * t472 + (t136 * t328 - t156 * t472) * t517) * t460 + (t32 * t324 - t67 * t473 + (t135 * t324 - t155 * t473) * t519) * t461 + (t31 * t320 - t66 * t474 + (t134 * t320 - t154 * t474) * t521) * t462 + (t29 * t307 - t59 * t475 + (t125 * t307 - t153 * t475) * t529) * t463 + (t307 ^ 2 + t320 ^ 2 + t324 ^ 2 + t328 ^ 2) * t232 * t390 ^ 2) * t351; (t24 * t548 + t28 * t484) * t569 + (t24 * t551 + t28 * t487) * t573 + t104 * t4 + t76 * t8 + (t23 * t549 + t27 * t485) * t570 + (t23 * t552 + t27 * t488) * t574 + t103 * t3 + t75 * t7 + (t22 * t550 + t26 * t486) * t571 + (t22 * t553 + t26 * t489) * t575 + t102 * t2 + t74 * t6 + (t21 * t555 + t25 * t490) * t572 + (t21 * t556 + t25 * t491) * t576 + t101 * t1 + t73 * t5 - t162 * t605 + t163 * t606 + (Icges(4,3) + m(4) * (t375 ^ 2 + t376 ^ 2)) * t350 + m(4) * ((g(1) * t376 + g(2) * t375) * t298 + (g(1) * t375 - g(2) * t376) * t299) + ((t104 * t136 + t156 * t76) * t517 + (t24 * t328 - t28 * t472) * t460 + (t103 * t135 + t155 * t75) * t519 + (t23 * t324 - t27 * t473) * t461 + (t102 * t134 + t154 * t74) * t521 + (t22 * t320 - t26 * t474) * t462 + (t101 * t125 + t153 * t73) * t529 + (t21 * t307 - t25 * t475) * t463) * t351;];
tauX  = t17;
