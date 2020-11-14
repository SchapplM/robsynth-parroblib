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
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
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

function tauX = P4PRRRR8V2G3A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_slag_vp2: xDP has to be [4x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [4 1]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_slag_vp2: xDDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G3A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:25:34
% EndTime: 2020-08-07 11:25:51
% DurationCPUTime: 17.18s
% Computational Cost: add. (79027->760), mult. (170417->1340), div. (9248->9), fcn. (162006->30), ass. (0->473)
t310 = sin(qJ(2,1));
t316 = cos(qJ(2,1));
t323 = pkin(7) + pkin(6);
t231 = pkin(2) * t310 - t316 * t323;
t289 = sin(pkin(4));
t291 = cos(pkin(4));
t309 = sin(qJ(3,1));
t432 = t291 * t309;
t172 = pkin(3) * t432 + t231 * t289;
t315 = cos(qJ(3,1));
t447 = t289 * t310;
t286 = t315 ^ 2;
t530 = pkin(3) * t286;
t136 = 0.1e1 / (pkin(2) * t432 + t172 * t315 + t447 * t530);
t308 = sin(qJ(2,2));
t314 = cos(qJ(2,2));
t230 = pkin(2) * t308 - t314 * t323;
t307 = sin(qJ(3,2));
t434 = t291 * t307;
t171 = pkin(3) * t434 + t230 * t289;
t313 = cos(qJ(3,2));
t449 = t289 * t308;
t285 = t313 ^ 2;
t531 = pkin(3) * t285;
t135 = 0.1e1 / (pkin(2) * t434 + t171 * t313 + t449 * t531);
t306 = sin(qJ(2,3));
t312 = cos(qJ(2,3));
t229 = pkin(2) * t306 - t312 * t323;
t305 = sin(qJ(3,3));
t436 = t291 * t305;
t170 = pkin(3) * t436 + t229 * t289;
t311 = cos(qJ(3,3));
t451 = t289 * t306;
t284 = t311 ^ 2;
t532 = pkin(3) * t284;
t134 = 0.1e1 / (pkin(2) * t436 + t170 * t311 + t451 * t532);
t293 = sin(qJ(2,4));
t295 = cos(qJ(2,4));
t227 = pkin(2) * t293 - t295 * t323;
t292 = sin(qJ(3,4));
t439 = t291 * t292;
t169 = pkin(3) * t439 + t227 * t289;
t294 = cos(qJ(3,4));
t455 = t289 * t293;
t281 = t294 ^ 2;
t533 = pkin(3) * t281;
t133 = 0.1e1 / (pkin(2) * t439 + t169 * t294 + t455 * t533);
t300 = legFrame(1,2);
t270 = sin(t300);
t274 = cos(t300);
t222 = g(1) * t270 + g(2) * t274;
t226 = g(1) * t274 - g(2) * t270;
t290 = cos(pkin(8));
t260 = g(3) * t290;
t288 = sin(pkin(8));
t378 = -t226 * t288 - t260;
t564 = t222 * t289 + t378 * t291;
t299 = legFrame(2,2);
t269 = sin(t299);
t273 = cos(t299);
t221 = g(1) * t269 + g(2) * t273;
t225 = g(1) * t273 - g(2) * t269;
t380 = -t225 * t288 - t260;
t563 = t221 * t289 + t380 * t291;
t298 = legFrame(3,2);
t268 = sin(t298);
t272 = cos(t298);
t220 = g(1) * t268 + g(2) * t272;
t224 = g(1) * t272 - g(2) * t268;
t382 = -t224 * t288 - t260;
t562 = t220 * t289 + t382 * t291;
t297 = legFrame(4,2);
t267 = sin(t297);
t271 = cos(t297);
t219 = g(1) * t267 + g(2) * t271;
t223 = g(1) * t271 - g(2) * t267;
t384 = -t223 * t288 - t260;
t561 = t219 * t289 + t384 * t291;
t228 = pkin(2) * t295 + t293 * t323;
t437 = t291 * t295;
t440 = t290 * t291;
t529 = pkin(3) * t294;
t125 = (t288 * t293 - t290 * t437) * t529 - t228 * t440 + t227 * t288;
t457 = t288 * t291;
t126 = (t288 * t437 + t290 * t293) * t529 + t228 * t457 + t290 * t227;
t320 = xDP(3);
t325 = xP(4);
t279 = sin(t325);
t280 = cos(t325);
t328 = koppelP(4,2);
t332 = koppelP(4,1);
t207 = t279 * t332 + t280 * t328;
t211 = -t279 * t328 + t280 * t332;
t319 = xDP(4);
t321 = xDP(2);
t322 = xDP(1);
t372 = (t211 * t319 + t321) * t267 - (-t207 * t319 + t322) * t271;
t338 = 0.1e1 / pkin(3);
t496 = t133 * t338;
t86 = (-t372 * t125 + t126 * t320) * t496;
t560 = -0.2e1 * t86;
t232 = pkin(2) * t312 + t306 * t323;
t430 = t291 * t312;
t528 = pkin(3) * t311;
t127 = (t288 * t306 - t290 * t430) * t528 - t232 * t440 + t229 * t288;
t130 = (t288 * t430 + t290 * t306) * t528 + t232 * t457 + t290 * t229;
t329 = koppelP(3,2);
t333 = koppelP(3,1);
t208 = t279 * t333 + t280 * t329;
t212 = -t279 * t329 + t280 * t333;
t371 = (t212 * t319 + t321) * t268 - (-t208 * t319 + t322) * t272;
t495 = t134 * t338;
t91 = (-t127 * t371 + t130 * t320) * t495;
t559 = -0.2e1 * t91;
t233 = pkin(2) * t314 + t308 * t323;
t429 = t291 * t314;
t527 = pkin(3) * t313;
t128 = (t288 * t308 - t290 * t429) * t527 - t233 * t440 + t230 * t288;
t131 = (t288 * t429 + t290 * t308) * t527 + t233 * t457 + t290 * t230;
t330 = koppelP(2,2);
t334 = koppelP(2,1);
t209 = t279 * t334 + t280 * t330;
t213 = -t279 * t330 + t280 * t334;
t370 = (t213 * t319 + t321) * t269 - (-t209 * t319 + t322) * t273;
t494 = t135 * t338;
t92 = (-t128 * t370 + t131 * t320) * t494;
t558 = -0.2e1 * t92;
t234 = pkin(2) * t316 + t310 * t323;
t428 = t291 * t316;
t526 = pkin(3) * t315;
t129 = (t288 * t310 - t290 * t428) * t526 - t234 * t440 + t231 * t288;
t132 = (t288 * t428 + t290 * t310) * t526 + t234 * t457 + t290 * t231;
t331 = koppelP(1,2);
t335 = koppelP(1,1);
t210 = t279 * t335 + t280 * t331;
t214 = -t279 * t331 + t280 * t335;
t369 = (t214 * t319 + t321) * t270 - (-t210 * t319 + t322) * t274;
t493 = t136 * t338;
t93 = (-t129 * t369 + t132 * t320) * t493;
t557 = -0.2e1 * t93;
t392 = t294 * mrSges(3,1) - mrSges(3,2) * t292;
t391 = t311 * mrSges(3,1) - mrSges(3,2) * t305;
t390 = t313 * mrSges(3,1) - mrSges(3,2) * t307;
t389 = t315 * mrSges(3,1) - mrSges(3,2) * t309;
t296 = Ifges(3,1) - Ifges(3,2);
t317 = pkin(2) * mrSges(3,2);
t544 = -2 * Ifges(3,4);
t548 = Ifges(3,4) + t286 * t544 + (-t296 * t309 + t317) * t315;
t547 = Ifges(3,4) + t285 * t544 + (-t296 * t307 + t317) * t313;
t546 = Ifges(3,4) + t284 * t544 + (-t296 * t305 + t317) * t311;
t545 = Ifges(3,4) + t281 * t544 + (-t296 * t292 + t317) * t294;
t318 = mrSges(3,1) * pkin(2);
t543 = pkin(3) * t86;
t542 = pkin(3) * t91;
t541 = pkin(3) * t92;
t540 = pkin(3) * t93;
t261 = mrSges(3,2) * pkin(6) - Ifges(3,6);
t539 = -t261 / 0.2e1;
t262 = mrSges(3,1) * pkin(6) - Ifges(3,5);
t538 = t262 / 0.2e1;
t537 = pkin(2) * t292;
t536 = pkin(2) * t305;
t535 = pkin(2) * t307;
t534 = pkin(2) * t309;
t525 = g(3) * t288;
t339 = pkin(2) ^ 2;
t257 = t323 ^ 2 + t339;
t337 = pkin(3) ^ 2;
t408 = t292 * t543;
t418 = 0.2e1 * pkin(2) * pkin(3);
t438 = t291 * t293;
t189 = t288 * t438 - t290 * t295;
t454 = t289 * t294;
t154 = t189 * t292 + t288 * t454;
t190 = t288 * t295 + t290 * t438;
t155 = t190 * t292 + t290 * t454;
t98 = (t154 * t320 + t372 * t155) * t133;
t524 = (-t323 * t408 + (t281 * t337 + t294 * t418 + t257) * t98) * t98;
t523 = mrSges(3,1) * t292;
t522 = mrSges(3,1) * t305;
t521 = mrSges(3,1) * t307;
t520 = mrSges(3,1) * t309;
t435 = t291 * t306;
t192 = t288 * t435 - t290 * t312;
t446 = t289 * t311;
t156 = t192 * t305 + t288 * t446;
t195 = t288 * t312 + t290 * t435;
t159 = t195 * t305 + t290 * t446;
t102 = (t156 * t320 + t159 * t371) * t134;
t407 = t305 * t542;
t515 = t102 * (-t323 * t407 + (t284 * t337 + t311 * t418 + t257) * t102);
t433 = t291 * t308;
t193 = t288 * t433 - t290 * t314;
t444 = t289 * t313;
t157 = t193 * t307 + t288 * t444;
t196 = t288 * t314 + t290 * t433;
t160 = t196 * t307 + t290 * t444;
t103 = (t157 * t320 + t160 * t370) * t135;
t406 = t307 * t541;
t514 = t103 * (-t323 * t406 + (t285 * t337 + t313 * t418 + t257) * t103);
t431 = t291 * t310;
t194 = t288 * t431 - t290 * t316;
t442 = t289 * t315;
t158 = t194 * t309 + t288 * t442;
t197 = t288 * t316 + t290 * t431;
t161 = t197 * t309 + t290 * t442;
t104 = (t158 * t320 + t161 * t369) * t136;
t405 = t309 * t540;
t513 = t104 * (-t323 * t405 + (t286 * t337 + t315 * t418 + t257) * t104);
t512 = t295 * t98;
t511 = t323 * t98;
t275 = m(3) * pkin(2) + mrSges(2,1);
t510 = t102 * t312;
t509 = t102 * t323;
t508 = t103 * t314;
t507 = t103 * t323;
t506 = t104 * t316;
t505 = t104 * t323;
t504 = t125 * t338;
t503 = t126 * t338;
t502 = t127 * t338;
t501 = t128 * t338;
t500 = t129 * t338;
t499 = t130 * t338;
t498 = t131 * t338;
t497 = t132 * t338;
t388 = mrSges(3,2) * t294 + t523;
t153 = t392 * t291 - t388 * t455;
t492 = t153 * t338;
t491 = t155 * t267;
t490 = t155 * t271;
t489 = t159 * t268;
t488 = t159 * t272;
t487 = t160 * t269;
t486 = t160 * t273;
t485 = t161 * t270;
t484 = t161 * t274;
t387 = mrSges(3,2) * t311 + t522;
t162 = t391 * t291 - t387 * t451;
t483 = t162 * t338;
t386 = mrSges(3,2) * t313 + t521;
t163 = t390 * t291 - t386 * t449;
t482 = t163 * t338;
t385 = mrSges(3,2) * t315 + t520;
t164 = t389 * t291 - t385 * t447;
t481 = t164 * t338;
t215 = t392 + t275;
t258 = m(3) * pkin(6) - mrSges(2,2) + mrSges(3,3);
t165 = t215 * t295 + t258 * t293;
t480 = t165 * t289;
t216 = t391 + t275;
t166 = t216 * t312 + t258 * t306;
t479 = t166 * t289;
t217 = t390 + t275;
t167 = t217 * t314 + t258 * t308;
t478 = t167 * t289;
t218 = t389 + t275;
t168 = t218 * t316 + t258 * t310;
t477 = t168 * t289;
t191 = -t261 * t294 - t292 * t262;
t476 = t191 * t338;
t199 = -t261 * t311 - t305 * t262;
t475 = t199 * t338;
t200 = -t261 * t313 - t307 * t262;
t474 = t200 * t338;
t201 = -t261 * t315 - t309 * t262;
t473 = t201 * t338;
t471 = t219 * t291;
t469 = t220 * t291;
t467 = t221 * t291;
t465 = t222 * t291;
t460 = mrSges(3,2) * t289 * t260;
t459 = t258 * t290;
t458 = t288 * t289;
t456 = t289 * t292;
t453 = t289 * t295;
t452 = t289 * t305;
t450 = t289 * t307;
t448 = t289 * t309;
t445 = t289 * t312;
t443 = t289 * t314;
t441 = t289 * t316;
t427 = t291 * t338;
t426 = t293 * t294;
t421 = t306 * t311;
t420 = t308 * t313;
t419 = t310 * t315;
t413 = -0.2e1 * t317;
t404 = t292 * t511;
t403 = t305 * t509;
t402 = t307 * t507;
t401 = t309 * t505;
t400 = t267 * t504;
t399 = t271 * t504;
t398 = t268 * t502;
t397 = t272 * t502;
t396 = t269 * t501;
t395 = t273 * t501;
t394 = t270 * t500;
t393 = t274 * t500;
t383 = t223 * t290 - t525;
t381 = t224 * t290 - t525;
t379 = t225 * t290 - t525;
t377 = t226 * t290 - t525;
t368 = t207 * t271 + t211 * t267;
t367 = t208 * t272 + t212 * t268;
t366 = t209 * t273 + t213 * t269;
t365 = t210 * t274 + t214 * t270;
t41 = t404 - t543;
t13 = (((t291 * t86 + t98 * t453) * t533 + ((-t408 + t511) * t293 + pkin(2) * t512) * t454 + t41 * t291) * t98 + (t86 * t453 + (t281 * t291 - t426 * t456 - t291) * t98) * t543) * t133;
t344 = Ifges(3,1) + Ifges(2,3) + (pkin(6) ^ 2 + t339) * m(3) + 0.2e1 * pkin(6) * mrSges(3,3);
t149 = -t296 * t281 + 0.2e1 * (Ifges(3,4) * t292 + t318) * t294 + t292 * t413 + t344;
t185 = pkin(3) * t426 + t227;
t17 = t133 * t427 * t524 + (-t291 * t404 + (-t185 * t456 + (pkin(2) * t294 + t533) * t291) * t86) / (t185 * t454 + (pkin(2) + t529) * t439) * t86;
t21 = (-t294 * t524 - (pkin(2) * t86 - t294 * t41) * t543) * t133;
t239 = t258 * t525;
t5 = -t21 * t480 - t149 * t13 - t191 * t17 + ((t539 * t292 + t538 * t294) * t86 + (t318 * t292 + t545) * t98) * t560 + (-t215 * t561 - t223 * t459 + t239) * t295 + t293 * (t383 * t215 - t258 * t561);
t343 = t293 * t561 + t383 * t295;
t97 = t98 ^ 2;
t9 = -t153 * t21 - t191 * t13 - Ifges(3,3) * t17 + (pkin(2) * t523 + t545) * t97 + ((t384 * t289 - t471) * mrSges(3,1) + t343 * mrSges(3,2)) * t294 + t292 * (t460 + (t223 * t458 + t471) * mrSges(3,2) + t343 * mrSges(3,1));
t364 = t155 * t5 - t9 * t504;
t43 = t403 - t542;
t14 = (((t102 * t445 + t291 * t91) * t532 + ((-t407 + t509) * t306 + pkin(2) * t510) * t446 + t43 * t291) * t102 + (t91 * t445 + (t284 * t291 - t421 * t452 - t291) * t102) * t542) * t134;
t186 = pkin(3) * t421 + t229;
t18 = t134 * t427 * t515 + (-t291 * t403 + (-t186 * t452 + t291 * (pkin(2) * t311 + t532)) * t91) / (t186 * t446 + (pkin(2) + t528) * t436) * t91;
t22 = (-t311 * t515 - (pkin(2) * t91 - t311 * t43) * t542) * t134;
t342 = t306 * t562 + t381 * t312;
t99 = t102 ^ 2;
t10 = -t162 * t22 - t199 * t14 - Ifges(3,3) * t18 + (pkin(2) * t522 + t546) * t99 + ((t382 * t289 - t469) * mrSges(3,1) + t342 * mrSges(3,2)) * t311 + t305 * (t460 + (t224 * t458 + t469) * mrSges(3,2) + t342 * mrSges(3,1));
t150 = -t296 * t284 + 0.2e1 * (Ifges(3,4) * t305 + t318) * t311 + t305 * t413 + t344;
t6 = -t22 * t479 - t150 * t14 - t199 * t18 + ((t539 * t305 + t538 * t311) * t91 + (t318 * t305 + t546) * t102) * t559 + (-t216 * t562 - t224 * t459 + t239) * t312 + t306 * (t381 * t216 - t258 * t562);
t363 = -t10 * t502 + t159 * t6;
t100 = t103 ^ 2;
t44 = t402 - t541;
t15 = (((t103 * t443 + t291 * t92) * t531 + ((-t406 + t507) * t308 + pkin(2) * t508) * t444 + t44 * t291) * t103 + (t92 * t443 + (t285 * t291 - t420 * t450 - t291) * t103) * t541) * t135;
t187 = pkin(3) * t420 + t230;
t19 = t135 * t427 * t514 + (-t291 * t402 + (-t187 * t450 + t291 * (pkin(2) * t313 + t531)) * t92) / (t187 * t444 + (pkin(2) + t527) * t434) * t92;
t23 = (-t313 * t514 - (pkin(2) * t92 - t313 * t44) * t541) * t135;
t341 = t308 * t563 + t379 * t314;
t11 = -t163 * t23 - t200 * t15 - Ifges(3,3) * t19 + (pkin(2) * t521 + t547) * t100 + ((t380 * t289 - t467) * mrSges(3,1) + t341 * mrSges(3,2)) * t313 + t307 * (t460 + (t225 * t458 + t467) * mrSges(3,2) + t341 * mrSges(3,1));
t151 = -t296 * t285 + 0.2e1 * (Ifges(3,4) * t307 + t318) * t313 + t307 * t413 + t344;
t7 = -t23 * t478 - t151 * t15 - t200 * t19 + ((t539 * t307 + t538 * t313) * t92 + (t318 * t307 + t547) * t103) * t558 + (-t217 * t563 - t225 * t459 + t239) * t314 + t308 * (t379 * t217 - t258 * t563);
t362 = -t11 * t501 + t160 * t7;
t101 = t104 ^ 2;
t45 = t401 - t540;
t16 = (((t104 * t441 + t291 * t93) * t530 + ((-t405 + t505) * t310 + pkin(2) * t506) * t442 + t45 * t291) * t104 + (t93 * t441 + (t286 * t291 - t419 * t448 - t291) * t104) * t540) * t136;
t188 = pkin(3) * t419 + t231;
t20 = t136 * t427 * t513 + (-t291 * t401 + (-t188 * t448 + t291 * (pkin(2) * t315 + t530)) * t93) / (t188 * t442 + (pkin(2) + t526) * t432) * t93;
t24 = (-t315 * t513 - (pkin(2) * t93 - t315 * t45) * t540) * t136;
t340 = t310 * t564 + t377 * t316;
t12 = -t164 * t24 - t201 * t16 - Ifges(3,3) * t20 + (pkin(2) * t520 + t548) * t101 + ((t378 * t289 - t465) * mrSges(3,1) + t340 * mrSges(3,2)) * t315 + t309 * (t460 + (t226 * t458 + t465) * mrSges(3,2) + t340 * mrSges(3,1));
t152 = -t296 * t286 + 0.2e1 * (Ifges(3,4) * t309 + t318) * t315 + t309 * t413 + t344;
t8 = -t24 * t477 - t152 * t16 - t201 * t20 + ((t539 * t309 + t538 * t315) * t93 + (t318 * t309 + t548) * t104) * t557 + (-t218 * t564 - t226 * t459 + t239) * t316 + t310 * (t377 * t218 - t258 * t564);
t361 = -t12 * t500 + t161 * t8;
t360 = pkin(3) * t456 - t227 * t291;
t359 = pkin(3) * t452 - t229 * t291;
t358 = pkin(3) * t450 - t230 * t291;
t357 = pkin(3) * t448 - t231 * t291;
t356 = Ifges(3,3) * t504 - t155 * t191;
t355 = Ifges(3,3) * t502 - t159 * t199;
t354 = Ifges(3,3) * t501 - t160 * t200;
t353 = Ifges(3,3) * t500 - t161 * t201;
t352 = t125 * t476 - t149 * t155;
t351 = t127 * t475 - t150 * t159;
t350 = t128 * t474 - t151 * t160;
t349 = t129 * t473 - t152 * t161;
t348 = t125 * t492 - t155 * t480;
t347 = t127 * t483 - t159 * t479;
t346 = t128 * t482 - t160 * t478;
t345 = t129 * t481 - t161 * t477;
t327 = mrSges(4,1);
t326 = mrSges(4,2);
t304 = xDDP(1);
t303 = xDDP(2);
t302 = xDDP(3);
t301 = xDDP(4);
t287 = t319 ^ 2;
t282 = m(1) + m(2) + m(3);
t206 = -t326 * t279 + t280 * t327;
t205 = t327 * t279 + t280 * t326;
t148 = -t210 * t301 - t214 * t287 + t304;
t147 = -t209 * t301 - t213 * t287 + t304;
t146 = -t208 * t301 - t212 * t287 + t304;
t145 = -t207 * t301 - t211 * t287 + t304;
t144 = -t210 * t287 + t214 * t301 + t303;
t143 = -t209 * t287 + t213 * t301 + t303;
t142 = -t208 * t287 + t212 * t301 + t303;
t141 = -t207 * t287 + t211 * t301 + t303;
t140 = t234 * t290 + t288 * t357;
t139 = t233 * t290 + t288 * t358;
t138 = t232 * t290 + t288 * t359;
t137 = t228 * t290 + t288 * t360;
t124 = -t197 * t530 - t234 * t288 * t315 + (pkin(2) * t448 + t315 * t357) * t290;
t123 = -t196 * t531 - t233 * t288 * t313 + (pkin(2) * t450 + t313 * t358) * t290;
t122 = -t195 * t532 - t232 * t288 * t311 + (pkin(2) * t452 + t311 * t359) * t290;
t121 = -t190 * t533 - t228 * t288 * t294 + (pkin(2) * t456 + t294 * t360) * t290;
t120 = -(t194 * t274 - t270 * t447) * t530 + (t140 * t274 + t172 * t270) * t315 + (t270 * t291 + t274 * t458) * t534;
t119 = (t194 * t270 + t274 * t447) * t530 + (-t140 * t270 + t172 * t274) * t315 + (-t270 * t458 + t274 * t291) * t534;
t118 = -(t193 * t273 - t269 * t449) * t531 + (t139 * t273 + t171 * t269) * t313 + (t269 * t291 + t273 * t458) * t535;
t117 = (t193 * t269 + t273 * t449) * t531 + (-t139 * t269 + t171 * t273) * t313 + (-t269 * t458 + t273 * t291) * t535;
t116 = -(t192 * t272 - t268 * t451) * t532 + (t138 * t272 + t170 * t268) * t311 + (t268 * t291 + t272 * t458) * t536;
t115 = (t192 * t268 + t272 * t451) * t532 + (-t138 * t268 + t170 * t272) * t311 + (-t268 * t458 + t272 * t291) * t536;
t114 = -(t189 * t271 - t267 * t455) * t533 + (t137 * t271 + t169 * t267) * t294 + (t267 * t291 + t271 * t458) * t537;
t113 = (t189 * t267 + t271 * t455) * t533 + (-t137 * t267 + t169 * t271) * t294 + (-t267 * t458 + t271 * t291) * t537;
t112 = t365 * t161 * t136;
t111 = t366 * t160 * t135;
t110 = t367 * t159 * t134;
t109 = t368 * t155 * t133;
t108 = t365 * t129 * t493;
t107 = t366 * t128 * t494;
t106 = t367 * t127 * t495;
t105 = t368 * t125 * t496;
t96 = (Ifges(3,3) * t497 + t124 * t164 + t158 * t201) * t136;
t95 = (Ifges(3,3) * t498 + t123 * t163 + t157 * t200) * t135;
t94 = (Ifges(3,3) * t499 + t122 * t162 + t156 * t199) * t134;
t90 = t93 ^ 2;
t89 = t92 ^ 2;
t88 = t91 ^ 2;
t87 = (Ifges(3,3) * t503 + t121 * t153 + t154 * t191) * t133;
t85 = t86 ^ 2;
t84 = (t124 * t282 + t132 * t481 + t158 * t477) * t136;
t83 = (t123 * t282 + t131 * t482 + t157 * t478) * t135;
t82 = (t122 * t282 + t130 * t483 + t156 * t479) * t134;
t81 = (t121 * t282 + t126 * t492 + t154 * t480) * t133;
t80 = (t124 * t477 + t132 * t473 + t152 * t158) * t136;
t79 = (t123 * t478 + t131 * t474 + t151 * t157) * t135;
t78 = (t122 * t479 + t130 * t475 + t150 * t156) * t134;
t77 = (t121 * t480 + t126 * t476 + t149 * t154) * t133;
t76 = (t119 * t214 - t120 * t210) * t136;
t75 = (t117 * t213 - t118 * t209) * t135;
t74 = (t115 * t212 - t116 * t208) * t134;
t73 = (t113 * t211 - t114 * t207) * t133;
t72 = (t120 * t164 + t274 * t353) * t136;
t71 = (t119 * t164 - t270 * t353) * t136;
t70 = (t118 * t163 + t273 * t354) * t135;
t69 = (t117 * t163 - t269 * t354) * t135;
t68 = (t116 * t162 + t272 * t355) * t134;
t67 = (t115 * t162 - t268 * t355) * t134;
t66 = (t114 * t153 + t271 * t356) * t133;
t65 = (t113 * t153 - t267 * t356) * t133;
t64 = (t120 * t282 + t274 * t345) * t136;
t63 = (t119 * t282 - t270 * t345) * t136;
t62 = (t118 * t282 + t273 * t346) * t135;
t61 = (t117 * t282 - t269 * t346) * t135;
t60 = (t116 * t282 + t272 * t347) * t134;
t59 = (t115 * t282 - t268 * t347) * t134;
t58 = (t114 * t282 + t271 * t348) * t133;
t57 = (t113 * t282 - t267 * t348) * t133;
t56 = (t120 * t477 + t274 * t349) * t136;
t55 = (t119 * t477 - t270 * t349) * t136;
t54 = (t118 * t478 + t273 * t350) * t135;
t53 = (t117 * t478 - t269 * t350) * t135;
t52 = (t116 * t479 + t272 * t351) * t134;
t51 = (t115 * t479 - t268 * t351) * t134;
t50 = (t114 * t480 + t271 * t352) * t133;
t49 = (t113 * t480 - t267 * t352) * t133;
t40 = -Ifges(3,3) * t108 + t112 * t201 + t164 * t76;
t39 = -Ifges(3,3) * t107 + t111 * t200 + t163 * t75;
t38 = -Ifges(3,3) * t106 + t110 * t199 + t162 * t74;
t37 = -t108 * t164 + t112 * t477 + t282 * t76;
t36 = -t107 * t163 + t111 * t478 + t282 * t75;
t35 = -t106 * t162 + t110 * t479 + t282 * t74;
t34 = -Ifges(3,3) * t105 + t109 * t191 + t153 * t73;
t33 = -t105 * t153 + t109 * t480 + t282 * t73;
t32 = -t108 * t201 + t112 * t152 + t76 * t477;
t31 = -t107 * t200 + t111 * t151 + t75 * t478;
t30 = -t106 * t199 + t110 * t150 + t74 * t479;
t29 = -t105 * t191 + t109 * t149 + t73 * t480;
t4 = -t164 * t20 - t90 * t291 * t385 + (-t168 * t16 + (-t101 * t275 - t389 * (t101 + t90)) * t310 + (t104 * t258 + t385 * t557) * t506) * t289 + (-t24 - t222) * t282;
t3 = -t163 * t19 - t89 * t291 * t386 + (-t167 * t15 + (-t100 * t275 - t390 * (t100 + t89)) * t308 + (t103 * t258 + t386 * t558) * t508) * t289 + (-t23 - t221) * t282;
t2 = -t162 * t18 - t88 * t291 * t387 + (-t166 * t14 + (-t275 * t99 - t391 * (t99 + t88)) * t306 + (t102 * t258 + t387 * t559) * t510) * t289 + (-t22 - t220) * t282;
t1 = -t153 * t17 - t85 * t291 * t388 + (-t165 * t13 + (-t275 * t97 - t392 * (t97 + t85)) * t293 + (t258 * t98 + t388 * t560) * t512) * t289 + (-t21 - t219) * t282;
t25 = [-t205 * t301 - t287 * t206 + (t304 - g(1)) * m(4) + ((t124 * t64 + t158 * t56 + t72 * t497) * t302 + (t119 * t64 - t72 * t394 + t56 * t485) * t144 + (t64 * t148 + t4) * t120 + ((-t161 * t56 + t72 * t500) * t148 - t361) * t274) * t136 + ((t123 * t62 + t157 * t54 + t70 * t498) * t302 + (t117 * t62 - t70 * t396 + t54 * t487) * t143 + (t62 * t147 + t3) * t118 + ((-t160 * t54 + t70 * t501) * t147 - t362) * t273) * t135 + ((t122 * t60 + t156 * t52 + t68 * t499) * t302 + (t115 * t60 - t68 * t398 + t52 * t489) * t142 + (t60 * t146 + t2) * t116 + ((-t159 * t52 + t68 * t502) * t146 - t363) * t272) * t134 + ((t121 * t58 + t154 * t50 + t66 * t503) * t302 + (t113 * t58 - t66 * t400 + t50 * t491) * t141 + (t58 * t145 + t1) * t114 + ((-t155 * t50 + t66 * t504) * t145 - t364) * t271) * t133; -t287 * t205 + t206 * t301 + (t303 - g(2)) * m(4) + ((t124 * t63 + t158 * t55 + t497 * t71) * t302 + (t120 * t63 + t393 * t71 - t484 * t55) * t148 + (t63 * t144 + t4) * t119 + ((t161 * t55 - t500 * t71) * t144 + t361) * t270) * t136 + ((t123 * t61 + t157 * t53 + t498 * t69) * t302 + (t118 * t61 + t395 * t69 - t486 * t53) * t147 + (t61 * t143 + t3) * t117 + ((t160 * t53 - t501 * t69) * t143 + t362) * t269) * t135 + ((t122 * t59 + t156 * t51 + t499 * t67) * t302 + (t116 * t59 + t397 * t67 - t488 * t51) * t146 + (t59 * t142 + t2) * t115 + ((t159 * t51 - t502 * t67) * t142 + t363) * t268) * t134 + ((t121 * t57 + t154 * t49 + t503 * t65) * t302 + (t114 * t57 + t399 * t65 - t49 * t490) * t145 + (t57 * t141 + t1) * t113 + ((t155 * t49 - t504 * t65) * t141 + t364) * t267) * t133; (t302 - g(3)) * m(4) + (t158 * t8 + t124 * t4 + (t124 * t84 + t158 * t80 + t497 * t96) * t302 + (t120 * t84 + t393 * t96 - t484 * t80) * t148 + (t119 * t84 - t394 * t96 + t485 * t80) * t144 + t12 * t497) * t136 + (t157 * t7 + t123 * t3 + (t123 * t83 + t157 * t79 + t498 * t95) * t302 + (t118 * t83 + t395 * t95 - t486 * t79) * t147 + (t117 * t83 - t396 * t95 + t487 * t79) * t143 + t11 * t498) * t135 + (t156 * t6 + t122 * t2 + (t122 * t82 + t156 * t78 + t499 * t94) * t302 + (t116 * t82 + t397 * t94 - t488 * t78) * t146 + (t115 * t82 - t398 * t94 + t489 * t78) * t142 + t10 * t499) * t134 + (t154 * t5 + t121 * t1 + (t121 * t81 + t154 * t77 + t503 * t87) * t302 + (t114 * t81 + t399 * t87 - t490 * t77) * t145 + (t113 * t81 - t400 * t87 + t491 * t77) * t141 + t9 * t503) * t133; -(-g(1) * t327 - g(2) * t326) * t279 + t280 * (g(1) * t326 - g(2) * t327) + Ifges(4,3) * t301 + t206 * t303 - t205 * t304 + t112 * t8 - t105 * t9 - t106 * t10 - t107 * t11 - t108 * t12 + t109 * t5 + t110 * t6 + t111 * t7 + t73 * t1 + t74 * t2 + t75 * t3 + t76 * t4 + ((t124 * t37 + t158 * t32 + t40 * t497) * t302 + (t120 * t37 - t32 * t484 + t393 * t40) * t148 + (t119 * t37 + t32 * t485 - t394 * t40) * t144) * t136 + ((t123 * t36 + t157 * t31 + t39 * t498) * t302 + (t118 * t36 - t31 * t486 + t39 * t395) * t147 + (t117 * t36 + t31 * t487 - t39 * t396) * t143) * t135 + ((t122 * t35 + t156 * t30 + t38 * t499) * t302 + (t116 * t35 - t30 * t488 + t38 * t397) * t146 + (t115 * t35 + t30 * t489 - t38 * t398) * t142) * t134 + ((t113 * t33 + t29 * t491 - t34 * t400) * t141 + (t121 * t33 + t154 * t29 + t34 * t503) * t302 + (t114 * t33 - t29 * t490 + t34 * t399) * t145) * t133;];
tauX  = t25;
