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
% Datum: 2020-08-07 17:26
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P4RRRRR2G1A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4RRRRR2G1A0_invdyn_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4RRRRR2G1A0_invdyn_para_pf_slag_vp2: xDP has to be [4x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [4 1]), ...
  'P4RRRRR2G1A0_invdyn_para_pf_slag_vp2: xDDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4RRRRR2G1A0_invdyn_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4RRRRR2G1A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4RRRRR2G1A0_invdyn_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4RRRRR2G1A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4RRRRR2G1A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P4RRRRR2G1A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4RRRRR2G1A0_invdyn_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4RRRRR2G1A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 17:24:27
% EndTime: 2020-08-07 17:24:37
% DurationCPUTime: 9.95s
% Computational Cost: add. (32883->702), mult. (29561->1194), div. (13244->23), fcn. (25218->98), ass. (0->529)
t582 = 2 * pkin(1);
t581 = -2 * pkin(1);
t580 = mrSges(3,1) * g(3);
t579 = mrSges(3,2) * g(3);
t305 = sin(qJ(2,4));
t578 = pkin(1) * t305;
t307 = cos(qJ(2,4));
t577 = pkin(1) * t307;
t315 = sin(qJ(2,3));
t576 = pkin(1) * t315;
t317 = sin(qJ(2,2));
t575 = pkin(1) * t317;
t319 = sin(qJ(2,1));
t574 = pkin(1) * t319;
t321 = cos(qJ(2,3));
t573 = pkin(1) * t321;
t323 = cos(qJ(2,2));
t572 = pkin(1) * t323;
t325 = cos(qJ(2,1));
t571 = pkin(1) * t325;
t306 = cos(qJ(3,4));
t276 = t306 ^ 2;
t570 = pkin(2) * t276;
t320 = cos(qJ(3,3));
t283 = t320 ^ 2;
t569 = pkin(2) * t283;
t322 = cos(qJ(3,2));
t287 = t322 ^ 2;
t568 = pkin(2) * t287;
t324 = cos(qJ(3,1));
t291 = t324 ^ 2;
t567 = pkin(2) * t291;
t566 = Ifges(3,1) + Ifges(2,3);
t300 = legFrame(4,3);
t258 = sin(t300);
t262 = cos(t300);
t175 = -g(1) * t258 + g(2) * t262;
t565 = mrSges(3,1) * t175;
t301 = legFrame(3,3);
t259 = sin(t301);
t263 = cos(t301);
t176 = -g(1) * t259 + g(2) * t263;
t564 = mrSges(3,1) * t176;
t302 = legFrame(2,3);
t260 = sin(t302);
t264 = cos(t302);
t177 = -g(1) * t260 + g(2) * t264;
t563 = mrSges(3,1) * t177;
t303 = legFrame(1,3);
t261 = sin(t303);
t265 = cos(t303);
t178 = -g(1) * t261 + g(2) * t265;
t562 = mrSges(3,1) * t178;
t561 = mrSges(3,2) * t175;
t560 = mrSges(3,2) * t176;
t559 = mrSges(3,2) * t177;
t558 = mrSges(3,2) * t178;
t179 = g(1) * t262 + g(2) * t258;
t557 = mrSges(3,2) * t179;
t180 = g(1) * t263 + g(2) * t259;
t556 = mrSges(3,2) * t180;
t181 = g(1) * t264 + g(2) * t260;
t555 = mrSges(3,2) * t181;
t182 = g(1) * t265 + g(2) * t261;
t554 = mrSges(3,2) * t182;
t304 = sin(qJ(3,4));
t553 = mrSges(3,2) * t304;
t314 = sin(qJ(3,3));
t552 = mrSges(3,2) * t314;
t316 = sin(qJ(3,2));
t551 = mrSges(3,2) * t316;
t318 = sin(qJ(3,1));
t550 = mrSges(3,2) * t318;
t549 = Ifges(3,4) * t304;
t548 = Ifges(3,4) * t314;
t547 = Ifges(3,4) * t316;
t546 = Ifges(3,4) * t318;
t278 = 0.1e1 / t306 ^ 2;
t545 = Ifges(3,6) * t278;
t285 = 0.1e1 / t320 ^ 2;
t544 = Ifges(3,6) * t285;
t289 = 0.1e1 / t322 ^ 2;
t543 = Ifges(3,6) * t289;
t293 = 0.1e1 / t324 ^ 2;
t542 = Ifges(3,6) * t293;
t327 = xDP(3);
t296 = t327 ^ 2;
t495 = t296 / pkin(2) ^ 2;
t405 = t495 / 0.2e1;
t346 = 1 / pkin(1);
t343 = 0.1e1 / pkin(2);
t509 = t278 * t343;
t433 = (pkin(2) * t306 + t577) * t509;
t275 = 0.1e1 / t305;
t512 = t275 * t304;
t371 = t433 * t512;
t358 = t346 * t371;
t354 = t327 * t358;
t331 = xP(4);
t272 = sin(t331);
t273 = cos(t331);
t334 = koppelP(4,2);
t338 = koppelP(4,1);
t171 = -t272 * t334 + t273 * t338;
t326 = xDP(4);
t328 = xDP(2);
t145 = t171 * t326 + t328;
t468 = qJ(2,4) + qJ(3,4);
t247 = qJ(1,4) + t468;
t214 = t300 + t247;
t469 = qJ(2,4) - qJ(3,4);
t248 = qJ(1,4) + t469;
t215 = t300 + t248;
t244 = qJ(1,4) + t300;
t133 = sin(t244) * t581 + (-sin(t215) - sin(t214)) * pkin(2);
t187 = 0.1e1 / (sin(t468) + sin(t469));
t476 = t343 * t346;
t440 = t187 * t476;
t403 = t133 * t440;
t97 = t145 * t403;
t167 = t272 * t338 + t273 * t334;
t329 = xDP(1);
t149 = -t167 * t326 + t329;
t134 = cos(t244) * t581 + (-cos(t214) - cos(t215)) * pkin(2);
t402 = t134 * t440;
t98 = t149 * t402;
t57 = t97 + t98 - t354;
t277 = 0.1e1 / t306;
t511 = t275 * t346;
t421 = t304 * t511;
t395 = t277 * t421;
t274 = qJ(1,4) + qJ(2,4);
t235 = t300 + t274;
t213 = cos(t235);
t428 = t213 * t511;
t212 = sin(t235);
t429 = t212 * t511;
t90 = t145 * t429 + t149 * t428 + t327 * t395;
t540 = (t278 * t405 + (t90 + t57 / 0.2e1) * t57) * t305;
t335 = koppelP(3,2);
t339 = koppelP(3,1);
t168 = t272 * t339 + t273 * t335;
t150 = -t168 * t326 + t329;
t249 = qJ(1,3) + t301;
t236 = qJ(2,3) + t249;
t225 = qJ(3,3) + t236;
t226 = -qJ(3,3) + t236;
t138 = cos(t249) * t581 + (-cos(t225) - cos(t226)) * pkin(2);
t470 = qJ(2,3) + qJ(3,3);
t471 = qJ(2,3) - qJ(3,3);
t190 = 0.1e1 / (sin(t470) + sin(t471));
t438 = t190 * t476;
t398 = t138 * t438;
t102 = t150 * t398;
t501 = t285 * t343;
t432 = (pkin(2) * t320 + t573) * t501;
t280 = 0.1e1 / t315;
t508 = t280 * t314;
t370 = t432 * t508;
t357 = t346 * t370;
t353 = t327 * t357;
t172 = -t272 * t335 + t273 * t339;
t146 = t172 * t326 + t328;
t135 = sin(t249) * t581 + (-sin(t226) - sin(t225)) * pkin(2);
t401 = t135 * t438;
t99 = t146 * t401;
t59 = t102 + t99 - t353;
t284 = 0.1e1 / t320;
t507 = t280 * t346;
t419 = t314 * t507;
t393 = t284 * t419;
t222 = cos(t236);
t424 = t222 * t507;
t219 = sin(t236);
t427 = t219 * t507;
t94 = t146 * t427 + t150 * t424 + t327 * t393;
t539 = (t285 * t405 + (t94 + t59 / 0.2e1) * t59) * t315;
t336 = koppelP(2,2);
t340 = koppelP(2,1);
t173 = -t272 * t336 + t273 * t340;
t147 = t173 * t326 + t328;
t472 = qJ(2,2) + qJ(3,2);
t268 = qJ(1,2) + t472;
t227 = t302 + t268;
t473 = qJ(2,2) - qJ(3,2);
t269 = qJ(1,2) + t473;
t228 = t302 + t269;
t250 = qJ(1,2) + t302;
t136 = sin(t250) * t581 + (-sin(t228) - sin(t227)) * pkin(2);
t191 = 0.1e1 / (sin(t472) + sin(t473));
t436 = t191 * t476;
t400 = t136 * t436;
t100 = t147 * t400;
t169 = t272 * t340 + t273 * t336;
t151 = -t169 * t326 + t329;
t139 = cos(t250) * t581 + (-cos(t227) - cos(t228)) * pkin(2);
t397 = t139 * t436;
t103 = t151 * t397;
t499 = t289 * t343;
t431 = (pkin(2) * t322 + t572) * t499;
t281 = 0.1e1 / t317;
t506 = t281 * t316;
t369 = t431 * t506;
t356 = t346 * t369;
t352 = t327 * t356;
t60 = t100 + t103 - t352;
t288 = 0.1e1 / t322;
t505 = t281 * t346;
t418 = t316 * t505;
t392 = t288 * t418;
t298 = qJ(2,2) + qJ(1,2);
t237 = t302 + t298;
t223 = cos(t237);
t423 = t223 * t505;
t220 = sin(t237);
t426 = t220 * t505;
t95 = t147 * t426 + t151 * t423 + t327 * t392;
t538 = (t289 * t405 + (t95 + t60 / 0.2e1) * t60) * t317;
t337 = koppelP(1,2);
t341 = koppelP(1,1);
t174 = -t272 * t337 + t273 * t341;
t148 = t174 * t326 + t328;
t251 = qJ(1,1) + t303;
t238 = qJ(2,1) + t251;
t229 = qJ(3,1) + t238;
t230 = -qJ(3,1) + t238;
t137 = sin(t251) * t581 + (-sin(t230) - sin(t229)) * pkin(2);
t474 = qJ(2,1) + qJ(3,1);
t475 = qJ(2,1) - qJ(3,1);
t192 = 0.1e1 / (sin(t474) + sin(t475));
t434 = t192 * t476;
t399 = t137 * t434;
t101 = t148 * t399;
t170 = t272 * t341 + t273 * t337;
t152 = -t170 * t326 + t329;
t140 = cos(t251) * t581 + (-cos(t229) - cos(t230)) * pkin(2);
t396 = t140 * t434;
t104 = t152 * t396;
t497 = t293 * t343;
t430 = (pkin(2) * t324 + t571) * t497;
t282 = 0.1e1 / t319;
t504 = t282 * t318;
t368 = t430 * t504;
t355 = t346 * t368;
t351 = t327 * t355;
t61 = t101 + t104 - t351;
t292 = 0.1e1 / t324;
t503 = t282 * t346;
t417 = t318 * t503;
t391 = t292 * t417;
t224 = cos(t238);
t422 = t224 * t503;
t221 = sin(t238);
t425 = t221 * t503;
t96 = t148 * t425 + t152 * t422 + t327 * t391;
t537 = (t293 * t405 + (t96 + t61 / 0.2e1) * t61) * t319;
t40 = t57 + t90;
t536 = t306 * t40;
t51 = t59 + t94;
t535 = t320 * t51;
t52 = t60 + t95;
t534 = t322 * t52;
t53 = t61 + t96;
t533 = t324 * t53;
t295 = t326 ^ 2;
t310 = xDDP(4);
t312 = xDDP(2);
t532 = (-t167 * t295 + t171 * t310 + t312) * t346;
t531 = (-t168 * t295 + t172 * t310 + t312) * t346;
t530 = (-t169 * t295 + t173 * t310 + t312) * t346;
t529 = (-t170 * t295 + t174 * t310 + t312) * t346;
t313 = xDDP(1);
t528 = (-t167 * t310 - t171 * t295 + t313) * t346;
t527 = (-t168 * t310 - t172 * t295 + t313) * t346;
t526 = (-t169 * t310 - t173 * t295 + t313) * t346;
t525 = (-t170 * t310 - t174 * t295 + t313) * t346;
t524 = t187 * t343;
t523 = t190 * t343;
t522 = t191 * t343;
t521 = t192 * t343;
t520 = t212 * t275;
t519 = t213 * t275;
t518 = t219 * t280;
t517 = t220 * t281;
t516 = t221 * t282;
t515 = t222 * t280;
t514 = t223 * t281;
t513 = t224 * t282;
t510 = t277 * t343;
t502 = t284 * t343;
t500 = t288 * t343;
t498 = t292 * t343;
t496 = t296 * t343;
t494 = t304 * t305;
t309 = mrSges(3,3) - mrSges(2,2);
t493 = t307 * t309;
t308 = Ifges(3,2) - Ifges(3,1);
t492 = t308 * t276;
t491 = t308 * t283;
t490 = t308 * t287;
t489 = t308 * t291;
t246 = cos(t274);
t488 = t309 * t246;
t297 = qJ(1,3) + qJ(2,3);
t255 = cos(t297);
t487 = t309 * t255;
t256 = cos(t298);
t486 = t309 * t256;
t299 = qJ(1,1) + qJ(2,1);
t257 = cos(t299);
t485 = t309 * t257;
t484 = t309 * t321;
t483 = t309 * t323;
t482 = t309 * t325;
t311 = xDDP(3);
t481 = t311 * t346;
t480 = t314 * t315;
t479 = t316 * t317;
t478 = t318 * t319;
t477 = t327 * t343;
t467 = 0.2e1 * t549;
t466 = 0.2e1 * t548;
t465 = 0.2e1 * t547;
t464 = 0.2e1 * t546;
t89 = t90 ^ 2;
t463 = t89 * t577;
t91 = t94 ^ 2;
t462 = t91 * t573;
t92 = t95 ^ 2;
t461 = t92 * t572;
t93 = t96 ^ 2;
t460 = t93 * t571;
t459 = t307 * t570;
t458 = t321 * t569;
t457 = t323 * t568;
t456 = t325 * t567;
t455 = Ifges(3,5) * t495;
t454 = Ifges(3,3) * t495;
t453 = t133 * t524;
t452 = t134 * t524;
t451 = t135 * t523;
t450 = t136 * t522;
t449 = t137 * t521;
t448 = t138 * t523;
t447 = t139 * t522;
t446 = t140 * t521;
t387 = t492 + t566;
t153 = t306 * t467 + t387;
t445 = t153 * t524;
t386 = t491 + t566;
t154 = t320 * t466 + t386;
t444 = t154 * t523;
t385 = t490 + t566;
t155 = t322 * t465 + t385;
t443 = t155 * t522;
t384 = t489 + t566;
t156 = t324 * t464 + t384;
t442 = t156 * t521;
t199 = Ifges(3,5) * t304 + Ifges(3,6) * t306;
t441 = t199 * t524;
t202 = Ifges(3,5) * t314 + Ifges(3,6) * t320;
t439 = t202 * t523;
t203 = Ifges(3,5) * t316 + Ifges(3,6) * t322;
t437 = t203 * t522;
t204 = Ifges(3,5) * t318 + Ifges(3,6) * t324;
t435 = t204 * t521;
t420 = t277 * t477;
t416 = t284 * t477;
t415 = t288 * t477;
t414 = t292 * t477;
t413 = t304 * t495;
t412 = t314 * t495;
t411 = t316 * t495;
t410 = t318 * t495;
t409 = t327 * t494;
t408 = t327 * t480;
t407 = t327 * t479;
t406 = t327 * t478;
t404 = 0.4e1 * Ifges(3,4) * t477;
t394 = t307 * t420;
t390 = t321 * t416;
t389 = t323 * t415;
t388 = t325 * t414;
t383 = -0.2e1 * t40 * t420;
t382 = -0.2e1 * t51 * t416;
t381 = -0.2e1 * t52 * t415;
t380 = -0.2e1 * t53 * t414;
t379 = t40 * t394;
t378 = t51 * t390;
t377 = t52 * t389;
t376 = t53 * t388;
t375 = mrSges(3,1) * t306 - t553;
t374 = mrSges(3,1) * t320 - t552;
t373 = mrSges(3,1) * t322 - t551;
t372 = mrSges(3,1) * t324 - t550;
t330 = (m(2) + m(3));
t345 = pkin(1) ^ 2;
t367 = (t330 * t345) + Ifges(1,3) + t566;
t366 = t404 * t536 + Ifges(3,4) * t383 + (t308 * t304 * t383 + t278 * t455) * t306;
t365 = t404 * t535 + Ifges(3,4) * t382 + (t308 * t314 * t382 + t285 * t455) * t320;
t364 = t404 * t534 + Ifges(3,4) * t381 + (t308 * t316 * t381 + t289 * t455) * t322;
t363 = t404 * t533 + Ifges(3,4) * t380 + (t308 * t318 * t380 + t293 * t455) * t324;
t362 = mrSges(2,1) + t375;
t361 = mrSges(2,1) + t374;
t360 = mrSges(2,1) + t373;
t359 = mrSges(2,1) + t372;
t350 = t324 * t291;
t349 = t322 * t287;
t348 = t320 * t283;
t347 = t306 * t276;
t342 = pkin(2) ^ 2;
t333 = mrSges(4,1);
t332 = mrSges(4,2);
t294 = 0.1e1 / t350;
t290 = 0.1e1 / t349;
t286 = 0.1e1 / t348;
t279 = 0.1e1 / t347;
t271 = qJ(1,1) + t475;
t270 = qJ(1,1) + t474;
t267 = qJ(1,3) + t471;
t266 = qJ(1,3) + t470;
t254 = sin(t299);
t253 = sin(t298);
t252 = sin(t297);
t245 = sin(t274);
t243 = mrSges(3,1) * t571;
t242 = mrSges(3,1) * t572;
t241 = mrSges(3,1) * t573;
t240 = mrSges(3,1) * t577;
t239 = pkin(1) * t330 + mrSges(1,1);
t234 = t309 * t574;
t233 = t309 * t575;
t232 = t309 * t576;
t231 = t309 * t578;
t195 = (-mrSges(2,1) + t550) * t571;
t194 = (-mrSges(2,1) + t551) * t572;
t193 = (-mrSges(2,1) + t552) * t573;
t188 = (-mrSges(2,1) + t553) * t577;
t166 = mrSges(2,1) * t182;
t165 = mrSges(2,1) * t181;
t164 = mrSges(2,1) * t180;
t163 = mrSges(2,1) * t179;
t162 = mrSges(3,1) * t182;
t161 = mrSges(3,1) * t181;
t160 = mrSges(3,1) * t180;
t159 = mrSges(3,1) * t179;
t158 = -t272 * t332 + t273 * t333;
t157 = t272 * t333 + t273 * t332;
t144 = t324 * (-mrSges(3,2) * t574 + Ifges(3,6)) - t318 * (mrSges(3,1) * t574 - Ifges(3,5));
t143 = t322 * (-mrSges(3,2) * t575 + Ifges(3,6)) - t316 * (mrSges(3,1) * t575 - Ifges(3,5));
t142 = t320 * (-mrSges(3,2) * t576 + Ifges(3,6)) - t314 * (mrSges(3,1) * t576 - Ifges(3,5));
t141 = t306 * (-mrSges(3,2) * t578 + Ifges(3,6)) - t304 * (mrSges(3,1) * t578 - Ifges(3,5));
t116 = (t243 + t464) * t324 - t195 + t234 + t384;
t115 = (t242 + t465) * t322 - t194 + t233 + t385;
t114 = (t241 + t466) * t320 - t193 + t232 + t386;
t113 = (t240 + t467) * t306 - t188 + t231 + t387;
t112 = t489 + 0.2e1 * (t243 + t546) * t324 - 0.2e1 * t195 + 0.2e1 * t234 + t367;
t111 = t490 + 0.2e1 * (t242 + t547) * t322 - 0.2e1 * t194 + 0.2e1 * t233 + t367;
t110 = t491 + 0.2e1 * (t241 + t548) * t320 - 0.2e1 * t193 + 0.2e1 * t232 + t367;
t109 = t492 + 0.2e1 * (t240 + t549) * t306 - 0.2e1 * t188 + 0.2e1 * t231 + t367;
t108 = (-t170 * t224 + t174 * t221) * t503;
t107 = (-t169 * t223 + t173 * t220) * t505;
t106 = (-t168 * t222 + t172 * t219) * t507;
t105 = (-t167 * t213 + t171 * t212) * t511;
t88 = t204 * t498 + (t116 * t292 - t156 * t430) * t417;
t87 = t203 * t500 + (t115 * t288 - t155 * t431) * t418;
t86 = t202 * t502 + (t114 * t284 - t154 * t432) * t419;
t85 = t199 * t510 + (t113 * t277 - t153 * t433) * t421;
t84 = (t116 * t513 + t140 * t442) * t346;
t83 = (t115 * t514 + t139 * t443) * t346;
t82 = (t114 * t515 + t138 * t444) * t346;
t81 = (t116 * t516 + t137 * t442) * t346;
t80 = (t115 * t517 + t136 * t443) * t346;
t79 = (t114 * t518 + t135 * t444) * t346;
t78 = (t137 * t174 - t140 * t170) * t434;
t77 = (t136 * t173 - t139 * t169) * t436;
t76 = (t135 * t172 - t138 * t168) * t438;
t75 = (t113 * t519 + t134 * t445) * t346;
t74 = (t113 * t520 + t133 * t445) * t346;
t73 = (t133 * t171 - t134 * t167) * t440;
t72 = (t112 * t513 + t116 * t446) * t346;
t71 = (t111 * t514 + t115 * t447) * t346;
t70 = (t110 * t515 + t114 * t448) * t346;
t69 = (t112 * t516 + t116 * t449) * t346;
t68 = (t111 * t517 + t115 * t450) * t346;
t67 = (t110 * t518 + t114 * t451) * t346;
t66 = t144 * t498 + (t112 * t292 - t116 * t430) * t417;
t65 = t143 * t500 + (t111 * t288 - t115 * t431) * t418;
t64 = t142 * t502 + (t110 * t284 - t114 * t432) * t419;
t63 = (t109 * t519 + t113 * t452) * t346;
t62 = (t109 * t520 + t113 * t453) * t346;
t58 = t141 * t510 + (t109 * t277 - t113 * t433) * t421;
t56 = t108 * t116 + t156 * t78;
t55 = t107 * t115 + t155 * t77;
t54 = t106 * t114 + t154 * t76;
t50 = t53 ^ 2;
t49 = t52 ^ 2;
t48 = t51 ^ 2;
t47 = t104 / 0.2e1 + t101 / 0.2e1 - t351 / 0.2e1 + t96;
t46 = t103 / 0.2e1 + t100 / 0.2e1 - t352 / 0.2e1 + t95;
t45 = t102 / 0.2e1 + t99 / 0.2e1 - t353 / 0.2e1 + t94;
t44 = t342 * t53 * t350;
t43 = t342 * t52 * t349;
t42 = t342 * t51 * t348;
t41 = t105 * t113 + t153 * t73;
t39 = t40 ^ 2;
t38 = t98 / 0.2e1 + t97 / 0.2e1 - t354 / 0.2e1 + t90;
t31 = t342 * t40 * t347;
t24 = t108 * t112 + t116 * t78;
t23 = t107 * t111 + t115 * t77;
t22 = t106 * t110 + t114 * t76;
t21 = t105 * t109 + t113 * t73;
t16 = ((-t324 * t571 * t96 - t53 * t567) * t292 * t96 - pkin(2) * t61 * t533 - t294 * t496) * t503;
t15 = ((-t322 * t572 * t95 - t52 * t568) * t288 * t95 - pkin(2) * t60 * t534 - t290 * t496) * t505;
t14 = ((-t320 * t573 * t94 - t51 * t569) * t284 * t94 - pkin(2) * t59 * t535 - t286 * t496) * t507;
t13 = ((-t306 * t577 * t90 - t40 * t570) * t277 * t90 - pkin(2) * t57 * t536 - t279 * t496) * t511;
t12 = (((-pkin(1) * t53 * t478 + t292 * t327) * t324 + pkin(1) * t388) * t294 * t477 + ((t44 + t47 * t456 * t582 + (-pkin(1) * t292 * t406 + t345 * t96) * t324) * t96 + (t44 + (t456 * t53 - t406) * pkin(1)) * t61) * t497) * t503;
t11 = (((-pkin(1) * t52 * t479 + t288 * t327) * t322 + pkin(1) * t389) * t290 * t477 + ((t43 + t46 * t457 * t582 + (-pkin(1) * t288 * t407 + t345 * t95) * t322) * t95 + (t43 + (t457 * t52 - t407) * pkin(1)) * t60) * t499) * t505;
t10 = (((-pkin(1) * t51 * t480 + t284 * t327) * t320 + pkin(1) * t390) * t286 * t477 + ((t42 + t45 * t458 * t582 + (-pkin(1) * t284 * t408 + t345 * t94) * t320) * t94 + (t42 + (t458 * t51 - t408) * pkin(1)) * t59) * t501) * t507;
t9 = (((-pkin(1) * t40 * t494 + t277 * t327) * t306 + pkin(1) * t394) * t279 * t477 + ((t31 + t38 * t459 * t582 + (-pkin(1) * t277 * t409 + t345 * t90) * t306) * t90 + (t31 + (t40 * t459 - t409) * pkin(1)) * t57) * t509) * t511;
t8 = -t116 * t16 - t156 * t12 + t166 * t254 + (t204 * t294 - t542) * t410 + (t254 * t372 - t485) * t182 + (-t309 * t254 - t257 * t359) * t178 + (t319 * t359 - t482) * t93 * pkin(1) + t363;
t7 = -t115 * t15 - t155 * t11 + t165 * t253 + (t203 * t290 - t543) * t411 + (t253 * t373 - t486) * t181 + (-t309 * t253 - t256 * t360) * t177 + (t317 * t360 - t483) * t92 * pkin(1) + t364;
t6 = -t114 * t14 - t154 * t10 + t164 * t252 + (t202 * t286 - t544) * t412 + (t252 * t374 - t487) * t180 + (-t309 * t252 - t255 * t361) * t176 + (t315 * t361 - t484) * t91 * pkin(1) + t365;
t5 = -t113 * t13 - t153 * t9 + t163 * t245 + (t199 * t279 - t545) * t413 + (t245 * t375 - t488) * t179 + (-t309 * t245 - t246 * t362) * t175 + (t305 * t362 - t493) * t89 * pkin(1) + t366;
t4 = -t112 * t16 - t116 * t12 + (-t554 - t562) * cos(t271) / 0.2e1 + (t162 - t558) * sin(t271) / 0.2e1 + (t554 - t562) * cos(t270) / 0.2e1 + (t162 + t558) * sin(t270) / 0.2e1 - mrSges(2,1) * t178 * t257 + (mrSges(1,2) * t178 + t182 * t239) * sin(qJ(1,1)) + (t144 * t294 - t542) * t410 - t182 * t485 + (-t178 * t309 + t166) * t254 + (mrSges(1,2) * t182 - t178 * t239) * cos(qJ(1,1)) + ((-mrSges(3,1) * t537 - mrSges(3,2) * t376) * t324 + (-mrSges(3,1) * t376 + mrSges(3,2) * t537) * t318 + (-mrSges(2,1) * t319 + t482) * t61 * t47) * t582 + t363;
t3 = -t111 * t15 - t115 * t11 + (-t555 - t563) * cos(t269) / 0.2e1 + (t161 - t559) * sin(t269) / 0.2e1 + (t555 - t563) * cos(t268) / 0.2e1 + (t161 + t559) * sin(t268) / 0.2e1 - mrSges(2,1) * t177 * t256 + (mrSges(1,2) * t177 + t181 * t239) * sin(qJ(1,2)) + (t143 * t290 - t543) * t411 - t181 * t486 + (-t177 * t309 + t165) * t253 + (mrSges(1,2) * t181 - t177 * t239) * cos(qJ(1,2)) + ((-mrSges(3,1) * t538 - mrSges(3,2) * t377) * t322 + (-mrSges(3,1) * t377 + mrSges(3,2) * t538) * t316 + (-mrSges(2,1) * t317 + t483) * t60 * t46) * t582 + t364;
t2 = -t110 * t14 - t114 * t10 + (-t556 - t564) * cos(t267) / 0.2e1 + (t160 - t560) * sin(t267) / 0.2e1 + (t556 - t564) * cos(t266) / 0.2e1 + (t160 + t560) * sin(t266) / 0.2e1 - mrSges(2,1) * t176 * t255 + (mrSges(1,2) * t176 + t180 * t239) * sin(qJ(1,3)) + (t142 * t286 - t544) * t412 - t180 * t487 + (-t176 * t309 + t164) * t252 + (mrSges(1,2) * t180 - t176 * t239) * cos(qJ(1,3)) + ((-mrSges(3,1) * t539 - mrSges(3,2) * t378) * t320 + (-mrSges(3,1) * t378 + mrSges(3,2) * t539) * t314 + (-mrSges(2,1) * t315 + t484) * t59 * t45) * t582 + t365;
t1 = -t109 * t13 - t113 * t9 + (-t557 - t565) * cos(t248) / 0.2e1 + (t159 - t561) * sin(t248) / 0.2e1 + (t557 - t565) * cos(t247) / 0.2e1 + (t159 + t561) * sin(t247) / 0.2e1 - mrSges(2,1) * t175 * t246 + (mrSges(1,2) * t175 + t179 * t239) * sin(qJ(1,4)) + (t141 * t279 - t545) * t413 - t179 * t488 + (-t175 * t309 + t163) * t245 + (mrSges(1,2) * t179 - t175 * t239) * cos(qJ(1,4)) + ((-mrSges(3,1) * t540 - mrSges(3,2) * t379) * t306 + (-mrSges(3,1) * t379 + mrSges(3,2) * t540) * t304 + (-mrSges(2,1) * t305 + t493) * t57 * t38) * t582 + t366;
t17 = [(t446 * t84 + t513 * t72) * t525 + (t449 * t84 + t516 * t72) * t529 + t4 * t422 + t8 * t396 + (t447 * t83 + t514 * t71) * t526 + (t450 * t83 + t517 * t71) * t530 + t3 * t423 + t7 * t397 + (t448 * t82 + t515 * t70) * t527 + (t451 * t82 + t518 * t70) * t531 + t2 * t424 + t6 * t398 + (t452 * t75 + t519 * t63) * t528 + (t453 * t75 + t520 * t63) * t532 + t1 * t428 + t5 * t402 - t157 * t310 - t295 * t158 + (t313 - g(1)) * m(4) + (-t84 * t368 + (t72 * t504 + (t140 * t435 + t144 * t513) * t343) * t292 - t83 * t369 + (t71 * t506 + (t139 * t437 + t143 * t514) * t343) * t288 - t82 * t370 + (t70 * t508 + (t138 * t439 + t142 * t515) * t343) * t284 - t75 * t371 + (t63 * t512 + (t134 * t441 + t141 * t519) * t343) * t277) * t481; (t446 * t81 + t513 * t69) * t525 + (t449 * t81 + t516 * t69) * t529 + t4 * t425 + t8 * t399 + (t447 * t80 + t514 * t68) * t526 + (t450 * t80 + t517 * t68) * t530 + t3 * t426 + t7 * t400 + (t448 * t79 + t515 * t67) * t527 + (t451 * t79 + t518 * t67) * t531 + t2 * t427 + t6 * t401 + (t452 * t74 + t519 * t62) * t528 + (t453 * t74 + t520 * t62) * t532 + t1 * t429 + t5 * t403 + t158 * t310 - t295 * t157 + (t312 - g(2)) * m(4) + (-t81 * t368 + (t69 * t504 + (t137 * t435 + t144 * t516) * t343) * t292 - t80 * t369 + (t68 * t506 + (t136 * t437 + t143 * t517) * t343) * t288 - t79 * t370 + (t67 * t508 + (t135 * t439 + t142 * t518) * t343) * t284 - t74 * t371 + (t62 * t512 + (t133 * t441 + t141 * t520) * t343) * t277) * t481; (-t141 * t13 - t199 * t9 + (mrSges(3,2) * t463 - t580) * t306 + (t306 * t308 * t39 + mrSges(3,1) * t463 + t279 * t454 + t579) * t304 + (t175 * t245 + t179 * t246) * (mrSges(3,1) * t304 + mrSges(3,2) * t306) + (-0.2e1 * t276 + 0.1e1) * t39 * Ifges(3,4)) * t510 + (-t202 * t10 - t142 * t14 + (mrSges(3,2) * t462 - t580) * t320 + (t308 * t320 * t48 + mrSges(3,1) * t462 + t286 * t454 + t579) * t314 + (t176 * t252 + t180 * t255) * (mrSges(3,1) * t314 + mrSges(3,2) * t320) + (-0.2e1 * t283 + 0.1e1) * t48 * Ifges(3,4)) * t502 + (-t204 * t12 - t144 * t16 + (mrSges(3,2) * t460 - t580) * t324 + (t308 * t324 * t50 + mrSges(3,1) * t460 + t294 * t454 + t579) * t318 + (t178 * t254 + t182 * t257) * (mrSges(3,1) * t318 + mrSges(3,2) * t324) + (-0.2e1 * t291 + 0.1e1) * t50 * Ifges(3,4)) * t498 + (-t203 * t11 - t143 * t15 + (mrSges(3,2) * t461 - t580) * t322 + (t308 * t322 * t49 + mrSges(3,1) * t461 + t290 * t454 + t579) * t316 + (t177 * t253 + t181 * t256) * (mrSges(3,1) * t316 + mrSges(3,2) * t322) + (-0.2e1 * t287 + 0.1e1) * t49 * Ifges(3,4)) * t500 - t6 * t357 - t5 * t358 - t8 * t355 - t7 * t356 + (t446 * t88 + t513 * t66) * t525 + (t447 * t87 + t514 * t65) * t526 + t4 * t391 + t3 * t392 + t2 * t393 + t1 * t395 + (t448 * t86 + t515 * t64) * t527 + (t452 * t85 + t519 * t58) * t528 + (t449 * t88 + t516 * t66) * t529 + (t450 * t87 + t517 * t65) * t530 + (t451 * t86 + t518 * t64) * t531 + (t453 * t85 + t520 * t58) * t532 - g(3) * m(4) + (m(4) + (t66 * t292 - t88 * t430 + (t144 * t292 - t204 * t430) * t498) * t417 + (t65 * t288 - t87 * t431 + (t143 * t288 - t203 * t431) * t500) * t418 + (t64 * t284 - t86 * t432 + (t142 * t284 - t202 * t432) * t502) * t419 + (t58 * t277 - t85 * t433 + (t141 * t277 - t199 * t433) * t510) * t421 + (t277 ^ 2 + t284 ^ 2 + t288 ^ 2 + t292 ^ 2) * Ifges(3,3) * t343 ^ 2) * t311; t273 * (g(1) * t332 - g(2) * t333) - (-g(1) * t333 - g(2) * t332) * t272 + t158 * t312 - t157 * t313 + Ifges(4,3) * t310 + t106 * t2 + t107 * t3 + t108 * t4 + t105 * t1 + (t24 * t513 + t446 * t56) * t525 + t76 * t6 + t77 * t7 + t78 * t8 + t73 * t5 + (t23 * t514 + t447 * t55) * t526 + (t22 * t515 + t448 * t54) * t527 + (t21 * t519 + t41 * t452) * t528 + (t24 * t516 + t449 * t56) * t529 + (t23 * t517 + t450 * t55) * t530 + (t22 * t518 + t451 * t54) * t531 + (t21 * t520 + t41 * t453) * t532 + ((t108 * t144 + t204 * t78) * t498 + (t24 * t292 - t430 * t56) * t417 + (t107 * t143 + t203 * t77) * t500 + (t23 * t288 - t431 * t55) * t418 + (t106 * t142 + t202 * t76) * t502 + (t22 * t284 - t432 * t54) * t419 + (t105 * t141 + t199 * t73) * t510 + (t21 * t277 - t41 * t433) * t421) * t311;];
tauX  = t17;
