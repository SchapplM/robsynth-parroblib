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
% Datum: 2020-08-07 11:22
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P4PRRRR8V2G2A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_slag_vp2: xDP has to be [4x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [4 1]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_slag_vp2: xDDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G2A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:16:39
% EndTime: 2020-08-07 11:16:56
% DurationCPUTime: 17.04s
% Computational Cost: add. (79027->756), mult. (170417->1342), div. (9248->9), fcn. (162006->30), ass. (0->471)
t309 = sin(qJ(2,1));
t315 = cos(qJ(2,1));
t322 = pkin(7) + pkin(6);
t231 = pkin(2) * t309 - t315 * t322;
t288 = sin(pkin(4));
t290 = cos(pkin(4));
t308 = sin(qJ(3,1));
t431 = t290 * t308;
t172 = pkin(3) * t431 + t231 * t288;
t314 = cos(qJ(3,1));
t450 = t288 * t309;
t285 = t314 ^ 2;
t528 = pkin(3) * t285;
t136 = 0.1e1 / (pkin(2) * t431 + t172 * t314 + t450 * t528);
t307 = sin(qJ(2,2));
t313 = cos(qJ(2,2));
t230 = pkin(2) * t307 - t313 * t322;
t306 = sin(qJ(3,2));
t433 = t290 * t306;
t171 = pkin(3) * t433 + t230 * t288;
t312 = cos(qJ(3,2));
t452 = t288 * t307;
t284 = t312 ^ 2;
t529 = pkin(3) * t284;
t135 = 0.1e1 / (pkin(2) * t433 + t171 * t312 + t452 * t529);
t305 = sin(qJ(2,3));
t311 = cos(qJ(2,3));
t229 = pkin(2) * t305 - t311 * t322;
t304 = sin(qJ(3,3));
t435 = t290 * t304;
t170 = pkin(3) * t435 + t229 * t288;
t310 = cos(qJ(3,3));
t454 = t288 * t305;
t283 = t310 ^ 2;
t530 = pkin(3) * t283;
t134 = 0.1e1 / (pkin(2) * t435 + t170 * t310 + t454 * t530);
t292 = sin(qJ(2,4));
t294 = cos(qJ(2,4));
t227 = pkin(2) * t292 - t294 * t322;
t291 = sin(qJ(3,4));
t438 = t290 * t291;
t169 = pkin(3) * t438 + t227 * t288;
t293 = cos(qJ(3,4));
t458 = t288 * t292;
t280 = t293 ^ 2;
t531 = pkin(3) * t280;
t133 = 0.1e1 / (pkin(2) * t438 + t169 * t293 + t458 * t531);
t299 = legFrame(1,2);
t269 = sin(t299);
t273 = cos(t299);
t222 = g(1) * t269 + g(2) * t273;
t226 = g(1) * t273 - g(2) * t269;
t287 = sin(pkin(8));
t259 = g(3) * t287;
t289 = cos(pkin(8));
t376 = t226 * t289 - t259;
t562 = t222 * t288 + t376 * t290;
t298 = legFrame(2,2);
t268 = sin(t298);
t272 = cos(t298);
t221 = g(1) * t268 + g(2) * t272;
t225 = g(1) * t272 - g(2) * t268;
t378 = t225 * t289 - t259;
t561 = t221 * t288 + t378 * t290;
t297 = legFrame(3,2);
t267 = sin(t297);
t271 = cos(t297);
t220 = g(1) * t267 + g(2) * t271;
t224 = g(1) * t271 - g(2) * t267;
t380 = t224 * t289 - t259;
t560 = t220 * t288 + t380 * t290;
t296 = legFrame(4,2);
t266 = sin(t296);
t270 = cos(t296);
t219 = g(1) * t266 + g(2) * t270;
t223 = g(1) * t270 - g(2) * t266;
t382 = t223 * t289 - t259;
t559 = t219 * t288 + t382 * t290;
t228 = pkin(2) * t294 + t292 * t322;
t436 = t290 * t294;
t443 = t289 * t290;
t527 = pkin(3) * t293;
t125 = (t287 * t292 - t289 * t436) * t527 - t228 * t443 + t227 * t287;
t461 = t287 * t290;
t126 = (t287 * t436 + t289 * t292) * t527 + t228 * t461 + t227 * t289;
t319 = xDP(3);
t337 = 0.1e1 / pkin(3);
t324 = xP(4);
t278 = sin(t324);
t279 = cos(t324);
t327 = koppelP(4,2);
t331 = koppelP(4,1);
t207 = t278 * t331 + t279 * t327;
t211 = -t278 * t327 + t279 * t331;
t318 = xDP(4);
t320 = xDP(2);
t321 = xDP(1);
t371 = (t211 * t318 + t320) * t266 - (-t207 * t318 + t321) * t270;
t86 = (t125 * t319 + t126 * t371) * t337 * t133;
t558 = -0.2e1 * t86;
t232 = pkin(2) * t311 + t305 * t322;
t429 = t290 * t311;
t526 = pkin(3) * t310;
t127 = (t287 * t305 - t289 * t429) * t526 - t232 * t443 + t229 * t287;
t130 = (t287 * t429 + t289 * t305) * t526 + t232 * t461 + t229 * t289;
t328 = koppelP(3,2);
t332 = koppelP(3,1);
t208 = t278 * t332 + t279 * t328;
t212 = -t278 * t328 + t279 * t332;
t370 = (t212 * t318 + t320) * t267 - (-t208 * t318 + t321) * t271;
t91 = (t127 * t319 + t130 * t370) * t337 * t134;
t557 = -0.2e1 * t91;
t233 = pkin(2) * t313 + t307 * t322;
t428 = t290 * t313;
t525 = pkin(3) * t312;
t128 = (t287 * t307 - t289 * t428) * t525 - t233 * t443 + t230 * t287;
t131 = (t287 * t428 + t289 * t307) * t525 + t233 * t461 + t230 * t289;
t329 = koppelP(2,2);
t333 = koppelP(2,1);
t209 = t278 * t333 + t279 * t329;
t213 = -t278 * t329 + t279 * t333;
t369 = (t213 * t318 + t320) * t268 - (-t209 * t318 + t321) * t272;
t92 = (t128 * t319 + t131 * t369) * t337 * t135;
t556 = -0.2e1 * t92;
t234 = pkin(2) * t315 + t309 * t322;
t427 = t290 * t315;
t524 = pkin(3) * t314;
t129 = (t287 * t309 - t289 * t427) * t524 - t234 * t443 + t231 * t287;
t132 = (t287 * t427 + t289 * t309) * t524 + t234 * t461 + t231 * t289;
t330 = koppelP(1,2);
t334 = koppelP(1,1);
t210 = t278 * t334 + t279 * t330;
t214 = -t278 * t330 + t279 * t334;
t368 = (t214 * t318 + t320) * t269 - (-t210 * t318 + t321) * t273;
t93 = (t129 * t319 + t132 * t368) * t337 * t136;
t555 = -0.2e1 * t93;
t391 = mrSges(3,1) * t293 - mrSges(3,2) * t291;
t390 = mrSges(3,1) * t310 - mrSges(3,2) * t304;
t389 = mrSges(3,1) * t312 - mrSges(3,2) * t306;
t388 = mrSges(3,1) * t314 - mrSges(3,2) * t308;
t295 = Ifges(3,1) - Ifges(3,2);
t316 = pkin(2) * mrSges(3,2);
t542 = -2 * Ifges(3,4);
t546 = Ifges(3,4) + t285 * t542 + (-t295 * t308 + t316) * t314;
t545 = Ifges(3,4) + t284 * t542 + (-t295 * t306 + t316) * t312;
t544 = Ifges(3,4) + t283 * t542 + (-t295 * t304 + t316) * t310;
t543 = Ifges(3,4) + t280 * t542 + (-t291 * t295 + t316) * t293;
t317 = mrSges(3,1) * pkin(2);
t541 = pkin(3) * t86;
t540 = pkin(3) * t91;
t539 = pkin(3) * t92;
t538 = pkin(3) * t93;
t260 = mrSges(3,2) * pkin(6) - Ifges(3,6);
t537 = -t260 / 0.2e1;
t261 = mrSges(3,1) * pkin(6) - Ifges(3,5);
t536 = t261 / 0.2e1;
t535 = pkin(2) * t291;
t534 = pkin(2) * t304;
t533 = pkin(2) * t306;
t532 = pkin(2) * t308;
t523 = g(3) * t289;
t338 = pkin(2) ^ 2;
t256 = t322 ^ 2 + t338;
t336 = pkin(3) ^ 2;
t407 = t291 * t541;
t417 = 0.2e1 * pkin(2) * pkin(3);
t437 = t290 * t292;
t189 = t287 * t437 - t289 * t294;
t457 = t288 * t293;
t154 = t189 * t291 + t287 * t457;
t190 = t287 * t294 + t289 * t437;
t442 = t289 * t293;
t155 = -t190 * t291 - t288 * t442;
t98 = (t154 * t371 + t155 * t319) * t133;
t522 = (-t322 * t407 + (t280 * t336 + t293 * t417 + t256) * t98) * t98;
t521 = mrSges(3,1) * t291;
t520 = mrSges(3,1) * t304;
t519 = mrSges(3,1) * t306;
t518 = mrSges(3,1) * t308;
t434 = t290 * t305;
t192 = t287 * t434 - t289 * t311;
t449 = t288 * t310;
t156 = t192 * t304 + t287 * t449;
t195 = t287 * t311 + t289 * t434;
t441 = t289 * t310;
t159 = -t195 * t304 - t288 * t441;
t102 = (t156 * t370 + t159 * t319) * t134;
t406 = t304 * t540;
t513 = t102 * (-t322 * t406 + (t283 * t336 + t310 * t417 + t256) * t102);
t432 = t290 * t307;
t193 = t287 * t432 - t289 * t313;
t447 = t288 * t312;
t157 = t193 * t306 + t287 * t447;
t196 = t287 * t313 + t289 * t432;
t440 = t289 * t312;
t160 = -t196 * t306 - t288 * t440;
t103 = (t157 * t369 + t160 * t319) * t135;
t405 = t306 * t539;
t512 = t103 * (-t322 * t405 + (t284 * t336 + t312 * t417 + t256) * t103);
t430 = t290 * t309;
t194 = t287 * t430 - t289 * t315;
t445 = t288 * t314;
t158 = t194 * t308 + t287 * t445;
t197 = t287 * t315 + t289 * t430;
t439 = t289 * t314;
t161 = -t197 * t308 - t288 * t439;
t104 = (t158 * t368 + t161 * t319) * t136;
t404 = t308 * t538;
t511 = t104 * (-t322 * t404 + (t285 * t336 + t314 * t417 + t256) * t104);
t510 = t294 * t98;
t509 = t322 * t98;
t274 = m(3) * pkin(2) + mrSges(2,1);
t508 = t102 * t311;
t507 = t102 * t322;
t506 = t103 * t313;
t505 = t103 * t322;
t504 = t104 * t315;
t503 = t104 * t322;
t502 = t125 * t337;
t501 = t126 * t337;
t500 = t127 * t337;
t499 = t128 * t337;
t498 = t129 * t337;
t497 = t130 * t337;
t496 = t131 * t337;
t495 = t132 * t337;
t387 = mrSges(3,2) * t293 + t521;
t153 = t290 * t391 - t387 * t458;
t494 = t153 * t337;
t493 = t154 * t266;
t492 = t154 * t270;
t491 = t156 * t267;
t490 = t156 * t271;
t489 = t157 * t268;
t488 = t157 * t272;
t487 = t158 * t269;
t486 = t158 * t273;
t386 = mrSges(3,2) * t310 + t520;
t162 = t290 * t390 - t386 * t454;
t485 = t162 * t337;
t385 = mrSges(3,2) * t312 + t519;
t163 = t290 * t389 - t385 * t452;
t484 = t163 * t337;
t384 = mrSges(3,2) * t314 + t518;
t164 = t290 * t388 - t384 * t450;
t483 = t164 * t337;
t215 = t391 + t274;
t257 = m(3) * pkin(6) - mrSges(2,2) + mrSges(3,3);
t165 = t215 * t294 + t257 * t292;
t482 = t165 * t288;
t216 = t390 + t274;
t166 = t216 * t311 + t257 * t305;
t481 = t166 * t288;
t217 = t389 + t274;
t167 = t217 * t313 + t257 * t307;
t480 = t167 * t288;
t218 = t388 + t274;
t168 = t218 * t315 + t257 * t309;
t479 = t168 * t288;
t191 = -t260 * t293 - t261 * t291;
t478 = t191 * t337;
t199 = -t260 * t310 - t261 * t304;
t477 = t199 * t337;
t200 = -t260 * t312 - t261 * t306;
t476 = t200 * t337;
t201 = -t260 * t314 - t261 * t308;
t475 = t201 * t337;
t473 = t219 * t290;
t471 = t220 * t290;
t469 = t221 * t290;
t467 = t222 * t290;
t462 = mrSges(3,2) * t259 * t288;
t460 = t288 * t289;
t459 = t288 * t291;
t456 = t288 * t294;
t455 = t288 * t304;
t453 = t288 * t306;
t451 = t288 * t308;
t448 = t288 * t311;
t446 = t288 * t313;
t444 = t288 * t315;
t426 = t290 * t337;
t424 = t292 * t293;
t420 = t305 * t310;
t419 = t307 * t312;
t418 = t309 * t314;
t412 = -0.2e1 * t316;
t403 = t291 * t509;
t402 = t304 * t507;
t401 = t306 * t505;
t400 = t308 * t503;
t399 = t266 * t501;
t398 = t270 * t501;
t397 = t267 * t497;
t396 = t271 * t497;
t395 = t268 * t496;
t394 = t272 * t496;
t393 = t269 * t495;
t392 = t273 * t495;
t383 = t223 * t287 + t523;
t381 = t224 * t287 + t523;
t379 = t225 * t287 + t523;
t377 = t226 * t287 + t523;
t41 = t403 - t541;
t13 = (((t290 * t86 + t456 * t98) * t531 + ((-t407 + t509) * t292 + pkin(2) * t510) * t457 + t41 * t290) * t98 + (t86 * t456 + (t280 * t290 - t424 * t459 - t290) * t98) * t541) * t133;
t343 = Ifges(3,1) + Ifges(2,3) + (pkin(6) ^ 2 + t338) * m(3) + 0.2e1 * pkin(6) * mrSges(3,3);
t149 = -t295 * t280 + 0.2e1 * (Ifges(3,4) * t291 + t317) * t293 + t291 * t412 + t343;
t185 = pkin(3) * t424 + t227;
t17 = t133 * t426 * t522 + (-t290 * t403 + (-t185 * t459 + (pkin(2) * t293 + t531) * t290) * t86) / (t185 * t457 + (pkin(2) + t527) * t438) * t86;
t21 = (-t293 * t522 - (pkin(2) * t86 - t293 * t41) * t541) * t133;
t5 = -t21 * t482 - t149 * t13 - t191 * t17 + ((t291 * t537 + t293 * t536) * t86 + (t317 * t291 + t543) * t98) * t558 + (-t215 * t559 - t257 * t383) * t294 + (t215 * t383 - t257 * t559) * t292;
t342 = t292 * t559 + t294 * t383;
t97 = t98 ^ 2;
t9 = -t153 * t21 - t191 * t13 - Ifges(3,3) * t17 + t97 * (pkin(2) * t521 + t543) + ((t288 * t382 - t473) * mrSges(3,1) + t342 * mrSges(3,2)) * t293 + (t462 + (-t223 * t460 + t473) * mrSges(3,2) + t342 * mrSges(3,1)) * t291;
t367 = -t154 * t5 - t501 * t9;
t43 = t402 - t540;
t14 = (((t102 * t448 + t290 * t91) * t530 + ((-t406 + t507) * t305 + pkin(2) * t508) * t449 + t43 * t290) * t102 + (t91 * t448 + (t283 * t290 - t420 * t455 - t290) * t102) * t540) * t134;
t186 = pkin(3) * t420 + t229;
t18 = t134 * t426 * t513 + (-t290 * t402 + (-t186 * t455 + (pkin(2) * t310 + t530) * t290) * t91) / (t186 * t449 + (pkin(2) + t526) * t435) * t91;
t22 = (-t310 * t513 - (pkin(2) * t91 - t310 * t43) * t540) * t134;
t341 = t305 * t560 + t311 * t381;
t99 = t102 ^ 2;
t10 = -t162 * t22 - t199 * t14 - Ifges(3,3) * t18 + t99 * (pkin(2) * t520 + t544) + ((t288 * t380 - t471) * mrSges(3,1) + t341 * mrSges(3,2)) * t310 + (t462 + (-t224 * t460 + t471) * mrSges(3,2) + t341 * mrSges(3,1)) * t304;
t150 = -t295 * t283 + 0.2e1 * (Ifges(3,4) * t304 + t317) * t310 + t304 * t412 + t343;
t6 = -t22 * t481 - t150 * t14 - t199 * t18 + ((t304 * t537 + t310 * t536) * t91 + (t317 * t304 + t544) * t102) * t557 + (-t216 * t560 - t257 * t381) * t311 + (t216 * t381 - t257 * t560) * t305;
t366 = -t10 * t497 - t156 * t6;
t100 = t103 ^ 2;
t44 = t401 - t539;
t15 = (((t103 * t446 + t290 * t92) * t529 + ((-t405 + t505) * t307 + pkin(2) * t506) * t447 + t44 * t290) * t103 + (t92 * t446 + (t284 * t290 - t419 * t453 - t290) * t103) * t539) * t135;
t187 = pkin(3) * t419 + t230;
t19 = t135 * t426 * t512 + (-t290 * t401 + (-t187 * t453 + (pkin(2) * t312 + t529) * t290) * t92) / (t187 * t447 + (pkin(2) + t525) * t433) * t92;
t23 = (-t312 * t512 - (pkin(2) * t92 - t312 * t44) * t539) * t135;
t340 = t307 * t561 + t313 * t379;
t11 = -t163 * t23 - t200 * t15 - Ifges(3,3) * t19 + t100 * (pkin(2) * t519 + t545) + ((t288 * t378 - t469) * mrSges(3,1) + t340 * mrSges(3,2)) * t312 + (t462 + (-t225 * t460 + t469) * mrSges(3,2) + t340 * mrSges(3,1)) * t306;
t151 = -t295 * t284 + 0.2e1 * (Ifges(3,4) * t306 + t317) * t312 + t306 * t412 + t343;
t7 = -t23 * t480 - t151 * t15 - t200 * t19 + ((t306 * t537 + t312 * t536) * t92 + (t317 * t306 + t545) * t103) * t556 + (-t217 * t561 - t257 * t379) * t313 + (t217 * t379 - t257 * t561) * t307;
t365 = -t11 * t496 - t157 * t7;
t101 = t104 ^ 2;
t45 = t400 - t538;
t16 = (((t104 * t444 + t290 * t93) * t528 + ((-t404 + t503) * t309 + pkin(2) * t504) * t445 + t45 * t290) * t104 + (t93 * t444 + (t285 * t290 - t418 * t451 - t290) * t104) * t538) * t136;
t188 = pkin(3) * t418 + t231;
t20 = t136 * t426 * t511 + (-t290 * t400 + (-t188 * t451 + (pkin(2) * t314 + t528) * t290) * t93) / (t188 * t445 + (pkin(2) + t524) * t431) * t93;
t24 = (-t314 * t511 - (pkin(2) * t93 - t314 * t45) * t538) * t136;
t339 = t309 * t562 + t315 * t377;
t12 = -t164 * t24 - t201 * t16 - Ifges(3,3) * t20 + t101 * (pkin(2) * t518 + t546) + ((t288 * t376 - t467) * mrSges(3,1) + t339 * mrSges(3,2)) * t314 + (t462 + (-t226 * t460 + t467) * mrSges(3,2) + t339 * mrSges(3,1)) * t308;
t152 = -t295 * t285 + 0.2e1 * (Ifges(3,4) * t308 + t317) * t314 + t308 * t412 + t343;
t8 = -t24 * t479 - t152 * t16 - t201 * t20 + ((t308 * t537 + t314 * t536) * t93 + (t317 * t308 + t546) * t104) * t555 + (-t218 * t562 - t257 * t377) * t315 + (t218 * t377 - t257 * t562) * t309;
t364 = -t12 * t495 - t158 * t8;
t363 = pkin(3) * t459 - t227 * t290;
t362 = pkin(3) * t455 - t229 * t290;
t361 = pkin(3) * t453 - t230 * t290;
t360 = pkin(3) * t451 - t231 * t290;
t359 = Ifges(3,3) * t501 + t154 * t191;
t358 = Ifges(3,3) * t497 + t156 * t199;
t357 = Ifges(3,3) * t496 + t157 * t200;
t356 = Ifges(3,3) * t495 + t158 * t201;
t355 = t126 * t478 + t149 * t154;
t354 = t130 * t477 + t150 * t156;
t353 = t131 * t476 + t151 * t157;
t352 = t132 * t475 + t152 * t158;
t351 = t133 * (t207 * t270 + t211 * t266);
t350 = t134 * (t208 * t271 + t212 * t267);
t349 = t135 * (t209 * t272 + t213 * t268);
t348 = t136 * (t210 * t273 + t214 * t269);
t347 = t126 * t494 + t154 * t482;
t346 = t130 * t485 + t156 * t481;
t345 = t131 * t484 + t157 * t480;
t344 = t132 * t483 + t158 * t479;
t326 = mrSges(4,1);
t325 = mrSges(4,2);
t303 = xDDP(1);
t302 = xDDP(2);
t301 = xDDP(3);
t300 = xDDP(4);
t286 = t318 ^ 2;
t281 = m(1) + m(2) + m(3);
t206 = -t278 * t325 + t279 * t326;
t205 = t278 * t326 + t279 * t325;
t148 = -t210 * t300 - t214 * t286 + t303;
t147 = -t209 * t300 - t213 * t286 + t303;
t146 = -t208 * t300 - t212 * t286 + t303;
t145 = -t207 * t300 - t211 * t286 + t303;
t144 = -t210 * t286 + t214 * t300 + t302;
t143 = -t209 * t286 + t213 * t300 + t302;
t142 = -t208 * t286 + t212 * t300 + t302;
t141 = -t207 * t286 + t211 * t300 + t302;
t140 = -t234 * t287 + t289 * t360;
t139 = -t233 * t287 + t289 * t361;
t138 = -t232 * t287 + t289 * t362;
t137 = -t228 * t287 + t289 * t363;
t124 = -t194 * t528 + t234 * t439 + (pkin(2) * t451 + t314 * t360) * t287;
t123 = -t193 * t529 + t233 * t440 + (pkin(2) * t453 + t312 * t361) * t287;
t122 = -t192 * t530 + t232 * t441 + (pkin(2) * t455 + t310 * t362) * t287;
t121 = -t189 * t531 + t228 * t442 + (pkin(2) * t459 + t293 * t363) * t287;
t120 = -(t197 * t269 - t273 * t450) * t528 + (t140 * t269 + t172 * t273) * t314 + (t269 * t460 + t273 * t290) * t532;
t119 = -(t196 * t268 - t272 * t452) * t529 + (t139 * t268 + t171 * t272) * t312 + (t268 * t460 + t272 * t290) * t533;
t118 = -(t195 * t267 - t271 * t454) * t530 + (t138 * t267 + t170 * t271) * t310 + (t267 * t460 + t271 * t290) * t534;
t117 = (t197 * t273 + t269 * t450) * t528 + (-t140 * t273 + t172 * t269) * t314 + (t269 * t290 - t273 * t460) * t532;
t116 = (t196 * t272 + t268 * t452) * t529 + (-t139 * t272 + t171 * t268) * t312 + (t268 * t290 - t272 * t460) * t533;
t115 = (t195 * t271 + t267 * t454) * t530 + (-t138 * t271 + t170 * t267) * t310 + (t267 * t290 - t271 * t460) * t534;
t114 = -(t190 * t266 - t270 * t458) * t531 + (t137 * t266 + t169 * t270) * t293 + (t266 * t460 + t270 * t290) * t535;
t113 = (t190 * t270 + t266 * t458) * t531 + (-t137 * t270 + t169 * t266) * t293 + (t266 * t290 - t270 * t460) * t535;
t112 = t158 * t348;
t111 = t157 * t349;
t110 = t156 * t350;
t109 = t154 * t351;
t108 = t348 * t495;
t107 = t349 * t496;
t106 = t350 * t497;
t105 = t351 * t501;
t96 = (Ifges(3,3) * t498 + t124 * t164 + t161 * t201) * t136;
t95 = (Ifges(3,3) * t499 + t123 * t163 + t160 * t200) * t135;
t94 = (Ifges(3,3) * t500 + t122 * t162 + t159 * t199) * t134;
t90 = t93 ^ 2;
t89 = t92 ^ 2;
t88 = t91 ^ 2;
t87 = (Ifges(3,3) * t502 + t121 * t153 + t155 * t191) * t133;
t85 = t86 ^ 2;
t84 = (t124 * t281 + t129 * t483 + t161 * t479) * t136;
t83 = (t123 * t281 + t128 * t484 + t160 * t480) * t135;
t82 = (t122 * t281 + t127 * t485 + t159 * t481) * t134;
t81 = (t121 * t281 + t125 * t494 + t155 * t482) * t133;
t80 = (t124 * t479 + t129 * t475 + t152 * t161) * t136;
t79 = (t123 * t480 + t128 * t476 + t151 * t160) * t135;
t78 = (t122 * t481 + t127 * t477 + t150 * t159) * t134;
t77 = (t121 * t482 + t125 * t478 + t149 * t155) * t133;
t76 = (-t117 * t210 + t120 * t214) * t136;
t75 = (-t116 * t209 + t119 * t213) * t135;
t74 = (-t115 * t208 + t118 * t212) * t134;
t73 = (-t113 * t207 + t114 * t211) * t133;
t72 = (t120 * t164 + t269 * t356) * t136;
t71 = (t119 * t163 + t268 * t357) * t135;
t70 = (t118 * t162 + t267 * t358) * t134;
t69 = (t117 * t164 - t273 * t356) * t136;
t68 = (t116 * t163 - t272 * t357) * t135;
t67 = (t115 * t162 - t271 * t358) * t134;
t66 = (t114 * t153 + t266 * t359) * t133;
t65 = (t113 * t153 - t270 * t359) * t133;
t64 = (t120 * t281 + t269 * t344) * t136;
t63 = (t119 * t281 + t268 * t345) * t135;
t62 = (t118 * t281 + t267 * t346) * t134;
t61 = (t117 * t281 - t273 * t344) * t136;
t60 = (t116 * t281 - t272 * t345) * t135;
t59 = (t115 * t281 - t271 * t346) * t134;
t58 = (t114 * t281 + t266 * t347) * t133;
t57 = (t113 * t281 - t270 * t347) * t133;
t56 = (t120 * t479 + t269 * t352) * t136;
t55 = (t119 * t480 + t268 * t353) * t135;
t54 = (t118 * t481 + t267 * t354) * t134;
t53 = (t117 * t479 - t273 * t352) * t136;
t52 = (t116 * t480 - t272 * t353) * t135;
t51 = (t115 * t481 - t271 * t354) * t134;
t50 = (t114 * t482 + t266 * t355) * t133;
t49 = (t113 * t482 - t270 * t355) * t133;
t40 = Ifges(3,3) * t108 + t112 * t201 + t164 * t76;
t39 = Ifges(3,3) * t107 + t111 * t200 + t163 * t75;
t38 = Ifges(3,3) * t106 + t110 * t199 + t162 * t74;
t37 = t108 * t164 + t112 * t479 + t281 * t76;
t36 = t107 * t163 + t111 * t480 + t281 * t75;
t35 = t106 * t162 + t110 * t481 + t281 * t74;
t34 = Ifges(3,3) * t105 + t109 * t191 + t153 * t73;
t33 = t105 * t153 + t109 * t482 + t281 * t73;
t32 = t108 * t201 + t112 * t152 + t479 * t76;
t31 = t107 * t200 + t111 * t151 + t480 * t75;
t30 = t106 * t199 + t110 * t150 + t481 * t74;
t29 = t105 * t191 + t109 * t149 + t482 * t73;
t4 = -t164 * t20 - t90 * t290 * t384 + (-t168 * t16 + (-t101 * t274 - t388 * (t101 + t90)) * t309 + (t104 * t257 + t384 * t555) * t504) * t288 + (-t24 - t222) * t281;
t3 = -t163 * t19 - t89 * t290 * t385 + (-t167 * t15 + (-t100 * t274 - t389 * (t100 + t89)) * t307 + (t103 * t257 + t385 * t556) * t506) * t288 + (-t23 - t221) * t281;
t2 = -t162 * t18 - t88 * t290 * t386 + (-t166 * t14 + (-t274 * t99 - t390 * (t99 + t88)) * t305 + (t102 * t257 + t386 * t557) * t508) * t288 + (-t22 - t220) * t281;
t1 = -t153 * t17 - t85 * t290 * t387 + (-t165 * t13 + (-t274 * t97 - t391 * (t97 + t85)) * t292 + (t257 * t98 + t387 * t558) * t510) * t288 + (-t21 - t219) * t281;
t25 = [-t205 * t300 - t286 * t206 + (t303 - g(1)) * m(4) + ((t124 * t61 + t161 * t53 + t498 * t69) * t301 + (t120 * t61 + t393 * t69 + t487 * t53) * t144 + (t148 * t61 + t4) * t117 + ((-t158 * t53 - t495 * t69) * t148 + t364) * t273) * t136 + ((t123 * t60 + t160 * t52 + t499 * t68) * t301 + (t119 * t60 + t395 * t68 + t489 * t52) * t143 + (t147 * t60 + t3) * t116 + ((-t157 * t52 - t496 * t68) * t147 + t365) * t272) * t135 + ((t122 * t59 + t159 * t51 + t500 * t67) * t301 + (t118 * t59 + t397 * t67 + t491 * t51) * t142 + (t146 * t59 + t2) * t115 + ((-t156 * t51 - t497 * t67) * t146 + t366) * t271) * t134 + ((t121 * t57 + t155 * t49 + t502 * t65) * t301 + (t114 * t57 + t399 * t65 + t49 * t493) * t141 + (t145 * t57 + t1) * t113 + ((-t154 * t49 - t501 * t65) * t145 + t367) * t270) * t133; -t286 * t205 + t206 * t300 + (t302 - g(2)) * m(4) + ((t124 * t64 + t161 * t56 + t498 * t72) * t301 + (t117 * t64 - t392 * t72 - t486 * t56) * t148 + (t144 * t64 + t4) * t120 + ((t158 * t56 + t495 * t72) * t144 - t364) * t269) * t136 + ((t123 * t63 + t160 * t55 + t499 * t71) * t301 + (t116 * t63 - t394 * t71 - t488 * t55) * t147 + (t143 * t63 + t3) * t119 + ((t157 * t55 + t496 * t71) * t143 - t365) * t268) * t135 + ((t122 * t62 + t159 * t54 + t500 * t70) * t301 + (t115 * t62 - t396 * t70 - t490 * t54) * t146 + (t142 * t62 + t2) * t118 + ((t156 * t54 + t497 * t70) * t142 - t366) * t267) * t134 + ((t121 * t58 + t155 * t50 + t502 * t66) * t301 + (t113 * t58 - t398 * t66 - t492 * t50) * t145 + (t141 * t58 + t1) * t114 + ((t154 * t50 + t501 * t66) * t141 - t367) * t266) * t133; (t301 - g(3)) * m(4) + (t124 * t4 + t161 * t8 + (t117 * t84 - t392 * t96 - t486 * t80) * t148 + (t120 * t84 + t393 * t96 + t487 * t80) * t144 + (t124 * t84 + t161 * t80 + t498 * t96) * t301 + t12 * t498) * t136 + (t123 * t3 + t160 * t7 + (t116 * t83 - t394 * t95 - t488 * t79) * t147 + (t119 * t83 + t395 * t95 + t489 * t79) * t143 + (t123 * t83 + t160 * t79 + t499 * t95) * t301 + t11 * t499) * t135 + (t122 * t2 + t159 * t6 + (t115 * t82 - t396 * t94 - t490 * t78) * t146 + (t118 * t82 + t397 * t94 + t491 * t78) * t142 + (t122 * t82 + t159 * t78 + t500 * t94) * t301 + t10 * t500) * t134 + (t121 * t1 + t155 * t5 + (t121 * t81 + t155 * t77 + t502 * t87) * t301 + (t113 * t81 - t398 * t87 - t492 * t77) * t145 + (t114 * t81 + t399 * t87 + t493 * t77) * t141 + t9 * t502) * t133; -(-g(1) * t326 - g(2) * t325) * t278 + t279 * (g(1) * t325 - g(2) * t326) + Ifges(4,3) * t300 + t206 * t302 - t205 * t303 + t110 * t6 + t111 * t7 + t112 * t8 + t106 * t10 + t107 * t11 + t108 * t12 + t109 * t5 + t105 * t9 + t75 * t3 + t76 * t4 + t73 * t1 + t74 * t2 + ((t124 * t37 + t161 * t32 + t40 * t498) * t301 + (t117 * t37 - t32 * t486 - t392 * t40) * t148 + (t120 * t37 + t32 * t487 + t393 * t40) * t144) * t136 + ((t123 * t36 + t160 * t31 + t39 * t499) * t301 + (t116 * t36 - t31 * t488 - t39 * t394) * t147 + (t119 * t36 + t31 * t489 + t39 * t395) * t143) * t135 + ((t122 * t35 + t159 * t30 + t38 * t500) * t301 + (t115 * t35 - t30 * t490 - t38 * t396) * t146 + (t118 * t35 + t30 * t491 + t38 * t397) * t142) * t134 + ((t121 * t33 + t155 * t29 + t34 * t502) * t301 + (t113 * t33 - t29 * t492 - t34 * t398) * t145 + (t114 * t33 + t29 * t493 + t34 * t399) * t141) * t133;];
tauX  = t25;
