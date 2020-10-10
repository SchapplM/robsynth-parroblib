% Calculate vector of inverse dynamics forces for parallel robot
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
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
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
% Datum: 2020-08-07 11:06
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P4PRRRR8V1G2A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V1G2A0_invdyn_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR8V1G2A0_invdyn_para_pf_slag_vp2: xDP has to be [4x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [4 1]), ...
  'P4PRRRR8V1G2A0_invdyn_para_pf_slag_vp2: xDDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V1G2A0_invdyn_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR8V1G2A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P4PRRRR8V1G2A0_invdyn_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V1G2A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRRR8V1G2A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P4PRRRR8V1G2A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V1G2A0_invdyn_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V1G2A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:02:49
% EndTime: 2020-08-07 11:03:01
% DurationCPUTime: 11.96s
% Computational Cost: add. (39875->623), mult. (112801->1170), div. (13536->9), fcn. (119718->30), ass. (0->468)
t285 = cos(qJ(3,1));
t286 = cos(qJ(2,1));
t280 = sin(qJ(2,1));
t405 = t280 * t285;
t215 = pkin(2) * t405 - pkin(5) * t286;
t258 = sin(pkin(3));
t279 = sin(qJ(3,1));
t260 = cos(pkin(3));
t515 = pkin(2) * t260;
t168 = t215 * t258 + t279 * t515;
t532 = 0.1e1 / t168;
t461 = t532 / t285;
t283 = cos(qJ(3,2));
t284 = cos(qJ(2,2));
t278 = sin(qJ(2,2));
t408 = t278 * t283;
t214 = pkin(2) * t408 - pkin(5) * t284;
t277 = sin(qJ(3,2));
t167 = t214 * t258 + t277 * t515;
t533 = 0.1e1 / t167;
t462 = t533 / t283;
t281 = cos(qJ(3,3));
t282 = cos(qJ(2,3));
t276 = sin(qJ(2,3));
t411 = t276 * t281;
t213 = pkin(2) * t411 - pkin(5) * t282;
t275 = sin(qJ(3,3));
t166 = t213 * t258 + t275 * t515;
t534 = 0.1e1 / t166;
t463 = t534 / t281;
t263 = cos(qJ(3,4));
t264 = cos(qJ(2,4));
t262 = sin(qJ(2,4));
t414 = t262 * t263;
t203 = pkin(2) * t414 - pkin(5) * t264;
t261 = sin(qJ(3,4));
t162 = t203 * t258 + t261 * t515;
t535 = 0.1e1 / t162;
t464 = t535 / t263;
t380 = t263 * mrSges(3,1) - mrSges(3,2) * t261;
t379 = t281 * mrSges(3,1) - mrSges(3,2) * t275;
t378 = t283 * mrSges(3,1) - mrSges(3,2) * t277;
t377 = t285 * mrSges(3,1) - mrSges(3,2) * t279;
t270 = legFrame(1,2);
t240 = sin(t270);
t244 = cos(t270);
t212 = g(1) * t244 - g(2) * t240;
t259 = cos(pkin(6));
t257 = sin(pkin(6));
t510 = g(3) * t257;
t173 = t212 * t259 - t510;
t208 = g(1) * t240 + g(2) * t244;
t361 = t173 * t260 + t208 * t258;
t269 = legFrame(2,2);
t239 = sin(t269);
t243 = cos(t269);
t211 = g(1) * t243 - g(2) * t239;
t172 = t211 * t259 - t510;
t207 = g(1) * t239 + g(2) * t243;
t362 = t172 * t260 + t207 * t258;
t268 = legFrame(3,2);
t238 = sin(t268);
t242 = cos(t268);
t210 = g(1) * t242 - g(2) * t238;
t171 = t210 * t259 - t510;
t206 = g(1) * t238 + g(2) * t242;
t363 = t171 * t260 + t206 * t258;
t267 = legFrame(4,2);
t237 = sin(t267);
t241 = cos(t267);
t209 = g(1) * t241 - g(2) * t237;
t170 = t209 * t259 - t510;
t205 = g(1) * t237 + g(2) * t241;
t364 = t170 * t260 + t205 * t258;
t291 = xP(4);
t245 = sin(t291);
t246 = cos(t291);
t297 = koppelP(1,2);
t301 = koppelP(1,1);
t198 = t245 * t301 + t246 * t297;
t202 = -t245 * t297 + t246 * t301;
t287 = xDP(4);
t256 = t287 ^ 2;
t271 = xDDP(4);
t274 = xDDP(1);
t144 = -t198 * t271 - t202 * t256 + t274;
t473 = t144 * t244;
t273 = xDDP(2);
t140 = -t198 * t256 + t202 * t271 + t273;
t477 = t140 * t240;
t540 = -t473 + t477;
t296 = koppelP(2,2);
t300 = koppelP(2,1);
t197 = t245 * t300 + t246 * t296;
t201 = -t245 * t296 + t246 * t300;
t143 = -t197 * t271 - t201 * t256 + t274;
t474 = t143 * t243;
t139 = -t197 * t256 + t201 * t271 + t273;
t478 = t139 * t239;
t539 = -t474 + t478;
t295 = koppelP(3,2);
t299 = koppelP(3,1);
t196 = t245 * t299 + t246 * t295;
t200 = -t245 * t295 + t246 * t299;
t142 = -t196 * t271 - t200 * t256 + t274;
t475 = t142 * t242;
t138 = -t196 * t256 + t200 * t271 + t273;
t479 = t138 * t238;
t538 = -t475 + t479;
t294 = koppelP(4,2);
t298 = koppelP(4,1);
t195 = t245 * t298 + t246 * t294;
t199 = -t245 * t294 + t246 * t298;
t141 = -t195 * t271 - t199 * t256 + t274;
t476 = t141 * t241;
t137 = -t195 * t256 + t199 * t271 + t273;
t480 = t137 * t237;
t537 = -t476 + t480;
t536 = 2 * Ifges(3,4);
t247 = t263 ^ 2;
t531 = 0.2e1 * t247;
t250 = t281 ^ 2;
t530 = 0.2e1 * t250;
t252 = t283 ^ 2;
t529 = 0.2e1 * t252;
t254 = t285 ^ 2;
t528 = 0.2e1 * t254;
t527 = Ifges(3,5) / 0.2e1;
t526 = -Ifges(3,6) / 0.2e1;
t424 = t260 * t264;
t425 = t260 * t262;
t514 = pkin(2) * t263;
t129 = -(-t257 * t262 + t259 * t424) * t514 - pkin(5) * (t257 * t264 + t259 * t425);
t130 = (t257 * t424 + t259 * t262) * t514 + (t257 * t425 - t259 * t264) * pkin(5);
t288 = xDP(3);
t304 = 0.1e1 / pkin(2);
t289 = xDP(2);
t290 = xDP(1);
t360 = (t199 * t287 + t289) * t237 - (-t195 * t287 + t290) * t241;
t74 = (t129 * t288 + t360 * t130) * t304 * t464;
t525 = pkin(2) * t74;
t420 = t260 * t282;
t423 = t260 * t276;
t513 = pkin(2) * t281;
t131 = -(-t257 * t276 + t259 * t420) * t513 - pkin(5) * (t257 * t282 + t259 * t423);
t134 = (t257 * t420 + t259 * t276) * t513 + (t257 * t423 - t259 * t282) * pkin(5);
t359 = (t200 * t287 + t289) * t238 - (-t196 * t287 + t290) * t242;
t79 = (t131 * t288 + t359 * t134) * t304 * t463;
t524 = pkin(2) * t79;
t419 = t260 * t284;
t422 = t260 * t278;
t512 = pkin(2) * t283;
t132 = -(-t257 * t278 + t259 * t419) * t512 - pkin(5) * (t257 * t284 + t259 * t422);
t135 = (t257 * t419 + t259 * t278) * t512 + (t257 * t422 - t259 * t284) * pkin(5);
t358 = (t201 * t287 + t289) * t239 - (-t197 * t287 + t290) * t243;
t80 = (t132 * t288 + t358 * t135) * t304 * t462;
t523 = pkin(2) * t80;
t418 = t260 * t286;
t421 = t260 * t280;
t511 = pkin(2) * t285;
t133 = -(-t257 * t280 + t259 * t418) * t511 - pkin(5) * (t257 * t286 + t259 * t421);
t136 = (t257 * t418 + t259 * t280) * t511 + (t257 * t421 - t259 * t286) * pkin(5);
t357 = (t202 * t287 + t289) * t240 - (-t198 * t287 + t290) * t244;
t81 = (t133 * t288 + t357 * t136) * t304 * t461;
t522 = pkin(2) * t81;
t324 = t258 * t263 + t261 * t425;
t416 = t261 * t264;
t145 = t324 * t257 - t259 * t416;
t146 = -t257 * t416 - t324 * t259;
t95 = (t360 * t145 + t146 * t288) * t464;
t521 = pkin(5) * t95;
t266 = mrSges(2,2) - mrSges(3,3);
t520 = t266 / 0.2e1;
t519 = pkin(2) * t247;
t518 = pkin(2) * t250;
t517 = pkin(2) * t252;
t516 = pkin(2) * t254;
t509 = g(3) * t259;
t302 = pkin(5) ^ 2;
t303 = pkin(2) ^ 2;
t499 = t261 * t74;
t404 = pkin(2) * t499;
t508 = (-pkin(5) * t404 + (t247 * t303 + t302) * t95) * t95;
t507 = Ifges(3,1) + Ifges(2,3);
t323 = t258 * t281 + t275 * t423;
t412 = t275 * t282;
t147 = t323 * t257 - t259 * t412;
t150 = -t257 * t412 - t323 * t259;
t102 = (t359 * t147 + t150 * t288) * t463;
t498 = t275 * t79;
t403 = pkin(2) * t498;
t502 = t102 * (-pkin(5) * t403 + (t250 * t303 + t302) * t102);
t322 = t258 * t283 + t277 * t422;
t409 = t277 * t284;
t148 = t322 * t257 - t259 * t409;
t151 = -t257 * t409 - t322 * t259;
t103 = (t358 * t148 + t151 * t288) * t462;
t497 = t277 * t80;
t402 = pkin(2) * t497;
t501 = t103 * (-pkin(5) * t402 + (t252 * t303 + t302) * t103);
t321 = t258 * t285 + t279 * t421;
t406 = t279 * t286;
t149 = t321 * t257 - t259 * t406;
t152 = -t257 * t406 - t321 * t259;
t104 = (t357 * t149 + t152 * t288) * t461;
t496 = t279 * t81;
t401 = pkin(2) * t496;
t500 = t104 * (-pkin(5) * t401 + (t254 * t303 + t302) * t104);
t495 = t102 * t275;
t494 = t103 * t277;
t493 = t104 * t279;
t204 = pkin(5) * t262 + t264 * t514;
t433 = t258 * t261;
t352 = pkin(2) * t433 - t203 * t260;
t122 = t204 * t259 + t352 * t257;
t492 = t122 * t535;
t216 = pkin(5) * t276 + t282 * t513;
t431 = t258 * t275;
t351 = pkin(2) * t431 - t213 * t260;
t126 = t216 * t259 + t351 * t257;
t491 = t126 * t534;
t217 = pkin(5) * t278 + t284 * t512;
t430 = t258 * t277;
t350 = pkin(2) * t430 - t214 * t260;
t127 = t217 * t259 + t350 * t257;
t490 = t127 * t533;
t218 = pkin(5) * t280 + t286 * t511;
t429 = t258 * t279;
t349 = pkin(2) * t429 - t215 * t260;
t128 = t218 * t259 + t349 * t257;
t489 = t128 * t532;
t488 = t129 * t304;
t487 = t130 * t304;
t486 = t131 * t304;
t485 = t132 * t304;
t484 = t133 * t304;
t483 = t134 * t304;
t482 = t135 * t304;
t481 = t136 * t304;
t376 = mrSges(3,1) * t261 + mrSges(3,2) * t263;
t153 = -t258 * t376 * t262 + t380 * t260;
t472 = t153 * t535;
t471 = t153 * t304;
t375 = mrSges(3,1) * t275 + mrSges(3,2) * t281;
t154 = -t258 * t375 * t276 + t379 * t260;
t470 = t154 * t534;
t469 = t154 * t304;
t374 = mrSges(3,1) * t277 + mrSges(3,2) * t283;
t155 = -t258 * t374 * t278 + t378 * t260;
t468 = t155 * t533;
t467 = t155 * t304;
t373 = mrSges(3,1) * t279 + mrSges(3,2) * t285;
t156 = -t258 * t373 * t280 + t377 * t260;
t466 = t156 * t532;
t465 = t156 * t304;
t249 = m(1) + m(2) + m(3);
t460 = t535 * t249;
t459 = t534 * t249;
t458 = t533 * t249;
t457 = t532 * t249;
t219 = mrSges(2,1) + t380;
t456 = (t219 * t264 - t262 * t266) * t258;
t220 = mrSges(2,1) + t379;
t451 = (t220 * t282 - t266 * t276) * t258;
t221 = mrSges(2,1) + t378;
t450 = (t221 * t284 - t266 * t278) * t258;
t222 = mrSges(2,1) + t377;
t449 = (t222 * t286 - t266 * t280) * t258;
t447 = t205 * t260;
t445 = t206 * t260;
t443 = t207 * t260;
t441 = t208 * t260;
t223 = Ifges(3,5) * t261 + Ifges(3,6) * t263;
t440 = t223 * t304;
t224 = Ifges(3,5) * t275 + Ifges(3,6) * t281;
t439 = t224 * t304;
t225 = Ifges(3,5) * t277 + Ifges(3,6) * t283;
t438 = t225 * t304;
t226 = Ifges(3,5) * t279 + Ifges(3,6) * t285;
t437 = t226 * t304;
t436 = mrSges(3,2) * t510 * t258;
t435 = t257 * t266;
t434 = t258 * t259;
t432 = t258 * t264;
t428 = t258 * t282;
t427 = t258 * t284;
t426 = t258 * t286;
t417 = t260 * t304;
t265 = Ifges(3,1) - Ifges(3,2);
t415 = t261 * t265;
t413 = t275 * t281;
t410 = t277 * t283;
t407 = t279 * t285;
t400 = t261 * t521;
t399 = pkin(5) * t495;
t398 = pkin(5) * t494;
t397 = pkin(5) * t493;
t396 = t237 * t464;
t395 = t241 * t464;
t394 = t238 * t463;
t393 = t242 * t463;
t392 = t239 * t462;
t391 = t243 * t462;
t390 = t240 * t461;
t389 = t244 * t461;
t388 = t535 * t456;
t387 = t534 * t451;
t386 = t533 * t450;
t385 = t532 * t449;
t384 = t258 * t414;
t383 = t258 * t411;
t382 = t258 * t408;
t381 = t258 * t405;
t372 = t209 * t257 + t509;
t371 = t210 * t257 + t509;
t370 = t211 * t257 + t509;
t369 = t212 * t257 + t509;
t121 = -t204 * t257 + t352 * t259;
t113 = t121 * t237 + t162 * t241;
t114 = -t121 * t241 + t162 * t237;
t368 = t113 * t137 + t114 * t141;
t123 = -t216 * t257 + t351 * t259;
t115 = t123 * t238 + t166 * t242;
t116 = -t123 * t242 + t166 * t238;
t367 = t115 * t138 + t116 * t142;
t124 = -t217 * t257 + t350 * t259;
t117 = t124 * t239 + t167 * t243;
t118 = -t124 * t243 + t167 * t239;
t366 = t117 * t139 + t118 * t143;
t125 = -t218 * t257 + t349 * t259;
t119 = t125 * t240 + t168 * t244;
t120 = -t125 * t244 + t168 * t240;
t365 = t119 * t140 + t120 * t144;
t41 = t400 - t525;
t13 = (((t260 * t74 + t95 * t432) * t519 - (t404 - t521) * t384 + t260 * t41) * t95 - (-t74 * t432 + (-t247 * t260 + t261 * t384 + t260) * t95) * t525) * t464;
t17 = (t417 * t508 + (-t74 * t203 * t433 + t260 * (t74 * t519 - t400)) * t74) * t464;
t185 = t261 * t263 * t536 - t247 * t265 + t507;
t21 = (t41 * t525 - t508) * t535;
t227 = t266 * t509;
t5 = -t21 * t456 - t185 * t13 - t223 * t17 + 0.2e1 * t74 * ((t95 * t415 + t74 * t527) * t263 + t499 * t526 + (t531 - 0.1e1) * t95 * Ifges(3,4)) + (t209 * t435 - t364 * t219 + t227) * t264 + (t372 * t219 + t364 * t266) * t262;
t308 = t364 * t262 + t372 * t264;
t94 = t95 ^ 2;
t9 = -t153 * t21 - t223 * t13 - Ifges(3,3) * t17 - t94 * (Ifges(3,4) * t531 + t263 * t415 - Ifges(3,4)) + ((t170 * t258 - t447) * mrSges(3,1) + t308 * mrSges(3,2)) * t263 + (t436 + (-t209 * t434 + t447) * mrSges(3,2) + t308 * mrSges(3,1)) * t261;
t356 = t145 * t5 + t9 * t487;
t43 = t399 - t524;
t14 = (((t102 * t428 + t260 * t79) * t518 - (-pkin(5) * t102 + t403) * t383 + t260 * t43) * t102 - (-t79 * t428 + (-t250 * t260 + t275 * t383 + t260) * t102) * t524) * t463;
t18 = (t417 * t502 + (-t79 * t213 * t431 + t260 * (t79 * t518 - t399)) * t79) * t463;
t22 = (t43 * t524 - t502) * t534;
t307 = t363 * t276 + t371 * t282;
t99 = t102 ^ 2;
t10 = -t154 * t22 - t224 * t14 - Ifges(3,3) * t18 - t99 * (Ifges(3,4) * t530 + t265 * t413 - Ifges(3,4)) + ((t171 * t258 - t445) * mrSges(3,1) + t307 * mrSges(3,2)) * t281 + (t436 + (-t210 * t434 + t445) * mrSges(3,2) + t307 * mrSges(3,1)) * t275;
t187 = -t250 * t265 + t413 * t536 + t507;
t6 = -t22 * t451 - t187 * t14 - t224 * t18 + 0.2e1 * t79 * ((t265 * t495 + t79 * t527) * t281 + t498 * t526 + (t530 - 0.1e1) * t102 * Ifges(3,4)) + (t210 * t435 - t363 * t220 + t227) * t282 + (t371 * t220 + t363 * t266) * t276;
t355 = t10 * t483 + t147 * t6;
t100 = t103 ^ 2;
t44 = t398 - t523;
t15 = (((t103 * t427 + t260 * t80) * t517 - (-pkin(5) * t103 + t402) * t382 + t260 * t44) * t103 + (t80 * t427 + (t252 * t260 - t277 * t382 - t260) * t103) * t523) * t462;
t19 = (t417 * t501 + (-t80 * t214 * t430 + t260 * (t80 * t517 - t398)) * t80) * t462;
t23 = (t44 * t523 - t501) * t533;
t306 = t362 * t278 + t370 * t284;
t11 = -t155 * t23 - t225 * t15 - Ifges(3,3) * t19 - t100 * (Ifges(3,4) * t529 + t265 * t410 - Ifges(3,4)) + ((t172 * t258 - t443) * mrSges(3,1) + t306 * mrSges(3,2)) * t283 + (t436 + (-t211 * t434 + t443) * mrSges(3,2) + t306 * mrSges(3,1)) * t277;
t188 = -t252 * t265 + t410 * t536 + t507;
t7 = -t23 * t450 - t188 * t15 - t225 * t19 + 0.2e1 * t80 * ((t265 * t494 + t80 * t527) * t283 + t497 * t526 + (t529 - 0.1e1) * t103 * Ifges(3,4)) + (t211 * t435 - t362 * t221 + t227) * t284 + (t370 * t221 + t362 * t266) * t278;
t354 = t11 * t482 + t148 * t7;
t101 = t104 ^ 2;
t45 = t397 - t522;
t16 = (((t104 * t426 + t260 * t81) * t516 - (-pkin(5) * t104 + t401) * t381 + t260 * t45) * t104 - (-t81 * t426 + (-t254 * t260 + t279 * t381 + t260) * t104) * t522) * t461;
t20 = (t417 * t500 + (-t81 * t215 * t429 + t260 * (t81 * t516 - t397)) * t81) * t461;
t24 = (t45 * t522 - t500) * t532;
t305 = t361 * t280 + t369 * t286;
t12 = -t156 * t24 - t226 * t16 - Ifges(3,3) * t20 - t101 * (Ifges(3,4) * t528 + t265 * t407 - Ifges(3,4)) + ((t173 * t258 - t441) * mrSges(3,1) + t305 * mrSges(3,2)) * t285 + (t436 + (-t212 * t434 + t441) * mrSges(3,2) + t305 * mrSges(3,1)) * t279;
t189 = -t254 * t265 + t407 * t536 + t507;
t8 = -t24 * t449 - t189 * t16 - t226 * t20 + 0.2e1 * t81 * ((t265 * t493 + t81 * t527) * t285 + t496 * t526 + (t528 - 0.1e1) * t104 * Ifges(3,4)) + (t212 * t435 - t361 * t222 + t227) * t286 + (t369 * t222 + t361 * t266) * t280;
t353 = t12 * t481 + t149 * t8;
t328 = t130 * t440 + t145 * t185;
t51 = t113 * t388 + t328 * t396;
t332 = Ifges(3,3) * t487 + t145 * t223;
t65 = t113 * t472 + t332 * t396;
t347 = t145 * t51 + t65 * t487;
t52 = t114 * t388 - t328 * t395;
t66 = t114 * t472 - t332 * t395;
t346 = t145 * t52 + t66 * t487;
t327 = t134 * t439 + t147 * t187;
t59 = t115 * t387 + t327 * t394;
t331 = Ifges(3,3) * t483 + t147 * t224;
t67 = t115 * t470 + t331 * t394;
t343 = t147 * t59 + t67 * t483;
t60 = t116 * t387 - t327 * t393;
t68 = t116 * t470 - t331 * t393;
t342 = t147 * t60 + t68 * t483;
t326 = t135 * t438 + t148 * t188;
t61 = t117 * t386 + t326 * t392;
t330 = Ifges(3,3) * t482 + t148 * t225;
t69 = t117 * t468 + t330 * t392;
t339 = t148 * t61 + t69 * t482;
t62 = t118 * t386 - t326 * t391;
t70 = t118 * t468 - t330 * t391;
t338 = t148 * t62 + t70 * t482;
t325 = t136 * t437 + t149 * t189;
t63 = t119 * t385 + t325 * t390;
t329 = Ifges(3,3) * t481 + t149 * t226;
t71 = t119 * t466 + t329 * t390;
t335 = t149 * t63 + t71 * t481;
t64 = t120 * t385 - t325 * t389;
t72 = t120 * t466 - t329 * t389;
t334 = t149 * t64 + t72 * t481;
t320 = t130 * t471 + t145 * t456;
t319 = t134 * t469 + t147 * t451;
t318 = t135 * t467 + t148 * t450;
t317 = t136 * t465 + t149 * t449;
t316 = (t195 * t241 + t199 * t237) * t464;
t315 = (t196 * t242 + t200 * t238) * t463;
t314 = (t197 * t243 + t201 * t239) * t462;
t313 = (t198 * t244 + t202 * t240) * t461;
t272 = xDDP(3);
t312 = t122 * t272 + t368;
t311 = t126 * t272 + t367;
t310 = t127 * t272 + t366;
t309 = t128 * t272 + t365;
t293 = mrSges(4,1);
t292 = mrSges(4,2);
t194 = -t292 * t245 + t246 * t293;
t193 = t293 * t245 + t246 * t292;
t112 = t149 * t313;
t111 = t148 * t314;
t110 = t147 * t315;
t109 = t145 * t316;
t108 = t313 * t481;
t107 = t314 * t482;
t106 = t315 * t483;
t105 = t316 * t487;
t98 = (t119 * t202 - t120 * t198) * t532;
t97 = (t117 * t201 - t118 * t197) * t533;
t96 = (t115 * t200 - t116 * t196) * t534;
t93 = t128 * t466 + (Ifges(3,3) * t484 + t152 * t226) * t461;
t92 = t127 * t468 + (Ifges(3,3) * t485 + t151 * t225) * t462;
t91 = t126 * t470 + (Ifges(3,3) * t486 + t150 * t224) * t463;
t90 = (t113 * t199 - t114 * t195) * t535;
t89 = t122 * t472 + (Ifges(3,3) * t488 + t146 * t223) * t464;
t88 = t128 * t385 + (t133 * t437 + t152 * t189) * t461;
t87 = t127 * t386 + (t132 * t438 + t151 * t188) * t462;
t86 = t126 * t387 + (t131 * t439 + t150 * t187) * t463;
t85 = t122 * t388 + (t129 * t440 + t146 * t185) * t464;
t84 = t128 * t457 + (t133 * t465 + t152 * t449) * t461;
t83 = t127 * t458 + (t132 * t467 + t151 * t450) * t462;
t82 = t126 * t459 + (t131 * t469 + t150 * t451) * t463;
t78 = t81 ^ 2;
t77 = t80 ^ 2;
t76 = t79 ^ 2;
t75 = t122 * t460 + (t129 * t471 + t146 * t456) * t464;
t73 = t74 ^ 2;
t40 = Ifges(3,3) * t108 + t112 * t226 + t156 * t98;
t39 = Ifges(3,3) * t107 + t111 * t225 + t155 * t97;
t38 = Ifges(3,3) * t106 + t110 * t224 + t154 * t96;
t37 = t108 * t226 + t112 * t189 + t98 * t449;
t36 = t107 * t225 + t111 * t188 + t97 * t450;
t35 = t106 * t224 + t110 * t187 + t96 * t451;
t34 = Ifges(3,3) * t105 + t109 * t223 + t153 * t90;
t33 = t108 * t156 + t112 * t449 + t249 * t98;
t32 = t107 * t155 + t111 * t450 + t249 * t97;
t31 = t106 * t154 + t110 * t451 + t249 * t96;
t30 = t105 * t223 + t109 * t185 + t90 * t456;
t29 = t105 * t153 + t109 * t456 + t249 * t90;
t4 = -t16 * t449 - t156 * t20 + ((-mrSges(2,1) * t101 - t377 * (t101 + t78)) * t280 - 0.2e1 * t286 * (t104 * t520 + t373 * t81) * t104) * t258 - t260 * t78 * t373 + (-t24 - t208) * t249;
t3 = -t15 * t450 - t155 * t19 + ((-mrSges(2,1) * t100 - t378 * (t100 + t77)) * t278 - 0.2e1 * t284 * (t103 * t520 + t374 * t80) * t103) * t258 - t260 * t77 * t374 + (-t23 - t207) * t249;
t2 = -t14 * t451 - t154 * t18 + ((-mrSges(2,1) * t99 - t379 * (t99 + t76)) * t276 - 0.2e1 * t282 * (t102 * t520 + t375 * t79) * t102) * t258 - t260 * t76 * t375 + (-t22 - t206) * t249;
t1 = -t13 * t456 - t153 * t17 + ((-mrSges(2,1) * t94 - t380 * (t94 + t73)) * t262 - 0.2e1 * t264 * (t376 * t74 + t95 * t520) * t95) * t258 - t260 * t73 * t376 + (-t21 - t205) * t249;
t25 = [-t193 * t271 - t256 * t194 + (t274 - g(1)) * m(4) + (t120 * t4 + t309 * (t120 * t457 - t317 * t389)) * t532 + (t118 * t3 + t310 * (t118 * t458 - t318 * t391)) * t533 + (t116 * t2 + t311 * (t116 * t459 - t319 * t393)) * t534 + (t114 * t1 + t312 * (t114 * t460 - t320 * t395)) * t535 + ((t152 * t64 + t72 * t484) * t272 + t334 * t477 + (-t334 * t144 - t353) * t244) * t461 + ((t151 * t62 + t70 * t485) * t272 + t338 * t478 + (-t338 * t143 - t354) * t243) * t462 + ((t150 * t60 + t68 * t486) * t272 + t342 * t479 + (-t342 * t142 - t355) * t242) * t463 + ((t146 * t52 + t66 * t488) * t272 + t346 * t480 + (-t346 * t141 - t356) * t241) * t464; -t256 * t193 + t194 * t271 + (t273 - g(2)) * m(4) + (t119 * t4 + t309 * (t119 * t457 + t317 * t390)) * t532 + (t117 * t3 + t310 * (t117 * t458 + t318 * t392)) * t533 + (t115 * t2 + t311 * (t115 * t459 + t319 * t394)) * t534 + (t113 * t1 + t312 * (t113 * t460 + t320 * t396)) * t535 + ((t152 * t63 + t71 * t484) * t272 - t335 * t473 + (t335 * t140 + t353) * t240) * t461 + ((t151 * t61 + t69 * t485) * t272 - t339 * t474 + (t339 * t139 + t354) * t239) * t462 + ((t150 * t59 + t67 * t486) * t272 - t343 * t475 + (t343 * t138 + t355) * t238) * t463 + ((t146 * t51 + t65 * t488) * t272 - t347 * t476 + (t347 * t137 + t356) * t237) * t464; -g(3) * m(4) + (t128 * t4 + t365 * t84) * t532 + (t127 * t3 + t366 * t83) * t533 + (t126 * t2 + t367 * t82) * t534 + (t122 * t1 + t368 * t75) * t535 + (t84 * t489 + t83 * t490 + t82 * t491 + t75 * t492 + m(4)) * t272 + ((t152 * t88 + t93 * t484) * t272 + t152 * t8 + t12 * t484 + t540 * (t149 * t88 + t93 * t481)) * t461 + ((t151 * t87 + t92 * t485) * t272 + t151 * t7 + t11 * t485 + t539 * (t148 * t87 + t92 * t482)) * t462 + ((t150 * t86 + t91 * t486) * t272 + t150 * t6 + t10 * t486 + t538 * (t147 * t86 + t91 * t483)) * t463 + ((t146 * t85 + t89 * t488) * t272 + t146 * t5 + t9 * t488 + t537 * (t145 * t85 + t89 * t487)) * t464; t537 * t464 * (t145 * t30 + t34 * t487) + t538 * t463 * (t147 * t35 + t38 * t483) + t539 * t462 * (t148 * t36 + t39 * t482) + t540 * t461 * (t149 * t37 + t40 * t481) - (-g(1) * t293 - g(2) * t292) * t245 + t246 * (g(1) * t292 - g(2) * t293) + Ifges(4,3) * t271 + t194 * t273 - t193 * t274 + (t30 * t146 * t464 + t35 * t150 * t463 + t36 * t151 * t462 + t37 * t152 * t461 + t29 * t492 + t31 * t491 + t32 * t490 + t33 * t489 + (t129 * t34 * t464 + t131 * t38 * t463 + t132 * t39 * t462 + t133 * t40 * t461) * t304) * t272 + t106 * t10 + t107 * t11 + t108 * t12 + t109 * t5 + t110 * t6 + t111 * t7 + t112 * t8 + t105 * t9 + t96 * t2 + t97 * t3 + t98 * t4 + t90 * t1 + t367 * t31 * t534 + t366 * t32 * t533 + t365 * t33 * t532 + t368 * t29 * t535;];
tauX  = t25;
