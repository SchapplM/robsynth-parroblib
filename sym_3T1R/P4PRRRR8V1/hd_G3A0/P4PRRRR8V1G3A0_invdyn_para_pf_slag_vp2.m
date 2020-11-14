% Calculate vector of inverse dynamics forces for parallel robot
% P4PRRRR8V1G3A0
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
% Datum: 2020-09-20 23:03
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P4PRRRR8V1G3A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V1G3A0_invdyn_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR8V1G3A0_invdyn_para_pf_slag_vp2: xDP has to be [4x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [4 1]), ...
  'P4PRRRR8V1G3A0_invdyn_para_pf_slag_vp2: xDDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V1G3A0_invdyn_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR8V1G3A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P4PRRRR8V1G3A0_invdyn_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V1G3A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRRR8V1G3A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P4PRRRR8V1G3A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V1G3A0_invdyn_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V1G3A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-09-20 22:59:19
% EndTime: 2020-09-20 22:59:31
% DurationCPUTime: 11.88s
% Computational Cost: add. (39875->619), mult. (112801->1168), div. (13536->9), fcn. (119718->30), ass. (0->466)
t284 = cos(qJ(3,1));
t285 = cos(qJ(2,1));
t279 = sin(qJ(2,1));
t404 = t279 * t284;
t215 = pkin(2) * t404 - pkin(5) * t285;
t257 = sin(pkin(3));
t278 = sin(qJ(3,1));
t259 = cos(pkin(3));
t513 = pkin(2) * t259;
t168 = t215 * t257 + t278 * t513;
t530 = 0.1e1 / t168;
t459 = t530 / t284;
t282 = cos(qJ(3,2));
t283 = cos(qJ(2,2));
t277 = sin(qJ(2,2));
t407 = t277 * t282;
t214 = pkin(2) * t407 - pkin(5) * t283;
t276 = sin(qJ(3,2));
t167 = t214 * t257 + t276 * t513;
t531 = 0.1e1 / t167;
t460 = t531 / t282;
t280 = cos(qJ(3,3));
t281 = cos(qJ(2,3));
t275 = sin(qJ(2,3));
t410 = t275 * t280;
t213 = pkin(2) * t410 - pkin(5) * t281;
t274 = sin(qJ(3,3));
t166 = t213 * t257 + t274 * t513;
t532 = 0.1e1 / t166;
t461 = t532 / t280;
t262 = cos(qJ(3,4));
t263 = cos(qJ(2,4));
t261 = sin(qJ(2,4));
t413 = t261 * t262;
t203 = pkin(2) * t413 - pkin(5) * t263;
t260 = sin(qJ(3,4));
t162 = t203 * t257 + t260 * t513;
t533 = 0.1e1 / t162;
t462 = t533 / t262;
t258 = cos(pkin(6));
t507 = g(3) * t258;
t379 = t262 * mrSges(3,1) - mrSges(3,2) * t260;
t378 = t280 * mrSges(3,1) - mrSges(3,2) * t274;
t377 = t282 * mrSges(3,1) - mrSges(3,2) * t276;
t376 = t284 * mrSges(3,1) - mrSges(3,2) * t278;
t269 = legFrame(1,2);
t239 = sin(t269);
t243 = cos(t269);
t212 = g(1) * t243 - g(2) * t239;
t256 = sin(pkin(6));
t173 = -t212 * t256 - t507;
t208 = g(1) * t239 + g(2) * t243;
t360 = t173 * t259 + t208 * t257;
t268 = legFrame(2,2);
t238 = sin(t268);
t242 = cos(t268);
t211 = g(1) * t242 - g(2) * t238;
t172 = -t211 * t256 - t507;
t207 = g(1) * t238 + g(2) * t242;
t361 = t172 * t259 + t207 * t257;
t267 = legFrame(3,2);
t237 = sin(t267);
t241 = cos(t267);
t210 = g(1) * t241 - g(2) * t237;
t171 = -t210 * t256 - t507;
t206 = g(1) * t237 + g(2) * t241;
t362 = t171 * t259 + t206 * t257;
t266 = legFrame(4,2);
t236 = sin(t266);
t240 = cos(t266);
t209 = g(1) * t240 - g(2) * t236;
t170 = -t209 * t256 - t507;
t205 = g(1) * t236 + g(2) * t240;
t363 = t170 * t259 + t205 * t257;
t290 = xP(4);
t244 = sin(t290);
t245 = cos(t290);
t296 = koppelP(1,2);
t300 = koppelP(1,1);
t198 = t244 * t300 + t245 * t296;
t202 = -t244 * t296 + t245 * t300;
t286 = xDP(4);
t255 = t286 ^ 2;
t270 = xDDP(4);
t273 = xDDP(1);
t144 = -t198 * t270 - t202 * t255 + t273;
t471 = t144 * t243;
t272 = xDDP(2);
t140 = -t198 * t255 + t202 * t270 + t272;
t475 = t140 * t239;
t538 = -t471 + t475;
t295 = koppelP(2,2);
t299 = koppelP(2,1);
t197 = t244 * t299 + t245 * t295;
t201 = -t244 * t295 + t245 * t299;
t143 = -t197 * t270 - t201 * t255 + t273;
t472 = t143 * t242;
t139 = -t197 * t255 + t201 * t270 + t272;
t476 = t139 * t238;
t537 = -t472 + t476;
t294 = koppelP(3,2);
t298 = koppelP(3,1);
t196 = t244 * t298 + t245 * t294;
t200 = -t244 * t294 + t245 * t298;
t142 = -t196 * t270 - t200 * t255 + t273;
t473 = t142 * t241;
t138 = -t196 * t255 + t200 * t270 + t272;
t477 = t138 * t237;
t536 = -t473 + t477;
t293 = koppelP(4,2);
t297 = koppelP(4,1);
t195 = t244 * t297 + t245 * t293;
t199 = -t244 * t293 + t245 * t297;
t141 = -t195 * t270 - t199 * t255 + t273;
t474 = t141 * t240;
t137 = -t195 * t255 + t199 * t270 + t272;
t478 = t137 * t236;
t535 = -t474 + t478;
t534 = 2 * Ifges(3,4);
t246 = t262 ^ 2;
t529 = 0.2e1 * t246;
t249 = t280 ^ 2;
t528 = 0.2e1 * t249;
t251 = t282 ^ 2;
t527 = 0.2e1 * t251;
t253 = t284 ^ 2;
t526 = 0.2e1 * t253;
t525 = Ifges(3,5) / 0.2e1;
t524 = -Ifges(3,6) / 0.2e1;
t423 = t259 * t263;
t424 = t259 * t261;
t512 = pkin(2) * t262;
t129 = (-t256 * t261 + t258 * t423) * t512 + pkin(5) * (t256 * t263 + t258 * t424);
t130 = (t256 * t423 + t258 * t261) * t512 + (t256 * t424 - t258 * t263) * pkin(5);
t287 = xDP(3);
t303 = 0.1e1 / pkin(2);
t288 = xDP(2);
t289 = xDP(1);
t359 = (t199 * t286 + t288) * t236 - (-t195 * t286 + t289) * t240;
t74 = (t359 * t129 + t130 * t287) * t303 * t462;
t523 = pkin(2) * t74;
t419 = t259 * t281;
t422 = t259 * t275;
t511 = pkin(2) * t280;
t131 = (-t256 * t275 + t258 * t419) * t511 + pkin(5) * (t256 * t281 + t258 * t422);
t134 = (t256 * t419 + t258 * t275) * t511 + (t256 * t422 - t258 * t281) * pkin(5);
t358 = (t200 * t286 + t288) * t237 - (-t196 * t286 + t289) * t241;
t79 = (t358 * t131 + t134 * t287) * t303 * t461;
t522 = pkin(2) * t79;
t418 = t259 * t283;
t421 = t259 * t277;
t510 = pkin(2) * t282;
t132 = (-t256 * t277 + t258 * t418) * t510 + pkin(5) * (t256 * t283 + t258 * t421);
t135 = (t256 * t418 + t258 * t277) * t510 + (t256 * t421 - t258 * t283) * pkin(5);
t357 = (t201 * t286 + t288) * t238 - (-t197 * t286 + t289) * t242;
t80 = (t357 * t132 + t135 * t287) * t303 * t460;
t521 = pkin(2) * t80;
t417 = t259 * t285;
t420 = t259 * t279;
t509 = pkin(2) * t284;
t133 = (-t256 * t279 + t258 * t417) * t509 + pkin(5) * (t256 * t285 + t258 * t420);
t136 = (t256 * t417 + t258 * t279) * t509 + (t256 * t420 - t258 * t285) * pkin(5);
t356 = (t202 * t286 + t288) * t239 - (-t198 * t286 + t289) * t243;
t81 = (t356 * t133 + t136 * t287) * t303 * t459;
t520 = pkin(2) * t81;
t323 = t257 * t262 + t260 * t424;
t415 = t260 * t263;
t145 = t323 * t256 - t258 * t415;
t146 = t256 * t415 + t323 * t258;
t95 = (t145 * t287 + t359 * t146) * t462;
t519 = pkin(5) * t95;
t265 = mrSges(2,2) - mrSges(3,3);
t518 = t265 / 0.2e1;
t517 = pkin(2) * t246;
t516 = pkin(2) * t249;
t515 = pkin(2) * t251;
t514 = pkin(2) * t253;
t508 = g(3) * t256;
t301 = pkin(5) ^ 2;
t302 = pkin(2) ^ 2;
t497 = t260 * t74;
t403 = pkin(2) * t497;
t506 = (-pkin(5) * t403 + (t246 * t302 + t301) * t95) * t95;
t505 = Ifges(3,1) + Ifges(2,3);
t322 = t257 * t280 + t274 * t422;
t411 = t274 * t281;
t147 = t322 * t256 - t258 * t411;
t150 = t256 * t411 + t322 * t258;
t102 = (t147 * t287 + t358 * t150) * t461;
t496 = t274 * t79;
t402 = pkin(2) * t496;
t500 = t102 * (-pkin(5) * t402 + (t249 * t302 + t301) * t102);
t321 = t257 * t282 + t276 * t421;
t408 = t276 * t283;
t148 = t321 * t256 - t258 * t408;
t151 = t256 * t408 + t321 * t258;
t103 = (t148 * t287 + t357 * t151) * t460;
t495 = t276 * t80;
t401 = pkin(2) * t495;
t499 = t103 * (-pkin(5) * t401 + (t251 * t302 + t301) * t103);
t320 = t257 * t284 + t278 * t420;
t405 = t278 * t285;
t149 = t320 * t256 - t258 * t405;
t152 = t256 * t405 + t320 * t258;
t104 = (t149 * t287 + t356 * t152) * t459;
t494 = t278 * t81;
t400 = pkin(2) * t494;
t498 = t104 * (-pkin(5) * t400 + (t253 * t302 + t301) * t104);
t493 = t102 * t274;
t492 = t103 * t276;
t491 = t104 * t278;
t204 = pkin(5) * t261 + t263 * t512;
t432 = t257 * t260;
t351 = pkin(2) * t432 - t203 * t259;
t122 = -t204 * t256 + t351 * t258;
t490 = t122 * t533;
t216 = pkin(5) * t275 + t281 * t511;
t430 = t257 * t274;
t350 = pkin(2) * t430 - t213 * t259;
t126 = -t216 * t256 + t350 * t258;
t489 = t126 * t532;
t217 = pkin(5) * t277 + t283 * t510;
t429 = t257 * t276;
t349 = pkin(2) * t429 - t214 * t259;
t127 = -t217 * t256 + t349 * t258;
t488 = t127 * t531;
t218 = pkin(5) * t279 + t285 * t509;
t428 = t257 * t278;
t348 = pkin(2) * t428 - t215 * t259;
t128 = -t218 * t256 + t348 * t258;
t487 = t128 * t530;
t486 = t129 * t303;
t485 = t130 * t303;
t484 = t131 * t303;
t483 = t132 * t303;
t482 = t133 * t303;
t481 = t134 * t303;
t480 = t135 * t303;
t479 = t136 * t303;
t375 = mrSges(3,1) * t260 + mrSges(3,2) * t262;
t153 = -t261 * t257 * t375 + t379 * t259;
t470 = t153 * t533;
t469 = t153 * t303;
t374 = mrSges(3,1) * t274 + mrSges(3,2) * t280;
t154 = -t275 * t257 * t374 + t378 * t259;
t468 = t154 * t532;
t467 = t154 * t303;
t373 = mrSges(3,1) * t276 + mrSges(3,2) * t282;
t155 = -t277 * t257 * t373 + t377 * t259;
t466 = t155 * t531;
t465 = t155 * t303;
t372 = mrSges(3,1) * t278 + mrSges(3,2) * t284;
t156 = -t279 * t257 * t372 + t376 * t259;
t464 = t156 * t530;
t463 = t156 * t303;
t248 = m(1) + m(2) + m(3);
t458 = t533 * t248;
t457 = t532 * t248;
t456 = t531 * t248;
t455 = t530 * t248;
t219 = mrSges(2,1) + t379;
t454 = (t219 * t263 - t261 * t265) * t257;
t220 = mrSges(2,1) + t378;
t449 = (t220 * t281 - t265 * t275) * t257;
t221 = mrSges(2,1) + t377;
t448 = (t221 * t283 - t265 * t277) * t257;
t222 = mrSges(2,1) + t376;
t447 = (t222 * t285 - t265 * t279) * t257;
t445 = t205 * t259;
t443 = t206 * t259;
t441 = t207 * t259;
t439 = t208 * t259;
t223 = Ifges(3,5) * t260 + Ifges(3,6) * t262;
t438 = t223 * t303;
t224 = Ifges(3,5) * t274 + Ifges(3,6) * t280;
t437 = t224 * t303;
t225 = Ifges(3,5) * t276 + Ifges(3,6) * t282;
t436 = t225 * t303;
t226 = Ifges(3,5) * t278 + Ifges(3,6) * t284;
t435 = t226 * t303;
t434 = mrSges(3,2) * t257 * t507;
t433 = t256 * t257;
t431 = t257 * t263;
t427 = t257 * t281;
t426 = t257 * t283;
t425 = t257 * t285;
t416 = t259 * t303;
t264 = Ifges(3,1) - Ifges(3,2);
t414 = t260 * t264;
t412 = t274 * t280;
t409 = t276 * t282;
t406 = t278 * t284;
t399 = t260 * t519;
t398 = pkin(5) * t493;
t397 = pkin(5) * t492;
t396 = pkin(5) * t491;
t395 = t236 * t462;
t394 = t240 * t462;
t393 = t237 * t461;
t392 = t241 * t461;
t391 = t238 * t460;
t390 = t242 * t460;
t389 = t239 * t459;
t388 = t243 * t459;
t387 = t533 * t454;
t386 = t532 * t449;
t385 = t531 * t448;
t384 = t530 * t447;
t383 = t257 * t413;
t382 = t257 * t410;
t381 = t257 * t407;
t380 = t257 * t404;
t371 = t209 * t258 - t508;
t370 = t210 * t258 - t508;
t369 = t211 * t258 - t508;
t368 = t212 * t258 - t508;
t121 = t204 * t258 + t351 * t256;
t113 = t121 * t240 + t162 * t236;
t114 = -t121 * t236 + t162 * t240;
t367 = t113 * t141 + t114 * t137;
t123 = t216 * t258 + t350 * t256;
t115 = t123 * t241 + t166 * t237;
t116 = -t123 * t237 + t166 * t241;
t366 = t115 * t142 + t116 * t138;
t124 = t217 * t258 + t349 * t256;
t117 = t124 * t242 + t167 * t238;
t118 = -t124 * t238 + t167 * t242;
t365 = t117 * t143 + t118 * t139;
t125 = t218 * t258 + t348 * t256;
t119 = t125 * t243 + t168 * t239;
t120 = -t125 * t239 + t168 * t243;
t364 = t119 * t144 + t120 * t140;
t41 = t399 - t523;
t13 = (((t259 * t74 + t95 * t431) * t517 - (t403 - t519) * t383 + t259 * t41) * t95 - (-t74 * t431 + (-t246 * t259 + t260 * t383 + t259) * t95) * t523) * t462;
t17 = (t416 * t506 + (-t74 * t203 * t432 + t259 * (t74 * t517 - t399)) * t74) * t462;
t185 = t260 * t262 * t534 - t246 * t264 + t505;
t21 = (t41 * t523 - t506) * t533;
t5 = -t21 * t454 - t185 * t13 - t223 * t17 + 0.2e1 * t74 * ((t95 * t414 + t74 * t525) * t262 + t497 * t524 + (t529 - 0.1e1) * t95 * Ifges(3,4)) + (-t363 * t219 + t371 * t265) * t263 + t261 * (t371 * t219 + t363 * t265);
t307 = t261 * t363 + t371 * t263;
t94 = t95 ^ 2;
t9 = -t153 * t21 - t223 * t13 - Ifges(3,3) * t17 - (Ifges(3,4) * t529 + t262 * t414 - Ifges(3,4)) * t94 + ((t170 * t257 - t445) * mrSges(3,1) + t307 * mrSges(3,2)) * t262 + t260 * (t434 + (t209 * t433 + t445) * mrSges(3,2) + t307 * mrSges(3,1));
t355 = -t146 * t5 - t9 * t486;
t43 = t398 - t522;
t14 = (((t102 * t427 + t259 * t79) * t516 - (-pkin(5) * t102 + t402) * t382 + t259 * t43) * t102 + (t79 * t427 + (t249 * t259 - t274 * t382 - t259) * t102) * t522) * t461;
t18 = (t416 * t500 + (-t79 * t213 * t430 + t259 * (t79 * t516 - t398)) * t79) * t461;
t22 = (t43 * t522 - t500) * t532;
t306 = t275 * t362 + t370 * t281;
t99 = t102 ^ 2;
t10 = -t154 * t22 - t224 * t14 - Ifges(3,3) * t18 - (Ifges(3,4) * t528 + t264 * t412 - Ifges(3,4)) * t99 + ((t171 * t257 - t443) * mrSges(3,1) + t306 * mrSges(3,2)) * t280 + t274 * (t434 + (t210 * t433 + t443) * mrSges(3,2) + t306 * mrSges(3,1));
t187 = -t249 * t264 + t412 * t534 + t505;
t6 = -t22 * t449 - t187 * t14 - t224 * t18 + 0.2e1 * t79 * ((t264 * t493 + t79 * t525) * t280 + t496 * t524 + (t528 - 0.1e1) * t102 * Ifges(3,4)) + (-t362 * t220 + t370 * t265) * t281 + t275 * (t370 * t220 + t362 * t265);
t354 = -t10 * t484 - t150 * t6;
t100 = t103 ^ 2;
t44 = t397 - t521;
t15 = (((t103 * t426 + t259 * t80) * t515 - (-pkin(5) * t103 + t401) * t381 + t259 * t44) * t103 + (t80 * t426 + (t251 * t259 - t276 * t381 - t259) * t103) * t521) * t460;
t19 = (t416 * t499 + (-t80 * t214 * t429 + t259 * (t80 * t515 - t397)) * t80) * t460;
t23 = (t44 * t521 - t499) * t531;
t305 = t277 * t361 + t369 * t283;
t11 = -t155 * t23 - t225 * t15 - Ifges(3,3) * t19 - (Ifges(3,4) * t527 + t264 * t409 - Ifges(3,4)) * t100 + ((t172 * t257 - t441) * mrSges(3,1) + t305 * mrSges(3,2)) * t282 + t276 * (t434 + (t211 * t433 + t441) * mrSges(3,2) + t305 * mrSges(3,1));
t188 = -t251 * t264 + t409 * t534 + t505;
t7 = -t23 * t448 - t188 * t15 - t225 * t19 + 0.2e1 * t80 * ((t264 * t492 + t80 * t525) * t282 + t495 * t524 + (t527 - 0.1e1) * t103 * Ifges(3,4)) + (-t361 * t221 + t369 * t265) * t283 + t277 * (t369 * t221 + t361 * t265);
t353 = -t11 * t483 - t151 * t7;
t101 = t104 ^ 2;
t45 = t396 - t520;
t16 = (((t104 * t425 + t259 * t81) * t514 - (-pkin(5) * t104 + t400) * t380 + t259 * t45) * t104 - (-t81 * t425 + (-t253 * t259 + t278 * t380 + t259) * t104) * t520) * t459;
t20 = (t416 * t498 + (-t81 * t215 * t428 + t259 * (t81 * t514 - t396)) * t81) * t459;
t24 = (t45 * t520 - t498) * t530;
t304 = t279 * t360 + t368 * t285;
t12 = -t156 * t24 - t226 * t16 - Ifges(3,3) * t20 - (Ifges(3,4) * t526 + t264 * t406 - Ifges(3,4)) * t101 + ((t173 * t257 - t439) * mrSges(3,1) + t304 * mrSges(3,2)) * t284 + t278 * (t434 + (t212 * t433 + t439) * mrSges(3,2) + t304 * mrSges(3,1));
t189 = -t253 * t264 + t406 * t534 + t505;
t8 = -t24 * t447 - t189 * t16 - t226 * t20 + 0.2e1 * t81 * ((t264 * t491 + t81 * t525) * t284 + t494 * t524 + (t526 - 0.1e1) * t104 * Ifges(3,4)) + (-t360 * t222 + t368 * t265) * t285 + t279 * (t368 * t222 + t360 * t265);
t352 = -t12 * t482 - t152 * t8;
t327 = t129 * t438 + t146 * t185;
t51 = t113 * t387 - t327 * t394;
t331 = Ifges(3,3) * t486 + t146 * t223;
t65 = t113 * t470 - t331 * t394;
t346 = t146 * t51 + t65 * t486;
t52 = t114 * t387 + t327 * t395;
t66 = t114 * t470 + t331 * t395;
t345 = t146 * t52 + t66 * t486;
t326 = t131 * t437 + t150 * t187;
t59 = t115 * t386 - t326 * t392;
t330 = Ifges(3,3) * t484 + t150 * t224;
t67 = t115 * t468 - t330 * t392;
t342 = t150 * t59 + t67 * t484;
t60 = t116 * t386 + t326 * t393;
t68 = t116 * t468 + t330 * t393;
t341 = t150 * t60 + t68 * t484;
t325 = t132 * t436 + t151 * t188;
t61 = t117 * t385 - t325 * t390;
t329 = Ifges(3,3) * t483 + t151 * t225;
t69 = t117 * t466 - t329 * t390;
t338 = t151 * t61 + t69 * t483;
t62 = t118 * t385 + t325 * t391;
t70 = t118 * t466 + t329 * t391;
t337 = t151 * t62 + t70 * t483;
t324 = t133 * t435 + t152 * t189;
t63 = t119 * t384 - t324 * t388;
t328 = Ifges(3,3) * t482 + t152 * t226;
t71 = t119 * t464 - t328 * t388;
t334 = t152 * t63 + t71 * t482;
t64 = t120 * t384 + t324 * t389;
t72 = t120 * t464 + t328 * t389;
t333 = t152 * t64 + t72 * t482;
t319 = t129 * t469 + t146 * t454;
t318 = t131 * t467 + t150 * t449;
t317 = t132 * t465 + t151 * t448;
t316 = t133 * t463 + t152 * t447;
t315 = (t195 * t240 + t199 * t236) * t462;
t314 = (t196 * t241 + t200 * t237) * t461;
t313 = (t197 * t242 + t201 * t238) * t460;
t312 = (t198 * t243 + t202 * t239) * t459;
t271 = xDDP(3);
t311 = t122 * t271 + t367;
t310 = t126 * t271 + t366;
t309 = t127 * t271 + t365;
t308 = t128 * t271 + t364;
t292 = mrSges(4,1);
t291 = mrSges(4,2);
t194 = -t244 * t291 + t245 * t292;
t193 = t244 * t292 + t245 * t291;
t112 = t152 * t312;
t111 = t151 * t313;
t110 = t150 * t314;
t109 = t146 * t315;
t108 = t312 * t482;
t107 = t313 * t483;
t106 = t314 * t484;
t105 = t315 * t486;
t98 = (-t119 * t198 + t120 * t202) * t530;
t97 = (-t117 * t197 + t118 * t201) * t531;
t96 = (-t115 * t196 + t116 * t200) * t532;
t93 = t128 * t464 + (Ifges(3,3) * t479 + t149 * t226) * t459;
t92 = t127 * t466 + (Ifges(3,3) * t480 + t148 * t225) * t460;
t91 = t126 * t468 + (Ifges(3,3) * t481 + t147 * t224) * t461;
t90 = (-t113 * t195 + t114 * t199) * t533;
t89 = t122 * t470 + (Ifges(3,3) * t485 + t145 * t223) * t462;
t88 = t128 * t384 + (t136 * t435 + t149 * t189) * t459;
t87 = t127 * t385 + (t135 * t436 + t148 * t188) * t460;
t86 = t126 * t386 + (t134 * t437 + t147 * t187) * t461;
t85 = t122 * t387 + (t130 * t438 + t145 * t185) * t462;
t84 = t128 * t455 + (t136 * t463 + t149 * t447) * t459;
t83 = t127 * t456 + (t135 * t465 + t148 * t448) * t460;
t82 = t126 * t457 + (t134 * t467 + t147 * t449) * t461;
t78 = t81 ^ 2;
t77 = t80 ^ 2;
t76 = t79 ^ 2;
t75 = t122 * t458 + (t130 * t469 + t145 * t454) * t462;
t73 = t74 ^ 2;
t40 = Ifges(3,3) * t108 + t112 * t226 + t156 * t98;
t39 = Ifges(3,3) * t107 + t111 * t225 + t155 * t97;
t38 = Ifges(3,3) * t106 + t110 * t224 + t154 * t96;
t37 = t108 * t226 + t112 * t189 + t98 * t447;
t36 = t107 * t225 + t111 * t188 + t97 * t448;
t35 = t106 * t224 + t110 * t187 + t96 * t449;
t34 = Ifges(3,3) * t105 + t109 * t223 + t153 * t90;
t33 = t108 * t156 + t112 * t447 + t248 * t98;
t32 = t107 * t155 + t111 * t448 + t248 * t97;
t31 = t106 * t154 + t110 * t449 + t248 * t96;
t30 = t105 * t223 + t109 * t185 + t90 * t454;
t29 = t105 * t153 + t109 * t454 + t248 * t90;
t4 = -t16 * t447 - t156 * t20 + ((-mrSges(2,1) * t101 - t376 * (t101 + t78)) * t279 - 0.2e1 * (t104 * t518 + t372 * t81) * t104 * t285) * t257 - t259 * t78 * t372 + (-t24 - t208) * t248;
t3 = -t15 * t448 - t155 * t19 + ((-mrSges(2,1) * t100 - t377 * (t100 + t77)) * t277 - 0.2e1 * (t103 * t518 + t373 * t80) * t103 * t283) * t257 - t259 * t77 * t373 + (-t23 - t207) * t248;
t2 = -t14 * t449 - t154 * t18 + ((-mrSges(2,1) * t99 - t378 * (t99 + t76)) * t275 - 0.2e1 * (t102 * t518 + t374 * t79) * t102 * t281) * t257 - t259 * t76 * t374 + (-t22 - t206) * t248;
t1 = -t13 * t454 - t153 * t17 + ((-mrSges(2,1) * t94 - t379 * (t94 + t73)) * t261 - 0.2e1 * (t375 * t74 + t95 * t518) * t95 * t263) * t257 - t259 * t73 * t375 + (-t21 - t205) * t248;
t25 = [-t193 * t270 - t255 * t194 + (t273 - g(1)) * m(4) + (t119 * t4 + t308 * (t119 * t455 - t316 * t388)) * t530 + (t117 * t3 + t309 * (t117 * t456 - t317 * t390)) * t531 + (t115 * t2 + t310 * (t115 * t457 - t318 * t392)) * t532 + (t113 * t1 + t311 * (t113 * t458 - t319 * t394)) * t533 + ((t149 * t63 + t71 * t479) * t271 + t334 * t475 + (-t334 * t144 + t352) * t243) * t459 + ((t148 * t61 + t69 * t480) * t271 + t338 * t476 + (-t338 * t143 + t353) * t242) * t460 + ((t147 * t59 + t67 * t481) * t271 + t342 * t477 + (-t342 * t142 + t354) * t241) * t461 + ((t145 * t51 + t65 * t485) * t271 + t346 * t478 + (-t346 * t141 + t355) * t240) * t462; -t255 * t193 + t194 * t270 + (t272 - g(2)) * m(4) + (t120 * t4 + t308 * (t120 * t455 + t316 * t389)) * t530 + (t118 * t3 + t309 * (t118 * t456 + t317 * t391)) * t531 + (t116 * t2 + t310 * (t116 * t457 + t318 * t393)) * t532 + (t114 * t1 + t311 * (t114 * t458 + t319 * t395)) * t533 + ((t149 * t64 + t72 * t479) * t271 - t333 * t471 + (t333 * t140 - t352) * t239) * t459 + ((t148 * t62 + t70 * t480) * t271 - t337 * t472 + (t337 * t139 - t353) * t238) * t460 + ((t147 * t60 + t68 * t481) * t271 - t341 * t473 + (t341 * t138 - t354) * t237) * t461 + ((t145 * t52 + t66 * t485) * t271 - t345 * t474 + (t345 * t137 - t355) * t236) * t462; -g(3) * m(4) + (t128 * t4 + t364 * t84) * t530 + (t127 * t3 + t365 * t83) * t531 + (t126 * t2 + t366 * t82) * t532 + (t122 * t1 + t367 * t75) * t533 + (t487 * t84 + t488 * t83 + t489 * t82 + t490 * t75 + m(4)) * t271 + ((t149 * t88 + t479 * t93) * t271 + t12 * t479 + t149 * t8 + t538 * (t152 * t88 + t93 * t482)) * t459 + ((t148 * t87 + t480 * t92) * t271 + t11 * t480 + t148 * t7 + t537 * (t151 * t87 + t92 * t483)) * t460 + ((t147 * t86 + t481 * t91) * t271 + t10 * t481 + t147 * t6 + t536 * (t150 * t86 + t91 * t484)) * t461 + ((t145 * t85 + t485 * t89) * t271 + t9 * t485 + t145 * t5 + t535 * (t146 * t85 + t89 * t486)) * t462; t538 * t459 * (t152 * t37 + t40 * t482) + t535 * t462 * (t146 * t30 + t34 * t486) + t536 * t461 * (t150 * t35 + t38 * t484) + t537 * t460 * (t151 * t36 + t39 * t483) + (t30 * t145 * t462 + t35 * t147 * t461 + t36 * t148 * t460 + t37 * t149 * t459 + t29 * t490 + t31 * t489 + t32 * t488 + t33 * t487 + (t130 * t34 * t462 + t134 * t38 * t461 + t135 * t39 * t460 + t136 * t40 * t459) * t303) * t271 - (-g(1) * t292 - g(2) * t291) * t244 + t245 * (g(1) * t291 - g(2) * t292) + t109 * t5 + t110 * t6 + t111 * t7 + t112 * t8 + t105 * t9 + t106 * t10 + t107 * t11 + t108 * t12 + t96 * t2 + t97 * t3 + t98 * t4 + t90 * t1 + t366 * t31 * t532 + t365 * t32 * t531 + t364 * t33 * t530 + t367 * t29 * t533 + t194 * t272 - t193 * t273 + Ifges(4,3) * t270;];
tauX  = t25;
