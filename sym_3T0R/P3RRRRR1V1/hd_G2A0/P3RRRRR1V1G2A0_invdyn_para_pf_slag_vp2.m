% Calculate vector of inverse dynamics forces for parallel robot
% P3RRRRR1V1G2A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1]';
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
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 03:36
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRRRR1V1G2A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR1V1G2A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRRRR1V1G2A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRRRR1V1G2A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR1V1G2A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR1V1G2A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRRRR1V1G2A0_invdyn_para_pf_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR1V1G2A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR1V1G2A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRRRR1V1G2A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR1V1G2A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR1V1G2A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:33:43
% EndTime: 2020-08-07 03:33:55
% DurationCPUTime: 12.87s
% Computational Cost: add. (44397->746), mult. (79572->1203), div. (18684->8), fcn. (84957->36), ass. (0->467)
t250 = qJ(2,3) + qJ(3,3);
t160 = cos(qJ(1,3) - t250) + cos(qJ(1,3) + t250);
t251 = qJ(2,2) + qJ(3,2);
t161 = cos(qJ(1,2) - t251) + cos(qJ(1,2) + t251);
t252 = qJ(2,1) + qJ(3,1);
t162 = cos(qJ(1,1) - t252) + cos(qJ(1,1) + t252);
t270 = cos(qJ(3,3));
t241 = t270 ^ 2;
t273 = cos(qJ(3,2));
t244 = t273 ^ 2;
t276 = cos(qJ(3,1));
t247 = t276 ^ 2;
t272 = cos(qJ(1,3));
t518 = 0.2e1 * t272 ^ 2;
t275 = cos(qJ(1,2));
t517 = 0.2e1 * t275 ^ 2;
t278 = cos(qJ(1,1));
t516 = 0.2e1 * t278 ^ 2;
t284 = xDP(1);
t481 = pkin(3) * t276;
t204 = pkin(2) + t481;
t268 = sin(qJ(2,1));
t267 = sin(qJ(3,1));
t277 = cos(qJ(2,1));
t406 = t267 * t277;
t358 = pkin(3) * t406;
t138 = t204 * t268 + t358;
t177 = t204 * t277;
t407 = t267 * t268;
t359 = pkin(3) * t407;
t141 = t177 - t359;
t257 = legFrame(1,2);
t224 = sin(t257);
t227 = cos(t257);
t269 = sin(qJ(1,1));
t425 = t227 * t269;
t123 = t138 * t224 - t141 * t425;
t240 = 0.1e1 / t267;
t290 = 0.1e1 / pkin(3);
t292 = 0.1e1 / pkin(2);
t400 = t290 * t292;
t334 = t240 * t400;
t328 = t123 * t334;
t102 = t284 * t328;
t283 = xDP(2);
t428 = t224 * t269;
t126 = t138 * t227 + t141 * t428;
t325 = t126 * t334;
t105 = t283 * t325;
t282 = xDP(3);
t510 = t141 * t278;
t313 = t334 * t510;
t300 = t282 * t313;
t69 = t102 + t105 - t300;
t422 = t240 * t292;
t331 = t422 / 0.2e1;
t322 = t162 * t331;
t304 = -t276 * t277 + t407;
t305 = -t268 * t276 - t406;
t120 = t224 * t305 - t304 * t425;
t346 = t120 * t422;
t119 = t227 * t305 + t304 * t428;
t347 = t119 * t422;
t81 = t282 * t322 + t283 * t347 + t284 * t346;
t51 = t69 + t81;
t493 = pkin(3) * t51;
t515 = -0.2e1 * t493;
t482 = pkin(3) * t273;
t203 = pkin(2) + t482;
t265 = sin(qJ(2,2));
t264 = sin(qJ(3,2));
t274 = cos(qJ(2,2));
t410 = t264 * t274;
t360 = pkin(3) * t410;
t137 = t203 * t265 + t360;
t176 = t203 * t274;
t411 = t264 * t265;
t361 = pkin(3) * t411;
t140 = t176 - t361;
t256 = legFrame(2,2);
t223 = sin(t256);
t226 = cos(t256);
t266 = sin(qJ(1,2));
t426 = t226 * t266;
t122 = t137 * t223 - t140 * t426;
t239 = 0.1e1 / t264;
t335 = t239 * t400;
t329 = t122 * t335;
t101 = t284 * t329;
t429 = t223 * t266;
t125 = t137 * t226 + t140 * t429;
t326 = t125 * t335;
t104 = t283 * t326;
t511 = t140 * t275;
t314 = t335 * t511;
t301 = t282 * t314;
t68 = t101 + t104 - t301;
t423 = t239 * t292;
t332 = t423 / 0.2e1;
t323 = t161 * t332;
t306 = -t273 * t274 + t411;
t307 = -t265 * t273 - t410;
t118 = t223 * t307 - t306 * t426;
t348 = t118 * t423;
t117 = t226 * t307 + t306 * t429;
t349 = t117 * t423;
t80 = t282 * t323 + t283 * t349 + t284 * t348;
t50 = t68 + t80;
t494 = pkin(3) * t50;
t514 = -0.2e1 * t494;
t483 = pkin(3) * t270;
t202 = pkin(2) + t483;
t262 = sin(qJ(2,3));
t261 = sin(qJ(3,3));
t271 = cos(qJ(2,3));
t414 = t261 * t271;
t362 = pkin(3) * t414;
t136 = t202 * t262 + t362;
t175 = t202 * t271;
t415 = t261 * t262;
t363 = pkin(3) * t415;
t139 = t175 - t363;
t255 = legFrame(3,2);
t222 = sin(t255);
t225 = cos(t255);
t263 = sin(qJ(1,3));
t427 = t225 * t263;
t121 = t136 * t222 - t139 * t427;
t238 = 0.1e1 / t261;
t336 = t238 * t400;
t330 = t121 * t336;
t100 = t284 * t330;
t430 = t222 * t263;
t124 = t136 * t225 + t139 * t430;
t327 = t124 * t336;
t103 = t283 * t327;
t512 = t139 * t272;
t315 = t336 * t512;
t302 = t282 * t315;
t67 = t100 + t103 - t302;
t424 = t238 * t292;
t333 = t424 / 0.2e1;
t324 = t160 * t333;
t308 = -t270 * t271 + t415;
t309 = -t262 * t270 - t414;
t116 = t222 * t309 - t308 * t427;
t350 = t116 * t424;
t115 = t225 * t309 + t308 * t430;
t351 = t115 * t424;
t79 = t282 * t324 + t283 * t351 + t284 * t350;
t49 = t67 + t79;
t495 = pkin(3) * t49;
t513 = -0.2e1 * t495;
t254 = (Ifges(2,4) - Ifges(3,4));
t509 = 2 * pkin(1);
t212 = cos(t250);
t154 = 0.1e1 / (pkin(2) * t271 + pkin(3) * t212 + pkin(1));
t97 = (-t263 * t282 + (-t222 * t283 + t225 * t284) * t272) * t154;
t508 = -0.2e1 * t97;
t213 = cos(t251);
t155 = 0.1e1 / (pkin(2) * t274 + pkin(3) * t213 + pkin(1));
t98 = (-t266 * t282 + (-t223 * t283 + t226 * t284) * t275) * t155;
t507 = -0.2e1 * t98;
t214 = cos(t252);
t156 = 0.1e1 / (pkin(2) * t277 + pkin(3) * t214 + pkin(1));
t99 = (-t269 * t282 + (-t224 * t283 + t227 * t284) * t278) * t156;
t506 = -0.2e1 * t99;
t253 = Ifges(3,2) - Ifges(3,1);
t291 = pkin(2) ^ 2;
t287 = m(3) * t291;
t188 = -Ifges(2,1) + Ifges(2,2) + t287 - t253;
t505 = -0.2e1 * t188;
t289 = pkin(3) ^ 2;
t504 = -0.2e1 * t289;
t285 = m(2) + m(3);
t281 = pkin(1) * mrSges(3,1);
t503 = pkin(2) * mrSges(3,1);
t502 = pkin(2) * mrSges(3,2);
t501 = pkin(2) * t79 ^ 2;
t500 = pkin(2) * t80 ^ 2;
t499 = pkin(2) * t81 ^ 2;
t498 = pkin(2) * t79;
t497 = pkin(2) * t80;
t496 = pkin(2) * t81;
t492 = t160 / 0.2e1;
t491 = t161 / 0.2e1;
t490 = t162 / 0.2e1;
t489 = pkin(1) * t262;
t488 = pkin(1) * t265;
t487 = pkin(1) * t268;
t486 = pkin(3) * t241;
t485 = pkin(3) * t244;
t484 = pkin(3) * t247;
t230 = t270 * pkin(2);
t231 = t273 * pkin(2);
t232 = t276 * pkin(2);
t480 = t49 * t67;
t479 = t50 * t68;
t478 = t51 * t69;
t477 = 2 * t254;
t476 = mrSges(3,1) * t270;
t475 = mrSges(3,1) * t273;
t474 = mrSges(3,1) * t276;
t473 = mrSges(3,2) * t261;
t472 = mrSges(3,2) * t264;
t471 = mrSges(3,2) * t267;
t470 = Ifges(3,4) * t261;
t469 = Ifges(3,4) * t264;
t468 = Ifges(3,4) * t267;
t464 = t241 * Ifges(3,4);
t463 = t241 * t49;
t242 = t271 ^ 2;
t94 = t97 ^ 2;
t462 = t242 * t94;
t461 = t244 * Ifges(3,4);
t245 = t274 ^ 2;
t95 = t98 ^ 2;
t460 = t245 * t95;
t459 = t247 * Ifges(3,4);
t248 = t277 ^ 2;
t96 = t99 ^ 2;
t458 = t248 * t96;
t216 = t261 * mrSges(3,1);
t457 = t261 * t49;
t217 = t264 * mrSges(3,1);
t456 = t264 * t50;
t218 = t267 * mrSges(3,1);
t455 = t267 * t51;
t454 = t271 * t94;
t453 = t274 * t95;
t452 = t277 * t96;
t451 = t50 * t244;
t450 = t51 * t247;
t228 = m(3) * pkin(2) + mrSges(2,1);
t190 = t216 + mrSges(2,2);
t191 = t217 + mrSges(2,2);
t192 = t218 + mrSges(2,2);
t449 = pkin(2) * mrSges(3,3) - Ifges(2,5);
t448 = t121 * t290;
t447 = t122 * t290;
t446 = t123 * t290;
t445 = t124 * t290;
t444 = t125 * t290;
t443 = t126 * t290;
t130 = (-Ifges(3,5) * t271 + Ifges(3,6) * t262) * t261 - t270 * (Ifges(3,5) * t262 + Ifges(3,6) * t271);
t442 = t130 * t290;
t131 = (-Ifges(3,5) * t274 + Ifges(3,6) * t265) * t264 - t273 * (Ifges(3,5) * t265 + Ifges(3,6) * t274);
t441 = t131 * t290;
t132 = (-Ifges(3,5) * t277 + Ifges(3,6) * t268) * t267 - t276 * (Ifges(3,5) * t268 + Ifges(3,6) * t277);
t440 = t132 * t290;
t439 = t154 * t263;
t438 = t154 * t272;
t437 = t155 * t266;
t436 = t155 * t275;
t435 = t156 * t269;
t434 = t156 * t278;
t206 = pkin(2) * t473;
t381 = pkin(2) * t476;
t157 = Ifges(3,3) - t206 + t381;
t433 = t157 * t290;
t207 = pkin(2) * t472;
t380 = pkin(2) * t475;
t158 = Ifges(3,3) - t207 + t380;
t432 = t158 * t290;
t208 = pkin(2) * t471;
t379 = pkin(2) * t474;
t159 = Ifges(3,3) - t208 + t379;
t431 = t159 * t290;
t421 = t253 * t241;
t420 = t253 * t244;
t419 = t253 * t247;
t418 = t253 * t261;
t417 = t253 * t264;
t416 = t253 * t267;
t413 = t262 * t271;
t412 = t263 * t272;
t409 = t265 * t274;
t408 = t266 * t275;
t405 = t268 * t277;
t404 = t269 * t278;
t280 = pkin(1) * mrSges(3,2);
t403 = t280 * t270;
t402 = t280 * t273;
t401 = t280 * t276;
t399 = t289 - t291;
t398 = -2 * t281;
t397 = t97 * t509;
t396 = t98 * t509;
t395 = t99 * t509;
t394 = -0.2e1 * pkin(1) * t228;
t393 = pkin(3) * t230;
t392 = pkin(3) * t231;
t391 = pkin(3) * t232;
t390 = 0.4e1 * t470;
t389 = 0.4e1 * t469;
t388 = 0.4e1 * t468;
t387 = 0.2e1 * t464;
t386 = -0.4e1 * t463;
t385 = 0.2e1 * t461;
t384 = 0.2e1 * t459;
t383 = -0.4e1 * t451;
t382 = -0.4e1 * t450;
t378 = t94 * t489;
t377 = t95 * t488;
t376 = t96 * t487;
t375 = t49 * t483;
t374 = t50 * t482;
t373 = t51 * t481;
t372 = pkin(2) * t216;
t371 = pkin(2) * t217;
t370 = pkin(2) * t218;
t369 = t79 * t230;
t368 = t80 * t231;
t367 = t81 * t232;
t366 = 0.2e1 * t421;
t365 = 0.2e1 * t420;
t364 = 0.2e1 * t419;
t357 = -t503 / 0.4e1;
t356 = t503 / 0.2e1;
t355 = t498 / 0.2e1;
t354 = t497 / 0.2e1;
t353 = t496 / 0.2e1;
t352 = Ifges(2,3) + Ifges(3,3) + t287;
t345 = t290 * t512;
t344 = t290 * t511;
t343 = t290 * t510;
t342 = t222 * t438;
t341 = t225 * t438;
t340 = t223 * t436;
t339 = t226 * t436;
t338 = t224 * t434;
t337 = t227 * t434;
t321 = mrSges(3,1) * t378 - t418 * t94;
t320 = mrSges(3,1) * t377 - t417 * t95;
t319 = mrSges(3,1) * t376 - t416 * t96;
t178 = t228 - t473;
t179 = t228 - t472;
t180 = t228 - t471;
t166 = g(1) * t225 - g(2) * t222;
t318 = g(3) * t272 + t166 * t263;
t167 = g(1) * t226 - g(2) * t223;
t317 = g(3) * t275 + t167 * t266;
t168 = g(1) * t227 - g(2) * t224;
t316 = g(3) * t278 + t168 * t269;
t148 = -t178 - t476;
t169 = mrSges(3,2) * t270 + t190;
t312 = -t148 * t271 - t169 * t262;
t149 = -t179 - t475;
t170 = mrSges(3,2) * t273 + t191;
t311 = -t149 * t274 - t170 * t265;
t150 = -t180 - t474;
t171 = mrSges(3,2) * t276 + t192;
t310 = -t150 * t277 - t171 * t268;
t303 = (pkin(1) ^ 2) * t285 + Ifges(2,1) + Ifges(3,2) + Ifges(1,3);
t299 = pkin(1) * t415 - pkin(3) + t486;
t298 = pkin(1) * t411 - pkin(3) + t485;
t297 = pkin(1) * t407 - pkin(3) + t484;
t193 = -0.2e1 * t206;
t236 = 0.2e1 * t503;
t296 = (t236 + t390) * t270 + t188 + t193 + t366;
t194 = -0.2e1 * t207;
t295 = (t236 + t389) * t273 + t188 + t194 + t365;
t195 = -0.2e1 * t208;
t294 = (t236 + t388) * t276 + t188 + t195 + t364;
t286 = -Ifges(3,4) / 0.2e1;
t260 = xDDP(1);
t259 = xDDP(2);
t258 = xDDP(3);
t237 = 2 * t281;
t235 = -0.2e1 * t502;
t234 = mrSges(2,3) + mrSges(3,3) + mrSges(1,2);
t233 = -t502 / 0.4e1;
t229 = Ifges(3,1) / 0.2e1 - Ifges(3,2) / 0.2e1;
t211 = sin(t252);
t210 = sin(t251);
t209 = sin(t250);
t205 = t234 * g(3);
t189 = pkin(1) * t285 + mrSges(1,1);
t184 = t189 * g(3);
t183 = pkin(1) - 0.2e1 * t359;
t182 = pkin(1) - 0.2e1 * t361;
t181 = pkin(1) - 0.2e1 * t363;
t174 = -pkin(3) + t232 + 0.2e1 * t484;
t173 = -pkin(3) + t231 + 0.2e1 * t485;
t172 = -pkin(3) + t230 + 0.2e1 * t486;
t165 = g(1) * t224 + g(2) * t227;
t164 = g(1) * t223 + g(2) * t226;
t163 = g(1) * t222 + g(2) * t225;
t144 = t195 + t352 + 0.2e1 * t379;
t143 = t194 + t352 + 0.2e1 * t380;
t142 = t193 + t352 + 0.2e1 * t381;
t129 = (-Ifges(3,5) * t276 + Ifges(3,6) * t267 + t449) * t268 - (Ifges(3,5) * t267 + Ifges(3,6) * t276 + Ifges(2,6)) * t277;
t128 = (-Ifges(3,5) * t273 + Ifges(3,6) * t264 + t449) * t265 - (Ifges(3,5) * t264 + Ifges(3,6) * t273 + Ifges(2,6)) * t274;
t127 = (-Ifges(3,5) * t270 + Ifges(3,6) * t261 + t449) * t262 - (Ifges(3,5) * t261 + Ifges(3,6) * t270 + Ifges(2,6)) * t271;
t93 = t96 * t384;
t92 = t95 * t385;
t91 = t94 * t387;
t87 = -t132 * t435 + (-Ifges(3,3) * t343 + t159 * t490) * t422;
t86 = -t131 * t437 + (-Ifges(3,3) * t344 + t158 * t491) * t423;
t85 = -t130 * t439 + (-Ifges(3,3) * t345 + t157 * t492) * t424;
t84 = -t129 * t435 + (t144 * t490 - t159 * t343) * t422;
t83 = -t128 * t437 + (t143 * t491 - t158 * t344) * t423;
t82 = -t127 * t439 + (t142 * t492 - t157 * t345) * t424;
t75 = t254 * t81;
t74 = t254 * t80;
t73 = t254 * t79;
t72 = t294 * t248 + (t180 * t509 + t237 * t276 + (0.4e1 * t459 + t235 * t276 + 0.2e1 * (-t253 * t276 - t503) * t267 + t477) * t268) * t277 - t419 + 0.2e1 * (-mrSges(3,2) * t487 - t468) * t276 - 0.2e1 * t192 * t487 + t303;
t71 = t295 * t245 + (t179 * t509 + t237 * t273 + (0.4e1 * t461 + t235 * t273 + 0.2e1 * (-t253 * t273 - t503) * t264 + t477) * t265) * t274 - t420 + 0.2e1 * (-mrSges(3,2) * t488 - t469) * t273 - 0.2e1 * t191 * t488 + t303;
t70 = t296 * t242 + (t178 * t509 + t237 * t270 + (0.4e1 * t464 + t235 * t270 + 0.2e1 * (-t253 * t270 - t503) * t261 + t477) * t262) * t271 - t421 + 0.2e1 * (-mrSges(3,2) * t489 - t470) * t270 - 0.2e1 * t190 * t489 + t303;
t66 = t132 * t337 + (Ifges(3,3) * t446 + t120 * t159) * t422;
t65 = t131 * t339 + (Ifges(3,3) * t447 + t118 * t158) * t423;
t64 = t130 * t341 + (Ifges(3,3) * t448 + t116 * t157) * t424;
t63 = -t132 * t338 + (Ifges(3,3) * t443 + t119 * t159) * t422;
t62 = -t131 * t340 + (Ifges(3,3) * t444 + t117 * t158) * t423;
t61 = -t130 * t342 + (Ifges(3,3) * t445 + t115 * t157) * t424;
t60 = t129 * t337 + (t120 * t144 + t123 * t431) * t422;
t59 = t128 * t339 + (t118 * t143 + t122 * t432) * t423;
t58 = t127 * t341 + (t116 * t142 + t121 * t433) * t424;
t57 = -t129 * t338 + (t119 * t144 + t126 * t431) * t422;
t56 = -t128 * t340 + (t117 * t143 + t125 * t432) * t423;
t55 = -t127 * t342 + (t115 * t142 + t124 * t433) * t424;
t54 = -t72 * t435 + (t129 * t490 - t132 * t343) * t422;
t53 = -t71 * t437 + (t128 * t491 - t131 * t344) * t423;
t52 = -t70 * t439 + (t127 * t492 - t130 * t345) * t424;
t48 = t102 / 0.2e1 + t105 / 0.2e1 - t300 / 0.2e1 + t81;
t47 = t101 / 0.2e1 + t104 / 0.2e1 - t301 / 0.2e1 + t80;
t46 = t100 / 0.2e1 + t103 / 0.2e1 - t302 / 0.2e1 + t79;
t45 = t51 * t416;
t44 = t50 * t417;
t43 = t49 * t418;
t42 = t72 * t337 + (t120 * t129 + t123 * t440) * t422;
t41 = t71 * t339 + (t118 * t128 + t122 * t441) * t423;
t40 = t70 * t341 + (t116 * t127 + t121 * t442) * t424;
t39 = -t72 * t338 + (t119 * t129 + t126 * t440) * t422;
t38 = -t71 * t340 + (t117 * t128 + t125 * t441) * t423;
t37 = -t70 * t342 + (t115 * t127 + t124 * t442) * t424;
t36 = Ifges(3,6) * t51 + t398 * t99;
t35 = Ifges(3,6) * t50 + t398 * t98;
t34 = Ifges(3,6) * t49 + t398 * t97;
t33 = mrSges(3,2) * t395 + Ifges(3,5) * t51;
t32 = mrSges(3,2) * t396 + Ifges(3,5) * t50;
t31 = mrSges(3,2) * t397 + Ifges(3,5) * t49;
t30 = t373 + t496;
t29 = t374 + t497;
t28 = t375 + t498;
t27 = t96 + (0.2e1 * t81 + t69) * t69;
t26 = t95 + (0.2e1 * t80 + t68) * t68;
t25 = t94 + (0.2e1 * t79 + t67) * t67;
t18 = (t211 * t515 - 0.2e1 * t268 * t496) * t99 * t156;
t17 = (t210 * t514 - 0.2e1 * t265 * t497) * t98 * t155;
t16 = (t209 * t513 - 0.2e1 * t262 * t498) * t97 * t154;
t15 = ((t232 + pkin(3)) * t478 + (-((t247 * t504 - 0.2e1 * t391 + t399) * t248 - t183 * t177 + pkin(3) * t297) * t96 + (t289 * t51 + t291 * t81 + 0.2e1 * t391 * t48) * t81) * t290) * t422;
t14 = ((t231 + pkin(3)) * t479 + (-((t244 * t504 - 0.2e1 * t392 + t399) * t245 - t182 * t176 + pkin(3) * t298) * t95 + (t289 * t50 + t291 * t80 + 0.2e1 * t392 * t47) * t80) * t290) * t423;
t13 = ((t230 + pkin(3)) * t480 + (-((t241 * t504 - 0.2e1 * t393 + t399) * t242 - t181 * t175 + pkin(3) * t299) * t94 + (t289 * t49 + t291 * t79 + 0.2e1 * t393 * t46) * t79) * t290) * t424;
t12 = -pkin(3) * t422 * t478 + ((-0.4e1 * ((t276 * t353 + (t247 - 0.1e1 / 0.2e1) * t493) * t405 + ((t353 + t373) * t248 - t30 / 0.2e1) * t267) * t404 + (t516 - 0.2e1) * (t174 * t248 + (-pkin(2) * t407 + t183 * t276) * t277 - t297) * t99 + t162 * ((t268 * t30 + t358 * t51) * t269 - t278 * (pkin(1) + t141) * t99)) * t99 + (((t367 + (0.2e1 * t247 - 0.1e1) * t493) * t248 - (0.2e1 * t373 + t496) * t267 * t405 + t493 - t493 * t247) * t516 + (t174 * t405 + ((pkin(2) + 0.2e1 * t481) * t248 - t204) * t267) * t404 * t506 - 0.2e1 * t367 + t515 + t162 * ((-t277 * t30 + t359 * t51) * t278 + t269 * t138 * t99)) * t81) * t331;
t11 = -pkin(3) * t423 * t479 + ((-0.4e1 * ((t273 * t354 + (t244 - 0.1e1 / 0.2e1) * t494) * t409 + ((t354 + t374) * t245 - t29 / 0.2e1) * t264) * t408 + (t517 - 0.2e1) * (t173 * t245 + (-pkin(2) * t411 + t182 * t273) * t274 - t298) * t98 + t161 * ((t265 * t29 + t360 * t50) * t266 - t275 * (pkin(1) + t140) * t98)) * t98 + (((t368 + (0.2e1 * t244 - 0.1e1) * t494) * t245 - (0.2e1 * t374 + t497) * t264 * t409 + t494 - t494 * t244) * t517 + (t173 * t409 + ((pkin(2) + 0.2e1 * t482) * t245 - t203) * t264) * t408 * t507 - 0.2e1 * t368 + t514 + t161 * ((-t274 * t29 + t361 * t50) * t275 + t266 * t137 * t98)) * t80) * t332;
t10 = -pkin(3) * t424 * t480 + ((-0.4e1 * ((t270 * t355 + (t241 - 0.1e1 / 0.2e1) * t495) * t413 + ((t355 + t375) * t242 - t28 / 0.2e1) * t261) * t412 + (t518 - 0.2e1) * (t172 * t242 + (-pkin(2) * t415 + t181 * t270) * t271 - t299) * t97 + t160 * ((t262 * t28 + t362 * t49) * t263 - t272 * (pkin(1) + t139) * t97)) * t97 + (((t369 + (0.2e1 * t241 - 0.1e1) * t495) * t242 - (0.2e1 * t375 + t498) * t261 * t413 + t495 - t495 * t241) * t518 + (t172 * t413 + ((pkin(2) + 0.2e1 * t483) * t242 - t202) * t261) * t412 * t508 - 0.2e1 * t369 + t513 + t160 * ((-t271 * t28 + t363 * t49) * t272 + t263 * t136 * t97)) * t79) * t333;
t9 = -t132 * t18 - t159 * t12 - Ifges(3,3) * t15 - 0.4e1 * (t459 + (t229 * t267 + t233) * t276 + t267 * t357 + t286) * t458 + (t281 * t267 + t401 + (t364 + (t388 + t503) * t276 - t208 - t253) * t268) * t452 + t93 + (mrSges(3,2) * t499 + t319) * t276 + (mrSges(3,1) * t499 - mrSges(3,2) * t376) * t267 - Ifges(3,4) * t96 + (mrSges(3,1) * t165 + mrSges(3,2) * t316) * t214 + t211 * (mrSges(3,1) * t316 - mrSges(3,2) * t165);
t8 = -t131 * t17 - t158 * t11 - Ifges(3,3) * t14 - 0.4e1 * (t461 + (t229 * t264 + t233) * t273 + t264 * t357 + t286) * t460 + (t281 * t264 + t402 + (t365 + (t389 + t503) * t273 - t207 - t253) * t265) * t453 + t92 + (mrSges(3,2) * t500 + t320) * t273 + (mrSges(3,1) * t500 - mrSges(3,2) * t377) * t264 - Ifges(3,4) * t95 + (mrSges(3,1) * t164 + mrSges(3,2) * t317) * t213 + t210 * (mrSges(3,1) * t317 - mrSges(3,2) * t164);
t7 = -t130 * t16 - t157 * t10 - Ifges(3,3) * t13 - 0.4e1 * (t464 + (t229 * t261 + t233) * t270 + t261 * t357 + t286) * t462 + (t281 * t261 + t403 + (t366 + (t390 + t503) * t270 - t206 - t253) * t262) * t454 + t91 + (mrSges(3,2) * t501 + t321) * t270 + (mrSges(3,1) * t501 - mrSges(3,2) * t378) * t261 - Ifges(3,4) * t94 + (mrSges(3,1) * t163 + mrSges(3,2) * t318) * t212 + t209 * (mrSges(3,1) * t318 - mrSges(3,2) * t163);
t6 = -t129 * t18 - t144 * t12 - t159 * t15 - 0.2e1 * (t384 + (-t416 - t502) * t276 - t370 + t254) * t458 + (pkin(1) * t192 + t268 * t294 + t401) * t452 + t93 + (-t27 * t502 + t319) * t276 + t180 * t376 - t27 * t370 + t96 * t254 + (-t150 * t316 - t165 * t171) * t268 + t277 * (-t150 * t165 + t171 * t316);
t5 = -t128 * t17 - t143 * t11 - t158 * t14 - 0.2e1 * (t385 + (-t417 - t502) * t273 - t371 + t254) * t460 + (pkin(1) * t191 + t265 * t295 + t402) * t453 + t92 + (-t26 * t502 + t320) * t273 + t179 * t377 - t26 * t371 + t95 * t254 + (-t149 * t317 - t164 * t170) * t265 + t274 * (-t149 * t164 + t170 * t317);
t4 = -t127 * t16 - t142 * t10 - t157 * t13 - 0.2e1 * (t387 + (-t418 - t502) * t270 - t372 + t254) * t462 + (pkin(1) * t190 + t262 * t296 + t403) * t454 + t91 + (-t25 * t502 + t321) * t270 + t178 * t378 - t25 * t372 + t94 * t254 + (-t148 * t318 - t163 * t169) * t262 + t271 * (-t148 * t163 + t169 * t318);
t3 = -t72 * t18 - t129 * t12 - t132 * t15 + 0.4e1 * t99 * ((-t48 * t502 - t45) * t276 - t48 * t370 + t75 + (-t69 + 0.2e1 * t450) * Ifges(3,4)) * t248 + (-(mrSges(2,2) * t395 - t449 * t81) * t81 + (t36 * t267 - t276 * t33) * t51 + (-0.8e1 * (Ifges(3,4) * t455 + t356 * t48) * t276 + 0.4e1 * t48 * t208 + t81 * t505 + (0.2e1 * t69 + t382) * t253) * t99 * t268) * t277 + Ifges(3,4) * t99 * t382 + (t36 * t51 * t268 + (-mrSges(3,2) * t496 - t45) * t506) * t276 + (t33 * t455 + (Ifges(2,6) * t81 + t394 * t99) * t81) * t268 + (-Ifges(3,4) * t69 - t370 * t81 + t75) * t506 + (t205 + (-t189 - t310) * t168) * t278 + (g(3) * t310 + t168 * t234 + t184) * t269;
t2 = -t71 * t17 - t128 * t11 - t131 * t14 + 0.4e1 * t98 * ((-t47 * t502 - t44) * t273 - t47 * t371 + t74 + (-t68 + 0.2e1 * t451) * Ifges(3,4)) * t245 + (-(mrSges(2,2) * t396 - t449 * t80) * t80 + (t264 * t35 - t273 * t32) * t50 + (-0.8e1 * (Ifges(3,4) * t456 + t356 * t47) * t273 + 0.4e1 * t47 * t207 + t80 * t505 + (0.2e1 * t68 + t383) * t253) * t98 * t265) * t274 + Ifges(3,4) * t98 * t383 + (t35 * t50 * t265 + (-mrSges(3,2) * t497 - t44) * t507) * t273 + (t32 * t456 + (Ifges(2,6) * t80 + t394 * t98) * t80) * t265 + (-Ifges(3,4) * t68 - t371 * t80 + t74) * t507 + (t205 + (-t189 - t311) * t167) * t275 + (g(3) * t311 + t167 * t234 + t184) * t266;
t1 = -t70 * t16 - t127 * t10 - t130 * t13 + 0.4e1 * t97 * ((-t46 * t502 - t43) * t270 - t46 * t372 + t73 + (-t67 + 0.2e1 * t463) * Ifges(3,4)) * t242 + (-(mrSges(2,2) * t397 - t449 * t79) * t79 + (t261 * t34 - t270 * t31) * t49 + (-0.8e1 * (Ifges(3,4) * t457 + t356 * t46) * t270 + 0.4e1 * t46 * t206 + t79 * t505 + (0.2e1 * t67 + t386) * t253) * t97 * t262) * t271 + Ifges(3,4) * t97 * t386 + (t34 * t49 * t262 + (-mrSges(3,2) * t498 - t43) * t508) * t270 + (t31 * t457 + (Ifges(2,6) * t79 + t394 * t97) * t79) * t262 + (-Ifges(3,4) * t67 - t372 * t79 + t73) * t508 + (t205 + (-t189 - t312) * t166) * t272 + (g(3) * t312 + t166 * t234 + t184) * t263;
t19 = [-g(1) * m(4) + t1 * t341 + t2 * t339 + t3 * t337 + t9 * t328 + t8 * t329 + t7 * t330 + t6 * t346 + t5 * t348 + t4 * t350 + (-t42 * t338 + (t119 * t60 + t443 * t66) * t422 - t41 * t340 + (t117 * t59 + t444 * t65) * t423 - t40 * t342 + (t115 * t58 + t445 * t64) * t424) * t259 + (-t42 * t435 + (-t343 * t66 + t490 * t60) * t422 - t41 * t437 + (-t344 * t65 + t491 * t59) * t423 - t40 * t439 + (-t345 * t64 + t492 * t58) * t424) * t258 + (t42 * t337 + (t120 * t60 + t446 * t66) * t422 + t41 * t339 + (t118 * t59 + t447 * t65) * t423 + t40 * t341 + (t116 * t58 + t448 * t64) * t424 + m(4)) * t260; -g(2) * m(4) - t1 * t342 - t2 * t340 - t3 * t338 + t9 * t325 + t8 * t326 + t7 * t327 + t6 * t347 + t5 * t349 + t4 * t351 + (t39 * t337 + (t120 * t57 + t446 * t63) * t422 + t38 * t339 + (t118 * t56 + t447 * t62) * t423 + t37 * t341 + (t116 * t55 + t448 * t61) * t424) * t260 + (-t39 * t435 + (-t343 * t63 + t490 * t57) * t422 - t38 * t437 + (-t344 * t62 + t491 * t56) * t423 - t37 * t439 + (-t345 * t61 + t492 * t55) * t424) * t258 + (-t39 * t338 + (t119 * t57 + t443 * t63) * t422 - t38 * t340 + (t117 * t56 + t444 * t62) * t423 - t37 * t342 + (t115 * t55 + t445 * t61) * t424 + m(4)) * t259; -g(3) * m(4) - t1 * t439 - t2 * t437 - t3 * t435 - t9 * t313 - t8 * t314 - t7 * t315 + t6 * t322 + t5 * t323 + t4 * t324 + (t54 * t337 + (t120 * t84 + t446 * t87) * t422 + t53 * t339 + (t118 * t83 + t447 * t86) * t423 + t52 * t341 + (t116 * t82 + t448 * t85) * t424) * t260 + (-t54 * t338 + (t119 * t84 + t443 * t87) * t422 - t53 * t340 + (t117 * t83 + t444 * t86) * t423 - t52 * t342 + (t115 * t82 + t445 * t85) * t424) * t259 + (-t54 * t435 + (-t343 * t87 + t490 * t84) * t422 - t53 * t437 + (-t344 * t86 + t491 * t83) * t423 - t52 * t439 + (-t345 * t85 + t492 * t82) * t424 + m(4)) * t258;];
tauX  = t19;
