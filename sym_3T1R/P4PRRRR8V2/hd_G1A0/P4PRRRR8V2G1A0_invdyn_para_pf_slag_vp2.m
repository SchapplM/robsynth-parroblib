% Calculate vector of inverse dynamics forces for parallel robot
% P4PRRRR8V2G1A0
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
% Datum: 2020-08-07 11:13
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P4PRRRR8V2G1A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G1A0_invdyn_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR8V2G1A0_invdyn_para_pf_slag_vp2: xDP has to be [4x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [4 1]), ...
  'P4PRRRR8V2G1A0_invdyn_para_pf_slag_vp2: xDDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G1A0_invdyn_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR8V2G1A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G1A0_invdyn_para_pf_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V2G1A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRRR8V2G1A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P4PRRRR8V2G1A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G1A0_invdyn_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G1A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:09:53
% EndTime: 2020-08-07 11:10:05
% DurationCPUTime: 11.98s
% Computational Cost: add. (58017->689), mult. (110081->1178), div. (5368->13), fcn. (112126->30), ass. (0->453)
t320 = cos(qJ(2,4));
t344 = pkin(7) + pkin(6);
t272 = t320 * t344;
t318 = sin(qJ(2,4));
t243 = pkin(2) * t318 - t272;
t310 = sin(pkin(4));
t319 = cos(qJ(3,4));
t312 = cos(pkin(4));
t317 = sin(qJ(3,4));
t445 = t312 * t317;
t460 = t310 * t318;
t302 = t319 ^ 2;
t518 = pkin(3) * t302;
t137 = 0.1e1 / ((pkin(3) * t445 + t243 * t310) * t319 + pkin(2) * t445 + t460 * t518);
t333 = cos(qJ(2,3));
t276 = t333 * t344;
t327 = sin(qJ(2,3));
t245 = pkin(2) * t327 - t276;
t332 = cos(qJ(3,3));
t326 = sin(qJ(3,3));
t443 = t312 * t326;
t456 = t310 * t327;
t305 = t332 ^ 2;
t517 = pkin(3) * t305;
t138 = 0.1e1 / ((pkin(3) * t443 + t245 * t310) * t332 + pkin(2) * t443 + t456 * t517);
t335 = cos(qJ(2,2));
t277 = t335 * t344;
t329 = sin(qJ(2,2));
t246 = pkin(2) * t329 - t277;
t334 = cos(qJ(3,2));
t328 = sin(qJ(3,2));
t441 = t312 * t328;
t454 = t310 * t329;
t306 = t334 ^ 2;
t516 = pkin(3) * t306;
t139 = 0.1e1 / ((pkin(3) * t441 + t246 * t310) * t334 + pkin(2) * t441 + t454 * t516);
t337 = cos(qJ(2,1));
t278 = t337 * t344;
t331 = sin(qJ(2,1));
t247 = pkin(2) * t331 - t278;
t336 = cos(qJ(3,1));
t330 = sin(qJ(3,1));
t439 = t312 * t330;
t452 = t310 * t331;
t307 = t336 ^ 2;
t515 = pkin(3) * t307;
t140 = 0.1e1 / ((pkin(3) * t439 + t247 * t310) * t336 + pkin(2) * t439 + t452 * t515);
t313 = legFrame(4,3);
t282 = sin(t313);
t286 = cos(t313);
t309 = sin(pkin(8));
t311 = cos(pkin(8));
t185 = -t282 * t309 + t286 * t311;
t189 = t282 * t311 + t286 * t309;
t271 = pkin(3) * t319 + pkin(2);
t207 = t271 * t318 - t272;
t434 = t318 * t344;
t469 = (t271 * t320 + t434) * t312;
t129 = -t185 * t469 + t189 * t207;
t130 = -t185 * t207 - t189 * t469;
t346 = xP(4);
t300 = sin(t346);
t301 = cos(t346);
t349 = koppelP(4,2);
t353 = koppelP(4,1);
t223 = -t300 * t349 + t301 * t353;
t341 = xDP(4);
t342 = xDP(2);
t173 = t223 * t341 + t342;
t219 = t300 * t353 + t301 * t349;
t343 = xDP(1);
t177 = -t219 * t341 + t343;
t239 = t271 * t445;
t459 = t310 * t319;
t153 = 0.1e1 / (t207 * t459 + t239);
t359 = 0.1e1 / pkin(3);
t477 = t153 * t359;
t102 = (t129 * t177 + t130 * t173) * t477;
t545 = -0.2e1 * t102;
t314 = legFrame(3,3);
t283 = sin(t314);
t287 = cos(t314);
t186 = -t283 * t309 + t287 * t311;
t190 = t283 * t311 + t287 * t309;
t273 = pkin(3) * t332 + pkin(2);
t214 = t273 * t327 - t276;
t429 = t327 * t344;
t468 = (t273 * t333 + t429) * t312;
t131 = -t186 * t468 + t190 * t214;
t134 = -t186 * t214 - t190 * t468;
t350 = koppelP(3,2);
t354 = koppelP(3,1);
t224 = -t300 * t350 + t301 * t354;
t174 = t224 * t341 + t342;
t220 = t300 * t354 + t301 * t350;
t178 = -t220 * t341 + t343;
t240 = t273 * t443;
t451 = t310 * t332;
t155 = 0.1e1 / (t214 * t451 + t240);
t476 = t155 * t359;
t106 = (t131 * t178 + t134 * t174) * t476;
t544 = -0.2e1 * t106;
t315 = legFrame(2,3);
t284 = sin(t315);
t288 = cos(t315);
t187 = -t284 * t309 + t288 * t311;
t191 = t284 * t311 + t288 * t309;
t274 = pkin(3) * t334 + pkin(2);
t215 = t274 * t329 - t277;
t427 = t329 * t344;
t467 = (t274 * t335 + t427) * t312;
t132 = -t187 * t467 + t191 * t215;
t135 = -t187 * t215 - t191 * t467;
t351 = koppelP(2,2);
t355 = koppelP(2,1);
t225 = -t300 * t351 + t301 * t355;
t175 = t225 * t341 + t342;
t221 = t300 * t355 + t301 * t351;
t179 = -t221 * t341 + t343;
t241 = t274 * t441;
t449 = t310 * t334;
t156 = 0.1e1 / (t215 * t449 + t241);
t475 = t156 * t359;
t107 = (t132 * t179 + t135 * t175) * t475;
t543 = -0.2e1 * t107;
t316 = legFrame(1,3);
t285 = sin(t316);
t289 = cos(t316);
t188 = -t285 * t309 + t289 * t311;
t192 = t285 * t311 + t289 * t309;
t275 = pkin(3) * t336 + pkin(2);
t216 = t275 * t331 - t278;
t425 = t331 * t344;
t466 = (t275 * t337 + t425) * t312;
t133 = -t188 * t466 + t192 * t216;
t136 = -t188 * t216 - t192 * t466;
t352 = koppelP(1,2);
t356 = koppelP(1,1);
t226 = -t300 * t352 + t301 * t356;
t176 = t226 * t341 + t342;
t222 = t300 * t356 + t301 * t352;
t180 = -t222 * t341 + t343;
t242 = t275 * t439;
t447 = t310 * t336;
t157 = 0.1e1 / (t216 * t447 + t242);
t474 = t157 * t359;
t108 = (t133 * t180 + t136 * t176) * t474;
t542 = -0.2e1 * t108;
t512 = g(3) * t312;
t385 = t319 * mrSges(3,1) - mrSges(3,2) * t317;
t384 = t332 * mrSges(3,1) - mrSges(3,2) * t326;
t383 = t334 * mrSges(3,1) - mrSges(3,2) * t328;
t382 = t336 * mrSges(3,1) - mrSges(3,2) * t330;
t250 = pkin(2) * t337 + t425;
t453 = t310 * t330;
t366 = pkin(3) * t453 - t247 * t312;
t537 = t250 * t311 + t366 * t309;
t249 = pkin(2) * t335 + t427;
t455 = t310 * t328;
t367 = pkin(3) * t455 - t246 * t312;
t536 = t249 * t311 + t367 * t309;
t248 = pkin(2) * t333 + t429;
t457 = t310 * t326;
t368 = pkin(3) * t457 - t245 * t312;
t535 = t248 * t311 + t368 * t309;
t244 = pkin(2) * t320 + t434;
t461 = t310 * t317;
t369 = pkin(3) * t461 - t243 * t312;
t534 = t244 * t311 + t369 * t309;
t234 = -g(1) * t285 + g(2) * t289;
t238 = g(1) * t289 + g(2) * t285;
t371 = t234 * t311 - t238 * t309;
t513 = g(3) * t310;
t533 = t371 * t312 + t513;
t233 = -g(1) * t284 + g(2) * t288;
t237 = g(1) * t288 + g(2) * t284;
t373 = t233 * t311 - t237 * t309;
t532 = t373 * t312 + t513;
t232 = -g(1) * t283 + g(2) * t287;
t236 = g(1) * t287 + g(2) * t283;
t375 = t232 * t311 - t236 * t309;
t531 = t375 * t312 + t513;
t231 = -g(1) * t282 + g(2) * t286;
t235 = g(1) * t286 + g(2) * t282;
t377 = t231 * t311 - t235 * t309;
t530 = t377 * t312 + t513;
t321 = Ifges(3,1) - Ifges(3,2);
t338 = pkin(2) * mrSges(3,2);
t525 = -2 * Ifges(3,4);
t529 = Ifges(3,4) + t307 * t525 + (-t321 * t330 + t338) * t336;
t528 = Ifges(3,4) + t306 * t525 + (-t321 * t328 + t338) * t334;
t527 = Ifges(3,4) + t305 * t525 + (-t321 * t326 + t338) * t332;
t526 = Ifges(3,4) + t302 * t525 + (-t317 * t321 + t338) * t319;
t340 = mrSges(3,1) * pkin(2);
t290 = mrSges(3,2) * pkin(6) - Ifges(3,6);
t524 = -t290 / 0.2e1;
t291 = mrSges(3,1) * pkin(6) - Ifges(3,5);
t523 = t291 / 0.2e1;
t522 = pkin(3) * t102;
t521 = pkin(3) * t106;
t520 = pkin(3) * t107;
t519 = pkin(3) * t108;
t303 = m(1) + m(2) + m(3);
t514 = g(3) * t303;
t360 = pkin(2) ^ 2;
t279 = t344 ^ 2 + t360;
t358 = pkin(3) ^ 2;
t405 = t317 * t522;
t423 = 0.2e1 * pkin(2) * pkin(3);
t444 = t312 * t318;
t121 = -t185 * t459 - (t185 * t444 + t189 * t320) * t317;
t122 = -t189 * t459 - (-t185 * t320 + t189 * t444) * t317;
t90 = (t121 * t177 + t122 * t173) * t137;
t511 = (-t344 * t405 + (t302 * t358 + t319 * t423 + t279) * t90) * t90;
t404 = t326 * t521;
t442 = t312 * t327;
t123 = -t186 * t451 - (t186 * t442 + t190 * t333) * t326;
t126 = -t190 * t451 - (-t186 * t333 + t190 * t442) * t326;
t94 = (t123 * t178 + t126 * t174) * t138;
t510 = (-t344 * t404 + (t305 * t358 + t332 * t423 + t279) * t94) * t94;
t403 = t328 * t520;
t440 = t312 * t329;
t124 = -t187 * t449 - (t187 * t440 + t191 * t335) * t328;
t127 = -t191 * t449 - (-t187 * t335 + t191 * t440) * t328;
t95 = (t124 * t179 + t127 * t175) * t139;
t509 = (-t344 * t403 + (t306 * t358 + t334 * t423 + t279) * t95) * t95;
t402 = t330 * t519;
t438 = t312 * t331;
t125 = -t188 * t447 - (t188 * t438 + t192 * t337) * t330;
t128 = -t192 * t447 - (-t188 * t337 + t192 * t438) * t330;
t96 = (t125 * t180 + t128 * t176) * t140;
t508 = (-t344 * t402 + (t307 * t358 + t336 * t423 + t279) * t96) * t96;
t506 = mrSges(3,1) * t317;
t505 = mrSges(3,1) * t326;
t504 = mrSges(3,1) * t328;
t503 = mrSges(3,1) * t330;
t502 = mrSges(3,2) * t310;
t497 = t320 * t90;
t496 = t333 * t94;
t495 = t335 * t95;
t494 = t337 * t96;
t493 = t344 * t90;
t492 = t344 * t94;
t491 = t344 * t95;
t490 = t344 * t96;
t252 = mrSges(3,2) * t319 + t506;
t98 = t102 ^ 2;
t489 = t98 * t252;
t296 = m(3) * pkin(2) + mrSges(2,1);
t103 = t106 ^ 2;
t256 = mrSges(3,2) * t332 + t505;
t488 = t103 * t256;
t104 = t107 ^ 2;
t257 = mrSges(3,2) * t334 + t504;
t487 = t104 * t257;
t105 = t108 ^ 2;
t258 = mrSges(3,2) * t336 + t503;
t486 = t105 * t258;
t308 = t341 ^ 2;
t322 = xDDP(4);
t325 = xDDP(1);
t149 = -t219 * t322 - t223 * t308 + t325;
t485 = t129 * t149;
t324 = xDDP(2);
t145 = -t219 * t308 + t223 * t322 + t324;
t484 = t130 * t145;
t150 = -t220 * t322 - t224 * t308 + t325;
t483 = t131 * t150;
t151 = -t221 * t322 - t225 * t308 + t325;
t482 = t132 * t151;
t152 = -t222 * t322 - t226 * t308 + t325;
t481 = t133 * t152;
t146 = -t220 * t308 + t224 * t322 + t324;
t480 = t134 * t146;
t147 = -t221 * t308 + t225 * t322 + t324;
t479 = t135 * t147;
t148 = -t222 * t308 + t226 * t322 + t324;
t478 = t136 * t148;
t227 = t385 + t296;
t280 = m(3) * pkin(6) - mrSges(2,2) + mrSges(3,3);
t169 = t227 * t320 + t280 * t318;
t473 = t169 * t310;
t228 = t384 + t296;
t170 = t228 * t333 + t280 * t327;
t472 = t170 * t310;
t229 = t383 + t296;
t171 = t229 * t335 + t280 * t329;
t471 = t171 * t310;
t230 = t382 + t296;
t172 = t230 * t337 + t280 * t331;
t470 = t172 * t310;
t458 = t310 * t320;
t450 = t310 * t333;
t448 = t310 * t335;
t446 = t310 * t337;
t437 = t312 * t359;
t435 = t318 * t319;
t430 = t327 * t332;
t428 = t329 * t334;
t426 = t331 * t336;
t424 = mrSges(3,2) * t512;
t418 = -0.2e1 * t338;
t165 = -t252 * t460 + t385 * t312;
t181 = pkin(3) * t435 + t243;
t397 = t317 * t493;
t18 = t137 * t437 * t511 + (-t312 * t397 + (-t181 * t461 + t312 * (pkin(2) * t319 + t518)) * t102) / (t181 * t459 + t239) * t102;
t53 = t397 - t522;
t26 = (-t319 * t511 - (pkin(2) * t102 - t319 * t53) * t522) * t137;
t89 = t90 ^ 2;
t413 = -t165 * t18 - t303 * t26 + ((-t296 * t89 - t385 * (t89 + t98)) * t318 + (t252 * t545 + t280 * t90) * t497) * t310;
t166 = -t256 * t456 + t384 * t312;
t182 = pkin(3) * t430 + t245;
t396 = t326 * t492;
t22 = t138 * t437 * t510 + (-t312 * t396 + (-t182 * t457 + t312 * (pkin(2) * t332 + t517)) * t106) / (t182 * t451 + t240) * t106;
t55 = t396 - t521;
t30 = (-t332 * t510 - (pkin(2) * t106 - t332 * t55) * t521) * t138;
t91 = t94 ^ 2;
t412 = -t166 * t22 - t303 * t30 + ((-t296 * t91 - t384 * (t91 + t103)) * t327 + (t256 * t544 + t280 * t94) * t496) * t310;
t167 = -t257 * t454 + t383 * t312;
t183 = pkin(3) * t428 + t246;
t395 = t328 * t491;
t23 = t139 * t437 * t509 + (-t312 * t395 + (-t183 * t455 + t312 * (pkin(2) * t334 + t516)) * t107) / (t183 * t449 + t241) * t107;
t56 = t395 - t520;
t31 = (-t334 * t509 - (pkin(2) * t107 - t334 * t56) * t520) * t139;
t92 = t95 ^ 2;
t411 = -t167 * t23 - t303 * t31 + ((-t296 * t92 - t383 * (t92 + t104)) * t329 + (t257 * t543 + t280 * t95) * t495) * t310;
t168 = -t258 * t452 + t382 * t312;
t184 = pkin(3) * t426 + t247;
t394 = t330 * t490;
t24 = t140 * t437 * t508 + (-t312 * t394 + (-t184 * t453 + t312 * (pkin(2) * t336 + t515)) * t108) / (t184 * t447 + t242) * t108;
t57 = t394 - t519;
t32 = (-t336 * t508 - (pkin(2) * t108 - t336 * t57) * t519) * t140;
t93 = t96 ^ 2;
t410 = -t168 * t24 - t303 * t32 + ((-t296 * t93 - t382 * (t93 + t105)) * t331 + (t258 * t542 + t280 * t96) * t494) * t310;
t409 = pkin(2) * t461;
t408 = pkin(2) * t457;
t407 = pkin(2) * t455;
t406 = pkin(2) * t453;
t401 = Ifges(3,3) * t477;
t400 = Ifges(3,3) * t476;
t399 = Ifges(3,3) * t475;
t398 = Ifges(3,3) * t474;
t393 = t165 * t477;
t195 = -t290 * t319 - t317 * t291;
t392 = t195 * t477;
t391 = t166 * t476;
t203 = -t290 * t332 - t326 * t291;
t390 = t203 * t476;
t389 = t167 * t475;
t204 = -t290 * t334 - t328 * t291;
t388 = t204 * t475;
t387 = t168 * t474;
t205 = -t290 * t336 - t330 * t291;
t386 = t205 * t474;
t376 = t231 * t309 + t235 * t311;
t374 = t232 * t309 + t236 * t311;
t372 = t233 * t309 + t237 * t311;
t370 = t234 * t309 + t238 * t311;
t365 = Ifges(3,1) + Ifges(2,3) + (pkin(6) ^ 2 + t360) * m(3) + 0.2e1 * pkin(6) * mrSges(3,3);
t364 = t530 * t318 + t376 * t320;
t363 = t531 * t327 + t374 * t333;
t362 = t532 * t329 + t372 * t335;
t361 = t533 * t331 + t370 * t337;
t348 = mrSges(4,1);
t347 = mrSges(4,2);
t323 = xDDP(3);
t218 = -t347 * t300 + t301 * t348;
t217 = t348 * t300 + t301 * t347;
t201 = t309 * t337 + t311 * t438;
t200 = t309 * t335 + t311 * t440;
t199 = t309 * t333 + t311 * t442;
t198 = t309 * t438 - t311 * t337;
t197 = t309 * t440 - t311 * t335;
t196 = t309 * t442 - t311 * t333;
t194 = t309 * t320 + t311 * t444;
t193 = t309 * t444 - t311 * t320;
t160 = -t321 * t307 + 0.2e1 * (Ifges(3,4) * t330 + t340) * t336 + t330 * t418 + t365;
t159 = -t321 * t306 + 0.2e1 * (Ifges(3,4) * t328 + t340) * t334 + t328 * t418 + t365;
t158 = -t321 * t305 + 0.2e1 * (Ifges(3,4) * t326 + t340) * t332 + t326 * t418 + t365;
t154 = -t321 * t302 + 0.2e1 * (Ifges(3,4) * t317 + t340) * t319 + t317 * t418 + t365;
t144 = t250 * t309 - t366 * t311;
t143 = t249 * t309 - t367 * t311;
t142 = t248 * t309 - t368 * t311;
t141 = t244 * t309 - t369 * t311;
t120 = (-t198 * t285 + t201 * t289) * t515 + (t144 * t289 + t537 * t285) * t336 - t188 * t406;
t119 = (-t197 * t284 + t200 * t288) * t516 + (t143 * t288 + t536 * t284) * t334 - t187 * t407;
t118 = (-t196 * t283 + t199 * t287) * t517 + (t142 * t287 + t535 * t283) * t332 - t186 * t408;
t117 = -(t198 * t289 + t201 * t285) * t515 + (-t285 * t144 + t537 * t289) * t336 + t192 * t406;
t116 = -(t197 * t288 + t200 * t284) * t516 + (-t284 * t143 + t536 * t288) * t334 + t191 * t407;
t115 = -(t196 * t287 + t199 * t283) * t517 + (-t283 * t142 + t535 * t287) * t332 + t190 * t408;
t114 = (-t193 * t282 + t194 * t286) * t518 + (t141 * t286 + t534 * t282) * t319 - t185 * t409;
t113 = -(t193 * t286 + t194 * t282) * t518 + (-t282 * t141 + t534 * t286) * t319 + t189 * t409;
t112 = (-t133 * t222 + t136 * t226) * t474;
t111 = (-t132 * t221 + t135 * t225) * t475;
t110 = (-t131 * t220 + t134 * t224) * t476;
t109 = (-t129 * t219 + t130 * t223) * t477;
t101 = (-t125 * t222 + t128 * t226) * t140;
t100 = (-t124 * t221 + t127 * t225) * t139;
t99 = (-t123 * t220 + t126 * t224) * t138;
t97 = (-t121 * t219 + t122 * t223) * t137;
t88 = (-t117 * t222 + t120 * t226) * t140;
t87 = (-t116 * t221 + t119 * t225) * t139;
t86 = (-t115 * t220 + t118 * t224) * t138;
t85 = (-t113 * t219 + t114 * t223) * t137;
t84 = t136 * t398 + (t120 * t168 + t128 * t205) * t140;
t83 = t135 * t399 + (t119 * t167 + t127 * t204) * t139;
t82 = t134 * t400 + (t118 * t166 + t126 * t203) * t138;
t81 = t133 * t398 + (t117 * t168 + t125 * t205) * t140;
t80 = t132 * t399 + (t116 * t167 + t124 * t204) * t139;
t79 = t131 * t400 + (t115 * t166 + t123 * t203) * t138;
t78 = t130 * t401 + (t114 * t165 + t122 * t195) * t137;
t77 = t129 * t401 + (t113 * t165 + t121 * t195) * t137;
t76 = t136 * t387 + (t120 * t303 + t128 * t470) * t140;
t75 = t135 * t389 + (t119 * t303 + t127 * t471) * t139;
t74 = t134 * t391 + (t118 * t303 + t126 * t472) * t138;
t73 = t133 * t387 + (t117 * t303 + t125 * t470) * t140;
t72 = t132 * t389 + (t116 * t303 + t124 * t471) * t139;
t71 = t131 * t391 + (t115 * t303 + t123 * t472) * t138;
t70 = t130 * t393 + (t114 * t303 + t122 * t473) * t137;
t69 = t129 * t393 + (t113 * t303 + t121 * t473) * t137;
t68 = t136 * t386 + (t120 * t470 + t128 * t160) * t140;
t67 = t135 * t388 + (t119 * t471 + t127 * t159) * t139;
t66 = t134 * t390 + (t118 * t472 + t126 * t158) * t138;
t65 = t133 * t386 + (t117 * t470 + t125 * t160) * t140;
t64 = t132 * t388 + (t116 * t471 + t124 * t159) * t139;
t63 = t131 * t390 + (t115 * t472 + t123 * t158) * t138;
t62 = t130 * t392 + (t114 * t473 + t122 * t154) * t137;
t61 = t129 * t392 + (t113 * t473 + t121 * t154) * t137;
t49 = t101 * t470 + t112 * t168 + t303 * t88;
t48 = t100 * t471 + t111 * t167 + t303 * t87;
t47 = t110 * t166 + t303 * t86 + t99 * t472;
t45 = t101 * t160 + t112 * t205 + t88 * t470;
t44 = t100 * t159 + t111 * t204 + t87 * t471;
t43 = t110 * t203 + t158 * t99 + t86 * t472;
t42 = t109 * t165 + t303 * t85 + t97 * t473;
t41 = t109 * t195 + t154 * t97 + t85 * t473;
t16 = (((t108 * t312 + t96 * t446) * t515 + ((-t402 + t490) * t331 + pkin(2) * t494) * t447 + t312 * t57) * t96 + (t108 * t446 + (t307 * t312 - t426 * t453 - t312) * t96) * t519) * t140;
t15 = (((t107 * t312 + t95 * t448) * t516 + ((-t403 + t491) * t329 + pkin(2) * t495) * t449 + t312 * t56) * t95 + (t107 * t448 + (t306 * t312 - t428 * t455 - t312) * t95) * t520) * t139;
t14 = (((t106 * t312 + t94 * t450) * t517 + ((-t404 + t492) * t327 + pkin(2) * t496) * t451 + t312 * t55) * t94 + (t106 * t450 + (t305 * t312 - t430 * t457 - t312) * t94) * t521) * t138;
t13 = (((t102 * t312 + t90 * t458) * t518 + ((-t405 + t493) * t318 + pkin(2) * t497) * t459 + t53 * t312) * t90 + (t102 * t458 + (t302 * t312 - t435 * t461 - t312) * t90) * t522) * t137;
t12 = -t168 * t32 - t205 * t16 - Ifges(3,3) * t24 + t93 * (pkin(2) * t503 + t529) + ((t371 * t310 - t512) * mrSges(3,1) + t361 * mrSges(3,2)) * t336 + (t361 * mrSges(3,1) - t371 * t502 + t424) * t330;
t11 = -t167 * t31 - t204 * t15 - Ifges(3,3) * t23 + t92 * (pkin(2) * t504 + t528) + ((t373 * t310 - t512) * mrSges(3,1) + t362 * mrSges(3,2)) * t334 + (t362 * mrSges(3,1) - t373 * t502 + t424) * t328;
t10 = -t166 * t30 - t203 * t14 - Ifges(3,3) * t22 + t91 * (pkin(2) * t505 + t527) + ((t375 * t310 - t512) * mrSges(3,1) + t363 * mrSges(3,2)) * t332 + (t363 * mrSges(3,1) - t375 * t502 + t424) * t326;
t9 = -t165 * t26 - t195 * t13 - Ifges(3,3) * t18 + t89 * (pkin(2) * t506 + t526) + ((t377 * t310 - t512) * mrSges(3,1) + t364 * mrSges(3,2)) * t319 + (t364 * mrSges(3,1) - t377 * t502 + t424) * t317;
t8 = -t32 * t470 - t160 * t16 - t205 * t24 + ((t524 * t330 + t523 * t336) * t108 + (t340 * t330 + t529) * t96) * t542 + (-t230 * t533 - t370 * t280) * t337 + (t370 * t230 - t280 * t533) * t331;
t7 = -t31 * t471 - t159 * t15 - t204 * t23 + ((t524 * t328 + t523 * t334) * t107 + (t340 * t328 + t528) * t95) * t543 + (-t229 * t532 - t372 * t280) * t335 + (t372 * t229 - t280 * t532) * t329;
t6 = -t30 * t472 - t158 * t14 - t203 * t22 + ((t524 * t326 + t523 * t332) * t106 + (t340 * t326 + t527) * t94) * t544 + (-t228 * t531 - t374 * t280) * t333 + (t374 * t228 - t280 * t531) * t327;
t5 = -t26 * t473 - t154 * t13 - t195 * t18 + ((t524 * t317 + t523 * t319) * t102 + (t340 * t317 + t526) * t90) * t545 + (-t227 * t530 - t376 * t280) * t320 + (t376 * t227 - t280 * t530) * t318;
t4 = -t16 * t470 - t312 * t486 + t410 - t514;
t3 = -t15 * t471 - t312 * t487 + t411 - t514;
t2 = -t14 * t472 - t312 * t488 + t412 - t514;
t1 = -t13 * t473 - t312 * t489 + t413 - t514;
t17 = [-t217 * t322 - t308 * t218 + (t325 - g(1)) * m(4) + (t69 + t71 + t72 + t73) * t323 + ((t120 * t73 + t128 * t65) * t148 + (t117 * t73 + t125 * t65) * t152 + t117 * t4 + t125 * t8) * t140 + ((t119 * t72 + t127 * t64) * t147 + (t116 * t72 + t124 * t64) * t151 + t116 * t3 + t124 * t7) * t139 + ((t115 * t71 + t123 * t63) * t150 + (t118 * t71 + t126 * t63) * t146 + t115 * t2 + t123 * t6) * t138 + ((t81 * t478 + (t152 * t81 + t12) * t133) * t157 + (t80 * t479 + (t151 * t80 + t11) * t132) * t156 + (t79 * t480 + (t150 * t79 + t10) * t131) * t155 + (t77 * t484 + (t149 * t77 + t9) * t129) * t153) * t359 + ((t113 * t69 + t121 * t61) * t149 + (t114 * t69 + t122 * t61) * t145 + t113 * t1 + t121 * t5) * t137; -t308 * t217 + t218 * t322 + (t324 - g(2)) * m(4) + (t70 + t74 + t75 + t76) * t323 + ((t120 * t76 + t128 * t68) * t148 + (t117 * t76 + t125 * t68) * t152 + t120 * t4 + t128 * t8) * t140 + ((t119 * t75 + t127 * t67) * t147 + (t116 * t75 + t124 * t67) * t151 + t119 * t3 + t127 * t7) * t139 + ((t115 * t74 + t123 * t66) * t150 + (t118 * t74 + t126 * t66) * t146 + t118 * t2 + t126 * t6) * t138 + ((t84 * t481 + (t148 * t84 + t12) * t136) * t157 + (t83 * t482 + (t147 * t83 + t11) * t135) * t156 + (t82 * t483 + (t146 * t82 + t10) * t134) * t155 + (t78 * t485 + (t145 * t78 + t9) * t130) * t153) * t359 + ((t113 * t70 + t121 * t62) * t149 + (t114 * t70 + t122 * t62) * t145 + t114 * t1 + t122 * t5) * t137; t411 + (-t169 * t13 - t170 * t14 - t171 * t15 - t172 * t16) * t310 + (-t486 - t487 - t488 - t489) * t312 + t76 * t148 + t69 * t149 + t71 * t150 + t72 * t151 + t73 * t152 + t70 * t145 + t74 * t146 + t75 * t147 + t410 + t412 + t413 + (t323 - g(3)) * (0.4e1 * t303 + m(4)); t109 * t9 + t110 * t10 + t111 * t11 + t112 * t12 + t97 * t5 + t99 * t6 + t100 * t7 + t101 * t8 - (-g(1) * t348 - g(2) * t347) * t300 + t301 * (g(1) * t347 - g(2) * t348) + t85 * t1 + t86 * t2 + t87 * t3 + t88 * t4 + Ifges(4,3) * t322 + t218 * t324 - t217 * t325 + ((t120 * t49 + t128 * t45) * t148 + (t117 * t49 + t125 * t45) * t152) * t140 + ((t116 * t48 + t124 * t44) * t151 + (t119 * t48 + t127 * t44) * t147) * t139 + ((t478 + t481) * (Ifges(3,3) * t112 + t101 * t205 + t168 * t88) * t157 + (t479 + t482) * (Ifges(3,3) * t111 + t100 * t204 + t167 * t87) * t156 + (t480 + t483) * (Ifges(3,3) * t110 + t166 * t86 + t203 * t99) * t155 + (t484 + t485) * (Ifges(3,3) * t109 + t165 * t85 + t195 * t97) * t153) * t359 + ((t115 * t47 + t123 * t43) * t150 + (t118 * t47 + t126 * t43) * t146) * t138 + ((t113 * t42 + t121 * t41) * t149 + (t114 * t42 + t122 * t41) * t145) * t137 + (t42 + t47 + t48 + t49) * t323;];
tauX  = t17;
