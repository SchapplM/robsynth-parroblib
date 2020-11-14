% Calculate vector of inverse dynamics forces for parallel robot
% P3RPRRR6V1G2A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:37
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRRR6V1G2A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G2A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:35:36
% EndTime: 2020-08-06 18:35:49
% DurationCPUTime: 12.88s
% Computational Cost: add. (51516->700), mult. (58989->1080), div. (8619->16), fcn. (39375->106), ass. (0->457)
t250 = sin(pkin(7));
t281 = -pkin(6) - pkin(5);
t440 = t250 * t281;
t143 = -pkin(1) + t440;
t262 = sin(qJ(1,3));
t123 = t143 * t262;
t251 = cos(pkin(7));
t268 = cos(qJ(1,3));
t438 = t268 * t281;
t443 = t250 * t268;
t562 = t123 - pkin(2) * t443 - (pkin(2) * t262 + t438) * t251;
t264 = sin(qJ(1,2));
t124 = t143 * t264;
t270 = cos(qJ(1,2));
t437 = t270 * t281;
t442 = t250 * t270;
t561 = t124 - pkin(2) * t442 - (pkin(2) * t264 + t437) * t251;
t266 = sin(qJ(1,1));
t125 = t143 * t266;
t272 = cos(qJ(1,1));
t436 = t272 * t281;
t441 = t250 * t272;
t560 = t125 - pkin(2) * t441 - (pkin(2) * t266 + t436) * t251;
t259 = xDDP(2);
t559 = t259 / 0.2e1;
t260 = xDDP(1);
t558 = t260 / 0.2e1;
t267 = cos(qJ(3,3));
t197 = t267 * pkin(3);
t164 = t197 + pkin(2);
t190 = t251 * pkin(1);
t126 = 0.1e1 / (t190 + t164);
t261 = sin(qJ(3,3));
t243 = 0.1e1 / t261;
t467 = t126 * t243;
t464 = t164 * t262;
t72 = (t438 + t464) * t251 - t123 + t164 * t443;
t386 = t72 * t467;
t278 = xDP(3);
t209 = pkin(7) + qJ(3,3);
t171 = sin(t209);
t213 = -pkin(7) + qJ(3,3);
t173 = sin(t213);
t220 = qJ(1,3) + pkin(7);
t175 = sin(t220);
t186 = cos(t220);
t198 = t268 * pkin(1);
t294 = 0.2e1 * qJ(3,3);
t227 = sin(t294);
t517 = 0.2e1 * t281;
t527 = -0.2e1 * pkin(2);
t550 = 0.2e1 * pkin(2);
t497 = (t175 * t517 + t186 * t527 - 0.2e1 * t198 + (-cos(qJ(1,3) - t213) - cos(qJ(1,3) + t209)) * pkin(3)) / (t261 * t550 + pkin(3) * t227 + (t171 + t173) * pkin(1));
t390 = t278 * t497;
t305 = pkin(3) ^ 2;
t306 = 0.1e1 / pkin(3);
t429 = t305 * t306;
t255 = legFrame(3,2);
t194 = cos(t255);
t280 = xDP(1);
t450 = t194 * t280;
t191 = sin(t255);
t279 = xDP(2);
t455 = t191 * t279;
t557 = t429 * (t390 / 0.4e1 - (t450 / 0.4e1 - t455 / 0.4e1) * t386);
t269 = cos(qJ(3,2));
t199 = t269 * pkin(3);
t165 = t199 + pkin(2);
t127 = 0.1e1 / (t190 + t165);
t263 = sin(qJ(3,2));
t244 = 0.1e1 / t263;
t466 = t127 * t244;
t463 = t165 * t264;
t74 = (t437 + t463) * t251 - t124 + t165 * t442;
t385 = t74 * t466;
t211 = pkin(7) + qJ(3,2);
t172 = sin(t211);
t215 = -pkin(7) + qJ(3,2);
t174 = sin(t215);
t223 = qJ(1,2) + pkin(7);
t176 = sin(t223);
t187 = cos(t223);
t200 = t270 * pkin(1);
t297 = 0.2e1 * qJ(3,2);
t230 = sin(t297);
t496 = (t176 * t517 + t187 * t527 - 0.2e1 * t200 + (-cos(qJ(1,2) - t215) - cos(qJ(1,2) + t211)) * pkin(3)) / (t263 * t550 + pkin(3) * t230 + (t172 + t174) * pkin(1));
t388 = t278 * t496;
t256 = legFrame(2,2);
t195 = cos(t256);
t448 = t195 * t280;
t192 = sin(t256);
t453 = t192 * t279;
t556 = t429 * (t388 / 0.4e1 - (t448 / 0.4e1 - t453 / 0.4e1) * t385);
t249 = t281 ^ 2;
t307 = pkin(2) ^ 2;
t284 = 0.3e1 * t307;
t518 = t284 + t249;
t308 = pkin(1) ^ 2;
t555 = t308 + t518;
t510 = Icges(3,2) / 0.2e1;
t554 = t510 - Icges(3,1) / 0.2e1;
t149 = t255 + t220;
t150 = -t255 + t220;
t106 = cos(t150) + cos(t149);
t258 = xDDP(3);
t461 = t175 * t258;
t553 = t106 * t558 - t461;
t151 = t256 + t223;
t152 = -t256 + t223;
t107 = cos(t152) + cos(t151);
t459 = t176 * t258;
t552 = t107 * t558 - t459;
t224 = qJ(1,1) + pkin(7);
t257 = legFrame(1,2);
t153 = t257 + t224;
t154 = -t257 + t224;
t108 = cos(t154) + cos(t153);
t177 = sin(t224);
t457 = t177 * t258;
t551 = t108 * t558 - t457;
t549 = 0.4e1 * pkin(3);
t515 = m(3) / 0.2e1;
t103 = -sin(t149) + sin(t150);
t479 = t103 / 0.2e1;
t104 = -sin(t151) + sin(t152);
t477 = t104 / 0.2e1;
t105 = -sin(t153) + sin(t154);
t475 = t105 / 0.2e1;
t473 = t106 / 0.2e1;
t471 = t107 / 0.2e1;
t469 = t108 / 0.2e1;
t435 = t279 / 0.2e1;
t434 = t280 / 0.2e1;
t544 = t105 * t559;
t543 = t104 * t559;
t542 = t103 * t559;
t483 = t306 * t74;
t375 = t244 * t483;
t541 = (t448 / 0.6e1 - t453 / 0.6e1) * t375;
t540 = (t448 / 0.3e1 - t453 / 0.3e1) * t375;
t539 = (t448 / 0.2e1 - t453 / 0.2e1) * t375;
t484 = t306 * t72;
t377 = t243 * t484;
t538 = (t450 / 0.6e1 - t455 / 0.6e1) * t377;
t537 = (t450 / 0.3e1 - t455 / 0.3e1) * t377;
t536 = (t450 / 0.2e1 - t455 / 0.2e1) * t377;
t503 = pkin(1) * t250;
t535 = -t281 + t503;
t514 = pkin(1) * pkin(2);
t534 = 0.8e1 * t281 * t514;
t533 = 0.4e1 * t281 * t308;
t290 = 0.2e1 * pkin(7);
t169 = t308 * cos(t290);
t532 = -t169 - t305 / 0.2e1;
t428 = t308 * sin(t290);
t398 = 0.2e1 * t428;
t410 = -0.4e1 * t190;
t531 = t281 * t410 + t398;
t528 = -0.4e1 * pkin(1);
t526 = 0.4e1 * pkin(2);
t41 = (t390 + (-t450 + t455) * t386) * t306;
t524 = t41 ^ 2;
t42 = (t388 + (-t448 + t453) * t385) * t306;
t523 = t42 ^ 2;
t193 = sin(t257);
t196 = cos(t257);
t271 = cos(qJ(3,1));
t201 = t271 * pkin(3);
t167 = t201 + pkin(2);
t128 = 0.1e1 / (t190 + t167);
t265 = sin(qJ(3,1));
t245 = 0.1e1 / t265;
t465 = t128 * t245;
t188 = cos(t224);
t202 = t272 * pkin(1);
t408 = pkin(7) + qJ(3,1);
t409 = -pkin(7) + qJ(3,1);
t300 = 0.2e1 * qJ(3,1);
t233 = sin(t300);
t499 = pkin(3) * t233;
t495 = (t177 * t517 + t188 * t527 - 0.2e1 * t202 + (-cos(qJ(1,1) - t409) - cos(qJ(1,1) + t408)) * pkin(3)) / (t265 * t550 + t499 + (sin(t408) + sin(t409)) * pkin(1));
t462 = t167 * t266;
t76 = (t436 + t462) * t251 - t125 + t167 * t441;
t43 = (t278 * t495 + (t193 * t279 - t196 * t280) * t76 * t465) * t306;
t522 = t43 ^ 2;
t141 = pkin(2) * t535;
t148 = t281 * t190;
t346 = -0.2e1 * t148 + t428;
t95 = 0.2e1 * t141 + t346;
t521 = -0.2e1 * t95;
t285 = 0.2e1 * t307;
t203 = t285 + t308;
t404 = pkin(2) * t190;
t109 = 0.4e1 * t404 + t169 + t203;
t520 = -0.2e1 * t305 - 0.2e1 * t109;
t519 = 0.2e1 * t109;
t516 = -0.6e1 * t305;
t289 = -0.2e1 * t308;
t513 = m(2) * rSges(2,1);
t512 = m(3) * rSges(3,2);
t511 = -Icges(3,5) / 0.4e1;
t509 = Icges(3,6) / 0.4e1;
t301 = rSges(3,2) ^ 2;
t302 = rSges(3,1) ^ 2;
t508 = (-t301 + t302) * t515 + t554;
t507 = t306 / 0.2e1;
t273 = rSges(3,3) + pkin(5);
t135 = rSges(3,1) * t267 - rSges(3,2) * t261;
t506 = m(3) * t135;
t136 = rSges(3,1) * t269 - rSges(3,2) * t263;
t505 = m(3) * t136;
t137 = rSges(3,1) * t271 - rSges(3,2) * t265;
t504 = m(3) * t137;
t502 = pkin(1) * t262;
t501 = pkin(1) * t264;
t500 = pkin(1) * t266;
t498 = pkin(3) * t281;
t88 = t103 * t126 * t435;
t90 = t106 * t126 * t434;
t494 = t88 + t90;
t89 = t104 * t127 * t435;
t91 = t107 * t127 * t434;
t493 = t89 + t91;
t492 = t261 * t41;
t460 = t175 * t278;
t363 = t126 * t460;
t69 = -t363 + t494;
t491 = t261 * t69;
t490 = t263 * t42;
t458 = t176 * t278;
t362 = t127 * t458;
t70 = -t362 + t493;
t489 = t263 * t70;
t488 = t265 * t43;
t73 = (t164 * t268 - t262 * t281) * t251 - t143 * t268 - t250 * t464;
t487 = t267 * t73;
t75 = (t165 * t270 - t264 * t281) * t251 - t143 * t270 - t250 * t463;
t486 = t269 * t75;
t77 = (t167 * t272 - t266 * t281) * t251 - t143 * t272 - t250 * t462;
t485 = t271 * t77;
t482 = t306 * t76;
t242 = cos(t300);
t481 = t95 * t242;
t468 = t109 * t242;
t456 = t191 * t261;
t454 = t192 * t263;
t452 = t193 * t265;
t451 = t194 * t261;
t449 = t195 * t263;
t447 = t196 * t265;
t446 = t243 * t267;
t445 = t244 * t269;
t444 = t245 * t271;
t439 = t258 * t306;
t292 = 0.4e1 * qJ(3,3);
t432 = t305 * cos(t292);
t295 = 0.4e1 * qJ(3,2);
t431 = t305 * cos(t295);
t298 = 0.4e1 * qJ(3,1);
t430 = t305 * cos(t298);
t427 = t141 - t148;
t407 = rSges(3,1) * t512;
t426 = -t407 / 0.2e1 + Icges(3,4) / 0.2e1;
t425 = t249 / 0.2e1 + t308;
t293 = 0.3e1 * qJ(3,3);
t235 = cos(t293);
t424 = t235 - t267;
t296 = 0.3e1 * qJ(3,2);
t238 = cos(t296);
t423 = t238 - t269;
t422 = t249 + t308;
t421 = t301 + t302;
t420 = 0.2e1 * pkin(1);
t418 = 0.2e1 * pkin(3);
t416 = -0.3e1 * pkin(1) * t305;
t415 = pkin(2) * t516;
t414 = -0.4e1 * pkin(2) * t308;
t413 = -0.8e1 * pkin(3) * (t305 / 0.4e1 + 0.3e1 / 0.2e1 * t307 + t425);
t155 = t190 + pkin(2);
t412 = -0.4e1 * pkin(3) * t155;
t411 = pkin(3) * t289;
t405 = pkin(1) * t498;
t129 = 0.3e1 / 0.8e1 * t305 + t307 / 0.2e1 + t425;
t403 = -0.16e2 * pkin(2) * t129;
t402 = pkin(3) * t492;
t401 = pkin(3) * t490;
t400 = pkin(3) * t488;
t399 = t305 * t517;
t166 = t201 + t550;
t397 = t166 * t190;
t396 = pkin(1) * t440;
t246 = t267 ^ 2;
t395 = pkin(3) * (t251 * t262 + t443) * t246;
t247 = t269 ^ 2;
t394 = pkin(3) * (t251 * t264 + t442) * t247;
t248 = t271 ^ 2;
t393 = pkin(3) * (t251 * t266 + t441) * t248;
t392 = t306 * t495;
t391 = t306 * t497;
t389 = t306 * t496;
t387 = g(3) * t273;
t384 = t191 * t484;
t383 = t192 * t483;
t382 = t193 * t482;
t381 = t194 * t484;
t380 = t195 * t483;
t379 = t196 * t482;
t378 = t73 * t446;
t376 = t75 * t445;
t374 = t77 * t444;
t373 = t258 * t487;
t372 = t258 * t486;
t371 = t258 * t485;
t368 = t155 * t515;
t367 = -t498 / 0.2e1;
t366 = -0.12e2 * pkin(3) * t514;
t365 = m(3) * t550;
t364 = pkin(2) * t418;
t361 = -0.4e1 * t405;
t360 = 0.4e1 * t405;
t359 = t498 * t526;
t358 = -0.2e1 * t307 + t532;
t355 = t535 * t418;
t352 = pkin(3) * t398;
t351 = t273 + t503;
t349 = t484 * t506;
t348 = t483 * t505;
t347 = t482 * t504;
t110 = -t351 * t512 + Icges(3,6);
t311 = -m(3) * rSges(3,1) * t351 + Icges(3,5);
t81 = t110 * t267 + t261 * t311;
t345 = t81 * t377;
t82 = t110 * t269 + t263 * t311;
t344 = t82 * t375;
t83 = t110 * t271 + t265 * t311;
t343 = t245 * t83 * t482;
t342 = t306 * t390;
t341 = t306 * t388;
t340 = rSges(3,1) * t368;
t339 = -0.4e1 * t404 + t358;
t338 = -0.2e1 * t169 - 0.4e1 * t307 + t289 - t305;
t337 = t190 / 0.2e1 + pkin(2) / 0.2e1;
t336 = t503 / 0.4e1 + t273 / 0.4e1;
t335 = rSges(3,1) * t261 + rSges(3,2) * t267;
t334 = rSges(3,1) * t263 + rSges(3,2) * t269;
t333 = rSges(3,1) * t265 + rSges(3,2) * t271;
t332 = t69 * t337;
t331 = t70 * t337;
t71 = (t105 * t435 + t108 * t434 - t177 * t278) * t128;
t330 = t71 * t337;
t329 = pkin(2) + t135;
t328 = pkin(2) + t136;
t327 = pkin(2) + t137;
t323 = (rSges(2,2) * m(2) - m(3) * t273) * t250;
t299 = 0.3e1 * qJ(3,1);
t241 = cos(t299);
t322 = pkin(3) * (-t271 * pkin(2) + t155 * t241);
t313 = (0.3e1 / 0.4e1 * t90 + 0.3e1 / 0.4e1 * t88 - 0.3e1 / 0.4e1 * t363) * t305 + t555 * t69;
t312 = (0.3e1 / 0.4e1 * t91 + 0.3e1 / 0.4e1 * t89 - 0.3e1 / 0.4e1 * t362) * t305 + t555 * t70;
t288 = 0.2e1 * t308;
t310 = Icges(1,3) + Icges(2,3) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + Icges(3,1) / 0.2e1 + t510 + (0.2e1 * t273 ^ 2 + t285 + t288 + t421) * t515 + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2 + t308) * m(2);
t304 = pkin(3) * t305;
t291 = -0.2e1 * pkin(7);
t283 = m(2) + m(3);
t239 = cos(t297);
t236 = cos(t294);
t232 = sin(t299);
t231 = sin(t298);
t229 = sin(t296);
t228 = sin(t295);
t226 = sin(t293);
t225 = sin(t292);
t222 = t296 - pkin(7);
t221 = t296 + pkin(7);
t219 = t293 - pkin(7);
t218 = t293 + pkin(7);
t217 = t291 + qJ(3,2);
t216 = t291 + qJ(3,3);
t214 = -pkin(7) + t297;
t212 = -pkin(7) + t294;
t210 = pkin(7) + t297;
t208 = pkin(7) + t294;
t207 = t290 + qJ(3,2);
t206 = t290 + qJ(3,3);
t185 = cos(t215);
t184 = cos(t214);
t183 = cos(t213);
t182 = cos(t212);
t181 = cos(t211);
t180 = cos(t210);
t179 = cos(t209);
t178 = cos(t208);
t168 = -Icges(3,4) + t407;
t162 = 0.2e1 * t215;
t161 = 0.2e1 * t213;
t160 = 0.2e1 * t211;
t159 = 0.2e1 * t209;
t145 = -0.2e1 * t396;
t142 = m(3) * t421 + Icges(3,3);
t134 = rSges(3,2) * t368;
t122 = (t302 / 0.2e1 - t301 / 0.2e1) * m(3) + t554;
t121 = g(1) * t196 - g(2) * t193;
t120 = g(1) * t195 - g(2) * t192;
t119 = g(1) * t194 - g(2) * t191;
t118 = g(1) * t193 + g(2) * t196;
t117 = g(1) * t192 + g(2) * t195;
t116 = g(1) * t191 + g(2) * t194;
t98 = t145 + t307 + 0.2e1 * t404 + t422;
t66 = t196 * t393 + (pkin(3) * t452 - t560 * t196) * t271 + t155 * t452;
t65 = -t193 * t393 + (pkin(3) * t447 + t560 * t193) * t271 + t155 * t447;
t64 = t195 * t394 + (pkin(3) * t454 - t561 * t195) * t269 + t155 * t454;
t63 = -t192 * t394 + (pkin(3) * t449 + t561 * t192) * t269 + t155 * t449;
t62 = t194 * t395 + (pkin(3) * t456 - t562 * t194) * t267 + t155 * t456;
t61 = -t191 * t395 + (pkin(3) * t451 + t562 * t191) * t267 + t155 * t451;
t58 = t242 * t508 + t137 * t365 + (-(-m(3) * t327 - t513) * t251 - t323) * t420 - t168 * t233 + t310;
t57 = t239 * t508 + t136 * t365 + (-(-m(3) * t328 - t513) * t251 - t323) * t420 - t168 * t230 + t310;
t56 = t236 * t508 + t135 * t365 + (-(-m(3) * t329 - t513) * t251 - t323) * t420 - t168 * t227 + t310;
t55 = t128 * t283 * t374 + t392 * t504;
t54 = t127 * t283 * t376 + t389 * t505;
t53 = t126 * t283 * t378 + t391 * t506;
t52 = (-t196 * t347 + t283 * t66) * t465;
t51 = (t193 * t347 + t283 * t65) * t465;
t50 = (-t195 * t348 + t283 * t64) * t466;
t49 = (t192 * t348 + t283 * t63) * t466;
t48 = (-t194 * t349 + t283 * t62) * t467;
t47 = (t191 * t349 + t283 * t61) * t467;
t46 = t142 * t392 + (-t177 * t83 + t374 * t504) * t128;
t45 = t142 * t389 + (-t176 * t82 + t376 * t505) * t127;
t44 = t142 * t391 + (-t175 * t81 + t378 * t506) * t126;
t31 = (t83 * t469 + (-t142 * t379 + t504 * t66) * t245) * t128;
t30 = (t82 * t471 + (-t142 * t380 + t505 * t64) * t244) * t127;
t29 = (t81 * t473 + (-t142 * t381 + t506 * t62) * t243) * t126;
t28 = (t83 * t475 + (t142 * t382 + t504 * t65) * t245) * t128;
t27 = (t82 * t477 + (t142 * t383 + t505 * t63) * t244) * t127;
t26 = (t81 * t479 + (t142 * t384 + t506 * t61) * t243) * t126;
t21 = t535 * t71 - t400;
t20 = t535 * t70 - t401;
t19 = t535 * t69 - t402;
t18 = (t21 - t400) * t71 * t128;
t17 = (t20 - t401) * t70 * t127;
t16 = (t19 - t402) * t69 * t126;
t15 = (pkin(3) * (-t155 - t201) * t522 * t245 + (-(t248 * t305 + t98) * t71 + (-t155 * t271 * t71 + t488 * t535) * t418) * t71 * t444) * t128;
t14 = (pkin(3) * (-t155 - t199) * t523 * t244 + (-(t247 * t305 + t98) * t70 + (-t155 * t269 * t70 + t490 * t535) * t418) * t70 * t445) * t127;
t13 = (pkin(3) * (-t155 - t197) * t524 * t243 + (-(t246 * t305 + t98) * t69 + (-t155 * t267 * t69 + t492 * t535) * t418) * t69 * t446) * t126;
t12 = ((-0.2e1 * t305 * t535 * t241 + (t166 * t535 + t346 - t481) * t418) * t43 + (-0.4e1 * (0.6e1 * t404 + t145 + t288 - t532 + t518) * t499 + t155 * t232 * t516 - t304 * t231 + 0.8e1 * (-pkin(2) * t169 + t281 * t428 - (0.3e1 / 0.4e1 * t305 + t284 + t422) * t190 + (t129 - t396) * t527) * t265) * t71) / (t468 + t430 / 0.2e1 - 0.2e1 * t397 - t308 + 0.2e1 * t322 + t358) * t71 * t507 + (t21 * t526 + t400 * t410 + (-0.2e1 * t481 + (-t241 + t271) * t355 + t531) * t71 + (-t305 * t231 + t232 * t412 + t233 * t520) * t43) / (0.4e1 * t322 + t338 - 0.4e1 * t397 + t430 + 0.2e1 * t468) * t43;
t11 = -m(2) * t118 - t283 * t15 + (-t137 * t12 - t333 * t522 - t118) * m(3);
t10 = -t83 * t18 - t15 * t504 - t142 * t12 + 0.2e1 * t71 ^ 2 * (t168 * t248 + (t122 * t265 + t134) * t271 + t265 * t340 + t426) - m(3) * (t118 * t137 + (-g(3) * t188 - t121 * t177) * t333);
t9 = -t58 * t18 - t83 * t12 - 0.4e1 * t43 * ((t71 * t265 * t508 + t43 * t511) * t271 + t488 * t509 + (t248 - 0.1e1 / 0.2e1) * t71 * t168 + ((t271 * t330 - t336 * t488) * rSges(3,2) + (t271 * t336 * t43 + t265 * t330) * rSges(3,1)) * m(3)) - m(1) * (g(3) * (-rSges(1,1) * t266 - rSges(1,2) * t272) + t121 * (rSges(1,1) * t272 - rSges(1,2) * t266)) - m(2) * (g(3) * (-rSges(2,1) * t177 - rSges(2,2) * t188 - t500) + t121 * (rSges(2,1) * t188 - rSges(2,2) * t177 + t202)) - m(3) * (-g(3) * t500 + t121 * t202 + (t121 * t327 + t387) * t188 + (-g(3) * t327 + t121 * t273) * t177);
t8 = (t403 * t489 + ((t312 + t556) * t174 + (t312 - t556) * t172) * t528 + ((-t341 / 0.3e1 + (-t458 + t540) * t127 + t493) * sin(t222) + (t341 / 0.3e1 + (-t458 - t540) * t127 + t493) * sin(t221)) * t416 + ((-t341 / 0.6e1 + (-t458 + t541) * t127 + t493) * sin(t214) + (t341 / 0.6e1 + (-t458 - t541) * t127 + t493) * sin(t210)) * t366 + (sin(t160) * t411 + t184 * t360) * (t341 / 0.2e1 + (-t458 - t539) * t127 + t493) + (sin(t162) * t411 + t180 * t361) * (-t341 / 0.2e1 + (-t458 + t539) * t127 + t493) + (t239 * t359 + t238 * t399 + t352 + (t269 * t367 + t427) * t549) * t42 + ((cos(t217) - cos(t207)) * t533 + (sin(t217) + sin(t207)) * t414 + (-t181 + t185) * t534 - t304 * t228 + t229 * t415 + t230 * t413) * t70) / (t203 * t239 + t431 / 0.2e1 + (cos(t162) / 0.2e1 + cos(t160) / 0.2e1 - 0.1e1) * t308 + t423 * t364 + ((t180 + t184) * t550 + (cos(t222) + cos(t221) - t181 - t185) * pkin(3)) * pkin(1) + t339) * t70 * t507 + (t20 * t526 + t401 * t410 + (t239 * t521 - t423 * t355 + t531) * t70 + (-t305 * t228 + t229 * t412 + t230 * t520) * t42) / (t239 * t519 + t431 + (t199 + t550) * t410 + (-t269 * pkin(2) + t155 * t238) * t549 + t338) * t42;
t7 = (t403 * t491 + ((t313 + t557) * t173 + (t313 - t557) * t171) * t528 + ((-t342 / 0.3e1 + (-t460 + t537) * t126 + t494) * sin(t219) + (t342 / 0.3e1 + (-t460 - t537) * t126 + t494) * sin(t218)) * t416 + ((-t342 / 0.6e1 + (-t460 + t538) * t126 + t494) * sin(t212) + (t342 / 0.6e1 + (-t460 - t538) * t126 + t494) * sin(t208)) * t366 + (sin(t159) * t411 + t182 * t360) * (t342 / 0.2e1 + (-t460 - t536) * t126 + t494) + (sin(t161) * t411 + t178 * t361) * (-t342 / 0.2e1 + (-t460 + t536) * t126 + t494) + (t236 * t359 + t235 * t399 + t352 + (t267 * t367 + t427) * t549) * t41 + ((cos(t216) - cos(t206)) * t533 + (sin(t216) + sin(t206)) * t414 + (-t179 + t183) * t534 - t304 * t225 + t226 * t415 + t227 * t413) * t69) / (t203 * t236 + t432 / 0.2e1 + (cos(t161) / 0.2e1 + cos(t159) / 0.2e1 - 0.1e1) * t308 + t424 * t364 + ((t178 + t182) * t550 + (cos(t219) + cos(t218) - t179 - t183) * pkin(3)) * pkin(1) + t339) * t69 * t507 + (t19 * t526 + t402 * t410 + (t236 * t521 - t424 * t355 + t531) * t69 + (-t305 * t225 + t226 * t412 + t227 * t520) * t41) / (t236 * t519 + t432 + (t197 + t550) * t410 + (-t267 * pkin(2) + t155 * t235) * t549 + t338) * t41;
t6 = -m(2) * t117 - t283 * t14 + (-t136 * t8 - t334 * t523 - t117) * m(3);
t5 = -m(2) * t116 - t283 * t13 + (-t135 * t7 - t335 * t524 - t116) * m(3);
t4 = -t82 * t17 - t14 * t505 - t142 * t8 + 0.2e1 * t70 ^ 2 * (t168 * t247 + (t122 * t263 + t134) * t269 + t263 * t340 + t426) - m(3) * (t117 * t136 + (-g(3) * t187 - t120 * t176) * t334);
t3 = -t81 * t16 - t13 * t506 - t142 * t7 + 0.2e1 * t69 ^ 2 * (t168 * t246 + (t122 * t261 + t134) * t267 + t261 * t340 + t426) - m(3) * (t116 * t135 + (-g(3) * t186 - t119 * t175) * t335);
t2 = -t57 * t17 - t82 * t8 - 0.4e1 * t42 * ((t42 * t511 + t489 * t508) * t269 + t490 * t509 + (t247 - 0.1e1 / 0.2e1) * t70 * t168 + ((t269 * t331 - t336 * t490) * rSges(3,2) + (t269 * t336 * t42 + t263 * t331) * rSges(3,1)) * m(3)) - m(1) * (g(3) * (-rSges(1,1) * t264 - rSges(1,2) * t270) + t120 * (rSges(1,1) * t270 - rSges(1,2) * t264)) - m(2) * (g(3) * (-rSges(2,1) * t176 - rSges(2,2) * t187 - t501) + t120 * (rSges(2,1) * t187 - rSges(2,2) * t176 + t200)) - m(3) * (-g(3) * t501 + t120 * t200 + (t120 * t328 + t387) * t187 + (-g(3) * t328 + t120 * t273) * t176);
t1 = -t56 * t16 - t81 * t7 - 0.4e1 * t41 * ((t41 * t511 + t491 * t508) * t267 + t492 * t509 + (t246 - 0.1e1 / 0.2e1) * t69 * t168 + ((t267 * t332 - t336 * t492) * rSges(3,2) + (t267 * t336 * t41 + t261 * t332) * rSges(3,1)) * m(3)) - m(1) * (g(3) * (-rSges(1,1) * t262 - rSges(1,2) * t268) + t119 * (rSges(1,1) * t268 - rSges(1,2) * t262)) - m(2) * (g(3) * (-rSges(2,1) * t175 - rSges(2,2) * t186 - t502) + t119 * (rSges(2,1) * t186 - rSges(2,2) * t175 + t198)) - m(3) * (-g(3) * t502 + t119 * t198 + (t119 * t329 + t387) * t186 + (-g(3) * t329 + t119 * t273) * t175);
t22 = [(-g(1) + t260) * m(4) + (t29 * t497 + t30 * t496 + t31 * t495) * t439 + (t9 * t469 + ((-t31 * t379 + t52 * t66) * t260 + (t31 * t382 + t52 * t65) * t259 + t52 * t371 + t66 * t11 - t10 * t379) * t245 + (t260 * t469 - t457 + t544) * (-t196 * t343 + t469 * t58) * t128) * t128 + (t2 * t471 + ((-t30 * t380 + t50 * t64) * t260 + (t30 * t383 + t50 * t63) * t259 + t50 * t372 + t64 * t6 - t4 * t380) * t244 + (t260 * t471 - t459 + t543) * (-t195 * t344 + t471 * t57) * t127) * t127 + (t1 * t473 + ((-t29 * t381 + t48 * t62) * t260 + (t29 * t384 + t48 * t61) * t259 + t48 * t373 + t62 * t5 - t3 * t381) * t243 + (t260 * t473 - t461 + t542) * (-t194 * t345 + t473 * t56) * t126) * t126; (-g(2) + t259) * m(4) + (t26 * t497 + t27 * t496 + t28 * t495) * t439 + (t9 * t475 + ((-t28 * t379 + t51 * t66) * t260 + (t28 * t382 + t51 * t65) * t259 + t51 * t371 + t65 * t11 + t10 * t382) * t245 + (t259 * t475 + t551) * (t193 * t343 + t475 * t58) * t128) * t128 + (t2 * t477 + ((-t27 * t380 + t49 * t64) * t260 + (t27 * t383 + t49 * t63) * t259 + t49 * t372 + t63 * t6 + t4 * t383) * t244 + (t259 * t477 + t552) * (t192 * t344 + t477 * t57) * t127) * t127 + (t1 * t479 + ((-t26 * t381 + t47 * t62) * t260 + (t26 * t384 + t47 * t61) * t259 + t47 * t373 + t61 * t5 + t3 * t384) * t243 + (t259 * t479 + t553) * (t191 * t345 + t479 * t56) * t126) * t126; (-g(3) + t258) * m(4) + (t10 * t495 + t3 * t497 + t4 * t496 + (t44 * t497 + t45 * t496 + t46 * t495) * t258) * t306 + (-t177 * t9 + (t544 + t551) * (-t128 * t177 * t58 + t392 * t83) + ((-t379 * t46 + t55 * t66) * t260 + (t382 * t46 + t55 * t65) * t259 + (t258 * t55 + t11) * t485) * t245) * t128 + (-t176 * t2 + (t543 + t552) * (-t127 * t176 * t57 + t389 * t82) + ((-t380 * t45 + t54 * t64) * t260 + (t383 * t45 + t54 * t63) * t259 + (t258 * t54 + t6) * t486) * t244) * t127 + (-t175 * t1 + (t542 + t553) * (-t126 * t175 * t56 + t391 * t81) + ((-t381 * t44 + t53 * t62) * t260 + (t384 * t44 + t53 * t61) * t259 + (t258 * t53 + t5) * t487) * t243) * t126;];
tauX  = t22;
