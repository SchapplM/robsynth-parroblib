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
% Datum: 2020-09-20 23:03
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P4PRRRR8V1G3A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V1G3A0_invdyn_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR8V1G3A0_invdyn_para_pf_slag_vp1: xDP has to be [4x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [4 1]), ...
  'P4PRRRR8V1G3A0_invdyn_para_pf_slag_vp1: xDDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V1G3A0_invdyn_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR8V1G3A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P4PRRRR8V1G3A0_invdyn_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V1G3A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRRR8V1G3A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P4PRRRR8V1G3A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V1G3A0_invdyn_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V1G3A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-09-20 22:58:42
% EndTime: 2020-09-20 22:58:56
% DurationCPUTime: 14.09s
% Computational Cost: add. (40997->641), mult. (116475->1232), div. (13536->11), fcn. (119654->38), ass. (0->497)
t289 = cos(qJ(3,3));
t290 = cos(qJ(2,3));
t284 = sin(qJ(2,3));
t439 = t284 * t289;
t221 = pkin(2) * t439 - pkin(5) * t290;
t268 = sin(pkin(3));
t283 = sin(qJ(3,3));
t270 = cos(pkin(3));
t542 = pkin(2) * t270;
t178 = t221 * t268 + t283 * t542;
t562 = 0.1e1 / t178;
t488 = t562 / t289;
t273 = cos(qJ(3,4));
t274 = cos(qJ(2,4));
t272 = sin(qJ(2,4));
t441 = t272 * t273;
t211 = pkin(2) * t441 - pkin(5) * t274;
t271 = sin(qJ(3,4));
t174 = t211 * t268 + t271 * t542;
t563 = 0.1e1 / t174;
t489 = t563 / t273;
t278 = legFrame(1,2);
t250 = sin(t278);
t254 = cos(t278);
t216 = g(1) * t250 + g(2) * t254;
t220 = g(1) * t254 - g(2) * t250;
t269 = cos(pkin(6));
t246 = g(3) * t269;
t267 = sin(pkin(6));
t388 = -t220 * t267 - t246;
t581 = t216 * t268 + t388 * t270;
t277 = legFrame(2,2);
t249 = sin(t277);
t253 = cos(t277);
t215 = g(1) * t249 + g(2) * t253;
t219 = g(1) * t253 - g(2) * t249;
t390 = -t219 * t267 - t246;
t580 = t215 * t268 + t390 * t270;
t276 = legFrame(3,2);
t248 = sin(t276);
t252 = cos(t276);
t214 = g(1) * t248 + g(2) * t252;
t218 = g(1) * t252 - g(2) * t248;
t392 = -t218 * t267 - t246;
t579 = t214 * t268 + t392 * t270;
t275 = legFrame(4,2);
t247 = sin(t275);
t251 = cos(t275);
t213 = g(1) * t247 + g(2) * t251;
t217 = g(1) * t251 - g(2) * t247;
t394 = -t217 * t267 - t246;
t578 = t213 * t268 + t394 * t270;
t302 = xP(4);
t255 = sin(t302);
t256 = cos(t302);
t312 = koppelP(1,2);
t316 = koppelP(1,1);
t206 = t255 * t316 + t256 * t312;
t210 = -t255 * t312 + t256 * t316;
t298 = xDP(4);
t266 = t298 ^ 2;
t279 = xDDP(4);
t282 = xDDP(1);
t148 = -t206 * t279 - t210 * t266 + t282;
t494 = t148 * t254;
t281 = xDDP(2);
t144 = -t206 * t266 + t210 * t279 + t281;
t498 = t144 * t250;
t573 = -t494 + t498;
t311 = koppelP(2,2);
t315 = koppelP(2,1);
t205 = t255 * t315 + t256 * t311;
t209 = -t255 * t311 + t256 * t315;
t147 = -t205 * t279 - t209 * t266 + t282;
t495 = t147 * t253;
t143 = -t205 * t266 + t209 * t279 + t281;
t499 = t143 * t249;
t572 = -t495 + t499;
t310 = koppelP(3,2);
t314 = koppelP(3,1);
t204 = t255 * t314 + t256 * t310;
t208 = -t255 * t310 + t256 * t314;
t146 = -t204 * t279 - t208 * t266 + t282;
t496 = t146 * t252;
t142 = -t204 * t266 + t208 * t279 + t281;
t500 = t142 * t248;
t571 = -t496 + t500;
t309 = koppelP(4,2);
t313 = koppelP(4,1);
t203 = t255 * t313 + t256 * t309;
t207 = -t255 * t309 + t256 * t313;
t145 = -t203 * t279 - t207 * t266 + t282;
t497 = t145 * t251;
t141 = -t203 * t266 + t207 * t279 + t281;
t501 = t141 * t247;
t570 = -t497 + t501;
t287 = sin(qJ(3,1));
t293 = cos(qJ(3,1));
t569 = -rSges(3,1) * t293 + rSges(3,2) * t287;
t285 = sin(qJ(3,2));
t291 = cos(qJ(3,2));
t568 = -rSges(3,1) * t291 + rSges(3,2) * t285;
t567 = -rSges(3,1) * t289 + rSges(3,2) * t283;
t566 = -rSges(3,1) * t273 + rSges(3,2) * t271;
t453 = t270 * t272;
t342 = t268 * t273 + t271 * t453;
t442 = t271 * t274;
t149 = t342 * t267 - t269 * t442;
t150 = t267 * t442 + t342 * t269;
t299 = xDP(3);
t300 = xDP(2);
t301 = xDP(1);
t378 = (t207 * t298 + t300) * t247 - (-t203 * t298 + t301) * t251;
t95 = (t149 * t299 + t378 * t150) * t489;
t94 = t95 ^ 2;
t451 = t270 * t284;
t341 = t268 * t289 + t283 * t451;
t440 = t283 * t290;
t151 = t341 * t267 - t269 * t440;
t154 = t267 * t440 + t341 * t269;
t377 = (t208 * t298 + t300) * t248 - (-t204 * t298 + t301) * t252;
t102 = (t151 * t299 + t377 * t154) * t488;
t99 = t102 ^ 2;
t286 = sin(qJ(2,2));
t449 = t270 * t286;
t340 = t268 * t291 + t285 * t449;
t292 = cos(qJ(2,2));
t438 = t285 * t292;
t152 = t340 * t267 - t269 * t438;
t155 = t267 * t438 + t340 * t269;
t376 = (t209 * t298 + t300) * t249 - (-t205 * t298 + t301) * t253;
t263 = 0.1e1 / t291;
t437 = t286 * t291;
t400 = t268 * t437;
t450 = t270 * t285;
t455 = t268 * t292;
t565 = 0.1e1 / (-pkin(5) * t455 + (t400 + t450) * pkin(2));
t487 = t565 * t263;
t103 = (t152 * t299 + t376 * t155) * t487;
t100 = t103 ^ 2;
t288 = sin(qJ(2,1));
t447 = t270 * t288;
t339 = t268 * t293 + t287 * t447;
t294 = cos(qJ(2,1));
t436 = t287 * t294;
t153 = t339 * t267 - t269 * t436;
t156 = t267 * t436 + t339 * t269;
t375 = (t210 * t298 + t300) * t250 - (-t206 * t298 + t301) * t254;
t265 = 0.1e1 / t293;
t435 = t288 * t293;
t399 = t268 * t435;
t448 = t270 * t287;
t454 = t268 * t294;
t564 = 0.1e1 / (-pkin(5) * t454 + (t399 + t448) * pkin(2));
t486 = t564 * t265;
t104 = (t153 * t299 + t375 * t156) * t486;
t101 = t104 ^ 2;
t244 = m(3) * rSges(3,1) * rSges(3,2) - Icges(3,4);
t561 = 0.2e1 * t244;
t297 = m(2) * rSges(2,1);
t560 = m(3) * rSges(3,3);
t452 = t270 * t274;
t541 = pkin(2) * t273;
t133 = (-t267 * t272 + t269 * t452) * t541 + pkin(5) * (t267 * t274 + t269 * t453);
t134 = (t267 * t452 + t269 * t272) * t541 + (t267 * t453 - t269 * t274) * pkin(5);
t321 = 0.1e1 / pkin(2);
t79 = (t378 * t133 + t134 * t299) * t321 * t489;
t559 = pkin(2) * t79;
t446 = t270 * t290;
t540 = pkin(2) * t289;
t135 = (-t267 * t284 + t269 * t446) * t540 + pkin(5) * (t267 * t290 + t269 * t451);
t138 = (t267 * t446 + t269 * t284) * t540 + (t267 * t451 - t269 * t290) * pkin(5);
t87 = (t377 * t135 + t138 * t299) * t321 * t488;
t558 = pkin(2) * t87;
t445 = t270 * t292;
t539 = pkin(2) * t291;
t136 = (-t267 * t286 + t269 * t445) * t539 + pkin(5) * (t267 * t292 + t269 * t449);
t139 = (t267 * t445 + t269 * t286) * t539 + (t267 * t449 - t269 * t292) * pkin(5);
t88 = (t376 * t136 + t139 * t299) * t321 * t487;
t557 = pkin(2) * t88;
t444 = t270 * t294;
t538 = pkin(2) * t293;
t137 = (-t267 * t288 + t269 * t444) * t538 + pkin(5) * (t267 * t294 + t269 * t447);
t140 = (t267 * t444 + t269 * t288) * t538 + (t267 * t447 - t269 * t294) * pkin(5);
t89 = (t375 * t137 + t140 * t299) * t321 * t486;
t556 = pkin(2) * t89;
t555 = pkin(5) * t95;
t317 = rSges(3,2) ^ 2;
t318 = rSges(3,1) ^ 2;
t233 = (-t317 + t318) * m(3) - Icges(3,1) + Icges(3,2);
t554 = t233 / 0.2e1;
t242 = rSges(3,2) * t560 - Icges(3,6);
t553 = -t242 / 0.4e1;
t243 = rSges(3,1) * t560 - Icges(3,5);
t552 = t243 / 0.4e1;
t527 = rSges(3,2) * t273;
t227 = rSges(3,1) * t271 + t527;
t157 = -t268 * t272 * t227 - t270 * t566;
t551 = m(3) * t157;
t523 = rSges(3,2) * t289;
t228 = rSges(3,1) * t283 + t523;
t158 = -t268 * t284 * t228 - t270 * t567;
t550 = m(3) * t158;
t522 = rSges(3,2) * t291;
t229 = rSges(3,1) * t285 + t522;
t159 = -t268 * t286 * t229 - t270 * t568;
t549 = m(3) * t159;
t521 = rSges(3,2) * t293;
t230 = rSges(3,1) * t287 + t521;
t160 = -t268 * t288 * t230 - t270 * t569;
t548 = m(3) * t160;
t547 = m(3) * t321;
t257 = t273 ^ 2;
t546 = pkin(2) * t257;
t260 = t289 ^ 2;
t545 = pkin(2) * t260;
t262 = t291 ^ 2;
t544 = pkin(2) * t262;
t264 = t293 ^ 2;
t543 = pkin(2) * t264;
t537 = pkin(5) * t102;
t536 = pkin(5) * t103;
t535 = pkin(5) * t104;
t534 = g(3) * t267;
t319 = pkin(5) ^ 2;
t320 = pkin(2) ^ 2;
t519 = t271 * t79;
t433 = pkin(2) * t519;
t533 = (-pkin(5) * t433 + (t257 * t320 + t319) * t95) * t95;
t518 = t283 * t87;
t432 = pkin(2) * t518;
t520 = t102 * (-pkin(5) * t432 + (t260 * t320 + t319) * t102);
t517 = t285 * t88;
t516 = t287 * t89;
t515 = t103 * t565;
t514 = t104 * t564;
t212 = pkin(5) * t272 + t274 * t541;
t461 = t268 * t271;
t370 = pkin(2) * t461 - t211 * t270;
t126 = -t212 * t267 + t370 * t269;
t513 = t126 * t563;
t224 = pkin(5) * t284 + t290 * t540;
t459 = t268 * t283;
t369 = pkin(2) * t459 - t221 * t270;
t130 = -t224 * t267 + t369 * t269;
t512 = t130 * t562;
t225 = pkin(5) * t286 + t292 * t539;
t222 = pkin(2) * t437 - pkin(5) * t292;
t458 = t268 * t285;
t368 = pkin(2) * t458 - t222 * t270;
t131 = -t225 * t267 + t368 * t269;
t179 = pkin(2) * t450 + t222 * t268;
t176 = 0.1e1 / t179;
t511 = t131 * t176;
t226 = pkin(5) * t288 + t294 * t538;
t223 = pkin(2) * t435 - pkin(5) * t294;
t457 = t268 * t287;
t367 = pkin(2) * t457 - t223 * t270;
t132 = -t226 * t267 + t367 * t269;
t180 = pkin(2) * t448 + t223 * t268;
t177 = 0.1e1 / t180;
t510 = t132 * t177;
t509 = t133 * t321;
t508 = t134 * t321;
t507 = t135 * t321;
t506 = t136 * t321;
t505 = t137 * t321;
t504 = t138 * t321;
t503 = t139 * t321;
t502 = t140 * t321;
t199 = -m(3) * t566 + t297;
t241 = m(2) * rSges(2,2) - t560;
t161 = t199 * t274 - t241 * t272;
t493 = t161 * t268;
t200 = -m(3) * t567 + t297;
t162 = t200 * t290 - t241 * t284;
t492 = t162 * t268;
t201 = -m(3) * t568 + t297;
t163 = t201 * t292 - t241 * t286;
t491 = t163 * t268;
t202 = -m(3) * t569 + t297;
t164 = t202 * t294 - t241 * t288;
t490 = t164 * t268;
t259 = m(1) + m(2) + m(3);
t485 = t563 * t259;
t484 = t562 * t259;
t483 = t176 * t259;
t482 = t177 * t259;
t193 = -t242 * t273 - t243 * t271;
t481 = t193 * t321;
t194 = -t242 * t289 - t243 * t283;
t480 = t194 * t321;
t195 = -t242 * t291 - t243 * t285;
t479 = t195 * t321;
t196 = -t242 * t293 - t243 * t287;
t478 = t196 * t321;
t476 = t213 * t270;
t474 = t214 * t270;
t472 = t215 * t270;
t470 = t216 * t270;
t465 = t233 * t271;
t434 = t317 + t318;
t238 = t434 * m(3) + Icges(3,3);
t464 = t238 * t321;
t463 = rSges(3,2) * t268 * t246;
t462 = t267 * t268;
t460 = t268 * t274;
t456 = t268 * t290;
t443 = t270 * t321;
t431 = pkin(2) * t517;
t430 = pkin(2) * t516;
t429 = t271 * t555;
t428 = t563 * t551;
t427 = t157 * t547;
t426 = t562 * t550;
t425 = t158 * t547;
t424 = t176 * t549;
t423 = t159 * t547;
t422 = t177 * t548;
t421 = t160 * t547;
t420 = t283 * t537;
t419 = t285 * t536;
t418 = t287 * t535;
t417 = t563 * t493;
t416 = t562 * t492;
t415 = t176 * t491;
t414 = t177 * t490;
t413 = t247 * t489;
t412 = t251 * t489;
t411 = t248 * t488;
t410 = t252 * t488;
t409 = t249 * t487;
t408 = t253 * t487;
t407 = t250 * t486;
t406 = t254 * t486;
t405 = t233 * t283 * t289;
t404 = t233 * t285 * t291;
t403 = t233 * t287 * t293;
t402 = t268 * t441;
t401 = t268 * t439;
t393 = -t217 * t269 + t534;
t391 = -t218 * t269 + t534;
t389 = -t219 * t269 + t534;
t387 = -t220 * t269 + t534;
t125 = t212 * t269 + t370 * t267;
t113 = t125 * t251 + t174 * t247;
t114 = -t125 * t247 + t174 * t251;
t386 = t113 * t145 + t114 * t141;
t127 = t224 * t269 + t369 * t267;
t115 = t127 * t252 + t178 * t248;
t116 = -t127 * t248 + t178 * t252;
t385 = t115 * t146 + t116 * t142;
t128 = t225 * t269 + t368 * t267;
t117 = t128 * t253 + t179 * t249;
t118 = -t128 * t249 + t179 * t253;
t384 = t117 * t147 + t118 * t143;
t129 = t226 * t269 + t367 * t267;
t119 = t129 * t254 + t180 * t250;
t120 = -t129 * t250 + t180 * t254;
t383 = t119 * t148 + t120 * t144;
t305 = 0.2e1 * qJ(3,4);
t326 = Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1 + (0.2e1 * rSges(3,3) ^ 2 + t434) * m(3) / 0.2e1;
t121 = cos(t305) * t554 - t244 * sin(t305) + t326;
t41 = t429 - t559;
t13 = (((t270 * t79 + t95 * t460) * t546 - (t433 - t555) * t402 + t270 * t41) * t95 - (-t79 * t460 + (-t257 * t270 + t271 * t402 + t270) * t95) * t559) * t489;
t17 = (t443 * t533 + (-t79 * t211 * t461 + t270 * (t79 * t546 - t429)) * t79) * t489;
t21 = (t41 * t559 - t533) * t563;
t5 = -t21 * t493 - t121 * t13 - t193 * t17 - 0.4e1 * ((t95 * t465 / 0.2e1 + t79 * t552) * t273 + t519 * t553 + (t257 - 0.1e1 / 0.2e1) * t95 * t244) * t79 + (-t199 * t578 - t393 * t241) * t274 - (t393 * t199 - t241 * t578) * t272;
t325 = t578 * t272 - t393 * t274;
t9 = -t21 * t551 - t193 * t13 - t238 * t17 + t94 * (t257 * t561 + t273 * t465 - t244) + (((t394 * t268 - t476) * rSges(3,1) + t325 * rSges(3,2)) * t273 + (t463 + (t217 * t462 + t476) * rSges(3,2) + t325 * rSges(3,1)) * t271) * m(3);
t374 = -t150 * t5 - t9 * t509;
t43 = t420 - t558;
t14 = (((t102 * t456 + t270 * t87) * t545 - (t432 - t537) * t401 + t270 * t43) * t102 - (-t87 * t456 + (-t260 * t270 + t283 * t401 + t270) * t102) * t558) * t488;
t18 = (t443 * t520 + (-t87 * t221 * t459 + t270 * (t87 * t545 - t420)) * t87) * t488;
t22 = (t43 * t558 - t520) * t562;
t324 = t579 * t284 - t391 * t290;
t10 = -t22 * t550 - t194 * t14 - t238 * t18 + t99 * (t260 * t561 - t244 + t405) + (((t392 * t268 - t474) * rSges(3,1) + t324 * rSges(3,2)) * t289 + (t463 + (t218 * t462 + t474) * rSges(3,2) + t324 * rSges(3,1)) * t283) * m(3);
t306 = 0.2e1 * qJ(3,3);
t122 = cos(t306) * t554 - t244 * sin(t306) + t326;
t6 = -t22 * t492 - t122 * t14 - t194 * t18 - 0.4e1 * ((t283 * t553 + t289 * t552) * t87 + (t405 / 0.2e1 + (t260 - 0.1e1 / 0.2e1) * t244) * t102) * t87 + (-t200 * t579 - t391 * t241) * t290 - (t391 * t200 - t241 * t579) * t284;
t373 = -t10 * t507 - t154 * t6;
t44 = t419 - t557;
t15 = (((t103 * t455 + t270 * t88) * t544 - (t431 - t536) * t400 + t270 * t44) * t515 + (t88 * t455 + (t262 * t270 - t285 * t400 - t270) * t103) * t565 * t557) * t263;
t27 = -pkin(5) * t431 + (t262 * t320 + t319) * t103;
t19 = (t27 * t443 * t515 + (-t88 * t222 * t458 + t270 * (t88 * t544 - t419)) * t176 * t88) * t263;
t23 = (t103 * t27 - t44 * t557) * t565;
t323 = t580 * t286 - t389 * t292;
t11 = t23 * t549 - t195 * t15 - t238 * t19 + t100 * (t262 * t561 - t244 + t404) + (((t390 * t268 - t472) * rSges(3,1) + t323 * rSges(3,2)) * t291 + (t463 + (t219 * t462 + t472) * rSges(3,2) + t323 * rSges(3,1)) * t285) * m(3);
t307 = 0.2e1 * qJ(3,2);
t123 = cos(t307) * t554 - t244 * sin(t307) + t326;
t7 = t23 * t491 - t123 * t15 - t195 * t19 - 0.4e1 * ((t285 * t553 + t291 * t552) * t88 + (t404 / 0.2e1 + (t262 - 0.1e1 / 0.2e1) * t244) * t103) * t88 + (-t201 * t580 - t389 * t241) * t292 - (t389 * t201 - t241 * t580) * t286;
t372 = -t11 * t506 - t155 * t7;
t45 = t418 - t556;
t16 = (((t104 * t454 + t270 * t89) * t543 - (t430 - t535) * t399 + t270 * t45) * t514 + (t89 * t454 + (t264 * t270 - t287 * t399 - t270) * t104) * t564 * t556) * t265;
t28 = -pkin(5) * t430 + (t264 * t320 + t319) * t104;
t20 = (t28 * t443 * t514 + (-t89 * t223 * t457 + t270 * (t89 * t543 - t418)) * t177 * t89) * t265;
t24 = (t104 * t28 - t45 * t556) * t564;
t322 = t581 * t288 - t387 * t294;
t12 = t24 * t548 - t196 * t16 - t238 * t20 + t101 * (t264 * t561 - t244 + t403) + (((t388 * t268 - t470) * rSges(3,1) + t322 * rSges(3,2)) * t293 + (t463 + (t220 * t462 + t470) * rSges(3,2) + t322 * rSges(3,1)) * t287) * m(3);
t308 = 0.2e1 * qJ(3,1);
t124 = cos(t308) * t554 - t244 * sin(t308) + t326;
t8 = t24 * t490 - t124 * t16 - t196 * t20 - 0.4e1 * ((t287 * t553 + t293 * t552) * t89 + (t403 / 0.2e1 + (t264 - 0.1e1 / 0.2e1) * t244) * t104) * t89 + (-t202 * t581 - t387 * t241) * t294 - (t387 * t202 - t241 * t581) * t288;
t371 = -t12 * t505 - t156 * t8;
t350 = t121 * t150 + t133 * t481;
t49 = t113 * t417 - t350 * t412;
t346 = t133 * t464 + t150 * t193;
t59 = t113 * t428 - t346 * t412;
t365 = t150 * t49 + t59 * t509;
t50 = t114 * t417 + t350 * t413;
t60 = t114 * t428 + t346 * t413;
t364 = t150 * t50 + t60 * t509;
t349 = t122 * t154 + t135 * t480;
t51 = t115 * t416 - t349 * t410;
t345 = t135 * t464 + t154 * t194;
t68 = t115 * t426 - t345 * t410;
t361 = t154 * t51 + t68 * t507;
t52 = t116 * t416 + t349 * t411;
t69 = t116 * t426 + t345 * t411;
t360 = t154 * t52 + t69 * t507;
t348 = t123 * t155 + t136 * t479;
t53 = t117 * t415 - t348 * t408;
t344 = t136 * t464 + t155 * t195;
t70 = t117 * t424 - t344 * t408;
t357 = t155 * t53 + t70 * t506;
t54 = t118 * t415 + t348 * t409;
t71 = t118 * t424 + t344 * t409;
t356 = t155 * t54 + t71 * t506;
t347 = t124 * t156 + t137 * t478;
t55 = t119 * t414 - t347 * t406;
t343 = t137 * t464 + t156 * t196;
t72 = t119 * t422 - t343 * t406;
t353 = t156 * t55 + t72 * t505;
t56 = t120 * t414 + t347 * t407;
t73 = t120 * t422 + t343 * t407;
t352 = t156 * t56 + t73 * t505;
t338 = (t203 * t251 + t207 * t247) * t489;
t337 = (t204 * t252 + t208 * t248) * t488;
t336 = (t205 * t253 + t209 * t249) * t487;
t335 = (t206 * t254 + t210 * t250) * t486;
t280 = xDDP(3);
t334 = t126 * t280 + t386;
t333 = t130 * t280 + t385;
t332 = t131 * t280 + t384;
t331 = t132 * t280 + t383;
t330 = t133 * t427 + t150 * t493;
t329 = t135 * t425 + t154 * t492;
t328 = t136 * t423 + t155 * t491;
t327 = t137 * t421 + t156 * t490;
t304 = rSges(4,1);
t303 = rSges(4,2);
t198 = -t255 * t303 + t256 * t304;
t197 = t255 * t304 + t256 * t303;
t112 = t156 * t335;
t111 = t155 * t336;
t110 = t154 * t337;
t109 = t150 * t338;
t108 = t335 * t505;
t107 = t336 * t506;
t106 = t337 * t507;
t105 = t338 * t509;
t98 = (-t119 * t206 + t120 * t210) * t177;
t97 = (-t117 * t205 + t118 * t209) * t176;
t96 = (-t115 * t204 + t116 * t208) * t562;
t93 = (-t113 * t203 + t114 * t207) * t563;
t92 = t132 * t422 + (t140 * t464 + t153 * t196) * t486;
t91 = t131 * t424 + (t139 * t464 + t152 * t195) * t487;
t90 = t130 * t426 + (t138 * t464 + t151 * t194) * t488;
t86 = t89 ^ 2;
t85 = t88 ^ 2;
t84 = t87 ^ 2;
t83 = t126 * t428 + (t134 * t464 + t149 * t193) * t489;
t82 = t132 * t482 + (t140 * t421 + t153 * t490) * t486;
t81 = t131 * t483 + (t139 * t423 + t152 * t491) * t487;
t80 = t130 * t484 + (t138 * t425 + t151 * t492) * t488;
t78 = t79 ^ 2;
t77 = t126 * t485 + (t134 * t427 + t149 * t493) * t489;
t76 = t132 * t414 + (t124 * t153 + t140 * t478) * t486;
t75 = t131 * t415 + (t123 * t152 + t139 * t479) * t487;
t74 = t130 * t416 + (t122 * t151 + t138 * t480) * t488;
t67 = t126 * t417 + (t121 * t149 + t134 * t481) * t489;
t40 = t108 * t238 + t112 * t196 + t98 * t548;
t39 = t107 * t238 + t111 * t195 + t97 * t549;
t38 = t106 * t238 + t110 * t194 + t96 * t550;
t37 = t108 * t548 + t112 * t490 + t259 * t98;
t36 = t107 * t549 + t111 * t491 + t259 * t97;
t35 = t106 * t550 + t110 * t492 + t259 * t96;
t34 = t105 * t238 + t109 * t193 + t93 * t551;
t33 = t105 * t551 + t109 * t493 + t259 * t93;
t32 = t108 * t196 + t112 * t124 + t98 * t490;
t31 = t107 * t195 + t111 * t123 + t97 * t491;
t30 = t106 * t194 + t110 * t122 + t96 * t492;
t29 = t105 * t193 + t109 * t121 + t93 * t493;
t4 = (-t16 * t164 + (-t241 * t294 - t288 * t297) * t101) * t268 + (t24 - t216) * t259 + (-t160 * t20 + (-0.2e1 * t294 * t104 * (rSges(3,1) * t516 + t89 * t521) + t569 * t288 * (t101 + t86)) * t268 - t86 * t270 * t230) * m(3);
t3 = (-t15 * t163 + (-t241 * t292 - t286 * t297) * t100) * t268 + (t23 - t215) * t259 + (-t159 * t19 + (-0.2e1 * t292 * t103 * (rSges(3,1) * t517 + t88 * t522) + t568 * t286 * (t100 + t85)) * t268 - t85 * t270 * t229) * m(3);
t2 = (-t14 * t162 + (-t241 * t290 - t284 * t297) * t99) * t268 + (-t22 - t214) * t259 + (-t158 * t18 + (-0.2e1 * t290 * t102 * (rSges(3,1) * t518 + t87 * t523) + t567 * t284 * (t99 + t84)) * t268 - t84 * t270 * t228) * m(3);
t1 = (-t13 * t161 + (-t241 * t274 - t272 * t297) * t94) * t268 + (-t21 - t213) * t259 + (-t157 * t17 + (-0.2e1 * t274 * t95 * (rSges(3,1) * t519 + t79 * t527) + t566 * t272 * (t94 + t78)) * t268 - t78 * t270 * t227) * m(3);
t25 = [(t119 * t4 + t331 * (t119 * t482 - t327 * t406)) * t177 + (t117 * t3 + t332 * (t117 * t483 - t328 * t408)) * t176 + (t115 * t2 + t333 * (t115 * t484 - t329 * t410)) * t562 + (t113 * t1 + t334 * (t113 * t485 - t330 * t412)) * t563 + (-t197 * t279 - t266 * t198 - g(1) + t282) * m(4) + ((t153 * t55 + t72 * t502) * t280 + t353 * t498 + (-t353 * t148 + t371) * t254) * t486 + ((t152 * t53 + t70 * t503) * t280 + t357 * t499 + (-t357 * t147 + t372) * t253) * t487 + ((t151 * t51 + t68 * t504) * t280 + t361 * t500 + (-t361 * t146 + t373) * t252) * t488 + ((t149 * t49 + t59 * t508) * t280 + t365 * t501 + (-t365 * t145 + t374) * t251) * t489; (t120 * t4 + t331 * (t120 * t482 + t327 * t407)) * t177 + (t118 * t3 + t332 * (t118 * t483 + t328 * t409)) * t176 + (t116 * t2 + t333 * (t116 * t484 + t329 * t411)) * t562 + (t114 * t1 + t334 * (t114 * t485 + t330 * t413)) * t563 + (-t266 * t197 + t198 * t279 - g(2) + t281) * m(4) + ((t153 * t56 + t73 * t502) * t280 - t352 * t494 + (t352 * t144 - t371) * t250) * t486 + ((t152 * t54 + t71 * t503) * t280 - t356 * t495 + (t356 * t143 - t372) * t249) * t487 + ((t151 * t52 + t69 * t504) * t280 - t360 * t496 + (t360 * t142 - t373) * t248) * t488 + ((t149 * t50 + t60 * t508) * t280 - t364 * t497 + (t364 * t141 - t374) * t247) * t489; -m(4) * g(3) + (t132 * t4 + t383 * t82) * t177 + (t131 * t3 + t384 * t81) * t176 + (t130 * t2 + t385 * t80) * t562 + (t126 * t1 + t386 * t77) * t563 + (t82 * t510 + t81 * t511 + t80 * t512 + t77 * t513 + m(4)) * t280 + (t12 * t502 + (t153 * t76 + t92 * t502) * t280 + t153 * t8 + t573 * (t156 * t76 + t92 * t505)) * t486 + (t11 * t503 + (t152 * t75 + t91 * t503) * t280 + t152 * t7 + t572 * (t155 * t75 + t91 * t506)) * t487 + (t10 * t504 + (t151 * t74 + t90 * t504) * t280 + t151 * t6 + t571 * (t154 * t74 + t90 * t507)) * t488 + (t9 * t508 + (t149 * t67 + t83 * t508) * t280 + t149 * t5 + t570 * (t150 * t67 + t83 * t509)) * t489; Icges(4,3) * t279 + t93 * t1 + t106 * t10 + t105 * t9 + t107 * t11 + t108 * t12 + t109 * t5 + t110 * t6 + t111 * t7 + t112 * t8 + t96 * t2 + t97 * t3 + t98 * t4 + t383 * t37 * t177 + t384 * t36 * t176 + t385 * t35 * t562 + t386 * t33 * t563 + ((t303 ^ 2 + t304 ^ 2) * t279 + (g(1) * t304 + g(2) * t303) * t255 + (g(1) * t303 - g(2) * t304) * t256 + t198 * t281 - t197 * t282) * m(4) + (t29 * t149 * t489 + t30 * t151 * t488 + t31 * t152 * t487 + t32 * t153 * t486 + t33 * t513 + t35 * t512 + t36 * t511 + t37 * t510 + (t134 * t34 * t489 + t138 * t38 * t488 + t139 * t39 * t487 + t140 * t40 * t486) * t321) * t280 + t570 * t489 * (t150 * t29 + t34 * t509) + t571 * t488 * (t154 * t30 + t38 * t507) + t572 * t487 * (t155 * t31 + t39 * t506) + t573 * t486 * (t156 * t32 + t40 * t505);];
tauX  = t25;
