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
% Datum: 2020-08-07 11:13
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P4PRRRR8V2G1A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G1A0_invdyn_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR8V2G1A0_invdyn_para_pf_slag_vp1: xDP has to be [4x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [4 1]), ...
  'P4PRRRR8V2G1A0_invdyn_para_pf_slag_vp1: xDDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G1A0_invdyn_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR8V2G1A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G1A0_invdyn_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V2G1A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRRR8V2G1A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P4PRRRR8V2G1A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G1A0_invdyn_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G1A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:09:19
% EndTime: 2020-08-07 11:09:32
% DurationCPUTime: 12.68s
% Computational Cost: add. (58571->717), mult. (111835->1233), div. (5368->13), fcn. (112126->30), ass. (0->476)
t316 = cos(qJ(2,4));
t340 = pkin(7) + pkin(6);
t271 = t316 * t340;
t314 = sin(qJ(2,4));
t236 = pkin(2) * t314 - t271;
t306 = sin(pkin(4));
t315 = cos(qJ(3,4));
t308 = cos(pkin(4));
t313 = sin(qJ(3,4));
t435 = t308 * t313;
t450 = t306 * t314;
t299 = t315 ^ 2;
t513 = pkin(3) * t299;
t137 = 0.1e1 / ((pkin(3) * t435 + t236 * t306) * t315 + pkin(2) * t435 + t450 * t513);
t328 = cos(qJ(2,3));
t276 = t328 * t340;
t322 = sin(qJ(2,3));
t238 = pkin(2) * t322 - t276;
t327 = cos(qJ(3,3));
t321 = sin(qJ(3,3));
t433 = t308 * t321;
t446 = t306 * t322;
t301 = t327 ^ 2;
t512 = pkin(3) * t301;
t138 = 0.1e1 / ((pkin(3) * t433 + t238 * t306) * t327 + pkin(2) * t433 + t446 * t512);
t330 = cos(qJ(2,2));
t277 = t330 * t340;
t324 = sin(qJ(2,2));
t239 = pkin(2) * t324 - t277;
t329 = cos(qJ(3,2));
t323 = sin(qJ(3,2));
t431 = t308 * t323;
t444 = t306 * t324;
t302 = t329 ^ 2;
t511 = pkin(3) * t302;
t139 = 0.1e1 / ((pkin(3) * t431 + t239 * t306) * t329 + pkin(2) * t431 + t444 * t511);
t332 = cos(qJ(2,1));
t278 = t332 * t340;
t326 = sin(qJ(2,1));
t240 = pkin(2) * t326 - t278;
t331 = cos(qJ(3,1));
t325 = sin(qJ(3,1));
t429 = t308 * t325;
t442 = t306 * t326;
t303 = t331 ^ 2;
t510 = pkin(3) * t303;
t140 = 0.1e1 / ((pkin(3) * t429 + t240 * t306) * t331 + pkin(2) * t429 + t442 * t510);
t507 = g(3) * t308;
t378 = rSges(3,1) * t331 - rSges(3,2) * t325;
t549 = t378 * m(3);
t379 = rSges(3,1) * t329 - rSges(3,2) * t323;
t548 = t379 * m(3);
t380 = rSges(3,1) * t327 - rSges(3,2) * t321;
t547 = t380 * m(3);
t381 = rSges(3,1) * t315 - rSges(3,2) * t313;
t546 = t381 * m(3);
t419 = t326 * t340;
t243 = pkin(2) * t332 + t419;
t305 = sin(pkin(8));
t307 = cos(pkin(8));
t443 = t306 * t325;
t362 = pkin(3) * t443 - t240 * t308;
t541 = t243 * t307 + t362 * t305;
t421 = t324 * t340;
t242 = pkin(2) * t330 + t421;
t445 = t306 * t323;
t363 = pkin(3) * t445 - t239 * t308;
t540 = t242 * t307 + t363 * t305;
t423 = t322 * t340;
t241 = pkin(2) * t328 + t423;
t447 = t306 * t321;
t364 = pkin(3) * t447 - t238 * t308;
t539 = t241 * t307 + t364 * t305;
t425 = t314 * t340;
t237 = pkin(2) * t316 + t425;
t451 = t306 * t313;
t365 = pkin(3) * t451 - t236 * t308;
t538 = t237 * t307 + t365 * t305;
t312 = legFrame(1,3);
t286 = sin(t312);
t290 = cos(t312);
t226 = -g(1) * t286 + g(2) * t290;
t230 = g(1) * t290 + g(2) * t286;
t367 = t226 * t307 - t230 * t305;
t508 = g(3) * t306;
t537 = t367 * t308 + t508;
t311 = legFrame(2,3);
t285 = sin(t311);
t289 = cos(t311);
t225 = -g(1) * t285 + g(2) * t289;
t229 = g(1) * t289 + g(2) * t285;
t369 = t225 * t307 - t229 * t305;
t536 = t369 * t308 + t508;
t310 = legFrame(3,3);
t284 = sin(t310);
t288 = cos(t310);
t224 = -g(1) * t284 + g(2) * t288;
t228 = g(1) * t288 + g(2) * t284;
t371 = t224 * t307 - t228 * t305;
t535 = t371 * t308 + t508;
t309 = legFrame(4,3);
t283 = sin(t309);
t287 = cos(t309);
t223 = -g(1) * t283 + g(2) * t287;
t227 = g(1) * t287 + g(2) * t283;
t373 = t223 * t307 - t227 * t305;
t534 = t373 * t308 + t508;
t527 = m(3) * rSges(3,1);
t409 = rSges(3,2) * t527;
t279 = -Icges(3,4) + t409;
t341 = pkin(2) * m(3);
t407 = t341 / 0.2e1;
t387 = rSges(3,1) * t407;
t533 = t279 * t299 + t313 * t387;
t532 = t279 * t301 + t321 * t387;
t531 = t279 * t302 + t323 * t387;
t530 = t279 * t303 + t325 * t387;
t529 = 0.2e1 * pkin(2);
t528 = -0.2e1 * t279;
t185 = -t283 * t305 + t287 * t307;
t189 = t283 * t307 + t287 * t305;
t270 = pkin(3) * t315 + pkin(2);
t206 = t270 * t314 - t271;
t464 = (t270 * t316 + t425) * t308;
t125 = -t185 * t464 + t189 * t206;
t126 = -t185 * t206 - t189 * t464;
t342 = xP(4);
t297 = sin(t342);
t298 = cos(t342);
t345 = koppelP(4,2);
t349 = koppelP(4,1);
t219 = -t297 * t345 + t298 * t349;
t337 = xDP(4);
t338 = xDP(2);
t173 = t219 * t337 + t338;
t215 = t297 * t349 + t298 * t345;
t339 = xDP(1);
t177 = -t215 * t337 + t339;
t231 = t270 * t435;
t449 = t306 * t315;
t153 = 0.1e1 / (t206 * t449 + t231);
t356 = 0.1e1 / pkin(3);
t472 = t153 * t356;
t98 = (t125 * t177 + t126 * t173) * t472;
t526 = pkin(3) * t98;
t353 = rSges(3,2) ^ 2;
t354 = rSges(3,1) ^ 2;
t256 = (-t353 + t354) * m(3) + Icges(3,2) - Icges(3,1);
t525 = t256 / 0.2e1;
t333 = pkin(6) + rSges(3,3);
t517 = m(3) * t333;
t262 = rSges(3,2) * t517 - Icges(3,6);
t524 = -t262 / 0.4e1;
t263 = rSges(3,1) * t517 - Icges(3,5);
t523 = t263 / 0.4e1;
t522 = -t279 / 0.2e1;
t249 = rSges(3,1) * t313 + rSges(3,2) * t315;
t165 = -t249 * t450 + t381 * t308;
t521 = m(3) * t165;
t253 = rSges(3,1) * t321 + rSges(3,2) * t327;
t166 = -t253 * t446 + t380 * t308;
t520 = m(3) * t166;
t254 = rSges(3,1) * t323 + rSges(3,2) * t329;
t167 = -t254 * t444 + t379 * t308;
t519 = m(3) * t167;
t255 = rSges(3,1) * t325 + rSges(3,2) * t331;
t168 = -t255 * t442 + t378 * t308;
t518 = m(3) * t168;
t186 = -t284 * t305 + t288 * t307;
t190 = t284 * t307 + t288 * t305;
t273 = pkin(3) * t327 + pkin(2);
t210 = t273 * t322 - t276;
t463 = (t273 * t328 + t423) * t308;
t127 = -t186 * t463 + t190 * t210;
t130 = -t186 * t210 - t190 * t463;
t346 = koppelP(3,2);
t350 = koppelP(3,1);
t220 = -t297 * t346 + t298 * t350;
t174 = t220 * t337 + t338;
t216 = t297 * t350 + t298 * t346;
t178 = -t216 * t337 + t339;
t232 = t273 * t433;
t441 = t306 * t327;
t154 = 0.1e1 / (t210 * t441 + t232);
t471 = t154 * t356;
t102 = (t127 * t178 + t130 * t174) * t471;
t516 = pkin(3) * t102;
t187 = -t285 * t305 + t289 * t307;
t191 = t285 * t307 + t289 * t305;
t274 = pkin(3) * t329 + pkin(2);
t211 = t274 * t324 - t277;
t462 = (t274 * t330 + t421) * t308;
t128 = -t187 * t462 + t191 * t211;
t131 = -t187 * t211 - t191 * t462;
t347 = koppelP(2,2);
t351 = koppelP(2,1);
t221 = -t297 * t347 + t298 * t351;
t175 = t221 * t337 + t338;
t217 = t297 * t351 + t298 * t347;
t179 = -t217 * t337 + t339;
t233 = t274 * t431;
t439 = t306 * t329;
t155 = 0.1e1 / (t211 * t439 + t233);
t470 = t155 * t356;
t103 = (t128 * t179 + t131 * t175) * t470;
t515 = pkin(3) * t103;
t188 = -t286 * t305 + t290 * t307;
t192 = t286 * t307 + t290 * t305;
t275 = pkin(3) * t331 + pkin(2);
t212 = t275 * t326 - t278;
t461 = (t275 * t332 + t419) * t308;
t129 = -t188 * t461 + t192 * t212;
t132 = -t188 * t212 - t192 * t461;
t348 = koppelP(1,2);
t352 = koppelP(1,1);
t222 = -t297 * t348 + t298 * t352;
t176 = t222 * t337 + t338;
t218 = t297 * t352 + t298 * t348;
t180 = -t218 * t337 + t339;
t234 = t275 * t429;
t437 = t306 * t331;
t156 = 0.1e1 / (t212 * t437 + t234);
t469 = t156 * t356;
t104 = (t129 * t180 + t132 * t176) * t469;
t514 = pkin(3) * t104;
t300 = m(1) + m(2) + m(3);
t509 = g(3) * t300;
t357 = pkin(2) ^ 2;
t280 = t340 ^ 2 + t357;
t355 = pkin(3) ^ 2;
t408 = t313 * t526;
t417 = pkin(3) * t529;
t434 = t308 * t314;
t117 = -t185 * t449 - (t185 * t434 + t189 * t316) * t313;
t118 = -t189 * t449 - (-t185 * t316 + t189 * t434) * t313;
t86 = (t117 * t177 + t118 * t173) * t137;
t506 = (-t340 * t408 + (t299 * t355 + t315 * t417 + t280) * t86) * t86;
t402 = t321 * t516;
t432 = t308 * t322;
t119 = -t186 * t441 - (t186 * t432 + t190 * t328) * t321;
t122 = -t190 * t441 - (-t186 * t328 + t190 * t432) * t321;
t90 = (t119 * t178 + t122 * t174) * t138;
t505 = (-t340 * t402 + (t301 * t355 + t327 * t417 + t280) * t90) * t90;
t401 = t323 * t515;
t430 = t308 * t324;
t120 = -t187 * t439 - (t187 * t430 + t191 * t330) * t323;
t123 = -t191 * t439 - (-t187 * t330 + t191 * t430) * t323;
t91 = (t120 * t179 + t123 * t175) * t139;
t504 = (-t340 * t401 + (t302 * t355 + t329 * t417 + t280) * t91) * t91;
t400 = t325 * t514;
t428 = t308 * t326;
t121 = -t188 * t437 - (t188 * t428 + t192 * t332) * t325;
t124 = -t192 * t437 - (-t188 * t332 + t192 * t428) * t325;
t92 = (t121 * t180 + t124 * t176) * t140;
t503 = (-t340 * t400 + (t303 * t355 + t331 * t417 + t280) * t92) * t92;
t486 = t340 * t86;
t399 = t313 * t486;
t57 = t399 - t526;
t22 = (-t315 * t506 - (pkin(2) * t98 - t315 * t57) * t526) * t137;
t258 = m(2) * rSges(2,2) - t517;
t272 = m(2) * rSges(2,1) + t341;
t410 = 0.2e1 * m(3);
t490 = t316 * t86;
t85 = t86 ^ 2;
t94 = t98 ^ 2;
t502 = -t300 * t22 + ((-t272 * t85 - (t85 + t94) * t546) * t314 - (t249 * t98 * t410 + t258 * t86) * t490) * t306;
t485 = t340 * t90;
t398 = t321 * t485;
t59 = t398 - t516;
t26 = (-t327 * t505 - (pkin(2) * t102 - t327 * t59) * t516) * t138;
t489 = t328 * t90;
t87 = t90 ^ 2;
t99 = t102 ^ 2;
t501 = -t300 * t26 + ((-t272 * t87 - (t87 + t99) * t547) * t322 - (t253 * t102 * t410 + t258 * t90) * t489) * t306;
t100 = t103 ^ 2;
t484 = t340 * t91;
t397 = t323 * t484;
t60 = t397 - t515;
t27 = (-t329 * t504 - (pkin(2) * t103 - t329 * t60) * t515) * t139;
t488 = t330 * t91;
t88 = t91 ^ 2;
t499 = -t300 * t27 + ((-t272 * t88 - (t88 + t100) * t548) * t324 - (t254 * t103 * t410 + t258 * t91) * t488) * t306;
t101 = t104 ^ 2;
t483 = t340 * t92;
t396 = t325 * t483;
t61 = t396 - t514;
t28 = (-t331 * t503 - (pkin(2) * t104 - t331 * t61) * t514) * t140;
t487 = t332 * t92;
t89 = t92 ^ 2;
t498 = -t300 * t28 + ((-t272 * t89 - (t89 + t101) * t549) * t326 - (t255 * t104 * t410 + t258 * t92) * t487) * t306;
t497 = rSges(3,2) * t306;
t426 = t314 * t315;
t181 = pkin(3) * t426 + t236;
t427 = t308 * t356;
t17 = t137 * t427 * t506 + (-t308 * t399 + (-t181 * t451 + (pkin(2) * t315 + t513) * t308) * t98) / (t181 * t449 + t231) * t98;
t496 = t165 * t17;
t424 = t322 * t327;
t182 = pkin(3) * t424 + t238;
t18 = t138 * t427 * t505 + (-t308 * t398 + (-t182 * t447 + t308 * (pkin(2) * t327 + t512)) * t102) / (t182 * t441 + t232) * t102;
t495 = t166 * t18;
t422 = t324 * t329;
t183 = pkin(3) * t422 + t239;
t19 = t139 * t427 * t504 + (-t308 * t397 + (-t183 * t445 + t308 * (pkin(2) * t329 + t511)) * t103) / (t183 * t439 + t233) * t103;
t494 = t167 * t19;
t420 = t326 * t331;
t184 = pkin(3) * t420 + t240;
t20 = t140 * t427 * t503 + (-t308 * t396 + (-t184 * t443 + t308 * (pkin(2) * t331 + t510)) * t104) / (t184 * t437 + t234) * t104;
t493 = t168 * t20;
t492 = t249 * t94;
t491 = t253 * t99;
t482 = t100 * t254;
t481 = t101 * t255;
t304 = t337 ^ 2;
t317 = xDDP(4);
t320 = xDDP(1);
t149 = -t215 * t317 - t219 * t304 + t320;
t480 = t125 * t149;
t319 = xDDP(2);
t145 = -t215 * t304 + t219 * t317 + t319;
t479 = t126 * t145;
t150 = -t216 * t317 - t220 * t304 + t320;
t478 = t127 * t150;
t151 = -t217 * t317 - t221 * t304 + t320;
t477 = t128 * t151;
t152 = -t218 * t317 - t222 * t304 + t320;
t476 = t129 * t152;
t146 = -t216 * t304 + t220 * t317 + t319;
t475 = t130 * t146;
t147 = -t217 * t304 + t221 * t317 + t319;
t474 = t131 * t147;
t148 = -t218 * t304 + t222 * t317 + t319;
t473 = t132 * t148;
t193 = t272 + t546;
t161 = t193 * t316 - t258 * t314;
t468 = t161 * t306;
t196 = t272 + t547;
t162 = t196 * t328 - t258 * t322;
t467 = t162 * t306;
t197 = t272 + t548;
t163 = t197 * t330 - t258 * t324;
t466 = t163 * t306;
t198 = t272 + t549;
t164 = t198 * t332 - t258 * t326;
t465 = t164 * t306;
t257 = (t353 + t354) * m(3) + Icges(3,3);
t456 = t257 * t356;
t448 = t306 * t316;
t440 = t306 * t328;
t438 = t306 * t330;
t436 = t306 * t332;
t418 = rSges(3,2) * t507;
t412 = -t409 / 0.2e1 + Icges(3,4) / 0.2e1;
t411 = -0.2e1 * rSges(3,2) * pkin(2);
t406 = pkin(2) * t451;
t405 = pkin(2) * t447;
t404 = pkin(2) * t445;
t403 = pkin(2) * t443;
t169 = -t262 * t315 - t263 * t313;
t395 = t169 * t472;
t394 = t153 * t456;
t170 = -t262 * t327 - t263 * t321;
t393 = t170 * t471;
t392 = t154 * t456;
t171 = -t262 * t329 - t263 * t323;
t391 = t171 * t470;
t390 = t155 * t456;
t172 = -t262 * t331 - t263 * t325;
t389 = t172 * t469;
t388 = t156 * t456;
t282 = rSges(3,2) * t407;
t386 = t472 * t521;
t385 = t471 * t520;
t384 = t470 * t519;
t383 = t469 * t518;
t382 = (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,1) + Icges(2,3);
t372 = t223 * t305 + t227 * t307;
t370 = t224 * t305 + t228 * t307;
t368 = t225 * t305 + t229 * t307;
t366 = t226 * t305 + t230 * t307;
t361 = t534 * t314 + t372 * t316;
t360 = t535 * t322 + t370 * t328;
t359 = t536 * t324 + t368 * t330;
t358 = t537 * t326 + t366 * t332;
t344 = rSges(4,1);
t343 = rSges(4,2);
t318 = xDDP(3);
t293 = t527 * t529;
t259 = t333 ^ 2 + t353 + t357;
t235 = (t354 / 0.2e1 - t353 / 0.2e1) * m(3) - Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1;
t214 = -t343 * t297 + t344 * t298;
t213 = t344 * t297 + t343 * t298;
t204 = t305 * t332 + t307 * t428;
t203 = t305 * t330 + t307 * t430;
t202 = t305 * t328 + t307 * t432;
t201 = t305 * t428 - t307 * t332;
t200 = t305 * t430 - t307 * t330;
t199 = t305 * t432 - t307 * t328;
t195 = t305 * t316 + t307 * t434;
t194 = t305 * t434 - t307 * t316;
t144 = t243 * t305 - t362 * t307;
t143 = t242 * t305 - t363 * t307;
t142 = t241 * t305 - t364 * t307;
t141 = t237 * t305 - t365 * t307;
t136 = t256 * t303 + (t325 * t528 + t293) * t331 + (t325 * t411 + t259) * m(3) + t382;
t135 = t256 * t302 + (t323 * t528 + t293) * t329 + (t323 * t411 + t259) * m(3) + t382;
t134 = t256 * t301 + (t321 * t528 + t293) * t327 + (t321 * t411 + t259) * m(3) + t382;
t133 = t256 * t299 + (t313 * t528 + t293) * t315 + (t313 * t411 + t259) * m(3) + t382;
t116 = (-t201 * t286 + t204 * t290) * t510 + (t144 * t290 + t541 * t286) * t331 - t188 * t403;
t115 = (-t200 * t285 + t203 * t289) * t511 + (t143 * t289 + t540 * t285) * t329 - t187 * t404;
t114 = (-t199 * t284 + t202 * t288) * t512 + (t142 * t288 + t539 * t284) * t327 - t186 * t405;
t113 = -(t201 * t290 + t204 * t286) * t510 + (-t286 * t144 + t541 * t290) * t331 + t192 * t403;
t112 = -(t200 * t289 + t203 * t285) * t511 + (-t285 * t143 + t540 * t289) * t329 + t191 * t404;
t111 = -(t199 * t288 + t202 * t284) * t512 + (-t284 * t142 + t539 * t288) * t327 + t190 * t405;
t110 = (-t194 * t283 + t195 * t287) * t513 + (t141 * t287 + t538 * t283) * t315 - t185 * t406;
t109 = -(t194 * t287 + t195 * t283) * t513 + (-t283 * t141 + t538 * t287) * t315 + t189 * t406;
t108 = (-t129 * t218 + t132 * t222) * t469;
t107 = (-t128 * t217 + t131 * t221) * t470;
t106 = (-t127 * t216 + t130 * t220) * t471;
t105 = (-t125 * t215 + t126 * t219) * t472;
t97 = (-t121 * t218 + t124 * t222) * t140;
t96 = (-t120 * t217 + t123 * t221) * t139;
t95 = (-t119 * t216 + t122 * t220) * t138;
t93 = (-t117 * t215 + t118 * t219) * t137;
t84 = (-t113 * t218 + t116 * t222) * t140;
t83 = (-t112 * t217 + t115 * t221) * t139;
t82 = (-t111 * t216 + t114 * t220) * t138;
t81 = (-t109 * t215 + t110 * t219) * t137;
t80 = t132 * t388 + (t116 * t518 + t124 * t172) * t140;
t79 = t131 * t390 + (t115 * t519 + t123 * t171) * t139;
t78 = t130 * t392 + (t114 * t520 + t122 * t170) * t138;
t77 = t129 * t388 + (t113 * t518 + t121 * t172) * t140;
t76 = t128 * t390 + (t112 * t519 + t120 * t171) * t139;
t75 = t127 * t392 + (t111 * t520 + t119 * t170) * t138;
t74 = t132 * t383 + (t116 * t300 + t124 * t465) * t140;
t73 = t131 * t384 + (t115 * t300 + t123 * t466) * t139;
t72 = t130 * t385 + (t114 * t300 + t122 * t467) * t138;
t71 = t129 * t383 + (t113 * t300 + t121 * t465) * t140;
t70 = t128 * t384 + (t112 * t300 + t120 * t466) * t139;
t69 = t127 * t385 + (t111 * t300 + t119 * t467) * t138;
t68 = t126 * t394 + (t110 * t521 + t118 * t169) * t137;
t67 = t125 * t394 + (t109 * t521 + t117 * t169) * t137;
t66 = t126 * t386 + (t110 * t300 + t118 * t468) * t137;
t65 = t125 * t386 + (t109 * t300 + t117 * t468) * t137;
t56 = t132 * t389 + (t116 * t465 + t124 * t136) * t140;
t55 = t131 * t391 + (t115 * t466 + t123 * t135) * t139;
t54 = t130 * t393 + (t114 * t467 + t122 * t134) * t138;
t53 = t129 * t389 + (t113 * t465 + t121 * t136) * t140;
t52 = t128 * t391 + (t112 * t466 + t120 * t135) * t139;
t51 = t127 * t393 + (t111 * t467 + t119 * t134) * t138;
t50 = t126 * t395 + (t110 * t468 + t118 * t133) * t137;
t49 = t125 * t395 + (t109 * t468 + t117 * t133) * t137;
t45 = t108 * t518 + t300 * t84 + t97 * t465;
t44 = t107 * t519 + t300 * t83 + t96 * t466;
t43 = t106 * t520 + t300 * t82 + t95 * t467;
t41 = t105 * t521 + t300 * t81 + t93 * t468;
t40 = t108 * t172 + t136 * t97 + t84 * t465;
t39 = t107 * t171 + t135 * t96 + t83 * t466;
t38 = t106 * t170 + t134 * t95 + t82 * t467;
t37 = t105 * t169 + t133 * t93 + t81 * t468;
t16 = (((t104 * t308 + t92 * t436) * t510 + ((-t400 + t483) * t326 + pkin(2) * t487) * t437 + t61 * t308) * t92 + (t104 * t436 + (t303 * t308 - t420 * t443 - t308) * t92) * t514) * t140;
t15 = (((t103 * t308 + t91 * t438) * t511 + ((-t401 + t484) * t324 + pkin(2) * t488) * t439 + t60 * t308) * t91 + (t103 * t438 + (t302 * t308 - t422 * t445 - t308) * t91) * t515) * t139;
t14 = (((t102 * t308 + t90 * t440) * t512 + ((-t402 + t485) * t322 + pkin(2) * t489) * t441 + t59 * t308) * t90 + (t102 * t440 + (t301 * t308 - t424 * t447 - t308) * t90) * t516) * t138;
t13 = (((t308 * t98 + t86 * t448) * t513 + ((-t408 + t486) * t314 + pkin(2) * t490) * t449 + t57 * t308) * t86 + (t98 * t448 + (t299 * t308 - t426 * t451 - t308) * t86) * t526) * t137;
t12 = -t28 * t518 - t172 * t16 - t257 * t20 + 0.2e1 * ((t235 * t325 + t282) * t331 + t412 + t530) * t89 + (((t367 * t306 - t507) * rSges(3,1) + t358 * rSges(3,2)) * t331 + t325 * (t358 * rSges(3,1) - t367 * t497 + t418)) * m(3);
t11 = -t27 * t519 - t171 * t15 - t257 * t19 + 0.2e1 * ((t235 * t323 + t282) * t329 + t412 + t531) * t88 + (((t369 * t306 - t507) * rSges(3,1) + t359 * rSges(3,2)) * t329 + t323 * (t359 * rSges(3,1) - t369 * t497 + t418)) * m(3);
t10 = -t26 * t520 - t170 * t14 - t257 * t18 + 0.2e1 * ((t235 * t321 + t282) * t327 + t412 + t532) * t87 + (((t371 * t306 - t507) * rSges(3,1) + t360 * rSges(3,2)) * t327 + t321 * (t360 * rSges(3,1) - t371 * t497 + t418)) * m(3);
t9 = -t22 * t521 - t169 * t13 - t257 * t17 + 0.2e1 * ((t235 * t313 + t282) * t315 + t412 + t533) * t85 + (((t373 * t306 - t507) * rSges(3,1) + t361 * rSges(3,2)) * t315 + t313 * (t361 * rSges(3,1) - t373 * t497 + t418)) * m(3);
t8 = -t28 * t465 - t136 * t16 - t172 * t20 - 0.4e1 * t104 * ((t524 * t325 + t523 * t331) * t104 + ((t325 * t525 + t282) * t331 + t522 + t530) * t92) + (-t198 * t537 + t366 * t258) * t332 - (-t366 * t198 - t258 * t537) * t326;
t7 = -t27 * t466 - t135 * t15 - t171 * t19 - 0.4e1 * t103 * ((t524 * t323 + t523 * t329) * t103 + ((t323 * t525 + t282) * t329 + t522 + t531) * t91) + (-t197 * t536 + t368 * t258) * t330 - (-t368 * t197 - t258 * t536) * t324;
t6 = -t26 * t467 - t134 * t14 - t170 * t18 - 0.4e1 * t102 * ((t524 * t321 + t523 * t327) * t102 + ((t321 * t525 + t282) * t327 + t522 + t532) * t90) + (-t196 * t535 + t370 * t258) * t328 - (-t370 * t196 - t258 * t535) * t322;
t5 = -t22 * t468 - t133 * t13 - t169 * t17 - 0.4e1 * t98 * ((t524 * t313 + t523 * t315) * t98 + ((t313 * t525 + t282) * t315 + t522 + t533) * t86) + (-t193 * t534 + t372 * t258) * t316 - (-t372 * t193 - t258 * t534) * t314;
t4 = -t16 * t465 - t509 + (-t308 * t481 - t493) * m(3) + t498;
t3 = -t15 * t466 - t509 + (-t308 * t482 - t494) * m(3) + t499;
t2 = -t14 * t467 - t509 + (-t308 * t491 - t495) * m(3) + t501;
t1 = -t13 * t468 - t509 + (-t308 * t492 - t496) * m(3) + t502;
t21 = [(t69 + t70 + t71 + t65) * t318 + (t121 * t8 + t113 * t4 + (t116 * t71 + t124 * t53) * t148 + (t113 * t71 + t121 * t53) * t152) * t140 + (t120 * t7 + t112 * t3 + (t115 * t70 + t123 * t52) * t147 + (t112 * t70 + t120 * t52) * t151) * t139 + (t119 * t6 + t111 * t2 + (t111 * t69 + t119 * t51) * t150 + (t114 * t69 + t122 * t51) * t146) * t138 + ((t77 * t473 + (t152 * t77 + t12) * t129) * t156 + (t76 * t474 + (t151 * t76 + t11) * t128) * t155 + (t75 * t475 + (t150 * t75 + t10) * t127) * t154 + (t67 * t479 + (t149 * t67 + t9) * t125) * t153) * t356 + (t117 * t5 + t109 * t1 + (t109 * t65 + t117 * t49) * t149 + (t110 * t65 + t118 * t49) * t145) * t137 + (-t213 * t317 - t304 * t214 - g(1) + t320) * m(4); (t66 + t72 + t73 + t74) * t318 + (t124 * t8 + t116 * t4 + (t116 * t74 + t124 * t56) * t148 + (t113 * t74 + t121 * t56) * t152) * t140 + (t123 * t7 + t115 * t3 + (t112 * t73 + t120 * t55) * t151 + (t115 * t73 + t123 * t55) * t147) * t139 + (t122 * t6 + t114 * t2 + (t111 * t72 + t119 * t54) * t150 + (t114 * t72 + t122 * t54) * t146) * t138 + ((t80 * t476 + (t148 * t80 + t12) * t132) * t156 + (t79 * t477 + (t147 * t79 + t11) * t131) * t155 + (t78 * t478 + (t146 * t78 + t10) * t130) * t154 + (t68 * t480 + (t145 * t68 + t9) * t126) * t153) * t356 + (t118 * t5 + t110 * t1 + (t109 * t66 + t117 * t50) * t149 + (t110 * t66 + t118 * t50) * t145) * t137 + (-t304 * t213 + t214 * t317 - g(2) + t319) * m(4); t66 * t145 + t72 * t146 + t73 * t147 + t74 * t148 + t65 * t149 + t69 * t150 + t70 * t151 + t71 * t152 + (-t161 * t13 - t162 * t14 - t163 * t15 - t164 * t16) * t306 + (-t496 - t495 - t494 - t493 + (-t481 - t482 - t491 - t492) * t308) * m(3) + t498 + t499 + t501 + t502 + (-t318 + g(3)) * (-0.4e1 * t300 - m(4)); Icges(4,3) * t317 + t81 * t1 + t106 * t10 + t105 * t9 + t107 * t11 + t108 * t12 + t82 * t2 + t83 * t3 + t84 * t4 + t93 * t5 + t95 * t6 + t96 * t7 + t97 * t8 + ((t116 * t45 + t124 * t40) * t148 + (t113 * t45 + t121 * t40) * t152) * t140 + ((t112 * t44 + t120 * t39) * t151 + (t115 * t44 + t123 * t39) * t147) * t139 + ((t473 + t476) * (t108 * t257 + t172 * t97 + t84 * t518) * t156 + (t474 + t477) * (t107 * t257 + t171 * t96 + t83 * t519) * t155 + (t475 + t478) * (t106 * t257 + t170 * t95 + t82 * t520) * t154 + (t479 + t480) * (t105 * t257 + t169 * t93 + t81 * t521) * t153) * t356 + ((t111 * t43 + t119 * t38) * t150 + (t114 * t43 + t122 * t38) * t146) * t138 + ((t109 * t41 + t117 * t37) * t149 + (t110 * t41 + t118 * t37) * t145) * t137 + (t214 * t319 - t213 * t320 + (g(1) * t344 + g(2) * t343) * t297 + (g(1) * t343 - g(2) * t344) * t298 + (t343 ^ 2 + t344 ^ 2) * t317) * m(4) + (t45 + t41 + t43 + t44) * t318;];
tauX  = t21;
