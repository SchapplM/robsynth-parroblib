% Analytische Jacobi-Matrix für parallelen Roboter
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorpose und aktiven Gelenkkoordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Input:
% xP [6x1]
%   Generalized platform coordinates
% qJ [3x6]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [6x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,alpha3,alpha4,d2,d3,d4,theta1]';
% koppelP [6x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% Jinv [6x6]
%   Analytische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-18 16:42
% Revision: 36f6366a01c4a552c0708fcd8ed3e0fb9da693e2 (2019-05-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P6PRRRRR6V2A2_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,6),zeros(10,1),zeros(6,3),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [3 6]), ...
  'P6PRRRRR6V2A2_Jinv: qJ has to be [3x6] (double)');
assert(isreal(xP) && all(size(xP) == [6 1]), ...
  'P6PRRRRR6V2A2_Jinv: xP has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'P6PRRRRR6V2A2_Jinv: pkin has to be [10x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [6 3]), ...
  'P6PRRRRR6V2A2_Jinv: legFrame has to be [6x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6PRRRRR6V2A2_Jinv: Koppelpunkt has to be [6x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-18 15:50:44
% EndTime: 2019-05-18 15:52:13
% DurationCPUTime: 225.70s
% Computational Cost: add. (5664->718), mult. (14400->1461), div. (36->12), fcn. (14388->74), ass. (0->556)
t338 = cos(pkin(5));
t325 = t338 ^ 2;
t574 = t325 - 0.1e1;
t339 = cos(pkin(4));
t434 = t574 * t339;
t334 = sin(pkin(5));
t573 = pkin(2) * t334;
t333 = sin(pkin(6));
t572 = pkin(3) * t333;
t571 = pkin(9) * t333;
t359 = cos(qJ(2,6));
t474 = t339 * t359;
t335 = sin(pkin(4));
t485 = t334 * t335;
t224 = t338 * t474 - t485;
t332 = sin(pkin(10));
t336 = cos(pkin(10));
t353 = sin(qJ(2,6));
t483 = t338 * t353;
t181 = t224 * t332 + t336 * t483;
t364 = legFrame(6,2);
t307 = sin(t364);
t570 = t181 * t307;
t361 = cos(qJ(2,5));
t473 = t339 * t361;
t225 = t338 * t473 - t485;
t355 = sin(qJ(2,5));
t482 = t338 * t355;
t182 = t225 * t332 + t336 * t482;
t365 = legFrame(5,2);
t308 = sin(t365);
t569 = t182 * t308;
t363 = cos(qJ(2,4));
t472 = t339 * t363;
t226 = t338 * t472 - t485;
t357 = sin(qJ(2,4));
t481 = t338 * t357;
t183 = t226 * t332 + t336 * t481;
t366 = legFrame(4,2);
t309 = sin(t366);
t568 = t183 * t309;
t184 = t224 * t336 - t332 * t483;
t567 = t184 * t307;
t185 = t225 * t336 - t332 * t482;
t566 = t185 * t308;
t186 = t226 * t336 - t332 * t481;
t565 = t186 * t309;
t377 = cos(qJ(2,3));
t468 = t339 * t377;
t239 = t338 * t468 - t485;
t371 = sin(qJ(2,3));
t480 = t338 * t371;
t187 = t239 * t332 + t336 * t480;
t367 = legFrame(3,2);
t310 = sin(t367);
t564 = t187 * t310;
t379 = cos(qJ(2,2));
t467 = t339 * t379;
t240 = t338 * t467 - t485;
t373 = sin(qJ(2,2));
t479 = t338 * t373;
t188 = t240 * t332 + t336 * t479;
t368 = legFrame(2,2);
t311 = sin(t368);
t563 = t188 * t311;
t381 = cos(qJ(2,1));
t466 = t339 * t381;
t241 = t338 * t466 - t485;
t375 = sin(qJ(2,1));
t478 = t338 * t375;
t189 = t241 * t332 + t336 * t478;
t369 = legFrame(1,2);
t312 = sin(t369);
t562 = t189 * t312;
t190 = t239 * t336 - t332 * t480;
t561 = t190 * t310;
t191 = t240 * t336 - t332 * t479;
t560 = t191 * t311;
t192 = t241 * t336 - t332 * t478;
t559 = t192 * t312;
t383 = xP(5);
t320 = sin(t383);
t323 = cos(t383);
t385 = koppelP(6,3);
t382 = xP(6);
t319 = sin(t382);
t322 = cos(t382);
t391 = koppelP(6,2);
t397 = koppelP(6,1);
t416 = -t319 * t391 + t322 * t397;
t199 = t320 * t385 + t323 * t416;
t558 = t199 * t353;
t557 = t199 * t359;
t386 = koppelP(5,3);
t392 = koppelP(5,2);
t398 = koppelP(5,1);
t415 = -t319 * t392 + t322 * t398;
t200 = t320 * t386 + t323 * t415;
t556 = t200 * t355;
t555 = t200 * t361;
t387 = koppelP(4,3);
t393 = koppelP(4,2);
t399 = koppelP(4,1);
t414 = -t319 * t393 + t322 * t399;
t201 = t320 * t387 + t323 * t414;
t554 = t201 * t357;
t553 = t201 * t363;
t388 = koppelP(3,3);
t394 = koppelP(3,2);
t400 = koppelP(3,1);
t413 = -t319 * t394 + t322 * t400;
t202 = t320 * t388 + t323 * t413;
t552 = t202 * t371;
t551 = t202 * t377;
t389 = koppelP(2,3);
t395 = koppelP(2,2);
t401 = koppelP(2,1);
t412 = -t319 * t395 + t322 * t401;
t203 = t320 * t389 + t323 * t412;
t550 = t203 * t373;
t549 = t203 * t379;
t390 = koppelP(1,3);
t396 = koppelP(1,2);
t402 = koppelP(1,1);
t411 = -t319 * t396 + t322 * t402;
t204 = t320 * t390 + t323 * t411;
t548 = t204 * t375;
t547 = t204 * t381;
t477 = t339 * t353;
t242 = t332 * t359 + t336 * t477;
t546 = t242 * t307;
t476 = t339 * t355;
t243 = t332 * t361 + t336 * t476;
t545 = t243 * t308;
t475 = t339 * t357;
t244 = t332 * t363 + t336 * t475;
t544 = t244 * t309;
t245 = -t332 * t477 + t336 * t359;
t543 = t245 * t307;
t246 = -t332 * t476 + t336 * t361;
t542 = t246 * t308;
t247 = -t332 * t475 + t336 * t363;
t541 = t247 * t309;
t471 = t339 * t371;
t249 = t332 * t377 + t336 * t471;
t540 = t249 * t310;
t470 = t339 * t373;
t250 = t332 * t379 + t336 * t470;
t539 = t250 * t311;
t469 = t339 * t375;
t251 = t332 * t381 + t336 * t469;
t538 = t251 * t312;
t252 = -t332 * t471 + t336 * t377;
t537 = t252 * t310;
t253 = -t332 * t470 + t336 * t379;
t536 = t253 * t311;
t254 = -t332 * t469 + t336 * t381;
t535 = t254 * t312;
t337 = cos(pkin(6));
t279 = pkin(9) * t337 + pkin(8);
t534 = t279 * t359;
t533 = t279 * t361;
t532 = t279 * t363;
t531 = t279 * t377;
t530 = t279 * t379;
t529 = t279 * t381;
t346 = legFrame(6,1);
t287 = sin(t346);
t528 = t287 * t307;
t347 = legFrame(5,1);
t288 = sin(t347);
t527 = t288 * t308;
t348 = legFrame(4,1);
t289 = sin(t348);
t526 = t289 * t309;
t349 = legFrame(3,1);
t290 = sin(t349);
t525 = t290 * t310;
t350 = legFrame(2,1);
t291 = sin(t350);
t524 = t291 * t311;
t351 = legFrame(1,1);
t292 = sin(t351);
t523 = t292 * t312;
t299 = cos(t346);
t522 = t299 * t307;
t300 = cos(t347);
t521 = t300 * t308;
t301 = cos(t348);
t520 = t301 * t309;
t302 = cos(t349);
t519 = t302 * t310;
t303 = cos(t350);
t518 = t303 * t311;
t304 = cos(t351);
t517 = t304 * t312;
t384 = xP(4);
t321 = sin(t384);
t516 = t321 * t391;
t515 = t321 * t392;
t514 = t321 * t393;
t513 = t321 * t394;
t512 = t321 * t395;
t511 = t321 * t396;
t510 = t321 * t397;
t509 = t321 * t398;
t508 = t321 * t399;
t507 = t321 * t400;
t506 = t321 * t401;
t505 = t321 * t402;
t504 = t323 * t385;
t503 = t323 * t386;
t502 = t323 * t387;
t501 = t323 * t388;
t500 = t323 * t389;
t499 = t323 * t390;
t324 = cos(t384);
t498 = t324 * t391;
t497 = t324 * t392;
t496 = t324 * t393;
t495 = t324 * t394;
t494 = t324 * t395;
t493 = t324 * t396;
t492 = t324 * t397;
t491 = t324 * t398;
t490 = t324 * t399;
t489 = t324 * t400;
t488 = t324 * t401;
t487 = t324 * t402;
t486 = t333 * t339;
t484 = t334 * t338;
t255 = t339 * t279 * t484;
t280 = t339 * pkin(2) * pkin(3);
t465 = t255 * t571 + t280;
t272 = pkin(3) * t486;
t432 = t325 * t272;
t463 = pkin(3) * t571;
t437 = -0.2e1 * t463;
t464 = 0.2e1 * pkin(9) * t432 + t339 * t437;
t462 = t199 * t485;
t461 = t199 * t477;
t460 = t199 * t474;
t459 = t200 * t485;
t458 = t200 * t476;
t457 = t200 * t473;
t456 = t201 * t485;
t455 = t201 * t475;
t454 = t201 * t472;
t453 = t202 * t485;
t452 = t202 * t471;
t451 = t202 * t468;
t450 = t203 * t485;
t449 = t203 * t470;
t448 = t203 * t467;
t447 = t204 * t485;
t446 = t204 * t469;
t445 = t204 * t466;
t403 = pkin(9) ^ 2;
t271 = t337 ^ 2 * t403 + pkin(3) ^ 2 - t403;
t444 = t271 * t485;
t443 = t359 * t484;
t442 = t361 * t484;
t441 = t363 * t484;
t440 = t377 * t484;
t439 = t379 * t484;
t438 = t381 * t484;
t436 = t279 * t574;
t435 = pkin(2) * pkin(9) * t486;
t433 = pkin(9) * (t337 + 0.1e1) * (t337 - 0.1e1) * t334;
t273 = t353 * t573;
t431 = t325 * t534 + t273 - t534;
t274 = t355 * t573;
t430 = t325 * t533 + t274 - t533;
t275 = t357 * t573;
t429 = t325 * t532 + t275 - t532;
t276 = t371 * t573;
t428 = t325 * t531 + t276 - t531;
t277 = t373 * t573;
t427 = t325 * t530 + t277 - t530;
t278 = t375 * t573;
t426 = t325 * t529 + t278 - t529;
t425 = t271 * t434;
t424 = t353 * t433;
t423 = t355 * t433;
t422 = t357 * t433;
t421 = t371 * t433;
t420 = t373 * t433;
t419 = t375 * t433;
t418 = 0.2e1 * t463 * t485;
t417 = pkin(3) * t255 - t435;
t340 = legFrame(6,3);
t281 = sin(t340);
t293 = cos(t340);
t227 = t281 * t336 + t293 * t332;
t233 = -t281 * t332 + t293 * t336;
t410 = (-t227 * t353 + t233 * t474) * t338 - t233 * t485;
t341 = legFrame(5,3);
t282 = sin(t341);
t294 = cos(t341);
t228 = t282 * t336 + t294 * t332;
t234 = -t282 * t332 + t294 * t336;
t409 = (-t228 * t355 + t234 * t473) * t338 - t234 * t485;
t342 = legFrame(4,3);
t283 = sin(t342);
t295 = cos(t342);
t229 = t283 * t336 + t295 * t332;
t235 = -t283 * t332 + t295 * t336;
t408 = (-t229 * t357 + t235 * t472) * t338 - t235 * t485;
t343 = legFrame(3,3);
t284 = sin(t343);
t296 = cos(t343);
t230 = t284 * t336 + t296 * t332;
t236 = -t284 * t332 + t296 * t336;
t407 = (-t230 * t371 + t236 * t468) * t338 - t236 * t485;
t344 = legFrame(2,3);
t285 = sin(t344);
t297 = cos(t344);
t231 = t285 * t336 + t297 * t332;
t237 = -t285 * t332 + t297 * t336;
t406 = (-t231 * t373 + t237 * t467) * t338 - t237 * t485;
t345 = legFrame(1,3);
t286 = sin(t345);
t298 = cos(t345);
t232 = t286 * t336 + t298 * t332;
t238 = -t286 * t332 + t298 * t336;
t405 = (-t232 * t375 + t238 * t466) * t338 - t238 * t485;
t404 = (t272 - t432) * pkin(9);
t380 = cos(qJ(3,1));
t378 = cos(qJ(3,2));
t376 = cos(qJ(3,3));
t374 = sin(qJ(3,1));
t372 = sin(qJ(3,2));
t370 = sin(qJ(3,3));
t362 = cos(qJ(3,4));
t360 = cos(qJ(3,5));
t358 = cos(qJ(3,6));
t356 = sin(qJ(3,4));
t354 = sin(qJ(3,5));
t352 = sin(qJ(3,6));
t331 = t380 ^ 2;
t330 = t378 ^ 2;
t329 = t376 ^ 2;
t328 = t362 ^ 2;
t327 = t360 ^ 2;
t326 = t358 ^ 2;
t318 = cos(t369);
t317 = cos(t368);
t316 = cos(t367);
t315 = cos(t366);
t314 = cos(t365);
t313 = cos(t364);
t267 = t319 * t402 + t322 * t396;
t266 = t319 * t401 + t322 * t395;
t265 = t319 * t400 + t322 * t394;
t264 = t319 * t399 + t322 * t393;
t263 = t319 * t398 + t322 * t392;
t262 = t319 * t397 + t322 * t391;
t261 = t375 * t418;
t260 = t373 * t418;
t259 = t371 * t418;
t258 = t357 * t418;
t257 = t355 * t418;
t256 = t353 * t418;
t220 = t375 * t444;
t219 = t373 * t444;
t218 = t371 * t444;
t214 = t357 * t444;
t213 = t355 * t444;
t212 = t353 * t444;
t210 = t335 * t438 - t434;
t209 = t335 * t439 - t434;
t208 = t335 * t440 - t434;
t207 = t335 * t441 - t434;
t206 = t335 * t442 - t434;
t205 = t335 * t443 - t434;
t198 = t320 * t411 - t499;
t197 = t320 * t412 - t500;
t196 = t320 * t413 - t501;
t195 = t320 * t414 - t502;
t194 = t320 * t415 - t503;
t193 = t320 * t416 - t504;
t177 = t232 * t381 + t238 * t469;
t176 = t231 * t379 + t237 * t470;
t175 = t230 * t377 + t236 * t471;
t171 = t229 * t363 + t235 * t475;
t170 = t228 * t361 + t234 * t476;
t169 = t227 * t359 + t233 * t477;
t168 = t232 * t523 - t238 * t304;
t167 = t231 * t524 - t237 * t303;
t166 = t230 * t525 - t236 * t302;
t165 = t229 * t526 - t235 * t301;
t164 = t228 * t527 - t234 * t300;
t163 = t227 * t528 - t233 * t299;
t162 = t232 * t304 + t238 * t523;
t161 = t231 * t303 + t237 * t524;
t160 = t230 * t302 + t236 * t525;
t159 = t229 * t301 + t235 * t526;
t158 = t228 * t300 + t234 * t527;
t157 = t227 * t299 + t233 * t528;
t156 = -t232 * t292 + t238 * t517;
t155 = t232 * t517 + t238 * t292;
t154 = -t231 * t291 + t237 * t518;
t153 = t231 * t518 + t237 * t291;
t152 = -t230 * t290 + t236 * t519;
t151 = t230 * t519 + t236 * t290;
t150 = -t229 * t289 + t235 * t520;
t149 = t229 * t520 + t235 * t289;
t148 = -t228 * t288 + t234 * t521;
t147 = t228 * t521 + t234 * t288;
t146 = -t227 * t287 + t233 * t522;
t145 = t227 * t522 + t233 * t287;
t144 = t198 * t324 - t267 * t321;
t143 = t197 * t324 - t266 * t321;
t142 = t196 * t324 - t265 * t321;
t141 = t195 * t324 - t264 * t321;
t140 = t194 * t324 - t263 * t321;
t139 = t193 * t324 - t262 * t321;
t138 = t198 * t321 + t267 * t324;
t137 = t197 * t321 + t266 * t324;
t136 = t196 * t321 + t265 * t324;
t135 = t195 * t321 + t264 * t324;
t134 = t194 * t321 + t263 * t324;
t133 = t193 * t321 + t262 * t324;
t132 = t255 + (t381 * t436 + t278) * t335;
t131 = t255 + (t379 * t436 + t277) * t335;
t130 = t255 + (t377 * t436 + t276) * t335;
t129 = t255 + (t363 * t436 + t275) * t335;
t128 = t255 + (t361 * t436 + t274) * t335;
t127 = t255 + (t359 * t436 + t273) * t335;
t126 = (-t324 * t499 + (t320 * t487 - t511) * t322 - t319 * (t320 * t493 + t505)) * t318 - t204 * t517;
t125 = (-t321 * t499 + (t320 * t505 + t493) * t322 + t319 * (-t320 * t511 + t487)) * t318 - t204 * t523;
t124 = (-t324 * t500 + (t320 * t488 - t512) * t322 - t319 * (t320 * t494 + t506)) * t317 - t203 * t518;
t123 = (-t321 * t500 + (t320 * t506 + t494) * t322 + t319 * (-t320 * t512 + t488)) * t317 - t203 * t524;
t122 = (-t324 * t501 + (t320 * t489 - t513) * t322 - t319 * (t320 * t495 + t507)) * t316 - t202 * t519;
t121 = (-t321 * t501 + (t320 * t507 + t495) * t322 + t319 * (-t320 * t513 + t489)) * t316 - t202 * t525;
t120 = (-t324 * t502 + (t320 * t490 - t514) * t322 - t319 * (t320 * t496 + t508)) * t315 - t201 * t520;
t119 = (-t321 * t502 + (t320 * t508 + t496) * t322 + t319 * (-t320 * t514 + t490)) * t315 - t201 * t526;
t118 = (-t324 * t503 + (t320 * t491 - t515) * t322 - t319 * (t320 * t497 + t509)) * t314 - t200 * t521;
t117 = (-t321 * t503 + (t320 * t509 + t497) * t322 + t319 * (-t320 * t515 + t491)) * t314 - t200 * t527;
t116 = (-t324 * t504 + (t320 * t492 - t516) * t322 - t319 * (t320 * t498 + t510)) * t313 - t199 * t522;
t115 = (-t321 * t504 + (t320 * t510 + t498) * t322 + t319 * (-t320 * t516 + t492)) * t313 - t199 * t528;
t114 = t162 * t469 + t168 * t381;
t113 = t161 * t470 + t167 * t379;
t112 = t160 * t471 + t166 * t377;
t111 = t155 * t381 + t156 * t469;
t110 = t153 * t379 + t154 * t470;
t109 = t151 * t377 + t152 * t471;
t108 = t159 * t475 + t165 * t363;
t107 = t158 * t476 + t164 * t361;
t106 = t157 * t477 + t163 * t359;
t105 = t149 * t363 + t150 * t475;
t104 = t147 * t361 + t148 * t476;
t103 = t145 * t359 + t146 * t477;
t102 = t126 * t469 - t292 * t547;
t101 = t125 * t469 + t304 * t547;
t100 = t124 * t470 - t291 * t549;
t99 = t123 * t470 + t303 * t549;
t98 = t122 * t471 - t290 * t551;
t97 = t121 * t471 + t302 * t551;
t96 = t125 * t381 - t304 * t446;
t95 = t126 * t381 + t292 * t446;
t94 = t123 * t379 - t303 * t449;
t93 = t124 * t379 + t291 * t449;
t92 = t121 * t377 - t302 * t452;
t91 = t122 * t377 + t290 * t452;
t90 = t120 * t475 - t289 * t553;
t89 = t119 * t475 + t301 * t553;
t88 = t118 * t476 - t288 * t555;
t87 = t117 * t476 + t300 * t555;
t86 = t116 * t477 - t287 * t557;
t85 = t115 * t477 + t299 * t557;
t84 = t119 * t363 - t301 * t455;
t83 = t117 * t361 - t300 * t458;
t82 = t115 * t359 - t299 * t461;
t81 = t120 * t363 + t289 * t455;
t80 = t118 * t361 + t288 * t458;
t79 = t116 * t359 + t287 * t461;
t78 = (t162 * t466 - t168 * t375) * t338 - t162 * t485;
t77 = (t161 * t467 - t167 * t373) * t338 - t161 * t485;
t76 = (t160 * t468 - t166 * t371) * t338 - t160 * t485;
t75 = (-t155 * t375 + t156 * t466) * t338 - t156 * t485;
t74 = (-t153 * t373 + t154 * t467) * t338 - t154 * t485;
t73 = (-t151 * t371 + t152 * t468) * t338 - t152 * t485;
t72 = (t159 * t472 - t165 * t357) * t338 - t159 * t485;
t71 = (t158 * t473 - t164 * t355) * t338 - t158 * t485;
t70 = (t157 * t474 - t163 * t353) * t338 - t157 * t485;
t69 = (-t149 * t357 + t150 * t472) * t338 - t150 * t485;
t68 = (-t147 * t355 + t148 * t473) * t338 - t148 * t485;
t67 = (-t145 * t353 + t146 * t474) * t338 - t146 * t485;
t66 = (t125 * t375 + t304 * t445) * t338 - t304 * t447;
t65 = (-t126 * t375 + t292 * t445) * t338 - t292 * t447;
t64 = (t123 * t373 + t303 * t448) * t338 - t303 * t450;
t63 = (-t124 * t373 + t291 * t448) * t338 - t291 * t450;
t62 = (t121 * t371 + t302 * t451) * t338 - t302 * t453;
t61 = (-t122 * t371 + t290 * t451) * t338 - t290 * t453;
t60 = (t119 * t357 + t301 * t454) * t338 - t301 * t456;
t59 = (t117 * t355 + t300 * t457) * t338 - t300 * t459;
t58 = (t115 * t353 + t299 * t460) * t338 - t299 * t462;
t57 = (-t120 * t357 + t289 * t454) * t338 - t289 * t456;
t56 = (-t118 * t355 + t288 * t457) * t338 - t288 * t459;
t55 = (-t116 * t353 + t287 * t460) * t338 - t287 * t462;
t54 = (t126 * t466 + t292 * t548) * t338 - t126 * t485;
t53 = (t125 * t466 - t304 * t548) * t338 - t125 * t485;
t52 = (t124 * t467 + t291 * t550) * t338 - t124 * t485;
t51 = (t123 * t467 - t303 * t550) * t338 - t123 * t485;
t50 = (t122 * t468 + t290 * t552) * t338 - t122 * t485;
t49 = (t121 * t468 - t302 * t552) * t338 - t121 * t485;
t48 = (t120 * t472 + t289 * t554) * t338 - t120 * t485;
t47 = (t119 * t472 - t301 * t554) * t338 - t119 * t485;
t46 = (t118 * t473 + t288 * t556) * t338 - t118 * t485;
t45 = (t117 * t473 - t300 * t556) * t338 - t117 * t485;
t44 = (t116 * t474 + t287 * t558) * t338 - t116 * t485;
t43 = (t115 * t474 - t299 * t558) * t338 - t115 * t485;
t42 = 0.1e1 / ((t210 * t437 + t220) * t331 + ((t210 * t271 + t261) * t374 - t435 + t132 * pkin(3)) * t380 + (t132 * t571 + t280) * t374 + (t210 * t572 - t335 * t419) * pkin(9));
t41 = 0.1e1 / ((t209 * t437 + t219) * t330 + ((t209 * t271 + t260) * t372 - t435 + t131 * pkin(3)) * t378 + (t131 * t571 + t280) * t372 + (t209 * t572 - t335 * t420) * pkin(9));
t40 = 0.1e1 / ((t208 * t437 + t218) * t329 + ((t208 * t271 + t259) * t370 - t435 + t130 * pkin(3)) * t376 + (t130 * t571 + t280) * t370 + (t208 * t572 - t335 * t421) * pkin(9));
t39 = 0.1e1 / ((t207 * t437 + t214) * t328 + ((t207 * t271 + t258) * t356 - t435 + t129 * pkin(3)) * t362 + (t129 * t571 + t280) * t356 + (t207 * t572 - t335 * t422) * pkin(9));
t38 = 0.1e1 / ((t206 * t437 + t213) * t327 + ((t206 * t271 + t257) * t354 - t435 + t128 * pkin(3)) * t360 + (t128 * t571 + t280) * t354 + (t206 * t572 - t335 * t423) * pkin(9));
t37 = 0.1e1 / ((t205 * t437 + t212) * t326 + ((t205 * t271 + t256) * t352 - t435 + t127 * pkin(3)) * t358 + (t127 * t571 + t280) * t352 + (t205 * t572 - t335 * t424) * pkin(9));
t36 = (t102 * t336 + t332 * t95) * t298 + (-t102 * t332 + t336 * t95) * t286;
t35 = (t101 * t336 + t332 * t96) * t298 + (-t101 * t332 + t336 * t96) * t286;
t34 = (t100 * t336 + t332 * t93) * t297 + (-t100 * t332 + t336 * t93) * t285;
t33 = (t332 * t94 + t336 * t99) * t297 + (-t332 * t99 + t336 * t94) * t285;
t32 = (t332 * t91 + t336 * t98) * t296 + (-t332 * t98 + t336 * t91) * t284;
t31 = (t332 * t92 + t336 * t97) * t296 + (-t332 * t97 + t336 * t92) * t284;
t30 = (t332 * t81 + t336 * t90) * t295 + (-t332 * t90 + t336 * t81) * t283;
t29 = (t332 * t84 + t336 * t89) * t295 + t283 * (-t332 * t89 + t336 * t84);
t28 = (t332 * t80 + t336 * t88) * t294 + (-t332 * t88 + t336 * t80) * t282;
t27 = (t332 * t83 + t336 * t87) * t294 + t282 * (-t332 * t87 + t336 * t83);
t26 = (t332 * t79 + t336 * t86) * t293 + (-t332 * t86 + t336 * t79) * t281;
t25 = (t332 * t82 + t336 * t85) * t293 + t281 * (-t332 * t85 + t336 * t82);
t24 = ((-t138 * t538 - t144 * t254) * t298 - t286 * (t138 * t535 - t144 * t251)) * t304 + ((-t138 * t254 + t144 * t538) * t298 + (t138 * t251 + t144 * t535) * t286) * t292;
t23 = ((-t137 * t539 - t143 * t253) * t297 - t285 * (t137 * t536 - t143 * t250)) * t303 + ((-t137 * t253 + t143 * t539) * t297 + (t137 * t250 + t143 * t536) * t285) * t291;
t22 = ((-t136 * t540 - t142 * t252) * t296 - t284 * (t136 * t537 - t142 * t249)) * t302 + ((-t136 * t252 + t142 * t540) * t296 + (t136 * t249 + t142 * t537) * t284) * t290;
t21 = ((-t135 * t544 - t141 * t247) * t295 - (t135 * t541 - t141 * t244) * t283) * t301 + ((-t135 * t247 + t141 * t544) * t295 + t283 * (t135 * t244 + t141 * t541)) * t289;
t20 = ((-t134 * t545 - t140 * t246) * t294 - (t134 * t542 - t140 * t243) * t282) * t300 + ((-t134 * t246 + t140 * t545) * t294 + t282 * (t134 * t243 + t140 * t542)) * t288;
t19 = ((-t133 * t546 - t139 * t245) * t293 - (t133 * t543 - t139 * t242) * t281) * t299 + ((-t133 * t245 + t139 * t546) * t293 + t281 * (t133 * t242 + t139 * t543)) * t287;
t18 = ((-t138 * t559 + t144 * t189) * t298 + t286 * (t138 * t562 + t144 * t192)) * t304 + ((t138 * t189 + t144 * t559) * t298 + t286 * (t138 * t192 - t144 * t562)) * t292;
t17 = ((-t137 * t560 + t143 * t188) * t297 + t285 * (t137 * t563 + t143 * t191)) * t303 + ((t137 * t188 + t143 * t560) * t297 + t285 * (t137 * t191 - t143 * t563)) * t291;
t16 = ((-t136 * t561 + t142 * t187) * t296 + t284 * (t136 * t564 + t142 * t190)) * t302 + ((t136 * t187 + t142 * t561) * t296 + t284 * (t136 * t190 - t142 * t564)) * t290;
t15 = ((-t135 * t565 + t141 * t183) * t295 + (t135 * t568 + t141 * t186) * t283) * t301 + ((t135 * t183 + t141 * t565) * t295 + (t135 * t186 - t141 * t568) * t283) * t289;
t14 = ((-t134 * t566 + t140 * t182) * t294 + (t134 * t569 + t140 * t185) * t282) * t300 + ((t134 * t182 + t140 * t566) * t294 + (t134 * t185 - t140 * t569) * t282) * t288;
t13 = ((-t133 * t567 + t139 * t181) * t293 + (t133 * t570 + t139 * t184) * t281) * t299 + ((t133 * t181 + t139 * t567) * t293 + (t133 * t184 - t139 * t570) * t281) * t287;
t12 = (t332 * t65 + t336 * t54) * t298 + (-t332 * t54 + t336 * t65) * t286;
t11 = (-t332 * t66 + t336 * t53) * t298 - (t332 * t53 + t336 * t66) * t286;
t10 = (t332 * t63 + t336 * t52) * t297 + (-t332 * t52 + t336 * t63) * t285;
t9 = (-t332 * t64 + t336 * t51) * t297 - t285 * (t332 * t51 + t336 * t64);
t8 = (t332 * t61 + t336 * t50) * t296 + (-t332 * t50 + t336 * t61) * t284;
t7 = (-t332 * t62 + t336 * t49) * t296 - t284 * (t332 * t49 + t336 * t62);
t6 = (t332 * t57 + t336 * t48) * t295 + (-t332 * t48 + t336 * t57) * t283;
t5 = (-t332 * t60 + t336 * t47) * t295 - t283 * (t332 * t47 + t336 * t60);
t4 = (t332 * t56 + t336 * t46) * t294 + (-t332 * t46 + t336 * t56) * t282;
t3 = (-t332 * t59 + t336 * t45) * t294 - t282 * (t332 * t45 + t336 * t59);
t2 = (t332 * t55 + t336 * t44) * t293 + (-t332 * t44 + t336 * t55) * t281;
t1 = (-t332 * t58 + t336 * t43) * t293 - t281 * (t332 * t43 + t336 * t58);
t172 = [((pkin(3) * t405 + t177 * t571) * t380 + (-t177 * pkin(3) + t405 * t571) * t374) * t318 / ((t220 + t464) * t331 + ((t261 - t425) * t374 + t417) * t380 + t465 * t374 + t404 + ((t271 * t374 * t438 + pkin(3) * t426) * t380 + (-t419 + (t426 * t374 + (-0.2e1 * t331 + 0.1e1) * pkin(3) * t438) * t333) * pkin(9)) * t335), ((pkin(3) * t78 + t114 * t571) * t380 + (-pkin(3) * t114 + t571 * t78) * t374) * t42, ((-pkin(3) * t75 - t111 * t571) * t380 - (-pkin(3) * t111 + t571 * t75) * t374) * t42, ((pkin(3) * t18 + t24 * t571) * t380 - (pkin(3) * t24 - t18 * t571) * t374) * t42, ((-pkin(3) * t12 - t36 * t571) * t380 + t374 * (pkin(3) * t36 - t12 * t571)) * t42, ((-pkin(3) * t11 - t35 * t571) * t380 + t374 * (pkin(3) * t35 - t11 * t571)) * t42; ((pkin(3) * t406 + t176 * t571) * t378 + (-t176 * pkin(3) + t406 * t571) * t372) * t317 / ((t219 + t464) * t330 + ((t260 - t425) * t372 + t417) * t378 + t465 * t372 + t404 + ((t271 * t372 * t439 + pkin(3) * t427) * t378 + (-t420 + (t427 * t372 + (-0.2e1 * t330 + 0.1e1) * pkin(3) * t439) * t333) * pkin(9)) * t335), ((pkin(3) * t77 + t113 * t571) * t378 + (-pkin(3) * t113 + t571 * t77) * t372) * t41, ((-pkin(3) * t74 - t110 * t571) * t378 - (-pkin(3) * t110 + t571 * t74) * t372) * t41, ((pkin(3) * t17 + t23 * t571) * t378 - (pkin(3) * t23 - t17 * t571) * t372) * t41, ((-pkin(3) * t10 - t34 * t571) * t378 + t372 * (pkin(3) * t34 - t10 * t571)) * t41, ((-pkin(3) * t9 - t33 * t571) * t378 + (pkin(3) * t33 - t571 * t9) * t372) * t41; ((pkin(3) * t407 + t175 * t571) * t376 + (-t175 * pkin(3) + t407 * t571) * t370) * t316 / ((t218 + t464) * t329 + ((t259 - t425) * t370 + t417) * t376 + t465 * t370 + t404 + ((t271 * t370 * t440 + pkin(3) * t428) * t376 + (-t421 + (t428 * t370 + (-0.2e1 * t329 + 0.1e1) * pkin(3) * t440) * t333) * pkin(9)) * t335), ((pkin(3) * t76 + t112 * t571) * t376 + (-pkin(3) * t112 + t571 * t76) * t370) * t40, ((-pkin(3) * t73 - t109 * t571) * t376 - (-pkin(3) * t109 + t571 * t73) * t370) * t40, ((pkin(3) * t16 + t22 * t571) * t376 - (pkin(3) * t22 - t16 * t571) * t370) * t40, ((-pkin(3) * t8 - t32 * t571) * t376 + t370 * (pkin(3) * t32 - t571 * t8)) * t40, ((-pkin(3) * t7 - t31 * t571) * t376 + (pkin(3) * t31 - t571 * t7) * t370) * t40; ((pkin(3) * t408 + t171 * t571) * t362 + (-t171 * pkin(3) + t408 * t571) * t356) * t315 / ((t214 + t464) * t328 + ((t258 - t425) * t356 + t417) * t362 + t465 * t356 + t404 + ((t271 * t356 * t441 + pkin(3) * t429) * t362 + (-t422 + (t429 * t356 + (-0.2e1 * t328 + 0.1e1) * pkin(3) * t441) * t333) * pkin(9)) * t335), ((pkin(3) * t72 + t108 * t571) * t362 + (-pkin(3) * t108 + t571 * t72) * t356) * t39, ((-pkin(3) * t69 - t105 * t571) * t362 - (-pkin(3) * t105 + t571 * t69) * t356) * t39, ((pkin(3) * t15 + t21 * t571) * t362 - (pkin(3) * t21 - t15 * t571) * t356) * t39, ((-pkin(3) * t6 - t30 * t571) * t362 + (pkin(3) * t30 - t571 * t6) * t356) * t39, ((-pkin(3) * t5 - t29 * t571) * t362 + (pkin(3) * t29 - t5 * t571) * t356) * t39; ((pkin(3) * t409 + t170 * t571) * t360 + (-t170 * pkin(3) + t409 * t571) * t354) * t314 / ((t213 + t464) * t327 + ((t257 - t425) * t354 + t417) * t360 + t465 * t354 + t404 + ((t271 * t354 * t442 + pkin(3) * t430) * t360 + (-t423 + (t430 * t354 + (-0.2e1 * t327 + 0.1e1) * pkin(3) * t442) * t333) * pkin(9)) * t335), ((pkin(3) * t71 + t107 * t571) * t360 + (-pkin(3) * t107 + t571 * t71) * t354) * t38, ((-pkin(3) * t68 - t104 * t571) * t360 - (-pkin(3) * t104 + t571 * t68) * t354) * t38, ((pkin(3) * t14 + t20 * t571) * t360 - (pkin(3) * t20 - t14 * t571) * t354) * t38, ((-pkin(3) * t4 - t28 * t571) * t360 + (pkin(3) * t28 - t4 * t571) * t354) * t38, ((-pkin(3) * t3 - t27 * t571) * t360 + (pkin(3) * t27 - t3 * t571) * t354) * t38; ((pkin(3) * t410 + t169 * t571) * t358 + (-t169 * pkin(3) + t410 * t571) * t352) * t313 / ((t212 + t464) * t326 + ((t256 - t425) * t352 + t417) * t358 + t465 * t352 + t404 + ((t271 * t352 * t443 + pkin(3) * t431) * t358 + (-t424 + (t431 * t352 + (-0.2e1 * t326 + 0.1e1) * pkin(3) * t443) * t333) * pkin(9)) * t335), ((pkin(3) * t70 + t106 * t571) * t358 + (-pkin(3) * t106 + t571 * t70) * t352) * t37, ((-pkin(3) * t67 - t103 * t571) * t358 - (-pkin(3) * t103 + t571 * t67) * t352) * t37, ((pkin(3) * t13 + t19 * t571) * t358 - (pkin(3) * t19 - t13 * t571) * t352) * t37, ((-pkin(3) * t2 - t26 * t571) * t358 + (pkin(3) * t26 - t2 * t571) * t352) * t37, ((-pkin(3) * t1 - t25 * t571) * t358 + (pkin(3) * t25 - t1 * t571) * t352) * t37;];
Jinv  = t172;
