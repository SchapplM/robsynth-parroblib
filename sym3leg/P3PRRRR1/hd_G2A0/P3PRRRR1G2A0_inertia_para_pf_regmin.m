% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3PRRRR1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3*3x12]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:16
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRRR1G2A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G2A0_inertia_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G2A0_inertia_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G2A0_inertia_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G2A0_inertia_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G2A0_inertia_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:16:33
% EndTime: 2020-03-09 21:16:35
% DurationCPUTime: 2.24s
% Computational Cost: add. (526->210), mult. (2592->540), div. (1215->17), fcn. (2760->18), ass. (0->270)
t367 = sin(qJ(3,1));
t339 = t367 ^ 2;
t373 = cos(qJ(3,1));
t391 = t373 ^ 2;
t355 = 0.1e1 / t391;
t368 = sin(qJ(2,1));
t375 = 0.1e1 / pkin(2);
t340 = 0.1e1 / t368;
t374 = cos(qJ(2,1));
t359 = t374 ^ 2;
t530 = t340 * t359;
t570 = t375 * (t339 * t355 * t530 - t368);
t365 = sin(qJ(3,2));
t334 = t365 ^ 2;
t371 = cos(qJ(3,2));
t387 = t371 ^ 2;
t349 = 0.1e1 / t387;
t366 = sin(qJ(2,2));
t335 = 0.1e1 / t366;
t372 = cos(qJ(2,2));
t353 = t372 ^ 2;
t536 = t335 * t353;
t569 = t375 * (t334 * t349 * t536 - t366);
t363 = sin(qJ(3,3));
t329 = t363 ^ 2;
t369 = cos(qJ(3,3));
t383 = t369 ^ 2;
t343 = 0.1e1 / t383;
t364 = sin(qJ(2,3));
t330 = 0.1e1 / t364;
t370 = cos(qJ(2,3));
t347 = t370 ^ 2;
t542 = t330 * t347;
t568 = t375 * (t329 * t343 * t542 - t364);
t354 = 0.1e1 / t373;
t512 = t354 * t367;
t567 = (t368 + t530) * t375 * t512;
t348 = 0.1e1 / t371;
t519 = t348 * t365;
t566 = (t366 + t536) * t375 * t519;
t342 = 0.1e1 / t369;
t526 = t342 * t363;
t565 = (t364 + t542) * t375 * t526;
t341 = 0.1e1 / t368 ^ 2;
t507 = t359 * t341;
t510 = t355 * t367;
t564 = (0.1e1 + t507) * t510;
t336 = 0.1e1 / t366 ^ 2;
t514 = t353 * t336;
t517 = t349 * t365;
t563 = (0.1e1 + t514) * t517;
t331 = 0.1e1 / t364 ^ 2;
t521 = t347 * t331;
t524 = t343 * t363;
t562 = (0.1e1 + t521) * t524;
t561 = 0.2e1 * t375;
t376 = 1 / (pkin(2) ^ 2);
t560 = 2 * t376;
t360 = legFrame(3,2);
t321 = sin(t360);
t504 = t364 * t369;
t324 = cos(t360);
t547 = t324 * t363;
t309 = t321 * t504 - t547;
t559 = t309 * t324;
t361 = legFrame(2,2);
t322 = sin(t361);
t502 = t366 * t371;
t325 = cos(t361);
t546 = t325 * t365;
t310 = t322 * t502 - t546;
t558 = t310 * t325;
t362 = legFrame(1,2);
t323 = sin(t362);
t500 = t368 * t373;
t326 = cos(t362);
t545 = t326 * t367;
t311 = t323 * t500 - t545;
t557 = t311 * t326;
t552 = t321 * t363;
t312 = t324 * t504 + t552;
t556 = t312 * t321;
t550 = t322 * t365;
t313 = t325 * t502 + t550;
t555 = t313 * t322;
t548 = t323 * t367;
t314 = t326 * t500 + t548;
t554 = t314 * t323;
t553 = t321 * t324;
t551 = t322 * t325;
t549 = t323 * t326;
t544 = t330 * t342;
t543 = t330 * t343;
t541 = t330 * t370;
t540 = t331 * t343;
t539 = t331 * t370;
t538 = t335 * t348;
t537 = t335 * t349;
t535 = t335 * t372;
t534 = t336 * t349;
t533 = t336 * t372;
t532 = t340 * t354;
t531 = t340 * t355;
t529 = t340 * t374;
t528 = t341 * t355;
t527 = t341 * t374;
t525 = t342 * t370;
t344 = t342 * t343;
t523 = t344 * t370;
t346 = t370 * t347;
t522 = t346 * t363;
t520 = t347 * t363;
t518 = t348 * t372;
t350 = t348 * t349;
t516 = t350 * t372;
t352 = t372 * t353;
t515 = t352 * t365;
t513 = t353 * t365;
t511 = t354 * t374;
t356 = t354 * t355;
t509 = t356 * t374;
t358 = t374 * t359;
t508 = t358 * t367;
t506 = t359 * t367;
t505 = t363 * t370;
t503 = t365 * t372;
t501 = t367 * t374;
t499 = t309 * t544;
t498 = t310 * t538;
t497 = t311 * t532;
t496 = t312 * t544;
t495 = t313 * t538;
t494 = t314 * t532;
t493 = t329 * t540;
t492 = t330 * t524;
t491 = t331 * t525;
t490 = t344 * t521;
t489 = 0.1e1 / t383 ^ 2 * t521;
t488 = t331 * t520;
t487 = t334 * t534;
t486 = t335 * t517;
t485 = t336 * t518;
t484 = t350 * t514;
t483 = 0.1e1 / t387 ^ 2 * t514;
t482 = t336 * t513;
t481 = t339 * t528;
t480 = t340 * t510;
t479 = t341 * t511;
t478 = t356 * t507;
t477 = 0.1e1 / t391 ^ 2 * t507;
t476 = t341 * t506;
t474 = t343 * t505;
t473 = t344 * t505;
t471 = t349 * t503;
t470 = t350 * t503;
t468 = t355 * t501;
t467 = t356 * t501;
t420 = t341 * t375 * t468;
t421 = t336 * t375 * t471;
t422 = t331 * t375 * t474;
t466 = t309 * t422 + t310 * t421 + t311 * t420;
t465 = t312 * t422 + t313 * t421 + t314 * t420;
t460 = t331 * t346 + t370;
t458 = t336 * t352 + t372;
t456 = t341 * t358 + t374;
t327 = t329 ^ 2;
t455 = t327 * t489;
t328 = t363 * t329;
t454 = t328 * t490;
t453 = t328 * t331 * t523;
t452 = t329 * t330 * t523;
t451 = t370 * t493;
t450 = t329 * t489;
t449 = t330 * t474;
t448 = t330 * t473;
t447 = t344 * t488;
t446 = t331 * t473;
t332 = t334 ^ 2;
t445 = t332 * t483;
t333 = t365 * t334;
t444 = t333 * t484;
t443 = t333 * t336 * t516;
t442 = t334 * t335 * t516;
t441 = t372 * t487;
t440 = t334 * t483;
t439 = t335 * t471;
t438 = t335 * t470;
t437 = t350 * t482;
t436 = t336 * t470;
t337 = t339 ^ 2;
t435 = t337 * t477;
t338 = t367 * t339;
t434 = t338 * t478;
t433 = t338 * t341 * t509;
t432 = t339 * t340 * t509;
t431 = t374 * t481;
t430 = t339 * t477;
t429 = t340 * t468;
t428 = t340 * t467;
t427 = t356 * t476;
t426 = t341 * t467;
t425 = t489 * t553;
t424 = t483 * t551;
t423 = t477 * t549;
t419 = t309 * t321 - t312 * t324;
t418 = t310 * t322 - t313 * t325;
t417 = t311 * t323 - t314 * t326;
t415 = t329 * t490 - t342;
t414 = t346 * t493 - t370;
t412 = t334 * t484 - t348;
t411 = t352 * t487 - t372;
t409 = t339 * t478 - t354;
t408 = t358 * t481 - t374;
t407 = t321 * t562;
t406 = t322 * t563;
t405 = t323 * t564;
t404 = t324 * t562;
t403 = t325 * t563;
t402 = t326 * t564;
t401 = t321 * t415;
t400 = t323 * t409;
t399 = t324 * t415;
t398 = t326 * t409;
t397 = t412 * t325;
t396 = t412 * t322;
t395 = t507 + t514 + t521;
t273 = t309 * t491 + t310 * t485 + t311 * t479;
t274 = t312 * t491 + t313 * t485 + t314 * t479;
t320 = t326 ^ 2;
t319 = t325 ^ 2;
t318 = t324 ^ 2;
t317 = t323 ^ 2;
t316 = t322 ^ 2;
t315 = t321 ^ 2;
t308 = t326 * t570;
t307 = t323 * t570;
t306 = t325 * t569;
t305 = t322 * t569;
t304 = t324 * t568;
t303 = t321 * t568;
t296 = t326 * t567;
t295 = t323 * t567;
t294 = t325 * t566;
t293 = t322 * t566;
t292 = t324 * t565;
t291 = t321 * t565;
t290 = (t324 * t544 + t325 * t538 + t326 * t532) * t376;
t289 = (-t321 * t544 - t322 * t538 - t323 * t532) * t376;
t288 = (-t343 * t553 - t349 * t551 - t355 * t549) * t376;
t287 = (t324 * t492 + t325 * t486 + t326 * t480) * t376;
t286 = (-t321 * t492 - t322 * t486 - t323 * t480) * t376;
t285 = (-t324 * t446 - t325 * t436 - t326 * t426) * t376;
t284 = (-t324 * t453 - t325 * t443 - t326 * t433) * t376;
t283 = (t321 * t446 + t322 * t436 + t323 * t426) * t376;
t282 = (t321 * t453 + t322 * t443 + t323 * t433) * t376;
t281 = (-t324 * t451 - t325 * t441 - t326 * t431) * t560;
t280 = (t321 * t451 + t322 * t441 + t323 * t431) * t560;
t279 = (-t329 * t425 - t334 * t424 - t339 * t423) * t376;
t278 = (-t327 * t425 - t332 * t424 - t337 * t423) * t376;
t277 = (t429 * t549 + t439 * t551 + t449 * t553) * t560;
t276 = (t432 * t549 + t442 * t551 + t452 * t553) * t560;
t275 = (-t434 * t549 - t444 * t551 - t454 * t553) * t560;
t272 = ((t323 * t506 + t314) * t531 + (t322 * t513 + t313) * t537 + (t321 * t520 + t312) * t543) * t375;
t271 = ((-t326 * t506 + t311) * t531 + (-t325 * t513 + t310) * t537 + (-t324 * t520 + t309) * t543) * t375;
t270 = t309 * t312 * t540 + t310 * t313 * t534 + t311 * t314 * t528;
t269 = ((-t314 * t374 - t323 * t508) * t528 + (-t313 * t372 - t322 * t515) * t534 + (-t312 * t370 - t321 * t522) * t540) * t375;
t268 = ((-t311 * t374 + t326 * t508) * t528 + (-t310 * t372 + t325 * t515) * t534 + (-t309 * t370 + t324 * t522) * t540) * t375;
t267 = (t417 * t428 + t418 * t438 + t419 * t448) * t375;
t266 = (-t417 * t427 - t418 * t437 - t419 * t447) * t375;
t1 = [t309 ^ 2 * t540 + t310 ^ 2 * t534 + t311 ^ 2 * t528, (t318 * t450 + t319 * t440 + t320 * t430) * t376, (t427 * t557 + t437 * t558 + t447 * t559) * t561, (-t428 * t557 - t438 * t558 - t448 * t559) * t561, (t318 * t455 + t319 * t445 + t320 * t435) * t376, (t318 * t454 + t319 * t444 + t320 * t434) * t560, (-t318 * t452 - t319 * t442 - t320 * t432) * t560, (-t318 * t449 - t319 * t439 - t320 * t429) * t560, (t318 * t343 + t319 * t349 + t320 * t355) * t376, t292 * t499 + t294 * t498 + t296 * t497 + (t309 * t404 + t310 * t403 + t311 * t402) * t375, -t304 * t499 - t306 * t498 - t308 * t497 + (-t309 * t399 - t310 * t397 - t311 * t398) * t375, 1; t270, t279, t266, t267, t278, t275, t276, t277, t288, -t291 * t499 - t293 * t498 - t295 * t497 + (t312 * t404 + t313 * t403 + t314 * t402) * t375, t303 * t499 + t305 * t498 + t307 * t497 + (-t312 * t399 - t313 * t397 - t314 * t398) * t375, 0; t273, t285, t268, t271, t284, t281, t287, t290, 0, ((-t311 * t527 + t456 * t545) * t354 + (-t310 * t533 + t458 * t546) * t348 + (-t309 * t539 + t460 * t547) * t342) * t375, (-t324 * t414 - t325 * t411 - t326 * t408) * t375 + t466, 0; t270, t279, t266, t267, t278, t275, t276, t277, t288, t292 * t496 + t294 * t495 + t296 * t494 + (-t309 * t407 - t310 * t406 - t311 * t405) * t375, -t304 * t496 - t306 * t495 - t308 * t494 + (t309 * t401 + t310 * t396 + t311 * t400) * t375, 0; t312 ^ 2 * t540 + t313 ^ 2 * t534 + t314 ^ 2 * t528, (t315 * t450 + t316 * t440 + t317 * t430) * t376, (-t427 * t554 - t437 * t555 - t447 * t556) * t561, (t428 * t554 + t438 * t555 + t448 * t556) * t561, (t315 * t455 + t316 * t445 + t317 * t435) * t376, (t315 * t454 + t316 * t444 + t317 * t434) * t560, (-t315 * t452 - t316 * t442 - t317 * t432) * t560, (-t315 * t449 - t316 * t439 - t317 * t429) * t560, (t315 * t343 + t316 * t349 + t317 * t355) * t376, -t291 * t496 - t293 * t495 - t295 * t494 + (-t312 * t407 - t313 * t406 - t314 * t405) * t375, t303 * t496 + t305 * t495 + t307 * t494 + (t312 * t401 + t313 * t396 + t314 * t400) * t375, 1; t274, t283, t269, t272, t282, t280, t286, t289, 0, ((-t314 * t527 - t456 * t548) * t354 + (-t313 * t533 - t458 * t550) * t348 + (-t312 * t539 - t460 * t552) * t342) * t375, (t321 * t414 + t322 * t411 + t323 * t408) * t375 + t465, 0; t273, t285, t268, t271, t284, t281, t287, t290, 0, -t273 * t375 + t292 * t541 + t294 * t535 + t296 * t529, -t304 * t541 - t306 * t535 - t308 * t529 + t466, 0; t274, t283, t269, t272, t282, t280, t286, t289, 0, -t274 * t375 - t291 * t541 - t293 * t535 - t295 * t529, t303 * t541 + t305 * t535 + t307 * t529 + t465, 0; t395, (t528 + t534 + t540) * t376, (-t342 * t521 - t348 * t514 - t354 * t507) * t561, (t330 * t525 + t335 * t518 + t340 * t511) * t561, (t481 + t487 + t493) * t376, (t331 * t526 + t336 * t519 + t341 * t512) * t560, 0, 0, 0, -0.2e1 * t395 * t375, (t342 * t488 + t348 * t482 + t354 * t476) * t561, 1;];
tau_reg  = t1;
