% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3PRRR1G3P3A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3*3x8]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:07
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRR1G3P3A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G3P3A0_inertia_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G3P3A0_inertia_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G3P3A0_inertia_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G3P3A0_inertia_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G3P3A0_inertia_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:06:55
% EndTime: 2020-03-09 21:06:57
% DurationCPUTime: 1.50s
% Computational Cost: add. (4676->189), mult. (2019->409), div. (981->9), fcn. (2856->24), ass. (0->211)
t367 = legFrame(3,2);
t358 = sin(t367);
t364 = pkin(7) + qJ(2,3);
t355 = qJ(3,3) + t364;
t337 = sin(t355);
t340 = cos(t355);
t343 = sin(t364);
t346 = cos(t364);
t325 = t337 * t346 - t343 * t340;
t502 = 0.1e1 / t325;
t508 = t358 * t502;
t368 = legFrame(2,2);
t359 = sin(t368);
t365 = pkin(7) + qJ(2,2);
t356 = qJ(3,2) + t365;
t338 = sin(t356);
t341 = cos(t356);
t344 = sin(t365);
t347 = cos(t365);
t326 = t338 * t347 - t344 * t341;
t501 = 0.1e1 / t326;
t507 = t359 * t501;
t369 = legFrame(1,2);
t360 = sin(t369);
t366 = pkin(7) + qJ(2,1);
t357 = qJ(3,1) + t366;
t339 = sin(t357);
t342 = cos(t357);
t345 = sin(t366);
t348 = cos(t366);
t327 = t339 * t348 - t345 * t342;
t500 = 0.1e1 / t327;
t506 = t360 * t500;
t505 = t500 * t339;
t504 = t501 * t338;
t503 = t502 * t337;
t499 = t502 ^ 2;
t328 = pkin(2) * t343 + pkin(3) * t337;
t498 = t502 * t328;
t377 = 0.1e1 / pkin(2);
t496 = t502 * t377;
t317 = 0.1e1 / t325 ^ 2;
t495 = t317 * t328;
t494 = t317 * t337;
t493 = t501 ^ 2;
t329 = pkin(2) * t344 + pkin(3) * t338;
t492 = t501 * t329;
t490 = t501 * t377;
t319 = 0.1e1 / t326 ^ 2;
t489 = t319 * t329;
t488 = t319 * t338;
t487 = t500 ^ 2;
t330 = pkin(2) * t345 + pkin(3) * t339;
t486 = t500 * t330;
t484 = t500 * t377;
t321 = 0.1e1 / t327 ^ 2;
t483 = t321 * t330;
t482 = t321 * t339;
t361 = cos(t367);
t480 = t502 * t361;
t362 = cos(t368);
t478 = t501 * t362;
t363 = cos(t369);
t476 = t500 * t363;
t331 = pkin(2) * t346 + pkin(3) * t340;
t475 = t331 * t340;
t474 = t331 * t358;
t332 = pkin(2) * t347 + pkin(3) * t341;
t473 = t332 * t341;
t472 = t332 * t359;
t333 = pkin(2) * t348 + pkin(3) * t342;
t471 = t333 * t342;
t470 = t333 * t360;
t469 = t340 * t361;
t468 = t341 * t362;
t467 = t342 * t363;
t466 = t358 * t361;
t465 = t359 * t362;
t464 = t360 * t363;
t376 = 0.1e1 / pkin(3);
t463 = t376 * t377;
t462 = t502 * t474;
t461 = t502 * t469;
t460 = t340 * t496;
t370 = sin(qJ(3,3));
t459 = t370 * t495;
t373 = cos(qJ(3,3));
t458 = t373 * t495;
t457 = t501 * t472;
t456 = t501 * t468;
t455 = t341 * t490;
t371 = sin(qJ(3,2));
t454 = t371 * t489;
t374 = cos(qJ(3,2));
t453 = t374 * t489;
t452 = t500 * t470;
t451 = t500 * t467;
t450 = t342 * t484;
t372 = sin(qJ(3,1));
t449 = t372 * t483;
t375 = cos(qJ(3,1));
t448 = t375 * t483;
t447 = t331 * t480;
t446 = t370 * t503;
t445 = t373 * t503;
t444 = t340 * t508;
t443 = t332 * t478;
t442 = t371 * t504;
t441 = t374 * t504;
t440 = t341 * t507;
t439 = t333 * t476;
t438 = t372 * t505;
t437 = t375 * t505;
t436 = t342 * t506;
t435 = t331 * t463;
t434 = t332 * t463;
t433 = t333 * t463;
t334 = t340 ^ 2;
t432 = t334 * t466;
t335 = t341 ^ 2;
t431 = t335 * t465;
t336 = t342 ^ 2;
t430 = t336 * t464;
t429 = t337 * t496;
t428 = t338 * t490;
t427 = t339 * t484;
t426 = t480 * t503;
t425 = t370 * t461;
t424 = t373 * t461;
t423 = t340 * t459;
t422 = t340 * t458;
t421 = t474 * t494;
t349 = t358 ^ 2;
t420 = t317 * t349 * t475;
t419 = t478 * t504;
t418 = t371 * t456;
t417 = t374 * t456;
t416 = t341 * t454;
t415 = t341 * t453;
t414 = t472 * t488;
t350 = t359 ^ 2;
t413 = t319 * t350 * t473;
t412 = t476 * t505;
t411 = t372 * t451;
t410 = t375 * t451;
t409 = t342 * t449;
t408 = t342 * t448;
t407 = t470 * t482;
t351 = t360 ^ 2;
t406 = t321 * t351 * t471;
t405 = t370 * t444;
t404 = t373 * t444;
t403 = t371 * t440;
t402 = t374 * t440;
t401 = t372 * t436;
t400 = t375 * t436;
t399 = t466 * t475;
t398 = t465 * t473;
t397 = t464 * t471;
t396 = t358 * t460;
t395 = t359 * t455;
t394 = t360 * t450;
t393 = t361 * t460;
t392 = t362 * t455;
t391 = t363 * t450;
t390 = t331 * t426;
t352 = t361 ^ 2;
t389 = t352 * t475 * t499;
t388 = t317 * t399;
t387 = t332 * t419;
t353 = t362 ^ 2;
t386 = t353 * t473 * t493;
t385 = t319 * t398;
t384 = t333 * t412;
t354 = t363 ^ 2;
t383 = t354 * t471 * t487;
t382 = t321 * t397;
t381 = t399 * t499;
t380 = t398 * t493;
t379 = t397 * t487;
t378 = 0.1e1 / pkin(2) ^ 2;
t315 = t464 + t465 + t466;
t314 = t463 * t486;
t313 = t463 * t492;
t312 = t463 * t498;
t311 = t433 * t506;
t310 = t434 * t507;
t309 = t435 * t508;
t308 = t433 * t476;
t307 = t434 * t478;
t306 = t435 * t480;
t305 = t314 - 0.2e1 * t427;
t304 = t314 - t427;
t303 = t313 - 0.2e1 * t428;
t302 = t313 - t428;
t301 = t312 - 0.2e1 * t429;
t300 = t312 - t429;
t299 = -t308 + t391;
t298 = -t308 + 0.2e1 * t391;
t297 = -t307 + t392;
t296 = -t307 + 0.2e1 * t392;
t295 = -t306 + t393;
t294 = -t306 + 0.2e1 * t393;
t293 = t311 - 0.2e1 * t394;
t292 = t311 - t394;
t291 = t310 - 0.2e1 * t395;
t290 = t310 - t395;
t289 = t309 - 0.2e1 * t396;
t288 = t309 - t396;
t287 = (t436 * t505 + t440 * t504 + t444 * t503) * t378;
t1 = [t351 + t350 + t349, (t317 * t334 * t352 + t319 * t335 * t353 + t321 * t336 * t354) * t378, 0, 0, (t295 * t461 + t297 * t456 + t299 * t451 + (-t295 * t447 - t297 * t443 - t299 * t439) * t376) * t377, t294 * t424 + t296 * t417 + t298 * t410 + (-t373 * t389 - t374 * t386 - t375 * t383) * t463, -t294 * t425 - t296 * t418 - t298 * t411 + (t370 * t389 + t371 * t386 + t372 * t383) * t463, 1; t315, (-t317 * t432 - t319 * t431 - t321 * t430) * t378, 0, 0, (t288 * t461 + t290 * t456 + t292 * t451 + (-t288 * t447 - t290 * t443 - t292 * t439) * t376) * t377, t289 * t424 + t291 * t417 + t293 * t410 + (t373 * t381 + t374 * t380 + t375 * t379) * t463, -t289 * t425 - t291 * t418 - t293 * t411 + (-t370 * t381 - t371 * t380 - t372 * t379) * t463, 0; 0, (-t467 * t482 - t468 * t488 - t469 * t494) * t378, 0, 0, (t300 * t461 + t302 * t456 + t304 * t451 + (-t300 * t447 - t302 * t443 - t304 * t439) * t376) * t377, t301 * t424 + t303 * t417 + t305 * t410 + (t373 * t390 + t374 * t387 + t375 * t384) * t463, -t301 * t425 - t303 * t418 - t305 * t411 + (-t370 * t390 - t371 * t387 - t372 * t384) * t463, 0; t315, (-t430 * t487 - t431 * t493 - t432 * t499) * t378, 0, 0, (-t295 * t444 - t297 * t440 - t299 * t436 + (t295 * t462 + t297 * t457 + t299 * t452) * t376) * t377, -t294 * t404 - t296 * t402 - t298 * t400 + (t373 * t388 + t374 * t385 + t375 * t382) * t463, t294 * t405 + t296 * t403 + t298 * t401 + (-t370 * t388 - t371 * t385 - t372 * t382) * t463, 0; t354 + t353 + t352, (t334 * t349 * t499 + t335 * t350 * t493 + t336 * t351 * t487) * t378, 0, 0, (-t288 * t444 - t290 * t440 - t292 * t436 + (t288 * t462 + t290 * t457 + t292 * t452) * t376) * t377, -t289 * t404 - t291 * t402 - t293 * t400 + (-t373 * t420 - t374 * t413 - t375 * t406) * t463, t289 * t405 + t291 * t403 + t293 * t401 + (t370 * t420 + t371 * t413 + t372 * t406) * t463, 1; 0, t287, 0, 0, (-t300 * t444 - t302 * t440 - t304 * t436 + (t300 * t462 + t302 * t457 + t304 * t452) * t376) * t377, -t301 * t404 - t303 * t402 - t305 * t400 + (-t373 * t421 - t374 * t414 - t375 * t407) * t463, t301 * t405 + t303 * t403 + t305 * t401 + (t370 * t421 + t371 * t414 + t372 * t407) * t463, 0; 0, (-t340 * t426 - t341 * t419 - t342 * t412) * t378, 0, 0, (-t295 * t503 - t297 * t504 - t299 * t505 + (t295 * t498 + t297 * t492 + t299 * t486) * t376) * t377, -t294 * t445 - t296 * t441 - t298 * t437 + (t361 * t422 + t362 * t415 + t363 * t408) * t463, t294 * t446 + t296 * t442 + t298 * t438 + (-t361 * t423 - t362 * t416 - t363 * t409) * t463, 0; 0, t287, 0, 0, (-t288 * t503 - t290 * t504 - t292 * t505 + (t288 * t498 + t290 * t492 + t292 * t486) * t376) * t377, -t289 * t445 - t291 * t441 - t293 * t437 + (-t358 * t422 - t359 * t415 - t360 * t408) * t463, t289 * t446 + t291 * t442 + t293 * t438 + (t358 * t423 + t359 * t416 + t360 * t409) * t463, 0; 0, (t337 ^ 2 * t499 + t338 ^ 2 * t493 + t339 ^ 2 * t487) * t378, 0, 0, (-t300 * t503 - t302 * t504 - t304 * t505 + (t300 * t498 + t302 * t492 + t304 * t486) * t376) * t377, -t301 * t445 - t303 * t441 - t305 * t437 + (-t337 * t458 - t338 * t453 - t339 * t448) * t463, t301 * t446 + t303 * t442 + t305 * t438 + (t337 * t459 + t338 * t454 + t339 * t449) * t463, 1;];
tau_reg  = t1;
