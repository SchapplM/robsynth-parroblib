% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3PRRRR1G1P3A0
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
% Datum: 2020-03-09 20:34
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRRR1G1P3A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G1P3A0_inertia_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G1P3A0_inertia_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G1P3A0_inertia_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G1P3A0_inertia_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G1P3A0_inertia_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 20:34:10
% EndTime: 2020-03-09 20:34:12
% DurationCPUTime: 1.22s
% Computational Cost: add. (608->161), mult. (2112->405), div. (840->17), fcn. (2742->18), ass. (0->213)
t364 = legFrame(3,3);
t337 = sin(t364);
t340 = cos(t364);
t373 = cos(qJ(3,3));
t367 = sin(qJ(3,3));
t374 = cos(qJ(2,3));
t464 = t367 * t374;
t319 = -t337 * t464 - t340 * t373;
t461 = t373 * t374;
t494 = t337 * t367;
t322 = t340 * t461 + t494;
t511 = t319 * t322;
t323 = -t373 * t337 + t340 * t464;
t510 = t319 * t323;
t509 = t319 * t337;
t365 = legFrame(2,3);
t338 = sin(t365);
t341 = cos(t365);
t375 = cos(qJ(3,2));
t369 = sin(qJ(3,2));
t376 = cos(qJ(2,2));
t463 = t369 * t376;
t320 = -t338 * t463 - t341 * t375;
t460 = t375 * t376;
t493 = t338 * t369;
t324 = t341 * t460 + t493;
t508 = t320 * t324;
t325 = -t375 * t338 + t341 * t463;
t507 = t320 * t325;
t506 = t320 * t338;
t366 = legFrame(1,3);
t339 = sin(t366);
t342 = cos(t366);
t377 = cos(qJ(3,1));
t371 = sin(qJ(3,1));
t378 = cos(qJ(2,1));
t462 = t371 * t378;
t321 = -t339 * t462 - t342 * t377;
t459 = t377 * t378;
t492 = t339 * t371;
t326 = t342 * t459 + t492;
t505 = t321 * t326;
t327 = -t377 * t339 + t342 * t462;
t504 = t321 * t327;
t503 = t321 * t339;
t491 = t340 * t367;
t328 = t337 * t461 - t491;
t502 = t323 * t328;
t501 = t323 * t340;
t490 = t341 * t369;
t329 = t338 * t460 - t490;
t500 = t325 * t329;
t499 = t325 * t341;
t489 = t342 * t371;
t330 = t339 * t459 - t489;
t498 = t327 * t330;
t497 = t327 * t342;
t352 = 0.1e1 / t373;
t356 = 0.1e1 / t375;
t360 = 0.1e1 / t377;
t381 = t373 ^ 2;
t353 = 0.1e1 / t381;
t384 = t375 ^ 2;
t357 = 0.1e1 / t384;
t387 = t377 ^ 2;
t361 = 0.1e1 / t387;
t379 = 1 / pkin(2);
t496 = 2 * t379;
t380 = 1 / (pkin(2) ^ 2);
t495 = 2 * t380;
t368 = sin(qJ(2,3));
t344 = 0.1e1 / t368;
t488 = t344 * t352;
t487 = t344 * t353;
t354 = t352 * t353;
t486 = t344 * t354;
t485 = t344 * t374;
t345 = 0.1e1 / t368 ^ 2;
t484 = t345 * t353;
t483 = t345 / t381 ^ 2;
t482 = t345 * t374;
t370 = sin(qJ(2,2));
t347 = 0.1e1 / t370;
t481 = t347 * t356;
t480 = t347 * t357;
t358 = t356 * t357;
t479 = t347 * t358;
t478 = t347 * t376;
t348 = 0.1e1 / t370 ^ 2;
t477 = t348 * t357;
t476 = t348 / t384 ^ 2;
t475 = t348 * t376;
t372 = sin(qJ(2,1));
t350 = 0.1e1 / t372;
t474 = t350 * t360;
t473 = t350 * t361;
t362 = t360 * t361;
t472 = t350 * t362;
t471 = t350 * t378;
t351 = 0.1e1 / t372 ^ 2;
t470 = t351 * t361;
t469 = t351 / t387 ^ 2;
t468 = t351 * t378;
t467 = t354 * t367;
t466 = t358 * t369;
t465 = t362 * t371;
t431 = t352 * t367 * t368;
t416 = t379 * t431;
t425 = t352 * t379 * t485;
t295 = t323 * t425 + t340 * t416;
t430 = t356 * t369 * t370;
t415 = t379 * t430;
t422 = t356 * t379 * t478;
t297 = t325 * t422 + t341 * t415;
t429 = t360 * t371 * t372;
t414 = t379 * t429;
t419 = t360 * t379 * t471;
t299 = t327 * t419 + t342 * t414;
t343 = t367 ^ 2;
t458 = t343 * t483;
t307 = t319 * t425;
t457 = (-t337 * t416 + t307) * t488;
t456 = t295 * t488;
t453 = t353 * t485;
t424 = t367 * t453;
t398 = -t319 * t424 - t337 * t368;
t455 = t398 * t379 * t488;
t334 = t368 * t379 * t340;
t413 = t323 * t424;
t454 = (-t379 * t413 + t334) * t488;
t452 = t344 * t467;
t451 = t345 * t467;
t450 = t354 * t482;
t346 = t369 ^ 2;
t449 = t346 * t476;
t309 = t320 * t422;
t448 = (-t338 * t415 + t309) * t481;
t447 = t297 * t481;
t444 = t357 * t478;
t421 = t369 * t444;
t396 = -t320 * t421 - t338 * t370;
t446 = t396 * t379 * t481;
t335 = t370 * t379 * t341;
t412 = t325 * t421;
t445 = (-t379 * t412 + t335) * t481;
t443 = t347 * t466;
t442 = t348 * t466;
t441 = t358 * t475;
t349 = t371 ^ 2;
t440 = t349 * t469;
t311 = t321 * t419;
t439 = (-t339 * t414 + t311) * t474;
t438 = t299 * t474;
t435 = t361 * t471;
t418 = t371 * t435;
t394 = -t321 * t418 - t339 * t372;
t437 = t394 * t379 * t474;
t336 = t372 * t379 * t342;
t411 = t327 * t418;
t436 = (-t379 * t411 + t336) * t474;
t434 = t350 * t465;
t433 = t351 * t465;
t432 = t362 * t468;
t428 = t483 * t510;
t427 = t476 * t507;
t426 = t469 * t504;
t423 = t367 * t450;
t420 = t369 * t441;
t417 = t371 * t432;
t410 = t319 * t328 + t322 * t323;
t409 = t320 * t329 + t324 * t325;
t408 = t321 * t330 + t326 * t327;
t407 = t344 * (-t319 * t340 + t323 * t337);
t406 = t347 * (-t320 * t341 + t325 * t338);
t405 = t350 * (-t321 * t342 + t327 * t339);
t404 = (t319 * t482 - t494) * t353;
t403 = (t320 * t475 - t493) * t357;
t402 = (t321 * t468 - t492) * t361;
t401 = (t323 * t482 + t491) * t353;
t400 = (t325 * t475 + t490) * t357;
t399 = (t327 * t468 + t489) * t361;
t397 = -t319 * t423 - t337 * t352;
t395 = -t320 * t420 - t338 * t356;
t393 = -t321 * t417 - t339 * t360;
t392 = -t323 * t423 + t340 * t352;
t391 = -t325 * t420 + t341 * t356;
t390 = -t327 * t417 + t342 * t360;
t318 = t327 ^ 2;
t317 = t325 ^ 2;
t316 = t323 ^ 2;
t315 = t321 ^ 2;
t314 = t320 ^ 2;
t313 = t319 ^ 2;
t306 = (-t337 * t340 * t353 - t338 * t341 * t357 - t339 * t342 * t361) * t380;
t293 = (-t323 * t353 - t325 * t357 - t327 * t361) * t379;
t292 = (-t319 * t353 - t320 * t357 - t321 * t361) * t379;
t291 = t328 * t488 + t329 * t481 + t330 * t474;
t290 = t322 * t488 + t324 * t481 + t326 * t474;
t289 = (t323 * t453 + t325 * t444 + t327 * t435) * t379;
t288 = (t319 * t453 + t320 * t444 + t321 * t435) * t379;
t287 = t334 + t335 + t336 + (-t411 - t412 - t413) * t379;
t286 = (t394 + t396 + t398) * t379;
t285 = t299 + t297 + t295;
t284 = t307 + t309 + t311 + (-t337 * t431 - t338 * t430 - t339 * t429) * t379;
t283 = t322 * t328 * t484 + t324 * t329 * t477 + t326 * t330 * t470;
t282 = (t426 + t427 + t428) * t380;
t281 = (t343 * t428 + t346 * t427 + t349 * t426) * t380;
t280 = (t433 * t504 + t442 * t507 + t451 * t510) * t495;
t279 = (t353 * t407 + t357 * t406 + t361 * t405) * t380;
t278 = (t405 * t465 + t406 * t466 + t407 * t467) * t380;
t277 = (-t408 * t472 - t409 * t479 - t410 * t486) * t379;
t276 = (t408 * t432 + t409 * t441 + t410 * t450) * t379;
t1 = [t322 ^ 2 * t484 + t324 ^ 2 * t477 + t326 ^ 2 * t470, (t313 * t483 + t314 * t476 + t315 * t469) * t380, (t432 * t505 + t441 * t508 + t450 * t511) * t496, (-t472 * t505 - t479 * t508 - t486 * t511) * t496, (t313 * t458 + t314 * t449 + t315 * t440) * t380, (t313 * t451 + t314 * t442 + t315 * t433) * t495, (t434 * t503 + t443 * t506 + t452 * t509) * t495, (t473 * t503 + t480 * t506 + t487 * t509) * t495, (t337 ^ 2 * t353 + t338 ^ 2 * t357 + t339 ^ 2 * t361) * t380, t322 * t457 + t324 * t448 + t326 * t439 + (t322 * t404 + t324 * t403 + t326 * t402) * t379, t322 * t455 + t324 * t446 + t326 * t437 + (t397 * t322 + t395 * t324 + t393 * t326) * t379, 1; t283, t282, t276, t277, t281, t280, t278, t279, t306, t322 * t456 + t324 * t447 + t326 * t438 + (t328 * t404 + t329 * t403 + t330 * t402) * t379, t322 * t454 + t324 * t445 + t326 * t436 + (t328 * t397 + t329 * t395 + t330 * t393) * t379, 0; t290, 0, t288, t292, 0, 0, 0, 0, 0, t284, t286, 0; t283, t282, t276, t277, t281, t280, t278, t279, t306, t328 * t457 + t329 * t448 + t330 * t439 + (t322 * t401 + t324 * t400 + t326 * t399) * t379, t328 * t455 + t329 * t446 + t330 * t437 + (t322 * t392 + t324 * t391 + t326 * t390) * t379, 0; t328 ^ 2 * t484 + t329 ^ 2 * t477 + t330 ^ 2 * t470, (t316 * t483 + t317 * t476 + t318 * t469) * t380, (t432 * t498 + t441 * t500 + t450 * t502) * t496, (-t472 * t498 - t479 * t500 - t486 * t502) * t496, (t316 * t458 + t317 * t449 + t318 * t440) * t380, (t316 * t451 + t317 * t442 + t318 * t433) * t495, (-t434 * t497 - t443 * t499 - t452 * t501) * t495, (-t473 * t497 - t480 * t499 - t487 * t501) * t495, (t340 ^ 2 * t353 + t341 ^ 2 * t357 + t342 ^ 2 * t361) * t380, t328 * t456 + t329 * t447 + t330 * t438 + (t328 * t401 + t329 * t400 + t330 * t399) * t379, t328 * t454 + t329 * t445 + t330 * t436 + (t328 * t392 + t329 * t391 + t330 * t390) * t379, 1; t291, 0, t289, t293, 0, 0, 0, 0, 0, t285, t287, 0; t290, 0, t288, t292, 0, 0, 0, 0, 0, t284, t286, 0; t291, 0, t289, t293, 0, 0, 0, 0, 0, t285, t287, 0; 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
tau_reg  = t1;
