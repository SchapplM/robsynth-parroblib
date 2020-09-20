% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3PRRR1G2P2A0
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
% Datum: 2020-03-09 21:18
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRR1G2P2A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G2P2A0_inertia_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G2P2A0_inertia_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G2P2A0_inertia_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G2P2A0_inertia_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G2P2A0_inertia_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:18:27
% EndTime: 2020-03-09 21:18:28
% DurationCPUTime: 1.56s
% Computational Cost: add. (4676->181), mult. (2019->384), div. (981->9), fcn. (2856->24), ass. (0->203)
t375 = legFrame(1,2);
t366 = sin(t375);
t369 = cos(t375);
t470 = t366 * t369;
t372 = pkin(7) + qJ(2,1);
t363 = qJ(3,1) + t372;
t348 = sin(t363);
t354 = sin(t372);
t336 = pkin(2) * t354 + pkin(3) * t348;
t479 = t336 * t348;
t351 = cos(t363);
t503 = cos(t372);
t330 = -t503 * t348 + t351 * t354;
t506 = 0.1e1 / t330 ^ 2;
t519 = t470 * t479 * t506;
t374 = legFrame(2,2);
t365 = sin(t374);
t368 = cos(t374);
t471 = t365 * t368;
t371 = pkin(7) + qJ(2,2);
t362 = qJ(3,2) + t371;
t347 = sin(t362);
t353 = sin(t371);
t335 = pkin(2) * t353 + pkin(3) * t347;
t483 = t335 * t347;
t350 = cos(t362);
t504 = cos(t371);
t329 = -t504 * t347 + t350 * t353;
t507 = 0.1e1 / t329 ^ 2;
t518 = t471 * t483 * t507;
t373 = legFrame(3,2);
t364 = sin(t373);
t367 = cos(t373);
t472 = t364 * t367;
t370 = pkin(7) + qJ(2,3);
t361 = qJ(3,3) + t370;
t346 = sin(t361);
t352 = sin(t370);
t334 = pkin(2) * t352 + pkin(3) * t346;
t487 = t334 * t346;
t349 = cos(t361);
t505 = cos(t370);
t328 = -t505 * t346 + t349 * t352;
t508 = 0.1e1 / t328 ^ 2;
t517 = t472 * t487 * t508;
t343 = t346 ^ 2;
t516 = t343 * t508;
t344 = t347 ^ 2;
t515 = t344 * t507;
t345 = t348 ^ 2;
t514 = t345 * t506;
t379 = cos(qJ(3,3));
t380 = cos(qJ(3,2));
t381 = cos(qJ(3,1));
t382 = 0.1e1 / pkin(3);
t383 = 0.1e1 / pkin(2);
t469 = t382 * t383;
t513 = (t379 * t517 + t380 * t518 + t381 * t519) * t469;
t376 = sin(qJ(3,3));
t377 = sin(qJ(3,2));
t378 = sin(qJ(3,1));
t512 = (-t376 * t517 - t377 * t518 - t378 * t519) * t469;
t511 = 0.1e1 / t328;
t510 = 0.1e1 / t329;
t509 = 0.1e1 / t330;
t502 = t511 ^ 2;
t501 = t511 * t364;
t500 = t510 ^ 2;
t499 = t510 * t365;
t498 = t509 ^ 2;
t497 = t509 * t366;
t337 = pkin(2) * t505 + pkin(3) * t349;
t496 = t511 * t337;
t495 = t511 * t349;
t494 = t508 * t349;
t338 = pkin(2) * t504 + pkin(3) * t350;
t493 = t510 * t338;
t492 = t510 * t350;
t491 = t507 * t350;
t339 = pkin(2) * t503 + pkin(3) * t351;
t490 = t509 * t339;
t489 = t509 * t351;
t488 = t506 * t351;
t486 = t334 * t364;
t485 = t334 * t367;
t484 = t334 * t382;
t482 = t335 * t365;
t481 = t335 * t368;
t480 = t335 * t382;
t478 = t336 * t366;
t477 = t336 * t369;
t476 = t336 * t382;
t475 = t346 * t367;
t474 = t347 * t368;
t473 = t348 * t369;
t468 = t511 * t485;
t467 = t346 * t501;
t466 = t383 * t501;
t465 = t511 * t469;
t464 = t510 * t481;
t463 = t347 * t499;
t462 = t383 * t499;
t461 = t510 * t469;
t460 = t509 * t477;
t459 = t348 * t497;
t458 = t383 * t497;
t457 = t509 * t469;
t456 = t511 * t486;
t455 = t511 * t475;
t454 = t376 * t495;
t453 = t379 * t495;
t452 = t337 * t494;
t451 = t508 * t475;
t450 = t510 * t482;
t449 = t510 * t474;
t448 = t377 * t492;
t447 = t380 * t492;
t446 = t338 * t491;
t445 = t507 * t474;
t444 = t509 * t478;
t443 = t509 * t473;
t442 = t378 * t489;
t441 = t381 * t489;
t440 = t339 * t488;
t439 = t506 * t473;
t438 = t383 * t495;
t437 = t383 * t492;
t436 = t383 * t489;
t435 = t487 * t502;
t434 = t511 * t467;
t433 = t376 * t467;
t432 = t379 * t467;
t431 = t483 * t500;
t430 = t510 * t463;
t429 = t377 * t463;
t428 = t380 * t463;
t427 = t479 * t498;
t426 = t509 * t459;
t425 = t378 * t459;
t424 = t381 * t459;
t423 = t376 * t455;
t422 = t379 * t455;
t421 = t486 * t494;
t420 = t337 * t451;
t419 = t377 * t449;
t418 = t380 * t449;
t417 = t482 * t491;
t416 = t338 * t445;
t415 = t378 * t443;
t414 = t381 * t443;
t413 = t478 * t488;
t412 = t339 * t439;
t408 = t383 * t455;
t407 = t383 * t449;
t406 = t383 * t443;
t405 = t376 * t435;
t404 = t379 * t435;
t403 = t468 * t495;
t402 = t337 * t434;
t401 = t377 * t431;
t400 = t380 * t431;
t399 = t464 * t492;
t398 = t338 * t430;
t397 = t378 * t427;
t396 = t381 * t427;
t395 = t460 * t489;
t394 = t339 * t426;
t384 = 0.1e1 / pkin(2) ^ 2;
t360 = t369 ^ 2;
t359 = t368 ^ 2;
t358 = t367 ^ 2;
t357 = t366 ^ 2;
t356 = t365 ^ 2;
t355 = t364 ^ 2;
t315 = t470 + t471 + t472;
t314 = t339 * t457;
t313 = t338 * t461;
t312 = t337 * t465;
t311 = t457 * t477;
t310 = t461 * t481;
t309 = t465 * t485;
t308 = t314 - t436;
t307 = t314 - 0.2e1 * t436;
t306 = t313 - t437;
t305 = t313 - 0.2e1 * t437;
t304 = t312 - t438;
t303 = t312 - 0.2e1 * t438;
t302 = (t348 - t476) * t458;
t301 = (0.2e1 * t348 - t476) * t458;
t300 = (t347 - t480) * t462;
t299 = (0.2e1 * t347 - t480) * t462;
t298 = (t346 - t484) * t466;
t297 = (0.2e1 * t346 - t484) * t466;
t296 = t311 - t406;
t295 = t311 - 0.2e1 * t406;
t294 = t310 - t407;
t293 = t310 - 0.2e1 * t407;
t292 = t309 - t408;
t291 = t309 - 0.2e1 * t408;
t290 = (t349 * t451 + t350 * t445 + t351 * t439) * t384;
t289 = (-t343 * t472 * t502 - t344 * t471 * t500 - t345 * t470 * t498) * t384;
t288 = (-t349 * t434 - t350 * t430 - t351 * t426) * t384;
t1 = [t357 + t356 + t355, (t358 * t516 + t359 * t515 + t360 * t514) * t384, 0, 0, (-t292 * t455 - t294 * t449 - t296 * t443 + (t292 * t468 + t294 * t464 + t296 * t460) * t382) * t383, -t291 * t422 - t293 * t418 - t295 * t414 + (-t358 * t404 - t359 * t400 - t360 * t396) * t469, t291 * t423 + t293 * t419 + t295 * t415 + (t358 * t405 + t359 * t401 + t360 * t397) * t469, 1; t315, t289, 0, 0, (-t298 * t455 - t300 * t449 - t302 * t443 + (t298 * t468 + t300 * t464 + t302 * t460) * t382) * t383, -t297 * t422 - t299 * t418 - t301 * t414 + t513, t297 * t423 + t299 * t419 + t301 * t415 + t512, 0; 0, t290, 0, 0, (-t304 * t455 - t306 * t449 - t308 * t443 + (t304 * t468 + t306 * t464 + t308 * t460) * t382) * t383, -t303 * t422 - t305 * t418 - t307 * t414 + (-t379 * t403 - t380 * t399 - t381 * t395) * t469, t303 * t423 + t305 * t419 + t307 * t415 + (t376 * t403 + t377 * t399 + t378 * t395) * t469, 0; t315, t289, 0, 0, (t292 * t467 + t294 * t463 + t296 * t459 + (-t292 * t456 - t294 * t450 - t296 * t444) * t382) * t383, t291 * t432 + t293 * t428 + t295 * t424 + t513, -t291 * t433 - t293 * t429 - t295 * t425 + t512, 0; t360 + t359 + t358, (t355 * t516 + t356 * t515 + t357 * t514) * t384, 0, 0, (t298 * t467 + t300 * t463 + t302 * t459 + (-t298 * t456 - t300 * t450 - t302 * t444) * t382) * t383, t297 * t432 + t299 * t428 + t301 * t424 + (-t355 * t404 - t356 * t400 - t357 * t396) * t469, -t297 * t433 - t299 * t429 - t301 * t425 + (t355 * t405 + t356 * t401 + t357 * t397) * t469, 1; 0, t288, 0, 0, (t304 * t467 + t306 * t463 + t308 * t459 + (-t304 * t456 - t306 * t450 - t308 * t444) * t382) * t383, t303 * t432 + t305 * t428 + t307 * t424 + (t379 * t421 + t380 * t417 + t381 * t413) * t469, -t303 * t433 - t305 * t429 - t307 * t425 + (-t376 * t421 - t377 * t417 - t378 * t413) * t469, 0; 0, t290, 0, 0, (-t292 * t495 - t294 * t492 - t296 * t489 + (t292 * t496 + t294 * t493 + t296 * t490) * t382) * t383, -t291 * t453 - t293 * t447 - t295 * t441 + (-t379 * t420 - t380 * t416 - t381 * t412) * t469, t291 * t454 + t293 * t448 + t295 * t442 + (t376 * t420 + t377 * t416 + t378 * t412) * t469, 0; 0, t288, 0, 0, (-t298 * t495 - t300 * t492 - t302 * t489 + (t298 * t496 + t300 * t493 + t302 * t490) * t382) * t383, -t297 * t453 - t299 * t447 - t301 * t441 + (t379 * t402 + t380 * t398 + t381 * t394) * t469, t297 * t454 + t299 * t448 + t301 * t442 + (-t376 * t402 - t377 * t398 - t378 * t394) * t469, 0; 0, (t349 ^ 2 * t508 + t350 ^ 2 * t507 + t351 ^ 2 * t506) * t384, 0, 0, (-t304 * t495 - t306 * t492 - t308 * t489 + (t304 * t496 + t306 * t493 + t308 * t490) * t382) * t383, -t303 * t453 - t305 * t447 - t307 * t441 + (-t379 * t452 - t380 * t446 - t381 * t440) * t469, t303 * t454 + t305 * t448 + t307 * t442 + (t376 * t452 + t377 * t446 + t378 * t440) * t469, 1;];
tau_reg  = t1;
