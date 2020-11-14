% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RRR1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [2x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3*(3+1)/2x10]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 15:38
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRR1G1A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(5,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRR1G1A0_inertia_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RRR1G1A0_inertia_para_pf_regmin: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRR1G1A0_inertia_para_pf_regmin: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRR1G1A0_inertia_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRR1G1A0_inertia_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 15:38:32
% EndTime: 2019-05-03 15:38:34
% DurationCPUTime: 1.33s
% Computational Cost: add. (3209->190), mult. (3643->388), div. (759->9), fcn. (4484->26), ass. (0->198)
t403 = qJ(1,3) + qJ(2,3);
t389 = sin(t403);
t392 = cos(t403);
t410 = sin(qJ(1,3));
t416 = cos(qJ(1,3));
t370 = t389 * t416 - t392 * t410;
t521 = 0.1e1 / t370;
t404 = qJ(1,2) + qJ(2,2);
t390 = sin(t404);
t393 = cos(t404);
t412 = sin(qJ(1,2));
t418 = cos(qJ(1,2));
t371 = t390 * t418 - t393 * t412;
t520 = 0.1e1 / t371;
t405 = qJ(1,1) + qJ(2,1);
t391 = sin(t405);
t394 = cos(t405);
t414 = sin(qJ(1,1));
t420 = cos(qJ(1,1));
t372 = t391 * t420 - t394 * t414;
t519 = 0.1e1 / t372;
t421 = xP(3);
t401 = sin(t421);
t402 = cos(t421);
t422 = koppelP(3,2);
t425 = koppelP(3,1);
t373 = t401 * t425 + t402 * t422;
t376 = -t401 * t422 + t402 * t425;
t406 = legFrame(3,3);
t395 = sin(t406);
t398 = cos(t406);
t322 = (t373 * t398 - t376 * t395) * t392 - (t373 * t395 + t376 * t398) * t389;
t379 = t410 * t425 - t416 * t422;
t380 = t410 * t422 + t416 * t425;
t304 = pkin(1) * ((t379 * t402 - t380 * t401) * t398 + t395 * (t379 * t401 + t380 * t402)) - t322 * pkin(2);
t518 = t304 * t521;
t423 = koppelP(2,2);
t426 = koppelP(2,1);
t374 = t401 * t426 + t402 * t423;
t377 = -t401 * t423 + t402 * t426;
t407 = legFrame(2,3);
t396 = sin(t407);
t399 = cos(t407);
t323 = (t374 * t399 - t377 * t396) * t393 - (t374 * t396 + t377 * t399) * t390;
t381 = t412 * t426 - t418 * t423;
t382 = t412 * t423 + t418 * t426;
t305 = pkin(1) * ((t381 * t402 - t382 * t401) * t399 + t396 * (t381 * t401 + t382 * t402)) - t323 * pkin(2);
t517 = t305 * t520;
t424 = koppelP(1,2);
t427 = koppelP(1,1);
t375 = t401 * t427 + t402 * t424;
t378 = -t401 * t424 + t402 * t427;
t408 = legFrame(1,3);
t397 = sin(t408);
t400 = cos(t408);
t324 = (t375 * t400 - t378 * t397) * t394 - (t375 * t397 + t378 * t400) * t391;
t383 = t414 * t427 - t420 * t424;
t384 = t414 * t424 + t420 * t427;
t306 = pkin(1) * ((t383 * t402 - t384 * t401) * t400 + t397 * (t383 * t401 + t384 * t402)) - t324 * pkin(2);
t516 = t306 * t519;
t356 = -t389 * t395 + t392 * t398;
t346 = pkin(1) * (-t395 * t410 + t398 * t416) + t356 * pkin(2);
t428 = 0.1e1 / pkin(2);
t429 = 0.1e1 / pkin(1);
t488 = t428 * t429;
t453 = t521 * t488;
t445 = t346 * t453;
t499 = t521 * t429;
t467 = t356 * t499;
t317 = -t445 + t467;
t515 = t317 * t521;
t358 = -t390 * t396 + t393 * t399;
t347 = pkin(1) * (-t396 * t412 + t399 * t418) + t358 * pkin(2);
t452 = t520 * t488;
t444 = t347 * t452;
t497 = t520 * t429;
t449 = t358 * t497;
t319 = -t444 + t449;
t514 = t319 * t520;
t360 = -t391 * t397 + t394 * t400;
t348 = pkin(1) * (-t397 * t414 + t400 * t420) + t360 * pkin(2);
t451 = t519 * t488;
t443 = t348 * t451;
t495 = t519 * t429;
t460 = t360 * t495;
t321 = -t443 + t460;
t513 = t321 * t519;
t512 = t322 * t521;
t511 = t323 * t520;
t510 = t324 * t519;
t355 = t389 * t398 + t392 * t395;
t343 = pkin(1) * (t395 * t416 + t398 * t410) + t355 * pkin(2);
t509 = t343 * t521;
t357 = t390 * t399 + t393 * t396;
t344 = pkin(1) * (t396 * t418 + t399 * t412) + t357 * pkin(2);
t508 = t344 * t520;
t359 = t391 * t400 + t394 * t397;
t345 = pkin(1) * (t397 * t420 + t400 * t414) + t359 * pkin(2);
t507 = t345 * t519;
t506 = t355 * t521;
t505 = t356 * t521;
t504 = t357 * t520;
t503 = t358 * t520;
t502 = t359 * t519;
t501 = t360 * t519;
t500 = t521 ^ 2;
t498 = t520 ^ 2;
t496 = t519 ^ 2;
t409 = sin(qJ(2,3));
t494 = t521 * t409;
t415 = cos(qJ(2,3));
t493 = t521 * t415;
t411 = sin(qJ(2,2));
t492 = t520 * t411;
t417 = cos(qJ(2,2));
t491 = t520 * t417;
t413 = sin(qJ(2,1));
t490 = t519 * t413;
t419 = cos(qJ(2,1));
t489 = t519 * t419;
t448 = t343 * t453;
t487 = t373 * t445 - t376 * t448;
t447 = t344 * t452;
t486 = t374 * t444 - t377 * t447;
t446 = t345 * t451;
t485 = t375 * t443 - t378 * t446;
t337 = t373 * t467;
t450 = t355 * t499;
t340 = t376 * t450;
t307 = -t337 + t340;
t338 = t374 * t449;
t464 = t357 * t497;
t341 = t377 * t464;
t308 = -t338 + t341;
t339 = t375 * t460;
t461 = t359 * t495;
t342 = t378 * t461;
t309 = -t339 + t342;
t484 = t307 * t518;
t483 = t308 * t517;
t482 = t309 * t516;
t313 = -t445 + 0.2e1 * t467;
t481 = t313 * t505;
t314 = -t444 + 0.2e1 * t449;
t480 = t314 * t503;
t315 = -t443 + 0.2e1 * t460;
t479 = t315 * t501;
t478 = t322 * t500;
t477 = t322 * t494;
t476 = t322 * t493;
t475 = t323 * t498;
t474 = t323 * t492;
t473 = t323 * t491;
t472 = t324 * t496;
t471 = t324 * t490;
t470 = t324 * t489;
t469 = t409 * t506;
t468 = t415 * t506;
t466 = t411 * t504;
t465 = t417 * t504;
t463 = t413 * t502;
t462 = t419 * t502;
t459 = t521 * t494;
t458 = t521 * t493;
t457 = t520 * t492;
t456 = t520 * t491;
t455 = t519 * t490;
t454 = t519 * t489;
t442 = t355 * t459;
t441 = t355 * t458;
t440 = t356 * t459;
t439 = t356 * t458;
t438 = t357 * t457;
t437 = t357 * t456;
t436 = t358 * t457;
t435 = t358 * t456;
t434 = t359 * t455;
t433 = t359 * t454;
t432 = t360 * t455;
t431 = t360 * t454;
t430 = 0.1e1 / pkin(1) ^ 2;
t385 = t401 ^ 2 + t402 ^ 2;
t366 = 0.1e1 / t372 ^ 2;
t364 = 0.1e1 / t371 ^ 2;
t362 = 0.1e1 / t370 ^ 2;
t320 = -t446 + t461;
t318 = -t447 + t464;
t316 = -t448 + t450;
t312 = -t446 + 0.2e1 * t461;
t311 = -t447 + 0.2e1 * t464;
t310 = -t448 + 0.2e1 * t450;
t303 = t309 + t485;
t302 = t308 + t486;
t301 = t307 + t487;
t300 = -0.2e1 * t339 + 0.2e1 * t342 + t485;
t299 = -0.2e1 * t338 + 0.2e1 * t341 + t486;
t298 = -0.2e1 * t337 + 0.2e1 * t340 + t487;
t1 = [(t356 ^ 2 * t500 + t358 ^ 2 * t498 + t360 ^ 2 * t496) * t430, 0, 0, (t317 * t505 + t319 * t503 + t321 * t501 + (-t346 * t515 - t347 * t514 - t348 * t513) * t428) * t429, t415 * t481 + t417 * t480 + t419 * t479 + (-t346 * t439 - t347 * t435 - t348 * t431) * t488, -t409 * t481 - t411 * t480 - t413 * t479 + (t346 * t440 + t347 * t436 + t348 * t432) * t488, 0, 0, 0, t385; (t355 * t356 * t362 + t357 * t358 * t364 + t359 * t360 * t366) * t430, 0, 0, (t317 * t506 + t319 * t504 + t321 * t502 + (-t317 * t509 - t319 * t508 - t321 * t507) * t428) * t429, t313 * t468 + t314 * t465 + t315 * t462 + (-t343 * t439 - t344 * t435 - t345 * t431) * t488, -t313 * t469 - t314 * t466 - t315 * t463 + (t343 * t440 + t344 * t436 + t345 * t432) * t488, 0, 0, 0, 0; (t355 ^ 2 * t362 + t357 ^ 2 * t364 + t359 ^ 2 * t366) * t430, 0, 0, (t316 * t506 + t318 * t504 + t320 * t502 + (-t316 * t509 - t318 * t508 - t320 * t507) * t428) * t429, t310 * t468 + t311 * t465 + t312 * t462 + (-t343 * t441 - t344 * t437 - t345 * t433) * t488, -t310 * t469 - t311 * t466 - t312 * t463 + (t343 * t442 + t344 * t438 + t345 * t434) * t488, 0, 0, 0, t385; (-t356 * t478 - t358 * t475 - t360 * t472) * t430, 0, 0, (-t317 * t512 - t319 * t511 - t321 * t510 + (-t304 * t515 - t305 * t514 - t306 * t513) * t428) * t429, -t313 * t476 - t314 * t473 - t315 * t470 + (-t304 * t439 - t305 * t435 - t306 * t431) * t488, t313 * t477 + t314 * t474 + t315 * t471 + (t304 * t440 + t305 * t436 + t306 * t432) * t488, 0, -t401, -t402, 0; (-t355 * t478 - t357 * t475 - t359 * t472) * t430, 0, 0, (-t316 * t512 - t318 * t511 - t320 * t510 + (-t316 * t518 - t318 * t517 - t320 * t516) * t428) * t429, -t310 * t476 - t311 * t473 - t312 * t470 + (-t304 * t441 - t305 * t437 - t306 * t433) * t488, t310 * t477 + t311 * t474 + t312 * t471 + (t304 * t442 + t305 * t438 + t306 * t434) * t488, 0, t402, -t401, 0; (-t307 * t512 - t308 * t511 - t309 * t510) * t429, 0, 0, (-t301 * t512 - t302 * t511 - t303 * t510 + (-t301 * t518 - t302 * t517 - t303 * t516) * t428) * t429, -t298 * t476 - t299 * t473 - t300 * t470 + (-t415 * t484 - t417 * t483 - t419 * t482) * t428, t298 * t477 + t299 * t474 + t300 * t471 + (t409 * t484 + t411 * t483 + t413 * t482) * t428, 1, 0, 0, 0;];
tau_reg  = t1;
