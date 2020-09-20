% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RPRR1G2A0
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-03-09 21:25
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RPRR1G2A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G2A0_inertia_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G2A0_inertia_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G2A0_inertia_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G2A0_inertia_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G2A0_inertia_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:25:03
% EndTime: 2020-03-09 21:25:05
% DurationCPUTime: 2.19s
% Computational Cost: add. (5700->262), mult. (3519->423), div. (669->7), fcn. (2937->62), ass. (0->229)
t484 = legFrame(1,2);
t552 = qJ(1,1) + pkin(7);
t462 = t484 + t552;
t453 = qJ(3,1) + t462;
t463 = -t484 + t552;
t454 = qJ(3,1) + t463;
t415 = sin(t453) + sin(t454);
t487 = sin(qJ(3,1));
t549 = pkin(7) + qJ(3,1);
t436 = pkin(1) * sin(t549) + t487 * pkin(2);
t432 = 0.1e1 / t436;
t559 = t432 / 0.2e1;
t616 = t415 * t559;
t418 = -cos(t454) + cos(t453);
t615 = t418 * t559;
t481 = cos(pkin(7));
t488 = cos(qJ(3,3));
t480 = sin(pkin(7));
t485 = sin(qJ(3,3));
t557 = t480 * t485;
t401 = -pkin(2) * t488 + (-t481 * t488 + t557) * pkin(1);
t614 = -t401 / 0.2e1;
t489 = cos(qJ(3,2));
t486 = sin(qJ(3,2));
t556 = t480 * t486;
t402 = -t489 * pkin(2) + (-t481 * t489 + t556) * pkin(1);
t613 = -t402 / 0.2e1;
t547 = pkin(7) + qJ(3,3);
t434 = pkin(1) * sin(t547) + t485 * pkin(2);
t429 = 0.1e1 / t434 ^ 2;
t612 = t429 / 0.4e1;
t548 = pkin(7) + qJ(3,2);
t435 = pkin(1) * sin(t548) + t486 * pkin(2);
t431 = 0.1e1 / t435 ^ 2;
t611 = t431 / 0.4e1;
t433 = 0.1e1 / t436 ^ 2;
t610 = t433 / 0.4e1;
t493 = 0.1e1 / pkin(3);
t554 = t493 / 0.2e1;
t466 = cos(qJ(1,1) + t549);
t412 = -pkin(3) * t466 - pkin(2) * cos(t552) - cos(qJ(1,1)) * pkin(1);
t579 = t412 * t432;
t400 = t493 * t579;
t421 = t466 * t432;
t385 = t421 + t400 / 0.2e1;
t467 = t481 * pkin(1) + pkin(2);
t490 = cos(qJ(3,1));
t591 = pkin(1) * t480;
t427 = t487 * t467 + t490 * t591;
t565 = t427 * t432;
t543 = -0.2e1 * t565;
t609 = t385 * t543 / 0.2e1;
t465 = cos(qJ(1,2) + t548);
t551 = qJ(1,2) + pkin(7);
t411 = -pkin(3) * t465 - pkin(2) * cos(t551) - cos(qJ(1,2)) * pkin(1);
t430 = 0.1e1 / t435;
t580 = t411 * t430;
t399 = t493 * t580;
t420 = t465 * t430;
t384 = t420 + t399 / 0.2e1;
t426 = t486 * t467 + t489 * t591;
t567 = t426 * t430;
t544 = -0.2e1 * t567;
t608 = t384 * t544 / 0.2e1;
t464 = cos(qJ(1,3) + t547);
t550 = qJ(1,3) + pkin(7);
t410 = -pkin(3) * t464 - pkin(2) * cos(t550) - cos(qJ(1,3)) * pkin(1);
t428 = 0.1e1 / t434;
t581 = t410 * t428;
t398 = t493 * t581;
t419 = t464 * t428;
t383 = t419 + t398 / 0.2e1;
t425 = t485 * t467 + t488 * t591;
t569 = t425 * t428;
t545 = -0.2e1 * t569;
t607 = t383 * t545 / 0.2e1;
t555 = t480 * t487;
t403 = t490 * pkin(2) + (t481 * t490 - t555) * pkin(1);
t606 = t403 * t615;
t605 = t403 * t616;
t483 = legFrame(2,2);
t460 = t483 + t551;
t451 = qJ(3,2) + t460;
t461 = -t483 + t551;
t452 = qJ(3,2) + t461;
t417 = -cos(t452) + cos(t451);
t574 = t417 * t430;
t604 = t574 * t613;
t414 = sin(t451) + sin(t452);
t577 = t414 * t430;
t603 = t577 * t613;
t482 = legFrame(3,2);
t458 = t482 + t550;
t449 = qJ(3,3) + t458;
t459 = -t482 + t550;
t450 = qJ(3,3) + t459;
t416 = -cos(t450) + cos(t449);
t575 = t416 * t428;
t602 = t575 * t614;
t413 = sin(t449) + sin(t450);
t578 = t413 * t428;
t601 = t578 * t614;
t558 = t433 * t466;
t600 = t558 / 0.2e1;
t560 = t431 * t465;
t599 = t560 / 0.2e1;
t562 = t429 * t464;
t598 = t562 / 0.2e1;
t596 = t574 / 0.2e1;
t595 = t575 / 0.2e1;
t593 = t577 / 0.2e1;
t592 = t578 / 0.2e1;
t468 = qJ(1,3) + t482;
t469 = qJ(1,3) - t482;
t377 = t413 * pkin(3) + (sin(t458) + sin(t459)) * pkin(2) + (sin(t468) + sin(t469)) * pkin(1);
t590 = t377 * t428;
t470 = qJ(1,2) + t483;
t471 = qJ(1,2) - t483;
t378 = t414 * pkin(3) + (sin(t460) + sin(t461)) * pkin(2) + (sin(t470) + sin(t471)) * pkin(1);
t589 = t378 * t430;
t472 = qJ(1,1) + t484;
t473 = qJ(1,1) - t484;
t379 = t415 * pkin(3) + (sin(t462) + sin(t463)) * pkin(2) + (sin(t472) + sin(t473)) * pkin(1);
t588 = t379 * t432;
t380 = -t416 * pkin(3) + (cos(t459) - cos(t458)) * pkin(2) + (cos(t469) - cos(t468)) * pkin(1);
t587 = t380 * t428;
t381 = -t417 * pkin(3) + (cos(t461) - cos(t460)) * pkin(2) + (cos(t471) - cos(t470)) * pkin(1);
t586 = t381 * t430;
t382 = -t418 * pkin(3) + (cos(t463) - cos(t462)) * pkin(2) + (cos(t473) - cos(t472)) * pkin(1);
t585 = t382 * t432;
t387 = t419 + t398;
t584 = t387 * t428;
t389 = t420 + t399;
t583 = t389 * t430;
t391 = t421 + t400;
t582 = t391 * t432;
t422 = -pkin(1) * t557 + t488 * t467;
t572 = t422 * t429;
t423 = -pkin(1) * t556 + t489 * t467;
t571 = t423 * t431;
t424 = -pkin(1) * t555 + t490 * t467;
t570 = t424 * t433;
t568 = t425 * t429;
t566 = t426 * t431;
t564 = t427 * t433;
t563 = t428 / 0.2e1;
t561 = t430 / 0.2e1;
t553 = t493 / 0.4e1;
t540 = t401 * t419;
t537 = t402 * t420;
t534 = t403 * t421;
t533 = t413 * t572;
t532 = t413 * t569;
t531 = t413 * t568;
t530 = t414 * t571;
t529 = t414 * t567;
t528 = t414 * t566;
t527 = t415 * t570;
t526 = t415 * t565;
t525 = t415 * t564;
t524 = t416 * t572;
t523 = t416 * t569;
t522 = t416 * t568;
t521 = t417 * t571;
t520 = t417 * t567;
t519 = t417 * t566;
t518 = t418 * t570;
t517 = t418 * t565;
t516 = t418 * t564;
t515 = t422 * t562;
t514 = t423 * t560;
t513 = t424 * t558;
t512 = t425 * t562;
t511 = t426 * t560;
t510 = t427 * t558;
t509 = t428 * t554;
t508 = t430 * t554;
t507 = t432 * t554;
t503 = t464 * t545;
t502 = t465 * t544;
t501 = t466 * t543;
t500 = t377 * t509;
t499 = t378 * t508;
t498 = t379 * t507;
t497 = t429 * t464 ^ 2 + t431 * t465 ^ 2 + t433 * t466 ^ 2;
t496 = t413 ^ 2 * t612 + t414 ^ 2 * t611 + t415 ^ 2 * t610;
t495 = t416 ^ 2 * t612 + t417 ^ 2 * t611 + t418 ^ 2 * t610;
t351 = t413 * t416 * t612 + t414 * t417 * t611 + t415 * t418 * t610;
t354 = t413 * t598 + t414 * t599 + t415 * t600;
t355 = t416 * t598 + t417 * t599 + t418 * t600;
t494 = pkin(1) ^ 2;
t479 = cos(t484);
t478 = cos(t483);
t477 = cos(t482);
t476 = sin(t484);
t475 = sin(t483);
t474 = sin(t482);
t396 = t417 * t561;
t395 = t416 * t563;
t393 = t414 * t561;
t392 = t413 * t563;
t390 = 0.2e1 * t421 + t400;
t388 = 0.2e1 * t420 + t399;
t386 = 0.2e1 * t419 + t398;
t376 = t382 * t507;
t375 = t381 * t508;
t374 = t380 * t509;
t373 = t615 + t376;
t372 = 0.2e1 * t615 + t376;
t371 = t396 + t375;
t370 = 0.2e1 * t396 + t375;
t369 = t395 + t374;
t368 = 0.2e1 * t395 + t374;
t367 = t616 - t498;
t366 = 0.2e1 * t616 - t498;
t365 = t393 - t499;
t364 = 0.2e1 * t393 - t499;
t363 = t392 - t500;
t362 = 0.2e1 * t392 - t500;
t361 = t615 + t376 / 0.2e1;
t360 = t396 + t375 / 0.2e1;
t359 = t395 + t374 / 0.2e1;
t358 = t616 - t498 / 0.2e1;
t357 = t393 - t499 / 0.2e1;
t356 = t392 - t500 / 0.2e1;
t353 = t494 * t355;
t352 = t494 * t354;
t350 = t494 * t351 + t474 * t477 + t475 * t478 + t476 * t479;
t1 = [t496, 0, 0, t474 ^ 2 + t475 ^ 2 + t476 ^ 2 + t494 * t496, t363 * t592 + t365 * t593 + t367 * t616 + (-t363 * t590 - t365 * t589 - t367 * t588) * t554, (-t377 * t533 - t378 * t530 - t379 * t527) * t553 + t362 * t601 + t364 * t603 + t366 * t605, (t377 * t531 + t378 * t528 + t379 * t525) * t553 - t356 * t532 - t357 * t529 - t358 * t526, 1; t351, 0, 0, t350, t369 * t592 + t371 * t593 + t373 * t616 + (-t369 * t590 - t371 * t589 - t373 * t588) * t554, (-t377 * t524 - t378 * t521 - t379 * t518) * t553 + t368 * t601 + t370 * t603 + t372 * t605, (t377 * t522 + t378 * t519 + t379 * t516) * t553 - t359 * t532 - t360 * t529 - t361 * t526, 0; t354, 0, 0, t352, t387 * t592 + t389 * t593 + t391 * t616 + (-t377 * t584 - t378 * t583 - t379 * t582) * t554, t386 * t601 + t388 * t603 + t390 * t605 + (-t377 * t515 - t378 * t514 - t379 * t513) * t554, t413 * t607 + t414 * t608 + t415 * t609 + (t377 * t512 + t378 * t511 + t379 * t510) * t554, 0; t351, 0, 0, t350, t363 * t595 + t365 * t596 + t367 * t615 + (t363 * t587 + t365 * t586 + t367 * t585) * t554, (t380 * t533 + t381 * t530 + t382 * t527) * t553 + t362 * t602 + t364 * t604 + t366 * t606, (-t380 * t531 - t381 * t528 - t382 * t525) * t553 - t356 * t523 - t357 * t520 - t358 * t517, 0; t495, 0, 0, t477 ^ 2 + t478 ^ 2 + t479 ^ 2 + t494 * t495, t369 * t595 + t371 * t596 + t373 * t615 + (t369 * t587 + t371 * t586 + t373 * t585) * t554, (t380 * t524 + t381 * t521 + t382 * t518) * t553 + t368 * t602 + t370 * t604 + t372 * t606, (-t380 * t522 - t381 * t519 - t382 * t516) * t553 - t359 * t523 - t360 * t520 - t361 * t517, 1; t355, 0, 0, t353, t387 * t595 + t389 * t596 + t391 * t615 + (t380 * t584 + t381 * t583 + t382 * t582) * t554, t386 * t602 + t388 * t604 + t390 * t606 + (t380 * t515 + t381 * t514 + t382 * t513) * t554, t416 * t607 + t417 * t608 + t418 * t609 + (-t380 * t512 - t381 * t511 - t382 * t510) * t554, 0; t354, 0, 0, t352, t363 * t419 + t365 * t420 + t367 * t421 + (t363 * t581 + t365 * t580 + t367 * t579) * t493, -t362 * t540 - t364 * t537 + t366 * t534 + (t410 * t533 + t411 * t530 + t412 * t527) * t554, t356 * t503 + t357 * t502 + t358 * t501 + (-t410 * t531 - t411 * t528 - t412 * t525) * t554, 0; t355, 0, 0, t353, t369 * t419 + t371 * t420 + t373 * t421 + (t369 * t581 + t371 * t580 + t373 * t579) * t493, -t368 * t540 - t370 * t537 + t372 * t534 + (t410 * t524 + t411 * t521 + t412 * t518) * t554, t359 * t503 + t360 * t502 + t361 * t501 + (-t410 * t522 - t411 * t519 - t412 * t516) * t554, 0; t497, 0, 0, t497 * t494, t387 * t419 + t389 * t420 + t391 * t421 + (t387 * t581 + t389 * t580 + t391 * t579) * t493, -t386 * t540 - t388 * t537 + t390 * t534 + (t410 * t515 + t411 * t514 + t412 * t513) * t493, t383 * t503 + t384 * t502 + t385 * t501 + (-t410 * t512 - t411 * t511 - t412 * t510) * t493, 1;];
tau_reg  = t1;
