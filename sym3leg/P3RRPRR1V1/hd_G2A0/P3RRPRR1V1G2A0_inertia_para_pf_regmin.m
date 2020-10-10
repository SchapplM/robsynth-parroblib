% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RRPRR1V1G2A0
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3*3x13]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:34
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRPRR1V1G2A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR1V1G2A0_inertia_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR1V1G2A0_inertia_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RRPRR1V1G2A0_inertia_para_pf_regmin: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR1V1G2A0_inertia_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR1V1G2A0_inertia_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:34:03
% EndTime: 2020-08-06 19:34:04
% DurationCPUTime: 1.88s
% Computational Cost: add. (1510->217), mult. (2250->545), div. (855->14), fcn. (2499->18), ass. (0->240)
t435 = legFrame(3,2);
t406 = sin(t435);
t409 = cos(t435);
t438 = sin(qJ(2,3));
t439 = sin(qJ(1,3));
t444 = cos(qJ(2,3));
t537 = t439 * t444;
t397 = -t406 * t537 + t438 * t409;
t400 = t438 * t406 + t409 * t537;
t580 = t397 * t400;
t436 = legFrame(2,2);
t407 = sin(t436);
t410 = cos(t436);
t440 = sin(qJ(2,2));
t441 = sin(qJ(1,2));
t446 = cos(qJ(2,2));
t534 = t441 * t446;
t398 = -t407 * t534 + t440 * t410;
t401 = t440 * t407 + t410 * t534;
t579 = t398 * t401;
t437 = legFrame(1,2);
t408 = sin(t437);
t411 = cos(t437);
t442 = sin(qJ(2,1));
t443 = sin(qJ(1,1));
t448 = cos(qJ(2,1));
t531 = t443 * t448;
t399 = -t408 * t531 + t442 * t411;
t402 = t442 * t408 + t411 * t531;
t578 = t399 * t402;
t577 = 2 * qJ(3,1);
t576 = 2 * qJ(3,2);
t575 = 2 * qJ(3,3);
t450 = pkin(1) + pkin(2);
t430 = 1 / t450;
t574 = 2 * t430;
t573 = pkin(1) * t430;
t572 = pkin(1) * t444;
t571 = pkin(1) * t446;
t570 = pkin(1) * t448;
t421 = 0.1e1 / t444;
t569 = t397 * t421;
t424 = 0.1e1 / t446;
t568 = t398 * t424;
t427 = 0.1e1 / t448;
t567 = t399 * t427;
t566 = t400 * t421;
t565 = t401 * t424;
t564 = t402 * t427;
t563 = t406 * t421;
t562 = t407 * t424;
t561 = t408 * t427;
t560 = t409 * t421;
t559 = t410 * t424;
t558 = t411 * t427;
t432 = pkin(3) + qJ(3,3);
t412 = 1 / t432;
t557 = t412 * t421;
t556 = t412 * t438;
t445 = cos(qJ(1,3));
t555 = t412 * t445;
t413 = 1 / t432 ^ 2;
t554 = t413 * t421;
t455 = t444 ^ 2;
t422 = 0.1e1 / t455;
t553 = t413 * t422;
t433 = pkin(3) + qJ(3,2);
t414 = 1 / t433;
t552 = t414 * t424;
t551 = t414 * t440;
t447 = cos(qJ(1,2));
t550 = t414 * t447;
t415 = 1 / t433 ^ 2;
t549 = t415 * t424;
t456 = t446 ^ 2;
t425 = 0.1e1 / t456;
t548 = t415 * t425;
t434 = pkin(3) + qJ(3,1);
t416 = 1 / t434;
t547 = t416 * t427;
t546 = t416 * t442;
t449 = cos(qJ(1,1));
t545 = t416 * t449;
t417 = 1 / t434 ^ 2;
t544 = t417 * t427;
t457 = t448 ^ 2;
t428 = 0.1e1 / t457;
t543 = t417 * t428;
t542 = t445 ^ 2 * t413;
t541 = t447 ^ 2 * t415;
t540 = t449 ^ 2 * t417;
t539 = t438 * t445;
t538 = t438 * t450;
t536 = t440 * t447;
t535 = t440 * t450;
t533 = t442 * t449;
t532 = t442 * t450;
t530 = t444 * t450;
t529 = t446 * t450;
t528 = t448 * t450;
t527 = t416 * t577;
t526 = t414 * t576;
t525 = t412 * t575;
t524 = t438 * t573;
t523 = t440 * t573;
t522 = t442 * t573;
t521 = qJ(3,1) * t546;
t520 = qJ(3,2) * t551;
t519 = qJ(3,3) * t556;
t518 = t397 * t557;
t517 = t398 * t552;
t516 = t399 * t547;
t515 = t400 * t557;
t514 = t401 * t552;
t513 = t402 * t547;
t512 = t422 * t556;
t418 = t438 ^ 2;
t511 = t418 * t553;
t510 = t438 * t554;
t509 = t445 * t554;
t508 = t413 * t539;
t507 = t425 * t551;
t419 = t440 ^ 2;
t506 = t419 * t548;
t505 = t440 * t549;
t504 = t447 * t549;
t503 = t415 * t536;
t502 = t428 * t546;
t420 = t442 ^ 2;
t501 = t420 * t543;
t500 = t442 * t544;
t499 = t449 * t544;
t498 = t417 * t533;
t497 = t553 * t580;
t496 = t548 * t579;
t495 = t543 * t578;
t494 = t406 * t512;
t493 = t407 * t507;
t492 = t408 * t502;
t491 = t409 * t512;
t490 = t410 * t507;
t489 = t411 * t502;
t488 = t539 * t557;
t487 = t418 * t509;
t486 = t536 * t552;
t485 = t419 * t504;
t484 = t533 * t547;
t483 = t420 * t499;
t482 = t499 * t577;
t481 = t504 * t576;
t480 = t509 * t575;
t479 = qJ(3,1) * t427 * t522;
t478 = qJ(3,2) * t424 * t523;
t477 = qJ(3,3) * t421 * t524;
t476 = t406 * t488;
t475 = t407 * t486;
t474 = t408 * t484;
t473 = t409 * t488;
t472 = t410 * t486;
t471 = t411 * t484;
t451 = qJ(3,3) ^ 2;
t454 = pkin(1) ^ 2;
t470 = t421 * t451 + t444 * t454;
t452 = qJ(3,2) ^ 2;
t469 = t424 * t452 + t446 * t454;
t453 = qJ(3,1) ^ 2;
t468 = t427 * t453 + t448 * t454;
t467 = -t432 * t445 + t439 * t530;
t466 = -t433 * t447 + t441 * t529;
t465 = -t434 * t449 + t443 * t528;
t464 = t412 * (t397 * t406 + t400 * t409);
t463 = t414 * (t398 * t407 + t401 * t410);
t462 = t416 * (t399 * t408 + t402 * t411);
t461 = t397 * t491 + t398 * t490 + t399 * t489;
t460 = t400 * t494 + t401 * t493 + t402 * t492;
t459 = t474 + t475 + t476;
t458 = t471 + t472 + t473;
t431 = 0.1e1 / t450 ^ 2;
t405 = t443 * t434 + t449 * t528;
t404 = t441 * t433 + t447 * t529;
t403 = t439 * t432 + t445 * t530;
t396 = t402 ^ 2;
t395 = t401 ^ 2;
t394 = t400 ^ 2;
t393 = t399 ^ 2;
t392 = t398 ^ 2;
t391 = t397 ^ 2;
t390 = (-t449 * t570 + t405) * t416;
t389 = (-t447 * t571 + t404) * t414;
t388 = (-t445 * t572 + t403) * t412;
t387 = t408 * t532 + t465 * t411;
t386 = t407 * t535 + t466 * t410;
t385 = t406 * t538 + t467 * t409;
t384 = -t466 * t407 + t410 * t535;
t383 = -t465 * t408 + t411 * t532;
t382 = -t467 * t406 + t409 * t538;
t381 = (t409 * t555 + t410 * t550 + t411 * t545) * t430;
t380 = (t406 * t555 + t407 * t550 + t408 * t545) * t430;
t379 = (t402 * t527 - t408 * t522) * t427;
t378 = (t401 * t526 - t407 * t523) * t424;
t377 = (t400 * t525 - t406 * t524) * t421;
t376 = (t399 * t527 - t411 * t522) * t427;
t375 = (t398 * t526 - t410 * t523) * t424;
t374 = (t397 * t525 - t409 * t524) * t421;
t373 = (-t402 * t521 + t408 * t573) * t427;
t372 = (-t399 * t521 + t411 * t573) * t427;
t371 = (-t401 * t520 + t407 * t573) * t424;
t370 = (-t398 * t520 + t410 * t573) * t424;
t369 = (-t400 * t519 + t406 * t573) * t421;
t368 = (-t397 * t519 + t409 * t573) * t421;
t367 = (-t405 * t570 + (t454 * t457 + t453) * t449) * t416;
t366 = (-t404 * t571 + (t454 * t456 + t452) * t447) * t414;
t365 = (-t403 * t572 + (t454 * t455 + t451) * t445) * t412;
t364 = (t406 * t409 * t422 + t407 * t410 * t425 + t408 * t411 * t428) * t431;
t363 = t458 * t430;
t362 = t459 * t430;
t361 = (-pkin(1) * t402 + t387) * t416;
t360 = (-pkin(1) * t401 + t386) * t414;
t359 = (-pkin(1) * t400 + t385) * t412;
t358 = (-pkin(1) * t399 + t383) * t416;
t357 = (-pkin(1) * t398 + t384) * t414;
t356 = (-pkin(1) * t397 + t382) * t412;
t355 = t400 * t509 + t401 * t504 + t402 * t499;
t354 = t397 * t509 + t398 * t504 + t399 * t499;
t353 = 0.2e1 * t400 * t508 + 0.2e1 * t401 * t503 + 0.2e1 * t402 * t498;
t352 = 0.2e1 * t397 * t508 + 0.2e1 * t398 * t503 + 0.2e1 * t399 * t498;
t351 = t400 * t487 + t401 * t485 + t402 * t483;
t350 = t397 * t487 + t398 * t485 + t399 * t483;
t349 = -t408 * t479 + (-t387 * t570 + t468 * t402) * t416;
t348 = -t411 * t479 + (-t383 * t570 + t468 * t399) * t416;
t347 = -t407 * t478 + (-t386 * t571 + t469 * t401) * t414;
t346 = -t410 * t478 + (-t384 * t571 + t469 * t398) * t414;
t345 = -t406 * t477 + (-t385 * t572 + t470 * t400) * t412;
t344 = -t409 * t477 + (-t382 * t572 + t470 * t397) * t412;
t343 = t495 + t496 + t497;
t342 = t418 * t497 + t419 * t496 + t420 * t495;
t341 = 0.2e1 * t500 * t578 + 0.2e1 * t505 * t579 + 0.2e1 * t510 * t580;
t340 = (t421 * t464 + t424 * t463 + t427 * t462) * t430;
t339 = (t438 * t422 * t464 + t440 * t425 * t463 + t442 * t428 * t462) * t430;
t1 = [t394 * t553 + t395 * t548 + t396 * t543, 0, 0, t394 * t511 + t395 * t506 + t396 * t501, 0.2e1 * t394 * t510 + 0.2e1 * t395 * t505 + 0.2e1 * t396 * t500, t460 * t574, (t406 * t515 + t407 * t514 + t408 * t513) * t574, (t406 ^ 2 * t422 + t407 ^ 2 * t425 + t408 ^ 2 * t428) * t431, 0, 0, t377 * t515 + t378 * t514 + t379 * t513 - t460 * t573, (t349 * t564 + t361 * t387) * t416 + (t347 * t565 + t360 * t386) * t414 + (t345 * t566 + t359 * t385) * t412 + (t369 * t563 + t371 * t562 + t373 * t561) * t573, 1; t343, 0, 0, t342, t341, t339, t340, t364, 0, 0, t374 * t515 + t375 * t514 + t376 * t513 + (-t397 * t494 - t398 * t493 - t399 * t492) * t573, (t348 * t564 + t358 * t387) * t416 + (t346 * t565 + t357 * t386) * t414 + (t344 * t566 + t356 * t385) * t412 + (t368 * t563 + t370 * t562 + t372 * t561) * t573, 0; t355, 0, 0, t351, t353, t362, t380, 0, 0, 0, t400 * t480 + t401 * t481 + t402 * t482 - t459 * t573, (t367 * t564 + t387 * t390) * t416 + (t366 * t565 + t386 * t389) * t414 + (t365 * t566 + t385 * t388) * t412 + (-qJ(3,1) * t474 - qJ(3,2) * t475 - qJ(3,3) * t476) * t573, 0; t343, 0, 0, t342, t341, t339, t340, t364, 0, 0, t377 * t518 + t378 * t517 + t379 * t516 + (-t400 * t491 - t401 * t490 - t402 * t489) * t573, (t349 * t567 + t361 * t383) * t416 + (t347 * t568 + t360 * t384) * t414 + (t345 * t569 + t359 * t382) * t412 + (t369 * t560 + t371 * t559 + t373 * t558) * t573, 0; t391 * t553 + t392 * t548 + t393 * t543, 0, 0, t391 * t511 + t392 * t506 + t393 * t501, 0.2e1 * t391 * t510 + 0.2e1 * t392 * t505 + 0.2e1 * t393 * t500, t461 * t574, (t409 * t518 + t410 * t517 + t411 * t516) * t574, (t409 ^ 2 * t422 + t410 ^ 2 * t425 + t411 ^ 2 * t428) * t431, 0, 0, t374 * t518 + t375 * t517 + t376 * t516 - t461 * t573, (t348 * t567 + t358 * t383) * t416 + (t346 * t568 + t357 * t384) * t414 + (t344 * t569 + t356 * t382) * t412 + (t368 * t560 + t370 * t559 + t372 * t558) * t573, 1; t354, 0, 0, t350, t352, t363, t381, 0, 0, 0, t397 * t480 + t398 * t481 + t399 * t482 - t458 * t573, (t367 * t567 + t383 * t390) * t416 + (t366 * t568 + t384 * t389) * t414 + (t365 * t569 + t382 * t388) * t412 + (-qJ(3,1) * t471 - qJ(3,2) * t472 - qJ(3,3) * t473) * t573, 0; t355, 0, 0, t351, t353, t362, t380, 0, 0, 0, t377 * t555 + t378 * t550 + t379 * t545, (t349 * t449 + t361 * t405) * t416 + (t347 * t447 + t360 * t404) * t414 + (t345 * t445 + t359 * t403) * t412, 0; t354, 0, 0, t350, t352, t363, t381, 0, 0, 0, t374 * t555 + t375 * t550 + t376 * t545, (t348 * t449 + t358 * t405) * t416 + (t346 * t447 + t357 * t404) * t414 + (t344 * t445 + t356 * t403) * t412, 0; t540 + t541 + t542, 0, 0, t418 * t542 + t419 * t541 + t420 * t540, 0.2e1 * t438 * t444 * t542 + 0.2e1 * t440 * t446 * t541 + 0.2e1 * t442 * t448 * t540, 0, 0, 0, 0, 0, 0.2e1 * qJ(3,1) * t540 + 0.2e1 * qJ(3,2) * t541 + 0.2e1 * qJ(3,3) * t542, (t367 * t449 + t390 * t405) * t416 + (t366 * t447 + t389 * t404) * t414 + (t365 * t445 + t388 * t403) * t412, 1;];
tau_reg  = t1;
