% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RPRR1G3A0
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
% Datum: 2020-03-09 21:27
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RPRR1G3A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G3A0_inertia_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G3A0_inertia_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G3A0_inertia_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G3A0_inertia_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G3A0_inertia_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:27:00
% EndTime: 2020-03-09 21:27:02
% DurationCPUTime: 2.14s
% Computational Cost: add. (5700->262), mult. (3519->416), div. (669->7), fcn. (2937->62), ass. (0->222)
t485 = legFrame(3,2);
t559 = qJ(1,3) + pkin(7);
t467 = t485 + t559;
t461 = qJ(3,3) + t467;
t516 = -t485 + t559;
t506 = qJ(3,3) + t516;
t425 = cos(t506) + cos(t461);
t488 = sin(qJ(3,3));
t556 = pkin(7) + qJ(3,3);
t440 = pkin(1) * sin(t556) + t488 * pkin(2);
t434 = 0.1e1 / t440;
t578 = t434 / 0.2e1;
t632 = t425 * t578;
t487 = legFrame(1,2);
t561 = qJ(1,1) + pkin(7);
t469 = t487 + t561;
t463 = qJ(3,1) + t469;
t518 = -t487 + t561;
t504 = qJ(3,1) + t518;
t427 = cos(t504) + cos(t463);
t490 = sin(qJ(3,1));
t558 = pkin(7) + qJ(3,1);
t442 = pkin(1) * sin(t558) + t490 * pkin(2);
t438 = 0.1e1 / t442;
t574 = t438 / 0.2e1;
t631 = t427 * t574;
t607 = sin(t463) - sin(t504);
t612 = -t607 * t438 / 0.2e1;
t609 = sin(t461) - sin(t506);
t610 = -t609 * t434 / 0.2e1;
t435 = 0.1e1 / t440 ^ 2;
t630 = t435 / 0.4e1;
t489 = sin(qJ(3,2));
t557 = pkin(7) + qJ(3,2);
t441 = pkin(1) * sin(t557) + t489 * pkin(2);
t437 = 0.1e1 / t441 ^ 2;
t629 = t437 / 0.4e1;
t439 = 0.1e1 / t442 ^ 2;
t628 = t439 / 0.4e1;
t496 = 0.1e1 / pkin(3);
t566 = t496 / 0.2e1;
t472 = sin(qJ(1,1) + t558);
t421 = pkin(3) * t472 + pkin(2) * sin(t561) + sin(qJ(1,1)) * pkin(1);
t594 = t421 * t438;
t409 = t496 * t594;
t570 = t472 * t438;
t394 = -t570 + t409 / 0.2e1;
t484 = cos(pkin(7));
t473 = t484 * pkin(1) + pkin(2);
t493 = cos(qJ(3,1));
t483 = sin(pkin(7));
t606 = pkin(1) * t483;
t433 = t490 * t473 + t493 * t606;
t580 = t433 * t438;
t627 = t394 * t580;
t471 = sin(qJ(1,2) + t557);
t560 = qJ(1,2) + pkin(7);
t420 = pkin(3) * t471 + pkin(2) * sin(t560) + sin(qJ(1,2)) * pkin(1);
t436 = 0.1e1 / t441;
t595 = t420 * t436;
t408 = t496 * t595;
t571 = t471 * t436;
t393 = -t571 + t408 / 0.2e1;
t492 = cos(qJ(3,2));
t432 = t489 * t473 + t492 * t606;
t582 = t432 * t436;
t626 = t393 * t582;
t470 = sin(qJ(1,3) + t556);
t419 = pkin(3) * t470 + pkin(2) * sin(t559) + sin(qJ(1,3)) * pkin(1);
t596 = t419 * t434;
t407 = t496 * t596;
t572 = t470 * t434;
t392 = -t572 + t407 / 0.2e1;
t491 = cos(qJ(3,3));
t431 = t488 * t473 + t491 * t606;
t584 = t431 * t434;
t625 = t392 * t584;
t567 = t483 * t490;
t412 = t493 * pkin(2) + (t484 * t493 - t567) * pkin(1);
t624 = t412 * t631;
t623 = t412 * t612;
t569 = t483 * t488;
t411 = t491 * pkin(2) + (t484 * t491 - t569) * pkin(1);
t622 = t411 * t632;
t621 = t411 * t610;
t568 = t483 * t489;
t410 = -t492 * pkin(2) + (-t484 * t492 + t568) * pkin(1);
t486 = legFrame(2,2);
t468 = t486 + t560;
t462 = qJ(3,2) + t468;
t517 = -t486 + t560;
t505 = qJ(3,2) + t517;
t426 = cos(t505) + cos(t462);
t589 = t426 * t436;
t620 = -t410 * t589 / 0.2e1;
t608 = sin(t462) - sin(t505);
t592 = t608 * t436;
t619 = t410 * t592 / 0.2e1;
t573 = t439 * t472;
t618 = -t573 / 0.2e1;
t575 = t437 * t471;
t617 = -t575 / 0.2e1;
t577 = t435 * t470;
t616 = -t577 / 0.2e1;
t614 = t589 / 0.2e1;
t611 = -t592 / 0.2e1;
t474 = qJ(1,3) + t485;
t562 = qJ(1,3) - t485;
t386 = t609 * pkin(3) + (sin(t467) - sin(t516)) * pkin(2) + (sin(t474) - sin(t562)) * pkin(1);
t605 = t386 * t434;
t475 = qJ(1,2) + t486;
t563 = qJ(1,2) - t486;
t387 = t608 * pkin(3) + (sin(t468) - sin(t517)) * pkin(2) + (sin(t475) - sin(t563)) * pkin(1);
t604 = t387 * t436;
t476 = qJ(1,1) + t487;
t564 = qJ(1,1) - t487;
t388 = t607 * pkin(3) + (sin(t469) - sin(t518)) * pkin(2) + (sin(t476) - sin(t564)) * pkin(1);
t603 = t388 * t438;
t389 = -t425 * pkin(3) + (-cos(t516) - cos(t467)) * pkin(2) + (-cos(t562) - cos(t474)) * pkin(1);
t602 = t389 * t434;
t390 = -t426 * pkin(3) + (-cos(t517) - cos(t468)) * pkin(2) + (-cos(t563) - cos(t475)) * pkin(1);
t601 = t390 * t436;
t391 = -t427 * pkin(3) + (-cos(t518) - cos(t469)) * pkin(2) + (-cos(t564) - cos(t476)) * pkin(1);
t600 = t391 * t438;
t395 = t407 - t572;
t599 = t395 * t434;
t397 = t408 - t571;
t598 = t397 * t436;
t399 = t409 - t570;
t597 = t399 * t438;
t428 = -pkin(1) * t569 + t491 * t473;
t587 = t428 * t435;
t429 = -pkin(1) * t568 + t492 * t473;
t586 = t429 * t437;
t430 = -pkin(1) * t567 + t493 * t473;
t585 = t430 * t439;
t583 = t431 * t435;
t581 = t432 * t437;
t579 = t433 * t439;
t576 = t436 / 0.2e1;
t565 = t496 / 0.4e1;
t552 = t410 * t571;
t549 = t411 * t572;
t546 = t412 * t570;
t545 = t609 * t587;
t544 = t609 * t584;
t543 = t609 * t583;
t542 = t608 * t586;
t541 = t608 * t582;
t540 = t608 * t581;
t539 = t607 * t585;
t538 = t607 * t580;
t537 = t607 * t579;
t536 = t425 * t587;
t535 = t425 * t584;
t534 = t425 * t583;
t533 = t426 * t586;
t532 = t426 * t582;
t531 = t426 * t581;
t530 = t427 * t585;
t529 = t427 * t580;
t528 = t427 * t579;
t527 = t428 * t577;
t526 = t429 * t575;
t525 = t430 * t573;
t524 = t431 * t577;
t523 = t432 * t575;
t522 = t433 * t573;
t521 = t434 * t566;
t520 = t436 * t566;
t519 = t438 * t566;
t512 = 0.2e1 * t431 * t572;
t511 = 0.2e1 * t432 * t571;
t510 = 0.2e1 * t433 * t570;
t509 = t386 * t521;
t508 = t387 * t520;
t507 = t388 * t519;
t500 = t435 * t470 ^ 2 + t437 * t471 ^ 2 + t439 * t472 ^ 2;
t499 = t607 ^ 2 * t628 + t608 ^ 2 * t629 + t609 ^ 2 * t630;
t498 = t425 ^ 2 * t630 + t426 ^ 2 * t629 + t427 ^ 2 * t628;
t360 = -t425 * t609 * t630 - t426 * t608 * t629 - t427 * t607 * t628;
t363 = -t607 * t618 - t608 * t617 - t609 * t616;
t364 = t425 * t616 + t426 * t617 + t427 * t618;
t497 = pkin(1) ^ 2;
t482 = cos(t487);
t481 = cos(t486);
t480 = cos(t485);
t479 = sin(t487);
t478 = sin(t486);
t477 = sin(t485);
t405 = t426 * t576;
t403 = t607 * t574;
t402 = t608 * t576;
t401 = t609 * t578;
t400 = t409 - 0.2e1 * t570;
t398 = t408 - 0.2e1 * t571;
t396 = t407 - 0.2e1 * t572;
t385 = t391 * t519;
t384 = t390 * t520;
t383 = t389 * t521;
t382 = t631 + t385;
t381 = 0.2e1 * t631 + t385;
t380 = t405 + t384;
t379 = 0.2e1 * t405 + t384;
t378 = t632 + t383;
t377 = 0.2e1 * t632 + t383;
t376 = -t403 + t507;
t375 = -0.2e1 * t403 + t507;
t374 = -t402 + t508;
t373 = -0.2e1 * t402 + t508;
t372 = -t401 + t509;
t371 = -0.2e1 * t401 + t509;
t370 = t631 + t385 / 0.2e1;
t369 = t405 + t384 / 0.2e1;
t368 = t632 + t383 / 0.2e1;
t367 = -t403 + t507 / 0.2e1;
t366 = -t402 + t508 / 0.2e1;
t365 = -t401 + t509 / 0.2e1;
t362 = t497 * t364;
t361 = t497 * t363;
t359 = t497 * t360 + t477 * t480 + t478 * t481 + t479 * t482;
t1 = [t498, 0, 0, t477 ^ 2 + t478 ^ 2 + t479 ^ 2 + t497 * t498, t378 * t632 + t380 * t614 + t382 * t631 + (t378 * t602 + t380 * t601 + t382 * t600) * t566, (t389 * t536 + t390 * t533 + t391 * t530) * t565 + t377 * t622 + t379 * t620 + t381 * t624, (-t389 * t534 - t390 * t531 - t391 * t528) * t565 - t368 * t535 - t369 * t532 - t370 * t529, 1; t360, 0, 0, t359, t372 * t632 + t374 * t614 + t376 * t631 + (t372 * t602 + t374 * t601 + t376 * t600) * t566, (-t389 * t545 - t390 * t542 - t391 * t539) * t565 + t371 * t622 + t373 * t620 + t375 * t624, (t389 * t543 + t390 * t540 + t391 * t537) * t565 - t365 * t535 - t366 * t532 - t367 * t529, 0; t364, 0, 0, t362, t395 * t632 + t397 * t614 + t399 * t631 + (t389 * t599 + t390 * t598 + t391 * t597) * t566, t396 * t622 + t398 * t620 + t400 * t624 + (-t389 * t527 - t390 * t526 - t391 * t525) * t566, -t425 * t625 - t426 * t626 - t427 * t627 + (t389 * t524 + t390 * t523 + t391 * t522) * t566, 0; t360, 0, 0, t359, t378 * t610 + t380 * t611 + t382 * t612 + (t378 * t605 + t380 * t604 + t382 * t603) * t566, (t386 * t536 + t387 * t533 + t388 * t530) * t565 + t377 * t621 + t379 * t619 + t381 * t623, (-t386 * t534 - t387 * t531 - t388 * t528) * t565 + t368 * t544 + t369 * t541 + t370 * t538, 0; t499, 0, 0, t480 ^ 2 + t481 ^ 2 + t482 ^ 2 + t497 * t499, t372 * t610 + t374 * t611 + t376 * t612 + (t372 * t605 + t374 * t604 + t376 * t603) * t566, (-t386 * t545 - t387 * t542 - t388 * t539) * t565 + t371 * t621 + t373 * t619 + t375 * t623, (t386 * t543 + t387 * t540 + t388 * t537) * t565 + t365 * t544 + t366 * t541 + t367 * t538, 1; t363, 0, 0, t361, t395 * t610 + t397 * t611 + t399 * t612 + (t386 * t599 + t387 * t598 + t388 * t597) * t566, t396 * t621 + t398 * t619 + t400 * t623 + (-t386 * t527 - t387 * t526 - t388 * t525) * t566, t609 * t625 + t608 * t626 + t607 * t627 + (t386 * t524 + t387 * t523 + t388 * t522) * t566, 0; t364, 0, 0, t362, -t378 * t572 - t380 * t571 - t382 * t570 + (t378 * t596 + t380 * t595 + t382 * t594) * t496, -t377 * t549 + t379 * t552 - t381 * t546 + (t419 * t536 + t420 * t533 + t421 * t530) * t566, t368 * t512 + t369 * t511 + t370 * t510 + (-t419 * t534 - t420 * t531 - t421 * t528) * t566, 0; t363, 0, 0, t361, -t372 * t572 - t374 * t571 - t376 * t570 + (t372 * t596 + t374 * t595 + t376 * t594) * t496, -t371 * t549 + t373 * t552 - t375 * t546 + (-t419 * t545 - t420 * t542 - t421 * t539) * t566, t365 * t512 + t366 * t511 + t367 * t510 + (t419 * t543 + t420 * t540 + t421 * t537) * t566, 0; t500, 0, 0, t500 * t497, -t395 * t572 - t397 * t571 - t399 * t570 + (t395 * t596 + t397 * t595 + t399 * t594) * t496, -t396 * t549 + t398 * t552 - t400 * t546 + (-t419 * t527 - t420 * t526 - t421 * t525) * t496, t392 * t512 + t393 * t511 + t394 * t510 + (t419 * t524 + t420 * t523 + t421 * t522) * t496, 1;];
tau_reg  = t1;
