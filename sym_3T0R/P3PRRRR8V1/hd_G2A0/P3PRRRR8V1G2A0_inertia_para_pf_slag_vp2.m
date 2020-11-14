% Calculate inertia matrix for parallel robot
% P3PRRRR8V1G2A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 17:03
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRRR8V1G2A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G2A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G2A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G2A0_inertia_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G2A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V1G2A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR8V1G2A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G2A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G2A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:02:51
% EndTime: 2020-08-06 17:02:52
% DurationCPUTime: 1.15s
% Computational Cost: add. (2667->199), mult. (7119->393), div. (648->7), fcn. (7497->22), ass. (0->193)
t564 = 2 * Ifges(3,4);
t459 = sin(qJ(2,3));
t465 = cos(qJ(2,3));
t464 = cos(qJ(3,3));
t559 = pkin(2) * t464;
t424 = -t465 * pkin(5) + t459 * t559;
t450 = sin(pkin(3));
t452 = cos(pkin(3));
t458 = sin(qJ(3,3));
t521 = t458 * t452;
t415 = pkin(2) * t521 + t424 * t450;
t563 = 0.1e1 / t415;
t461 = sin(qJ(2,2));
t467 = cos(qJ(2,2));
t466 = cos(qJ(3,2));
t558 = pkin(2) * t466;
t425 = -t467 * pkin(5) + t461 * t558;
t460 = sin(qJ(3,2));
t519 = t460 * t452;
t416 = pkin(2) * t519 + t425 * t450;
t562 = 0.1e1 / t416;
t463 = sin(qJ(2,1));
t469 = cos(qJ(2,1));
t468 = cos(qJ(3,1));
t557 = pkin(2) * t468;
t426 = -t469 * pkin(5) + t463 * t557;
t462 = sin(qJ(3,1));
t517 = t462 * t452;
t417 = pkin(2) * t517 + t426 * t450;
t561 = 0.1e1 / t417;
t560 = pkin(2) * t450;
t556 = Ifges(3,1) + Ifges(2,3);
t470 = 0.1e1 / pkin(2);
t555 = Ifges(3,3) * t470;
t427 = pkin(5) * t459 + t465 * t559;
t449 = sin(pkin(6));
t451 = cos(pkin(6));
t485 = -t424 * t452 + t458 * t560;
t388 = -t449 * t427 + t485 * t451;
t455 = legFrame(3,2);
t439 = sin(t455);
t442 = cos(t455);
t382 = t388 * t439 + t415 * t442;
t554 = t382 * t563;
t383 = -t388 * t442 + t415 * t439;
t553 = t383 * t563;
t428 = pkin(5) * t461 + t467 * t558;
t484 = -t425 * t452 + t460 * t560;
t389 = -t449 * t428 + t484 * t451;
t456 = legFrame(2,2);
t440 = sin(t456);
t443 = cos(t456);
t384 = t389 * t440 + t416 * t443;
t552 = t384 * t562;
t385 = -t389 * t443 + t416 * t440;
t551 = t385 * t562;
t429 = pkin(5) * t463 + t469 * t557;
t483 = -t426 * t452 + t462 * t560;
t390 = -t449 * t429 + t483 * t451;
t457 = legFrame(1,2);
t441 = sin(t457);
t444 = cos(t457);
t386 = t390 * t441 + t417 * t444;
t550 = t386 * t561;
t387 = -t390 * t444 + t417 * t441;
t549 = t387 * t561;
t391 = t451 * t427 + t485 * t449;
t548 = t391 * t563;
t392 = t451 * t428 + t484 * t449;
t547 = t392 * t562;
t393 = t451 * t429 + t483 * t449;
t546 = t393 * t561;
t500 = t464 * mrSges(3,1) - t458 * mrSges(3,2);
t406 = t500 * t452 - t459 * (t458 * mrSges(3,1) + t464 * mrSges(3,2)) * t450;
t545 = t406 * t563;
t544 = t406 * t470;
t499 = t466 * mrSges(3,1) - t460 * mrSges(3,2);
t407 = t499 * t452 - t461 * (t460 * mrSges(3,1) + t466 * mrSges(3,2)) * t450;
t543 = t407 * t562;
t542 = t407 * t470;
t498 = t468 * mrSges(3,1) - t462 * mrSges(3,2);
t408 = t498 * t452 - t463 * (t462 * mrSges(3,1) + t468 * mrSges(3,2)) * t450;
t541 = t408 * t561;
t540 = t408 * t470;
t539 = t563 / t464;
t538 = t562 / t466;
t537 = t561 / t468;
t445 = m(1) + m(2) + m(3);
t536 = t563 * t445;
t535 = t562 * t445;
t534 = t561 * t445;
t454 = mrSges(2,2) - mrSges(3,3);
t533 = ((mrSges(2,1) + t500) * t465 - t459 * t454) * t450;
t532 = ((mrSges(2,1) + t499) * t467 - t461 * t454) * t450;
t531 = ((mrSges(2,1) + t498) * t469 - t463 * t454) * t450;
t430 = Ifges(3,5) * t458 + Ifges(3,6) * t464;
t530 = t430 * t470;
t431 = Ifges(3,5) * t460 + Ifges(3,6) * t466;
t529 = t431 * t470;
t432 = Ifges(3,5) * t462 + Ifges(3,6) * t468;
t528 = t432 * t470;
t527 = t452 * t459;
t526 = t452 * t461;
t525 = t452 * t463;
t524 = t452 * t465;
t523 = t452 * t467;
t522 = t452 * t469;
t520 = t458 * t465;
t518 = t460 * t467;
t516 = t462 * t469;
t394 = -(-t449 * t459 + t451 * t524) * t559 - pkin(5) * (t449 * t465 + t451 * t527);
t515 = t394 * t539;
t395 = -(-t449 * t461 + t451 * t523) * t558 - pkin(5) * (t449 * t467 + t451 * t526);
t514 = t395 * t538;
t396 = -(-t449 * t463 + t451 * t522) * t557 - pkin(5) * (t449 * t469 + t451 * t525);
t513 = t396 * t537;
t476 = t450 * t464 + t459 * t521;
t403 = -t449 * t520 - t476 * t451;
t512 = t403 * t539;
t475 = t450 * t466 + t461 * t519;
t404 = -t449 * t518 - t475 * t451;
t511 = t404 * t538;
t474 = t450 * t468 + t463 * t517;
t405 = -t449 * t516 - t474 * t451;
t510 = t405 * t537;
t509 = t439 * t539;
t508 = t442 * t539;
t507 = t440 * t538;
t506 = t443 * t538;
t505 = t441 * t537;
t504 = t444 * t537;
t503 = t563 * t533;
t502 = t562 * t532;
t501 = t561 * t531;
t397 = (t449 * t524 + t451 * t459) * t559 + (t449 * t527 - t451 * t465) * pkin(5);
t497 = t397 * t509;
t496 = t397 * t508;
t398 = (t449 * t523 + t451 * t461) * t558 + (t449 * t526 - t451 * t467) * pkin(5);
t495 = t398 * t507;
t494 = t398 * t506;
t399 = (t449 * t522 + t451 * t463) * t557 + (t449 * t525 - t451 * t469) * pkin(5);
t493 = t399 * t505;
t492 = t399 * t504;
t400 = t476 * t449 - t451 * t520;
t491 = t400 * t509;
t490 = t400 * t508;
t401 = t475 * t449 - t451 * t518;
t489 = t401 * t507;
t488 = t401 * t506;
t402 = t474 * t449 - t451 * t516;
t487 = t402 * t505;
t486 = t402 * t504;
t482 = t397 * t555 + t400 * t430;
t481 = t398 * t555 + t401 * t431;
t480 = t399 * t555 + t402 * t432;
t453 = -Ifges(3,1) + Ifges(3,2);
t421 = (t453 * t464 + t458 * t564) * t464 + t556;
t479 = t397 * t530 + t400 * t421;
t422 = (t453 * t466 + t460 * t564) * t466 + t556;
t478 = t398 * t529 + t401 * t422;
t423 = (t453 * t468 + t462 * t564) * t468 + t556;
t477 = t399 * t528 + t402 * t423;
t473 = t397 * t544 + t400 * t533;
t472 = t398 * t542 + t401 * t532;
t471 = t399 * t540 + t402 * t531;
t381 = t393 * t541 + (t396 * t555 + t405 * t432) * t537;
t380 = t392 * t543 + (t395 * t555 + t404 * t431) * t538;
t379 = t391 * t545 + (t394 * t555 + t403 * t430) * t539;
t378 = t393 * t501 + (t396 * t528 + t405 * t423) * t537;
t377 = t392 * t502 + (t395 * t529 + t404 * t422) * t538;
t376 = t391 * t503 + (t394 * t530 + t403 * t421) * t539;
t375 = t393 * t534 + (t396 * t540 + t405 * t531) * t537;
t374 = t392 * t535 + (t395 * t542 + t404 * t532) * t538;
t373 = t391 * t536 + (t394 * t544 + t403 * t533) * t539;
t372 = t387 * t541 - t480 * t504;
t371 = t386 * t541 + t480 * t505;
t370 = t385 * t543 - t481 * t506;
t369 = t384 * t543 + t481 * t507;
t368 = t383 * t545 - t482 * t508;
t367 = t382 * t545 + t482 * t509;
t366 = t387 * t501 - t477 * t504;
t365 = t386 * t501 + t477 * t505;
t364 = t385 * t502 - t478 * t506;
t363 = t384 * t502 + t478 * t507;
t362 = t383 * t503 - t479 * t508;
t361 = t382 * t503 + t479 * t509;
t360 = t387 * t534 - t471 * t504;
t359 = t386 * t534 + t471 * t505;
t358 = t385 * t535 - t472 * t506;
t357 = t384 * t535 + t472 * t507;
t356 = t383 * t536 - t473 * t508;
t355 = t382 * t536 + t473 * t509;
t1 = [-t362 * t490 - t364 * t488 - t366 * t486 + t356 * t553 + t358 * t551 + t360 * t549 + m(4) + (-t368 * t496 - t370 * t494 - t372 * t492) * t470, t362 * t491 + t364 * t489 + t366 * t487 + t356 * t554 + t358 * t552 + t360 * t550 + (t368 * t497 + t370 * t495 + t372 * t493) * t470, t362 * t512 + t364 * t511 + t366 * t510 + t356 * t548 + t358 * t547 + t360 * t546 + (t368 * t515 + t370 * t514 + t372 * t513) * t470; -t361 * t490 - t363 * t488 - t365 * t486 + t355 * t553 + t357 * t551 + t359 * t549 + (-t367 * t496 - t369 * t494 - t371 * t492) * t470, t361 * t491 + t363 * t489 + t365 * t487 + t355 * t554 + t357 * t552 + t359 * t550 + m(4) + (t367 * t497 + t369 * t495 + t371 * t493) * t470, t361 * t512 + t363 * t511 + t365 * t510 + t355 * t548 + t357 * t547 + t359 * t546 + (t367 * t515 + t369 * t514 + t371 * t513) * t470; -t376 * t490 - t377 * t488 - t378 * t486 + t373 * t553 + t374 * t551 + t375 * t549 + (-t379 * t496 - t380 * t494 - t381 * t492) * t470, t376 * t491 + t377 * t489 + t378 * t487 + t373 * t554 + t374 * t552 + t375 * t550 + (t379 * t497 + t380 * t495 + t381 * t493) * t470, t376 * t512 + t377 * t511 + t378 * t510 + t373 * t548 + t374 * t547 + t375 * t546 + m(4) + (t379 * t515 + t380 * t514 + t381 * t513) * t470;];
MX  = t1;
