% Calculate inertia matrix for parallel robot
% P4PRRRR1G3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% MX [4x4]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-02 19:06
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P4PRRRR1G3A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G3A0_inertia_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G3A0_inertia_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G3A0_inertia_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR1G3A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRRR1G3A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P4PRRRR1G3A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G3A0_inertia_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G3A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-02 19:00:48
% EndTime: 2020-03-02 19:00:51
% DurationCPUTime: 3.48s
% Computational Cost: add. (2450->268), mult. (4775->496), div. (1184->13), fcn. (3280->34), ass. (0->212)
t573 = m(3) / 0.2e1;
t572 = Icges(3,2) / 0.2e1;
t461 = cos(qJ(3,4));
t447 = 0.1e1 / t461;
t571 = t447 ^ 2;
t473 = cos(qJ(3,3));
t453 = 0.1e1 / t473;
t570 = t453 ^ 2;
t475 = cos(qJ(3,2));
t455 = 0.1e1 / t475;
t569 = t455 ^ 2;
t477 = cos(qJ(3,1));
t457 = 0.1e1 / t477;
t568 = t457 ^ 2;
t567 = rSges(3,3) * m(3);
t497 = rSges(3,2) ^ 2;
t498 = rSges(3,1) ^ 2;
t566 = (-t497 + t498) * t573 - Icges(3,1) / 0.2e1 + t572;
t463 = legFrame(4,2);
t436 = sin(t463);
t440 = cos(t463);
t460 = sin(qJ(2,4));
t462 = cos(qJ(2,4));
t406 = -t436 * t462 + t440 * t460;
t446 = 0.1e1 / t460;
t565 = t406 * t446;
t407 = t460 * t436 + t440 * t462;
t564 = t407 * t446;
t464 = legFrame(3,2);
t437 = sin(t464);
t441 = cos(t464);
t468 = sin(qJ(2,3));
t474 = cos(qJ(2,3));
t410 = -t437 * t474 + t441 * t468;
t450 = 0.1e1 / t468;
t563 = t410 * t450;
t411 = t468 * t437 + t441 * t474;
t562 = t411 * t450;
t465 = legFrame(2,2);
t438 = sin(t465);
t442 = cos(t465);
t470 = sin(qJ(2,2));
t476 = cos(qJ(2,2));
t412 = -t438 * t476 + t442 * t470;
t451 = 0.1e1 / t470;
t561 = t412 * t451;
t413 = t470 * t438 + t442 * t476;
t560 = t413 * t451;
t466 = legFrame(1,2);
t439 = sin(t466);
t443 = cos(t466);
t472 = sin(qJ(2,1));
t478 = cos(qJ(2,1));
t414 = -t439 * t478 + t443 * t472;
t452 = 0.1e1 / t472;
t559 = t414 * t452;
t415 = t472 * t439 + t443 * t478;
t558 = t415 * t452;
t459 = sin(qJ(3,4));
t424 = t459 * rSges(3,1) + t461 * rSges(3,2);
t557 = t424 * t447;
t467 = sin(qJ(3,3));
t425 = t467 * rSges(3,1) + t473 * rSges(3,2);
t556 = t425 * t453;
t469 = sin(qJ(3,2));
t426 = t469 * rSges(3,1) + t475 * rSges(3,2);
t555 = t426 * t455;
t471 = sin(qJ(3,1));
t427 = t471 * rSges(3,1) + t477 * rSges(3,2);
t554 = t427 * t457;
t553 = t446 * t447;
t552 = t446 * t459;
t499 = 0.1e1 / pkin(2);
t551 = t447 * t499;
t550 = 0.1e1 / t461 ^ 2 * t462;
t549 = t450 * t453;
t548 = t450 * t467;
t547 = t451 * t455;
t546 = t451 * t469;
t545 = t452 * t457;
t544 = t452 * t471;
t543 = t453 * t499;
t542 = 0.1e1 / t473 ^ 2 * t474;
t541 = t455 * t499;
t540 = 0.1e1 / t475 ^ 2 * t476;
t539 = t457 * t499;
t538 = 0.1e1 / t477 ^ 2 * t478;
t537 = t497 + t498;
t536 = m(3) * t424 * t460;
t535 = m(3) * t425 * t468;
t534 = m(3) * t426 * t470;
t533 = m(3) * t427 * t472;
t433 = -rSges(3,2) * t567 + Icges(3,6);
t434 = rSges(3,1) * t567 - Icges(3,5);
t402 = t433 * t461 - t434 * t459;
t532 = t402 * t446 * t571;
t403 = t433 * t473 - t434 * t467;
t531 = t403 * t450 * t570;
t404 = t433 * t475 - t434 * t469;
t530 = t404 * t451 * t569;
t405 = t433 * t477 - t434 * t471;
t529 = t405 * t452 * t568;
t528 = t436 * t553;
t527 = t436 * t551;
t526 = t437 * t549;
t525 = t437 * t543;
t524 = t438 * t547;
t523 = t438 * t541;
t522 = t439 * t545;
t521 = t439 * t539;
t520 = t440 * t553;
t519 = t440 * t551;
t518 = t441 * t549;
t517 = t441 * t543;
t516 = t442 * t547;
t515 = t442 * t541;
t514 = t443 * t545;
t513 = t443 * t539;
t512 = t447 * t552;
t511 = t499 * t550;
t510 = t453 * t548;
t509 = t455 * t546;
t508 = t457 * t544;
t507 = t499 * t542;
t506 = t499 * t540;
t505 = t499 * t538;
t504 = t550 * t552;
t503 = t542 * t548;
t502 = t540 * t546;
t501 = t538 * t544;
t500 = Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + (0.2e1 * rSges(3,3) ^ 2 + t537) * t573 + t572 + Icges(3,1) / 0.2e1;
t496 = koppelP(1,1);
t495 = koppelP(2,1);
t494 = koppelP(3,1);
t493 = koppelP(4,1);
t492 = koppelP(1,2);
t491 = koppelP(2,2);
t490 = koppelP(3,2);
t489 = koppelP(4,2);
t488 = 0.2e1 * qJ(3,1);
t487 = 0.2e1 * qJ(3,2);
t486 = 0.2e1 * qJ(3,3);
t485 = 0.2e1 * qJ(3,4);
t484 = rSges(4,1);
t483 = rSges(4,2);
t482 = xP(4);
t481 = m(2) * rSges(2,1);
t449 = m(1) + m(2) + m(3);
t445 = cos(t482);
t444 = sin(t482);
t435 = -m(3) * rSges(3,1) * rSges(3,2) + Icges(3,4);
t432 = m(2) * rSges(2,2) - t567;
t431 = t537 * m(3) + Icges(3,3);
t423 = -t444 * t492 + t445 * t496;
t422 = -t444 * t491 + t445 * t495;
t421 = -t444 * t490 + t445 * t494;
t420 = -t444 * t489 + t445 * t493;
t419 = -t444 * t496 - t445 * t492;
t418 = -t444 * t495 - t445 * t491;
t417 = -t444 * t494 - t445 * t490;
t416 = -t444 * t493 - t445 * t489;
t409 = m(4) * (-t444 * t483 + t445 * t484);
t408 = m(4) * (-t444 * t484 - t445 * t483);
t401 = (t481 + (rSges(3,1) * t477 - rSges(3,2) * t471) * m(3)) * t478 - t472 * t432;
t400 = (t481 + (rSges(3,1) * t475 - rSges(3,2) * t469) * m(3)) * t476 - t470 * t432;
t399 = (t481 + (rSges(3,1) * t473 - rSges(3,2) * t467) * m(3)) * t474 - t468 * t432;
t398 = (t481 + (rSges(3,1) * t461 - rSges(3,2) * t459) * m(3)) * t462 - t460 * t432;
t397 = cos(t488) * t566 + t435 * sin(t488) + t500;
t396 = cos(t487) * t566 + t435 * sin(t487) + t500;
t395 = cos(t486) * t566 + t435 * sin(t486) + t500;
t394 = cos(t485) * t566 + t435 * sin(t485) + t500;
t393 = (-t419 * t443 + t423 * t439) * t452 * t539;
t392 = (-t418 * t442 + t422 * t438) * t451 * t541;
t391 = (-t417 * t441 + t421 * t437) * t450 * t543;
t390 = (-t416 * t440 + t420 * t436) * t446 * t551;
t389 = (-t401 * t513 + t415 * t449) * t452;
t388 = (t401 * t521 + t414 * t449) * t452;
t387 = (-t400 * t515 + t413 * t449) * t451;
t386 = (t400 * t523 + t412 * t449) * t451;
t385 = (-t399 * t517 + t411 * t449) * t450;
t384 = (t399 * t525 + t410 * t449) * t450;
t383 = (-t398 * t519 + t407 * t449) * t446;
t382 = (t398 * t527 + t406 * t449) * t446;
t381 = (t414 * t423 + t415 * t419) * t452;
t380 = (t412 * t422 + t413 * t418) * t451;
t379 = (t410 * t421 + t411 * t417) * t450;
t378 = (t406 * t420 + t407 * t416) * t446;
t377 = -t533 * t539 + (-t401 * t505 + t449 * t457) * t544;
t376 = -t534 * t541 + (-t400 * t506 + t449 * t455) * t546;
t375 = -t535 * t543 + (-t399 * t507 + t449 * t453) * t548;
t374 = -t536 * t551 + (-t398 * t511 + t447 * t449) * t552;
t373 = (-t397 * t513 + t401 * t415) * t452;
t372 = (t397 * t521 + t401 * t414) * t452;
t371 = (-t396 * t515 + t400 * t413) * t451;
t370 = (t396 * t523 + t400 * t412) * t451;
t369 = (-t395 * t517 + t399 * t411) * t450;
t368 = (t395 * t525 + t399 * t410) * t450;
t367 = (-t394 * t519 + t398 * t407) * t446;
t366 = (t394 * t527 + t398 * t406) * t446;
t365 = t405 * t539 + (-t397 * t505 + t401 * t457) * t544;
t364 = t404 * t541 + (-t396 * t506 + t400 * t455) * t546;
t363 = t403 * t543 + (-t395 * t507 + t399 * t453) * t548;
t362 = t402 * t551 + (-t394 * t511 + t398 * t447) * t552;
t361 = t381 * t449 + t393 * t401;
t360 = t380 * t449 + t392 * t400;
t359 = t379 * t449 + t391 * t399;
t358 = t378 * t449 + t390 * t398;
t357 = t381 * t401 + t393 * t397;
t356 = t380 * t400 + t392 * t396;
t355 = t379 * t399 + t391 * t395;
t354 = t378 * t398 + t390 * t394;
t1 = [t383 * t564 + t385 * t562 + t387 * t560 + t389 * t558 + m(4) + (-t367 * t520 - t369 * t518 - t371 * t516 - t373 * t514) * t499, t383 * t565 + t385 * t563 + t387 * t561 + t389 * t559 + (t367 * t528 + t369 * t526 + t371 * t524 + t373 * t522) * t499, t383 * t512 + t385 * t510 + t387 * t509 + t389 * t508 + (-t373 * t501 - t371 * t502 - t369 * t503 - t367 * t504 + (-t440 * t532 - t441 * t531 - t442 * t530 - t443 * t529) * t499 + (-t407 * t557 - t411 * t556 - t413 * t555 - t415 * t554) * m(3)) * t499, t367 * t390 + t369 * t391 + t371 * t392 + t373 * t393 + t383 * t378 + t385 * t379 + t387 * t380 + t389 * t381 + t408; t382 * t564 + t384 * t562 + t386 * t560 + t388 * t558 + (-t366 * t520 - t368 * t518 - t370 * t516 - t372 * t514) * t499, t382 * t565 + t384 * t563 + t386 * t561 + t388 * t559 + m(4) + (t366 * t528 + t368 * t526 + t370 * t524 + t372 * t522) * t499, t382 * t512 + t384 * t510 + t386 * t509 + t388 * t508 + (-t372 * t501 - t370 * t502 - t368 * t503 - t366 * t504 + (t436 * t532 + t437 * t531 + t438 * t530 + t439 * t529) * t499 + (-t406 * t557 - t410 * t556 - t412 * t555 - t414 * t554) * m(3)) * t499, t366 * t390 + t368 * t391 + t370 * t392 + t372 * t393 + t382 * t378 + t384 * t379 + t386 * t380 + t388 * t381 + t409; t374 * t564 + t375 * t562 + t376 * t560 + t377 * t558 + (-t362 * t520 - t363 * t518 - t364 * t516 - t365 * t514) * t499, t374 * t565 + t375 * t563 + t376 * t561 + t377 * t559 + (t362 * t528 + t363 * t526 + t364 * t524 + t365 * t522) * t499, t374 * t512 + t375 * t510 + t376 * t509 + t377 * t508 + m(4) + (-t362 * t504 - t363 * t503 - t364 * t502 - t365 * t501 + ((-t405 * t501 + t431 * t457) * t457 + (-t404 * t502 + t431 * t455) * t455 + (-t403 * t503 + t431 * t453) * t453 + (-t402 * t504 + t431 * t447) * t447) * t499 + (-t424 * t571 * t459 - t425 * t570 * t467 - t426 * t569 * t469 - t427 * t568 * t471) * m(3)) * t499, t362 * t390 + t363 * t391 + t364 * t392 + t365 * t393 + t374 * t378 + t375 * t379 + t376 * t380 + t377 * t381; t358 * t564 + t359 * t562 + t360 * t560 + t361 * t558 + t408 + (-t354 * t520 - t355 * t518 - t356 * t516 - t357 * t514) * t499, t358 * t565 + t359 * t563 + t360 * t561 + t361 * t559 + t409 + (t354 * t528 + t355 * t526 + t356 * t524 + t357 * t522) * t499, t358 * t512 + t359 * t510 + t360 * t509 + t361 * t508 + (-t357 * t501 + (-t381 * t533 + t393 * t405) * t457 - t356 * t502 + (-t380 * t534 + t392 * t404) * t455 - t355 * t503 + (-t379 * t535 + t391 * t403) * t453 - t354 * t504 + (-t378 * t536 + t390 * t402) * t447) * t499, t361 * t381 + t357 * t393 + t360 * t380 + t356 * t392 + t359 * t379 + t355 * t391 + t358 * t378 + t354 * t390 + Icges(4,3) + m(4) * (t483 ^ 2 + t484 ^ 2);];
MX  = t1;
