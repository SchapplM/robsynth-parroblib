% Calculate inertia matrix for parallel robot
% P3PRRRR8V1G3A0
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
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 17:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRRR8V1G3A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G3A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G3A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G3A0_inertia_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G3A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V1G3A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR8V1G3A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G3A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G3A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:15:55
% EndTime: 2020-08-06 17:15:56
% DurationCPUTime: 1.38s
% Computational Cost: add. (3072->211), mult. (8334->417), div. (648->7), fcn. (7470->28), ass. (0->206)
t607 = m(3) / 0.2e1;
t606 = Icges(3,2) / 0.2e1;
t489 = sin(qJ(2,3));
t495 = cos(qJ(2,3));
t494 = cos(qJ(3,3));
t598 = pkin(2) * t494;
t454 = -t495 * pkin(5) + t489 * t598;
t482 = sin(pkin(3));
t484 = cos(pkin(3));
t488 = sin(qJ(3,3));
t567 = t488 * t484;
t448 = pkin(2) * t567 + t454 * t482;
t605 = 0.1e1 / t448;
t491 = sin(qJ(2,2));
t497 = cos(qJ(2,2));
t496 = cos(qJ(3,2));
t597 = pkin(2) * t496;
t455 = -t497 * pkin(5) + t491 * t597;
t490 = sin(qJ(3,2));
t565 = t490 * t484;
t449 = pkin(2) * t565 + t455 * t482;
t604 = 0.1e1 / t449;
t493 = sin(qJ(2,1));
t499 = cos(qJ(2,1));
t498 = cos(qJ(3,1));
t596 = pkin(2) * t498;
t456 = -t499 * pkin(5) + t493 * t596;
t492 = sin(qJ(3,1));
t563 = t492 * t484;
t450 = pkin(2) * t563 + t456 * t482;
t603 = 0.1e1 / t450;
t602 = m(3) * rSges(3,3);
t506 = rSges(3,2) ^ 2;
t507 = rSges(3,1) ^ 2;
t601 = (-t506 + t507) * t607 - Icges(3,1) / 0.2e1 + t606;
t508 = 0.1e1 / pkin(2);
t600 = m(3) * t508;
t599 = pkin(2) * t482;
t457 = pkin(5) * t489 + t495 * t598;
t481 = sin(pkin(6));
t483 = cos(pkin(6));
t523 = -t454 * t484 + t488 * t599;
t418 = t483 * t457 + t523 * t481;
t485 = legFrame(3,2);
t471 = sin(t485);
t474 = cos(t485);
t409 = t418 * t474 + t471 * t448;
t595 = t409 * t605;
t410 = -t418 * t471 + t474 * t448;
t594 = t410 * t605;
t458 = pkin(5) * t491 + t497 * t597;
t522 = -t455 * t484 + t490 * t599;
t419 = t483 * t458 + t522 * t481;
t486 = legFrame(2,2);
t472 = sin(t486);
t475 = cos(t486);
t411 = t419 * t475 + t472 * t449;
t593 = t411 * t604;
t412 = -t419 * t472 + t475 * t449;
t592 = t412 * t604;
t459 = pkin(5) * t493 + t499 * t596;
t521 = -t456 * t484 + t492 * t599;
t420 = t483 * t459 + t521 * t481;
t487 = legFrame(1,2);
t473 = sin(t487);
t476 = cos(t487);
t413 = t420 * t476 + t473 * t450;
t591 = t413 * t603;
t414 = -t420 * t473 + t476 * t450;
t590 = t414 * t603;
t421 = -t481 * t457 + t523 * t483;
t589 = t421 * t605;
t422 = -t481 * t458 + t522 * t483;
t588 = t422 * t604;
t423 = -t481 * t459 + t521 * t483;
t587 = t423 * t603;
t467 = m(2) * rSges(2,2) - t602;
t502 = m(2) * rSges(2,1);
t527 = rSges(3,1) * t494 - rSges(3,2) * t488;
t586 = ((t527 * m(3) + t502) * t495 - t467 * t489) * t482;
t526 = rSges(3,1) * t496 - rSges(3,2) * t490;
t585 = ((t526 * m(3) + t502) * t497 - t467 * t491) * t482;
t525 = rSges(3,1) * t498 - rSges(3,2) * t492;
t584 = ((t525 * m(3) + t502) * t499 - t467 * t493) * t482;
t583 = t605 / t494;
t582 = t604 / t496;
t581 = t603 / t498;
t477 = m(1) + m(2) + m(3);
t580 = t605 * t477;
t579 = t604 * t477;
t578 = t603 * t477;
t468 = -rSges(3,2) * t602 + Icges(3,6);
t469 = rSges(3,1) * t602 - Icges(3,5);
t451 = t468 * t494 - t469 * t488;
t577 = t451 * t508;
t452 = t468 * t496 - t469 * t490;
t576 = t452 * t508;
t453 = t468 * t498 - t469 * t492;
t575 = t453 * t508;
t561 = t506 + t507;
t574 = (t561 * m(3) + Icges(3,3)) * t508;
t573 = t484 * t489;
t572 = t484 * t491;
t571 = t484 * t493;
t570 = t484 * t495;
t569 = t484 * t497;
t568 = t484 * t499;
t566 = t488 * t495;
t564 = t490 * t497;
t562 = t492 * t499;
t436 = t527 * t484 - t482 * t489 * (t488 * rSges(3,1) + t494 * rSges(3,2));
t560 = m(3) * t436 * t605;
t559 = t436 * t600;
t437 = t526 * t484 - t482 * t491 * (t490 * rSges(3,1) + t496 * rSges(3,2));
t558 = m(3) * t437 * t604;
t557 = t437 * t600;
t438 = t525 * t484 - t482 * t493 * (t492 * rSges(3,1) + t498 * rSges(3,2));
t556 = m(3) * t438 * t603;
t555 = t438 * t600;
t427 = (t481 * t570 + t483 * t489) * t598 + (t481 * t573 - t483 * t495) * pkin(5);
t554 = t427 * t583;
t428 = (t481 * t569 + t483 * t491) * t597 + (t481 * t572 - t483 * t497) * pkin(5);
t553 = t428 * t582;
t429 = (t481 * t568 + t483 * t493) * t596 + (t481 * t571 - t483 * t499) * pkin(5);
t552 = t429 * t581;
t514 = t482 * t494 + t489 * t567;
t430 = t514 * t481 - t483 * t566;
t551 = t430 * t583;
t513 = t482 * t496 + t491 * t565;
t431 = t513 * t481 - t483 * t564;
t550 = t431 * t582;
t512 = t482 * t498 + t493 * t563;
t432 = t512 * t481 - t483 * t562;
t549 = t432 * t581;
t548 = t605 * t586;
t547 = t604 * t585;
t546 = t603 * t584;
t545 = t471 * t583;
t544 = t474 * t583;
t543 = t472 * t582;
t542 = t475 * t582;
t541 = t473 * t581;
t540 = t476 * t581;
t424 = (-t481 * t489 + t483 * t570) * t598 + pkin(5) * (t481 * t495 + t483 * t573);
t539 = t424 * t545;
t538 = t424 * t544;
t425 = (-t481 * t491 + t483 * t569) * t597 + pkin(5) * (t481 * t497 + t483 * t572);
t537 = t425 * t543;
t536 = t425 * t542;
t426 = (-t481 * t493 + t483 * t568) * t596 + pkin(5) * (t481 * t499 + t483 * t571);
t535 = t426 * t541;
t534 = t426 * t540;
t433 = t481 * t566 + t514 * t483;
t533 = t433 * t545;
t532 = t433 * t544;
t434 = t481 * t564 + t513 * t483;
t531 = t434 * t543;
t530 = t434 * t542;
t435 = t481 * t562 + t512 * t483;
t529 = t435 * t541;
t528 = t435 * t540;
t524 = Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + (0.2e1 * rSges(3,3) ^ 2 + t561) * t607 + t606 + Icges(3,1) / 0.2e1;
t470 = -m(3) * rSges(3,1) * rSges(3,2) + Icges(3,4);
t503 = 0.2e1 * qJ(3,3);
t415 = cos(t503) * t601 + t470 * sin(t503) + t524;
t520 = t415 * t433 + t424 * t577;
t504 = 0.2e1 * qJ(3,2);
t416 = cos(t504) * t601 + t470 * sin(t504) + t524;
t519 = t416 * t434 + t425 * t576;
t505 = 0.2e1 * qJ(3,1);
t417 = cos(t505) * t601 + t470 * sin(t505) + t524;
t518 = t417 * t435 + t426 * t575;
t517 = t424 * t574 + t433 * t451;
t516 = t425 * t574 + t434 * t452;
t515 = t426 * t574 + t435 * t453;
t511 = t424 * t559 + t433 * t586;
t510 = t425 * t557 + t434 * t585;
t509 = t426 * t555 + t435 * t584;
t408 = t423 * t556 + (t429 * t574 + t432 * t453) * t581;
t407 = t422 * t558 + (t428 * t574 + t431 * t452) * t582;
t406 = t421 * t560 + (t427 * t574 + t430 * t451) * t583;
t405 = t423 * t578 + (t429 * t555 + t432 * t584) * t581;
t404 = t422 * t579 + (t428 * t557 + t431 * t585) * t582;
t403 = t421 * t580 + (t427 * t559 + t430 * t586) * t583;
t402 = t423 * t546 + (t417 * t432 + t429 * t575) * t581;
t401 = t422 * t547 + (t416 * t431 + t428 * t576) * t582;
t400 = t421 * t548 + (t415 * t430 + t427 * t577) * t583;
t399 = t414 * t556 + t515 * t541;
t398 = t413 * t556 - t515 * t540;
t397 = t412 * t558 + t516 * t543;
t396 = t411 * t558 - t516 * t542;
t395 = t410 * t560 + t517 * t545;
t394 = t409 * t560 - t517 * t544;
t393 = t414 * t578 + t509 * t541;
t392 = t413 * t578 - t509 * t540;
t391 = t412 * t579 + t510 * t543;
t390 = t411 * t579 - t510 * t542;
t389 = t410 * t580 + t511 * t545;
t388 = t409 * t580 - t511 * t544;
t387 = t414 * t546 + t518 * t541;
t386 = t413 * t546 - t518 * t540;
t385 = t412 * t547 + t519 * t543;
t384 = t411 * t547 - t519 * t542;
t383 = t410 * t548 + t520 * t545;
t382 = t409 * t548 - t520 * t544;
t1 = [-t382 * t532 - t384 * t530 - t386 * t528 + t388 * t595 + t390 * t593 + t392 * t591 + m(4) + (-t394 * t538 - t396 * t536 - t398 * t534) * t508, t382 * t533 + t384 * t531 + t386 * t529 + t388 * t594 + t390 * t592 + t392 * t590 + (t394 * t539 + t396 * t537 + t398 * t535) * t508, t382 * t551 + t384 * t550 + t386 * t549 + t388 * t589 + t390 * t588 + t392 * t587 + (t394 * t554 + t396 * t553 + t398 * t552) * t508; -t383 * t532 - t385 * t530 - t387 * t528 + t389 * t595 + t391 * t593 + t393 * t591 + (-t395 * t538 - t397 * t536 - t399 * t534) * t508, t383 * t533 + t385 * t531 + t387 * t529 + t389 * t594 + t391 * t592 + t393 * t590 + m(4) + (t395 * t539 + t397 * t537 + t399 * t535) * t508, t383 * t551 + t385 * t550 + t387 * t549 + t389 * t589 + t391 * t588 + t393 * t587 + (t395 * t554 + t397 * t553 + t399 * t552) * t508; -t400 * t532 - t401 * t530 - t402 * t528 + t403 * t595 + t404 * t593 + t405 * t591 + (-t406 * t538 - t407 * t536 - t408 * t534) * t508, t400 * t533 + t401 * t531 + t402 * t529 + t403 * t594 + t404 * t592 + t405 * t590 + (t406 * t539 + t407 * t537 + t408 * t535) * t508, t400 * t551 + t401 * t550 + t402 * t549 + t403 * t589 + t404 * t588 + t405 * t587 + m(4) + (t406 * t554 + t407 * t553 + t408 * t552) * t508;];
MX  = t1;
