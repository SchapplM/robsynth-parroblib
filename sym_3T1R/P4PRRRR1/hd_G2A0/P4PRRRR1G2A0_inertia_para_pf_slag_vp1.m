% Calculate inertia matrix for parallel robot
% P4PRRRR1G2A0
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
% Datum: 2020-08-07 11:00
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P4PRRRR1G2A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G2A0_inertia_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G2A0_inertia_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G2A0_inertia_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR1G2A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRRR1G2A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P4PRRRR1G2A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G2A0_inertia_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G2A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 10:58:28
% EndTime: 2020-08-07 10:58:30
% DurationCPUTime: 1.96s
% Computational Cost: add. (3190->296), mult. (7047->526), div. (1884->13), fcn. (5376->34), ass. (0->244)
t673 = m(3) / 0.2e1;
t672 = Icges(3,2) / 0.2e1;
t671 = rSges(3,3) * m(3);
t581 = rSges(3,2) ^ 2;
t582 = rSges(3,1) ^ 2;
t670 = (-t581 + t582) * t673 - Icges(3,1) / 0.2e1 + t672;
t543 = sin(qJ(3,4));
t545 = cos(qJ(3,4));
t669 = m(3) * (rSges(3,1) * t543 + rSges(3,2) * t545);
t551 = sin(qJ(3,3));
t557 = cos(qJ(3,3));
t668 = m(3) * (rSges(3,1) * t551 + rSges(3,2) * t557);
t553 = sin(qJ(3,2));
t559 = cos(qJ(3,2));
t667 = m(3) * (rSges(3,1) * t553 + rSges(3,2) * t559);
t555 = sin(qJ(3,1));
t561 = cos(qJ(3,1));
t666 = m(3) * (rSges(3,1) * t555 + rSges(3,2) * t561);
t547 = legFrame(4,2);
t520 = sin(t547);
t583 = 0.1e1 / pkin(2);
t665 = t520 * t583;
t548 = legFrame(3,2);
t521 = sin(t548);
t664 = t521 * t583;
t549 = legFrame(2,2);
t522 = sin(t549);
t663 = t522 * t583;
t550 = legFrame(1,2);
t523 = sin(t550);
t662 = t523 * t583;
t524 = cos(t547);
t661 = t524 * t583;
t525 = cos(t548);
t660 = t525 * t583;
t526 = cos(t549);
t659 = t526 * t583;
t527 = cos(t550);
t658 = t527 * t583;
t544 = sin(qJ(2,4));
t530 = 0.1e1 / t544;
t531 = 0.1e1 / t545;
t657 = t530 * t531;
t546 = cos(qJ(2,4));
t656 = t530 * t546;
t655 = t531 * t583;
t552 = sin(qJ(2,3));
t534 = 0.1e1 / t552;
t537 = 0.1e1 / t557;
t654 = t534 * t537;
t558 = cos(qJ(2,3));
t653 = t534 * t558;
t554 = sin(qJ(2,2));
t535 = 0.1e1 / t554;
t539 = 0.1e1 / t559;
t652 = t535 * t539;
t560 = cos(qJ(2,2));
t651 = t535 * t560;
t556 = sin(qJ(2,1));
t536 = 0.1e1 / t556;
t541 = 0.1e1 / t561;
t650 = t536 * t541;
t562 = cos(qJ(2,1));
t649 = t536 * t562;
t648 = t537 * t583;
t647 = t539 * t583;
t646 = t541 * t583;
t645 = t544 * t545;
t644 = t552 * t557;
t643 = t554 * t559;
t642 = t556 * t561;
t641 = t581 + t582;
t640 = t531 * t669;
t639 = t544 * t669;
t638 = t537 * t668;
t637 = t552 * t668;
t636 = t539 * t667;
t635 = t554 * t667;
t634 = t541 * t666;
t633 = t556 * t666;
t486 = t520 * t645 - t524 * t543;
t632 = t486 * t657;
t487 = t520 * t543 + t524 * t645;
t631 = t487 * t657;
t492 = t521 * t644 - t525 * t551;
t630 = t492 * t654;
t493 = t522 * t643 - t526 * t553;
t629 = t493 * t652;
t494 = t523 * t642 - t527 * t555;
t628 = t494 * t650;
t495 = t521 * t551 + t525 * t644;
t627 = t495 * t654;
t496 = t522 * t553 + t526 * t643;
t626 = t496 * t652;
t497 = t523 * t555 + t527 * t642;
t625 = t497 * t650;
t533 = m(1) + m(2) + m(3);
t624 = t533 * t657;
t623 = t533 * t654;
t622 = t533 * t652;
t621 = t533 * t650;
t620 = 0.1e1 / t545 ^ 2 * t543 * t656;
t619 = 0.1e1 / t557 ^ 2 * t551 * t653;
t618 = 0.1e1 / t559 ^ 2 * t553 * t651;
t617 = 0.1e1 / t561 ^ 2 * t555 * t649;
t616 = Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + (0.2e1 * rSges(3,3) ^ 2 + t641) * t673 + t672 + Icges(3,1) / 0.2e1;
t566 = xP(4);
t528 = sin(t566);
t529 = cos(t566);
t573 = koppelP(4,2);
t577 = koppelP(4,1);
t500 = -t528 * t577 - t529 * t573;
t504 = -t528 * t573 + t529 * t577;
t615 = t500 * t524 - t504 * t520;
t574 = koppelP(3,2);
t578 = koppelP(3,1);
t501 = -t528 * t578 - t529 * t574;
t505 = -t528 * t574 + t529 * t578;
t614 = t501 * t525 - t505 * t521;
t575 = koppelP(2,2);
t579 = koppelP(2,1);
t502 = -t528 * t579 - t529 * t575;
t506 = -t528 * t575 + t529 * t579;
t613 = t502 * t526 - t506 * t522;
t576 = koppelP(1,2);
t580 = koppelP(1,1);
t503 = -t528 * t580 - t529 * t576;
t507 = -t528 * t576 + t529 * t580;
t612 = t503 * t527 - t507 * t523;
t458 = (t486 * t500 + t487 * t504) * t657;
t462 = t615 * t583 * t620;
t519 = -m(3) * rSges(3,1) * rSges(3,2) + Icges(3,4);
t569 = 0.2e1 * qJ(3,4);
t466 = cos(t569) * t670 + t519 * sin(t569) + t616;
t471 = t615 * t655;
t516 = m(2) * rSges(2,2) - t671;
t565 = m(2) * rSges(2,1);
t482 = (t565 + (rSges(3,1) * t545 - rSges(3,2) * t543) * m(3)) * t546 - t544 * t516;
t517 = -rSges(3,2) * t671 + Icges(3,6);
t518 = rSges(3,1) * t671 - Icges(3,5);
t488 = t517 * t545 - t518 * t543;
t418 = t458 * t482 + t462 * t466 - t471 * t488;
t515 = t641 * m(3) + Icges(3,3);
t423 = -t458 * t639 + t462 * t488 - t471 * t515;
t611 = t418 * t620 - t423 * t531;
t459 = (t492 * t501 + t495 * t505) * t654;
t463 = t614 * t583 * t619;
t570 = 0.2e1 * qJ(3,3);
t467 = cos(t570) * t670 + t519 * sin(t570) + t616;
t472 = t614 * t648;
t483 = (t565 + (rSges(3,1) * t557 - rSges(3,2) * t551) * m(3)) * t558 - t552 * t516;
t489 = t517 * t557 - t518 * t551;
t419 = t459 * t483 + t463 * t467 - t472 * t489;
t427 = -t459 * t637 + t463 * t489 - t472 * t515;
t610 = t419 * t619 - t427 * t537;
t460 = (t493 * t502 + t496 * t506) * t652;
t464 = t613 * t583 * t618;
t571 = 0.2e1 * qJ(3,2);
t468 = cos(t571) * t670 + t519 * sin(t571) + t616;
t473 = t613 * t647;
t484 = (t565 + (rSges(3,1) * t559 - rSges(3,2) * t553) * m(3)) * t560 - t554 * t516;
t490 = t517 * t559 - t518 * t553;
t420 = t460 * t484 + t464 * t468 - t473 * t490;
t428 = -t460 * t635 + t464 * t490 - t473 * t515;
t609 = t420 * t618 - t428 * t539;
t461 = (t494 * t503 + t497 * t507) * t650;
t465 = t612 * t583 * t617;
t572 = 0.2e1 * qJ(3,1);
t469 = cos(t572) * t670 + t519 * sin(t572) + t616;
t474 = t612 * t646;
t485 = (t565 + (rSges(3,1) * t561 - rSges(3,2) * t555) * m(3)) * t562 - t556 * t516;
t491 = t517 * t561 - t518 * t555;
t421 = t461 * t485 + t465 * t469 - t474 * t491;
t429 = -t461 * t633 + t465 * t491 - t474 * t515;
t608 = t421 * t617 - t429 * t541;
t595 = t466 * t620 - t488 * t531;
t430 = t482 * t632 + t595 * t661;
t591 = t488 * t620 - t515 * t531;
t450 = -t486 * t640 + t591 * t661;
t607 = t430 * t620 - t450 * t531;
t431 = t482 * t631 - t595 * t665;
t451 = -t487 * t640 - t591 * t665;
t606 = t431 * t620 - t451 * t531;
t594 = t467 * t619 - t489 * t537;
t432 = t483 * t630 + t594 * t660;
t590 = t489 * t619 - t515 * t537;
t452 = -t492 * t638 + t590 * t660;
t605 = t432 * t619 - t452 * t537;
t593 = t468 * t618 - t490 * t539;
t433 = t484 * t629 + t593 * t659;
t589 = t490 * t618 - t515 * t539;
t453 = -t493 * t636 + t589 * t659;
t604 = t433 * t618 - t453 * t539;
t592 = t469 * t617 - t491 * t541;
t434 = t485 * t628 + t592 * t658;
t588 = t491 * t617 - t515 * t541;
t454 = -t494 * t634 + t588 * t658;
t603 = t434 * t617 - t454 * t541;
t435 = t483 * t627 - t594 * t664;
t455 = -t495 * t638 - t590 * t664;
t602 = t435 * t619 - t455 * t537;
t436 = t484 * t626 - t593 * t663;
t456 = -t496 * t636 - t589 * t663;
t601 = t436 * t618 - t456 * t539;
t437 = t485 * t625 - t592 * t662;
t457 = -t497 * t634 - t588 * t662;
t600 = t437 * t617 - t457 * t541;
t446 = (-t466 * t655 + t482 * t546) * t530;
t478 = -t488 * t530 * t655 - t546 * t669;
t599 = t446 * t620 - t478 * t531;
t447 = (-t467 * t648 + t483 * t558) * t534;
t479 = -t489 * t534 * t648 - t558 * t668;
t598 = t447 * t619 - t479 * t537;
t448 = (-t468 * t647 + t484 * t560) * t535;
t480 = -t490 * t535 * t647 - t560 * t667;
t597 = t448 * t618 - t480 * t539;
t449 = (-t469 * t646 + t485 * t562) * t536;
t481 = -t491 * t536 * t646 - t562 * t666;
t596 = t449 * t617 - t481 * t541;
t587 = t482 * t620 + t531 * t639;
t586 = t483 * t619 + t537 * t637;
t585 = t484 * t618 + t539 * t635;
t584 = t485 * t617 + t541 * t633;
t568 = rSges(4,1);
t567 = rSges(4,2);
t499 = m(4) * (-t528 * t567 + t529 * t568);
t498 = m(4) * (-t528 * t568 - t529 * t567);
t477 = (-t485 * t646 + t533 * t562) * t536;
t476 = (-t484 * t647 + t533 * t560) * t535;
t475 = (-t483 * t648 + t533 * t558) * t534;
t470 = (-t482 * t655 + t533 * t546) * t530;
t445 = t497 * t621 - t584 * t662;
t444 = t496 * t622 - t585 * t663;
t443 = t495 * t623 - t586 * t664;
t442 = t494 * t621 + t584 * t658;
t441 = t493 * t622 + t585 * t659;
t440 = t492 * t623 + t586 * t660;
t439 = t487 * t624 - t587 * t665;
t438 = t486 * t624 + t587 * t661;
t426 = t461 * t533 + t465 * t485 + t474 * t633;
t425 = t460 * t533 + t464 * t484 + t473 * t635;
t424 = t459 * t533 + t463 * t483 + t472 * t637;
t422 = t458 * t533 + t462 * t482 + t471 * t639;
t1 = [t438 * t632 + t440 * t630 + t441 * t629 + t442 * t628 + m(4) + (t607 * t524 + t605 * t525 + t604 * t526 + t603 * t527) * t583, t438 * t631 + t440 * t627 + t441 * t626 + t442 * t625 + (-t607 * t520 - t605 * t521 - t604 * t522 - t603 * t523) * t583, t438 * t656 + t440 * t653 + t441 * t651 + t442 * t649 + (-t430 * t657 - t432 * t654 - t433 * t652 - t434 * t650) * t583, t430 * t462 + t432 * t463 + t433 * t464 + t434 * t465 + t438 * t458 + t440 * t459 + t441 * t460 + t442 * t461 - t450 * t471 - t452 * t472 - t453 * t473 - t454 * t474 + t498; t439 * t632 + t443 * t630 + t444 * t629 + t445 * t628 + (t606 * t524 + t602 * t525 + t601 * t526 + t600 * t527) * t583, t439 * t631 + t443 * t627 + t444 * t626 + t445 * t625 + m(4) + (-t606 * t520 - t602 * t521 - t601 * t522 - t600 * t523) * t583, t439 * t656 + t443 * t653 + t444 * t651 + t445 * t649 + (-t431 * t657 - t435 * t654 - t436 * t652 - t437 * t650) * t583, t431 * t462 + t435 * t463 + t436 * t464 + t437 * t465 + t439 * t458 + t443 * t459 + t444 * t460 + t445 * t461 - t451 * t471 - t455 * t472 - t456 * t473 - t457 * t474 + t499; t470 * t632 + t475 * t630 + t476 * t629 + t477 * t628 + (t599 * t524 + t598 * t525 + t597 * t526 + t596 * t527) * t583, t470 * t631 + t475 * t627 + t476 * t626 + t477 * t625 + (-t599 * t520 - t598 * t521 - t597 * t522 - t596 * t523) * t583, t470 * t656 + t475 * t653 + t476 * t651 + t477 * t649 + m(4) + (-t446 * t657 - t447 * t654 - t448 * t652 - t449 * t650) * t583, t446 * t462 + t447 * t463 + t448 * t464 + t449 * t465 + t458 * t470 + t459 * t475 + t460 * t476 + t461 * t477 - t471 * t478 - t472 * t479 - t473 * t480 - t474 * t481; t422 * t632 + t424 * t630 + t425 * t629 + t426 * t628 + t498 + (t611 * t524 + t610 * t525 + t609 * t526 + t608 * t527) * t583, t422 * t631 + t424 * t627 + t425 * t626 + t426 * t625 + t499 + (-t611 * t520 - t610 * t521 - t609 * t522 - t608 * t523) * t583, t422 * t656 + t424 * t653 + t425 * t651 + t426 * t649 + (-t418 * t657 - t419 * t654 - t420 * t652 - t421 * t650) * t583, t426 * t461 + t421 * t465 - t429 * t474 + t425 * t460 + t420 * t464 - t428 * t473 + t424 * t459 + t419 * t463 - t427 * t472 + t422 * t458 + t418 * t462 - t423 * t471 + Icges(4,3) + m(4) * (t567 ^ 2 + t568 ^ 2);];
MX  = t1;
