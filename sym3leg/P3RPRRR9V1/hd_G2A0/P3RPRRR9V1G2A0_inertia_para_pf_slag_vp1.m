% Calculate inertia matrix for parallel robot
% P3RPRRR9V1G2A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-08-06 18:53
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RPRRR9V1G2A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G2A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G2A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G2A0_inertia_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR9V1G2A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR9V1G2A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRRR9V1G2A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G2A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G2A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:51:23
% EndTime: 2020-08-06 18:51:25
% DurationCPUTime: 1.46s
% Computational Cost: add. (4302->275), mult. (6591->480), div. (396->10), fcn. (3540->40), ass. (0->224)
t616 = 2 * pkin(1);
t686 = Icges(3,2) / 0.2e1;
t685 = 2 * pkin(3);
t684 = 4 * rSges(2,3);
t562 = cos(pkin(7));
t683 = 0.2e1 * t562 ^ 2;
t682 = m(2) / 0.2e1;
t681 = rSges(2,2) * m(2);
t588 = rSges(3,2) ^ 2;
t590 = rSges(3,1) ^ 2;
t680 = (-t588 + t590) * m(3) / 0.2e1 - Icges(3,1) / 0.2e1 + t686;
t668 = pkin(5) + qJ(2,3);
t545 = rSges(3,3) + t668;
t679 = m(3) * t545;
t669 = pkin(5) + qJ(2,2);
t546 = rSges(3,3) + t669;
t678 = m(3) * t546;
t670 = pkin(5) + qJ(2,1);
t547 = rSges(3,3) + t670;
t677 = m(3) * t547;
t572 = cos(qJ(3,3));
t558 = t572 ^ 2;
t676 = t558 * pkin(3);
t574 = cos(qJ(3,2));
t559 = t574 ^ 2;
t675 = t559 * pkin(3);
t576 = cos(qJ(3,1));
t560 = t576 ^ 2;
t674 = t560 * pkin(3);
t673 = t572 * pkin(2);
t672 = t574 * pkin(2);
t671 = t576 * pkin(2);
t567 = sin(qJ(1,3));
t548 = pkin(6) + t668;
t573 = cos(qJ(1,3));
t606 = pkin(1) * t567 - t573 * t548;
t561 = sin(pkin(7));
t566 = sin(qJ(3,3));
t622 = t561 * t566;
t475 = t606 * t622 + (t558 - 0.1e1) * t567 * pkin(3);
t484 = pkin(1) * t566 + (-pkin(3) + t673 + 0.2e1 * t676) * t561;
t586 = -pkin(3) / 0.2e1;
t501 = t676 + t673 / 0.2e1 + t586;
t563 = legFrame(3,2);
t536 = sin(t563);
t539 = cos(t563);
t609 = t567 * t622;
t597 = pkin(2) * t609 + (t609 * t685 - t606) * t572;
t532 = pkin(1) * t561;
t619 = t572 * (-t566 * pkin(3) + t532);
t634 = t536 * t567;
t587 = pkin(2) / 0.2e1;
t643 = (t572 * pkin(3) + t587) * t566;
t457 = (-t501 * t634 + t539 * t643) * t683 + (t539 * t484 + t597 * t536) * t562 + t475 * t536 + t539 * t619;
t497 = 0.1e1 / (t562 * t572 - t622);
t667 = t457 * t497;
t628 = t539 * t567;
t458 = (t501 * t628 + t536 * t643) * t683 + (t536 * t484 - t597 * t539) * t562 - t475 * t539 + t536 * t619;
t666 = t458 * t497;
t569 = sin(qJ(1,2));
t549 = pkin(6) + t669;
t575 = cos(qJ(1,2));
t605 = pkin(1) * t569 - t575 * t549;
t568 = sin(qJ(3,2));
t621 = t561 * t568;
t476 = t605 * t621 + (t559 - 0.1e1) * t569 * pkin(3);
t485 = pkin(1) * t568 + (-pkin(3) + t672 + 0.2e1 * t675) * t561;
t502 = t675 + t672 / 0.2e1 + t586;
t564 = legFrame(2,2);
t537 = sin(t564);
t540 = cos(t564);
t608 = t569 * t621;
t596 = pkin(2) * t608 + (t608 * t685 - t605) * t574;
t618 = t574 * (-t568 * pkin(3) + t532);
t632 = t537 * t569;
t642 = (t574 * pkin(3) + t587) * t568;
t459 = (-t502 * t632 + t540 * t642) * t683 + (t540 * t485 + t596 * t537) * t562 + t476 * t537 + t540 * t618;
t498 = 0.1e1 / (t562 * t574 - t621);
t665 = t459 * t498;
t626 = t540 * t569;
t460 = (t502 * t626 + t537 * t642) * t683 + (t537 * t485 - t596 * t540) * t562 - t476 * t540 + t537 * t618;
t664 = t460 * t498;
t571 = sin(qJ(1,1));
t550 = pkin(6) + t670;
t577 = cos(qJ(1,1));
t604 = pkin(1) * t571 - t577 * t550;
t570 = sin(qJ(3,1));
t620 = t561 * t570;
t477 = t604 * t620 + (t560 - 0.1e1) * t571 * pkin(3);
t486 = pkin(1) * t570 + (-pkin(3) + t671 + 0.2e1 * t674) * t561;
t503 = t674 + t671 / 0.2e1 + t586;
t565 = legFrame(1,2);
t538 = sin(t565);
t541 = cos(t565);
t607 = t571 * t620;
t595 = pkin(2) * t607 + (t607 * t685 - t604) * t576;
t617 = t576 * (-t570 * pkin(3) + t532);
t630 = t538 * t571;
t641 = (t576 * pkin(3) + t587) * t570;
t461 = (-t503 * t630 + t541 * t641) * t683 + (t541 * t486 + t595 * t538) * t562 + t477 * t538 + t541 * t617;
t499 = 0.1e1 / (t562 * t576 - t620);
t663 = t461 * t499;
t624 = t541 * t571;
t462 = (t503 * t624 + t538 * t641) * t683 + (t538 * t486 - t595 * t541) * t562 - t477 * t541 + t538 * t617;
t662 = t462 * t499;
t555 = pkin(7) + qJ(3,3);
t526 = sin(t555);
t529 = cos(t555);
t601 = rSges(3,1) * t529 - rSges(3,2) * t526;
t508 = (m(2) * rSges(2,1) + pkin(2) * m(3)) * t562;
t522 = t561 * t681;
t582 = m(2) + m(3);
t602 = -pkin(1) * t582 - t508 + t522;
t472 = -t601 * m(3) + t602;
t516 = t562 * pkin(2) + pkin(1);
t481 = t567 * t548 + (pkin(3) * t529 + t516) * t573;
t533 = 0.1e1 / t548;
t469 = (t472 * t573 + t481 * t582) * t533;
t661 = t469 * t497;
t556 = pkin(7) + qJ(3,2);
t527 = sin(t556);
t530 = cos(t556);
t600 = rSges(3,1) * t530 - rSges(3,2) * t527;
t473 = -t600 * m(3) + t602;
t482 = t569 * t549 + (pkin(3) * t530 + t516) * t575;
t534 = 0.1e1 / t549;
t470 = (t473 * t575 + t482 * t582) * t534;
t660 = t470 * t498;
t557 = pkin(7) + qJ(3,1);
t528 = sin(t557);
t531 = cos(t557);
t599 = rSges(3,1) * t531 - rSges(3,2) * t528;
t474 = -t599 * m(3) + t602;
t483 = t571 * t550 + (pkin(3) * t531 + t516) * t577;
t535 = 0.1e1 / t550;
t471 = (t474 * t577 + t483 * t582) * t535;
t659 = t471 * t499;
t487 = t539 * t526 - t529 * t634;
t523 = 0.1e1 / t529;
t658 = t487 * t523;
t657 = t487 * t533;
t488 = t536 * t526 + t529 * t628;
t656 = t488 * t523;
t655 = t488 * t533;
t489 = t540 * t527 - t530 * t632;
t524 = 0.1e1 / t530;
t654 = t489 * t524;
t653 = t489 * t534;
t490 = t537 * t527 + t530 * t626;
t652 = t490 * t524;
t651 = t490 * t534;
t491 = t541 * t528 - t531 * t630;
t525 = 0.1e1 / t531;
t650 = t491 * t525;
t649 = t491 * t535;
t492 = t538 * t528 + t531 * t624;
t648 = t492 * t525;
t647 = t492 * t535;
t646 = t497 * t582;
t645 = t498 * t582;
t644 = t499 * t582;
t640 = t523 * t536;
t639 = t523 * t539;
t638 = t524 * t537;
t637 = t524 * t540;
t636 = t525 * t538;
t635 = t525 * t541;
t592 = 1 / pkin(3);
t633 = t536 * t592;
t631 = t537 * t592;
t629 = t538 * t592;
t627 = t539 * t592;
t625 = t540 * t592;
t623 = t541 * t592;
t615 = t472 * t497 * t533;
t614 = t473 * t498 * t534;
t613 = t474 * t499 * t535;
t478 = (-rSges(3,2) * t679 + Icges(3,6)) * t529 - (rSges(3,1) * t679 - Icges(3,5)) * t526;
t612 = t478 * t573 * t592;
t479 = (-rSges(3,2) * t678 + Icges(3,6)) * t530 - (rSges(3,1) * t678 - Icges(3,5)) * t527;
t611 = t479 * t575 * t592;
t480 = (-rSges(3,2) * t677 + Icges(3,6)) * t531 - (rSges(3,1) * t677 - Icges(3,5)) * t528;
t610 = t480 * t577 * t592;
t584 = 2 * pkin(1) ^ 2;
t589 = rSges(2,2) ^ 2;
t591 = rSges(2,1) ^ 2;
t603 = (2 * rSges(2,3) ^ 2) + t584 + t589 + t591;
t593 = pkin(2) ^ 2;
t598 = t584 / 0.2e1 + t588 / 0.2e1 + t590 / 0.2e1 + t593 / 0.2e1;
t585 = 0.2e1 * pkin(7);
t594 = Icges(1,3) + (m(3) * t593 + (-t589 + t591) * m(2) - Icges(2,1) + Icges(2,2)) * cos(t585) / 0.2e1 + (-rSges(2,1) * t681 + Icges(2,4)) * sin(t585) + t508 * t616 + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) - 0.2e1 * pkin(1) * t522 + t686 + Icges(2,2) / 0.2e1 + Icges(3,1) / 0.2e1 + Icges(2,1) / 0.2e1;
t554 = t585 + qJ(3,1);
t553 = t585 + qJ(3,2);
t552 = t585 + qJ(3,3);
t521 = -m(3) * rSges(3,1) * rSges(3,2) + Icges(3,4);
t519 = 0.2e1 * t557;
t518 = 0.2e1 * t556;
t517 = 0.2e1 * t555;
t511 = (t588 + t590) * m(3) + Icges(3,3);
t468 = (t480 * t647 + t511 * t629) * t525;
t467 = (t480 * t649 + t511 * t623) * t525;
t466 = (t479 * t651 + t511 * t631) * t524;
t465 = (t479 * t653 + t511 * t625) * t524;
t464 = (t478 * t655 + t511 * t633) * t523;
t463 = (t478 * t657 + t511 * t627) * t523;
t456 = (t547 ^ 2 + t599 * t616 + ((-sin(t554) - t570) * rSges(3,2) + (cos(t554) + t576) * rSges(3,1)) * pkin(2) + t598) * m(3) + t521 * sin(t519) + (((t684 + 2 * qJ(2,1)) * qJ(2,1)) + t603) * t682 + t594 + cos(t519) * t680;
t455 = (t546 ^ 2 + t600 * t616 + ((-sin(t553) - t568) * rSges(3,2) + (cos(t553) + t574) * rSges(3,1)) * pkin(2) + t598) * m(3) + t521 * sin(t518) + (((t684 + 2 * qJ(2,2)) * qJ(2,2)) + t603) * t682 + t594 + cos(t518) * t680;
t454 = (t545 ^ 2 + t601 * t616 + ((-sin(t552) - t566) * rSges(3,2) + (cos(t552) + t572) * rSges(3,1)) * pkin(2) + t598) * m(3) + t521 * sin(t517) + (((t684 + 2 * qJ(2,3)) * qJ(2,3)) + t603) * t682 + t594 + cos(t517) * t680;
t453 = (t456 * t577 + t474 * t483) * t535;
t452 = (t455 * t575 + t473 * t482) * t534;
t451 = (t454 * t573 + t472 * t481) * t533;
t450 = (t462 * t644 + t474 * t648) * t535;
t449 = (t461 * t644 + t474 * t650) * t535;
t448 = (t460 * t645 + t473 * t652) * t534;
t447 = (t459 * t645 + t473 * t654) * t534;
t446 = (t458 * t646 + t472 * t656) * t533;
t445 = (t457 * t646 + t472 * t658) * t533;
t444 = t462 * t613 + (t456 * t647 + t480 * t629) * t525;
t443 = t461 * t613 + (t456 * t649 + t480 * t623) * t525;
t442 = t460 * t614 + (t455 * t651 + t479 * t631) * t524;
t441 = t459 * t614 + (t455 * t653 + t479 * t625) * t524;
t440 = t458 * t615 + (t454 * t655 + t478 * t633) * t523;
t439 = t457 * t615 + (t454 * t657 + t478 * t627) * t523;
t1 = [m(4) + (t444 * t648 + t450 * t662) * t535 + (t442 * t652 + t448 * t664) * t534 + (t440 * t656 + t446 * t666) * t533 + (t464 * t640 + t466 * t638 + t468 * t636) * t592, (t444 * t650 + t450 * t663) * t535 + (t442 * t654 + t448 * t665) * t534 + (t440 * t658 + t446 * t667) * t533 + (t464 * t639 + t466 * t637 + t468 * t635) * t592, (t444 * t577 + t450 * t483) * t535 + (t442 * t575 + t448 * t482) * t534 + (t440 * t573 + t446 * t481) * t533; (t443 * t648 + t449 * t662) * t535 + (t441 * t652 + t447 * t664) * t534 + (t439 * t656 + t445 * t666) * t533 + (t463 * t640 + t465 * t638 + t467 * t636) * t592, m(4) + (t443 * t650 + t449 * t663) * t535 + (t441 * t654 + t447 * t665) * t534 + (t439 * t658 + t445 * t667) * t533 + (t463 * t639 + t465 * t637 + t467 * t635) * t592, (t443 * t577 + t449 * t483) * t535 + (t441 * t575 + t447 * t482) * t534 + (t439 * t573 + t445 * t481) * t533; (t462 * t659 + (t453 * t492 + t538 * t610) * t525) * t535 + (t460 * t660 + (t452 * t490 + t537 * t611) * t524) * t534 + (t458 * t661 + (t451 * t488 + t536 * t612) * t523) * t533, (t461 * t659 + (t453 * t491 + t541 * t610) * t525) * t535 + (t459 * t660 + (t452 * t489 + t540 * t611) * t524) * t534 + (t457 * t661 + (t451 * t487 + t539 * t612) * t523) * t533, m(4) + (t453 * t577 + t471 * t483) * t535 + (t452 * t575 + t470 * t482) * t534 + (t451 * t573 + t469 * t481) * t533;];
MX  = t1;
