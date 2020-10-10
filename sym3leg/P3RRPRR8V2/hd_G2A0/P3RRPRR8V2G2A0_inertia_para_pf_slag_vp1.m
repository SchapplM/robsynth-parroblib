% Calculate inertia matrix for parallel robot
% P3RRPRR8V2G2A0
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2020-08-06 21:13
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRPRR8V2G2A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:11:43
% EndTime: 2020-08-06 21:11:45
% DurationCPUTime: 1.77s
% Computational Cost: add. (5187->302), mult. (8661->475), div. (348->9), fcn. (4200->44), ass. (0->227)
t701 = rSges(2,3) + pkin(5);
t700 = 2 * pkin(1);
t699 = m(2) / 0.2e1;
t698 = Icges(2,2) / 0.2e1;
t697 = Icges(3,2) / 0.2e1;
t611 = pkin(2) ^ 2;
t686 = t611 / 0.2e1;
t696 = -2 * pkin(1);
t592 = cos(qJ(2,3));
t694 = 0.2e1 * t592 ^ 2;
t594 = cos(qJ(2,2));
t693 = 0.2e1 * t594 ^ 2;
t596 = cos(qJ(2,1));
t692 = 0.2e1 * t596 ^ 2;
t691 = m(2) * rSges(2,1);
t690 = m(2) * rSges(2,2);
t689 = m(3) * rSges(3,1);
t607 = rSges(2,2) ^ 2;
t609 = rSges(2,1) ^ 2;
t688 = m(3) * t686 + (-t607 + t609) * t699 + t698 - Icges(2,1) / 0.2e1;
t606 = rSges(3,2) ^ 2;
t608 = rSges(3,1) ^ 2;
t687 = (-t606 + t608) * m(3) / 0.2e1 - Icges(3,1) / 0.2e1 + t697;
t573 = qJ(2,3) + pkin(7);
t541 = sin(t573);
t544 = cos(t573);
t556 = t592 * pkin(2);
t509 = -t544 * rSges(3,1) + t541 * rSges(3,2) - pkin(1) - t556;
t685 = m(3) * t509;
t574 = qJ(2,2) + pkin(7);
t542 = sin(t574);
t545 = cos(t574);
t557 = t594 * pkin(2);
t510 = -t545 * rSges(3,1) + t542 * rSges(3,2) - pkin(1) - t557;
t684 = m(3) * t510;
t576 = qJ(2,1) + pkin(7);
t543 = sin(t576);
t546 = cos(t576);
t558 = t596 * pkin(2);
t511 = -t546 * rSges(3,1) + t543 * rSges(3,2) - pkin(1) - t558;
t683 = m(3) * t511;
t672 = pkin(5) + qJ(3,3);
t563 = rSges(3,3) + t672;
t682 = m(3) * t563;
t673 = pkin(5) + qJ(3,2);
t564 = rSges(3,3) + t673;
t681 = m(3) * t564;
t674 = pkin(5) + qJ(3,1);
t565 = rSges(3,3) + t674;
t680 = m(3) * t565;
t587 = sin(qJ(1,3));
t566 = pkin(6) + t672;
t593 = cos(qJ(1,3));
t618 = pkin(1) * t587 - t593 * t566;
t582 = cos(pkin(7));
t570 = t582 ^ 2;
t632 = pkin(3) * (t570 - 0.1e1);
t581 = sin(pkin(7));
t586 = sin(qJ(2,3));
t650 = t581 * t586;
t679 = pkin(3) * (t587 * t632 + t618 * t650);
t589 = sin(qJ(1,2));
t567 = pkin(6) + t673;
t595 = cos(qJ(1,2));
t617 = pkin(1) * t589 - t595 * t567;
t588 = sin(qJ(2,2));
t649 = t581 * t588;
t678 = pkin(3) * (t589 * t632 + t617 * t649);
t591 = sin(qJ(1,1));
t568 = pkin(6) + t674;
t597 = cos(qJ(1,1));
t616 = pkin(1) * t591 - t597 * t568;
t590 = sin(qJ(2,1));
t648 = t581 * t590;
t677 = pkin(3) * (t591 * t632 + t616 * t648);
t676 = t581 * pkin(3);
t675 = t582 * pkin(3);
t533 = pkin(2) + t675;
t635 = pkin(3) * t650;
t512 = 0.1e1 / (t533 * t592 - t635);
t547 = 0.1e1 / t566;
t671 = t512 * t547;
t634 = pkin(3) * t649;
t513 = 0.1e1 / (t533 * t594 - t634);
t548 = 0.1e1 / t567;
t670 = t513 * t548;
t633 = pkin(3) * t648;
t514 = 0.1e1 / (t533 * t596 - t633);
t549 = 0.1e1 / t568;
t669 = t514 * t549;
t621 = pkin(3) * t544 + t556;
t517 = 0.1e1 / t621;
t583 = legFrame(3,2);
t550 = sin(t583);
t668 = t517 * t550;
t553 = cos(t583);
t667 = t517 * t553;
t620 = pkin(3) * t545 + t557;
t518 = 0.1e1 / t620;
t584 = legFrame(2,2);
t551 = sin(t584);
t666 = t518 * t551;
t554 = cos(t584);
t665 = t518 * t554;
t619 = pkin(3) * t546 + t558;
t519 = 0.1e1 / t619;
t585 = legFrame(1,2);
t552 = sin(t585);
t664 = t519 * t552;
t555 = cos(t585);
t663 = t519 * t555;
t662 = t533 * t553;
t661 = t533 * t554;
t660 = t533 * t555;
t659 = t550 * t533;
t658 = t550 * t587;
t657 = t551 * t533;
t656 = t551 * t589;
t655 = t552 * t533;
t654 = t552 * t591;
t653 = t553 * t587;
t652 = t554 * t589;
t651 = t555 * t591;
t647 = t607 + t609;
t645 = pkin(2) * t675;
t644 = m(3) * t671;
t643 = m(3) * t670;
t642 = m(3) * t669;
t641 = t553 * t676;
t640 = t554 * t676;
t639 = t555 * t676;
t638 = t550 * t676;
t637 = t551 * t676;
t636 = t552 * t676;
t526 = rSges(3,2) * t682 - Icges(3,6);
t529 = rSges(3,1) * t682 - Icges(3,5);
t615 = -t701 * t691 + Icges(2,5);
t622 = t701 * t690 - Icges(2,6);
t484 = (-pkin(2) * t682 + t526 * t581 - t529 * t582 + t615) * t586 - (t526 * t582 + t529 * t581 + t622) * t592;
t631 = t484 * t671;
t527 = rSges(3,2) * t681 - Icges(3,6);
t530 = rSges(3,1) * t681 - Icges(3,5);
t485 = (-pkin(2) * t681 + t527 * t581 - t530 * t582 + t615) * t588 - (t527 * t582 + t530 * t581 + t622) * t594;
t630 = t485 * t670;
t528 = rSges(3,2) * t680 - Icges(3,6);
t531 = rSges(3,1) * t680 - Icges(3,5);
t486 = (-pkin(2) * t680 + t528 * t581 - t531 * t582 + t615) * t590 - (t528 * t582 + t531 * t581 + t622) * t596;
t629 = t486 * t669;
t628 = t484 * t668;
t627 = t485 * t666;
t626 = t486 * t664;
t625 = t484 * t667;
t624 = t485 * t665;
t623 = t486 * t663;
t602 = 2 * pkin(1) ^ 2;
t614 = t602 / 0.2e1 + t606 / 0.2e1 + t608 / 0.2e1 + t686;
t610 = pkin(3) ^ 2;
t613 = 0.2e1 * t570 * t610 - t610 + t611 + 0.2e1 * t645;
t532 = t582 * pkin(2) * t689;
t612 = Icges(1,3) + ((2 * pkin(5) ^ 2) + t602 + ((4 * pkin(5) + 2 * rSges(2,3)) * rSges(2,3)) + t647) * t699 + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + t532 + t697 + t698 + Icges(3,1) / 0.2e1 + Icges(2,1) / 0.2e1;
t605 = 0.2e1 * qJ(2,1);
t604 = 0.2e1 * qJ(2,2);
t603 = 0.2e1 * qJ(2,3);
t575 = t605 + pkin(7);
t572 = t603 + pkin(7);
t571 = pkin(7) + t604;
t540 = pkin(1) * t676;
t539 = -rSges(2,1) * t690 + Icges(2,4);
t538 = -rSges(3,2) * t689 + Icges(3,4);
t537 = pkin(2) * m(3) + t691;
t536 = 0.2e1 * t576;
t535 = 0.2e1 * t574;
t534 = 0.2e1 * t573;
t523 = pkin(1) * t590 - t676;
t522 = pkin(1) * t588 - t676;
t521 = pkin(1) * t586 - t676;
t515 = t645 + t686 + (t570 - 0.1e1 / 0.2e1) * t610;
t508 = -0.2e1 * t591 * t633 + t616;
t507 = -0.2e1 * t589 * t634 + t617;
t506 = -0.2e1 * t587 * t635 + t618;
t505 = t591 * t568 + (pkin(1) + t619) * t597;
t504 = t589 * t567 + (pkin(1) + t620) * t595;
t503 = t587 * t566 + (pkin(1) + t621) * t593;
t502 = t590 * t613 + t540;
t501 = t588 * t613 + t540;
t500 = t586 * t613 + t540;
t499 = 0.2e1 * t532 + t647 * m(2) + Icges(2,3) + Icges(3,3) + (-0.2e1 * t581 * rSges(3,2) * pkin(2) + t606 + t608 + t611) * m(3);
t495 = (t533 * t651 + t636) * t596 + t590 * (-t591 * t639 + t655);
t494 = (t533 * t652 + t637) * t594 + t588 * (-t589 * t640 + t657);
t493 = (t533 * t653 + t638) * t592 + t586 * (-t587 * t641 + t659);
t492 = (-t533 * t654 + t639) * t596 + (t591 * t636 + t660) * t590;
t491 = (-t533 * t656 + t640) * t594 + (t589 * t637 + t661) * t588;
t490 = (-t533 * t658 + t641) * t592 + (t587 * t638 + t662) * t586;
t489 = (t511 * t597 + t505) * t549 * m(3);
t488 = (t510 * t595 + t504) * t548 * m(3);
t487 = (t509 * t593 + t503) * t547 * m(3);
t483 = (t515 * t651 + t533 * t636) * t692 + (t552 * t502 + t508 * t660) * t596 - t555 * t677 + t523 * t655;
t482 = (-t515 * t654 + t533 * t639) * t692 + (t502 * t555 - t508 * t655) * t596 + t552 * t677 + t523 * t660;
t481 = (t515 * t652 + t533 * t637) * t693 + (t551 * t501 + t507 * t661) * t594 - t554 * t678 + t522 * t657;
t480 = (-t515 * t656 + t533 * t640) * t693 + (t501 * t554 - t507 * t657) * t594 + t551 * t678 + t522 * t661;
t479 = (t515 * t653 + t533 * t638) * t694 + (t550 * t500 + t506 * t662) * t592 - t553 * t679 + t521 * t659;
t478 = (-t515 * t658 + t533 * t641) * t694 + (t500 * t553 - t506 * t659) * t592 + t550 * t679 + t521 * t662;
t477 = cos(t536) * t687 + cos(t605) * t688 + (t537 * t596 - t590 * t690) * t700 + t538 * sin(t536) + t539 * sin(t605) + (t565 ^ 2 + (pkin(2) * cos(t575) + t546 * t700) * rSges(3,1) + (t543 * t696 + (-sin(t575) - t581) * pkin(2)) * rSges(3,2) + t614) * m(3) + t612;
t476 = cos(t535) * t687 + cos(t604) * t688 + (t537 * t594 - t588 * t690) * t700 + t538 * sin(t535) + t539 * sin(t604) + (t564 ^ 2 + (pkin(2) * cos(t571) + t545 * t700) * rSges(3,1) + (t542 * t696 + (-sin(t571) - t581) * pkin(2)) * rSges(3,2) + t614) * m(3) + t612;
t475 = cos(t534) * t687 + cos(t603) * t688 + (t537 * t592 - t586 * t690) * t700 + t538 * sin(t534) + t539 * sin(t603) + (t563 ^ 2 + (pkin(2) * cos(t572) + t544 * t700) * rSges(3,1) + (t541 * t696 + (-sin(t572) - t581) * pkin(2)) * rSges(3,2) + t614) * m(3) + t612;
t474 = t495 * t629 + t499 * t664;
t473 = t494 * t630 + t499 * t666;
t472 = t493 * t631 + t499 * t668;
t471 = t492 * t629 + t499 * t663;
t470 = t491 * t630 + t499 * t665;
t469 = t490 * t631 + t499 * t667;
t468 = (t477 * t597 + t505 * t683) * t549;
t467 = (t476 * t595 + t504 * t684) * t548;
t466 = (t475 * t593 + t503 * t685) * t547;
t465 = (t495 * t511 + t483) * t642;
t464 = (t494 * t510 + t481) * t643;
t463 = (t493 * t509 + t479) * t644;
t462 = (t492 * t511 + t482) * t642;
t461 = (t491 * t510 + t480) * t643;
t460 = (t490 * t509 + t478) * t644;
t459 = t626 + (t477 * t495 + t483 * t683) * t669;
t458 = t627 + (t476 * t494 + t481 * t684) * t670;
t457 = t628 + (t475 * t493 + t479 * t685) * t671;
t456 = t623 + (t477 * t492 + t482 * t683) * t669;
t455 = t624 + (t476 * t491 + t480 * t684) * t670;
t454 = t625 + (t475 * t490 + t478 * t685) * t671;
t1 = [t472 * t668 + t473 * t666 + t474 * t664 + m(4) + (t459 * t495 + t465 * t483) * t669 + (t458 * t494 + t464 * t481) * t670 + (t457 * t493 + t463 * t479) * t671, t472 * t667 + t473 * t665 + t474 * t663 + (t459 * t492 + t465 * t482) * t669 + (t458 * t491 + t464 * t480) * t670 + (t457 * t490 + t463 * t478) * t671, (t459 * t597 + t465 * t505) * t549 + (t458 * t595 + t464 * t504) * t548 + (t457 * t593 + t463 * t503) * t547; t469 * t668 + t470 * t666 + t471 * t664 + (t456 * t495 + t462 * t483) * t669 + (t455 * t494 + t461 * t481) * t670 + (t454 * t493 + t460 * t479) * t671, t469 * t667 + t470 * t665 + t471 * t663 + m(4) + (t456 * t492 + t462 * t482) * t669 + (t455 * t491 + t461 * t480) * t670 + (t454 * t490 + t460 * t478) * t671, (t456 * t597 + t462 * t505) * t549 + (t455 * t595 + t461 * t504) * t548 + (t454 * t593 + t460 * t503) * t547; (t597 * t626 + (t468 * t495 + t483 * t489) * t514) * t549 + (t595 * t627 + (t467 * t494 + t481 * t488) * t513) * t548 + (t593 * t628 + (t466 * t493 + t479 * t487) * t512) * t547, (t597 * t623 + (t468 * t492 + t482 * t489) * t514) * t549 + (t595 * t624 + (t467 * t491 + t480 * t488) * t513) * t548 + (t593 * t625 + (t466 * t490 + t478 * t487) * t512) * t547, m(4) + (t468 * t597 + t489 * t505) * t549 + (t467 * t595 + t488 * t504) * t548 + (t466 * t593 + t487 * t503) * t547;];
MX  = t1;
