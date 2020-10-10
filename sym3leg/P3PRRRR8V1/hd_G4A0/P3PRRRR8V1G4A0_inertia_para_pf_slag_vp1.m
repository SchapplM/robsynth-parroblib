% Calculate inertia matrix for parallel robot
% P3PRRRR8V1G4A0
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
% Datum: 2020-08-06 17:27
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRRR8V1G4A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G4A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G4A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G4A0_inertia_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G4A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V1G4A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR8V1G4A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G4A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G4A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:24:48
% EndTime: 2020-08-06 17:24:50
% DurationCPUTime: 2.23s
% Computational Cost: add. (6132->280), mult. (15858->594), div. (648->7), fcn. (16866->40), ass. (0->258)
t715 = m(3) / 0.2e1;
t714 = Icges(3,2) / 0.2e1;
t603 = sin(qJ(2,3));
t609 = cos(qJ(2,3));
t608 = cos(qJ(3,3));
t702 = pkin(2) * t608;
t547 = -pkin(5) * t609 + t603 * t702;
t590 = sin(pkin(3));
t592 = cos(pkin(3));
t602 = sin(qJ(3,3));
t705 = pkin(2) * t602;
t529 = t547 * t590 + t592 * t705;
t713 = 0.1e1 / t529;
t605 = sin(qJ(2,2));
t611 = cos(qJ(2,2));
t610 = cos(qJ(3,2));
t701 = pkin(2) * t610;
t548 = -pkin(5) * t611 + t605 * t701;
t604 = sin(qJ(3,2));
t704 = pkin(2) * t604;
t530 = t548 * t590 + t592 * t704;
t712 = 0.1e1 / t530;
t607 = sin(qJ(2,1));
t613 = cos(qJ(2,1));
t612 = cos(qJ(3,1));
t700 = pkin(2) * t612;
t549 = -pkin(5) * t613 + t607 * t700;
t606 = sin(qJ(3,1));
t703 = pkin(2) * t606;
t531 = t549 * t590 + t592 * t703;
t711 = 0.1e1 / t531;
t710 = m(3) * rSges(3,3);
t620 = rSges(3,2) ^ 2;
t621 = rSges(3,1) ^ 2;
t709 = (-t620 + t621) * t715 - Icges(3,1) / 0.2e1 + t714;
t635 = rSges(3,1) * t608 - rSges(3,2) * t602;
t708 = m(3) * (t635 * t592 - t590 * t603 * (rSges(3,1) * t602 + rSges(3,2) * t608));
t634 = rSges(3,1) * t610 - rSges(3,2) * t604;
t707 = m(3) * (t634 * t592 - t590 * t605 * (rSges(3,1) * t604 + rSges(3,2) * t610));
t633 = rSges(3,1) * t612 - rSges(3,2) * t606;
t706 = m(3) * (t633 * t592 - t590 * t607 * (rSges(3,1) * t606 + rSges(3,2) * t612));
t593 = legFrame(3,3);
t567 = sin(t593);
t573 = cos(t593);
t589 = sin(pkin(6));
t591 = cos(pkin(6));
t535 = t567 * t591 + t573 * t589;
t538 = -t567 * t589 + t573 * t591;
t596 = legFrame(3,1);
t570 = sin(t596);
t576 = cos(t596);
t599 = legFrame(3,2);
t579 = sin(t599);
t678 = t576 * t579;
t494 = -t535 * t570 + t538 * t678;
t669 = t592 * t603;
t541 = t589 * t609 + t591 * t669;
t544 = -t589 * t669 + t591 * t609;
t514 = t541 * t573 + t544 * t567;
t628 = t541 * t567 - t544 * t573;
t672 = t590 * t608;
t472 = (t514 * t678 - t628 * t570) * t602 + t494 * t672;
t586 = 0.1e1 / t608;
t699 = t472 * t586;
t594 = legFrame(2,3);
t568 = sin(t594);
t574 = cos(t594);
t536 = t568 * t591 + t574 * t589;
t539 = -t568 * t589 + t574 * t591;
t597 = legFrame(2,1);
t571 = sin(t597);
t577 = cos(t597);
t600 = legFrame(2,2);
t580 = sin(t600);
t677 = t577 * t580;
t496 = -t536 * t571 + t539 * t677;
t668 = t592 * t605;
t542 = t589 * t611 + t591 * t668;
t545 = -t589 * t668 + t591 * t611;
t515 = t542 * t574 + t545 * t568;
t627 = t542 * t568 - t545 * t574;
t671 = t590 * t610;
t473 = (t515 * t677 - t627 * t571) * t604 + t496 * t671;
t587 = 0.1e1 / t610;
t698 = t473 * t587;
t595 = legFrame(1,3);
t569 = sin(t595);
t575 = cos(t595);
t537 = t569 * t591 + t575 * t589;
t540 = -t569 * t589 + t575 * t591;
t598 = legFrame(1,1);
t572 = sin(t598);
t578 = cos(t598);
t601 = legFrame(1,2);
t581 = sin(t601);
t676 = t578 * t581;
t498 = -t537 * t572 + t540 * t676;
t667 = t592 * t607;
t543 = t589 * t613 + t591 * t667;
t546 = -t589 * t667 + t591 * t613;
t516 = t543 * t575 + t546 * t569;
t626 = t543 * t569 - t546 * t575;
t670 = t590 * t612;
t474 = (t516 * t676 - t626 * t572) * t606 + t498 * t670;
t588 = 0.1e1 / t612;
t697 = t474 * t588;
t681 = t570 * t579;
t499 = t535 * t576 + t538 * t681;
t475 = (-t514 * t681 - t628 * t576) * t602 - t499 * t672;
t696 = t475 * t586;
t680 = t571 * t580;
t500 = t536 * t577 + t539 * t680;
t476 = (-t515 * t680 - t627 * t577) * t604 - t500 * t671;
t695 = t476 * t587;
t679 = t572 * t581;
t501 = t537 * t578 + t540 * t679;
t477 = (-t516 * t679 - t626 * t578) * t606 - t501 * t670;
t694 = t477 * t588;
t566 = -m(3) * rSges(3,1) * rSges(3,2) + Icges(3,4);
t617 = 0.2e1 * qJ(3,3);
t663 = t620 + t621;
t629 = Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + (0.2e1 * rSges(3,3) ^ 2 + t663) * t715 + t714 + Icges(3,1) / 0.2e1;
t505 = cos(t617) * t709 + t566 * sin(t617) + t629;
t693 = t505 * t586;
t618 = 0.2e1 * qJ(3,2);
t506 = cos(t618) * t709 + t566 * sin(t618) + t629;
t692 = t506 * t587;
t619 = 0.2e1 * qJ(3,1);
t507 = cos(t619) * t709 + t566 * sin(t619) + t629;
t691 = t507 * t588;
t563 = m(2) * rSges(2,2) - t710;
t616 = m(2) * rSges(2,1);
t690 = ((t635 * m(3) + t616) * t609 - t603 * t563) * t590;
t689 = ((t634 * m(3) + t616) * t611 - t605 * t563) * t590;
t688 = ((t633 * m(3) + t616) * t613 - t607 * t563) * t590;
t687 = t713 * t586;
t686 = t712 * t587;
t685 = t711 * t588;
t564 = -rSges(3,2) * t710 + Icges(3,6);
t565 = rSges(3,1) * t710 - Icges(3,5);
t532 = t564 * t608 - t565 * t602;
t684 = t532 * t586;
t533 = t564 * t610 - t565 * t604;
t683 = t533 * t587;
t534 = t564 * t612 - t565 * t606;
t682 = t534 * t588;
t582 = cos(t599);
t675 = t582 * t586;
t583 = cos(t600);
t674 = t583 * t587;
t584 = cos(t601);
t673 = t584 * t588;
t666 = t592 * t609;
t665 = t592 * t611;
t664 = t592 * t613;
t493 = t535 * t678 + t538 * t570;
t466 = (-t493 * t603 + t494 * t666) * t702 + pkin(5) * (t609 * t493 + t494 * t669);
t662 = t466 * t687;
t502 = t535 * t681 - t538 * t576;
t467 = -(t499 * t666 - t502 * t603) * t702 - pkin(5) * (t499 * t669 + t502 * t609);
t661 = t467 * t687;
t495 = t536 * t677 + t539 * t571;
t468 = (-t495 * t605 + t496 * t665) * t701 + pkin(5) * (t495 * t611 + t496 * t668);
t660 = t468 * t686;
t503 = t536 * t680 - t539 * t577;
t469 = -(t500 * t665 - t503 * t605) * t701 - pkin(5) * (t500 * t668 + t503 * t611);
t659 = t469 * t686;
t497 = t537 * t676 + t540 * t572;
t470 = (-t497 * t607 + t498 * t664) * t700 + pkin(5) * (t497 * t613 + t498 * t667);
t658 = t470 * t685;
t504 = t537 * t679 - t540 * t578;
t471 = -(t501 * t664 - t504 * t607) * t700 - pkin(5) * (t501 * t667 + t504 * t613);
t657 = t471 * t685;
t490 = t514 * t602 + t538 * t672;
t656 = t490 * t675;
t491 = t515 * t604 + t539 * t671;
t655 = t491 * t674;
t492 = t516 * t606 + t540 * t670;
t654 = t492 * t673;
t653 = t586 * t690;
t652 = t587 * t689;
t651 = t588 * t688;
t622 = 0.1e1 / pkin(2);
t650 = t622 * t687;
t649 = t622 * t686;
t648 = t622 * t685;
t647 = ((-t535 * t603 + t538 * t666) * t702 + pkin(5) * (t535 * t609 + t538 * t669)) * t713 * t675;
t646 = ((-t536 * t605 + t539 * t665) * t701 + pkin(5) * (t536 * t611 + t539 * t668)) * t712 * t674;
t645 = ((-t537 * t607 + t540 * t664) * t700 + pkin(5) * (t537 * t613 + t540 * t667)) * t711 * t673;
t644 = t532 * t650;
t559 = t663 * m(3) + Icges(3,3);
t643 = t559 * t650;
t642 = t533 * t649;
t641 = t559 * t649;
t640 = t534 * t648;
t639 = t559 * t648;
t638 = t650 * t708;
t637 = t649 * t707;
t636 = t648 * t706;
t632 = t622 * t647;
t631 = t622 * t646;
t630 = t622 * t645;
t625 = -t547 * t592 + t590 * t705;
t624 = -t548 * t592 + t590 * t704;
t623 = -t549 * t592 + t590 * t703;
t585 = m(1) + m(2) + m(3);
t552 = pkin(5) * t607 + t613 * t700;
t551 = pkin(5) * t605 + t611 * t701;
t550 = pkin(5) * t603 + t609 * t702;
t513 = -t552 * t589 + t623 * t591;
t512 = -t551 * t589 + t624 * t591;
t511 = -t550 * t589 + t625 * t591;
t510 = t552 * t591 + t623 * t589;
t509 = t551 * t591 + t624 * t589;
t508 = t550 * t591 + t625 * t589;
t486 = -t510 * t569 + t513 * t575;
t485 = -t509 * t568 + t512 * t574;
t484 = -t508 * t567 + t511 * t573;
t483 = (t623 * t537 + t540 * t552) * t584 + t581 * t531;
t482 = (t624 * t536 + t539 * t551) * t583 + t580 * t530;
t481 = (t625 * t535 + t538 * t550) * t582 + t579 * t529;
t480 = -t531 * t584 + (t510 * t575 + t513 * t569) * t581;
t479 = -t530 * t583 + (t509 * t574 + t512 * t568) * t580;
t478 = -t529 * t582 + (t508 * t573 + t511 * t567) * t579;
t465 = -t480 * t578 - t486 * t572;
t464 = t480 * t572 - t486 * t578;
t463 = -t479 * t577 - t485 * t571;
t462 = t479 * t571 - t485 * t577;
t461 = -t478 * t576 - t484 * t570;
t460 = t478 * t570 - t484 * t576;
t459 = -t559 * t630 + (t483 * t706 - t534 * t654) * t711;
t458 = -t559 * t631 + (t482 * t707 - t533 * t655) * t712;
t457 = -t559 * t632 + (t481 * t708 - t532 * t656) * t713;
t456 = -t630 * t706 + (-t492 * t584 * t651 + t483 * t585) * t711;
t455 = -t631 * t707 + (-t491 * t583 * t652 + t482 * t585) * t712;
t454 = -t632 * t708 + (-t490 * t582 * t653 + t481 * t585) * t713;
t453 = -t534 * t630 + (t483 * t688 - t507 * t654) * t711;
t452 = -t533 * t631 + (t482 * t689 - t506 * t655) * t712;
t451 = -t532 * t632 + (t481 * t690 - t505 * t656) * t713;
t450 = t470 * t639 + (t465 * t706 + t474 * t682) * t711;
t449 = t471 * t639 + (t464 * t706 + t477 * t682) * t711;
t448 = t468 * t641 + (t463 * t707 + t473 * t683) * t712;
t447 = t469 * t641 + (t462 * t707 + t476 * t683) * t712;
t446 = t466 * t643 + (t461 * t708 + t472 * t684) * t713;
t445 = t467 * t643 + (t460 * t708 + t475 * t684) * t713;
t444 = t470 * t636 + (t465 * t585 + t474 * t651) * t711;
t443 = t471 * t636 + (t464 * t585 + t477 * t651) * t711;
t442 = t468 * t637 + (t463 * t585 + t473 * t652) * t712;
t441 = t469 * t637 + (t462 * t585 + t476 * t652) * t712;
t440 = t466 * t638 + (t461 * t585 + t472 * t653) * t713;
t439 = t467 * t638 + (t460 * t585 + t475 * t653) * t713;
t438 = t470 * t640 + (t465 * t688 + t474 * t691) * t711;
t437 = t471 * t640 + (t464 * t688 + t477 * t691) * t711;
t436 = t468 * t642 + (t463 * t689 + t473 * t692) * t712;
t435 = t469 * t642 + (t462 * t689 + t476 * t692) * t712;
t434 = t466 * t644 + (t461 * t690 + t472 * t693) * t713;
t433 = t467 * t644 + (t460 * t690 + t475 * t693) * t713;
t1 = [m(4) + (-t453 * t654 + t456 * t483) * t711 + (-t452 * t655 + t455 * t482) * t712 + (-t451 * t656 + t454 * t481) * t713 + (-t457 * t647 - t458 * t646 - t459 * t645) * t622, (t453 * t694 + t456 * t464) * t711 + (t452 * t695 + t455 * t462) * t712 + (t451 * t696 + t454 * t460) * t713 + (t457 * t661 + t458 * t659 + t459 * t657) * t622, (t453 * t697 + t456 * t465) * t711 + (t452 * t698 + t455 * t463) * t712 + (t451 * t699 + t454 * t461) * t713 + (t457 * t662 + t458 * t660 + t459 * t658) * t622; (-t437 * t654 + t443 * t483) * t711 + (-t435 * t655 + t441 * t482) * t712 + (-t433 * t656 + t439 * t481) * t713 + (-t445 * t647 - t447 * t646 - t449 * t645) * t622, m(4) + (t437 * t694 + t443 * t464) * t711 + (t435 * t695 + t441 * t462) * t712 + (t433 * t696 + t439 * t460) * t713 + (t445 * t661 + t447 * t659 + t449 * t657) * t622, (t437 * t697 + t443 * t465) * t711 + (t435 * t698 + t441 * t463) * t712 + (t433 * t699 + t439 * t461) * t713 + (t445 * t662 + t447 * t660 + t449 * t658) * t622; (-t438 * t654 + t444 * t483) * t711 + (-t436 * t655 + t442 * t482) * t712 + (-t434 * t656 + t440 * t481) * t713 + (-t446 * t647 - t448 * t646 - t450 * t645) * t622, (t438 * t694 + t444 * t464) * t711 + (t436 * t695 + t442 * t462) * t712 + (t434 * t696 + t440 * t460) * t713 + (t446 * t661 + t448 * t659 + t450 * t657) * t622, m(4) + (t438 * t697 + t444 * t465) * t711 + (t436 * t698 + t442 * t463) * t712 + (t434 * t699 + t440 * t461) * t713 + (t446 * t662 + t448 * t660 + t450 * t658) * t622;];
MX  = t1;
