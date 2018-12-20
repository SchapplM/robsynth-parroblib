% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRP2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
% m [4x1]
%   mass of all robot links (including platform)
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
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2018-12-20 17:39
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taucX = P3PRP2A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRP2A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRP2A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRP2A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PRP2A0_coriolisvec_para_pf_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRP2A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRP2A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRP2A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRP2A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRP2A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:38:51
% EndTime: 2018-12-20 17:38:55
% DurationCPUTime: 4.04s
% Computational Cost: add. (26341->375), mult. (50050->590), div. (2250->3), fcn. (23653->14), ass. (0->258)
t618 = (qJ(3,1) ^ 2);
t628 = (pkin(2) ^ 2);
t720 = -t628 - 1;
t579 = -t618 - t720;
t607 = cos(qJ(2,1));
t593 = t607 ^ 2;
t604 = sin(qJ(2,1));
t672 = t604 * t607;
t659 = pkin(2) * t672;
t726 = 2 * qJ(3,1);
t548 = 0.1e1 / (t579 * t593 + t659 * t726 - t618 + t720);
t601 = legFrame(1,3);
t585 = cos(t601);
t582 = sin(t601);
t707 = qJ(3,1) * t585;
t664 = pkin(2) * t707;
t635 = t582 * t618 - t664;
t708 = qJ(3,1) * t582;
t573 = pkin(2) * t708;
t669 = t618 * t585 + t573;
t533 = t635 * t607 - t604 * (t585 + t669);
t613 = xP(3);
t587 = sin(t613);
t588 = cos(t613);
t621 = koppelP(1,2);
t624 = koppelP(1,1);
t566 = -t587 * t621 + t588 * t624;
t610 = xDP(3);
t611 = xDP(2);
t551 = t566 * t610 + t611;
t691 = t533 * t551;
t530 = t669 * t607 - t604 * (-t582 - t635);
t563 = t587 * t624 + t588 * t621;
t612 = xDP(1);
t554 = -t563 * t610 + t612;
t694 = t530 * t554;
t730 = (t691 + t694) * t548;
t617 = (qJ(3,2) ^ 2);
t578 = -t617 - t720;
t606 = cos(qJ(2,2));
t592 = t606 ^ 2;
t603 = sin(qJ(2,2));
t673 = t603 * t606;
t660 = pkin(2) * t673;
t725 = 2 * qJ(3,2);
t547 = 0.1e1 / (t578 * t592 + t660 * t725 - t617 + t720);
t600 = legFrame(2,3);
t584 = cos(t600);
t581 = sin(t600);
t704 = qJ(3,2) * t584;
t663 = pkin(2) * t704;
t634 = t581 * t617 - t663;
t705 = qJ(3,2) * t581;
t572 = pkin(2) * t705;
t670 = t617 * t584 + t572;
t532 = t634 * t606 - t603 * (t584 + t670);
t620 = koppelP(2,2);
t623 = koppelP(2,1);
t565 = -t587 * t620 + t588 * t623;
t550 = t565 * t610 + t611;
t692 = t532 * t550;
t529 = t670 * t606 - t603 * (-t581 - t634);
t562 = t587 * t623 + t588 * t620;
t553 = -t562 * t610 + t612;
t695 = t529 * t553;
t729 = (t692 + t695) * t547;
t616 = (qJ(3,3) ^ 2);
t577 = -t616 - t720;
t605 = cos(qJ(2,3));
t591 = t605 ^ 2;
t602 = sin(qJ(2,3));
t674 = t602 * t605;
t661 = pkin(2) * t674;
t724 = 2 * qJ(3,3);
t546 = 0.1e1 / (t577 * t591 + t661 * t724 - t616 + t720);
t599 = legFrame(3,3);
t583 = cos(t599);
t580 = sin(t599);
t701 = qJ(3,3) * t583;
t662 = pkin(2) * t701;
t633 = t580 * t616 - t662;
t702 = qJ(3,3) * t580;
t571 = pkin(2) * t702;
t671 = t616 * t583 + t571;
t531 = t633 * t605 - t602 * (t583 + t671);
t619 = koppelP(3,2);
t622 = koppelP(3,1);
t564 = -t587 * t619 + t588 * t622;
t549 = t564 * t610 + t611;
t693 = t531 * t549;
t528 = t671 * t605 - t602 * (-t580 - t633);
t561 = t587 * t622 + t588 * t619;
t552 = -t561 * t610 + t612;
t696 = t528 * t552;
t728 = (t693 + t696) * t546;
t727 = -2 * m(3);
t594 = t610 ^ 2;
t723 = 0.2e1 * t728;
t722 = 0.2e1 * t729;
t721 = 0.2e1 * t730;
t719 = m(3) * t546;
t718 = m(3) * t547;
t717 = m(3) * t548;
t716 = m(3) * t605;
t715 = m(3) * t606;
t714 = m(3) * t607;
t713 = m(4) * t594;
t712 = (rSges(3,3) + qJ(3,3)) * m(3);
t711 = (rSges(3,3) + qJ(3,2)) * m(3);
t710 = (rSges(3,3) + qJ(3,1)) * m(3);
t608 = pkin(2) + rSges(3,1);
t709 = t608 * m(3);
t706 = qJ(3,1) * t607;
t703 = qJ(3,2) * t606;
t700 = qJ(3,3) * t605;
t534 = 0.2e1 * t580 * t700 - t602 * (pkin(2) * t580 + t701);
t535 = -0.2e1 * t583 * t700 + t602 * (pkin(2) * t583 - t702);
t510 = (t534 * t552 + t535 * t549) * t546;
t699 = t510 * t546;
t536 = 0.2e1 * t581 * t703 - t603 * (pkin(2) * t581 + t704);
t537 = -0.2e1 * t584 * t703 + t603 * (pkin(2) * t584 - t705);
t511 = (t536 * t553 + t537 * t550) * t547;
t698 = t511 * t547;
t539 = -0.2e1 * t585 * t706 + t604 * (pkin(2) * t585 - t708);
t689 = t539 * t548;
t538 = 0.2e1 * t582 * t706 - t604 * (pkin(2) * t582 + t707);
t690 = t538 * t548;
t512 = t551 * t689 + t554 * t690;
t697 = t512 * t548;
t627 = pkin(2) * t628;
t498 = t616 * t510;
t644 = pkin(2) * t723 - t628 * t510 - t498;
t668 = -t628 + t720;
t447 = (t644 * t700 + ((t627 + (1 + t616) * pkin(2)) * t510 - t728 + t668 * t728) * t602) * t699;
t504 = pkin(2) * t510;
t468 = t504 - t728;
t450 = ((t468 - t728) * t674 + (-t591 * t510 + t510 - t644) * qJ(3,3)) * t699;
t677 = t591 * qJ(3,3);
t453 = ((0.2e1 * (t504 + (-t696 / 0.2e1 - t693 / 0.2e1) * t546) * t677 - (pkin(2) * t468 - t498) * t674 - qJ(3,3) * t728) * t546 + (-qJ(3,3) + t661 - t677) * t546 * t728) * t510;
t645 = (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,2) + Icges(2,3);
t649 = (rSges(3,3) ^ 2) + t628 + (0.2e1 * pkin(2) + rSges(3,1)) * rSges(3,1);
t540 = ((rSges(3,3) * t724) + t616 + t649) * m(3) + t645;
t570 = m(2) * rSges(2,1) + t709;
t609 = m(2) * rSges(2,2);
t655 = t609 - t712;
t543 = t570 * t605 - t602 * t655;
t438 = -t543 * t447 + t450 * t709 - t540 * t453;
t688 = t546 * t438;
t589 = m(1) + m(2) + m(3);
t441 = -t589 * t447 + t450 * t716 - t543 * t453;
t687 = t546 * t441;
t444 = (t447 * t605 + t453 * t608 - t450) * m(3);
t686 = t546 * t444;
t685 = t546 * t594;
t499 = t617 * t511;
t643 = pkin(2) * t722 - t628 * t511 - t499;
t448 = (t643 * t703 + ((t627 + (1 + t617) * pkin(2)) * t511 - t729 + t668 * t729) * t603) * t698;
t505 = pkin(2) * t511;
t469 = t505 - t729;
t451 = ((t469 - t729) * t673 + (-t592 * t511 + t511 - t643) * qJ(3,2)) * t698;
t676 = t592 * qJ(3,2);
t454 = ((0.2e1 * (t505 + (-t695 / 0.2e1 - t692 / 0.2e1) * t547) * t676 - (pkin(2) * t469 - t499) * t673 - qJ(3,2) * t729) * t547 + (-qJ(3,2) + t660 - t676) * t547 * t729) * t511;
t541 = ((rSges(3,3) * t725) + t617 + t649) * m(3) + t645;
t654 = t609 - t711;
t544 = t570 * t606 - t603 * t654;
t439 = -t544 * t448 + t451 * t709 - t541 * t454;
t684 = t547 * t439;
t442 = -t589 * t448 + t451 * t715 - t544 * t454;
t683 = t547 * t442;
t445 = (t448 * t606 + t454 * t608 - t451) * m(3);
t682 = t547 * t445;
t681 = t547 * t594;
t500 = t618 * t512;
t642 = pkin(2) * t721 - t628 * t512 - t500;
t449 = (t642 * t706 + ((t627 + (1 + t618) * pkin(2)) * t512 - t730 + t668 * t730) * t604) * t697;
t506 = pkin(2) * t512;
t470 = t506 - t730;
t452 = ((t470 - t730) * t672 + (-t593 * t512 + t512 - t642) * qJ(3,1)) * t697;
t675 = t593 * qJ(3,1);
t455 = ((0.2e1 * (t506 + (-t694 / 0.2e1 - t691 / 0.2e1) * t548) * t675 - (pkin(2) * t470 - t500) * t672 - qJ(3,1) * t730) * t548 + (-qJ(3,1) + t659 - t675) * t548 * t730) * t512;
t653 = t609 - t710;
t545 = t570 * t607 - t604 * t653;
t443 = -t589 * t449 + t452 * t714 - t545 * t455;
t680 = t548 * t443;
t446 = (t449 * t607 + t455 * t608 - t452) * m(3);
t679 = t548 * t446;
t678 = t548 * t594;
t667 = t510 ^ 2 * t712;
t666 = t511 ^ 2 * t711;
t665 = t512 ^ 2 * t710;
t465 = (t570 * t510 + t727 * t728) * t602 + t510 * t605 * t655;
t658 = t465 * t699;
t466 = (t570 * t511 + t727 * t729) * t603 + t511 * t606 * t654;
t657 = t466 * t698;
t467 = (t570 * t512 + t727 * t730) * t604 + t512 * t607 * t653;
t656 = t467 * t697;
t652 = t546 * t667;
t651 = t547 * t666;
t650 = t548 * t665;
t648 = t510 * t712 * t723;
t647 = t511 * t711 * t722;
t646 = t512 * t710 * t721;
t638 = t546 * t648;
t637 = t547 * t647;
t636 = t548 * t646;
t615 = rSges(4,1);
t614 = rSges(4,2);
t560 = t579 * t582 + 0.2e1 * t664;
t559 = t578 * t581 + 0.2e1 * t663;
t558 = t577 * t580 + 0.2e1 * t662;
t557 = t585 * t579 - 0.2e1 * t573;
t556 = t584 * t578 - 0.2e1 * t572;
t555 = t583 * t577 - 0.2e1 * t571;
t542 = ((rSges(3,3) * t726) + t618 + t649) * m(3) + t645;
t527 = -t557 * t672 + t560 * t593 - t582 * t628 - t582 - t664;
t526 = -t556 * t673 + t559 * t592 - t581 * t628 - t581 - t663;
t525 = -t555 * t674 + t558 * t591 - t580 * t628 - t580 - t662;
t524 = t557 * t593 + t560 * t672 - t628 * t585 + t573 - t585;
t523 = t556 * t592 + t559 * t673 - t628 * t584 + t572 - t584;
t522 = t555 * t591 + t558 * t674 - t628 * t583 + t571 - t583;
t515 = (-t538 * t563 + t539 * t566) * t548;
t514 = (-t536 * t562 + t537 * t565) * t547;
t513 = (-t534 * t561 + t535 * t564) * t546;
t497 = (-t530 * t563 + t533 * t566) * t548;
t496 = (-t529 * t562 + t532 * t565) * t547;
t495 = (-t528 * t561 + t531 * t564) * t546;
t491 = (-t524 * t563 + t527 * t566) * t548;
t490 = (-t523 * t562 + t526 * t565) * t547;
t489 = (-t522 * t561 + t525 * t564) * t546;
t488 = (-t527 * t607 - t539 * t608 + t533) * t717;
t487 = (-t526 * t606 - t537 * t608 + t532) * t718;
t486 = (-t525 * t605 - t535 * t608 + t531) * t719;
t485 = (-t524 * t607 - t538 * t608 + t530) * t717;
t484 = (-t523 * t606 - t536 * t608 + t529) * t718;
t483 = (-t522 * t605 - t534 * t608 + t528) * t719;
t482 = (t527 * t589 - t533 * t714 + t539 * t545) * t548;
t481 = (t526 * t589 - t532 * t715 + t537 * t544) * t547;
t480 = (t525 * t589 - t531 * t716 + t535 * t543) * t546;
t479 = (t524 * t589 - t530 * t714 + t538 * t545) * t548;
t478 = (t523 * t589 - t529 * t715 + t536 * t544) * t547;
t477 = (t522 * t589 - t528 * t716 + t534 * t543) * t546;
t476 = (t527 * t545 - t533 * t709 + t539 * t542) * t548;
t475 = (t526 * t544 - t532 * t709 + t537 * t541) * t547;
t474 = (t525 * t543 - t531 * t709 + t535 * t540) * t546;
t473 = (t524 * t545 - t530 * t709 + t538 * t542) * t548;
t472 = (t523 * t544 - t529 * t709 + t536 * t541) * t547;
t471 = (t522 * t543 - t528 * t709 + t534 * t540) * t546;
t464 = (-t491 * t607 - t515 * t608 + t497) * m(3);
t463 = (-t490 * t606 - t514 * t608 + t496) * m(3);
t462 = (-t489 * t605 - t513 * t608 + t495) * m(3);
t461 = t491 * t589 - t497 * t714 + t515 * t545;
t460 = t490 * t589 - t496 * t715 + t514 * t544;
t459 = t489 * t589 - t495 * t716 + t513 * t543;
t458 = t491 * t545 - t497 * t709 + t515 * t542;
t457 = t490 * t544 - t496 * t709 + t514 * t541;
t456 = t489 * t543 - t495 * t709 + t513 * t540;
t440 = -t545 * t449 + t452 * t709 - t542 * t455;
t1 = [(-(t473 * t538 + t479 * t524 + t485 * t530) * t566 - (t473 * t539 + t479 * t527 + t485 * t533) * t563) * t678 + t524 * t680 + t440 * t690 + t530 * t679 - t524 * t656 + t538 * t636 - t530 * t650 + (-(t472 * t536 + t478 * t523 + t484 * t529) * t565 - (t472 * t537 + t478 * t526 + t484 * t532) * t562) * t681 + t523 * t683 + t536 * t684 + t529 * t682 - t523 * t657 + t536 * t637 - t529 * t651 + (-(t471 * t534 + t477 * t522 + t483 * t528) * t564 - (t471 * t535 + t477 * t525 + t483 * t531) * t561) * t685 + t522 * t687 + t534 * t688 + t528 * t686 - t522 * t658 + t534 * t638 - t528 * t652 - (-t587 * t614 + t588 * t615) * t713; (-(t476 * t538 + t482 * t524 + t488 * t530) * t566 - (t476 * t539 + t482 * t527 + t488 * t533) * t563) * t678 + t527 * t680 + t440 * t689 + t533 * t679 - t527 * t656 + t539 * t636 - t533 * t650 + (-(t475 * t536 + t481 * t523 + t487 * t529) * t565 - (t475 * t537 + t481 * t526 + t487 * t532) * t562) * t681 + t526 * t683 + t537 * t684 + t532 * t682 - t526 * t657 + t537 * t637 - t532 * t651 + (-(t474 * t534 + t480 * t522 + t486 * t528) * t564 - (t474 * t535 + t480 * t525 + t486 * t531) * t561) * t685 + t525 * t687 + t535 * t688 + t531 * t686 - t525 * t658 + t535 * t638 - t531 * t652 - (t587 * t615 + t588 * t614) * t713; (-(t458 * t538 + t461 * t524 + t464 * t530) * t566 - (t458 * t539 + t461 * t527 + t464 * t533) * t563) * t678 + (-(t457 * t536 + t460 * t523 + t463 * t529) * t565 - (t457 * t537 + t460 * t526 + t463 * t532) * t562) * t681 + (-(t456 * t534 + t459 * t522 + t462 * t528) * t564 - (t456 * t535 + t459 * t525 + t462 * t531) * t561) * t685 + (t440 + t646) * t515 + (t439 + t647) * t514 + (t438 + t648) * t513 + (t446 - t665) * t497 + (t445 - t666) * t496 + (t444 - t667) * t495 + (-t467 * t512 + t443) * t491 + (-t466 * t511 + t442) * t490 + (-t465 * t510 + t441) * t489;];
taucX  = t1;
