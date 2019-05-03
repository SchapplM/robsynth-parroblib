% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRP1A0
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
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:42
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRP1A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRP1A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRP1A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRP1A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PRP1A0_coriolisvec_para_pf_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRP1A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRP1A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRP1A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRP1A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRP1A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:41:34
% EndTime: 2019-05-03 14:41:38
% DurationCPUTime: 4.08s
% Computational Cost: add. (26341->375), mult. (50050->587), div. (2250->3), fcn. (23653->14), ass. (0->263)
t616 = (pkin(2) ^ 2);
t583 = 1 + t616;
t606 = (qJ(3,1) ^ 2);
t567 = -t606 + t583;
t595 = cos(qJ(2,1));
t581 = t595 ^ 2;
t592 = sin(qJ(2,1));
t669 = t592 * t595;
t655 = pkin(2) * t669;
t724 = 2 * qJ(3,1);
t545 = 0.1e1 / (t567 * t581 + t655 * t724 - t583 - t606);
t589 = legFrame(1,3);
t570 = sin(t589);
t573 = cos(t589);
t707 = qJ(3,1) * t570;
t661 = pkin(2) * t707;
t627 = t573 * t606 - t661;
t706 = qJ(3,1) * t573;
t656 = pkin(2) * t706;
t628 = -t570 * t606 - t656;
t530 = t627 * t595 - t592 * (t570 - t628);
t601 = xP(3);
t575 = sin(t601);
t576 = cos(t601);
t609 = koppelP(1,2);
t612 = koppelP(1,1);
t560 = -t575 * t609 + t576 * t612;
t598 = xDP(3);
t599 = xDP(2);
t548 = t560 * t598 + t599;
t684 = t530 * t548;
t527 = t628 * t595 + t592 * (-t573 - t627);
t557 = t575 * t612 + t576 * t609;
t600 = xDP(1);
t551 = -t557 * t598 + t600;
t687 = t527 * t551;
t728 = (t684 + t687) * t545;
t605 = (qJ(3,2) ^ 2);
t566 = -t605 + t583;
t594 = cos(qJ(2,2));
t580 = t594 ^ 2;
t591 = sin(qJ(2,2));
t670 = t591 * t594;
t653 = pkin(2) * t670;
t723 = 2 * qJ(3,2);
t544 = 0.1e1 / (t566 * t580 + t653 * t723 - t583 - t605);
t588 = legFrame(2,3);
t569 = sin(t588);
t572 = cos(t588);
t704 = qJ(3,2) * t569;
t660 = pkin(2) * t704;
t624 = t572 * t605 - t660;
t703 = qJ(3,2) * t572;
t657 = pkin(2) * t703;
t625 = -t569 * t605 - t657;
t529 = t624 * t594 - t591 * (t569 - t625);
t608 = koppelP(2,2);
t611 = koppelP(2,1);
t559 = -t575 * t608 + t576 * t611;
t547 = t559 * t598 + t599;
t685 = t529 * t547;
t526 = t625 * t594 + t591 * (-t572 - t624);
t556 = t575 * t611 + t576 * t608;
t550 = -t556 * t598 + t600;
t688 = t526 * t550;
t727 = (t685 + t688) * t544;
t604 = (qJ(3,3) ^ 2);
t565 = -t604 + t583;
t593 = cos(qJ(2,3));
t579 = t593 ^ 2;
t590 = sin(qJ(2,3));
t671 = t590 * t593;
t654 = pkin(2) * t671;
t722 = 2 * qJ(3,3);
t543 = 0.1e1 / (t565 * t579 + t654 * t722 - t583 - t604);
t587 = legFrame(3,3);
t568 = sin(t587);
t571 = cos(t587);
t701 = qJ(3,3) * t568;
t659 = pkin(2) * t701;
t621 = t571 * t604 - t659;
t700 = qJ(3,3) * t571;
t658 = pkin(2) * t700;
t622 = -t568 * t604 - t658;
t528 = t621 * t593 - t590 * (t568 - t622);
t607 = koppelP(3,2);
t610 = koppelP(3,1);
t558 = -t575 * t607 + t576 * t610;
t546 = t558 * t598 + t599;
t686 = t528 * t546;
t525 = t622 * t593 + t590 * (-t571 - t621);
t555 = t575 * t610 + t576 * t607;
t549 = -t555 * t598 + t600;
t689 = t525 * t549;
t726 = (t686 + t689) * t543;
t725 = -2 * m(3);
t582 = t598 ^ 2;
t721 = 0.2e1 * t726;
t720 = 0.2e1 * t727;
t719 = 0.2e1 * t728;
t718 = m(3) * t543;
t717 = m(3) * t544;
t716 = m(3) * t545;
t715 = m(3) * t593;
t714 = m(3) * t594;
t713 = m(3) * t595;
t712 = m(4) * t582;
t711 = (rSges(3,3) + qJ(3,3)) * m(3);
t710 = (rSges(3,3) + qJ(3,2)) * m(3);
t709 = (rSges(3,3) + qJ(3,1)) * m(3);
t596 = pkin(2) + rSges(3,1);
t708 = t596 * m(3);
t705 = qJ(3,1) * t595;
t702 = qJ(3,2) * t594;
t699 = qJ(3,3) * t593;
t665 = -0.2e1 * t699;
t532 = t571 * t665 + t590 * (pkin(2) * t571 + t701);
t682 = t532 * t543;
t531 = t568 * t665 + t590 * (pkin(2) * t568 - t700);
t683 = t531 * t543;
t507 = t546 * t683 + t549 * t682;
t615 = pkin(2) * t616;
t495 = t604 * t507;
t638 = pkin(2) * t721 - t507 * t616 - t495;
t668 = -t616 - t583;
t692 = t507 * t543;
t444 = (t638 * t699 + ((t615 + (1 + t604) * pkin(2)) * t507 - t726 + t668 * t726) * t590) * t692;
t501 = pkin(2) * t507;
t465 = t501 - t726;
t447 = ((t465 - t726) * t671 + (-t507 * t579 + t507 - t638) * qJ(3,3)) * t692;
t674 = t579 * qJ(3,3);
t450 = ((0.2e1 * (t501 + (-t689 / 0.2e1 - t686 / 0.2e1) * t543) * t674 - (pkin(2) * t465 - t495) * t671 - qJ(3,3) * t726) * t543 + (-qJ(3,3) + t654 - t674) * t543 * t726) * t507;
t564 = m(2) * rSges(2,1) + t708;
t597 = m(2) * rSges(2,2);
t649 = t597 - t711;
t540 = t564 * t593 - t590 * t649;
t577 = m(1) + m(2) + m(3);
t438 = -t444 * t577 + t447 * t715 - t450 * t540;
t698 = t438 * t543;
t666 = -0.2e1 * t702;
t534 = t572 * t666 + t591 * (pkin(2) * t572 + t704);
t680 = t534 * t544;
t533 = t569 * t666 + t591 * (pkin(2) * t569 - t703);
t681 = t533 * t544;
t508 = t547 * t681 + t550 * t680;
t496 = t605 * t508;
t637 = pkin(2) * t720 - t508 * t616 - t496;
t691 = t508 * t544;
t445 = (t637 * t702 + ((t615 + (1 + t605) * pkin(2)) * t508 - t727 + t668 * t727) * t591) * t691;
t502 = pkin(2) * t508;
t466 = t502 - t727;
t448 = ((t466 - t727) * t670 + (-t508 * t580 + t508 - t637) * qJ(3,2)) * t691;
t673 = t580 * qJ(3,2);
t451 = ((0.2e1 * (t502 + (-t688 / 0.2e1 - t685 / 0.2e1) * t544) * t673 - (pkin(2) * t466 - t496) * t670 - qJ(3,2) * t727) * t544 + (-qJ(3,2) + t653 - t673) * t544 * t727) * t508;
t648 = t597 - t710;
t541 = t564 * t594 - t591 * t648;
t439 = -t445 * t577 + t448 * t714 - t451 * t541;
t697 = t439 * t544;
t667 = -0.2e1 * t705;
t536 = t573 * t667 + t592 * (pkin(2) * t573 + t707);
t678 = t536 * t545;
t535 = t570 * t667 + t592 * (pkin(2) * t570 - t706);
t679 = t535 * t545;
t509 = t548 * t679 + t551 * t678;
t497 = t606 * t509;
t636 = pkin(2) * t719 - t509 * t616 - t497;
t690 = t509 * t545;
t446 = (t636 * t705 + ((t615 + (1 + t606) * pkin(2)) * t509 - t728 + t668 * t728) * t592) * t690;
t503 = pkin(2) * t509;
t467 = t503 - t728;
t449 = ((t467 - t728) * t669 + (-t509 * t581 + t509 - t636) * qJ(3,1)) * t690;
t672 = t581 * qJ(3,1);
t452 = ((0.2e1 * (t503 + (-t687 / 0.2e1 - t684 / 0.2e1) * t545) * t672 - (pkin(2) * t467 - t497) * t669 - qJ(3,1) * t728) * t545 + (-qJ(3,1) + t655 - t672) * t545 * t728) * t509;
t647 = t597 - t709;
t542 = t564 * t595 - t592 * t647;
t440 = -t446 * t577 + t449 * t713 - t452 * t542;
t696 = t440 * t545;
t441 = (t444 * t593 + t450 * t596 - t447) * m(3);
t695 = t441 * t543;
t442 = (t445 * t594 + t451 * t596 - t448) * m(3);
t694 = t442 * t544;
t443 = (t446 * t595 + t452 * t596 - t449) * m(3);
t693 = t443 * t545;
t677 = t543 * t582;
t676 = t544 * t582;
t675 = t545 * t582;
t664 = t507 ^ 2 * t711;
t663 = t508 ^ 2 * t710;
t662 = t509 ^ 2 * t709;
t462 = (t507 * t564 + t725 * t726) * t590 + t649 * t593 * t507;
t652 = t462 * t692;
t463 = (t508 * t564 + t725 * t727) * t591 + t648 * t594 * t508;
t651 = t463 * t691;
t464 = (t509 * t564 + t725 * t728) * t592 + t647 * t595 * t509;
t650 = t464 * t690;
t646 = t543 * t664;
t645 = t544 * t663;
t644 = t545 * t662;
t643 = (rSges(3,3) ^ 2) + t616 + (0.2e1 * pkin(2) + rSges(3,1)) * rSges(3,1);
t642 = t507 * t711 * t721;
t641 = t508 * t710 * t720;
t640 = t509 * t709 * t719;
t639 = (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,2) + Icges(2,3);
t632 = t543 * t642;
t631 = t544 * t641;
t630 = t545 * t640;
t629 = -t567 * t570 + 0.2e1 * t656;
t626 = -t566 * t569 + 0.2e1 * t657;
t623 = -t565 * t568 + 0.2e1 * t658;
t603 = rSges(4,1);
t602 = rSges(4,2);
t554 = t567 * t573 + 0.2e1 * t661;
t553 = t566 * t572 + 0.2e1 * t660;
t552 = t565 * t571 + 0.2e1 * t659;
t539 = ((rSges(3,3) * t724) + t606 + t643) * m(3) + t639;
t538 = ((rSges(3,3) * t723) + t605 + t643) * m(3) + t639;
t537 = ((rSges(3,3) * t722) + t604 + t643) * m(3) + t639;
t524 = -t554 * t669 + t616 * t570 + t629 * t581 + t570 - t656;
t523 = t554 * t581 - t573 * t616 + t629 * t669 - t573 - t661;
t522 = -t553 * t670 + t616 * t569 + t626 * t580 + t569 - t657;
t521 = t553 * t580 - t572 * t616 + t626 * t670 - t572 - t660;
t520 = -t552 * t671 + t616 * t568 + t623 * t579 + t568 - t658;
t519 = t552 * t579 - t571 * t616 + t623 * t671 - t571 - t659;
t512 = (t535 * t560 - t536 * t557) * t545;
t511 = (t533 * t559 - t534 * t556) * t544;
t510 = (t531 * t558 - t532 * t555) * t543;
t494 = (-t527 * t557 + t530 * t560) * t545;
t493 = (-t526 * t556 + t529 * t559) * t544;
t492 = (-t525 * t555 + t528 * t558) * t543;
t488 = (t523 * t560 - t524 * t557) * t545;
t487 = (t521 * t559 - t522 * t556) * t544;
t486 = (t519 * t558 - t520 * t555) * t543;
t485 = (-t524 * t595 - t536 * t596 + t527) * t716;
t484 = (-t523 * t595 - t535 * t596 + t530) * t716;
t483 = (-t522 * t594 - t534 * t596 + t526) * t717;
t482 = (-t521 * t594 - t533 * t596 + t529) * t717;
t481 = (-t520 * t593 - t532 * t596 + t525) * t718;
t480 = (-t519 * t593 - t531 * t596 + t528) * t718;
t479 = (t524 * t577 - t527 * t713 + t536 * t542) * t545;
t478 = (t523 * t577 - t530 * t713 + t535 * t542) * t545;
t477 = (t522 * t577 - t526 * t714 + t534 * t541) * t544;
t476 = (t521 * t577 - t529 * t714 + t533 * t541) * t544;
t475 = (t520 * t577 - t525 * t715 + t532 * t540) * t543;
t474 = (t519 * t577 - t528 * t715 + t531 * t540) * t543;
t473 = (t524 * t542 - t527 * t708 + t536 * t539) * t545;
t472 = (t523 * t542 - t530 * t708 + t535 * t539) * t545;
t471 = (t522 * t541 - t526 * t708 + t534 * t538) * t544;
t470 = (t521 * t541 - t529 * t708 + t533 * t538) * t544;
t469 = (t520 * t540 - t525 * t708 + t532 * t537) * t543;
t468 = (t519 * t540 - t528 * t708 + t531 * t537) * t543;
t461 = (-t488 * t595 - t512 * t596 + t494) * m(3);
t460 = (-t487 * t594 - t511 * t596 + t493) * m(3);
t459 = (-t486 * t593 - t510 * t596 + t492) * m(3);
t458 = t488 * t577 - t494 * t713 + t512 * t542;
t457 = t487 * t577 - t493 * t714 + t511 * t541;
t456 = t486 * t577 - t492 * t715 + t510 * t540;
t455 = t488 * t542 - t494 * t708 + t512 * t539;
t454 = t487 * t541 - t493 * t708 + t511 * t538;
t453 = t486 * t540 - t492 * t708 + t510 * t537;
t437 = -t446 * t542 + t449 * t708 - t452 * t539;
t436 = -t445 * t541 + t448 * t708 - t451 * t538;
t435 = -t444 * t540 + t447 * t708 - t450 * t537;
t1 = [(-(t473 * t536 + t479 * t524 + t485 * t527) * t560 - (t473 * t535 + t479 * t523 + t485 * t530) * t557) * t675 + t524 * t696 + t437 * t678 + t527 * t693 - t524 * t650 + t536 * t630 - t527 * t644 + (-(t471 * t534 + t477 * t522 + t483 * t526) * t559 - (t471 * t533 + t477 * t521 + t483 * t529) * t556) * t676 + t522 * t697 + t436 * t680 + t526 * t694 - t522 * t651 + t534 * t631 - t526 * t645 + (-(t469 * t532 + t475 * t520 + t481 * t525) * t558 - (t469 * t531 + t475 * t519 + t481 * t528) * t555) * t677 + t520 * t698 + t435 * t682 + t525 * t695 - t520 * t652 + t532 * t632 - t525 * t646 - (-t575 * t602 + t576 * t603) * t712; (-(t472 * t536 + t478 * t524 + t484 * t527) * t560 - (t472 * t535 + t478 * t523 + t484 * t530) * t557) * t675 + t523 * t696 + t437 * t679 + t530 * t693 - t523 * t650 + t535 * t630 - t530 * t644 + (-(t470 * t534 + t476 * t522 + t482 * t526) * t559 - (t470 * t533 + t476 * t521 + t482 * t529) * t556) * t676 + t521 * t697 + t436 * t681 + t529 * t694 - t521 * t651 + t533 * t631 - t529 * t645 + (-(t468 * t532 + t474 * t520 + t480 * t525) * t558 - (t468 * t531 + t474 * t519 + t480 * t528) * t555) * t677 + t519 * t698 + t435 * t683 + t528 * t695 - t519 * t652 + t531 * t632 - t528 * t646 - (t575 * t603 + t576 * t602) * t712; (-(t455 * t536 + t458 * t524 + t461 * t527) * t560 - (t455 * t535 + t458 * t523 + t461 * t530) * t557) * t675 + (-(t454 * t534 + t457 * t522 + t460 * t526) * t559 - (t454 * t533 + t457 * t521 + t460 * t529) * t556) * t676 + (-(t453 * t532 + t456 * t520 + t459 * t525) * t558 - (t453 * t531 + t456 * t519 + t459 * t528) * t555) * t677 + (t437 + t640) * t512 + (t436 + t641) * t511 + (t435 + t642) * t510 + (t443 - t662) * t494 + (t442 - t663) * t493 + (t441 - t664) * t492 + (-t464 * t509 + t440) * t488 + (-t463 * t508 + t439) * t487 + (-t462 * t507 + t438) * t486;];
taucX  = t1;
