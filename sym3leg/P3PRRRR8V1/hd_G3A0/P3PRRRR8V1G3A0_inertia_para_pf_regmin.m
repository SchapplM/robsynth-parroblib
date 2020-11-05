% Calculate minimal parameter regressor of inertia matrix for parallel robot
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3*3x12]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 17:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRRR8V1G3A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G3A0_inertia_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G3A0_inertia_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G3A0_inertia_para_pf_regmin: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G3A0_inertia_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G3A0_inertia_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:16:29
% EndTime: 2020-08-06 17:16:37
% DurationCPUTime: 7.26s
% Computational Cost: add. (5621->307), mult. (18459->808), div. (1596->14), fcn. (19769->22), ass. (0->371)
t548 = cos(qJ(3,3));
t529 = 0.1e1 / t548;
t550 = cos(qJ(3,2));
t531 = 0.1e1 / t550;
t552 = cos(qJ(3,1));
t533 = 0.1e1 / t552;
t554 = 0.1e1 / pkin(2);
t535 = sin(pkin(6));
t537 = cos(pkin(6));
t547 = sin(qJ(2,1));
t553 = cos(qJ(2,1));
t538 = cos(pkin(3));
t730 = t538 * t553;
t733 = t538 * t547;
t777 = pkin(2) * t552;
t465 = (-t535 * t547 + t537 * t730) * t777 + pkin(5) * (t535 * t553 + t537 * t733);
t468 = (t535 * t730 + t537 * t547) * t777 + (t535 * t733 - t537 * t553) * pkin(5);
t536 = sin(pkin(3));
t546 = sin(qJ(3,1));
t722 = t546 * t547;
t598 = t536 * t552 + t538 * t722;
t721 = t546 * t553;
t477 = t535 * t598 - t537 * t721;
t480 = t535 * t721 + t537 * t598;
t624 = t465 * t477 + t468 * t480;
t541 = legFrame(1,2);
t525 = cos(t541);
t720 = t547 * t552;
t507 = pkin(2) * t720 - pkin(5) * t553;
t734 = t538 * t546;
t682 = pkin(2) * t734 + t507 * t536;
t806 = 0.1e1 / t682 ^ 2;
t751 = t806 * t525;
t585 = t624 * t751;
t545 = sin(qJ(2,2));
t551 = cos(qJ(2,2));
t731 = t538 * t551;
t735 = t538 * t545;
t778 = pkin(2) * t550;
t464 = (-t535 * t545 + t537 * t731) * t778 + pkin(5) * (t535 * t551 + t537 * t735);
t467 = (t535 * t731 + t537 * t545) * t778 + (t535 * t735 - t537 * t551) * pkin(5);
t544 = sin(qJ(3,2));
t725 = t544 * t545;
t599 = t536 * t550 + t538 * t725;
t724 = t544 * t551;
t476 = t535 * t599 - t537 * t724;
t479 = t535 * t724 + t537 * t599;
t625 = t464 * t476 + t467 * t479;
t540 = legFrame(2,2);
t524 = cos(t540);
t723 = t545 * t550;
t506 = pkin(2) * t723 - pkin(5) * t551;
t736 = t538 * t544;
t683 = pkin(2) * t736 + t506 * t536;
t804 = 0.1e1 / t683 ^ 2;
t755 = t804 * t524;
t587 = t625 * t755;
t543 = sin(qJ(2,3));
t549 = cos(qJ(2,3));
t732 = t538 * t549;
t737 = t538 * t543;
t779 = pkin(2) * t548;
t463 = pkin(5) * (t535 * t549 + t537 * t737) + (-t535 * t543 + t537 * t732) * t779;
t466 = (t535 * t732 + t537 * t543) * t779 + (t535 * t737 - t537 * t549) * pkin(5);
t542 = sin(qJ(3,3));
t728 = t542 * t543;
t600 = t536 * t548 + t538 * t728;
t727 = t542 * t549;
t475 = t535 * t600 - t537 * t727;
t478 = t535 * t727 + t537 * t600;
t626 = t463 * t475 + t466 * t478;
t539 = legFrame(3,2);
t523 = cos(t539);
t726 = t543 * t548;
t505 = pkin(2) * t726 - pkin(5) * t549;
t738 = t538 * t542;
t684 = pkin(2) * t738 + t505 * t536;
t802 = 0.1e1 / t684 ^ 2;
t758 = t802 * t523;
t589 = t626 * t758;
t814 = (-t529 * t589 - t531 * t587 - t533 * t585) * t554;
t534 = 0.1e1 / t552 ^ 2;
t785 = t534 * t546;
t532 = 0.1e1 / t550 ^ 2;
t786 = t532 * t544;
t530 = 0.1e1 / t548 ^ 2;
t787 = t530 * t542;
t813 = (-t585 * t785 - t587 * t786 - t589 * t787) * t554;
t812 = t463 * t466;
t687 = t536 * t727;
t748 = (t536 * t726 + t738) * t554;
t801 = 0.1e1 / t684;
t759 = t801 * t529;
t811 = t759 * (t463 * t748 + t478 * t687);
t690 = (-t536 * t728 + t538 * t548) * t529 * t554;
t742 = t536 * t549;
t714 = t478 * t742;
t810 = t801 * (t463 * t690 + t714);
t460 = t463 ^ 2;
t555 = 0.1e1 / pkin(2) ^ 2;
t749 = t806 * t534;
t678 = t465 * t468 * t749;
t753 = t804 * t532;
t680 = t464 * t467 * t753;
t757 = t802 * t530;
t707 = t523 * t757;
t807 = (-t524 * t680 - t525 * t678 - t707 * t812) * t555;
t805 = 0.1e1 / t682;
t803 = 0.1e1 / t683;
t518 = t524 ^ 2;
t800 = 0.2e1 * t518;
t519 = t525 ^ 2;
t799 = 0.2e1 * t519;
t798 = t464 * t479;
t797 = t465 * t480;
t796 = t467 * t476;
t795 = t468 * t477;
t794 = t476 * t479;
t793 = t477 * t480;
t508 = pkin(5) * t543 + t549 * t779;
t780 = pkin(2) * t536;
t622 = -t505 * t538 + t542 * t780;
t453 = t537 * t508 + t622 * t535;
t520 = sin(t539);
t447 = t453 * t523 + t520 * t684;
t448 = -t453 * t520 + t523 * t684;
t792 = t478 * (t447 * t520 - t448 * t523);
t509 = pkin(5) * t545 + t551 * t778;
t621 = -t506 * t538 + t544 * t780;
t454 = t537 * t509 + t535 * t621;
t521 = sin(t540);
t449 = t454 * t524 + t521 * t683;
t450 = -t454 * t521 + t524 * t683;
t791 = t479 * (t449 * t521 - t450 * t524);
t510 = pkin(5) * t547 + t553 * t777;
t620 = -t507 * t538 + t546 * t780;
t455 = t537 * t510 + t535 * t620;
t522 = sin(t541);
t451 = t455 * t525 + t522 * t682;
t452 = -t455 * t522 + t525 * t682;
t790 = t480 * (t451 * t522 - t452 * t525);
t689 = (-t536 * t725 + t538 * t550) * t531 * t554;
t741 = t536 * t551;
t713 = t479 * t741;
t789 = t803 * (t464 * t689 + t713);
t688 = (-t536 * t722 + t538 * t552) * t533 * t554;
t740 = t536 * t553;
t712 = t480 * t740;
t788 = t805 * (t465 * t688 + t712);
t685 = t536 * t721;
t746 = (t536 * t720 + t734) * t554;
t752 = t805 * t533;
t784 = (t465 * t746 + t480 * t685) * t752;
t686 = t536 * t724;
t747 = (t536 * t723 + t736) * t554;
t756 = t803 * t531;
t783 = (t464 * t747 + t479 * t686) * t756;
t782 = 0.2e1 * t536;
t781 = 0.2e1 * t554;
t776 = t447 * t801;
t775 = t448 * t801;
t774 = t449 * t803;
t773 = t450 * t803;
t772 = t451 * t805;
t771 = t452 * t805;
t456 = -t535 * t508 + t622 * t537;
t770 = t456 * t478;
t769 = t456 * t801;
t768 = t456 * t802;
t457 = -t535 * t509 + t537 * t621;
t767 = t457 * t479;
t766 = t457 * t803;
t765 = t457 * t804;
t458 = -t535 * t510 + t537 * t620;
t764 = t458 * t480;
t763 = t458 * t805;
t762 = t458 * t806;
t760 = t478 * t802;
t754 = t804 * t531;
t750 = t806 * t533;
t745 = t529 * t542;
t744 = t531 * t544;
t743 = t533 * t546;
t739 = t536 * t554;
t729 = t538 * t554;
t719 = -0.2e1 * t463;
t718 = t475 * t757;
t717 = t475 * t742;
t716 = t476 * t741;
t715 = t477 * t740;
t711 = t801 * t759;
t514 = t520 ^ 2;
t710 = t514 * t757;
t517 = t523 ^ 2;
t709 = t517 * t757;
t708 = t520 * t758;
t526 = t542 ^ 2;
t706 = t526 * t757;
t705 = t802 * t745;
t704 = t803 * t756;
t515 = t521 ^ 2;
t703 = t515 * t753;
t702 = t518 * t753;
t701 = t521 * t755;
t527 = t544 ^ 2;
t700 = t527 * t753;
t699 = t804 * t744;
t698 = t544 * t753;
t697 = t805 * t752;
t516 = t522 ^ 2;
t696 = t516 * t749;
t695 = t519 * t749;
t694 = t522 * t751;
t528 = t546 ^ 2;
t693 = t528 * t749;
t692 = t806 * t743;
t691 = t546 * t749;
t681 = t463 * t514 * t760;
t679 = t754 * t798;
t677 = t750 * t797;
t472 = t478 ^ 2;
t675 = t472 * t706;
t674 = t472 * t705;
t473 = t479 ^ 2;
t673 = t473 * t700;
t672 = t473 * t699;
t474 = t480 ^ 2;
t671 = t474 * t693;
t670 = t474 * t692;
t669 = t478 * t718;
t668 = t753 * t794;
t667 = t749 * t793;
t666 = t478 * t711;
t664 = t479 * t704;
t663 = t480 * t697;
t662 = t543 * t711;
t661 = t549 * t711;
t660 = t801 * t543 * t739;
t659 = t529 * t708;
t658 = t520 * t707;
t657 = t545 * t704;
t656 = t551 * t704;
t655 = t803 * t545 * t739;
t654 = t531 * t701;
t653 = t532 * t701;
t652 = t547 * t697;
t651 = t553 * t697;
t650 = t805 * t547 * t739;
t649 = t533 * t694;
t648 = t534 * t694;
t647 = t456 * t475 * t711;
t646 = t457 * t476 * t704;
t645 = t458 * t477 * t697;
t644 = t463 * t660;
t643 = t463 * t658;
t642 = t698 * t798;
t641 = t464 * t655;
t640 = t691 * t797;
t639 = t465 * t650;
t638 = t472 * t658;
t637 = t473 * t653;
t636 = t474 * t648;
t635 = t526 * t669;
t634 = t475 * t478 * t705;
t633 = t527 * t668;
t632 = t699 * t794;
t631 = t528 * t667;
t630 = t692 * t793;
t623 = t517 * t719 * t760;
t619 = t447 * t523 * t666;
t618 = t448 * t520 * t666;
t617 = t449 * t524 * t664;
t616 = t450 * t521 * t664;
t615 = t451 * t525 * t663;
t614 = t452 * t522 * t663;
t613 = t654 * t798;
t612 = t649 * t797;
t611 = t660 * t745;
t610 = t655 * t744;
t609 = t650 * t743;
t608 = -t447 * t475 + t523 * t770;
t607 = t448 * t475 + t520 * t770;
t606 = -t449 * t476 + t524 * t767;
t605 = t450 * t476 + t521 * t767;
t604 = -t451 * t477 + t525 * t764;
t603 = t452 * t477 + t522 * t764;
t597 = t478 * t719 * t708;
t596 = t463 * t611;
t595 = t521 * t524 * t642;
t594 = t464 * t610;
t593 = t522 * t525 * t640;
t592 = t465 * t609;
t591 = t802 * t520 * t626;
t588 = t804 * t521 * t625;
t586 = t806 * t522 * t624;
t580 = (t466 * t690 + t717) * t801;
t579 = (t467 * t689 + t716) * t803;
t578 = (t468 * t688 + t715) * t805;
t577 = (t463 * t729 + t714) * t759;
t576 = (t464 * t729 + t713) * t756;
t575 = (t465 * t729 + t712) * t752;
t574 = t523 * t810;
t573 = t520 * t810;
t572 = t524 * t789;
t571 = t521 * t789;
t570 = t525 * t788;
t569 = t522 * t788;
t568 = (-t466 * t748 - t475 * t687) * t759;
t567 = (-t467 * t747 - t476 * t686) * t756;
t566 = (-t468 * t746 - t477 * t685) * t752;
t461 = t464 ^ 2;
t462 = t465 ^ 2;
t565 = -t461 * t653 - t462 * t648;
t563 = t523 * t811;
t562 = t520 * t811;
t561 = t521 * t783;
t560 = t524 * t783;
t559 = t522 * t784;
t558 = t525 * t784;
t471 = t477 ^ 2;
t470 = t476 ^ 2;
t469 = t475 ^ 2;
t446 = (t468 * t729 + t715) * t752;
t445 = (t467 * t729 + t716) * t756;
t444 = (t466 * t729 + t717) * t759;
t443 = t525 * t575;
t442 = t522 * t575;
t441 = t524 * t576;
t440 = t521 * t576;
t439 = t523 * t577;
t438 = t520 * t577;
t437 = -t636 - t637 - t638;
t436 = -t526 * t638 - t527 * t637 - t528 * t636;
t435 = -t446 * t546 - t468 * t650;
t434 = -t445 * t544 - t467 * t655;
t433 = -t444 * t542 - t466 * t660;
t432 = -0.2e1 * t472 * t542 * t659 - 0.2e1 * t473 * t544 * t654 - 0.2e1 * t474 * t546 * t649;
t431 = t446 * t552 - t468 * t609;
t430 = t445 * t550 - t467 * t610;
t429 = t444 * t548 - t466 * t611;
t428 = t443 * t546 + t525 * t639;
t427 = -t442 * t546 - t522 * t639;
t426 = t441 * t544 + t524 * t641;
t425 = -t440 * t544 - t521 * t641;
t424 = t439 * t542 + t523 * t644;
t423 = -t438 * t542 - t520 * t644;
t422 = -t443 * t552 + t525 * t592;
t421 = t442 * t552 - t522 * t592;
t420 = -t441 * t550 + t524 * t594;
t419 = t440 * t550 - t521 * t594;
t418 = -t439 * t548 + t523 * t596;
t417 = t438 * t548 - t520 * t596;
t416 = -t523 * t669 - t524 * t668 - t525 * t667;
t415 = t520 * t669 + t521 * t668 + t522 * t667;
t414 = -t523 * t635 - t524 * t633 - t525 * t631;
t413 = t520 * t635 + t521 * t633 + t522 * t631;
t412 = -0.2e1 * t523 * t634 - 0.2e1 * t524 * t632 - 0.2e1 * t525 * t630;
t411 = 0.2e1 * t520 * t634 + 0.2e1 * t521 * t632 + 0.2e1 * t522 * t630;
t410 = (t520 * t757 * t812 + t521 * t680 + t522 * t678) * t555;
t409 = t448 * t768 + t450 * t765 + t452 * t762;
t408 = t447 * t768 + t449 * t765 + t451 * t762;
t407 = t447 * t448 * t802 + t449 * t450 * t804 + t451 * t452 * t806;
t406 = (t529 * t591 + t531 * t588 + t533 * t586) * t554;
t405 = (t586 * t785 + t588 * t786 + t591 * t787) * t554;
t404 = (t603 * t651 + t605 * t656 + t607 * t661) * t536;
t403 = (-t603 * t652 - t605 * t657 - t607 * t662) * t536;
t402 = (-t604 * t651 - t606 * t656 - t608 * t661) * t536;
t401 = (t604 * t652 + t606 * t657 + t608 * t662) * t536;
t400 = (t651 * t790 + t656 * t791 + t661 * t792) * t536;
t399 = (-t652 * t790 - t657 * t791 - t662 * t792) * t536;
t1 = [t447 ^ 2 * t802 + t449 ^ 2 * t804 + t451 ^ 2 * t806, t472 * t709 + t473 * t702 + t474 * t695, (-t549 * t619 - t551 * t617 - t553 * t615) * t782, (t543 * t619 + t545 * t617 + t547 * t615) * t782, t517 * t675 + t518 * t673 + t519 * t671, 0.2e1 * t517 * t674 + 0.2e1 * t518 * t672 + 0.2e1 * t519 * t670, (-t623 * t787 + t640 * t799 + t642 * t800) * t554, (-t529 * t623 + t677 * t799 + t679 * t800) * t554, (t460 * t709 + t461 * t702 + t462 * t695) * t555, (t422 - t570) * t772 + (t420 - t572) * t774 + (t418 - t574) * t776, (t428 + t558) * t772 + (t426 + t560) * t774 + (t424 + t563) * t776, 1; t407, t437, t400, t399, t436, t432, (t597 * t787 - 0.2e1 * t593 - 0.2e1 * t595) * t554, (t529 * t597 - 0.2e1 * t612 - 0.2e1 * t613) * t554, (-t463 * t643 + t565) * t555, (t451 * t421 - t452 * t570) * t805 + (t449 * t419 - t450 * t572) * t803 + (t447 * t417 - t448 * t574) * t801, (t451 * t427 + t452 * t558) * t805 + (t449 * t425 + t450 * t560) * t803 + (t447 * t423 + t448 * t563) * t801, 0; t408, t416, t402, t401, t414, t412, t813, t814, t807, (t451 * t431 - t458 * t570) * t805 + (t449 * t430 - t457 * t572) * t803 + (t447 * t429 - t456 * t574) * t801, (t451 * t435 + t458 * t558) * t805 + (t449 * t434 + t457 * t560) * t803 + (t447 * t433 + t456 * t563) * t801, 0; t407, t437, t400, t399, t436, t432, (-t478 * t542 * t643 - t593 - t595) * t781, (-t463 * t478 * t659 - t612 - t613) * t781, (-t460 * t658 + t565) * t555, (t452 * t422 + t451 * t569) * t805 + (t450 * t420 + t449 * t571) * t803 + (t448 * t418 + t447 * t573) * t801, (t452 * t428 - t451 * t559) * t805 + (t450 * t426 - t449 * t561) * t803 + (t448 * t424 - t447 * t562) * t801, 0; t448 ^ 2 * t802 + t450 ^ 2 * t804 + t452 ^ 2 * t806, t472 * t710 + t473 * t703 + t474 * t696, (t549 * t618 + t551 * t616 + t553 * t614) * t782, (-t543 * t618 - t545 * t616 - t547 * t614) * t782, t514 * t675 + t515 * t673 + t516 * t671, 0.2e1 * t514 * t674 + 0.2e1 * t515 * t672 + 0.2e1 * t516 * t670, (t515 * t642 + t516 * t640 + t681 * t787) * t781, (t515 * t679 + t516 * t677 + t529 * t681) * t781, (t460 * t710 + t461 * t703 + t462 * t696) * t555, (t421 + t569) * t771 + (t419 + t571) * t773 + (t417 + t573) * t775, (t427 - t559) * t771 + (t425 - t561) * t773 + (t423 - t562) * t775, 1; t409, t415, t404, t403, t413, t411, t405, t406, t410, (t452 * t431 + t458 * t569) * t805 + (t450 * t430 + t457 * t571) * t803 + (t448 * t429 + t456 * t573) * t801, (t452 * t435 - t458 * t559) * t805 + (t450 * t434 - t457 * t561) * t803 + (t448 * t433 - t456 * t562) * t801, 0; t408, t416, t402, t401, t414, t412, t813, t814, t807, (t458 * t422 + t451 * t578) * t805 + (t457 * t420 + t449 * t579) * t803 + (t456 * t418 + t447 * t580) * t801, (t458 * t428 + t451 * t566) * t805 + (t457 * t426 + t449 * t567) * t803 + (t456 * t424 + t447 * t568) * t801, 0; t409, t415, t404, t403, t413, t411, t405, t406, t410, (t458 * t421 + t452 * t578) * t805 + (t457 * t419 + t450 * t579) * t803 + (t456 * t417 + t448 * t580) * t801, (t458 * t427 + t452 * t566) * t805 + (t457 * t425 + t450 * t567) * t803 + (t456 * t423 + t448 * t568) * t801, 0; t456 ^ 2 * t802 + t457 ^ 2 * t804 + t458 ^ 2 * t806, t469 * t757 + t470 * t753 + t471 * t749, (t549 * t647 + t551 * t646 + t553 * t645) * t782, (-t543 * t647 - t545 * t646 - t547 * t645) * t782, t469 * t706 + t470 * t700 + t471 * t693, 0.2e1 * t469 * t705 + 0.2e1 * t470 * t699 + 0.2e1 * t471 * t692, (t466 * t542 * t718 + t691 * t795 + t698 * t796) * t781, (t466 * t475 * t529 * t802 + t750 * t795 + t754 * t796) * t781, (t466 ^ 2 * t757 + t467 ^ 2 * t753 + t468 ^ 2 * t749) * t555, (t431 + t578) * t763 + (t430 + t579) * t766 + (t429 + t580) * t769, (t435 + t566) * t763 + (t434 + t567) * t766 + (t433 + t568) * t769, 1;];
tau_reg  = t1;