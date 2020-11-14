% Calculate inertia matrix for parallel robot
% P3PRRRR8V2G4A0
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2020-08-06 18:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRRR8V2G4A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G4A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G4A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G4A0_inertia_para_pf_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G4A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V2G4A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR8V2G4A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G4A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G4A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:14:55
% EndTime: 2020-08-06 18:14:57
% DurationCPUTime: 2.76s
% Computational Cost: add. (9030->336), mult. (18036->654), div. (432->10), fcn. (19692->34), ass. (0->299)
t668 = sin(qJ(3,3));
t797 = pkin(2) * t668;
t674 = cos(qJ(3,3));
t651 = t674 ^ 2;
t796 = pkin(3) * t651;
t676 = cos(qJ(3,2));
t652 = t676 ^ 2;
t795 = pkin(3) * t652;
t678 = cos(qJ(3,1));
t653 = t678 ^ 2;
t794 = pkin(3) * t653;
t655 = sin(pkin(4));
t793 = pkin(3) * t655;
t792 = pkin(3) * t674;
t791 = pkin(3) * t676;
t790 = pkin(3) * t678;
t670 = sin(qJ(3,2));
t789 = t670 * pkin(2);
t672 = sin(qJ(3,1));
t788 = t672 * pkin(2);
t787 = m(3) * pkin(2) + mrSges(2,1);
t675 = cos(qJ(2,3));
t669 = sin(qJ(2,3));
t681 = pkin(7) + pkin(6);
t720 = t669 * t681;
t613 = pkin(2) * t675 + t720;
t654 = sin(pkin(8));
t656 = cos(pkin(8));
t619 = t681 * t675;
t610 = pkin(2) * t669 - t619;
t657 = cos(pkin(4));
t688 = -t610 * t657 + t668 * t793;
t559 = -t613 * t656 - t688 * t654;
t562 = t613 * t654 - t688 * t656;
t726 = t657 * t668;
t583 = pkin(3) * t726 + t610 * t655;
t658 = legFrame(3,3);
t623 = sin(t658);
t629 = cos(t658);
t665 = legFrame(3,2);
t640 = sin(t665);
t643 = cos(t665);
t532 = t583 * t643 + (t559 * t629 + t562 * t623) * t640;
t535 = -t559 * t623 + t562 * t629;
t725 = t657 * t669;
t592 = t654 * t725 - t656 * t675;
t595 = t654 * t675 + t656 * t725;
t691 = t592 * t629 + t595 * t623;
t732 = t655 * t669;
t544 = t691 * t640 + t643 * t732;
t556 = -t592 * t623 + t595 * t629;
t589 = t623 * t656 + t629 * t654;
t565 = -t589 * t655 * t640 + t643 * t657;
t661 = legFrame(3,1);
t626 = sin(t661);
t632 = cos(t661);
t586 = -t623 * t654 + t629 * t656;
t747 = t586 * t655;
t511 = -(t544 * t626 - t556 * t632) * t796 + (-t532 * t626 + t535 * t632) * t674 - (t565 * t626 + t632 * t747) * t797;
t553 = 0.1e1 / (pkin(2) * t726 + t583 * t674 + t732 * t796);
t786 = t511 * t553;
t512 = (t544 * t632 + t556 * t626) * t796 + (t532 * t632 + t535 * t626) * t674 + (t565 * t632 - t626 * t747) * t797;
t785 = t512 * t553;
t677 = cos(qJ(2,2));
t671 = sin(qJ(2,2));
t719 = t671 * t681;
t614 = pkin(2) * t677 + t719;
t620 = t681 * t677;
t611 = pkin(2) * t671 - t620;
t687 = -t611 * t657 + t670 * t793;
t560 = -t614 * t656 - t687 * t654;
t563 = t614 * t654 - t687 * t656;
t724 = t657 * t670;
t584 = pkin(3) * t724 + t611 * t655;
t659 = legFrame(2,3);
t624 = sin(t659);
t630 = cos(t659);
t666 = legFrame(2,2);
t641 = sin(t666);
t644 = cos(t666);
t533 = t584 * t644 + (t560 * t630 + t563 * t624) * t641;
t536 = -t560 * t624 + t563 * t630;
t723 = t657 * t671;
t593 = t654 * t723 - t656 * t677;
t596 = t654 * t677 + t656 * t723;
t690 = t593 * t630 + t596 * t624;
t731 = t655 * t671;
t545 = t690 * t641 + t644 * t731;
t557 = -t593 * t624 + t596 * t630;
t590 = t624 * t656 + t630 * t654;
t566 = -t590 * t655 * t641 + t644 * t657;
t662 = legFrame(2,1);
t627 = sin(t662);
t633 = cos(t662);
t587 = -t624 * t654 + t630 * t656;
t746 = t587 * t655;
t513 = -(t545 * t627 - t557 * t633) * t795 + (-t533 * t627 + t536 * t633) * t676 - (t566 * t627 + t633 * t746) * t789;
t554 = 0.1e1 / (pkin(2) * t724 + t584 * t676 + t731 * t795);
t784 = t513 * t554;
t514 = (t545 * t633 + t557 * t627) * t795 + (t533 * t633 + t536 * t627) * t676 + (t566 * t633 - t627 * t746) * t789;
t783 = t514 * t554;
t679 = cos(qJ(2,1));
t673 = sin(qJ(2,1));
t718 = t673 * t681;
t615 = pkin(2) * t679 + t718;
t621 = t681 * t679;
t612 = pkin(2) * t673 - t621;
t686 = -t612 * t657 + t672 * t793;
t561 = -t615 * t656 - t686 * t654;
t564 = t615 * t654 - t686 * t656;
t722 = t657 * t672;
t585 = pkin(3) * t722 + t612 * t655;
t660 = legFrame(1,3);
t625 = sin(t660);
t631 = cos(t660);
t667 = legFrame(1,2);
t642 = sin(t667);
t645 = cos(t667);
t534 = t585 * t645 + (t561 * t631 + t564 * t625) * t642;
t537 = -t561 * t625 + t564 * t631;
t721 = t657 * t673;
t594 = t654 * t721 - t656 * t679;
t597 = t654 * t679 + t656 * t721;
t689 = t594 * t631 + t597 * t625;
t730 = t655 * t673;
t546 = t689 * t642 + t645 * t730;
t558 = -t594 * t625 + t597 * t631;
t591 = t625 * t656 + t631 * t654;
t567 = -t591 * t655 * t642 + t645 * t657;
t663 = legFrame(1,1);
t628 = sin(t663);
t634 = cos(t663);
t588 = -t625 * t654 + t631 * t656;
t745 = t588 * t655;
t515 = -(t546 * t628 - t558 * t634) * t794 + (-t534 * t628 + t537 * t634) * t678 - (t567 * t628 + t634 * t745) * t788;
t555 = 0.1e1 / (pkin(2) * t722 + t585 * t678 + t730 * t794);
t782 = t515 * t555;
t516 = (t546 * t634 + t558 * t628) * t794 + (t534 * t634 + t537 * t628) * t678 + (t567 * t634 - t628 * t745) * t788;
t781 = t516 * t555;
t738 = t632 * t640;
t547 = t586 * t738 - t589 * t626;
t729 = t655 * t674;
t520 = (t556 * t738 - t626 * t691) * t668 + t547 * t729;
t616 = pkin(2) + t792;
t607 = t616 * t726;
t568 = 0.1e1 / (t607 + (t669 * t792 + t610) * t729);
t780 = t520 * t568;
t737 = t633 * t641;
t548 = t587 * t737 - t590 * t627;
t728 = t655 * t676;
t521 = (t557 * t737 - t627 * t690) * t670 + t548 * t728;
t617 = pkin(2) + t791;
t608 = t617 * t724;
t569 = 0.1e1 / (t608 + (t671 * t791 + t611) * t728);
t779 = t521 * t569;
t736 = t634 * t642;
t549 = t588 * t736 - t591 * t628;
t727 = t655 * t678;
t522 = (t558 * t736 - t628 * t689) * t672 + t549 * t727;
t618 = pkin(2) + t790;
t609 = t618 * t722;
t570 = 0.1e1 / (t609 + (t673 * t790 + t612) * t727);
t778 = t522 * t570;
t741 = t626 * t640;
t550 = t586 * t741 + t589 * t632;
t523 = (-t556 * t741 - t691 * t632) * t668 - t550 * t729;
t777 = t523 * t568;
t740 = t627 * t641;
t551 = t587 * t740 + t590 * t633;
t524 = (-t557 * t740 - t690 * t633) * t670 - t551 * t728;
t776 = t524 * t569;
t739 = t628 * t642;
t552 = t588 * t739 + t591 * t634;
t525 = (-t558 * t739 - t689 * t634) * t672 - t552 * t727;
t775 = t525 * t570;
t604 = t616 * t669 - t619;
t744 = (t616 * t675 + t720) * t657;
t526 = -t550 * t744 + (-t586 * t632 + t589 * t741) * t604;
t571 = 0.1e1 / (t604 * t729 + t607);
t774 = t526 * t571;
t605 = t617 * t671 - t620;
t743 = (t617 * t677 + t719) * t657;
t527 = -t551 * t743 + (-t587 * t633 + t590 * t740) * t605;
t572 = 0.1e1 / (t605 * t728 + t608);
t773 = t527 * t572;
t606 = t618 * t673 - t621;
t742 = (t618 * t679 + t718) * t657;
t528 = -t552 * t742 + (-t588 * t634 + t591 * t739) * t606;
t573 = 0.1e1 / (t606 * t727 + t609);
t772 = t528 * t573;
t529 = t547 * t744 - (t586 * t626 + t589 * t738) * t604;
t771 = t529 * t571;
t530 = t548 * t743 - (t587 * t627 + t590 * t737) * t605;
t770 = t530 * t572;
t531 = t549 * t742 - (t588 * t628 + t591 * t736) * t606;
t769 = t531 * t573;
t538 = t586 * t729 + t668 * (t586 * t725 + t589 * t675);
t768 = t538 * t643;
t539 = t587 * t728 + t670 * (t587 * t723 + t590 * t677);
t767 = t539 * t644;
t540 = t588 * t727 + t672 * (t588 * t721 + t591 * t679);
t766 = t540 * t645;
t698 = t674 * mrSges(3,1) - mrSges(3,2) * t668;
t577 = t698 * t657 - (mrSges(3,1) * t668 + mrSges(3,2) * t674) * t732;
t765 = t553 * t577;
t649 = m(1) + m(2) + m(3);
t764 = t553 * t649;
t697 = t676 * mrSges(3,1) - mrSges(3,2) * t670;
t578 = t697 * t657 - (mrSges(3,1) * t670 + mrSges(3,2) * t676) * t731;
t763 = t554 * t578;
t762 = t554 * t649;
t696 = t678 * mrSges(3,1) - mrSges(3,2) * t672;
t579 = t696 * t657 - (mrSges(3,1) * t672 + mrSges(3,2) * t678) * t730;
t761 = t555 * t579;
t760 = t555 * t649;
t664 = -Ifges(3,1) + Ifges(3,2);
t680 = mrSges(3,1) * pkin(2);
t692 = (2 * pkin(6) * mrSges(3,3)) + Ifges(3,1) + Ifges(2,3) + (pkin(2) ^ 2 + (pkin(6) ^ 2)) * m(3);
t717 = -0.2e1 * mrSges(3,2) * pkin(2);
t574 = t664 * t651 + 0.2e1 * (Ifges(3,4) * t668 + t680) * t674 + t668 * t717 + t692;
t759 = t568 * t574;
t635 = -pkin(6) * mrSges(3,2) + Ifges(3,6);
t636 = mrSges(3,1) * pkin(6) - Ifges(3,5);
t598 = t635 * t674 - t636 * t668;
t758 = t568 * t598;
t575 = t664 * t652 + 0.2e1 * (Ifges(3,4) * t670 + t680) * t676 + t670 * t717 + t692;
t757 = t569 * t575;
t599 = t635 * t676 - t636 * t670;
t756 = t569 * t599;
t576 = t664 * t653 + 0.2e1 * (Ifges(3,4) * t672 + t680) * t678 + t672 * t717 + t692;
t755 = t570 * t576;
t600 = t635 * t678 - t636 * t672;
t754 = t570 * t600;
t685 = 0.1e1 / pkin(3);
t753 = t571 * t685;
t752 = t572 * t685;
t751 = t573 * t685;
t622 = m(3) * pkin(6) - mrSges(2,2) + mrSges(3,3);
t580 = (t698 + t787) * t675 + t622 * t669;
t750 = t580 * t655;
t581 = (t697 + t787) * t677 + t622 * t671;
t749 = t581 * t655;
t582 = (t696 + t787) * t679 + t622 * t673;
t748 = t582 * t655;
t735 = t643 * t655;
t734 = t644 * t655;
t733 = t645 * t655;
t716 = Ifges(3,3) * t753;
t715 = Ifges(3,3) * t752;
t714 = Ifges(3,3) * t751;
t713 = (-t586 * t744 + t589 * t604) * t571 * t643;
t712 = (-t587 * t743 + t590 * t605) * t572 * t644;
t711 = (-t588 * t742 + t591 * t606) * t573 * t645;
t710 = t553 * t750;
t709 = t554 * t749;
t708 = t555 * t748;
t707 = t568 * t750;
t706 = t569 * t749;
t705 = t570 * t748;
t704 = t577 * t753;
t703 = t598 * t753;
t702 = t578 * t752;
t701 = t599 * t752;
t700 = t579 * t751;
t699 = t600 * t751;
t695 = t685 * t713;
t694 = t685 * t712;
t693 = t685 * t711;
t519 = -((-t588 * t679 + t591 * t721) * t645 - t642 * t730) * t794 + ((t588 * t615 + t686 * t591) * t645 + t585 * t642) * t678 + (t591 * t733 + t642 * t657) * t788;
t518 = -((-t587 * t677 + t590 * t723) * t644 - t641 * t731) * t795 + ((t587 * t614 + t687 * t590) * t644 + t584 * t641) * t676 + (t590 * t734 + t641 * t657) * t789;
t517 = -((-t586 * t675 + t589 * t725) * t643 - t640 * t732) * t796 + ((t586 * t613 + t688 * t589) * t643 + t583 * t640) * t674 + (t589 * t735 + t640 * t657) * t797;
t510 = Ifges(3,3) * t693 + (t519 * t579 - t600 * t766) * t555;
t509 = Ifges(3,3) * t694 + (t518 * t578 - t599 * t767) * t554;
t508 = Ifges(3,3) * t695 + (t517 * t577 - t598 * t768) * t553;
t507 = t579 * t693 + (-t540 * t582 * t733 + t519 * t649) * t555;
t506 = t578 * t694 + (-t539 * t581 * t734 + t518 * t649) * t554;
t505 = t577 * t695 + (-t538 * t580 * t735 + t517 * t649) * t553;
t504 = t600 * t693 + (t519 * t748 - t576 * t766) * t555;
t503 = t599 * t694 + (t518 * t749 - t575 * t767) * t554;
t502 = t598 * t695 + (t517 * t750 - t574 * t768) * t553;
t501 = t516 * t761 + t522 * t754 + t531 * t714;
t500 = t515 * t761 + t525 * t754 + t528 * t714;
t499 = t514 * t763 + t521 * t756 + t530 * t715;
t498 = t513 * t763 + t524 * t756 + t527 * t715;
t497 = t512 * t765 + t520 * t758 + t529 * t716;
t496 = t511 * t765 + t523 * t758 + t526 * t716;
t495 = t516 * t760 + t522 * t705 + t531 * t700;
t494 = t515 * t760 + t525 * t705 + t528 * t700;
t493 = t514 * t762 + t521 * t706 + t530 * t702;
t492 = t513 * t762 + t524 * t706 + t527 * t702;
t491 = t512 * t764 + t520 * t707 + t529 * t704;
t490 = t511 * t764 + t523 * t707 + t526 * t704;
t489 = t516 * t708 + t522 * t755 + t531 * t699;
t488 = t515 * t708 + t525 * t755 + t528 * t699;
t487 = t514 * t709 + t521 * t757 + t530 * t701;
t486 = t513 * t709 + t524 * t757 + t527 * t701;
t485 = t512 * t710 + t520 * t759 + t529 * t703;
t484 = t511 * t710 + t523 * t759 + t526 * t703;
t1 = [m(4) + (-t504 * t766 + t507 * t519) * t555 + (-t503 * t767 + t506 * t518) * t554 + (-t502 * t768 + t505 * t517) * t553 + (t508 * t713 + t509 * t712 + t510 * t711) * t685, t502 * t777 + t503 * t776 + t504 * t775 + t505 * t786 + t506 * t784 + t507 * t782 + (t508 * t774 + t509 * t773 + t510 * t772) * t685, t502 * t780 + t503 * t779 + t504 * t778 + t505 * t785 + t506 * t783 + t507 * t781 + (t508 * t771 + t509 * t770 + t510 * t769) * t685; (-t488 * t766 + t494 * t519) * t555 + (-t486 * t767 + t492 * t518) * t554 + (-t484 * t768 + t490 * t517) * t553 + (t496 * t713 + t498 * t712 + t500 * t711) * t685, t484 * t777 + t486 * t776 + t488 * t775 + t490 * t786 + t492 * t784 + t494 * t782 + m(4) + (t496 * t774 + t498 * t773 + t500 * t772) * t685, t484 * t780 + t486 * t779 + t488 * t778 + t490 * t785 + t492 * t783 + t494 * t781 + (t496 * t771 + t498 * t770 + t500 * t769) * t685; (-t489 * t766 + t495 * t519) * t555 + (-t487 * t767 + t493 * t518) * t554 + (-t485 * t768 + t491 * t517) * t553 + (t497 * t713 + t499 * t712 + t501 * t711) * t685, t485 * t777 + t487 * t776 + t489 * t775 + t491 * t786 + t493 * t784 + t495 * t782 + (t497 * t774 + t499 * t773 + t501 * t772) * t685, t485 * t780 + t487 * t779 + t489 * t778 + t491 * t785 + t493 * t783 + t495 * t781 + m(4) + (t497 * t771 + t499 * t770 + t501 * t769) * t685;];
MX  = t1;
