% Calculate inertia matrix for parallel robot
% P4PRRRR8V1G3A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
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
% MX [4x4]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-09-20 23:03
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P4PRRRR8V1G3A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V1G3A0_inertia_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V1G3A0_inertia_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P4PRRRR8V1G3A0_inertia_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V1G3A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRRR8V1G3A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P4PRRRR8V1G3A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V1G3A0_inertia_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V1G3A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-09-20 22:59:31
% EndTime: 2020-09-20 22:59:34
% DurationCPUTime: 2.88s
% Computational Cost: add. (8124->390), mult. (21960->699), div. (1920->9), fcn. (22920->30), ass. (0->308)
t870 = 2 * Ifges(3,4);
t704 = sin(qJ(2,4));
t706 = cos(qJ(2,4));
t705 = cos(qJ(3,4));
t864 = pkin(2) * t705;
t663 = -t706 * pkin(5) + t704 * t864;
t700 = sin(pkin(3));
t702 = cos(pkin(3));
t703 = sin(qJ(3,4));
t810 = t703 * t702;
t640 = pkin(2) * t810 + t663 * t700;
t869 = 0.1e1 / t640;
t714 = sin(qJ(2,3));
t720 = cos(qJ(2,3));
t719 = cos(qJ(3,3));
t863 = pkin(2) * t719;
t665 = -t720 * pkin(5) + t714 * t863;
t713 = sin(qJ(3,3));
t808 = t713 * t702;
t644 = pkin(2) * t808 + t665 * t700;
t868 = 0.1e1 / t644;
t716 = sin(qJ(2,2));
t722 = cos(qJ(2,2));
t721 = cos(qJ(3,2));
t862 = pkin(2) * t721;
t666 = -t722 * pkin(5) + t716 * t862;
t715 = sin(qJ(3,2));
t806 = t715 * t702;
t645 = pkin(2) * t806 + t666 * t700;
t867 = 0.1e1 / t645;
t718 = sin(qJ(2,1));
t724 = cos(qJ(2,1));
t723 = cos(qJ(3,1));
t861 = pkin(2) * t723;
t667 = -t724 * pkin(5) + t718 * t861;
t717 = sin(qJ(3,1));
t804 = t717 * t702;
t646 = pkin(2) * t804 + t667 * t700;
t866 = 0.1e1 / t646;
t865 = pkin(2) * t700;
t860 = Ifges(3,1) + Ifges(2,3);
t736 = 0.1e1 / pkin(2);
t859 = Ifges(3,3) * t736;
t664 = pkin(5) * t704 + t706 * t864;
t699 = sin(pkin(6));
t701 = cos(pkin(6));
t760 = -t663 * t702 + t703 * t865;
t607 = t701 * t664 + t699 * t760;
t709 = legFrame(4,2);
t684 = sin(t709);
t688 = cos(t709);
t599 = t607 * t688 + t640 * t684;
t858 = t599 * t869;
t600 = -t607 * t684 + t640 * t688;
t857 = t600 * t869;
t668 = pkin(5) * t714 + t720 * t863;
t759 = -t665 * t702 + t713 * t865;
t609 = t701 * t668 + t699 * t759;
t710 = legFrame(3,2);
t685 = sin(t710);
t689 = cos(t710);
t601 = t609 * t689 + t644 * t685;
t856 = t601 * t868;
t602 = -t609 * t685 + t644 * t689;
t855 = t602 * t868;
t669 = pkin(5) * t716 + t722 * t862;
t758 = -t666 * t702 + t715 * t865;
t610 = t701 * t669 + t699 * t758;
t711 = legFrame(2,2);
t686 = sin(t711);
t690 = cos(t711);
t603 = t610 * t690 + t645 * t686;
t854 = t603 * t867;
t604 = -t610 * t686 + t645 * t690;
t853 = t604 * t867;
t670 = pkin(5) * t718 + t724 * t861;
t757 = -t667 * t702 + t717 * t865;
t611 = t701 * t670 + t699 * t757;
t712 = legFrame(1,2);
t687 = sin(t712);
t691 = cos(t712);
t605 = t611 * t691 + t646 * t687;
t852 = t605 * t866;
t606 = -t611 * t687 + t646 * t691;
t851 = t606 * t866;
t608 = -t699 * t664 + t701 * t760;
t850 = t608 * t869;
t612 = -t699 * t668 + t701 * t759;
t849 = t612 * t868;
t613 = -t699 * t669 + t701 * t758;
t848 = t613 * t867;
t614 = -t699 * t670 + t701 * t757;
t847 = t614 * t866;
t817 = t702 * t706;
t818 = t702 * t704;
t615 = (-t699 * t704 + t701 * t817) * t864 + pkin(5) * (t699 * t706 + t701 * t818);
t846 = t615 * t736;
t813 = t702 * t720;
t816 = t702 * t714;
t617 = (-t699 * t714 + t701 * t813) * t863 + pkin(5) * (t699 * t720 + t701 * t816);
t845 = t617 * t736;
t812 = t702 * t722;
t815 = t702 * t716;
t618 = (-t699 * t716 + t701 * t812) * t862 + pkin(5) * (t699 * t722 + t701 * t815);
t844 = t618 * t736;
t811 = t702 * t724;
t814 = t702 * t718;
t619 = (-t699 * t718 + t701 * t811) * t861 + pkin(5) * (t699 * t724 + t701 * t814);
t843 = t619 * t736;
t782 = t705 * mrSges(3,1) - t703 * mrSges(3,2);
t631 = t782 * t702 - t704 * t700 * (t703 * mrSges(3,1) + t705 * mrSges(3,2));
t842 = t631 * t869;
t841 = t631 * t736;
t781 = t719 * mrSges(3,1) - t713 * mrSges(3,2);
t632 = t781 * t702 - t714 * t700 * (t713 * mrSges(3,1) + t719 * mrSges(3,2));
t840 = t632 * t868;
t839 = t632 * t736;
t780 = t721 * mrSges(3,1) - t715 * mrSges(3,2);
t633 = t780 * t702 - t716 * t700 * (t715 * mrSges(3,1) + t721 * mrSges(3,2));
t838 = t633 * t867;
t837 = t633 * t736;
t779 = t723 * mrSges(3,1) - t717 * mrSges(3,2);
t634 = t779 * t702 - t718 * t700 * (t717 * mrSges(3,1) + t723 * mrSges(3,2));
t836 = t634 * t866;
t835 = t634 * t736;
t834 = t869 / t705;
t833 = t868 / t719;
t832 = t867 / t721;
t831 = t866 / t723;
t695 = m(1) + m(2) + m(3);
t830 = t869 * t695;
t829 = t868 * t695;
t828 = t867 * t695;
t827 = t866 * t695;
t708 = mrSges(2,2) - mrSges(3,3);
t826 = ((mrSges(2,1) + t782) * t706 - t704 * t708) * t700;
t825 = ((mrSges(2,1) + t781) * t720 - t714 * t708) * t700;
t824 = ((mrSges(2,1) + t780) * t722 - t716 * t708) * t700;
t823 = ((mrSges(2,1) + t779) * t724 - t718 * t708) * t700;
t671 = Ifges(3,5) * t703 + Ifges(3,6) * t705;
t822 = t671 * t736;
t672 = Ifges(3,5) * t713 + Ifges(3,6) * t719;
t821 = t672 * t736;
t673 = Ifges(3,5) * t715 + Ifges(3,6) * t721;
t820 = t673 * t736;
t674 = Ifges(3,5) * t717 + Ifges(3,6) * t723;
t819 = t674 * t736;
t809 = t703 * t706;
t807 = t713 * t720;
t805 = t715 * t722;
t803 = t717 * t724;
t616 = (t699 * t817 + t701 * t704) * t864 + (t699 * t818 - t701 * t706) * pkin(5);
t802 = t616 * t834;
t620 = (t699 * t813 + t701 * t714) * t863 + (t699 * t816 - t701 * t720) * pkin(5);
t801 = t620 * t833;
t621 = (t699 * t812 + t701 * t716) * t862 + (t699 * t815 - t701 * t722) * pkin(5);
t800 = t621 * t832;
t622 = (t699 * t811 + t701 * t718) * t861 + (t699 * t814 - t701 * t724) * pkin(5);
t799 = t622 * t831;
t748 = t700 * t705 + t704 * t810;
t623 = t699 * t748 - t701 * t809;
t798 = t623 * t834;
t747 = t700 * t719 + t714 * t808;
t625 = t699 * t747 - t701 * t807;
t797 = t625 * t833;
t746 = t700 * t721 + t716 * t806;
t626 = t699 * t746 - t701 * t805;
t796 = t626 * t832;
t745 = t700 * t723 + t718 * t804;
t627 = t699 * t745 - t701 * t803;
t795 = t627 * t831;
t794 = t684 * t834;
t793 = t688 * t834;
t792 = t685 * t833;
t791 = t689 * t833;
t790 = t686 * t832;
t789 = t690 * t832;
t788 = t687 * t831;
t787 = t691 * t831;
t786 = t869 * t826;
t785 = t868 * t825;
t784 = t867 * t824;
t783 = t866 * t823;
t725 = xP(4);
t692 = sin(t725);
t693 = cos(t725);
t726 = mrSges(4,2);
t727 = mrSges(4,1);
t778 = -t692 * t726 + t693 * t727;
t777 = t615 * t794;
t776 = t615 * t793;
t775 = t617 * t792;
t774 = t617 * t791;
t773 = t618 * t790;
t772 = t618 * t789;
t771 = t619 * t788;
t770 = t619 * t787;
t624 = t699 * t809 + t701 * t748;
t769 = t624 * t794;
t768 = t624 * t793;
t628 = t699 * t807 + t701 * t747;
t767 = t628 * t792;
t766 = t628 * t791;
t629 = t699 * t805 + t701 * t746;
t765 = t629 * t790;
t764 = t629 * t789;
t630 = t699 * t803 + t701 * t745;
t763 = t630 * t788;
t762 = t630 * t787;
t761 = -t692 * t727 - t693 * t726;
t756 = Ifges(3,3) * t846 + t624 * t671;
t755 = Ifges(3,3) * t845 + t628 * t672;
t754 = Ifges(3,3) * t844 + t629 * t673;
t753 = Ifges(3,3) * t843 + t630 * t674;
t707 = -Ifges(3,1) + Ifges(3,2);
t651 = (t703 * t870 + t707 * t705) * t705 + t860;
t752 = t615 * t822 + t624 * t651;
t652 = (t707 * t719 + t713 * t870) * t719 + t860;
t751 = t617 * t821 + t628 * t652;
t653 = (t707 * t721 + t715 * t870) * t721 + t860;
t750 = t618 * t820 + t629 * t653;
t654 = (t707 * t723 + t717 * t870) * t723 + t860;
t749 = t619 * t819 + t630 * t654;
t744 = t615 * t841 + t624 * t826;
t743 = t617 * t839 + t628 * t825;
t742 = t618 * t837 + t629 * t824;
t741 = t619 * t835 + t630 * t823;
t728 = koppelP(4,2);
t732 = koppelP(4,1);
t655 = -t692 * t732 - t693 * t728;
t659 = -t692 * t728 + t693 * t732;
t740 = (-t655 * t688 + t659 * t684) * t834;
t729 = koppelP(3,2);
t733 = koppelP(3,1);
t656 = -t692 * t733 - t693 * t729;
t660 = -t692 * t729 + t693 * t733;
t739 = (-t656 * t689 + t660 * t685) * t833;
t730 = koppelP(2,2);
t734 = koppelP(2,1);
t657 = -t692 * t734 - t693 * t730;
t661 = -t692 * t730 + t693 * t734;
t738 = (-t657 * t690 + t661 * t686) * t832;
t731 = koppelP(1,2);
t735 = koppelP(1,1);
t658 = -t692 * t735 - t693 * t731;
t662 = -t692 * t731 + t693 * t735;
t737 = (-t658 * t691 + t662 * t687) * t831;
t598 = t630 * t737;
t597 = t629 * t738;
t596 = t628 * t739;
t595 = t624 * t740;
t594 = t737 * t843;
t593 = t738 * t844;
t592 = t739 * t845;
t591 = t740 * t846;
t590 = (t605 * t658 + t606 * t662) * t866;
t589 = (t603 * t657 + t604 * t661) * t867;
t588 = (t601 * t656 + t602 * t660) * t868;
t587 = t614 * t836 + (t622 * t859 + t627 * t674) * t831;
t586 = t613 * t838 + (t621 * t859 + t626 * t673) * t832;
t585 = t612 * t840 + (t620 * t859 + t625 * t672) * t833;
t584 = (t599 * t655 + t600 * t659) * t869;
t583 = t608 * t842 + (t616 * t859 + t623 * t671) * t834;
t582 = t614 * t783 + (t622 * t819 + t627 * t654) * t831;
t581 = t613 * t784 + (t621 * t820 + t626 * t653) * t832;
t580 = t612 * t785 + (t620 * t821 + t625 * t652) * t833;
t579 = t608 * t786 + (t616 * t822 + t623 * t651) * t834;
t578 = t614 * t827 + (t622 * t835 + t627 * t823) * t831;
t577 = t613 * t828 + (t621 * t837 + t626 * t824) * t832;
t576 = t612 * t829 + (t620 * t839 + t625 * t825) * t833;
t575 = t608 * t830 + (t616 * t841 + t623 * t826) * t834;
t574 = t606 * t836 + t753 * t788;
t573 = t605 * t836 - t753 * t787;
t572 = t604 * t838 + t754 * t790;
t571 = t603 * t838 - t754 * t789;
t570 = t602 * t840 + t755 * t792;
t569 = t601 * t840 - t755 * t791;
t568 = t600 * t842 + t756 * t794;
t567 = t599 * t842 - t756 * t793;
t566 = t606 * t783 + t749 * t788;
t565 = t605 * t783 - t749 * t787;
t564 = t604 * t784 + t750 * t790;
t563 = t603 * t784 - t750 * t789;
t562 = t602 * t785 + t751 * t792;
t561 = t601 * t785 - t751 * t791;
t560 = t606 * t827 + t741 * t788;
t559 = t605 * t827 - t741 * t787;
t558 = t604 * t828 + t742 * t790;
t557 = t603 * t828 - t742 * t789;
t556 = t602 * t829 + t743 * t792;
t555 = t601 * t829 - t743 * t791;
t554 = t600 * t786 + t752 * t794;
t553 = t599 * t786 - t752 * t793;
t552 = t600 * t830 + t744 * t794;
t551 = t599 * t830 - t744 * t793;
t550 = t594 * Ifges(3,3) + t590 * t634 + t598 * t674;
t549 = t593 * Ifges(3,3) + t589 * t633 + t597 * t673;
t548 = t592 * Ifges(3,3) + t588 * t632 + t596 * t672;
t547 = t590 * t823 + t594 * t674 + t598 * t654;
t546 = t589 * t824 + t593 * t673 + t597 * t653;
t545 = t588 * t825 + t592 * t672 + t596 * t652;
t544 = t591 * Ifges(3,3) + t584 * t631 + t595 * t671;
t543 = t590 * t695 + t594 * t634 + t598 * t823;
t542 = t589 * t695 + t593 * t633 + t597 * t824;
t541 = t588 * t695 + t592 * t632 + t596 * t825;
t540 = t584 * t826 + t591 * t671 + t595 * t651;
t539 = t584 * t695 + t591 * t631 + t595 * t826;
t1 = [t559 * t852 - t565 * t762 + t557 * t854 - t563 * t764 + t555 * t856 - t561 * t766 + t551 * t858 - t553 * t768 + m(4) + (-t567 * t776 - t569 * t774 - t571 * t772 - t573 * t770) * t736, t559 * t851 + t565 * t763 + t557 * t853 + t563 * t765 + t555 * t855 + t561 * t767 + t551 * t857 + t553 * t769 + (t567 * t777 + t569 * t775 + t571 * t773 + t573 * t771) * t736, t553 * t798 + t561 * t797 + t563 * t796 + t565 * t795 + t551 * t850 + t555 * t849 + t557 * t848 + t559 * t847 + (t567 * t802 + t569 * t801 + t571 * t800 + t573 * t799) * t736, t551 * t584 + t553 * t595 + t555 * t588 + t557 * t589 + t559 * t590 + t561 * t596 + t563 * t597 + t565 * t598 + t567 * t591 + t569 * t592 + t571 * t593 + t573 * t594 + t761; t560 * t852 - t566 * t762 + t558 * t854 - t564 * t764 + t556 * t856 - t562 * t766 + t552 * t858 - t554 * t768 + (-t568 * t776 - t570 * t774 - t572 * t772 - t574 * t770) * t736, t560 * t851 + t566 * t763 + t558 * t853 + t564 * t765 + t556 * t855 + t562 * t767 + t552 * t857 + t554 * t769 + m(4) + (t568 * t777 + t570 * t775 + t572 * t773 + t574 * t771) * t736, t554 * t798 + t562 * t797 + t564 * t796 + t566 * t795 + t552 * t850 + t556 * t849 + t558 * t848 + t560 * t847 + (t568 * t802 + t570 * t801 + t572 * t800 + t574 * t799) * t736, t552 * t584 + t554 * t595 + t556 * t588 + t558 * t589 + t560 * t590 + t562 * t596 + t564 * t597 + t566 * t598 + t568 * t591 + t570 * t592 + t572 * t593 + t574 * t594 + t778; t578 * t852 - t582 * t762 + t577 * t854 - t581 * t764 + t576 * t856 - t580 * t766 + t575 * t858 - t579 * t768 + (-t583 * t776 - t585 * t774 - t586 * t772 - t587 * t770) * t736, t578 * t851 + t582 * t763 + t577 * t853 + t581 * t765 + t576 * t855 + t580 * t767 + t575 * t857 + t579 * t769 + (t583 * t777 + t585 * t775 + t586 * t773 + t587 * t771) * t736, t579 * t798 + t580 * t797 + t581 * t796 + t582 * t795 + t575 * t850 + t576 * t849 + t577 * t848 + t578 * t847 + m(4) + (t583 * t802 + t585 * t801 + t586 * t800 + t587 * t799) * t736, t575 * t584 + t576 * t588 + t577 * t589 + t578 * t590 + t579 * t595 + t580 * t596 + t581 * t597 + t582 * t598 + t583 * t591 + t585 * t592 + t586 * t593 + t587 * t594; t543 * t852 - t547 * t762 + t542 * t854 - t546 * t764 + t541 * t856 - t545 * t766 + t539 * t858 - t540 * t768 + (-t544 * t776 - t548 * t774 - t549 * t772 - t550 * t770) * t736 + t761, t543 * t851 + t547 * t763 + t542 * t853 + t546 * t765 + t541 * t855 + t545 * t767 + t539 * t857 + t540 * t769 + (t544 * t777 + t548 * t775 + t549 * t773 + t550 * t771) * t736 + t778, t540 * t798 + t545 * t797 + t546 * t796 + t547 * t795 + t539 * t850 + t541 * t849 + t542 * t848 + t543 * t847 + (t544 * t802 + t548 * t801 + t549 * t800 + t550 * t799) * t736, t539 * t584 + t540 * t595 + t541 * t588 + t542 * t589 + t543 * t590 + t544 * t591 + t545 * t596 + t546 * t597 + t547 * t598 + t548 * t592 + t549 * t593 + t550 * t594 + Ifges(4,3);];
MX  = t1;
