% Calculate inertia matrix for parallel robot
% P4PRRRR8V1G2A0
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
% Datum: 2020-08-07 11:06
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P4PRRRR8V1G2A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V1G2A0_inertia_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V1G2A0_inertia_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P4PRRRR8V1G2A0_inertia_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V1G2A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRRR8V1G2A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P4PRRRR8V1G2A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V1G2A0_inertia_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V1G2A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:03:01
% EndTime: 2020-08-07 11:03:04
% DurationCPUTime: 2.92s
% Computational Cost: add. (8124->390), mult. (21960->699), div. (1920->9), fcn. (22920->30), ass. (0->308)
t872 = 2 * Ifges(3,4);
t706 = sin(qJ(2,4));
t708 = cos(qJ(2,4));
t707 = cos(qJ(3,4));
t866 = pkin(2) * t707;
t665 = -t708 * pkin(5) + t706 * t866;
t702 = sin(pkin(3));
t704 = cos(pkin(3));
t705 = sin(qJ(3,4));
t812 = t705 * t704;
t642 = pkin(2) * t812 + t665 * t702;
t871 = 0.1e1 / t642;
t716 = sin(qJ(2,3));
t722 = cos(qJ(2,3));
t721 = cos(qJ(3,3));
t865 = pkin(2) * t721;
t667 = -t722 * pkin(5) + t716 * t865;
t715 = sin(qJ(3,3));
t810 = t715 * t704;
t646 = pkin(2) * t810 + t667 * t702;
t870 = 0.1e1 / t646;
t718 = sin(qJ(2,2));
t724 = cos(qJ(2,2));
t723 = cos(qJ(3,2));
t864 = pkin(2) * t723;
t668 = -t724 * pkin(5) + t718 * t864;
t717 = sin(qJ(3,2));
t808 = t717 * t704;
t647 = pkin(2) * t808 + t668 * t702;
t869 = 0.1e1 / t647;
t720 = sin(qJ(2,1));
t726 = cos(qJ(2,1));
t725 = cos(qJ(3,1));
t863 = pkin(2) * t725;
t669 = -t726 * pkin(5) + t720 * t863;
t719 = sin(qJ(3,1));
t806 = t719 * t704;
t648 = pkin(2) * t806 + t669 * t702;
t868 = 0.1e1 / t648;
t867 = pkin(2) * t702;
t862 = Ifges(3,1) + Ifges(2,3);
t738 = 0.1e1 / pkin(2);
t861 = Ifges(3,3) * t738;
t666 = pkin(5) * t706 + t708 * t866;
t701 = sin(pkin(6));
t703 = cos(pkin(6));
t762 = -t665 * t704 + t705 * t867;
t609 = -t701 * t666 + t762 * t703;
t711 = legFrame(4,2);
t686 = sin(t711);
t690 = cos(t711);
t601 = t609 * t686 + t690 * t642;
t860 = t601 * t871;
t602 = -t609 * t690 + t686 * t642;
t859 = t602 * t871;
t670 = pkin(5) * t716 + t722 * t865;
t761 = -t667 * t704 + t715 * t867;
t611 = -t701 * t670 + t761 * t703;
t712 = legFrame(3,2);
t687 = sin(t712);
t691 = cos(t712);
t603 = t611 * t687 + t691 * t646;
t858 = t603 * t870;
t604 = -t611 * t691 + t687 * t646;
t857 = t604 * t870;
t671 = pkin(5) * t718 + t724 * t864;
t760 = -t668 * t704 + t717 * t867;
t612 = -t701 * t671 + t760 * t703;
t713 = legFrame(2,2);
t688 = sin(t713);
t692 = cos(t713);
t605 = t612 * t688 + t692 * t647;
t856 = t605 * t869;
t606 = -t612 * t692 + t688 * t647;
t855 = t606 * t869;
t672 = pkin(5) * t720 + t726 * t863;
t759 = -t669 * t704 + t719 * t867;
t613 = -t701 * t672 + t759 * t703;
t714 = legFrame(1,2);
t689 = sin(t714);
t693 = cos(t714);
t607 = t613 * t689 + t693 * t648;
t854 = t607 * t868;
t608 = -t613 * t693 + t689 * t648;
t853 = t608 * t868;
t610 = t703 * t666 + t762 * t701;
t852 = t610 * t871;
t614 = t703 * t670 + t761 * t701;
t851 = t614 * t870;
t615 = t703 * t671 + t760 * t701;
t850 = t615 * t869;
t616 = t703 * t672 + t759 * t701;
t849 = t616 * t868;
t819 = t704 * t708;
t820 = t704 * t706;
t618 = (t701 * t819 + t703 * t706) * t866 + (t701 * t820 - t703 * t708) * pkin(5);
t848 = t618 * t738;
t815 = t704 * t722;
t818 = t704 * t716;
t622 = (t701 * t815 + t703 * t716) * t865 + (t701 * t818 - t703 * t722) * pkin(5);
t847 = t622 * t738;
t814 = t704 * t724;
t817 = t704 * t718;
t623 = (t701 * t814 + t703 * t718) * t864 + (t701 * t817 - t703 * t724) * pkin(5);
t846 = t623 * t738;
t813 = t704 * t726;
t816 = t704 * t720;
t624 = (t701 * t813 + t703 * t720) * t863 + (t701 * t816 - t703 * t726) * pkin(5);
t845 = t624 * t738;
t784 = t707 * mrSges(3,1) - t705 * mrSges(3,2);
t633 = t784 * t704 - t702 * (t705 * mrSges(3,1) + t707 * mrSges(3,2)) * t706;
t844 = t633 * t871;
t843 = t633 * t738;
t783 = t721 * mrSges(3,1) - t715 * mrSges(3,2);
t634 = t783 * t704 - t702 * (t715 * mrSges(3,1) + t721 * mrSges(3,2)) * t716;
t842 = t634 * t870;
t841 = t634 * t738;
t782 = t723 * mrSges(3,1) - t717 * mrSges(3,2);
t635 = t782 * t704 - t702 * (t717 * mrSges(3,1) + t723 * mrSges(3,2)) * t718;
t840 = t635 * t869;
t839 = t635 * t738;
t781 = t725 * mrSges(3,1) - t719 * mrSges(3,2);
t636 = t781 * t704 - t702 * (t719 * mrSges(3,1) + t725 * mrSges(3,2)) * t720;
t838 = t636 * t868;
t837 = t636 * t738;
t836 = t871 / t707;
t835 = t870 / t721;
t834 = t869 / t723;
t833 = t868 / t725;
t697 = m(1) + m(2) + m(3);
t832 = t871 * t697;
t831 = t870 * t697;
t830 = t869 * t697;
t829 = t868 * t697;
t710 = mrSges(2,2) - mrSges(3,3);
t828 = ((mrSges(2,1) + t784) * t708 - t706 * t710) * t702;
t827 = ((mrSges(2,1) + t783) * t722 - t716 * t710) * t702;
t826 = ((mrSges(2,1) + t782) * t724 - t718 * t710) * t702;
t825 = ((mrSges(2,1) + t781) * t726 - t720 * t710) * t702;
t673 = Ifges(3,5) * t705 + Ifges(3,6) * t707;
t824 = t673 * t738;
t674 = Ifges(3,5) * t715 + Ifges(3,6) * t721;
t823 = t674 * t738;
t675 = Ifges(3,5) * t717 + Ifges(3,6) * t723;
t822 = t675 * t738;
t676 = Ifges(3,5) * t719 + Ifges(3,6) * t725;
t821 = t676 * t738;
t811 = t705 * t708;
t809 = t715 * t722;
t807 = t717 * t724;
t805 = t719 * t726;
t617 = -(-t701 * t706 + t703 * t819) * t866 - pkin(5) * (t701 * t708 + t703 * t820);
t804 = t617 * t836;
t619 = -(-t701 * t716 + t703 * t815) * t865 - pkin(5) * (t701 * t722 + t703 * t818);
t803 = t619 * t835;
t620 = -(-t701 * t718 + t703 * t814) * t864 - pkin(5) * (t701 * t724 + t703 * t817);
t802 = t620 * t834;
t621 = -(-t701 * t720 + t703 * t813) * t863 - pkin(5) * (t701 * t726 + t703 * t816);
t801 = t621 * t833;
t750 = t702 * t707 + t706 * t812;
t626 = -t701 * t811 - t750 * t703;
t800 = t626 * t836;
t749 = t702 * t721 + t716 * t810;
t630 = -t701 * t809 - t749 * t703;
t799 = t630 * t835;
t748 = t702 * t723 + t718 * t808;
t631 = -t701 * t807 - t748 * t703;
t798 = t631 * t834;
t747 = t702 * t725 + t720 * t806;
t632 = -t701 * t805 - t747 * t703;
t797 = t632 * t833;
t796 = t686 * t836;
t795 = t690 * t836;
t794 = t687 * t835;
t793 = t691 * t835;
t792 = t688 * t834;
t791 = t692 * t834;
t790 = t689 * t833;
t789 = t693 * t833;
t788 = t871 * t828;
t787 = t870 * t827;
t786 = t869 * t826;
t785 = t868 * t825;
t727 = xP(4);
t694 = sin(t727);
t695 = cos(t727);
t728 = mrSges(4,2);
t729 = mrSges(4,1);
t780 = -t694 * t728 + t695 * t729;
t779 = t618 * t796;
t778 = t618 * t795;
t777 = t622 * t794;
t776 = t622 * t793;
t775 = t623 * t792;
t774 = t623 * t791;
t773 = t624 * t790;
t772 = t624 * t789;
t625 = t750 * t701 - t703 * t811;
t771 = t625 * t796;
t770 = t625 * t795;
t627 = t749 * t701 - t703 * t809;
t769 = t627 * t794;
t768 = t627 * t793;
t628 = t748 * t701 - t703 * t807;
t767 = t628 * t792;
t766 = t628 * t791;
t629 = t747 * t701 - t703 * t805;
t765 = t629 * t790;
t764 = t629 * t789;
t763 = -t694 * t729 - t695 * t728;
t758 = Ifges(3,3) * t848 + t625 * t673;
t757 = Ifges(3,3) * t847 + t627 * t674;
t756 = Ifges(3,3) * t846 + t628 * t675;
t755 = Ifges(3,3) * t845 + t629 * t676;
t709 = -Ifges(3,1) + Ifges(3,2);
t653 = (t705 * t872 + t709 * t707) * t707 + t862;
t754 = t618 * t824 + t625 * t653;
t654 = (t709 * t721 + t715 * t872) * t721 + t862;
t753 = t622 * t823 + t627 * t654;
t655 = (t709 * t723 + t717 * t872) * t723 + t862;
t752 = t623 * t822 + t628 * t655;
t656 = (t709 * t725 + t719 * t872) * t725 + t862;
t751 = t624 * t821 + t629 * t656;
t746 = t618 * t843 + t625 * t828;
t745 = t622 * t841 + t627 * t827;
t744 = t623 * t839 + t628 * t826;
t743 = t624 * t837 + t629 * t825;
t730 = koppelP(4,2);
t734 = koppelP(4,1);
t657 = -t694 * t734 - t695 * t730;
t661 = -t694 * t730 + t695 * t734;
t742 = (-t657 * t690 + t661 * t686) * t836;
t731 = koppelP(3,2);
t735 = koppelP(3,1);
t658 = -t694 * t735 - t695 * t731;
t662 = -t694 * t731 + t695 * t735;
t741 = (-t658 * t691 + t662 * t687) * t835;
t732 = koppelP(2,2);
t736 = koppelP(2,1);
t659 = -t694 * t736 - t695 * t732;
t663 = -t694 * t732 + t695 * t736;
t740 = (-t659 * t692 + t663 * t688) * t834;
t733 = koppelP(1,2);
t737 = koppelP(1,1);
t660 = -t694 * t737 - t695 * t733;
t664 = -t694 * t733 + t695 * t737;
t739 = (-t660 * t693 + t664 * t689) * t833;
t600 = t629 * t739;
t599 = t628 * t740;
t598 = t627 * t741;
t597 = t625 * t742;
t596 = t739 * t845;
t595 = t740 * t846;
t594 = t741 * t847;
t593 = t742 * t848;
t592 = (t607 * t664 + t608 * t660) * t868;
t591 = (t605 * t663 + t606 * t659) * t869;
t590 = (t603 * t662 + t604 * t658) * t870;
t589 = t616 * t838 + (t621 * t861 + t632 * t676) * t833;
t588 = t615 * t840 + (t620 * t861 + t631 * t675) * t834;
t587 = t614 * t842 + (t619 * t861 + t630 * t674) * t835;
t586 = (t601 * t661 + t602 * t657) * t871;
t585 = t610 * t844 + (t617 * t861 + t626 * t673) * t836;
t584 = t616 * t785 + (t621 * t821 + t632 * t656) * t833;
t583 = t615 * t786 + (t620 * t822 + t631 * t655) * t834;
t582 = t614 * t787 + (t619 * t823 + t630 * t654) * t835;
t581 = t610 * t788 + (t617 * t824 + t626 * t653) * t836;
t580 = t616 * t829 + (t621 * t837 + t632 * t825) * t833;
t579 = t615 * t830 + (t620 * t839 + t631 * t826) * t834;
t578 = t614 * t831 + (t619 * t841 + t630 * t827) * t835;
t577 = t610 * t832 + (t617 * t843 + t626 * t828) * t836;
t576 = t608 * t838 - t755 * t789;
t575 = t607 * t838 + t755 * t790;
t574 = t606 * t840 - t756 * t791;
t573 = t605 * t840 + t756 * t792;
t572 = t604 * t842 - t757 * t793;
t571 = t603 * t842 + t757 * t794;
t570 = t602 * t844 - t758 * t795;
t569 = t601 * t844 + t758 * t796;
t568 = t608 * t785 - t751 * t789;
t567 = t607 * t785 + t751 * t790;
t566 = t606 * t786 - t752 * t791;
t565 = t605 * t786 + t752 * t792;
t564 = t604 * t787 - t753 * t793;
t563 = t603 * t787 + t753 * t794;
t562 = t608 * t829 - t743 * t789;
t561 = t607 * t829 + t743 * t790;
t560 = t606 * t830 - t744 * t791;
t559 = t605 * t830 + t744 * t792;
t558 = t604 * t831 - t745 * t793;
t557 = t603 * t831 + t745 * t794;
t556 = t602 * t788 - t754 * t795;
t555 = t601 * t788 + t754 * t796;
t554 = t602 * t832 - t746 * t795;
t553 = t601 * t832 + t746 * t796;
t552 = t596 * Ifges(3,3) + t592 * t636 + t600 * t676;
t551 = t595 * Ifges(3,3) + t591 * t635 + t599 * t675;
t550 = t594 * Ifges(3,3) + t590 * t634 + t598 * t674;
t549 = t592 * t825 + t596 * t676 + t600 * t656;
t548 = t591 * t826 + t595 * t675 + t599 * t655;
t547 = t590 * t827 + t594 * t674 + t598 * t654;
t546 = t593 * Ifges(3,3) + t586 * t633 + t597 * t673;
t545 = t592 * t697 + t596 * t636 + t600 * t825;
t544 = t591 * t697 + t595 * t635 + t599 * t826;
t543 = t590 * t697 + t594 * t634 + t598 * t827;
t542 = t586 * t828 + t593 * t673 + t597 * t653;
t541 = t586 * t697 + t593 * t633 + t597 * t828;
t1 = [t562 * t853 - t568 * t764 + t560 * t855 - t566 * t766 + t558 * t857 - t564 * t768 + t554 * t859 - t556 * t770 + m(4) + (-t570 * t778 - t572 * t776 - t574 * t774 - t576 * t772) * t738, t562 * t854 + t568 * t765 + t560 * t856 + t566 * t767 + t558 * t858 + t564 * t769 + t554 * t860 + t556 * t771 + (t570 * t779 + t572 * t777 + t574 * t775 + t576 * t773) * t738, t556 * t800 + t564 * t799 + t566 * t798 + t568 * t797 + t554 * t852 + t558 * t851 + t560 * t850 + t562 * t849 + (t570 * t804 + t572 * t803 + t574 * t802 + t576 * t801) * t738, t554 * t586 + t556 * t597 + t558 * t590 + t560 * t591 + t562 * t592 + t564 * t598 + t566 * t599 + t568 * t600 + t570 * t593 + t572 * t594 + t574 * t595 + t576 * t596 + t763; t561 * t853 - t567 * t764 + t559 * t855 - t565 * t766 + t557 * t857 - t563 * t768 + t553 * t859 - t555 * t770 + (-t569 * t778 - t571 * t776 - t573 * t774 - t575 * t772) * t738, t561 * t854 + t567 * t765 + t559 * t856 + t565 * t767 + t557 * t858 + t563 * t769 + t553 * t860 + t555 * t771 + m(4) + (t569 * t779 + t571 * t777 + t573 * t775 + t575 * t773) * t738, t555 * t800 + t563 * t799 + t565 * t798 + t567 * t797 + t553 * t852 + t557 * t851 + t559 * t850 + t561 * t849 + (t569 * t804 + t571 * t803 + t573 * t802 + t575 * t801) * t738, t553 * t586 + t555 * t597 + t557 * t590 + t559 * t591 + t561 * t592 + t563 * t598 + t565 * t599 + t567 * t600 + t569 * t593 + t571 * t594 + t573 * t595 + t575 * t596 + t780; t580 * t853 - t584 * t764 + t579 * t855 - t583 * t766 + t578 * t857 - t582 * t768 + t577 * t859 - t581 * t770 + (-t585 * t778 - t587 * t776 - t588 * t774 - t589 * t772) * t738, t580 * t854 + t584 * t765 + t579 * t856 + t583 * t767 + t578 * t858 + t582 * t769 + t577 * t860 + t581 * t771 + (t585 * t779 + t587 * t777 + t588 * t775 + t589 * t773) * t738, t581 * t800 + t582 * t799 + t583 * t798 + t584 * t797 + t577 * t852 + t578 * t851 + t579 * t850 + t580 * t849 + m(4) + (t585 * t804 + t587 * t803 + t588 * t802 + t589 * t801) * t738, t577 * t586 + t578 * t590 + t579 * t591 + t580 * t592 + t581 * t597 + t582 * t598 + t583 * t599 + t584 * t600 + t585 * t593 + t587 * t594 + t588 * t595 + t589 * t596; t545 * t853 - t549 * t764 + t544 * t855 - t548 * t766 + t543 * t857 - t547 * t768 + t541 * t859 - t542 * t770 + (-t546 * t778 - t550 * t776 - t551 * t774 - t552 * t772) * t738 + t763, t545 * t854 + t549 * t765 + t544 * t856 + t548 * t767 + t543 * t858 + t547 * t769 + t541 * t860 + t542 * t771 + (t546 * t779 + t550 * t777 + t551 * t775 + t552 * t773) * t738 + t780, t542 * t800 + t547 * t799 + t548 * t798 + t549 * t797 + t541 * t852 + t543 * t851 + t544 * t850 + t545 * t849 + (t546 * t804 + t550 * t803 + t551 * t802 + t552 * t801) * t738, t541 * t586 + t542 * t597 + t543 * t590 + t544 * t591 + t545 * t592 + t546 * t593 + t547 * t598 + t548 * t599 + t549 * t600 + t550 * t594 + t551 * t595 + t552 * t596 + Ifges(4,3);];
MX  = t1;
