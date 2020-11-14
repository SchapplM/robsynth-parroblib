% Calculate inertia matrix for parallel robot
% P4RRRRR2G1A0
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
% Datum: 2020-08-07 17:26
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P4RRRRR2G1A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4RRRRR2G1A0_inertia_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4RRRRR2G1A0_inertia_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4RRRRR2G1A0_inertia_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4RRRRR2G1A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4RRRRR2G1A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P4RRRRR2G1A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4RRRRR2G1A0_inertia_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4RRRRR2G1A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 17:24:37
% EndTime: 2020-08-07 17:24:39
% DurationCPUTime: 1.63s
% Computational Cost: add. (5900->328), mult. (5288->530), div. (1508->18), fcn. (4076->58), ass. (0->264)
t863 = -2 * pkin(1);
t725 = cos(qJ(3,4));
t708 = 0.1e1 / t725;
t862 = t708 ^ 2;
t735 = cos(qJ(3,3));
t713 = 0.1e1 / t735;
t861 = t713 ^ 2;
t737 = cos(qJ(3,2));
t715 = 0.1e1 / t737;
t860 = t715 ^ 2;
t739 = cos(qJ(3,1));
t717 = 0.1e1 / t739;
t859 = t717 ^ 2;
t724 = sin(qJ(2,4));
t858 = pkin(1) * t724;
t730 = sin(qJ(2,3));
t857 = pkin(1) * t730;
t732 = sin(qJ(2,2));
t856 = pkin(1) * t732;
t734 = sin(qJ(2,1));
t855 = pkin(1) * t734;
t726 = cos(qJ(2,4));
t854 = t726 * pkin(1);
t736 = cos(qJ(2,3));
t853 = t736 * pkin(1);
t738 = cos(qJ(2,2));
t852 = t738 * pkin(1);
t740 = cos(qJ(2,1));
t851 = t740 * pkin(1);
t850 = Ifges(3,1) + Ifges(2,3);
t723 = sin(qJ(3,4));
t849 = Ifges(3,4) * t723;
t729 = sin(qJ(3,3));
t848 = Ifges(3,4) * t729;
t731 = sin(qJ(3,2));
t847 = Ifges(3,4) * t731;
t733 = sin(qJ(3,1));
t846 = Ifges(3,4) * t733;
t700 = qJ(1,4) + legFrame(4,3);
t684 = qJ(2,4) + t700;
t669 = qJ(3,4) + t684;
t670 = -qJ(3,4) + t684;
t631 = sin(t700) * t863 + (-sin(t670) - sin(t669)) * pkin(2);
t655 = 0.1e1 / (sin(qJ(2,4) + qJ(3,4)) + sin(qJ(2,4) - qJ(3,4)));
t845 = t631 * t655;
t632 = cos(t700) * t863 + (-cos(t669) - cos(t670)) * pkin(2);
t844 = t632 * t655;
t701 = qJ(1,3) + legFrame(3,3);
t687 = qJ(2,3) + t701;
t677 = qJ(3,3) + t687;
t678 = -qJ(3,3) + t687;
t633 = sin(t701) * t863 + (-sin(t678) - sin(t677)) * pkin(2);
t656 = 0.1e1 / (sin(qJ(2,3) + qJ(3,3)) + sin(qJ(2,3) - qJ(3,3)));
t843 = t633 * t656;
t702 = qJ(1,2) + legFrame(2,3);
t688 = qJ(2,2) + t702;
t679 = qJ(3,2) + t688;
t680 = -qJ(3,2) + t688;
t634 = sin(t702) * t863 + (-sin(t680) - sin(t679)) * pkin(2);
t657 = 0.1e1 / (sin(qJ(2,2) + qJ(3,2)) + sin(qJ(2,2) - qJ(3,2)));
t842 = t634 * t657;
t703 = qJ(1,1) + legFrame(1,3);
t689 = qJ(2,1) + t703;
t681 = qJ(3,1) + t689;
t682 = -qJ(3,1) + t689;
t635 = sin(t703) * t863 + (-sin(t682) - sin(t681)) * pkin(2);
t658 = 0.1e1 / (sin(qJ(2,1) + qJ(3,1)) + sin(qJ(2,1) - qJ(3,1)));
t841 = t635 * t658;
t636 = cos(t701) * t863 + (-cos(t677) - cos(t678)) * pkin(2);
t840 = t636 * t656;
t637 = cos(t702) * t863 + (-cos(t679) - cos(t680)) * pkin(2);
t839 = t637 * t657;
t638 = cos(t703) * t863 + (-cos(t681) - cos(t682)) * pkin(2);
t838 = t638 * t658;
t639 = t725 * (-mrSges(3,2) * t858 + Ifges(3,6)) - t723 * (mrSges(3,1) * t858 - Ifges(3,5));
t837 = t639 * t708;
t640 = t735 * (-mrSges(3,2) * t857 + Ifges(3,6)) - t729 * (mrSges(3,1) * t857 - Ifges(3,5));
t836 = t640 * t713;
t641 = t737 * (-mrSges(3,2) * t856 + Ifges(3,6)) - t731 * (mrSges(3,1) * t856 - Ifges(3,5));
t835 = t641 * t715;
t642 = t739 * (-mrSges(3,2) * t855 + Ifges(3,6)) - t733 * (mrSges(3,1) * t855 - Ifges(3,5));
t834 = t642 * t717;
t752 = 0.1e1 / pkin(2);
t833 = t655 * t752;
t832 = t656 * t752;
t831 = t657 * t752;
t830 = t658 * t752;
t754 = t725 ^ 2;
t829 = (t725 * pkin(2) + t854) / t754;
t755 = t735 ^ 2;
t828 = (t735 * pkin(2) + t853) / t755;
t756 = t737 ^ 2;
t827 = (t737 * pkin(2) + t852) / t756;
t757 = t739 ^ 2;
t826 = (t739 * pkin(2) + t851) / t757;
t667 = sin(t684);
t707 = 0.1e1 / t724;
t825 = t667 * t707;
t668 = cos(t684);
t824 = t668 * t707;
t671 = sin(t687);
t710 = 0.1e1 / t730;
t823 = t671 * t710;
t672 = sin(t688);
t711 = 0.1e1 / t732;
t822 = t672 * t711;
t673 = sin(t689);
t712 = 0.1e1 / t734;
t821 = t673 * t712;
t674 = cos(t687);
t820 = t674 * t710;
t675 = cos(t688);
t819 = t675 * t711;
t676 = cos(t689);
t818 = t676 * t712;
t817 = t707 * t723;
t753 = 1 / pkin(1);
t816 = t707 * t753;
t815 = t708 * t752;
t814 = t710 * t729;
t813 = t710 * t753;
t812 = t711 * t731;
t811 = t711 * t753;
t810 = t712 * t733;
t809 = t712 * t753;
t808 = t713 * t752;
t807 = t715 * t752;
t806 = t717 * t752;
t805 = t752 * t753;
t804 = 0.2e1 * t849;
t803 = 0.2e1 * t848;
t802 = 0.2e1 * t847;
t801 = 0.2e1 * t846;
t727 = Ifges(3,2) - Ifges(3,1);
t686 = t727 * t754;
t800 = t686 + t850;
t693 = t727 * t755;
t799 = t693 + t850;
t694 = t727 * t756;
t798 = t694 + t850;
t695 = t727 * t757;
t797 = t695 + t850;
t696 = mrSges(3,1) * t854;
t728 = mrSges(2,2) - mrSges(3,3);
t761 = (-(t723 * mrSges(3,2) - mrSges(2,1)) * t726 - t728 * t724) * pkin(1);
t627 = (t696 + t804) * t725 + t761 + t800;
t796 = t627 * t833;
t697 = mrSges(3,1) * t853;
t760 = (-(t729 * mrSges(3,2) - mrSges(2,1)) * t736 - t728 * t730) * pkin(1);
t628 = (t697 + t803) * t735 + t760 + t799;
t795 = t628 * t832;
t698 = mrSges(3,1) * t852;
t759 = (-(t731 * mrSges(3,2) - mrSges(2,1)) * t738 - t728 * t732) * pkin(1);
t629 = (t698 + t802) * t737 + t759 + t798;
t794 = t629 * t831;
t699 = mrSges(3,1) * t851;
t758 = (-(t733 * mrSges(3,2) - mrSges(2,1)) * t740 - t728 * t734) * pkin(1);
t630 = (t699 + t801) * t739 + t758 + t797;
t793 = t630 * t830;
t643 = t725 * t804 + t800;
t792 = t643 * t833;
t644 = t735 * t803 + t799;
t791 = t644 * t832;
t645 = t737 * t802 + t798;
t790 = t645 * t831;
t646 = t739 * t801 + t797;
t789 = t646 * t830;
t659 = Ifges(3,5) * t723 + Ifges(3,6) * t725;
t788 = t655 * t659 * t708;
t661 = Ifges(3,5) * t729 + Ifges(3,6) * t735;
t787 = t656 * t661 * t713;
t662 = Ifges(3,5) * t731 + Ifges(3,6) * t737;
t786 = t657 * t662 * t715;
t663 = Ifges(3,5) * t733 + Ifges(3,6) * t739;
t785 = t658 * t663 * t717;
t784 = t723 * t829;
t783 = t752 * t829;
t782 = t729 * t828;
t781 = t752 * t828;
t780 = t731 * t827;
t779 = t752 * t827;
t778 = t733 * t826;
t777 = t752 * t826;
t776 = t708 * t817;
t775 = t723 * t816;
t774 = t659 * t815;
t773 = t713 * t814;
t772 = t729 * t813;
t771 = t715 * t812;
t770 = t731 * t811;
t769 = t717 * t810;
t768 = t733 * t809;
t767 = t661 * t808;
t766 = t662 * t807;
t765 = t663 * t806;
t764 = Ifges(1,3) + ((m(2) + m(3)) * pkin(1) ^ 2) + t850;
t741 = xP(4);
t705 = sin(t741);
t706 = cos(t741);
t742 = mrSges(4,2);
t743 = mrSges(4,1);
t763 = -t705 * t742 + t706 * t743;
t762 = -t705 * t743 - t706 * t742;
t751 = koppelP(1,1);
t750 = koppelP(2,1);
t749 = koppelP(3,1);
t748 = koppelP(4,1);
t747 = koppelP(1,2);
t746 = koppelP(2,2);
t745 = koppelP(3,2);
t744 = koppelP(4,2);
t654 = -t705 * t747 + t706 * t751;
t653 = -t705 * t746 + t706 * t750;
t652 = -t705 * t745 + t706 * t749;
t651 = -t705 * t744 + t706 * t748;
t650 = -t705 * t751 - t706 * t747;
t649 = -t705 * t750 - t706 * t746;
t648 = -t705 * t749 - t706 * t745;
t647 = -t705 * t748 - t706 * t744;
t626 = t695 + 0.2e1 * (t699 + t846) * t739 + 0.2e1 * t758 + t764;
t625 = t694 + 0.2e1 * (t698 + t847) * t737 + 0.2e1 * t759 + t764;
t624 = t693 + 0.2e1 * (t697 + t848) * t735 + 0.2e1 * t760 + t764;
t623 = t686 + 0.2e1 * (t696 + t849) * t725 + 0.2e1 * t761 + t764;
t622 = (t650 * t676 + t654 * t673) * t809;
t621 = (t649 * t675 + t653 * t672) * t811;
t620 = (t648 * t674 + t652 * t671) * t813;
t619 = (t647 * t668 + t651 * t667) * t816;
t618 = t765 + (t630 * t717 - t646 * t777) * t768;
t617 = t766 + (t629 * t715 - t645 * t779) * t770;
t616 = t767 + (t628 * t713 - t644 * t781) * t772;
t615 = t774 + (t627 * t708 - t643 * t783) * t775;
t614 = (t630 * t818 + t638 * t789) * t753;
t613 = (t629 * t819 + t637 * t790) * t753;
t612 = (t628 * t820 + t636 * t791) * t753;
t611 = (t630 * t821 + t635 * t789) * t753;
t610 = (t629 * t822 + t634 * t790) * t753;
t609 = (t628 * t823 + t633 * t791) * t753;
t608 = (t635 * t654 + t638 * t650) * t658 * t805;
t607 = (t634 * t653 + t637 * t649) * t657 * t805;
t606 = (t633 * t652 + t636 * t648) * t656 * t805;
t605 = (t627 * t824 + t632 * t792) * t753;
t604 = (t627 * t825 + t631 * t792) * t753;
t603 = (t631 * t651 + t632 * t647) * t655 * t805;
t602 = (t626 * t818 + t638 * t793) * t753;
t601 = (t625 * t819 + t637 * t794) * t753;
t600 = (t624 * t820 + t636 * t795) * t753;
t599 = (t626 * t821 + t635 * t793) * t753;
t598 = (t625 * t822 + t634 * t794) * t753;
t597 = (t624 * t823 + t633 * t795) * t753;
t596 = t642 * t806 + (t626 * t717 - t630 * t777) * t768;
t595 = t641 * t807 + (t625 * t715 - t629 * t779) * t770;
t594 = t640 * t808 + (t624 * t713 - t628 * t781) * t772;
t593 = (t623 * t824 + t632 * t796) * t753;
t592 = (t623 * t825 + t631 * t796) * t753;
t591 = t639 * t815 + (t623 * t708 - t627 * t783) * t775;
t590 = t608 * t646 + t622 * t630;
t589 = t607 * t645 + t621 * t629;
t588 = t606 * t644 + t620 * t628;
t587 = t603 * t643 + t619 * t627;
t586 = t608 * t630 + t622 * t626;
t585 = t607 * t629 + t621 * t625;
t584 = t606 * t628 + t620 * t624;
t583 = t603 * t627 + t619 * t623;
t1 = [m(4) + (t593 * t824 + t600 * t820 + t601 * t819 + t602 * t818 + (t605 * t844 + t612 * t840 + t613 * t839 + t614 * t838) * t752) * t753, (t593 * t825 + t600 * t823 + t601 * t822 + t602 * t821 + (t605 * t845 + t612 * t843 + t613 * t842 + t614 * t841) * t752) * t753, (t593 * t776 + t600 * t773 + t601 * t771 + t602 * t769 + ((t632 * t788 + t636 * t787 + t637 * t786 + t638 * t785) * t752 + (-t614 * t778 + t676 * t834) * t712 + (-t613 * t780 + t675 * t835) * t711 + (-t612 * t782 + t674 * t836) * t710 + (-t605 * t784 + t668 * t837) * t707) * t752) * t753, t593 * t619 + t600 * t620 + t601 * t621 + t602 * t622 + t605 * t603 + t612 * t606 + t613 * t607 + t614 * t608 + t762; (t592 * t824 + t597 * t820 + t598 * t819 + t599 * t818 + (t604 * t844 + t609 * t840 + t610 * t839 + t611 * t838) * t752) * t753, m(4) + (t592 * t825 + t597 * t823 + t598 * t822 + t599 * t821 + (t604 * t845 + t609 * t843 + t610 * t842 + t611 * t841) * t752) * t753, (t592 * t776 + t597 * t773 + t598 * t771 + t599 * t769 + ((t631 * t788 + t633 * t787 + t634 * t786 + t635 * t785) * t752 + (-t611 * t778 + t673 * t834) * t712 + (-t610 * t780 + t672 * t835) * t711 + (-t609 * t782 + t671 * t836) * t710 + (-t604 * t784 + t667 * t837) * t707) * t752) * t753, t592 * t619 + t597 * t620 + t598 * t621 + t599 * t622 + t604 * t603 + t609 * t606 + t610 * t607 + t611 * t608 + t763; (t591 * t824 + t594 * t820 + t595 * t819 + t596 * t818 + (t615 * t844 + t616 * t840 + t617 * t839 + t618 * t838) * t752) * t753, (t591 * t825 + t594 * t823 + t595 * t822 + t596 * t821 + (t615 * t845 + t616 * t843 + t617 * t842 + t618 * t841) * t752) * t753, m(4) + (t591 * t776 + t594 * t773 + t595 * t771 + t596 * t769) * t753 + ((t859 + t860 + t861 + t862) * t752 * Ifges(3,3) + ((t642 * t859 + (-t618 - t765) * t826) * t810 + (t641 * t860 + (-t617 - t766) * t827) * t812 + (t640 * t861 + (-t616 - t767) * t828) * t814 + (t639 * t862 + (-t615 - t774) * t829) * t817) * t753) * t752, t591 * t619 + t594 * t620 + t595 * t621 + t596 * t622 + t615 * t603 + t616 * t606 + t617 * t607 + t618 * t608; (t583 * t824 + t584 * t820 + t585 * t819 + t586 * t818 + (t587 * t844 + t588 * t840 + t589 * t839 + t590 * t838) * t752) * t753 + t762, (t583 * t825 + t584 * t823 + t585 * t822 + t586 * t821 + (t587 * t845 + t588 * t843 + t589 * t842 + t590 * t841) * t752) * t753 + t763, (t583 * t776 + t584 * t773 + t585 * t771 + t586 * t769) * t753 + ((t608 * t663 + t622 * t642) * t717 + (t607 * t662 + t621 * t641) * t715 + (t606 * t661 + t620 * t640) * t713 + (t603 * t659 + t619 * t639) * t708 + (-t587 * t707 * t784 - t588 * t710 * t782 - t589 * t711 * t780 - t590 * t712 * t778) * t753) * t752, t583 * t619 + t584 * t620 + t585 * t621 + t586 * t622 + t587 * t603 + t588 * t606 + t589 * t607 + t590 * t608 + Ifges(4,3);];
MX  = t1;
