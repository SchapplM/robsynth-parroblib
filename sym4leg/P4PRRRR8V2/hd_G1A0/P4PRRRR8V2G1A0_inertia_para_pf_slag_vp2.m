% Calculate inertia matrix for parallel robot
% P4PRRRR8V2G1A0
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
% MX [4x4]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:13
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P4PRRRR8V2G1A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G1A0_inertia_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G1A0_inertia_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G1A0_inertia_para_pf_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V2G1A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRRR8V2G1A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P4PRRRR8V2G1A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G1A0_inertia_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G1A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:10:05
% EndTime: 2020-08-07 11:10:08
% DurationCPUTime: 3.06s
% Computational Cost: add. (13175->378), mult. (25235->680), div. (896->9), fcn. (25640->30), ass. (0->290)
t781 = cos(qJ(2,1));
t775 = sin(qJ(2,1));
t783 = pkin(7) + pkin(6);
t834 = t775 * t783;
t724 = pkin(2) * t781 + t834;
t757 = sin(pkin(8));
t759 = cos(pkin(8));
t733 = t783 * t781;
t721 = pkin(2) * t775 - t733;
t760 = cos(pkin(4));
t758 = sin(pkin(4));
t774 = sin(qJ(3,1));
t848 = t758 * t774;
t799 = pkin(3) * t848 - t721 * t760;
t891 = t724 * t759 + t757 * t799;
t779 = cos(qJ(2,2));
t773 = sin(qJ(2,2));
t836 = t773 * t783;
t723 = pkin(2) * t779 + t836;
t732 = t783 * t779;
t720 = pkin(2) * t773 - t732;
t772 = sin(qJ(3,2));
t850 = t758 * t772;
t800 = pkin(3) * t850 - t720 * t760;
t890 = t723 * t759 + t757 * t800;
t777 = cos(qJ(2,3));
t771 = sin(qJ(2,3));
t838 = t771 * t783;
t722 = pkin(2) * t777 + t838;
t731 = t783 * t777;
t719 = pkin(2) * t771 - t731;
t770 = sin(qJ(3,3));
t852 = t758 * t770;
t801 = pkin(3) * t852 - t719 * t760;
t889 = t722 * t759 + t757 * t801;
t768 = cos(qJ(2,4));
t766 = sin(qJ(2,4));
t841 = t766 * t783;
t718 = pkin(2) * t768 + t841;
t727 = t783 * t768;
t717 = pkin(2) * t766 - t727;
t765 = sin(qJ(3,4));
t854 = t758 * t765;
t802 = pkin(3) * t854 - t717 * t760;
t888 = t718 * t759 + t757 * t802;
t767 = cos(qJ(3,4));
t751 = t767 ^ 2;
t887 = pkin(3) * t751;
t776 = cos(qJ(3,3));
t754 = t776 ^ 2;
t886 = pkin(3) * t754;
t778 = cos(qJ(3,2));
t755 = t778 ^ 2;
t885 = pkin(3) * t755;
t780 = cos(qJ(3,1));
t756 = t780 ^ 2;
t884 = pkin(3) * t756;
t883 = m(3) * pkin(2) + mrSges(2,1);
t761 = legFrame(4,3);
t736 = sin(t761);
t740 = cos(t761);
t681 = -t757 * t736 + t740 * t759;
t685 = t759 * t736 + t740 * t757;
t726 = t767 * pkin(3) + pkin(2);
t702 = t766 * t726 - t727;
t862 = (t726 * t768 + t841) * t760;
t649 = -t681 * t862 + t702 * t685;
t840 = t767 * t758;
t842 = t765 * t760;
t665 = 0.1e1 / (t702 * t840 + t726 * t842);
t882 = t649 * t665;
t650 = -t702 * t681 - t685 * t862;
t881 = t650 * t665;
t762 = legFrame(3,3);
t737 = sin(t762);
t741 = cos(t762);
t682 = -t757 * t737 + t741 * t759;
t686 = t759 * t737 + t741 * t757;
t728 = t776 * pkin(3) + pkin(2);
t706 = t771 * t728 - t731;
t861 = (t728 * t777 + t838) * t760;
t651 = -t682 * t861 + t706 * t686;
t833 = t776 * t758;
t839 = t770 * t760;
t667 = 0.1e1 / (t706 * t833 + t728 * t839);
t880 = t651 * t667;
t763 = legFrame(2,3);
t738 = sin(t763);
t742 = cos(t763);
t683 = -t757 * t738 + t742 * t759;
t687 = t759 * t738 + t742 * t757;
t729 = t778 * pkin(3) + pkin(2);
t707 = t773 * t729 - t732;
t860 = (t729 * t779 + t836) * t760;
t652 = -t683 * t860 + t707 * t687;
t832 = t778 * t758;
t837 = t772 * t760;
t668 = 0.1e1 / (t707 * t832 + t729 * t837);
t879 = t652 * t668;
t764 = legFrame(1,3);
t739 = sin(t764);
t743 = cos(t764);
t684 = -t757 * t739 + t743 * t759;
t688 = t759 * t739 + t743 * t757;
t730 = t780 * pkin(3) + pkin(2);
t708 = t775 * t730 - t733;
t859 = (t730 * t781 + t834) * t760;
t653 = -t684 * t859 + t708 * t688;
t831 = t780 * t758;
t835 = t774 * t760;
t669 = 0.1e1 / (t708 * t831 + t730 * t835);
t878 = t653 * t669;
t654 = -t706 * t682 - t686 * t861;
t877 = t654 * t667;
t655 = -t707 * t683 - t687 * t860;
t876 = t655 * t668;
t656 = -t708 * t684 - t688 * t859;
t875 = t656 * t669;
t853 = t758 * t766;
t657 = 0.1e1 / (t853 * t887 + (pkin(3) * t842 + t717 * t758) * t767 + pkin(2) * t842);
t752 = m(1) + m(2) + m(3);
t874 = t657 * t752;
t851 = t758 * t771;
t658 = 0.1e1 / (t851 * t886 + (pkin(3) * t839 + t719 * t758) * t776 + pkin(2) * t839);
t873 = t658 * t752;
t849 = t758 * t773;
t659 = 0.1e1 / (t849 * t885 + (pkin(3) * t837 + t720 * t758) * t778 + pkin(2) * t837);
t872 = t659 * t752;
t847 = t758 * t775;
t660 = 0.1e1 / (t847 * t884 + (pkin(3) * t835 + t721 * t758) * t780 + pkin(2) * t835);
t871 = t660 * t752;
t798 = 0.1e1 / pkin(3);
t870 = t665 * t798;
t869 = t667 * t798;
t868 = t668 * t798;
t867 = t669 * t798;
t734 = m(3) * pkin(6) - mrSges(2,2) + mrSges(3,3);
t809 = t767 * mrSges(3,1) - mrSges(3,2) * t765;
t866 = ((t809 + t883) * t768 + t766 * t734) * t758;
t808 = t776 * mrSges(3,1) - mrSges(3,2) * t770;
t865 = ((t808 + t883) * t777 + t771 * t734) * t758;
t807 = t778 * mrSges(3,1) - mrSges(3,2) * t772;
t864 = ((t807 + t883) * t779 + t773 * t734) * t758;
t806 = t780 * mrSges(3,1) - mrSges(3,2) * t774;
t863 = ((t806 + t883) * t781 + t775 * t734) * t758;
t846 = t760 * t766;
t845 = t760 * t771;
t844 = t760 * t773;
t843 = t760 * t775;
t830 = -0.2e1 * pkin(2) * mrSges(3,2);
t829 = pkin(2) * t854;
t828 = pkin(2) * t852;
t827 = pkin(2) * t850;
t826 = pkin(2) * t848;
t825 = Ifges(3,3) * t870;
t824 = Ifges(3,3) * t869;
t823 = Ifges(3,3) * t868;
t822 = Ifges(3,3) * t867;
t821 = t657 * t866;
t820 = t658 * t865;
t819 = t659 * t864;
t818 = t660 * t863;
t673 = t809 * t760 - (t765 * mrSges(3,1) + t767 * mrSges(3,2)) * t853;
t817 = t673 * t870;
t744 = -mrSges(3,2) * pkin(6) + Ifges(3,6);
t745 = mrSges(3,1) * pkin(6) - Ifges(3,5);
t691 = t744 * t767 - t765 * t745;
t816 = t691 * t870;
t674 = t808 * t760 - (t770 * mrSges(3,1) + t776 * mrSges(3,2)) * t851;
t815 = t674 * t869;
t698 = t744 * t776 - t770 * t745;
t814 = t698 * t869;
t675 = t807 * t760 - (t772 * mrSges(3,1) + t778 * mrSges(3,2)) * t849;
t813 = t675 * t868;
t699 = t744 * t778 - t772 * t745;
t812 = t699 * t868;
t676 = t806 * t760 - (t774 * mrSges(3,1) + t780 * mrSges(3,2)) * t847;
t811 = t676 * t867;
t700 = t744 * t780 - t774 * t745;
t810 = t700 * t867;
t661 = t757 * t718 - t759 * t802;
t689 = t757 * t846 - t759 * t768;
t690 = t757 * t768 + t759 * t846;
t617 = -(t689 * t740 + t736 * t690) * t887 + (-t736 * t661 + t888 * t740) * t767 + t685 * t829;
t618 = (-t736 * t689 + t690 * t740) * t887 + (t661 * t740 + t888 * t736) * t767 - t681 * t829;
t787 = xP(4);
t749 = sin(t787);
t750 = cos(t787);
t790 = koppelP(4,2);
t794 = koppelP(4,1);
t709 = -t749 * t794 - t750 * t790;
t713 = -t749 * t790 + t750 * t794;
t586 = (t617 * t709 + t618 * t713) * t657;
t641 = -t681 * t840 - (t681 * t846 + t768 * t685) * t765;
t642 = -t685 * t840 - (-t768 * t681 + t685 * t846) * t765;
t601 = (t641 * t709 + t642 * t713) * t657;
t605 = (t649 * t709 + t650 * t713) * t870;
t550 = t586 * t752 + t601 * t866 + t605 * t673;
t662 = t757 * t722 - t759 * t801;
t692 = t757 * t845 - t759 * t777;
t695 = t757 * t777 + t759 * t845;
t619 = -(t692 * t741 + t737 * t695) * t886 + (-t737 * t662 + t889 * t741) * t776 + t686 * t828;
t622 = (-t737 * t692 + t695 * t741) * t886 + (t662 * t741 + t889 * t737) * t776 - t682 * t828;
t791 = koppelP(3,2);
t795 = koppelP(3,1);
t710 = -t749 * t795 - t750 * t791;
t714 = -t749 * t791 + t750 * t795;
t590 = (t619 * t710 + t622 * t714) * t658;
t643 = -t682 * t833 - (t682 * t845 + t777 * t686) * t770;
t646 = -t686 * t833 - (-t777 * t682 + t686 * t845) * t770;
t602 = (t643 * t710 + t646 * t714) * t658;
t606 = (t651 * t710 + t654 * t714) * t869;
t555 = t590 * t752 + t602 * t865 + t606 * t674;
t663 = t757 * t723 - t759 * t800;
t693 = t757 * t844 - t759 * t779;
t696 = t757 * t779 + t759 * t844;
t620 = -(t693 * t742 + t738 * t696) * t885 + (-t738 * t663 + t890 * t742) * t778 + t687 * t827;
t623 = (-t738 * t693 + t696 * t742) * t885 + (t663 * t742 + t890 * t738) * t778 - t683 * t827;
t792 = koppelP(2,2);
t796 = koppelP(2,1);
t711 = -t749 * t796 - t750 * t792;
t715 = -t749 * t792 + t750 * t796;
t591 = (t620 * t711 + t623 * t715) * t659;
t644 = -t683 * t832 - (t683 * t844 + t779 * t687) * t772;
t647 = -t687 * t832 - (-t779 * t683 + t687 * t844) * t772;
t603 = (t644 * t711 + t647 * t715) * t659;
t607 = (t652 * t711 + t655 * t715) * t868;
t556 = t591 * t752 + t603 * t864 + t607 * t675;
t664 = t757 * t724 - t759 * t799;
t694 = t757 * t843 - t759 * t781;
t697 = t757 * t781 + t759 * t843;
t621 = -(t694 * t743 + t739 * t697) * t884 + (-t739 * t664 + t891 * t743) * t780 + t688 * t826;
t624 = (-t739 * t694 + t697 * t743) * t884 + (t664 * t743 + t891 * t739) * t780 - t684 * t826;
t793 = koppelP(1,2);
t797 = koppelP(1,1);
t712 = -t749 * t797 - t750 * t793;
t716 = -t749 * t793 + t750 * t797;
t592 = (t621 * t712 + t624 * t716) * t660;
t645 = -t684 * t831 - (t684 * t843 + t781 * t688) * t774;
t648 = -t688 * t831 - (-t781 * t684 + t688 * t843) * t774;
t604 = (t645 * t712 + t648 * t716) * t660;
t608 = (t653 * t712 + t656 * t716) * t867;
t557 = t592 * t752 + t604 * t863 + t608 * t676;
t569 = t617 * t874 + t641 * t821 + t649 * t817;
t570 = t618 * t874 + t642 * t821 + t650 * t817;
t571 = t619 * t873 + t643 * t820 + t651 * t815;
t572 = t620 * t872 + t644 * t819 + t652 * t813;
t573 = t621 * t871 + t645 * t818 + t653 * t811;
t574 = t622 * t873 + t646 * t820 + t654 * t815;
t575 = t623 * t872 + t647 * t819 + t655 * t813;
t576 = t624 * t871 + t648 * t818 + t656 * t811;
t788 = mrSges(4,2);
t789 = mrSges(4,1);
t805 = -t749 * t788 + t750 * t789;
t804 = 0.2e1 * pkin(6) * mrSges(3,3) + Ifges(3,1) + Ifges(2,3) + (pkin(2) ^ 2 + pkin(6) ^ 2) * m(3);
t803 = -t749 * t789 - t750 * t788;
t782 = mrSges(3,1) * pkin(2);
t769 = -Ifges(3,1) + Ifges(3,2);
t672 = t769 * t756 + 0.2e1 * (Ifges(3,4) * t774 + t782) * t780 + t774 * t830 + t804;
t671 = t769 * t755 + 0.2e1 * (Ifges(3,4) * t772 + t782) * t778 + t772 * t830 + t804;
t670 = t769 * t754 + 0.2e1 * (Ifges(3,4) * t770 + t782) * t776 + t770 * t830 + t804;
t666 = t769 * t751 + 0.2e1 * (Ifges(3,4) * t765 + t782) * t767 + t765 * t830 + t804;
t584 = t656 * t822 + (t624 * t676 + t648 * t700) * t660;
t583 = t655 * t823 + (t623 * t675 + t647 * t699) * t659;
t582 = t654 * t824 + (t622 * t674 + t646 * t698) * t658;
t581 = t653 * t822 + (t621 * t676 + t645 * t700) * t660;
t580 = t652 * t823 + (t620 * t675 + t644 * t699) * t659;
t579 = t651 * t824 + (t619 * t674 + t643 * t698) * t658;
t578 = t650 * t825 + (t618 * t673 + t642 * t691) * t657;
t577 = t649 * t825 + (t617 * t673 + t641 * t691) * t657;
t568 = t656 * t810 + (t624 * t863 + t648 * t672) * t660;
t567 = t655 * t812 + (t623 * t864 + t647 * t671) * t659;
t566 = t654 * t814 + (t622 * t865 + t646 * t670) * t658;
t565 = t653 * t810 + (t621 * t863 + t645 * t672) * t660;
t564 = t652 * t812 + (t620 * t864 + t644 * t671) * t659;
t563 = t651 * t814 + (t619 * t865 + t643 * t670) * t658;
t562 = t650 * t816 + (t618 * t866 + t642 * t666) * t657;
t561 = t649 * t816 + (t617 * t866 + t641 * t666) * t657;
t560 = t608 * Ifges(3,3) + t592 * t676 + t604 * t700;
t559 = t607 * Ifges(3,3) + t591 * t675 + t603 * t699;
t558 = t606 * Ifges(3,3) + t590 * t674 + t602 * t698;
t554 = t605 * Ifges(3,3) + t586 * t673 + t601 * t691;
t553 = t592 * t863 + t604 * t672 + t608 * t700;
t552 = t591 * t864 + t603 * t671 + t607 * t699;
t551 = t590 * t865 + t602 * t670 + t606 * t698;
t549 = t586 * t866 + t601 * t666 + t605 * t691;
t548 = t576 + t575 + t574 + t570;
t547 = t573 + t572 + t571 + t569;
t546 = t557 + t556 + t555 + t550;
t1 = [m(4) + (t565 * t645 + t573 * t621) * t660 + (t564 * t644 + t572 * t620) * t659 + (t563 * t643 + t571 * t619) * t658 + (t561 * t641 + t569 * t617) * t657 + (t577 * t882 + t579 * t880 + t580 * t879 + t581 * t878) * t798, (t565 * t648 + t573 * t624) * t660 + (t564 * t647 + t572 * t623) * t659 + (t563 * t646 + t571 * t622) * t658 + (t561 * t642 + t569 * t618) * t657 + (t577 * t881 + t579 * t877 + t580 * t876 + t581 * t875) * t798, t547, t561 * t601 + t563 * t602 + t564 * t603 + t565 * t604 + t569 * t586 + t571 * t590 + t572 * t591 + t573 * t592 + t577 * t605 + t579 * t606 + t580 * t607 + t581 * t608 + t803; (t568 * t645 + t576 * t621) * t660 + (t567 * t644 + t575 * t620) * t659 + (t566 * t643 + t574 * t619) * t658 + (t562 * t641 + t570 * t617) * t657 + (t578 * t882 + t582 * t880 + t583 * t879 + t584 * t878) * t798, m(4) + (t568 * t648 + t576 * t624) * t660 + (t567 * t647 + t575 * t623) * t659 + (t566 * t646 + t574 * t622) * t658 + (t562 * t642 + t570 * t618) * t657 + (t578 * t881 + t582 * t877 + t583 * t876 + t584 * t875) * t798, t548, t562 * t601 + t566 * t602 + t567 * t603 + t568 * t604 + t570 * t586 + t574 * t590 + t575 * t591 + t576 * t592 + t578 * t605 + t582 * t606 + t583 * t607 + t584 * t608 + t805; t547, t548, 0.4e1 * m(1) + 0.4e1 * m(2) + 0.4e1 * m(3) + m(4), t546; (t553 * t645 + t557 * t621) * t660 + (t552 * t644 + t556 * t620) * t659 + (t551 * t643 + t555 * t619) * t658 + (t549 * t641 + t550 * t617) * t657 + (t554 * t882 + t558 * t880 + t559 * t879 + t560 * t878) * t798 + t803, (t553 * t648 + t557 * t624) * t660 + (t552 * t647 + t556 * t623) * t659 + (t551 * t646 + t555 * t622) * t658 + (t549 * t642 + t550 * t618) * t657 + (t554 * t881 + t558 * t877 + t559 * t876 + t560 * t875) * t798 + t805, t546, t549 * t601 + t550 * t586 + t551 * t602 + t552 * t603 + t553 * t604 + t554 * t605 + t555 * t590 + t556 * t591 + t557 * t592 + t558 * t606 + t559 * t607 + t560 * t608 + Ifges(4,3);];
MX  = t1;
