% Calculate inertia matrix for parallel robot
% P4PRRRR8V2G2A0
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
% MX [4x4]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:22
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P4PRRRR8V2G2A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G2A0_inertia_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G2A0_inertia_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G2A0_inertia_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V2G2A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRRR8V2G2A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P4PRRRR8V2G2A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G2A0_inertia_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G2A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:16:11
% EndTime: 2020-08-07 11:16:14
% DurationCPUTime: 3.63s
% Computational Cost: add. (16254->467), mult. (34191->873), div. (1280->5), fcn. (30920->30), ass. (0->313)
t907 = (m(3) * rSges(3,2));
t908 = -2 * rSges(3,1) * t907 + 2 * Icges(3,4);
t782 = pkin(2) * m(3);
t758 = cos(pkin(4));
t759 = sin(qJ(3,4));
t761 = cos(qJ(3,4));
t820 = rSges(3,1) * t761 - rSges(3,2) * t759;
t756 = sin(pkin(4));
t760 = sin(qJ(2,4));
t859 = t756 * t760;
t685 = t820 * t758 - (t759 * rSges(3,1) + t761 * rSges(3,2)) * t859;
t906 = m(3) * t685;
t767 = sin(qJ(3,3));
t773 = cos(qJ(3,3));
t819 = rSges(3,1) * t773 - rSges(3,2) * t767;
t768 = sin(qJ(2,3));
t857 = t756 * t768;
t694 = t819 * t758 - (t767 * rSges(3,1) + t773 * rSges(3,2)) * t857;
t905 = m(3) * t694;
t769 = sin(qJ(3,2));
t775 = cos(qJ(3,2));
t818 = rSges(3,1) * t775 - rSges(3,2) * t769;
t770 = sin(qJ(2,2));
t855 = t756 * t770;
t695 = t818 * t758 - (t769 * rSges(3,1) + t775 * rSges(3,2)) * t855;
t904 = m(3) * t695;
t771 = sin(qJ(3,1));
t777 = cos(qJ(3,1));
t817 = rSges(3,1) * t777 - rSges(3,2) * t771;
t772 = sin(qJ(2,1));
t853 = t756 * t772;
t696 = t817 * t758 - (t771 * rSges(3,1) + t777 * rSges(3,2)) * t853;
t903 = m(3) * t696;
t796 = 0.1e1 / pkin(3);
t902 = m(3) * t796;
t901 = pkin(2) * t759;
t900 = pkin(2) * t767;
t899 = pkin(2) * t769;
t898 = pkin(2) * t771;
t750 = t761 ^ 2;
t897 = pkin(3) * t750;
t752 = t773 ^ 2;
t896 = pkin(3) * t752;
t753 = t775 ^ 2;
t895 = pkin(3) * t753;
t754 = t777 ^ 2;
t894 = pkin(3) * t754;
t893 = pkin(3) * t761;
t892 = pkin(3) * t773;
t891 = pkin(3) * t775;
t890 = pkin(3) * t777;
t779 = pkin(6) + rSges(3,3);
t889 = t779 * m(3);
t762 = cos(qJ(2,4));
t781 = pkin(7) + pkin(6);
t723 = pkin(2) * t760 - t781 * t762;
t724 = pkin(2) * t762 + t760 * t781;
t755 = sin(pkin(8));
t757 = cos(pkin(8));
t846 = t758 * t762;
t852 = t757 * t758;
t665 = (t755 * t760 - t757 * t846) * t893 - t724 * t852 + t723 * t755;
t839 = t759 * t758;
t697 = pkin(3) * t839 + t723 * t756;
t673 = 0.1e1 / (pkin(2) * t839 + t697 * t761 + t859 * t897);
t888 = t665 * t673;
t862 = t755 * t758;
t666 = (t755 * t846 + t757 * t760) * t893 + t724 * t862 + t723 * t757;
t887 = t666 * t673;
t774 = cos(qJ(2,3));
t725 = pkin(2) * t768 - t781 * t774;
t728 = pkin(2) * t774 + t768 * t781;
t842 = t758 * t774;
t667 = (t755 * t768 - t757 * t842) * t892 - t728 * t852 + t725 * t755;
t838 = t767 * t758;
t698 = pkin(3) * t838 + t725 * t756;
t674 = 0.1e1 / (pkin(2) * t838 + t698 * t773 + t857 * t896);
t886 = t667 * t674;
t776 = cos(qJ(2,2));
t726 = pkin(2) * t770 - t781 * t776;
t729 = pkin(2) * t776 + t770 * t781;
t841 = t758 * t776;
t668 = (t755 * t770 - t757 * t841) * t891 - t729 * t852 + t726 * t755;
t837 = t769 * t758;
t699 = pkin(3) * t837 + t726 * t756;
t675 = 0.1e1 / (pkin(2) * t837 + t699 * t775 + t855 * t895);
t885 = t668 * t675;
t778 = cos(qJ(2,1));
t727 = pkin(2) * t772 - t781 * t778;
t730 = pkin(2) * t778 + t772 * t781;
t840 = t758 * t778;
t669 = (t755 * t772 - t757 * t840) * t890 - t730 * t852 + t727 * t755;
t836 = t771 * t758;
t700 = pkin(3) * t836 + t727 * t756;
t676 = 0.1e1 / (pkin(2) * t836 + t700 * t777 + t853 * t894);
t884 = t669 * t676;
t670 = (t755 * t842 + t757 * t768) * t892 + t728 * t862 + t725 * t757;
t883 = t670 * t674;
t671 = (t755 * t841 + t757 * t770) * t891 + t729 * t862 + t726 * t757;
t882 = t671 * t675;
t672 = (t755 * t840 + t757 * t772) * t890 + t730 * t862 + t727 * t757;
t881 = t672 * t676;
t734 = m(2) * rSges(2,2) - t889;
t835 = m(2) * rSges(2,1) + t782;
t880 = ((t820 * m(3) + t835) * t762 - t760 * t734) * t756;
t879 = ((t819 * m(3) + t835) * t774 - t768 * t734) * t756;
t878 = ((t818 * m(3) + t835) * t776 - t770 * t734) * t756;
t877 = ((t817 * m(3) + t835) * t778 - t772 * t734) * t756;
t847 = t758 * t760;
t705 = t755 * t847 - t757 * t762;
t863 = t755 * t756;
t686 = t759 * t705 + t761 * t863;
t763 = legFrame(4,2);
t739 = sin(t763);
t876 = t686 * t739;
t743 = cos(t763);
t875 = t686 * t743;
t845 = t758 * t768;
t707 = t755 * t845 - t757 * t774;
t688 = t767 * t707 + t773 * t863;
t764 = legFrame(3,2);
t740 = sin(t764);
t874 = t688 * t740;
t744 = cos(t764);
t873 = t688 * t744;
t844 = t758 * t770;
t708 = t755 * t844 - t757 * t776;
t689 = t769 * t708 + t775 * t863;
t765 = legFrame(2,2);
t741 = sin(t765);
t872 = t689 * t741;
t745 = cos(t765);
t871 = t689 * t745;
t843 = t758 * t772;
t709 = t755 * t843 - t757 * t778;
t690 = t771 * t709 + t777 * t863;
t766 = legFrame(1,2);
t742 = sin(t766);
t870 = t690 * t742;
t746 = cos(t766);
t869 = t690 * t746;
t735 = -rSges(3,2) * t889 + Icges(3,6);
t736 = rSges(3,1) * t889 - Icges(3,5);
t701 = t735 * t761 - t759 * t736;
t868 = t701 * t796;
t702 = t735 * t773 - t767 * t736;
t867 = t702 * t796;
t703 = t735 * t775 - t769 * t736;
t866 = t703 * t796;
t704 = t735 * t777 - t771 * t736;
t865 = t704 * t796;
t794 = rSges(3,2) ^ 2;
t795 = rSges(3,1) ^ 2;
t733 = (t794 + t795) * m(3) + Icges(3,3);
t864 = t733 * t796;
t861 = t756 * t757;
t860 = t756 * t759;
t858 = t756 * t767;
t856 = t756 * t769;
t854 = t756 * t771;
t851 = t757 * t761;
t850 = t757 * t773;
t849 = t757 * t775;
t848 = t757 * t777;
t834 = t685 * t902;
t833 = t694 * t902;
t832 = t695 * t902;
t831 = t696 * t902;
t830 = -0.2e1 * pkin(2) * t907;
t829 = t739 * t887;
t828 = t743 * t887;
t827 = t740 * t883;
t826 = t744 * t883;
t825 = t741 * t882;
t824 = t745 * t882;
t823 = t742 * t881;
t822 = t746 * t881;
t821 = Icges(3,1) + Icges(2,3) + (pkin(2) ^ 2 + t779 ^ 2 + t794) * m(3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2);
t816 = pkin(3) * t860 - t723 * t758;
t815 = pkin(3) * t858 - t725 * t758;
t814 = pkin(3) * t856 - t726 * t758;
t813 = pkin(3) * t854 - t727 * t758;
t731 = (-t794 + t795) * m(3) + Icges(3,2) - Icges(3,1);
t747 = 0.2e1 * rSges(3,1) * t782;
t661 = t731 * t750 + (t759 * t908 + t747) * t761 + t759 * t830 + t821;
t812 = t661 * t686 + t666 * t868;
t662 = t731 * t752 + (t767 * t908 + t747) * t773 + t767 * t830 + t821;
t811 = t662 * t688 + t670 * t867;
t663 = t731 * t753 + (t769 * t908 + t747) * t775 + t769 * t830 + t821;
t810 = t663 * t689 + t671 * t866;
t664 = t731 * t754 + (t771 * t908 + t747) * t777 + t771 * t830 + t821;
t809 = t664 * t690 + t672 * t865;
t808 = t666 * t864 + t686 * t701;
t807 = t670 * t864 + t688 * t702;
t806 = t671 * t864 + t689 * t703;
t805 = t672 * t864 + t690 * t704;
t783 = xP(4);
t748 = sin(t783);
t749 = cos(t783);
t786 = koppelP(4,2);
t790 = koppelP(4,1);
t715 = -t748 * t790 - t749 * t786;
t719 = -t748 * t786 + t749 * t790;
t804 = t673 * (-t715 * t743 + t719 * t739);
t787 = koppelP(3,2);
t791 = koppelP(3,1);
t716 = -t748 * t791 - t749 * t787;
t720 = -t748 * t787 + t749 * t791;
t803 = t674 * (-t716 * t744 + t720 * t740);
t788 = koppelP(2,2);
t792 = koppelP(2,1);
t717 = -t748 * t792 - t749 * t788;
t721 = -t748 * t788 + t749 * t792;
t802 = t675 * (-t717 * t745 + t721 * t741);
t789 = koppelP(1,2);
t793 = koppelP(1,1);
t718 = -t748 * t793 - t749 * t789;
t722 = -t748 * t789 + t749 * t793;
t801 = t676 * (-t718 * t746 + t722 * t742);
t800 = t666 * t834 + t686 * t880;
t799 = t670 * t833 + t688 * t879;
t798 = t671 * t832 + t689 * t878;
t797 = t672 * t831 + t690 * t877;
t785 = rSges(4,1);
t784 = rSges(4,2);
t751 = m(1) + m(2) + m(3);
t714 = m(4) * (-t748 * t784 + t749 * t785);
t713 = m(4) * (-t748 * t785 - t749 * t784);
t712 = t755 * t778 + t757 * t843;
t711 = t755 * t776 + t757 * t844;
t710 = t755 * t774 + t757 * t845;
t706 = t755 * t762 + t757 * t847;
t693 = -t771 * t712 - t756 * t848;
t692 = -t769 * t711 - t756 * t849;
t691 = -t767 * t710 - t756 * t850;
t687 = -t759 * t706 - t756 * t851;
t680 = -t730 * t755 + t813 * t757;
t679 = -t729 * t755 + t814 * t757;
t678 = -t728 * t755 + t815 * t757;
t677 = -t724 * t755 + t816 * t757;
t660 = -t709 * t894 + t730 * t848 + (pkin(2) * t854 + t813 * t777) * t755;
t659 = -t708 * t895 + t729 * t849 + (pkin(2) * t856 + t814 * t775) * t755;
t658 = -t707 * t896 + t728 * t850 + (pkin(2) * t858 + t815 * t773) * t755;
t657 = -t705 * t897 + t724 * t851 + (pkin(2) * t860 + t816 * t761) * t755;
t656 = -(t712 * t742 - t746 * t853) * t894 + (t680 * t742 + t746 * t700) * t777 + (t742 * t861 + t758 * t746) * t898;
t655 = -(t711 * t741 - t745 * t855) * t895 + (t679 * t741 + t745 * t699) * t775 + (t741 * t861 + t758 * t745) * t899;
t654 = -(t710 * t740 - t744 * t857) * t896 + (t678 * t740 + t744 * t698) * t773 + (t740 * t861 + t758 * t744) * t900;
t653 = (t712 * t746 + t742 * t853) * t894 + (-t680 * t746 + t742 * t700) * t777 + (t758 * t742 - t746 * t861) * t898;
t652 = (t711 * t745 + t741 * t855) * t895 + (-t679 * t745 + t741 * t699) * t775 + (t758 * t741 - t745 * t861) * t899;
t651 = (t710 * t744 + t740 * t857) * t896 + (-t678 * t744 + t740 * t698) * t773 + (t758 * t740 - t744 * t861) * t900;
t650 = -(t706 * t739 - t743 * t859) * t897 + (t677 * t739 + t743 * t697) * t761 + (t739 * t861 + t758 * t743) * t901;
t649 = (t706 * t743 + t739 * t859) * t897 + (-t677 * t743 + t739 * t697) * t761 + (t758 * t739 - t743 * t861) * t901;
t648 = t690 * t801;
t647 = t689 * t802;
t646 = t688 * t803;
t645 = t686 * t804;
t644 = t796 * t672 * t801;
t643 = t796 * t671 * t802;
t642 = t796 * t670 * t803;
t641 = t796 * t666 * t804;
t640 = (t660 * t903 + t669 * t864 + t693 * t704) * t676;
t639 = (t659 * t904 + t668 * t864 + t692 * t703) * t675;
t638 = (t658 * t905 + t667 * t864 + t691 * t702) * t674;
t637 = (t660 * t751 + t669 * t831 + t693 * t877) * t676;
t636 = (t659 * t751 + t668 * t832 + t692 * t878) * t675;
t635 = (t658 * t751 + t667 * t833 + t691 * t879) * t674;
t634 = (t657 * t906 + t665 * t864 + t687 * t701) * t673;
t633 = (t657 * t751 + t665 * t834 + t687 * t880) * t673;
t632 = (t653 * t718 + t656 * t722) * t676;
t631 = (t652 * t717 + t655 * t721) * t675;
t630 = (t651 * t716 + t654 * t720) * t674;
t629 = (t649 * t715 + t650 * t719) * t673;
t628 = (t660 * t877 + t664 * t693 + t669 * t865) * t676;
t627 = (t659 * t878 + t663 * t692 + t668 * t866) * t675;
t626 = (t658 * t879 + t662 * t691 + t667 * t867) * t674;
t625 = (t657 * t880 + t661 * t687 + t665 * t868) * t673;
t624 = (t656 * t903 + t805 * t742) * t676;
t623 = (t655 * t904 + t806 * t741) * t675;
t622 = (t654 * t905 + t807 * t740) * t674;
t621 = (t653 * t903 - t805 * t746) * t676;
t620 = (t652 * t904 - t806 * t745) * t675;
t619 = (t651 * t905 - t807 * t744) * t674;
t618 = (t656 * t751 + t797 * t742) * t676;
t617 = (t655 * t751 + t798 * t741) * t675;
t616 = (t654 * t751 + t799 * t740) * t674;
t615 = (t653 * t751 - t797 * t746) * t676;
t614 = (t652 * t751 - t798 * t745) * t675;
t613 = (t651 * t751 - t799 * t744) * t674;
t612 = (t650 * t906 + t808 * t739) * t673;
t611 = (t649 * t906 - t808 * t743) * t673;
t610 = (t650 * t751 + t800 * t739) * t673;
t609 = (t649 * t751 - t800 * t743) * t673;
t608 = (t656 * t877 + t809 * t742) * t676;
t607 = (t655 * t878 + t810 * t741) * t675;
t606 = (t654 * t879 + t811 * t740) * t674;
t605 = (t653 * t877 - t809 * t746) * t676;
t604 = (t652 * t878 - t810 * t745) * t675;
t603 = (t651 * t879 - t811 * t744) * t674;
t602 = (t650 * t880 + t812 * t739) * t673;
t601 = (t649 * t880 - t812 * t743) * t673;
t600 = t632 * t903 + t644 * t733 + t648 * t704;
t599 = t631 * t904 + t643 * t733 + t647 * t703;
t598 = t630 * t905 + t642 * t733 + t646 * t702;
t597 = t632 * t751 + t644 * t903 + t648 * t877;
t596 = t631 * t751 + t643 * t904 + t647 * t878;
t595 = t630 * t751 + t642 * t905 + t646 * t879;
t594 = t629 * t906 + t641 * t733 + t645 * t701;
t593 = t629 * t751 + t641 * t906 + t645 * t880;
t592 = t632 * t877 + t644 * t704 + t648 * t664;
t591 = t631 * t878 + t643 * t703 + t647 * t663;
t590 = t630 * t879 + t642 * t702 + t646 * t662;
t589 = t629 * t880 + t641 * t701 + t645 * t661;
t1 = [m(4) + (-t605 * t869 + t615 * t653) * t676 + (-t604 * t871 + t614 * t652) * t675 + (-t603 * t873 + t613 * t651) * t674 + (-t601 * t875 + t609 * t649) * t673 + (-t611 * t828 - t619 * t826 - t620 * t824 - t621 * t822) * t796, (t605 * t870 + t615 * t656) * t676 + (t604 * t872 + t614 * t655) * t675 + (t603 * t874 + t613 * t654) * t674 + (t601 * t876 + t609 * t650) * t673 + (t611 * t829 + t619 * t827 + t620 * t825 + t621 * t823) * t796, (t605 * t693 + t615 * t660) * t676 + (t604 * t692 + t614 * t659) * t675 + (t603 * t691 + t613 * t658) * t674 + (t601 * t687 + t609 * t657) * t673 + (t611 * t888 + t619 * t886 + t620 * t885 + t621 * t884) * t796, t601 * t645 + t603 * t646 + t604 * t647 + t605 * t648 + t609 * t629 + t611 * t641 + t613 * t630 + t614 * t631 + t615 * t632 + t619 * t642 + t620 * t643 + t621 * t644 + t713; (-t608 * t869 + t618 * t653) * t676 + (-t607 * t871 + t617 * t652) * t675 + (-t606 * t873 + t616 * t651) * t674 + (-t602 * t875 + t610 * t649) * t673 + (-t612 * t828 - t622 * t826 - t623 * t824 - t624 * t822) * t796, m(4) + (t608 * t870 + t618 * t656) * t676 + (t607 * t872 + t617 * t655) * t675 + (t606 * t874 + t616 * t654) * t674 + (t602 * t876 + t610 * t650) * t673 + (t612 * t829 + t622 * t827 + t623 * t825 + t624 * t823) * t796, (t608 * t693 + t618 * t660) * t676 + (t607 * t692 + t617 * t659) * t675 + (t606 * t691 + t616 * t658) * t674 + (t602 * t687 + t610 * t657) * t673 + (t612 * t888 + t622 * t886 + t623 * t885 + t624 * t884) * t796, t602 * t645 + t606 * t646 + t607 * t647 + t608 * t648 + t610 * t629 + t612 * t641 + t616 * t630 + t617 * t631 + t618 * t632 + t622 * t642 + t623 * t643 + t624 * t644 + t714; (-t628 * t869 + t637 * t653) * t676 + (-t627 * t871 + t636 * t652) * t675 + (-t626 * t873 + t635 * t651) * t674 + (-t625 * t875 + t633 * t649) * t673 + (-t634 * t828 - t638 * t826 - t639 * t824 - t640 * t822) * t796, (t628 * t870 + t637 * t656) * t676 + (t627 * t872 + t636 * t655) * t675 + (t626 * t874 + t635 * t654) * t674 + (t625 * t876 + t633 * t650) * t673 + (t634 * t829 + t638 * t827 + t639 * t825 + t640 * t823) * t796, m(4) + (t628 * t693 + t637 * t660) * t676 + (t627 * t692 + t636 * t659) * t675 + (t626 * t691 + t635 * t658) * t674 + (t625 * t687 + t633 * t657) * t673 + (t634 * t888 + t638 * t886 + t639 * t885 + t640 * t884) * t796, t625 * t645 + t626 * t646 + t627 * t647 + t628 * t648 + t633 * t629 + t635 * t630 + t636 * t631 + t637 * t632 + t634 * t641 + t638 * t642 + t639 * t643 + t640 * t644; t713 + (-t592 * t869 + t597 * t653) * t676 + (-t591 * t871 + t596 * t652) * t675 + (-t590 * t873 + t595 * t651) * t674 + (-t589 * t875 + t593 * t649) * t673 + (-t594 * t828 - t598 * t826 - t599 * t824 - t600 * t822) * t796, t714 + (t592 * t870 + t597 * t656) * t676 + (t591 * t872 + t596 * t655) * t675 + (t590 * t874 + t595 * t654) * t674 + (t589 * t876 + t593 * t650) * t673 + (t594 * t829 + t598 * t827 + t599 * t825 + t600 * t823) * t796, (t592 * t693 + t597 * t660) * t676 + (t591 * t692 + t596 * t659) * t675 + (t590 * t691 + t595 * t658) * t674 + (t589 * t687 + t593 * t657) * t673 + (t594 * t888 + t598 * t886 + t599 * t885 + t600 * t884) * t796, t597 * t632 + t592 * t648 + t600 * t644 + t596 * t631 + t591 * t647 + t599 * t643 + t595 * t630 + t590 * t646 + t598 * t642 + t593 * t629 + t589 * t645 + t594 * t641 + Icges(4,3) + m(4) * (t784 ^ 2 + t785 ^ 2);];
MX  = t1;
