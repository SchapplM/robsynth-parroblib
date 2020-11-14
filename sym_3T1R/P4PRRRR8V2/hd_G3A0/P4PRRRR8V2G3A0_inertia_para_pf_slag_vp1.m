% Calculate inertia matrix for parallel robot
% P4PRRRR8V2G3A0
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
% Datum: 2020-08-07 11:30
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P4PRRRR8V2G3A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G3A0_inertia_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G3A0_inertia_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G3A0_inertia_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V2G3A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRRR8V2G3A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P4PRRRR8V2G3A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G3A0_inertia_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G3A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:25:08
% EndTime: 2020-08-07 11:25:11
% DurationCPUTime: 3.65s
% Computational Cost: add. (16254->467), mult. (34191->872), div. (1280->5), fcn. (30920->30), ass. (0->312)
t904 = (m(3) * rSges(3,2));
t905 = -2 * rSges(3,1) * t904 + 2 * Icges(3,4);
t780 = pkin(2) * m(3);
t756 = cos(pkin(4));
t757 = sin(qJ(3,4));
t759 = cos(qJ(3,4));
t818 = rSges(3,1) * t759 - rSges(3,2) * t757;
t754 = sin(pkin(4));
t758 = sin(qJ(2,4));
t857 = t754 * t758;
t683 = t818 * t756 - (t757 * rSges(3,1) + t759 * rSges(3,2)) * t857;
t903 = m(3) * t683;
t765 = sin(qJ(3,3));
t771 = cos(qJ(3,3));
t817 = rSges(3,1) * t771 - rSges(3,2) * t765;
t766 = sin(qJ(2,3));
t854 = t754 * t766;
t692 = t817 * t756 - (t765 * rSges(3,1) + t771 * rSges(3,2)) * t854;
t902 = m(3) * t692;
t767 = sin(qJ(3,2));
t773 = cos(qJ(3,2));
t816 = rSges(3,1) * t773 - rSges(3,2) * t767;
t768 = sin(qJ(2,2));
t852 = t754 * t768;
t693 = t816 * t756 - (t767 * rSges(3,1) + t773 * rSges(3,2)) * t852;
t901 = m(3) * t693;
t769 = sin(qJ(3,1));
t775 = cos(qJ(3,1));
t815 = rSges(3,1) * t775 - rSges(3,2) * t769;
t770 = sin(qJ(2,1));
t850 = t754 * t770;
t694 = t815 * t756 - (t769 * rSges(3,1) + t775 * rSges(3,2)) * t850;
t900 = m(3) * t694;
t794 = 0.1e1 / pkin(3);
t899 = m(3) * t794;
t898 = pkin(2) * t757;
t897 = pkin(2) * t765;
t896 = pkin(2) * t767;
t895 = pkin(2) * t769;
t748 = t759 ^ 2;
t894 = pkin(3) * t748;
t750 = t771 ^ 2;
t893 = pkin(3) * t750;
t751 = t773 ^ 2;
t892 = pkin(3) * t751;
t752 = t775 ^ 2;
t891 = pkin(3) * t752;
t890 = pkin(3) * t759;
t889 = pkin(3) * t771;
t888 = pkin(3) * t773;
t887 = pkin(3) * t775;
t777 = pkin(6) + rSges(3,3);
t886 = t777 * m(3);
t760 = cos(qJ(2,4));
t779 = pkin(7) + pkin(6);
t721 = pkin(2) * t758 - t779 * t760;
t722 = pkin(2) * t760 + t758 * t779;
t753 = sin(pkin(8));
t755 = cos(pkin(8));
t844 = t756 * t760;
t846 = t755 * t756;
t663 = (t753 * t758 - t755 * t844) * t890 - t722 * t846 + t721 * t753;
t837 = t757 * t756;
t695 = pkin(3) * t837 + t721 * t754;
t671 = 0.1e1 / (pkin(2) * t837 + t695 * t759 + t857 * t894);
t885 = t663 * t671;
t859 = t753 * t756;
t664 = (t753 * t844 + t755 * t758) * t890 + t722 * t859 + t755 * t721;
t884 = t664 * t671;
t772 = cos(qJ(2,3));
t723 = pkin(2) * t766 - t779 * t772;
t726 = pkin(2) * t772 + t766 * t779;
t840 = t756 * t772;
t665 = (t753 * t766 - t755 * t840) * t889 - t726 * t846 + t723 * t753;
t836 = t765 * t756;
t696 = pkin(3) * t836 + t723 * t754;
t672 = 0.1e1 / (pkin(2) * t836 + t696 * t771 + t854 * t893);
t883 = t665 * t672;
t774 = cos(qJ(2,2));
t724 = pkin(2) * t768 - t779 * t774;
t727 = pkin(2) * t774 + t768 * t779;
t839 = t756 * t774;
t666 = (t753 * t768 - t755 * t839) * t888 - t727 * t846 + t724 * t753;
t835 = t767 * t756;
t697 = pkin(3) * t835 + t724 * t754;
t673 = 0.1e1 / (pkin(2) * t835 + t697 * t773 + t852 * t892);
t882 = t666 * t673;
t776 = cos(qJ(2,1));
t725 = pkin(2) * t770 - t779 * t776;
t728 = pkin(2) * t776 + t770 * t779;
t838 = t756 * t776;
t667 = (t753 * t770 - t755 * t838) * t887 - t728 * t846 + t725 * t753;
t834 = t769 * t756;
t698 = pkin(3) * t834 + t725 * t754;
t674 = 0.1e1 / (pkin(2) * t834 + t698 * t775 + t850 * t891);
t881 = t667 * t674;
t668 = (t753 * t840 + t755 * t766) * t889 + t726 * t859 + t755 * t723;
t880 = t668 * t672;
t669 = (t753 * t839 + t755 * t768) * t888 + t727 * t859 + t755 * t724;
t879 = t669 * t673;
t670 = (t753 * t838 + t755 * t770) * t887 + t728 * t859 + t755 * t725;
t878 = t670 * t674;
t732 = m(2) * rSges(2,2) - t886;
t833 = m(2) * rSges(2,1) + t780;
t877 = ((t818 * m(3) + t833) * t760 - t758 * t732) * t754;
t876 = ((t817 * m(3) + t833) * t772 - t766 * t732) * t754;
t875 = ((t816 * m(3) + t833) * t774 - t768 * t732) * t754;
t874 = ((t815 * m(3) + t833) * t776 - t770 * t732) * t754;
t845 = t756 * t758;
t704 = t753 * t760 + t755 * t845;
t856 = t754 * t759;
t685 = t757 * t704 + t755 * t856;
t761 = legFrame(4,2);
t737 = sin(t761);
t873 = t685 * t737;
t741 = cos(t761);
t872 = t685 * t741;
t843 = t756 * t766;
t708 = t753 * t772 + t755 * t843;
t849 = t754 * t771;
t689 = t765 * t708 + t755 * t849;
t762 = legFrame(3,2);
t738 = sin(t762);
t871 = t689 * t738;
t742 = cos(t762);
t870 = t689 * t742;
t842 = t756 * t768;
t709 = t753 * t774 + t755 * t842;
t848 = t754 * t773;
t690 = t767 * t709 + t755 * t848;
t763 = legFrame(2,2);
t739 = sin(t763);
t869 = t690 * t739;
t743 = cos(t763);
t868 = t690 * t743;
t841 = t756 * t770;
t710 = t753 * t776 + t755 * t841;
t847 = t754 * t775;
t691 = t769 * t710 + t755 * t847;
t764 = legFrame(1,2);
t740 = sin(t764);
t867 = t691 * t740;
t744 = cos(t764);
t866 = t691 * t744;
t733 = -rSges(3,2) * t886 + Icges(3,6);
t734 = rSges(3,1) * t886 - Icges(3,5);
t699 = t733 * t759 - t757 * t734;
t865 = t699 * t794;
t700 = t733 * t771 - t765 * t734;
t864 = t700 * t794;
t701 = t733 * t773 - t767 * t734;
t863 = t701 * t794;
t702 = t733 * t775 - t769 * t734;
t862 = t702 * t794;
t792 = rSges(3,2) ^ 2;
t793 = rSges(3,1) ^ 2;
t731 = (t792 + t793) * m(3) + Icges(3,3);
t861 = t731 * t794;
t860 = t753 * t754;
t858 = t754 * t757;
t855 = t754 * t765;
t853 = t754 * t767;
t851 = t754 * t769;
t832 = t683 * t899;
t831 = t692 * t899;
t830 = t693 * t899;
t829 = t694 * t899;
t828 = -0.2e1 * pkin(2) * t904;
t827 = t737 * t885;
t826 = t741 * t885;
t825 = t738 * t883;
t824 = t742 * t883;
t823 = t739 * t882;
t822 = t743 * t882;
t821 = t740 * t881;
t820 = t744 * t881;
t819 = Icges(3,1) + Icges(2,3) + (pkin(2) ^ 2 + t777 ^ 2 + t792) * m(3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2);
t781 = xP(4);
t746 = sin(t781);
t747 = cos(t781);
t784 = koppelP(4,2);
t788 = koppelP(4,1);
t713 = -t746 * t788 - t747 * t784;
t717 = -t746 * t784 + t747 * t788;
t814 = t713 * t741 - t717 * t737;
t785 = koppelP(3,2);
t789 = koppelP(3,1);
t714 = -t746 * t789 - t747 * t785;
t718 = -t746 * t785 + t747 * t789;
t813 = t714 * t742 - t718 * t738;
t786 = koppelP(2,2);
t790 = koppelP(2,1);
t715 = -t746 * t790 - t747 * t786;
t719 = -t746 * t786 + t747 * t790;
t812 = t715 * t743 - t719 * t739;
t787 = koppelP(1,2);
t791 = koppelP(1,1);
t716 = -t746 * t791 - t747 * t787;
t720 = -t746 * t787 + t747 * t791;
t811 = t716 * t744 - t720 * t740;
t810 = pkin(3) * t858 - t721 * t756;
t809 = pkin(3) * t855 - t723 * t756;
t808 = pkin(3) * t853 - t724 * t756;
t807 = pkin(3) * t851 - t725 * t756;
t729 = (-t792 + t793) * m(3) + Icges(3,2) - Icges(3,1);
t745 = 0.2e1 * rSges(3,1) * t780;
t659 = t729 * t748 + (t757 * t905 + t745) * t759 + t757 * t828 + t819;
t806 = -t659 * t685 + t663 * t865;
t660 = t729 * t750 + (t765 * t905 + t745) * t771 + t765 * t828 + t819;
t805 = -t660 * t689 + t665 * t864;
t661 = t729 * t751 + (t767 * t905 + t745) * t773 + t767 * t828 + t819;
t804 = -t661 * t690 + t666 * t863;
t662 = t729 * t752 + (t769 * t905 + t745) * t775 + t769 * t828 + t819;
t803 = -t662 * t691 + t667 * t862;
t802 = t663 * t861 - t685 * t699;
t801 = t665 * t861 - t689 * t700;
t800 = t666 * t861 - t690 * t701;
t799 = t667 * t861 - t691 * t702;
t798 = t663 * t832 - t685 * t877;
t797 = t665 * t831 - t689 * t876;
t796 = t666 * t830 - t690 * t875;
t795 = t667 * t829 - t691 * t874;
t783 = rSges(4,1);
t782 = rSges(4,2);
t749 = m(1) + m(2) + m(3);
t712 = m(4) * (-t746 * t782 + t747 * t783);
t711 = m(4) * (-t746 * t783 - t747 * t782);
t707 = t753 * t841 - t755 * t776;
t706 = t753 * t842 - t755 * t774;
t705 = t753 * t843 - t755 * t772;
t703 = t753 * t845 - t755 * t760;
t688 = t769 * t707 + t753 * t847;
t687 = t767 * t706 + t753 * t848;
t686 = t765 * t705 + t753 * t849;
t684 = t757 * t703 + t753 * t856;
t678 = t728 * t755 + t807 * t753;
t677 = t727 * t755 + t808 * t753;
t676 = t726 * t755 + t809 * t753;
t675 = t722 * t755 + t810 * t753;
t658 = -t710 * t891 - t728 * t753 * t775 + (pkin(2) * t851 + t807 * t775) * t755;
t657 = -t709 * t892 - t727 * t753 * t773 + (pkin(2) * t853 + t808 * t773) * t755;
t656 = -t708 * t893 - t726 * t753 * t771 + (pkin(2) * t855 + t809 * t771) * t755;
t655 = -t704 * t894 - t722 * t753 * t759 + (pkin(2) * t858 + t810 * t759) * t755;
t654 = -(t707 * t744 - t740 * t850) * t891 + (t678 * t744 + t698 * t740) * t775 + (t756 * t740 + t744 * t860) * t895;
t653 = (t707 * t740 + t744 * t850) * t891 + (-t678 * t740 + t698 * t744) * t775 + (-t740 * t860 + t756 * t744) * t895;
t652 = -(t706 * t743 - t739 * t852) * t892 + (t677 * t743 + t697 * t739) * t773 + (t756 * t739 + t743 * t860) * t896;
t651 = (t706 * t739 + t743 * t852) * t892 + (-t677 * t739 + t697 * t743) * t773 + (-t739 * t860 + t756 * t743) * t896;
t650 = -(t705 * t742 - t738 * t854) * t893 + (t676 * t742 + t696 * t738) * t771 + (t756 * t738 + t742 * t860) * t897;
t649 = (t705 * t738 + t742 * t854) * t893 + (-t676 * t738 + t696 * t742) * t771 + (-t738 * t860 + t756 * t742) * t897;
t648 = -(t703 * t741 - t737 * t857) * t894 + (t675 * t741 + t695 * t737) * t759 + (t756 * t737 + t741 * t860) * t898;
t647 = (t703 * t737 + t741 * t857) * t894 + (-t675 * t737 + t695 * t741) * t759 + (-t737 * t860 + t756 * t741) * t898;
t646 = t811 * t691 * t674;
t645 = t812 * t690 * t673;
t644 = t813 * t689 * t672;
t643 = t814 * t685 * t671;
t642 = t811 * t794 * t881;
t641 = t812 * t794 * t882;
t640 = t813 * t794 * t883;
t639 = t814 * t794 * t885;
t638 = (t658 * t900 + t670 * t861 + t688 * t702) * t674;
t637 = (t657 * t901 + t669 * t861 + t687 * t701) * t673;
t636 = (t656 * t902 + t668 * t861 + t686 * t700) * t672;
t635 = (t658 * t749 + t670 * t829 + t688 * t874) * t674;
t634 = (t657 * t749 + t669 * t830 + t687 * t875) * t673;
t633 = (t656 * t749 + t668 * t831 + t686 * t876) * t672;
t632 = (t655 * t903 + t664 * t861 + t684 * t699) * t671;
t631 = (t655 * t749 + t664 * t832 + t684 * t877) * t671;
t630 = (t653 * t720 + t654 * t716) * t674;
t629 = (t651 * t719 + t652 * t715) * t673;
t628 = (t649 * t718 + t650 * t714) * t672;
t627 = (t647 * t717 + t648 * t713) * t671;
t626 = (t658 * t874 + t662 * t688 + t670 * t862) * t674;
t625 = (t657 * t875 + t661 * t687 + t669 * t863) * t673;
t624 = (t656 * t876 + t660 * t686 + t668 * t864) * t672;
t623 = (t655 * t877 + t659 * t684 + t664 * t865) * t671;
t622 = (t654 * t900 + t799 * t744) * t674;
t621 = (t653 * t900 - t799 * t740) * t674;
t620 = (t652 * t901 + t800 * t743) * t673;
t619 = (t651 * t901 - t800 * t739) * t673;
t618 = (t650 * t902 + t801 * t742) * t672;
t617 = (t649 * t902 - t801 * t738) * t672;
t616 = (t654 * t749 + t795 * t744) * t674;
t615 = (t653 * t749 - t795 * t740) * t674;
t614 = (t652 * t749 + t796 * t743) * t673;
t613 = (t651 * t749 - t796 * t739) * t673;
t612 = (t650 * t749 + t797 * t742) * t672;
t611 = (t649 * t749 - t797 * t738) * t672;
t610 = (t648 * t903 + t802 * t741) * t671;
t609 = (t647 * t903 - t802 * t737) * t671;
t608 = (t648 * t749 + t798 * t741) * t671;
t607 = (t647 * t749 - t798 * t737) * t671;
t606 = (t654 * t874 + t803 * t744) * t674;
t605 = (t653 * t874 - t803 * t740) * t674;
t604 = (t652 * t875 + t804 * t743) * t673;
t603 = (t651 * t875 - t804 * t739) * t673;
t602 = (t650 * t876 + t805 * t742) * t672;
t601 = (t649 * t876 - t805 * t738) * t672;
t600 = (t648 * t877 + t806 * t741) * t671;
t599 = (t647 * t877 - t806 * t737) * t671;
t598 = t630 * t900 + t642 * t731 - t646 * t702;
t597 = t629 * t901 + t641 * t731 - t645 * t701;
t596 = t628 * t902 + t640 * t731 - t644 * t700;
t595 = t630 * t749 + t642 * t900 - t646 * t874;
t594 = t629 * t749 + t641 * t901 - t645 * t875;
t593 = t628 * t749 + t640 * t902 - t644 * t876;
t592 = t627 * t903 + t639 * t731 - t643 * t699;
t591 = t627 * t749 + t639 * t903 - t643 * t877;
t590 = t630 * t874 + t642 * t702 - t646 * t662;
t589 = t629 * t875 + t641 * t701 - t645 * t661;
t588 = t628 * t876 + t640 * t700 - t644 * t660;
t587 = t627 * t877 + t639 * t699 - t643 * t659;
t1 = [m(4) + (-t606 * t866 + t616 * t654) * t674 + (-t604 * t868 + t614 * t652) * t673 + (-t602 * t870 + t612 * t650) * t672 + (-t600 * t872 + t608 * t648) * t671 + (t610 * t826 + t618 * t824 + t620 * t822 + t622 * t820) * t794, (t606 * t867 + t616 * t653) * t674 + (t604 * t869 + t614 * t651) * t673 + (t602 * t871 + t612 * t649) * t672 + (t600 * t873 + t608 * t647) * t671 + (-t610 * t827 - t618 * t825 - t620 * t823 - t622 * t821) * t794, (t606 * t688 + t616 * t658) * t674 + (t604 * t687 + t614 * t657) * t673 + (t602 * t686 + t612 * t656) * t672 + (t600 * t684 + t608 * t655) * t671 + (t610 * t884 + t618 * t880 + t620 * t879 + t622 * t878) * t794, -t600 * t643 - t602 * t644 - t604 * t645 - t606 * t646 + t608 * t627 + t610 * t639 + t612 * t628 + t614 * t629 + t616 * t630 + t618 * t640 + t620 * t641 + t622 * t642 + t711; (-t605 * t866 + t615 * t654) * t674 + (-t603 * t868 + t613 * t652) * t673 + (-t601 * t870 + t611 * t650) * t672 + (-t599 * t872 + t607 * t648) * t671 + (t609 * t826 + t617 * t824 + t619 * t822 + t621 * t820) * t794, m(4) + (t605 * t867 + t615 * t653) * t674 + (t603 * t869 + t613 * t651) * t673 + (t601 * t871 + t611 * t649) * t672 + (t599 * t873 + t607 * t647) * t671 + (-t609 * t827 - t617 * t825 - t619 * t823 - t621 * t821) * t794, (t605 * t688 + t615 * t658) * t674 + (t603 * t687 + t613 * t657) * t673 + (t601 * t686 + t611 * t656) * t672 + (t599 * t684 + t607 * t655) * t671 + (t609 * t884 + t617 * t880 + t619 * t879 + t621 * t878) * t794, -t599 * t643 - t601 * t644 - t603 * t645 - t605 * t646 + t607 * t627 + t609 * t639 + t611 * t628 + t613 * t629 + t615 * t630 + t617 * t640 + t619 * t641 + t621 * t642 + t712; (-t626 * t866 + t635 * t654) * t674 + (-t625 * t868 + t634 * t652) * t673 + (-t624 * t870 + t633 * t650) * t672 + (-t623 * t872 + t631 * t648) * t671 + (t632 * t826 + t636 * t824 + t637 * t822 + t638 * t820) * t794, (t626 * t867 + t635 * t653) * t674 + (t625 * t869 + t634 * t651) * t673 + (t624 * t871 + t633 * t649) * t672 + (t623 * t873 + t631 * t647) * t671 + (-t632 * t827 - t636 * t825 - t637 * t823 - t638 * t821) * t794, m(4) + (t626 * t688 + t635 * t658) * t674 + (t625 * t687 + t634 * t657) * t673 + (t624 * t686 + t633 * t656) * t672 + (t623 * t684 + t631 * t655) * t671 + (t632 * t884 + t636 * t880 + t637 * t879 + t638 * t878) * t794, -t623 * t643 - t624 * t644 - t625 * t645 - t626 * t646 + t631 * t627 + t633 * t628 + t634 * t629 + t635 * t630 + t632 * t639 + t636 * t640 + t637 * t641 + t638 * t642; t711 + (-t590 * t866 + t595 * t654) * t674 + (-t589 * t868 + t594 * t652) * t673 + (-t588 * t870 + t593 * t650) * t672 + (-t587 * t872 + t591 * t648) * t671 + (t592 * t826 + t596 * t824 + t597 * t822 + t598 * t820) * t794, t712 + (t590 * t867 + t595 * t653) * t674 + (t589 * t869 + t594 * t651) * t673 + (t588 * t871 + t593 * t649) * t672 + (t587 * t873 + t591 * t647) * t671 + (-t592 * t827 - t596 * t825 - t597 * t823 - t598 * t821) * t794, (t590 * t688 + t595 * t658) * t674 + (t589 * t687 + t594 * t657) * t673 + (t588 * t686 + t593 * t656) * t672 + (t587 * t684 + t591 * t655) * t671 + (t592 * t884 + t596 * t880 + t597 * t879 + t598 * t878) * t794, t595 * t630 - t590 * t646 + t598 * t642 + t594 * t629 - t589 * t645 + t597 * t641 + t593 * t628 - t588 * t644 + t596 * t640 + t591 * t627 - t587 * t643 + t592 * t639 + Icges(4,3) + m(4) * (t782 ^ 2 + t783 ^ 2);];
MX  = t1;
