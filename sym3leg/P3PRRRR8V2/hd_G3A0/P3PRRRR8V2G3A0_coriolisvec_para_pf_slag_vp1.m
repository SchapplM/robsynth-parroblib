% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRRR8V2G3A0
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
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:05
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:04:15
% EndTime: 2020-08-06 18:04:21
% DurationCPUTime: 5.89s
% Computational Cost: add. (33855->311), mult. (76560->590), div. (4716->7), fcn. (74370->22), ass. (0->221)
t762 = sin(qJ(2,1));
t768 = cos(qJ(2,1));
t775 = pkin(7) + pkin(6);
t705 = pkin(2) * t762 - t768 * t775;
t751 = sin(pkin(4));
t753 = cos(pkin(4));
t761 = sin(qJ(3,1));
t824 = t761 * t753;
t689 = pkin(3) * t824 + t751 * t705;
t767 = cos(qJ(3,1));
t840 = t751 * t762;
t749 = t767 ^ 2;
t864 = pkin(3) * t749;
t671 = 0.1e1 / (pkin(2) * t824 + t689 * t767 + t840 * t864);
t760 = sin(qJ(2,2));
t766 = cos(qJ(2,2));
t704 = pkin(2) * t760 - t766 * t775;
t759 = sin(qJ(3,2));
t826 = t759 * t753;
t688 = pkin(3) * t826 + t751 * t704;
t765 = cos(qJ(3,2));
t842 = t751 * t760;
t748 = t765 ^ 2;
t865 = pkin(3) * t748;
t670 = 0.1e1 / (pkin(2) * t826 + t688 * t765 + t842 * t865);
t758 = sin(qJ(2,3));
t764 = cos(qJ(2,3));
t703 = pkin(2) * t758 - t764 * t775;
t757 = sin(qJ(3,3));
t828 = t757 * t753;
t687 = pkin(3) * t828 + t751 * t703;
t763 = cos(qJ(3,3));
t844 = t751 * t758;
t747 = t763 ^ 2;
t866 = pkin(3) * t747;
t669 = 0.1e1 / (pkin(2) * t828 + t687 * t763 + t844 * t866);
t706 = pkin(2) * t764 + t758 * t775;
t750 = sin(pkin(8));
t752 = cos(pkin(8));
t832 = t753 * t764;
t836 = t752 * t753;
t863 = t763 * pkin(3);
t660 = (t750 * t758 - t752 * t832) * t863 - t706 * t836 + t703 * t750;
t835 = t753 * t758;
t699 = t750 * t764 + t752 * t835;
t839 = t751 * t763;
t681 = t757 * t699 + t752 * t839;
t696 = t750 * t835 - t752 * t764;
t678 = t757 * t696 + t750 * t839;
t772 = xDP(3);
t754 = legFrame(3,2);
t735 = sin(t754);
t738 = cos(t754);
t773 = xDP(2);
t774 = xDP(1);
t787 = t735 * t773 - t738 * t774;
t648 = (t678 * t772 + t681 * t787) * t669;
t856 = t648 * t775;
t804 = t757 * t856;
t846 = t750 * t753;
t663 = (t750 * t832 + t752 * t758) * t863 + t706 * t846 + t703 * t752;
t780 = 0.1e1 / pkin(3);
t642 = (-t660 * t787 + t663 * t772) * t780 * t669;
t869 = pkin(3) * t642;
t633 = t804 - t869;
t807 = t757 * t869;
t822 = t764 * t648;
t827 = t758 * t763;
t845 = t751 * t757;
t615 = (((t753 * t642 + t751 * t822) * t866 + ((-t807 + t856) * t758 + pkin(2) * t822) * t839 + t633 * t753) * t648 + (t764 * t751 * t642 + (t747 * t753 - t827 * t845 - t753) * t648) * t869) * t669;
t693 = pkin(3) * t827 + t703;
t829 = t753 * t780;
t781 = pkin(2) ^ 2;
t732 = t775 ^ 2 + t781;
t779 = pkin(3) ^ 2;
t883 = 0.2e1 * pkin(2);
t819 = pkin(3) * t883;
t859 = (-t775 * t807 + (t747 * t779 + t763 * t819 + t732) * t648) * t648;
t618 = t669 * t829 * t859 + (-t753 * t804 + (-t693 * t845 + (pkin(2) * t763 + t866) * t753) * t642) / (t693 * t839 + (pkin(2) + t863) * t828) * t642;
t624 = (-t763 * t859 - (pkin(2) * t642 - t633 * t763) * t869) * t669;
t645 = t648 ^ 2;
t769 = pkin(6) + rSges(3,3);
t860 = t769 * m(3);
t719 = rSges(3,2) * t860 - Icges(3,6);
t720 = rSges(3,1) * t860 - Icges(3,5);
t690 = t719 * t763 + t720 * t757;
t777 = rSges(3,2) ^ 2;
t778 = rSges(3,1) ^ 2;
t702 = (t778 / 0.2e1 - t777 / 0.2e1) * m(3) - Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1;
t716 = -(t777 + t778) * m(3) - Icges(3,3);
t776 = pkin(2) * m(3);
t808 = t776 / 0.2e1;
t734 = rSges(3,2) * t808;
t881 = m(3) * rSges(3,1);
t809 = rSges(3,2) * t881;
t812 = -t809 / 0.2e1 + Icges(3,4) / 0.2e1;
t712 = t757 * rSges(3,1) + t763 * rSges(3,2);
t793 = rSges(3,1) * t763 - rSges(3,2) * t757;
t876 = m(3) * (-t712 * t844 + t753 * t793);
t731 = -Icges(3,4) + t809;
t801 = rSges(3,1) * t808;
t886 = t731 * t747 + t757 * t801;
t790 = (t690 * t615 + t716 * t618 - t624 * t876 + 0.2e1 * t645 * ((t702 * t757 + t734) * t763 + t812 + t886)) * t780;
t715 = (-t777 + t778) * m(3) + Icges(3,2) - Icges(3,1);
t718 = t769 ^ 2 + t777 + t781;
t742 = t881 * t883;
t794 = -(rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) - Icges(3,1) - Icges(2,3);
t811 = rSges(3,2) * t883;
t717 = m(2) * rSges(2,2) - t860;
t727 = m(2) * rSges(2,1) + t776;
t896 = m(3) * t793;
t853 = ((t727 + t896) * t764 - t758 * t717) * t751;
t877 = -t731 / 0.2e1;
t878 = t720 / 0.4e1;
t879 = -t719 / 0.4e1;
t880 = t715 / 0.2e1;
t882 = -0.2e1 * t731;
t800 = 0.4e1 * ((t879 * t757 + t878 * t763) * t642 + ((t757 * t880 + t734) * t763 + t877 + t886) * t648) * t642 + t624 * t853 - (-t715 * t747 - (t757 * t882 + t742) * t763 + (t757 * t811 - t718) * m(3) + t794) * t615 - t690 * t618;
t901 = t660 * t790 + t800 * t681;
t707 = pkin(2) * t766 + t760 * t775;
t831 = t753 * t766;
t862 = t765 * pkin(3);
t661 = (t750 * t760 - t752 * t831) * t862 - t707 * t836 + t704 * t750;
t834 = t753 * t760;
t700 = t750 * t766 + t752 * t834;
t838 = t751 * t765;
t682 = t759 * t700 + t752 * t838;
t697 = t750 * t834 - t752 * t766;
t679 = t759 * t697 + t750 * t838;
t755 = legFrame(2,2);
t736 = sin(t755);
t739 = cos(t755);
t786 = t736 * t773 - t739 * t774;
t649 = (t679 * t772 + t682 * t786) * t670;
t855 = t649 * t775;
t803 = t759 * t855;
t664 = (t750 * t831 + t752 * t760) * t862 + t707 * t846 + t704 * t752;
t643 = (-t661 * t786 + t664 * t772) * t780 * t670;
t868 = pkin(3) * t643;
t634 = t803 - t868;
t806 = t759 * t868;
t821 = t766 * t649;
t825 = t760 * t765;
t843 = t751 * t759;
t616 = (((t753 * t643 + t751 * t821) * t865 + ((-t806 + t855) * t760 + pkin(2) * t821) * t838 + t634 * t753) * t649 + (t766 * t751 * t643 + (t748 * t753 - t825 * t843 - t753) * t649) * t868) * t670;
t694 = pkin(3) * t825 + t704;
t858 = (-t775 * t806 + (t748 * t779 + t765 * t819 + t732) * t649) * t649;
t619 = t670 * t829 * t858 + (-t753 * t803 + (-t694 * t843 + (pkin(2) * t765 + t865) * t753) * t643) / (t694 * t838 + (pkin(2) + t862) * t826) * t643;
t625 = (-t765 * t858 - (pkin(2) * t643 - t634 * t765) * t868) * t670;
t646 = t649 ^ 2;
t691 = t719 * t765 + t720 * t759;
t713 = t759 * rSges(3,1) + t765 * rSges(3,2);
t792 = rSges(3,1) * t765 - rSges(3,2) * t759;
t875 = m(3) * (-t713 * t842 + t753 * t792);
t885 = t731 * t748 + t759 * t801;
t789 = (t691 * t616 + t716 * t619 - t625 * t875 + 0.2e1 * t646 * ((t702 * t759 + t734) * t765 + t812 + t885)) * t780;
t897 = m(3) * t792;
t852 = ((t727 + t897) * t766 - t760 * t717) * t751;
t799 = 0.4e1 * ((t879 * t759 + t878 * t765) * t643 + ((t759 * t880 + t734) * t765 + t877 + t885) * t649) * t643 + t625 * t852 - (-t715 * t748 - (t759 * t882 + t742) * t765 + (t759 * t811 - t718) * m(3) + t794) * t616 - t691 * t619;
t900 = t661 * t789 + t799 * t682;
t708 = pkin(2) * t768 + t762 * t775;
t830 = t753 * t768;
t861 = t767 * pkin(3);
t662 = (t750 * t762 - t752 * t830) * t861 - t708 * t836 + t705 * t750;
t833 = t753 * t762;
t701 = t750 * t768 + t752 * t833;
t837 = t751 * t767;
t683 = t761 * t701 + t752 * t837;
t698 = t750 * t833 - t752 * t768;
t680 = t761 * t698 + t750 * t837;
t756 = legFrame(1,2);
t737 = sin(t756);
t740 = cos(t756);
t785 = t737 * t773 - t740 * t774;
t650 = (t680 * t772 + t683 * t785) * t671;
t854 = t650 * t775;
t802 = t761 * t854;
t665 = (t750 * t830 + t752 * t762) * t861 + t708 * t846 + t705 * t752;
t644 = (-t662 * t785 + t665 * t772) * t780 * t671;
t867 = pkin(3) * t644;
t635 = t802 - t867;
t805 = t761 * t867;
t820 = t768 * t650;
t823 = t762 * t767;
t841 = t751 * t761;
t617 = (((t753 * t644 + t751 * t820) * t864 + ((-t805 + t854) * t762 + pkin(2) * t820) * t837 + t635 * t753) * t650 + (t644 * t751 * t768 + (t749 * t753 - t823 * t841 - t753) * t650) * t867) * t671;
t695 = pkin(3) * t823 + t705;
t857 = (-t775 * t805 + (t749 * t779 + t767 * t819 + t732) * t650) * t650;
t620 = t671 * t829 * t857 + (-t753 * t802 + (-t695 * t841 + (pkin(2) * t767 + t864) * t753) * t644) / (t695 * t837 + (pkin(2) + t861) * t824) * t644;
t626 = (-t767 * t857 - (pkin(2) * t644 - t635 * t767) * t867) * t671;
t647 = t650 ^ 2;
t692 = t719 * t767 + t720 * t761;
t714 = t761 * rSges(3,1) + t767 * rSges(3,2);
t791 = rSges(3,1) * t767 - rSges(3,2) * t761;
t874 = m(3) * (-t714 * t840 + t753 * t791);
t884 = t731 * t749 + t761 * t801;
t788 = (t692 * t617 + t716 * t620 - t626 * t874 + 0.2e1 * t647 * ((t702 * t761 + t734) * t767 + t812 + t884)) * t780;
t898 = m(3) * t791;
t851 = ((t727 + t898) * t768 - t762 * t717) * t751;
t798 = 0.4e1 * ((t879 * t761 + t878 * t767) * t644 + ((t761 * t880 + t734) * t767 + t877 + t884) * t650) * t644 + t626 * t851 - (-t715 * t749 - (t761 * t882 + t742) * t767 + (t761 * t811 - t718) * m(3) + t794) * t617 - t692 * t620;
t899 = t662 * t788 + t683 * t798;
t873 = m(3) * t753;
t872 = pkin(2) * t757;
t871 = pkin(2) * t759;
t870 = pkin(2) * t761;
t847 = t750 * t751;
t639 = t642 ^ 2;
t746 = -m(1) - m(2) - m(3);
t810 = 0.2e1 * m(3);
t818 = -t615 * t853 - t618 * t876 + t746 * t624 + ((-t645 * t727 - (t645 + t639) * t896) * t758 - (t642 * t712 * t810 + t648 * t717) * t822) * t751 - t639 * t712 * t873;
t640 = t643 ^ 2;
t817 = -t616 * t852 - t619 * t875 + t746 * t625 + ((-t646 * t727 - (t646 + t640) * t897) * t760 - (t643 * t713 * t810 + t649 * t717) * t821) * t751 - t640 * t713 * t873;
t641 = t644 ^ 2;
t816 = -t617 * t851 - t620 * t874 + t746 * t626 + ((-t647 * t727 - (t647 + t641) * t898) * t762 - (t644 * t714 * t810 + t650 * t717) * t820) * t751 - t641 * t714 * t873;
t784 = pkin(3) * t845 - t703 * t753;
t783 = pkin(3) * t843 - t704 * t753;
t782 = pkin(3) * t841 - t705 * t753;
t674 = t752 * t708 + t750 * t782;
t673 = t752 * t707 + t750 * t783;
t672 = t752 * t706 + t750 * t784;
t1 = [(t816 * (-(t698 * t740 - t737 * t840) * t864 + (t674 * t740 + t737 * t689) * t767 + (t753 * t737 + t740 * t847) * t870) + t899 * t740) * t671 + (t817 * (-(t697 * t739 - t736 * t842) * t865 + (t673 * t739 + t736 * t688) * t765 + (t753 * t736 + t739 * t847) * t871) + t900 * t739) * t670 + (t818 * (-(t696 * t738 - t735 * t844) * t866 + (t672 * t738 + t735 * t687) * t763 + (t753 * t735 + t738 * t847) * t872) + t901 * t738) * t669; (t816 * ((t698 * t737 + t740 * t840) * t864 + (-t674 * t737 + t740 * t689) * t767 + (-t737 * t847 + t740 * t753) * t870) - t899 * t737) * t671 + (t817 * ((t697 * t736 + t739 * t842) * t865 + (-t673 * t736 + t739 * t688) * t765 + (-t736 * t847 + t739 * t753) * t871) - t900 * t736) * t670 + (t818 * ((t696 * t735 + t738 * t844) * t866 + (-t672 * t735 + t738 * t687) * t763 + (-t735 * t847 + t738 * t753) * t872) - t901 * t735) * t669; (-t798 * t680 + t665 * t788 + t816 * (-t701 * t864 - t708 * t750 * t767 + (pkin(2) * t841 + t767 * t782) * t752)) * t671 + (-t799 * t679 + t664 * t789 + t817 * (-t700 * t865 - t707 * t750 * t765 + (pkin(2) * t843 + t765 * t783) * t752)) * t670 + (-t800 * t678 + t663 * t790 + t818 * (-t699 * t866 - t706 * t750 * t763 + (pkin(2) * t845 + t763 * t784) * t752)) * t669;];
taucX  = t1;
