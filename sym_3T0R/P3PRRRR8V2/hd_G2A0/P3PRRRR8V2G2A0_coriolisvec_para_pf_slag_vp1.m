% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRRR8V2G2A0
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
% Datum: 2020-08-06 17:50
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:49:09
% EndTime: 2020-08-06 17:49:15
% DurationCPUTime: 5.97s
% Computational Cost: add. (33855->311), mult. (76560->590), div. (4716->7), fcn. (74370->22), ass. (0->227)
t770 = sin(qJ(2,1));
t776 = cos(qJ(2,1));
t783 = pkin(7) + pkin(6);
t713 = pkin(2) * t770 - t776 * t783;
t759 = sin(pkin(4));
t761 = cos(pkin(4));
t769 = sin(qJ(3,1));
t829 = t769 * t761;
t697 = pkin(3) * t829 + t713 * t759;
t775 = cos(qJ(3,1));
t851 = t759 * t770;
t757 = t775 ^ 2;
t878 = pkin(3) * t757;
t679 = 0.1e1 / (pkin(2) * t829 + t697 * t775 + t851 * t878);
t768 = sin(qJ(2,2));
t774 = cos(qJ(2,2));
t712 = pkin(2) * t768 - t774 * t783;
t767 = sin(qJ(3,2));
t831 = t767 * t761;
t696 = pkin(3) * t831 + t712 * t759;
t773 = cos(qJ(3,2));
t853 = t759 * t768;
t756 = t773 ^ 2;
t879 = pkin(3) * t756;
t678 = 0.1e1 / (pkin(2) * t831 + t696 * t773 + t853 * t879);
t766 = sin(qJ(2,3));
t772 = cos(qJ(2,3));
t711 = pkin(2) * t766 - t772 * t783;
t765 = sin(qJ(3,3));
t833 = t765 * t761;
t695 = pkin(3) * t833 + t711 * t759;
t771 = cos(qJ(3,3));
t855 = t759 * t766;
t755 = t771 ^ 2;
t880 = pkin(3) * t755;
t677 = 0.1e1 / (pkin(2) * t833 + t695 * t771 + t855 * t880);
t716 = pkin(2) * t776 + t770 * t783;
t758 = sin(pkin(8));
t760 = cos(pkin(8));
t835 = t761 * t776;
t858 = t758 * t761;
t875 = t775 * pkin(3);
t673 = (t758 * t835 + t760 * t770) * t875 + t716 * t858 + t713 * t760;
t838 = t761 * t770;
t706 = t758 * t838 - t760 * t776;
t846 = t759 * t775;
t688 = t769 * t706 + t758 * t846;
t709 = t758 * t776 + t760 * t838;
t841 = t760 * t775;
t691 = -t769 * t709 - t759 * t841;
t780 = xDP(3);
t764 = legFrame(1,2);
t745 = sin(t764);
t748 = cos(t764);
t781 = xDP(2);
t782 = xDP(1);
t793 = t745 * t781 - t748 * t782;
t658 = (t688 * t793 + t691 * t780) * t679;
t865 = t658 * t783;
t810 = t769 * t865;
t844 = t760 * t761;
t670 = (t758 * t770 - t760 * t835) * t875 - t716 * t844 + t758 * t713;
t788 = 0.1e1 / pkin(3);
t652 = (t670 * t780 + t673 * t793) * t788 * t679;
t881 = pkin(3) * t652;
t643 = t810 - t881;
t813 = t769 * t881;
t828 = t770 * t775;
t845 = t759 * t776;
t852 = t759 * t769;
t866 = t658 * t776;
t625 = (((t761 * t652 + t658 * t845) * t878 + ((-t813 + t865) * t770 + pkin(2) * t866) * t846 + t643 * t761) * t658 + (t652 * t845 + (t757 * t761 - t828 * t852 - t761) * t658) * t881) * t679;
t703 = pkin(3) * t828 + t713;
t834 = t761 * t788;
t789 = pkin(2) ^ 2;
t740 = t783 ^ 2 + t789;
t787 = pkin(3) ^ 2;
t897 = 0.2e1 * pkin(2);
t827 = pkin(3) * t897;
t871 = (-t783 * t813 + (t757 * t787 + t775 * t827 + t740) * t658) * t658;
t628 = t679 * t834 * t871 + (-t761 * t810 + (-t703 * t852 + t761 * (pkin(2) * t775 + t878)) * t652) / (t703 * t846 + (pkin(2) + t875) * t829) * t652;
t634 = (-t775 * t871 - (pkin(2) * t652 - t643 * t775) * t881) * t679;
t655 = t658 ^ 2;
t777 = pkin(6) + rSges(3,3);
t874 = t777 * m(3);
t727 = rSges(3,2) * t874 - Icges(3,6);
t728 = rSges(3,1) * t874 - Icges(3,5);
t700 = t727 * t775 + t728 * t769;
t785 = rSges(3,2) ^ 2;
t786 = rSges(3,1) ^ 2;
t710 = (t786 / 0.2e1 - t785 / 0.2e1) * m(3) - Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1;
t724 = -(t785 + t786) * m(3) - Icges(3,3);
t784 = pkin(2) * m(3);
t816 = t784 / 0.2e1;
t742 = rSges(3,2) * t816;
t895 = m(3) * rSges(3,1);
t817 = rSges(3,2) * t895;
t820 = -t817 / 0.2e1 + Icges(3,4) / 0.2e1;
t722 = t769 * rSges(3,1) + t775 * rSges(3,2);
t799 = rSges(3,1) * t775 - rSges(3,2) * t769;
t888 = m(3) * (-t722 * t851 + t761 * t799);
t739 = -Icges(3,4) + t817;
t809 = rSges(3,1) * t816;
t898 = t739 * t757 + t769 * t809;
t796 = (t700 * t625 + t724 * t628 - t634 * t888 + 0.2e1 * t655 * ((t710 * t769 + t742) * t775 + t820 + t898)) * t788;
t723 = (-t785 + t786) * m(3) + Icges(3,2) - Icges(3,1);
t726 = t777 ^ 2 + t785 + t789;
t750 = t895 * t897;
t802 = -(rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) - Icges(3,1) - Icges(2,3);
t819 = rSges(3,2) * t897;
t725 = m(2) * rSges(2,2) - t874;
t735 = m(2) * rSges(2,1) + t784;
t912 = m(3) * t799;
t862 = ((t735 + t912) * t776 - t770 * t725) * t759;
t891 = -t739 / 0.2e1;
t892 = t728 / 0.4e1;
t893 = -t727 / 0.4e1;
t894 = t723 / 0.2e1;
t896 = -0.2e1 * t739;
t806 = 0.4e1 * ((t893 * t769 + t892 * t775) * t652 + ((t769 * t894 + t742) * t775 + t891 + t898) * t658) * t652 + t634 * t862 - (-t723 * t757 - (t769 * t896 + t750) * t775 + (t769 * t819 - t726) * m(3) + t802) * t625 - t700 * t628;
t915 = -t673 * t796 + t688 * t806;
t715 = pkin(2) * t774 + t768 * t783;
t836 = t761 * t774;
t876 = t773 * pkin(3);
t672 = (t758 * t836 + t760 * t768) * t876 + t715 * t858 + t712 * t760;
t839 = t761 * t768;
t705 = t758 * t839 - t760 * t774;
t848 = t759 * t773;
t687 = t767 * t705 + t758 * t848;
t708 = t758 * t774 + t760 * t839;
t842 = t760 * t773;
t690 = -t767 * t708 - t759 * t842;
t763 = legFrame(2,2);
t744 = sin(t763);
t747 = cos(t763);
t794 = t744 * t781 - t747 * t782;
t657 = (t687 * t794 + t690 * t780) * t678;
t867 = t657 * t783;
t811 = t767 * t867;
t669 = (t758 * t768 - t760 * t836) * t876 - t715 * t844 + t758 * t712;
t651 = (t669 * t780 + t672 * t794) * t788 * t678;
t882 = pkin(3) * t651;
t642 = t811 - t882;
t814 = t767 * t882;
t830 = t768 * t773;
t847 = t759 * t774;
t854 = t759 * t767;
t868 = t657 * t774;
t624 = (((t761 * t651 + t657 * t847) * t879 + ((-t814 + t867) * t768 + pkin(2) * t868) * t848 + t642 * t761) * t657 + (t651 * t847 + (t756 * t761 - t830 * t854 - t761) * t657) * t882) * t678;
t702 = pkin(3) * t830 + t712;
t872 = (-t783 * t814 + (t756 * t787 + t773 * t827 + t740) * t657) * t657;
t627 = t678 * t834 * t872 + (-t761 * t811 + (-t702 * t854 + t761 * (pkin(2) * t773 + t879)) * t651) / (t702 * t848 + (pkin(2) + t876) * t831) * t651;
t633 = (-t773 * t872 - (pkin(2) * t651 - t642 * t773) * t882) * t678;
t654 = t657 ^ 2;
t699 = t727 * t773 + t728 * t767;
t721 = t767 * rSges(3,1) + t773 * rSges(3,2);
t800 = rSges(3,1) * t773 - rSges(3,2) * t767;
t889 = m(3) * (-t721 * t853 + t761 * t800);
t899 = t739 * t756 + t767 * t809;
t797 = (t699 * t624 + t724 * t627 - t633 * t889 + 0.2e1 * t654 * ((t710 * t767 + t742) * t773 + t820 + t899)) * t788;
t911 = m(3) * t800;
t863 = ((t735 + t911) * t774 - t768 * t725) * t759;
t807 = 0.4e1 * ((t893 * t767 + t892 * t773) * t651 + ((t767 * t894 + t742) * t773 + t891 + t899) * t657) * t651 + t633 * t863 - (-t723 * t756 - (t767 * t896 + t750) * t773 + (t767 * t819 - t726) * m(3) + t802) * t624 - t699 * t627;
t914 = -t672 * t797 + t687 * t807;
t714 = pkin(2) * t772 + t766 * t783;
t837 = t761 * t772;
t877 = t771 * pkin(3);
t671 = (t758 * t837 + t760 * t766) * t877 + t714 * t858 + t711 * t760;
t840 = t761 * t766;
t704 = t758 * t840 - t760 * t772;
t850 = t759 * t771;
t686 = t765 * t704 + t758 * t850;
t707 = t758 * t772 + t760 * t840;
t843 = t760 * t771;
t689 = -t765 * t707 - t759 * t843;
t762 = legFrame(3,2);
t743 = sin(t762);
t746 = cos(t762);
t795 = t743 * t781 - t746 * t782;
t656 = (t686 * t795 + t689 * t780) * t677;
t869 = t656 * t783;
t812 = t765 * t869;
t668 = (t758 * t766 - t760 * t837) * t877 - t714 * t844 + t758 * t711;
t650 = (t668 * t780 + t671 * t795) * t788 * t677;
t883 = pkin(3) * t650;
t641 = t812 - t883;
t815 = t765 * t883;
t832 = t766 * t771;
t849 = t759 * t772;
t856 = t759 * t765;
t870 = t656 * t772;
t623 = (((t761 * t650 + t656 * t849) * t880 + ((-t815 + t869) * t766 + pkin(2) * t870) * t850 + t641 * t761) * t656 + (t650 * t849 + (t755 * t761 - t832 * t856 - t761) * t656) * t883) * t677;
t701 = pkin(3) * t832 + t711;
t873 = (-t783 * t815 + (t755 * t787 + t771 * t827 + t740) * t656) * t656;
t626 = t677 * t834 * t873 + (-t761 * t812 + (-t701 * t856 + t761 * (pkin(2) * t771 + t880)) * t650) / (t701 * t850 + (pkin(2) + t877) * t833) * t650;
t632 = (-t771 * t873 - (pkin(2) * t650 - t641 * t771) * t883) * t677;
t653 = t656 ^ 2;
t698 = t727 * t771 + t728 * t765;
t720 = t765 * rSges(3,1) + t771 * rSges(3,2);
t801 = rSges(3,1) * t771 - rSges(3,2) * t765;
t890 = m(3) * (-t720 * t855 + t761 * t801);
t900 = t739 * t755 + t765 * t809;
t798 = (t698 * t623 + t724 * t626 - t632 * t890 + 0.2e1 * t653 * ((t710 * t765 + t742) * t771 + t820 + t900)) * t788;
t910 = m(3) * t801;
t864 = ((t735 + t910) * t772 - t766 * t725) * t759;
t808 = 0.4e1 * ((t893 * t765 + t892 * t771) * t650 + ((t765 * t894 + t742) * t771 + t891 + t900) * t656) * t650 + t632 * t864 - (-t723 * t755 - (t765 * t896 + t750) * t771 + (t765 * t819 - t726) * m(3) + t802) * t623 - t698 * t626;
t913 = -t671 * t798 + t686 * t808;
t887 = m(3) * t761;
t886 = pkin(2) * t765;
t885 = pkin(2) * t767;
t884 = pkin(2) * t769;
t857 = t759 * t760;
t647 = t650 ^ 2;
t754 = -m(1) - m(2) - m(3);
t818 = 0.2e1 * m(3);
t826 = -t623 * t864 - t626 * t890 + t754 * t632 + ((-t653 * t735 - (t653 + t647) * t910) * t766 - (t650 * t720 * t818 + t656 * t725) * t870) * t759 - t647 * t720 * t887;
t648 = t651 ^ 2;
t825 = -t624 * t863 - t627 * t889 + t754 * t633 + ((-t654 * t735 - (t654 + t648) * t911) * t768 - (t651 * t721 * t818 + t657 * t725) * t868) * t759 - t648 * t721 * t887;
t649 = t652 ^ 2;
t824 = -t625 * t862 - t628 * t888 + t754 * t634 + ((-t655 * t735 - (t655 + t649) * t912) * t770 - (t652 * t722 * t818 + t658 * t725) * t866) * t759 - t649 * t722 * t887;
t792 = pkin(3) * t856 - t711 * t761;
t791 = pkin(3) * t854 - t712 * t761;
t790 = pkin(3) * t852 - t713 * t761;
t682 = -t758 * t716 + t760 * t790;
t681 = -t758 * t715 + t760 * t791;
t680 = -t758 * t714 + t760 * t792;
t1 = [(t824 * ((t709 * t748 + t745 * t851) * t878 + (-t682 * t748 + t697 * t745) * t775 + (t761 * t745 - t748 * t857) * t884) + t915 * t748) * t679 + (t825 * ((t708 * t747 + t744 * t853) * t879 + (-t681 * t747 + t696 * t744) * t773 + (t761 * t744 - t747 * t857) * t885) + t914 * t747) * t678 + (t826 * ((t707 * t746 + t743 * t855) * t880 + (-t680 * t746 + t695 * t743) * t771 + (t761 * t743 - t746 * t857) * t886) + t913 * t746) * t677; (t824 * (-(t709 * t745 - t748 * t851) * t878 + (t682 * t745 + t748 * t697) * t775 + (t745 * t857 + t748 * t761) * t884) - t915 * t745) * t679 + (t825 * (-(t708 * t744 - t747 * t853) * t879 + (t681 * t744 + t747 * t696) * t773 + (t744 * t857 + t747 * t761) * t885) - t914 * t744) * t678 + (t826 * (-(t707 * t743 - t746 * t855) * t880 + (t680 * t743 + t746 * t695) * t771 + (t743 * t857 + t746 * t761) * t886) - t913 * t743) * t677; (-t806 * t691 + t670 * t796 + t824 * (-t706 * t878 + t716 * t841 + (pkin(2) * t852 + t775 * t790) * t758)) * t679 + (-t807 * t690 + t669 * t797 + t825 * (-t705 * t879 + t715 * t842 + (pkin(2) * t854 + t773 * t791) * t758)) * t678 + (-t808 * t689 + t668 * t798 + t826 * (-t704 * t880 + t714 * t843 + (pkin(2) * t856 + t771 * t792) * t758)) * t677;];
taucX  = t1;
