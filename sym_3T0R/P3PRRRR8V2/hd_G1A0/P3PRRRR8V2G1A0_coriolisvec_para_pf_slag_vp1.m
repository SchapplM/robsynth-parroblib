% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRRR8V2G1A0
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
% Datum: 2020-08-06 17:36
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G1A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:35:16
% EndTime: 2020-08-06 17:35:20
% DurationCPUTime: 4.33s
% Computational Cost: add. (22602->311), mult. (44376->571), div. (2541->10), fcn. (46854->22), ass. (0->244)
t787 = cos(qJ(2,3));
t797 = pkin(7) + pkin(6);
t751 = t787 * t797;
t781 = sin(qJ(2,3));
t723 = pkin(2) * t781 - t751;
t774 = sin(pkin(4));
t786 = cos(qJ(3,3));
t776 = cos(pkin(4));
t780 = sin(qJ(3,3));
t854 = t776 * t780;
t862 = t774 * t781;
t770 = t786 ^ 2;
t891 = pkin(3) * t770;
t680 = 0.1e1 / ((pkin(3) * t854 + t723 * t774) * t786 + pkin(2) * t854 + t862 * t891);
t789 = cos(qJ(2,2));
t752 = t789 * t797;
t783 = sin(qJ(2,2));
t724 = pkin(2) * t783 - t752;
t788 = cos(qJ(3,2));
t782 = sin(qJ(3,2));
t852 = t776 * t782;
t860 = t774 * t783;
t771 = t788 ^ 2;
t890 = pkin(3) * t771;
t681 = 0.1e1 / ((pkin(3) * t852 + t724 * t774) * t788 + pkin(2) * t852 + t860 * t890);
t791 = cos(qJ(2,1));
t753 = t791 * t797;
t785 = sin(qJ(2,1));
t725 = pkin(2) * t785 - t753;
t790 = cos(qJ(3,1));
t784 = sin(qJ(3,1));
t850 = t776 * t784;
t858 = t774 * t785;
t772 = t790 ^ 2;
t889 = pkin(3) * t772;
t682 = 0.1e1 / ((pkin(3) * t850 + t725 * t774) * t790 + pkin(2) * t850 + t858 * t889);
t810 = rSges(3,1) * t790 - rSges(3,2) * t784;
t917 = t810 * m(3);
t811 = rSges(3,1) * t788 - rSges(3,2) * t782;
t916 = t811 * m(3);
t812 = rSges(3,1) * t786 - rSges(3,2) * t780;
t915 = t812 * m(3);
t777 = legFrame(3,3);
t758 = sin(t777);
t761 = cos(t777);
t773 = sin(pkin(8));
t775 = cos(pkin(8));
t701 = -t773 * t758 + t761 * t775;
t704 = t775 * t758 + t761 * t773;
t853 = t776 * t781;
t857 = t774 * t786;
t665 = -t701 * t857 - (t701 * t853 + t787 * t704) * t780;
t668 = -t704 * t857 - (-t787 * t701 + t704 * t853) * t780;
t795 = xDP(2);
t796 = xDP(1);
t650 = (t665 * t796 + t668 * t795) * t680;
t881 = t650 * t797;
t820 = t780 * t881;
t748 = t786 * pkin(3) + pkin(2);
t716 = t781 * t748 - t751;
t846 = t781 * t797;
t872 = (t748 * t787 + t846) * t776;
t671 = -t716 * t701 - t704 * t872;
t674 = -t701 * t872 + t716 * t704;
t719 = t748 * t854;
t686 = 0.1e1 / (t716 * t857 + t719);
t802 = 0.1e1 / pkin(3);
t656 = (t671 * t795 + t674 * t796) * t802 * t686;
t894 = pkin(3) * t656;
t641 = t820 - t894;
t803 = pkin(2) ^ 2;
t755 = t797 ^ 2 + t803;
t801 = pkin(3) ^ 2;
t823 = t780 * t894;
t902 = 0.2e1 * pkin(2);
t838 = pkin(3) * t902;
t884 = (-t797 * t823 + (t770 * t801 + t786 * t838 + t755) * t650) * t650;
t629 = (-t786 * t884 - (pkin(2) * t656 - t641 * t786) * t894) * t680;
t647 = t650 ^ 2;
t653 = t656 ^ 2;
t732 = t780 * rSges(3,1) + t786 * rSges(3,2);
t792 = pkin(6) + rSges(3,3);
t888 = t792 * m(3);
t737 = m(2) * rSges(2,2) - t888;
t798 = pkin(2) * m(3);
t747 = m(2) * rSges(2,1) + t798;
t769 = -m(1) - m(2) - m(3);
t829 = 0.2e1 * m(3);
t841 = t787 * t650;
t911 = t769 * t629 + ((-t647 * t747 - (t647 + t653) * t915) * t781 - (t732 * t656 * t829 + t737 * t650) * t841) * t774;
t778 = legFrame(2,3);
t759 = sin(t778);
t762 = cos(t778);
t702 = -t773 * t759 + t762 * t775;
t705 = t775 * t759 + t762 * t773;
t851 = t776 * t783;
t856 = t774 * t788;
t666 = -t702 * t856 - (t702 * t851 + t789 * t705) * t782;
t669 = -t705 * t856 - (-t789 * t702 + t705 * t851) * t782;
t651 = (t666 * t796 + t669 * t795) * t681;
t880 = t651 * t797;
t819 = t782 * t880;
t749 = t788 * pkin(3) + pkin(2);
t717 = t783 * t749 - t752;
t844 = t783 * t797;
t871 = (t749 * t789 + t844) * t776;
t672 = -t717 * t702 - t705 * t871;
t675 = -t702 * t871 + t717 * t705;
t720 = t749 * t852;
t687 = 0.1e1 / (t717 * t856 + t720);
t657 = (t672 * t795 + t675 * t796) * t802 * t687;
t893 = pkin(3) * t657;
t642 = t819 - t893;
t822 = t782 * t893;
t883 = (-t797 * t822 + (t771 * t801 + t788 * t838 + t755) * t651) * t651;
t630 = (-t788 * t883 - (pkin(2) * t657 - t642 * t788) * t893) * t681;
t648 = t651 ^ 2;
t654 = t657 ^ 2;
t733 = t782 * rSges(3,1) + t788 * rSges(3,2);
t840 = t789 * t651;
t910 = t769 * t630 + ((-t648 * t747 - (t648 + t654) * t916) * t783 - (t733 * t657 * t829 + t737 * t651) * t840) * t774;
t779 = legFrame(1,3);
t760 = sin(t779);
t763 = cos(t779);
t703 = -t773 * t760 + t763 * t775;
t706 = t775 * t760 + t763 * t773;
t849 = t776 * t785;
t855 = t774 * t790;
t667 = -t703 * t855 - (t703 * t849 + t791 * t706) * t784;
t670 = -t706 * t855 - (-t791 * t703 + t706 * t849) * t784;
t652 = (t667 * t796 + t670 * t795) * t682;
t879 = t652 * t797;
t818 = t784 * t879;
t750 = t790 * pkin(3) + pkin(2);
t718 = t785 * t750 - t753;
t842 = t785 * t797;
t870 = (t750 * t791 + t842) * t776;
t673 = -t718 * t703 - t706 * t870;
t676 = -t703 * t870 + t718 * t706;
t721 = t750 * t850;
t688 = 0.1e1 / (t718 * t855 + t721);
t658 = (t673 * t795 + t676 * t796) * t802 * t688;
t892 = pkin(3) * t658;
t643 = t818 - t892;
t821 = t784 * t892;
t882 = (-t797 * t821 + (t772 * t801 + t790 * t838 + t755) * t652) * t652;
t631 = (-t790 * t882 - (pkin(2) * t658 - t643 * t790) * t892) * t682;
t649 = t652 ^ 2;
t655 = t658 ^ 2;
t734 = t784 * rSges(3,1) + t790 * rSges(3,2);
t839 = t791 * t652;
t909 = t769 * t631 + ((-t649 * t747 - (t649 + t655) * t917) * t785 - (t734 * t658 * t829 + t737 * t652) * t839) * t774;
t728 = pkin(2) * t791 + t842;
t859 = t774 * t784;
t804 = pkin(3) * t859 - t725 * t776;
t908 = t728 * t775 + t804 * t773;
t727 = pkin(2) * t789 + t844;
t861 = t774 * t782;
t805 = pkin(3) * t861 - t724 * t776;
t907 = t727 * t775 + t805 * t773;
t726 = pkin(2) * t787 + t846;
t863 = t774 * t780;
t806 = pkin(3) * t863 - t723 * t776;
t906 = t726 * t775 + t806 * t773;
t900 = m(3) * rSges(3,1);
t828 = rSges(3,2) * t900;
t754 = -Icges(3,4) + t828;
t827 = t798 / 0.2e1;
t817 = rSges(3,1) * t827;
t905 = t754 * t770 + t780 * t817;
t904 = t754 * t771 + t782 * t817;
t903 = t754 * t772 + t784 * t817;
t901 = -0.2e1 * t754;
t799 = rSges(3,2) ^ 2;
t800 = rSges(3,1) ^ 2;
t735 = (-t799 + t800) * m(3) + Icges(3,2) - Icges(3,1);
t899 = t735 / 0.2e1;
t739 = rSges(3,2) * t888 - Icges(3,6);
t898 = -t739 / 0.4e1;
t740 = rSges(3,1) * t888 - Icges(3,5);
t897 = t740 / 0.4e1;
t896 = -t754 / 0.2e1;
t895 = m(3) * t776;
t847 = t781 * t786;
t698 = pkin(3) * t847 + t723;
t848 = t776 * t802;
t620 = t680 * t848 * t884 + (-t776 * t820 + (-t698 * t863 + (pkin(2) * t786 + t891) * t776) * t656) / (t698 * t857 + t719) * t656;
t692 = -t732 * t862 + t812 * t776;
t887 = t620 * t692;
t845 = t783 * t788;
t699 = pkin(3) * t845 + t724;
t621 = t681 * t848 * t883 + (-t776 * t819 + (-t699 * t861 + (pkin(2) * t788 + t890) * t776) * t657) / (t699 * t856 + t720) * t657;
t693 = -t733 * t860 + t811 * t776;
t886 = t621 * t693;
t843 = t785 * t790;
t700 = pkin(3) * t843 + t725;
t622 = t682 * t848 * t882 + (-t776 * t818 + (-t700 * t859 + (pkin(2) * t790 + t889) * t776) * t658) / (t700 * t855 + t721) * t658;
t694 = -t734 * t858 + t810 * t776;
t885 = t622 * t694;
t878 = t653 * t732;
t877 = t654 * t733;
t876 = t655 * t734;
t689 = (t747 + t915) * t787 - t737 * t781;
t875 = t689 * t774;
t690 = (t747 + t916) * t789 - t737 * t783;
t874 = t690 * t774;
t691 = (t747 + t917) * t791 - t737 * t785;
t873 = t691 * t774;
t617 = (((t776 * t656 + t774 * t841) * t891 + ((-t823 + t881) * t781 + pkin(2) * t841) * t857 + t776 * t641) * t650 + (t787 * t774 * t656 + (t770 * t776 - t847 * t863 - t776) * t650) * t894) * t680;
t837 = -m(3) * t887 - t617 * t875 - t878 * t895 + t911;
t618 = (((t776 * t657 + t774 * t840) * t890 + ((-t822 + t880) * t783 + pkin(2) * t840) * t856 + t776 * t642) * t651 + (t789 * t774 * t657 + (t771 * t776 - t845 * t861 - t776) * t651) * t893) * t681;
t836 = -m(3) * t886 - t618 * t874 - t877 * t895 + t910;
t619 = (((t776 * t658 + t774 * t839) * t889 + ((-t821 + t879) * t785 + pkin(2) * t839) * t855 + t776 * t643) * t652 + (t658 * t774 * t791 + (t772 * t776 - t843 * t859 - t776) * t652) * t892) * t682;
t835 = -m(3) * t885 - t619 * t873 - t876 * t895 + t909;
t831 = -t828 / 0.2e1 + Icges(3,4) / 0.2e1;
t830 = rSges(3,2) * t902;
t826 = pkin(2) * t863;
t825 = pkin(2) * t861;
t824 = pkin(2) * t859;
t757 = rSges(3,2) * t827;
t695 = t739 * t786 + t740 * t780;
t738 = t792 ^ 2 + t799 + t803;
t765 = t900 * t902;
t813 = -(rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) - Icges(3,1) - Icges(2,3);
t816 = -0.4e1 * ((t898 * t780 + t897 * t786) * t656 + ((t780 * t899 + t757) * t786 + t896 + t905) * t650) * t656 - t629 * t875 + (-t735 * t770 - (t780 * t901 + t765) * t786 + (t780 * t830 - t738) * m(3) + t813) * t617 + t695 * t620;
t696 = t739 * t788 + t740 * t782;
t815 = -0.4e1 * ((t898 * t782 + t897 * t788) * t657 + ((t782 * t899 + t757) * t788 + t896 + t904) * t651) * t657 - t630 * t874 + (-t735 * t771 - (t782 * t901 + t765) * t788 + (t782 * t830 - t738) * m(3) + t813) * t618 + t696 * t621;
t697 = t739 * t790 + t740 * t784;
t814 = -0.4e1 * ((t898 * t784 + t897 * t790) * t658 + ((t784 * t899 + t757) * t790 + t896 + t903) * t652) * t658 - t631 * t873 + (-t735 * t772 - (t784 * t901 + t765) * t790 + (t784 * t830 - t738) * m(3) + t813) * t619 + t697 * t622;
t722 = (t800 / 0.2e1 - t799 / 0.2e1) * m(3) - Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1;
t736 = -(t799 + t800) * m(3) - Icges(3,3);
t809 = (-t692 * m(3) * t629 + t695 * t617 + t736 * t620 + 0.2e1 * t647 * ((t722 * t780 + t757) * t786 + t831 + t905)) * t686;
t808 = (-t693 * m(3) * t630 + t696 * t618 + t736 * t621 + 0.2e1 * t648 * ((t722 * t782 + t757) * t788 + t831 + t904)) * t687;
t807 = (-t694 * m(3) * t631 + t697 * t619 + t736 * t622 + 0.2e1 * t649 * ((t722 * t784 + t757) * t790 + t831 + t903)) * t688;
t712 = t773 * t791 + t775 * t849;
t711 = t773 * t789 + t775 * t851;
t710 = t773 * t787 + t775 * t853;
t709 = t773 * t849 - t775 * t791;
t708 = t773 * t851 - t775 * t789;
t707 = t773 * t853 - t775 * t787;
t685 = t773 * t728 - t804 * t775;
t684 = t773 * t727 - t805 * t775;
t683 = t773 * t726 - t806 * t775;
t1 = [(t814 * t667 + t835 * (-(t709 * t763 + t760 * t712) * t889 + (-t685 * t760 + t908 * t763) * t790 + t706 * t824)) * t682 + (t815 * t666 + t836 * (-(t708 * t762 + t759 * t711) * t890 + (-t684 * t759 + t907 * t762) * t788 + t705 * t825)) * t681 + (t816 * t665 + t837 * (-(t707 * t761 + t758 * t710) * t891 + (-t683 * t758 + t906 * t761) * t786 + t704 * t826)) * t680 + (t674 * t809 + t675 * t808 + t676 * t807) * t802; (t814 * t670 + t835 * ((-t760 * t709 + t712 * t763) * t889 + (t685 * t763 + t908 * t760) * t790 - t703 * t824)) * t682 + (t815 * t669 + t836 * ((-t759 * t708 + t711 * t762) * t890 + (t684 * t762 + t907 * t759) * t788 - t702 * t825)) * t681 + (t816 * t668 + t837 * ((-t758 * t707 + t710 * t761) * t891 + (t683 * t761 + t906 * t758) * t786 - t701 * t826)) * t680 + (t671 * t809 + t672 * t808 + t673 * t807) * t802; (-t617 * t689 - t618 * t690 - t619 * t691) * t774 + (-t887 - t886 - t885 + (-t876 - t877 - t878) * t776) * m(3) + t909 + t910 + t911;];
taucX  = t1;
