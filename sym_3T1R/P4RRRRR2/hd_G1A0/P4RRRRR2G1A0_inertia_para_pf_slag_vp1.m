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
% Datum: 2020-08-07 17:26
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P4RRRRR2G1A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4RRRRR2G1A0_inertia_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4RRRRR2G1A0_inertia_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4RRRRR2G1A0_inertia_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4RRRRR2G1A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4RRRRR2G1A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P4RRRRR2G1A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4RRRRR2G1A0_inertia_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4RRRRR2G1A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 17:24:15
% EndTime: 2020-08-07 17:24:17
% DurationCPUTime: 2.15s
% Computational Cost: add. (8470->337), mult. (12475->550), div. (1508->18), fcn. (3820->66), ass. (0->255)
t937 = m(3) / 0.2e1;
t947 = Icges(3,2) / 0.2e1;
t936 = m(3) * rSges(3,2);
t785 = -rSges(3,1) * t936 + Icges(3,4);
t834 = 2 * qJ(3,4);
t846 = rSges(3,2) ^ 2;
t847 = rSges(3,1) ^ 2;
t934 = (-t846 + t847) * t937 - Icges(3,1) / 0.2e1 + t947;
t946 = cos(t834) * t934 + t785 * sin(t834);
t835 = 2 * qJ(3,3);
t945 = cos(t835) * t934 + t785 * sin(t835);
t836 = 2 * qJ(3,2);
t944 = cos(t836) * t934 + t785 * sin(t836);
t837 = 2 * qJ(3,1);
t943 = cos(t837) * t934 + t785 * sin(t837);
t942 = -2 * pkin(1);
t811 = cos(qJ(3,4));
t794 = 0.1e1 / t811;
t941 = t794 ^ 2;
t819 = cos(qJ(3,3));
t799 = 0.1e1 / t819;
t940 = t799 ^ 2;
t821 = cos(qJ(3,2));
t801 = 0.1e1 / t821;
t939 = t801 ^ 2;
t823 = cos(qJ(3,1));
t803 = 0.1e1 / t823;
t938 = t803 ^ 2;
t935 = m(3) * rSges(3,3);
t784 = rSges(3,1) * t935 - Icges(3,5);
t809 = sin(qJ(3,4));
t810 = sin(qJ(2,4));
t891 = m(3) * rSges(3,1) * pkin(1);
t710 = (Icges(3,6) + (-pkin(1) * t810 - rSges(3,3)) * t936) * t811 - t809 * (t810 * t891 + t784);
t933 = t710 * t794;
t813 = sin(qJ(3,3));
t814 = sin(qJ(2,3));
t711 = (Icges(3,6) + (-pkin(1) * t814 - rSges(3,3)) * t936) * t819 - t813 * (t814 * t891 + t784);
t932 = t711 * t799;
t815 = sin(qJ(3,2));
t816 = sin(qJ(2,2));
t712 = (Icges(3,6) + (-pkin(1) * t816 - rSges(3,3)) * t936) * t821 - t815 * (t816 * t891 + t784);
t931 = t712 * t801;
t817 = sin(qJ(3,1));
t818 = sin(qJ(2,1));
t713 = (Icges(3,6) + (-pkin(1) * t818 - rSges(3,3)) * t936) * t823 - t817 * (t818 * t891 + t784);
t930 = t713 * t803;
t786 = qJ(1,4) + legFrame(4,3);
t778 = qJ(2,4) + t786;
t763 = qJ(3,4) + t778;
t764 = -qJ(3,4) + t778;
t714 = sin(t786) * t942 + (-sin(t764) - sin(t763)) * pkin(2);
t748 = 0.1e1 / (sin(qJ(2,4) + qJ(3,4)) + sin(qJ(2,4) - qJ(3,4)));
t929 = t714 * t748;
t715 = cos(t786) * t942 + (-cos(t763) - cos(t764)) * pkin(2);
t928 = t715 * t748;
t787 = qJ(1,3) + legFrame(3,3);
t779 = qJ(2,3) + t787;
t772 = qJ(3,3) + t779;
t773 = -qJ(3,3) + t779;
t716 = sin(t787) * t942 + (-sin(t773) - sin(t772)) * pkin(2);
t750 = 0.1e1 / (sin(qJ(2,3) + qJ(3,3)) + sin(qJ(2,3) - qJ(3,3)));
t927 = t716 * t750;
t788 = qJ(1,2) + legFrame(2,3);
t780 = qJ(2,2) + t788;
t774 = qJ(3,2) + t780;
t775 = -qJ(3,2) + t780;
t717 = sin(t788) * t942 + (-sin(t775) - sin(t774)) * pkin(2);
t751 = 0.1e1 / (sin(qJ(2,2) + qJ(3,2)) + sin(qJ(2,2) - qJ(3,2)));
t926 = t717 * t751;
t789 = qJ(1,1) + legFrame(1,3);
t781 = qJ(2,1) + t789;
t776 = qJ(3,1) + t781;
t777 = -qJ(3,1) + t781;
t718 = sin(t789) * t942 + (-sin(t777) - sin(t776)) * pkin(2);
t752 = 0.1e1 / (sin(qJ(2,1) + qJ(3,1)) + sin(qJ(2,1) - qJ(3,1)));
t925 = t718 * t752;
t719 = cos(t787) * t942 + (-cos(t772) - cos(t773)) * pkin(2);
t924 = t719 * t750;
t720 = cos(t788) * t942 + (-cos(t774) - cos(t775)) * pkin(2);
t923 = t720 * t751;
t721 = cos(t789) * t942 + (-cos(t776) - cos(t777)) * pkin(2);
t922 = t721 * t752;
t848 = 0.1e1 / pkin(2);
t921 = t748 * t848;
t920 = t750 * t848;
t919 = t751 * t848;
t918 = t752 * t848;
t812 = cos(qJ(2,4));
t917 = (t812 * pkin(1) + t811 * pkin(2)) / t811 ^ 2;
t820 = cos(qJ(2,3));
t916 = (t820 * pkin(1) + t819 * pkin(2)) / t819 ^ 2;
t822 = cos(qJ(2,2));
t915 = (t822 * pkin(1) + t821 * pkin(2)) / t821 ^ 2;
t824 = cos(qJ(2,1));
t914 = (t824 * pkin(1) + t823 * pkin(2)) / t823 ^ 2;
t761 = sin(t778);
t793 = 0.1e1 / t810;
t913 = t761 * t793;
t762 = cos(t778);
t912 = t762 * t793;
t766 = sin(t779);
t796 = 0.1e1 / t814;
t911 = t766 * t796;
t767 = sin(t780);
t797 = 0.1e1 / t816;
t910 = t767 * t797;
t768 = sin(t781);
t798 = 0.1e1 / t818;
t909 = t768 * t798;
t769 = cos(t779);
t908 = t769 * t796;
t770 = cos(t780);
t907 = t770 * t797;
t771 = cos(t781);
t906 = t771 * t798;
t905 = t793 * t809;
t850 = 1 / pkin(1);
t904 = t793 * t850;
t903 = t794 * t848;
t902 = t796 * t813;
t901 = t796 * t850;
t900 = t797 * t815;
t899 = t797 * t850;
t898 = t798 * t817;
t897 = t798 * t850;
t896 = t799 * t848;
t895 = t801 * t848;
t894 = t803 * t848;
t893 = t848 * t850;
t892 = t846 + t847;
t890 = t714 * t921;
t889 = t715 * t921;
t888 = t716 * t920;
t887 = t717 * t919;
t886 = t718 * t918;
t885 = t719 * t920;
t884 = t720 * t919;
t883 = t721 * t918;
t783 = -rSges(3,2) * t935 + Icges(3,6);
t722 = t783 * t811 - t784 * t809;
t882 = t722 * t748 * t794;
t723 = t783 * t819 - t784 * t813;
t881 = t723 * t750 * t799;
t724 = t783 * t821 - t784 * t815;
t880 = t724 * t751 * t801;
t725 = t783 * t823 - t784 * t817;
t879 = t725 * t752 * t803;
t878 = t809 * t917;
t877 = t848 * t917;
t876 = t813 * t916;
t875 = t848 * t916;
t874 = t815 * t915;
t873 = t848 * t915;
t872 = t817 * t914;
t871 = t848 * t914;
t870 = t794 * t905;
t869 = t809 * t904;
t868 = t722 * t903;
t867 = t799 * t902;
t866 = t813 * t901;
t865 = t801 * t900;
t864 = t815 * t899;
t863 = t803 * t898;
t862 = t817 * t897;
t861 = t723 * t896;
t860 = t724 * t895;
t859 = t725 * t894;
t858 = 0.2e1 * rSges(3,3) ^ 2 + t892;
t857 = Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + t947 + Icges(3,1) / 0.2e1;
t856 = t858 * t937 + t857;
t782 = m(2) * rSges(2,2) - t935;
t828 = m(2) * rSges(2,1);
t855 = (-(-t828 + (-rSges(3,1) * t811 + rSges(3,2) * t809) * m(3)) * t812 - t782 * t810) * pkin(1);
t854 = (-(-t828 + (-rSges(3,1) * t819 + rSges(3,2) * t813) * m(3)) * t820 - t782 * t814) * pkin(1);
t853 = (-(-t828 + (-rSges(3,1) * t821 + rSges(3,2) * t815) * m(3)) * t822 - t782 * t816) * pkin(1);
t852 = (-(-t828 + (-rSges(3,1) * t823 + rSges(3,2) * t817) * m(3)) * t824 - t782 * t818) * pkin(1);
t706 = t856 + t946;
t707 = t856 + t945;
t708 = t856 + t944;
t709 = t856 + t943;
t849 = pkin(1) ^ 2;
t851 = Icges(1,3) + ((2 * t849) + t858) * t937 + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + m(2) * t849 + t857;
t845 = koppelP(1,1);
t844 = koppelP(2,1);
t843 = koppelP(3,1);
t842 = koppelP(4,1);
t841 = koppelP(1,2);
t840 = koppelP(2,2);
t839 = koppelP(3,2);
t838 = koppelP(4,2);
t833 = rSges(4,1);
t832 = rSges(4,2);
t831 = xP(4);
t791 = cos(t831);
t790 = sin(t831);
t743 = -t790 * t841 + t791 * t845;
t742 = -t790 * t840 + t791 * t844;
t741 = -t790 * t839 + t791 * t843;
t740 = -t790 * t838 + t791 * t842;
t739 = -t790 * t845 - t791 * t841;
t738 = -t790 * t844 - t791 * t840;
t737 = -t790 * t843 - t791 * t839;
t736 = -t790 * t842 - t791 * t838;
t727 = m(4) * (-t790 * t832 + t791 * t833);
t726 = m(4) * (-t790 * t833 - t791 * t832);
t705 = (t739 * t771 + t743 * t768) * t897;
t704 = (t738 * t770 + t742 * t767) * t899;
t703 = (t737 * t769 + t741 * t766) * t901;
t702 = (t736 * t762 + t740 * t761) * t904;
t701 = t852 + t709;
t700 = t853 + t708;
t699 = t854 + t707;
t698 = t855 + t706;
t697 = t851 + 0.2e1 * t852 + t943;
t696 = t851 + 0.2e1 * t853 + t944;
t695 = t851 + 0.2e1 * t854 + t945;
t694 = t851 + 0.2e1 * t855 + t946;
t693 = (t718 * t743 + t721 * t739) * t752 * t893;
t692 = (t717 * t742 + t720 * t738) * t751 * t893;
t691 = (t716 * t741 + t719 * t737) * t750 * t893;
t690 = (t714 * t740 + t715 * t736) * t748 * t893;
t689 = t859 + (t701 * t803 - t709 * t871) * t862;
t688 = t860 + (t700 * t801 - t708 * t873) * t864;
t687 = t861 + (t699 * t799 - t707 * t875) * t866;
t686 = (t701 * t906 + t709 * t883) * t850;
t685 = (t700 * t907 + t708 * t884) * t850;
t684 = (t699 * t908 + t707 * t885) * t850;
t683 = (t701 * t909 + t709 * t886) * t850;
t682 = (t700 * t910 + t708 * t887) * t850;
t681 = (t699 * t911 + t707 * t888) * t850;
t680 = t868 + (t698 * t794 - t706 * t877) * t869;
t679 = (t698 * t912 + t706 * t889) * t850;
t678 = (t698 * t913 + t706 * t890) * t850;
t677 = (t697 * t906 + t701 * t883) * t850;
t676 = (t696 * t907 + t700 * t884) * t850;
t675 = (t695 * t908 + t699 * t885) * t850;
t674 = (t697 * t909 + t701 * t886) * t850;
t673 = (t696 * t910 + t700 * t887) * t850;
t672 = (t695 * t911 + t699 * t888) * t850;
t671 = (t694 * t912 + t698 * t889) * t850;
t670 = (t694 * t913 + t698 * t890) * t850;
t669 = t713 * t894 + (t697 * t803 - t701 * t871) * t862;
t668 = t712 * t895 + (t696 * t801 - t700 * t873) * t864;
t667 = t711 * t896 + (t695 * t799 - t699 * t875) * t866;
t666 = t710 * t903 + (t694 * t794 - t698 * t877) * t869;
t665 = t693 * t709 + t705 * t701;
t664 = t692 * t708 + t704 * t700;
t663 = t691 * t707 + t703 * t699;
t662 = t690 * t706 + t702 * t698;
t661 = t693 * t701 + t705 * t697;
t660 = t692 * t700 + t704 * t696;
t659 = t691 * t699 + t703 * t695;
t658 = t690 * t698 + t702 * t694;
t1 = [m(4) + (t671 * t912 + t675 * t908 + t676 * t907 + t677 * t906 + (t679 * t928 + t684 * t924 + t685 * t923 + t686 * t922) * t848) * t850, (t671 * t913 + t675 * t911 + t676 * t910 + t677 * t909 + (t679 * t929 + t684 * t927 + t685 * t926 + t686 * t925) * t848) * t850, (t671 * t870 + t675 * t867 + t676 * t865 + t677 * t863 + ((t715 * t882 + t719 * t881 + t720 * t880 + t721 * t879) * t848 + (-t686 * t872 + t771 * t930) * t798 + (-t685 * t874 + t770 * t931) * t797 + (-t684 * t876 + t769 * t932) * t796 + (-t679 * t878 + t762 * t933) * t793) * t848) * t850, t671 * t702 + t675 * t703 + t676 * t704 + t677 * t705 + t679 * t690 + t684 * t691 + t685 * t692 + t686 * t693 + t726; (t670 * t912 + t672 * t908 + t673 * t907 + t674 * t906 + (t678 * t928 + t681 * t924 + t682 * t923 + t683 * t922) * t848) * t850, m(4) + (t670 * t913 + t672 * t911 + t673 * t910 + t674 * t909 + (t678 * t929 + t681 * t927 + t682 * t926 + t683 * t925) * t848) * t850, (t670 * t870 + t672 * t867 + t673 * t865 + t674 * t863 + ((t714 * t882 + t716 * t881 + t717 * t880 + t718 * t879) * t848 + (-t683 * t872 + t768 * t930) * t798 + (-t682 * t874 + t767 * t931) * t797 + (-t681 * t876 + t766 * t932) * t796 + (-t678 * t878 + t761 * t933) * t793) * t848) * t850, t670 * t702 + t672 * t703 + t673 * t704 + t674 * t705 + t678 * t690 + t681 * t691 + t682 * t692 + t683 * t693 + t727; (t666 * t912 + t667 * t908 + t668 * t907 + t669 * t906 + (t680 * t928 + t687 * t924 + t688 * t923 + t689 * t922) * t848) * t850, (t666 * t913 + t667 * t911 + t668 * t910 + t669 * t909 + (t680 * t929 + t687 * t927 + t688 * t926 + t689 * t925) * t848) * t850, m(4) + (t666 * t870 + t667 * t867 + t668 * t865 + t669 * t863) * t850 + ((t938 + t939 + t940 + t941) * t848 * (t892 * m(3) + Icges(3,3)) + ((t713 * t938 + (-t689 - t859) * t914) * t898 + (t712 * t939 + (-t688 - t860) * t915) * t900 + (t711 * t940 + (-t687 - t861) * t916) * t902 + (t710 * t941 + (-t680 - t868) * t917) * t905) * t850) * t848, t666 * t702 + t667 * t703 + t668 * t704 + t669 * t705 + t680 * t690 + t687 * t691 + t688 * t692 + t689 * t693; t726 + (t658 * t912 + t659 * t908 + t660 * t907 + t661 * t906 + (t662 * t928 + t663 * t924 + t664 * t923 + t665 * t922) * t848) * t850, t727 + (t658 * t913 + t659 * t911 + t660 * t910 + t661 * t909 + (t662 * t929 + t663 * t927 + t664 * t926 + t665 * t925) * t848) * t850, (t658 * t870 + t659 * t867 + t660 * t865 + t661 * t863) * t850 + ((t693 * t725 + t705 * t713) * t803 + (t692 * t724 + t704 * t712) * t801 + (t691 * t723 + t703 * t711) * t799 + (t690 * t722 + t702 * t710) * t794 + (-t662 * t793 * t878 - t663 * t796 * t876 - t664 * t797 * t874 - t665 * t798 * t872) * t850) * t848, t661 * t705 + t665 * t693 + t660 * t704 + t664 * t692 + t659 * t703 + t663 * t691 + t658 * t702 + t662 * t690 + Icges(4,3) + m(4) * (t832 ^ 2 + t833 ^ 2);];
MX  = t1;
