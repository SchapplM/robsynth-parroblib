% Calculate vector of centrifugal and coriolis load on the joints for
% P3RRRRR2G2A0
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
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:10
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRRRR2G2A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G2A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRRRR2G2A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G2A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G2A0_coriolisvec_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR2G2A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR2G2A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRRRR2G2A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G2A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G2A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:08:50
% EndTime: 2020-03-09 21:08:53
% DurationCPUTime: 3.20s
% Computational Cost: add. (11961->264), mult. (23346->496), div. (9408->14), fcn. (26376->39), ass. (0->244)
t784 = legFrame(1,2);
t764 = sin(t784);
t767 = cos(t784);
t800 = cos(qJ(3,1));
t791 = sin(qJ(3,1));
t801 = cos(qJ(2,1));
t924 = t801 * pkin(1);
t855 = t791 * t924;
t792 = sin(qJ(2,1));
t793 = sin(qJ(1,1));
t802 = cos(qJ(1,1));
t743 = t802 * t792 + t793 * t801;
t777 = t800 ^ 2;
t858 = t743 * t777 * pkin(2);
t896 = t767 * t791;
t927 = pkin(1) * t793;
t706 = t764 * t858 + (-pkin(2) * t896 + t764 * t927) * t800 - t767 * t855;
t804 = xDP(2);
t779 = 0.1e1 / t800 ^ 2;
t770 = 0.1e1 / t792;
t808 = 0.1e1 / pkin(2);
t810 = 0.1e1 / pkin(1);
t873 = t808 * t810;
t837 = t770 * t873;
t827 = t779 * t837;
t700 = t706 * t804 * t827;
t899 = t764 * t791;
t707 = -t767 * t858 + (-pkin(2) * t899 - t767 * t927) * t800 - t764 * t855;
t805 = xDP(1);
t701 = t707 * t805 * t827;
t803 = xDP(3);
t836 = t803 * t873;
t871 = qJ(2,1) + qJ(3,1);
t872 = qJ(2,1) - qJ(3,1);
t936 = -0.2e1 * pkin(1);
t905 = (t802 * t936 + (-cos(qJ(1,1) + t872) - cos(qJ(1,1) + t871)) * pkin(2)) / (sin(t871) + sin(t872));
t719 = t836 * t905;
t689 = t701 + t700 + t719;
t902 = t743 * t800;
t724 = -t764 * t902 + t896;
t725 = t767 * t902 + t899;
t761 = cos(qJ(1,1) + qJ(2,1));
t778 = 0.1e1 / t800;
t893 = t770 * t810;
t838 = t778 * t893;
t874 = t803 * t810;
t695 = t761 * t770 * t874 + (t724 * t804 + t725 * t805) * t838;
t861 = -t689 - t695;
t683 = t861 ^ 2;
t783 = legFrame(2,2);
t763 = sin(t783);
t766 = cos(t783);
t797 = cos(qJ(3,2));
t788 = sin(qJ(3,2));
t798 = cos(qJ(2,2));
t925 = t798 * pkin(1);
t856 = t788 * t925;
t789 = sin(qJ(2,2));
t790 = sin(qJ(1,2));
t799 = cos(qJ(1,2));
t742 = t799 * t789 + t790 * t798;
t774 = t797 ^ 2;
t859 = t742 * t774 * pkin(2);
t897 = t766 * t788;
t929 = pkin(1) * t790;
t704 = t763 * t859 + (-pkin(2) * t897 + t763 * t929) * t797 - t766 * t856;
t776 = 0.1e1 / t797 ^ 2;
t769 = 0.1e1 / t789;
t839 = t769 * t873;
t828 = t776 * t839;
t698 = t704 * t804 * t828;
t900 = t763 * t788;
t705 = -t766 * t859 + (-pkin(2) * t900 - t766 * t929) * t797 - t763 * t856;
t699 = t705 * t805 * t828;
t869 = qJ(2,2) + qJ(3,2);
t870 = qJ(2,2) - qJ(3,2);
t906 = (t799 * t936 + (-cos(qJ(1,2) + t870) - cos(qJ(1,2) + t869)) * pkin(2)) / (sin(t869) + sin(t870));
t718 = t836 * t906;
t688 = t699 + t698 + t718;
t903 = t742 * t797;
t722 = -t763 * t903 + t897;
t723 = t766 * t903 + t900;
t760 = cos(qJ(1,2) + qJ(2,2));
t775 = 0.1e1 / t797;
t894 = t769 * t810;
t840 = t775 * t894;
t694 = t760 * t769 * t874 + (t722 * t804 + t723 * t805) * t840;
t862 = -t688 - t694;
t682 = t862 ^ 2;
t782 = legFrame(3,2);
t762 = sin(t782);
t765 = cos(t782);
t794 = cos(qJ(3,3));
t785 = sin(qJ(3,3));
t795 = cos(qJ(2,3));
t926 = t795 * pkin(1);
t857 = t785 * t926;
t786 = sin(qJ(2,3));
t787 = sin(qJ(1,3));
t796 = cos(qJ(1,3));
t741 = t796 * t786 + t787 * t795;
t771 = t794 ^ 2;
t860 = t741 * t771 * pkin(2);
t898 = t765 * t785;
t931 = pkin(1) * t787;
t702 = t762 * t860 + (-pkin(2) * t898 + t762 * t931) * t794 - t765 * t857;
t773 = 0.1e1 / t794 ^ 2;
t768 = 0.1e1 / t786;
t841 = t768 * t873;
t829 = t773 * t841;
t696 = t702 * t804 * t829;
t901 = t762 * t785;
t703 = -t765 * t860 + (-pkin(2) * t901 - t765 * t931) * t794 - t762 * t857;
t697 = t703 * t805 * t829;
t867 = qJ(2,3) + qJ(3,3);
t868 = qJ(2,3) - qJ(3,3);
t907 = (t796 * t936 + (-cos(qJ(1,3) + t868) - cos(qJ(1,3) + t867)) * pkin(2)) / (sin(t867) + sin(t868));
t717 = t836 * t907;
t687 = t697 + t696 + t717;
t904 = t741 * t794;
t720 = -t762 * t904 + t898;
t721 = t765 * t904 + t901;
t759 = cos(qJ(1,3) + qJ(2,3));
t772 = 0.1e1 / t794;
t895 = t768 * t810;
t842 = t772 * t895;
t693 = t759 * t768 * t874 + (t720 * t804 + t721 * t805) * t842;
t863 = -t687 - t693;
t681 = t863 ^ 2;
t940 = 0.2e1 * pkin(1);
t690 = t693 ^ 2;
t939 = t690 * t926;
t691 = t694 ^ 2;
t938 = t691 * t925;
t692 = t695 ^ 2;
t937 = t692 * t924;
t935 = -0.2e1 * t794;
t934 = -0.2e1 * t797;
t933 = -0.2e1 * t800;
t932 = pkin(1) * t786;
t930 = pkin(1) * t789;
t928 = pkin(1) * t792;
t923 = -Ifges(3,1) - Ifges(2,3);
t922 = Ifges(3,4) * t785;
t921 = Ifges(3,4) * t788;
t920 = Ifges(3,4) * t791;
t732 = (t762 * t805 + t765 * t804) * t808 * t772;
t729 = t732 ^ 2;
t919 = (t729 / 0.2e1 + (t693 + t687 / 0.2e1) * t687) * t786;
t733 = (t763 * t805 + t766 * t804) * t808 * t775;
t730 = t733 ^ 2;
t918 = (t730 / 0.2e1 + (t694 + t688 / 0.2e1) * t688) * t789;
t734 = (t764 * t805 + t767 * t804) * t808 * t778;
t731 = t734 ^ 2;
t917 = (t731 / 0.2e1 + (t695 + t689 / 0.2e1) * t689) * t792;
t916 = t863 * t732;
t915 = t863 * t771;
t914 = t862 * t733;
t913 = t862 * t774;
t912 = t861 * t734;
t911 = t861 * t777;
t910 = t729 * t785;
t909 = t730 * t788;
t908 = t731 * t791;
t780 = Ifges(3,2) - Ifges(3,1);
t892 = t780 * t771;
t891 = t780 * t774;
t890 = t780 * t777;
t889 = t780 * t785;
t888 = t780 * t788;
t887 = t780 * t791;
t781 = mrSges(2,2) - mrSges(3,3);
t886 = t781 * t795;
t885 = t781 * t798;
t884 = t781 * t801;
t883 = t785 * t786;
t882 = t788 * t789;
t881 = t791 * t792;
t880 = t794 * t795;
t879 = t795 * t732;
t878 = t797 * t798;
t877 = t798 * t733;
t876 = t800 * t801;
t875 = t801 * t734;
t678 = t697 / 0.2e1 + t696 / 0.2e1 + t717 / 0.2e1 + t693;
t807 = pkin(2) ^ 2;
t845 = t732 * t883;
t654 = (-t807 * t915 + (pkin(1) * t693 + (0.2e1 * t678 * t880 - t845) * pkin(2)) * pkin(1)) * t772 * t693 * t841 + (-pkin(2) * t915 + (-t863 * t880 - t845) * pkin(1)) * t687 * t842 + ((pkin(1) * t863 * t883 + pkin(2) * t732) * t794 + pkin(1) * t879) * t773 * t732 * t895;
t663 = (-t939 + (-t794 * t681 - t729 * t772) * pkin(2)) * t895;
t756 = mrSges(3,1) * t926;
t753 = t785 * mrSges(3,2) - mrSges(2,1);
t816 = (t753 * t795 + t781 * t786) * pkin(1);
t826 = -t892 + t923;
t708 = -(t756 + 0.2e1 * t922) * t794 + t816 + t826;
t750 = -Ifges(3,5) * t785 - Ifges(3,6) * t794;
t813 = (t729 * Ifges(3,5) + 0.2e1 * t889 * t916) * t794 - Ifges(3,6) * t910 - (0.4e1 * t771 - 0.2e1) * Ifges(3,4) * t916;
t848 = t772 * t910;
t866 = t708 * t663 + (t922 * t935 + t826) * t654 - t750 * t848 + (t886 + (mrSges(3,1) * t794 - t753) * t786) * t690 * pkin(1) + t813;
t679 = t699 / 0.2e1 + t698 / 0.2e1 + t718 / 0.2e1 + t694;
t844 = t733 * t882;
t655 = (-t807 * t913 + (pkin(1) * t694 + (0.2e1 * t679 * t878 - t844) * pkin(2)) * pkin(1)) * t775 * t694 * t839 + (-pkin(2) * t913 + (-t862 * t878 - t844) * pkin(1)) * t688 * t840 + ((pkin(1) * t862 * t882 + pkin(2) * t733) * t797 + pkin(1) * t877) * t776 * t733 * t894;
t664 = (-t938 + (-t797 * t682 - t730 * t775) * pkin(2)) * t894;
t757 = mrSges(3,1) * t925;
t754 = t788 * mrSges(3,2) - mrSges(2,1);
t815 = (t754 * t798 + t781 * t789) * pkin(1);
t825 = -t891 + t923;
t709 = -(t757 + 0.2e1 * t921) * t797 + t815 + t825;
t751 = -Ifges(3,5) * t788 - Ifges(3,6) * t797;
t812 = (t730 * Ifges(3,5) + 0.2e1 * t888 * t914) * t797 - Ifges(3,6) * t909 - (0.4e1 * t774 - 0.2e1) * Ifges(3,4) * t914;
t847 = t775 * t909;
t865 = t709 * t664 + (t921 * t934 + t825) * t655 - t751 * t847 + (t885 + (mrSges(3,1) * t797 - t754) * t789) * t691 * pkin(1) + t812;
t680 = t701 / 0.2e1 + t700 / 0.2e1 + t719 / 0.2e1 + t695;
t843 = t734 * t881;
t656 = (-t807 * t911 + (pkin(1) * t695 + (0.2e1 * t680 * t876 - t843) * pkin(2)) * pkin(1)) * t778 * t695 * t837 + (-pkin(2) * t911 + (-t861 * t876 - t843) * pkin(1)) * t689 * t838 + ((pkin(1) * t861 * t881 + pkin(2) * t734) * t800 + pkin(1) * t875) * t779 * t734 * t893;
t665 = (-t937 + (-t800 * t683 - t731 * t778) * pkin(2)) * t893;
t758 = mrSges(3,1) * t924;
t755 = t791 * mrSges(3,2) - mrSges(2,1);
t814 = (t755 * t801 + t781 * t792) * pkin(1);
t824 = -t890 + t923;
t710 = -(t758 + 0.2e1 * t920) * t800 + t814 + t824;
t752 = -Ifges(3,5) * t791 - Ifges(3,6) * t800;
t811 = (t731 * Ifges(3,5) + 0.2e1 * t887 * t912) * t800 - Ifges(3,6) * t908 - (0.4e1 * t777 - 0.2e1) * Ifges(3,4) * t912;
t846 = t778 * t908;
t864 = t710 * t665 + (t920 * t933 + t824) * t656 - t752 * t846 + (t884 + (mrSges(3,1) * t800 - t755) * t792) * t692 * pkin(1) + t811;
t851 = t863 * t879;
t850 = t862 * t877;
t849 = t861 * t875;
t735 = -(-mrSges(3,2) * t932 + Ifges(3,6)) * t794 + t785 * (mrSges(3,1) * t932 - Ifges(3,5));
t817 = -(m(2) + m(3)) * pkin(1) ^ 2 - Ifges(1,3) + t923;
t835 = (t813 + ((-mrSges(3,1) * t919 + mrSges(3,2) * t851) * t794 + (mrSges(3,1) * t851 + mrSges(3,2) * t919) * t785 + (-mrSges(2,1) * t786 - t886) * t687 * t678) * t940 + (-t892 + (t756 + t922) * t935 + 0.2e1 * t816 + t817) * t663 + t708 * t654 - t735 * t848) * t768;
t736 = -(-mrSges(3,2) * t930 + Ifges(3,6)) * t797 + t788 * (mrSges(3,1) * t930 - Ifges(3,5));
t834 = (t812 + ((-mrSges(3,1) * t918 + mrSges(3,2) * t850) * t797 + (mrSges(3,1) * t850 + mrSges(3,2) * t918) * t788 + (-mrSges(2,1) * t789 - t885) * t688 * t679) * t940 + (-t891 + (t757 + t921) * t934 + 0.2e1 * t815 + t817) * t664 + t709 * t655 - t736 * t847) * t769;
t737 = -(-mrSges(3,2) * t928 + Ifges(3,6)) * t800 + t791 * (mrSges(3,1) * t928 - Ifges(3,5));
t833 = (t811 + ((-mrSges(3,1) * t917 + mrSges(3,2) * t849) * t800 + (mrSges(3,1) * t849 + mrSges(3,2) * t917) * t791 + (-mrSges(2,1) * t792 - t884) * t689 * t680) * t940 + (-t890 + (t758 + t920) * t933 + 0.2e1 * t814 + t817) * t665 + t710 * t656 - t737 * t846) * t770;
t832 = (Ifges(3,3) * t848 + t750 * t654 + t735 * t663 + (mrSges(3,2) * t939 + t681 * t889) * t794 + mrSges(3,1) * t690 * t857 + (-0.2e1 * t771 + 0.1e1) * t681 * Ifges(3,4)) * t772;
t831 = (Ifges(3,3) * t847 + t751 * t655 + t736 * t664 + (mrSges(3,2) * t938 + t682 * t888) * t797 + mrSges(3,1) * t691 * t856 + (-0.2e1 * t774 + 0.1e1) * t682 * Ifges(3,4)) * t775;
t830 = (Ifges(3,3) * t846 + t752 * t656 + t737 * t665 + (mrSges(3,2) * t937 + t683 * t887) * t800 + mrSges(3,1) * t692 * t855 + (-0.2e1 * t777 + 0.1e1) * t683 * Ifges(3,4)) * t778;
t823 = t772 * t835;
t822 = t775 * t834;
t821 = t778 * t833;
t820 = t866 * t773 * t768;
t819 = t865 * t776 * t769;
t818 = t864 * t779 * t770;
t1 = [(t721 * t823 + t723 * t822 + t725 * t821) * t810 + (t764 * t830 + t763 * t831 + t762 * t832 + (t703 * t820 + t705 * t819 + t707 * t818) * t810) * t808; (t720 * t823 + t722 * t822 + t724 * t821) * t810 + (t767 * t830 + t766 * t831 + t765 * t832 + (t702 * t820 + t704 * t819 + t706 * t818) * t810) * t808; (t761 * t833 + t760 * t834 + t759 * t835 + (t864 * t905 + t865 * t906 + t866 * t907) * t808) * t810;];
taucX  = t1;
