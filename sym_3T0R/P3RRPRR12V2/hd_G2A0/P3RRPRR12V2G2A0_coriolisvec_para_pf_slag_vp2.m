% Calculate vector of centrifugal and coriolis load on the joints for
% P3RRPRR12V2G2A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2020-08-06 19:23
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR12V2G2A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G2A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V2G2A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G2A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G2A0_coriolisvec_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V2G2A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR12V2G2A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR12V2G2A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G2A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G2A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:20:58
% EndTime: 2020-08-06 19:21:06
% DurationCPUTime: 8.64s
% Computational Cost: add. (79992->459), mult. (94779->718), div. (11403->6), fcn. (72603->18), ass. (0->295)
t812 = sin(qJ(2,3));
t813 = sin(qJ(1,3));
t819 = cos(qJ(1,3));
t829 = pkin(5) - pkin(6);
t856 = pkin(1) * t813 - t829 * t819;
t923 = t813 * qJ(3,3);
t726 = t856 * t812 + t923;
t863 = t812 * t923;
t729 = t856 + 0.2e1 * t863;
t830 = pkin(2) + pkin(3);
t928 = (qJ(3,3) + t830) * (-qJ(3,3) + t830);
t871 = t812 * t928;
t976 = pkin(1) * qJ(3,3);
t736 = -t871 + t976;
t972 = pkin(1) * t812;
t756 = qJ(3,3) + t972;
t809 = legFrame(3,2);
t772 = sin(t809);
t775 = cos(t809);
t818 = cos(qJ(2,3));
t800 = t818 ^ 2;
t870 = t813 * t928;
t985 = -0.2e1 * t830;
t883 = qJ(3,3) * t985;
t936 = t775 * t830;
t942 = t772 * t830;
t943 = t772 * qJ(3,3);
t693 = (-t772 * t870 + t775 * t883) * t800 + (-t729 * t942 - t775 * t736) * t818 - t726 * t943 + t756 * t936;
t937 = t775 * qJ(3,3);
t694 = (t772 * t883 + t775 * t870) * t800 + (t729 * t936 - t772 * t736) * t818 + t726 * t937 + t756 * t942;
t759 = t813 * t829;
t769 = t812 * qJ(3,3);
t872 = t800 * t928;
t909 = t830 * t818;
t717 = t819 * t872 + ((0.2e1 * t769 + pkin(1)) * t819 + t759) * t909 + qJ(3,3) * (t756 * t819 + t812 * t759);
t826 = xDP(3);
t827 = xDP(2);
t828 = xDP(1);
t753 = t769 + pkin(1);
t992 = t753 + t909;
t739 = 0.1e1 / t992;
t833 = 0.1e1 / qJ(3,3);
t946 = t739 * t833;
t672 = (t693 * t827 + t694 * t828 + t717 * t826) * t946;
t730 = t856 + t863;
t922 = t813 * t830;
t924 = t812 * t830;
t702 = (-t772 * t922 - t937) * t800 + (-t730 * t772 + t775 * t924) * t818 + t775 * t756;
t703 = (t775 * t922 - t943) * t800 + (t730 * t775 + t772 * t924) * t818 + t772 * t756;
t949 = (t992 * t819 + t759) * t818;
t684 = (t702 * t827 + t703 * t828 + t826 * t949) * t946;
t912 = t830 * t684;
t667 = t672 - t912;
t814 = sin(qJ(2,2));
t815 = sin(qJ(1,2));
t821 = cos(qJ(1,2));
t855 = pkin(1) * t815 - t829 * t821;
t919 = t815 * qJ(3,2);
t727 = t855 * t814 + t919;
t862 = t814 * t919;
t731 = t855 + 0.2e1 * t862;
t927 = (qJ(3,2) + t830) * (-qJ(3,2) + t830);
t868 = t814 * t927;
t977 = pkin(1) * qJ(3,2);
t737 = -t868 + t977;
t971 = pkin(1) * t814;
t757 = qJ(3,2) + t971;
t810 = legFrame(2,2);
t773 = sin(t810);
t776 = cos(t810);
t820 = cos(qJ(2,2));
t801 = t820 ^ 2;
t867 = t815 * t927;
t884 = qJ(3,2) * t985;
t934 = t776 * t830;
t940 = t773 * t830;
t941 = t773 * qJ(3,2);
t695 = (-t773 * t867 + t776 * t884) * t801 + (-t731 * t940 - t776 * t737) * t820 - t727 * t941 + t757 * t934;
t935 = t776 * qJ(3,2);
t696 = (t773 * t884 + t776 * t867) * t801 + (t731 * t934 - t773 * t737) * t820 + t727 * t935 + t757 * t940;
t760 = t815 * t829;
t770 = t814 * qJ(3,2);
t869 = t801 * t927;
t908 = t830 * t820;
t718 = t821 * t869 + ((0.2e1 * t770 + pkin(1)) * t821 + t760) * t908 + qJ(3,2) * (t757 * t821 + t814 * t760);
t754 = t770 + pkin(1);
t991 = t754 + t908;
t740 = 0.1e1 / t991;
t835 = 0.1e1 / qJ(3,2);
t945 = t740 * t835;
t673 = (t695 * t827 + t696 * t828 + t718 * t826) * t945;
t732 = t855 + t862;
t918 = t815 * t830;
t920 = t814 * t830;
t704 = (-t773 * t918 - t935) * t801 + (-t732 * t773 + t776 * t920) * t820 + t776 * t757;
t705 = (t776 * t918 - t941) * t801 + (t732 * t776 + t773 * t920) * t820 + t773 * t757;
t948 = (t991 * t821 + t760) * t820;
t685 = (t704 * t827 + t705 * t828 + t826 * t948) * t945;
t911 = t830 * t685;
t668 = t673 - t911;
t816 = sin(qJ(2,1));
t817 = sin(qJ(1,1));
t823 = cos(qJ(1,1));
t854 = pkin(1) * t817 - t829 * t823;
t915 = t817 * qJ(3,1);
t728 = t854 * t816 + t915;
t861 = t816 * t915;
t733 = t854 + 0.2e1 * t861;
t926 = (qJ(3,1) + t830) * (-qJ(3,1) + t830);
t865 = t816 * t926;
t978 = pkin(1) * qJ(3,1);
t738 = -t865 + t978;
t970 = pkin(1) * t816;
t758 = qJ(3,1) + t970;
t811 = legFrame(1,2);
t774 = sin(t811);
t777 = cos(t811);
t822 = cos(qJ(2,1));
t802 = t822 ^ 2;
t864 = t817 * t926;
t885 = qJ(3,1) * t985;
t932 = t777 * t830;
t938 = t774 * t830;
t939 = t774 * qJ(3,1);
t697 = (-t774 * t864 + t777 * t885) * t802 + (-t733 * t938 - t777 * t738) * t822 - t728 * t939 + t758 * t932;
t933 = t777 * qJ(3,1);
t698 = (t774 * t885 + t777 * t864) * t802 + (t733 * t932 - t774 * t738) * t822 + t728 * t933 + t758 * t938;
t761 = t817 * t829;
t771 = t816 * qJ(3,1);
t866 = t802 * t926;
t907 = t830 * t822;
t719 = t823 * t866 + ((0.2e1 * t771 + pkin(1)) * t823 + t761) * t907 + qJ(3,1) * (t758 * t823 + t816 * t761);
t755 = t771 + pkin(1);
t990 = t755 + t907;
t741 = 0.1e1 / t990;
t837 = 0.1e1 / qJ(3,1);
t944 = t741 * t837;
t674 = (t697 * t827 + t698 * t828 + t719 * t826) * t944;
t734 = t854 + t861;
t914 = t817 * t830;
t916 = t816 * t830;
t706 = (-t774 * t914 - t933) * t802 + (-t734 * t774 + t777 * t916) * t822 + t777 * t758;
t707 = (t777 * t914 - t939) * t802 + (t734 * t777 + t774 * t916) * t822 + t774 * t758;
t947 = (t990 * t823 + t761) * t822;
t686 = (t706 * t827 + t707 * t828 + t826 * t947) * t944;
t910 = t830 * t686;
t666 = t674 - t910;
t879 = pkin(2) * mrSges(3,3) + Ifges(2,4) - Ifges(3,5);
t996 = m(3) * pkin(2) + mrSges(3,1);
t1003 = qJ(3,1) * t996 + t879;
t1002 = qJ(3,2) * t996 + t879;
t1001 = qJ(3,3) * t996 + t879;
t887 = 0.2e1 * pkin(1);
t1000 = 0.2e1 * t830;
t836 = qJ(3,1) ^ 2;
t840 = pkin(2) ^ 2;
t999 = (t836 - t840) * m(3);
t834 = qJ(3,2) ^ 2;
t998 = (t834 - t840) * m(3);
t832 = qJ(3,3) ^ 2;
t997 = (t832 - t840) * m(3);
t768 = qJ(3,1) * m(3) + mrSges(3,3);
t766 = m(3) * qJ(3,3) + mrSges(3,3);
t767 = m(3) * qJ(3,2) + mrSges(3,3);
t906 = t832 * t684;
t657 = t830 * t667 - t906;
t664 = t829 * t667;
t714 = (-t813 * t826 + (-t772 * t827 + t775 * t828) * t819) * t739;
t888 = -pkin(1) ^ 2 - pkin(5) ^ 2;
t989 = -0.2e1 * pkin(5);
t853 = -t888 + (t989 + pkin(6)) * pkin(6);
t749 = t832 + t853;
t804 = t830 ^ 2;
t880 = t714 * t976;
t894 = 0.2e1 * t880 + t664;
t913 = t829 * t830;
t925 = t812 * t829;
t956 = t714 * t829;
t957 = t714 * t812;
t958 = t714 * t800;
t875 = t714 * t925;
t961 = (t875 - t912) * t818;
t981 = t800 - 0.1e1;
t988 = -0.2e1 * t800;
t645 = ((-(t804 - 0.3e1 * t832) * t909 * t958 + (t829 * (-t1000 * t684 + t672) * qJ(3,3) + (-0.3e1 * (-t832 / 0.3e1 + t804) * t769 + (t832 - t804) * t887) * t714) * t800 + (-t906 * t925 + ((-0.4e1 * t880 - t664) * t812 - t714 * (0.3e1 * t832 + t853)) * t830) * t818 - (t749 * t957 + t894) * qJ(3,3)) * t714 + ((t830 * t657 + t871 * t956) * t818 + t657 * pkin(1) + (t657 * t812 + (t988 + 0.1e1) * t714 * t913) * qJ(3,3)) * t684 + ((pkin(1) * t684 - t961) * t830 + (t812 * t912 + t981 * t956) * qJ(3,3)) * t672) * t946;
t860 = t818 * qJ(3,3) * t684;
t648 = (-t818 * (t829 * t860 + t894 * t812 + (t1000 * t753 * t818 + t749 + t872) * t714) * t714 + (-qJ(3,3) * t800 * t956 + (t667 + t875) * t909 + t667 * t753) * t684 + (t684 * t753 - t961) * t672) * t946;
t654 = (0.2e1 * t667 * t812 + 0.2e1 * t860 + t956) * t714 * t739;
t678 = t1001 * t684;
t681 = t684 ^ 2;
t765 = mrSges(2,1) + t996;
t742 = pkin(2) * mrSges(3,2) + t765 * pkin(5) - Ifges(3,4) - Ifges(2,5);
t762 = -mrSges(2,2) + t766;
t968 = Ifges(3,6) - Ifges(2,6);
t850 = -qJ(3,3) * mrSges(3,2) + t968;
t720 = -(t762 * pkin(5) - t850) * t818 + t812 * t742;
t983 = m(3) * pkin(5);
t735 = (-mrSges(2,1) / 0.2e1 - mrSges(3,1) / 0.2e1) * pkin(5) + Ifges(3,4) / 0.2e1 + Ifges(2,5) / 0.2e1 + (-t983 / 0.2e1 - mrSges(3,2) / 0.2e1) * pkin(2);
t752 = t765 * t887;
t778 = mrSges(3,2) + t983;
t965 = mrSges(3,3) * qJ(3,3);
t792 = 0.2e1 * t965;
t793 = -0.2e1 * t965;
t982 = mrSges(3,1) * pkin(2);
t798 = 0.2e1 * t982;
t799 = -0.2e1 * t982;
t825 = mrSges(2,2) * pkin(5);
t969 = Ifges(2,1) + Ifges(3,1);
t845 = t888 * m(2) + (mrSges(3,2) + mrSges(2,3)) * t989 - Ifges(1,3) - t969;
t858 = Ifges(2,2) + Ifges(3,3) - t969;
t848 = -t858 + t997;
t931 = t778 * t812;
t964 = t672 * t766;
t975 = pkin(1) * t762;
t903 = (-(t793 + t798 - t848) * t800 - 0.2e1 * t762 * t972 - (t832 - t888) * m(3) + t793 + t845 + (-0.2e1 * t1001 * t812 - t752) * t818) * t654 + t720 * t648 - t645 * t931 + 0.4e1 * (t678 - t964 / 0.2e1) * t958 + 0.2e1 * (((t792 + t799 + t848) * t684 + t672 * t996) * t957 + (t672 * t778 + t735 * t684 + t714 * t975) * t684) * t818 + ((-t766 * pkin(5) + t825 + t850) * t681 + (m(3) * t672 - t684 * t765) * t714 * t887) * t812 - 0.2e1 * t714 * (t678 - t964);
t995 = t819 * t903;
t905 = t834 * t685;
t658 = t830 * t668 - t905;
t665 = t829 * t668;
t715 = (-t815 * t826 + (-t773 * t827 + t776 * t828) * t821) * t740;
t750 = t834 + t853;
t881 = t715 * t977;
t893 = 0.2e1 * t881 + t665;
t921 = t814 * t829;
t953 = t715 * t829;
t954 = t715 * t814;
t955 = t715 * t801;
t874 = t715 * t921;
t960 = (t874 - t911) * t820;
t980 = t801 - 0.1e1;
t987 = -0.2e1 * t801;
t646 = ((-(t804 - 0.3e1 * t834) * t908 * t955 + (t829 * (-t1000 * t685 + t673) * qJ(3,2) + (-0.3e1 * (-t834 / 0.3e1 + t804) * t770 + (t834 - t804) * t887) * t715) * t801 + (-t905 * t921 + ((-0.4e1 * t881 - t665) * t814 - t715 * (0.3e1 * t834 + t853)) * t830) * t820 - (t750 * t954 + t893) * qJ(3,2)) * t715 + ((t830 * t658 + t868 * t953) * t820 + t658 * pkin(1) + (t658 * t814 + (t987 + 0.1e1) * t715 * t913) * qJ(3,2)) * t685 + ((pkin(1) * t685 - t960) * t830 + (t814 * t911 + t980 * t953) * qJ(3,2)) * t673) * t945;
t859 = t820 * qJ(3,2) * t685;
t649 = (-t820 * (t829 * t859 + t893 * t814 + (t1000 * t754 * t820 + t750 + t869) * t715) * t715 + (-qJ(3,2) * t801 * t953 + (t668 + t874) * t908 + t668 * t754) * t685 + (t685 * t754 - t960) * t673) * t945;
t655 = (0.2e1 * t668 * t814 + 0.2e1 * t859 + t953) * t715 * t740;
t679 = t1002 * t685;
t682 = t685 ^ 2;
t763 = -mrSges(2,2) + t767;
t851 = -qJ(3,2) * mrSges(3,2) + t968;
t721 = -(t763 * pkin(5) - t851) * t820 + t814 * t742;
t966 = mrSges(3,3) * qJ(3,2);
t794 = 0.2e1 * t966;
t795 = -0.2e1 * t966;
t847 = -t858 + t998;
t930 = t778 * t814;
t963 = t673 * t767;
t974 = pkin(1) * t763;
t902 = (-(t795 + t798 - t847) * t801 - 0.2e1 * t763 * t971 - (t834 - t888) * m(3) + t795 + t845 + (-0.2e1 * t1002 * t814 - t752) * t820) * t655 + t721 * t649 - t646 * t930 + 0.4e1 * (t679 - t963 / 0.2e1) * t955 + 0.2e1 * (((t794 + t799 + t847) * t685 + t673 * t996) * t954 + (t673 * t778 + t735 * t685 + t715 * t974) * t685) * t820 + ((-t767 * pkin(5) + t825 + t851) * t682 + (m(3) * t673 - t685 * t765) * t715 * t887) * t814 - 0.2e1 * t715 * (t679 - t963);
t994 = t821 * t902;
t904 = t836 * t686;
t659 = t830 * t666 - t904;
t663 = t829 * t666;
t716 = (-t817 * t826 + (-t774 * t827 + t777 * t828) * t823) * t741;
t751 = t836 + t853;
t882 = t716 * t978;
t892 = 0.2e1 * t882 + t663;
t917 = t816 * t829;
t950 = t716 * t829;
t951 = t716 * t816;
t952 = t716 * t802;
t873 = t716 * t917;
t959 = (t873 - t910) * t822;
t979 = t802 - 0.1e1;
t986 = -0.2e1 * t802;
t647 = ((-(t804 - 0.3e1 * t836) * t907 * t952 + (t829 * (-t1000 * t686 + t674) * qJ(3,1) + (-0.3e1 * (-t836 / 0.3e1 + t804) * t771 + (t836 - t804) * t887) * t716) * t802 + (-t904 * t917 + ((-0.4e1 * t882 - t663) * t816 - t716 * (0.3e1 * t836 + t853)) * t830) * t822 - (t751 * t951 + t892) * qJ(3,1)) * t716 + ((t830 * t659 + t865 * t950) * t822 + t659 * pkin(1) + (t659 * t816 + (t986 + 0.1e1) * t716 * t913) * qJ(3,1)) * t686 + ((pkin(1) * t686 - t959) * t830 + (t816 * t910 + t979 * t950) * qJ(3,1)) * t674) * t944;
t876 = t686 * t822 * qJ(3,1);
t650 = (-t822 * (t829 * t876 + t892 * t816 + (t1000 * t755 * t822 + t751 + t866) * t716) * t716 + (-qJ(3,1) * t802 * t950 + (t666 + t873) * t907 + t666 * t755) * t686 + (t686 * t755 - t959) * t674) * t944;
t656 = (0.2e1 * t666 * t816 + 0.2e1 * t876 + t950) * t716 * t741;
t680 = t1003 * t686;
t683 = t686 ^ 2;
t764 = -mrSges(2,2) + t768;
t852 = -qJ(3,1) * mrSges(3,2) + t968;
t722 = -(t764 * pkin(5) - t852) * t822 + t816 * t742;
t967 = mrSges(3,3) * qJ(3,1);
t796 = 0.2e1 * t967;
t797 = -0.2e1 * t967;
t846 = -t858 + t999;
t929 = t778 * t816;
t962 = t674 * t768;
t973 = pkin(1) * t764;
t901 = (-(t797 + t798 - t846) * t802 - 0.2e1 * t764 * t970 - (t836 - t888) * m(3) + t797 + t845 + (-0.2e1 * t1003 * t816 - t752) * t822) * t656 + t722 * t650 - t647 * t929 + 0.4e1 * (t680 - t962 / 0.2e1) * t952 + 0.2e1 * (((t796 + t799 + t846) * t686 + t674 * t996) * t951 + (t674 * t778 + t735 * t686 + t716 * t973) * t686) * t822 + ((-t768 * pkin(5) + t825 + t852) * t683 + (m(3) * t674 - t686 * t765) * t716 * t887) * t816 - 0.2e1 * t716 * (t680 - t962);
t993 = t823 * t901;
t984 = m(3) * pkin(1);
t711 = t714 ^ 2;
t849 = t799 - t858;
t878 = t799 - Ifges(3,2) - Ifges(2,3);
t900 = t720 * t654 + (-m(3) * (t832 + t840) + t793 + t878) * t648 + t996 * t645 + 0.2e1 * t684 * t964 + (t1001 * t988 - ((t792 + t849 + t997) * t812 + t975) * t818 + t765 * t972 + t1001) * t711;
t712 = t715 ^ 2;
t899 = t721 * t655 + (-m(3) * (t834 + t840) + t795 + t878) * t649 + t996 * t646 + 0.2e1 * t685 * t963 + (t1002 * t987 - ((t794 + t849 + t998) * t814 + t974) * t820 + t765 * t971 + t1002) * t712;
t713 = t716 ^ 2;
t898 = t722 * t656 + (-m(3) * (t836 + t840) + t797 + t878) * t650 + t996 * t647 + 0.2e1 * t686 * t962 + (t1003 * t986 - ((t796 + t849 + t999) * t816 + t973) * t822 + t765 * t970 + t1003) * t713;
t897 = -m(3) * t645 + t996 * t648 - t654 * t931 - t681 * t766 + (t981 * t766 + (-t818 * t996 - t984) * t812) * t711;
t896 = -m(3) * t646 + t996 * t649 - t655 * t930 - t682 * t767 + (t980 * t767 + (-t820 * t996 - t984) * t814) * t712;
t895 = -m(3) * t647 + t996 * t650 - t656 * t929 - t683 * t768 + (t979 * t768 + (-t822 * t996 - t984) * t816) * t713;
t1 = [(t777 * t993 + (t895 * t698 + t898 * t707) * t837) * t741 + (t776 * t994 + (t896 * t696 + t899 * t705) * t835) * t740 + (t775 * t995 + (t897 * t694 + t900 * t703) * t833) * t739; (-t774 * t993 + (t895 * t697 + t898 * t706) * t837) * t741 + (-t773 * t994 + (t896 * t695 + t899 * t704) * t835) * t740 + (-t772 * t995 + (t897 * t693 + t900 * t702) * t833) * t739; (-t901 * t817 + (t719 * t895 + t898 * t947) * t837) * t741 + (-t902 * t815 + (t718 * t896 + t899 * t948) * t835) * t740 + (-t903 * t813 + (t717 * t897 + t900 * t949) * t833) * t739;];
taucX  = t1;
