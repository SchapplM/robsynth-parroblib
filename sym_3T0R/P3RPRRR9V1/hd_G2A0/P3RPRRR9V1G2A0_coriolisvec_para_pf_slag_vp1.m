% Calculate vector of centrifugal and coriolis load on the joints for
% P3RPRRR9V1G2A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-08-06 18:53
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRRR9V1G2A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G2A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR9V1G2A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G2A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G2A0_coriolisvec_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR9V1G2A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR9V1G2A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRRR9V1G2A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G2A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G2A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:51:25
% EndTime: 2020-08-06 18:51:28
% DurationCPUTime: 3.49s
% Computational Cost: add. (18744->360), mult. (24129->591), div. (3921->10), fcn. (20250->40), ass. (0->285)
t1037 = 2 * pkin(1);
t865 = cos(pkin(7));
t1021 = 0.2e1 * t865 ^ 2;
t875 = cos(qJ(3,3));
t860 = t875 ^ 2;
t870 = sin(qJ(1,3));
t998 = pkin(5) + qJ(2,3);
t850 = pkin(6) + t998;
t876 = cos(qJ(1,3));
t932 = pkin(1) * t870 - t876 * t850;
t864 = sin(pkin(7));
t869 = sin(qJ(3,3));
t972 = t864 * t869;
t744 = t932 * t972 + (t860 - 0.1e1) * t870 * pkin(3);
t1003 = t875 * pkin(2);
t1006 = t860 * pkin(3);
t759 = pkin(1) * t869 + (-pkin(3) + t1003 + 0.2e1 * t1006) * t864;
t897 = -pkin(3) / 0.2e1;
t779 = t1006 + t1003 / 0.2e1 + t897;
t866 = legFrame(3,2);
t835 = sin(t866);
t838 = cos(t866);
t1023 = 0.2e1 * pkin(3);
t944 = t870 * t972;
t914 = pkin(2) * t944 + (t944 * t1023 - t932) * t875;
t831 = pkin(1) * t864;
t969 = t875 * (-t869 * pkin(3) + t831);
t978 = t835 * t870;
t898 = pkin(2) / 0.2e1;
t985 = (t875 * pkin(3) + t898) * t869;
t723 = (-t779 * t978 + t838 * t985) * t1021 + (t838 * t759 + t914 * t835) * t865 + t744 * t835 + t838 * t969;
t887 = xDP(2);
t772 = 0.1e1 / (t865 * t875 - t972);
t832 = 0.1e1 / t850;
t988 = t772 * t832;
t717 = t723 * t887 * t988;
t975 = t838 * t870;
t724 = (t779 * t975 + t835 * t985) * t1021 + (t835 * t759 - t914 * t838) * t865 - t744 * t838 + t835 * t969;
t888 = xDP(1);
t718 = t724 * t888 * t988;
t857 = pkin(7) + qJ(3,3);
t827 = cos(t857);
t1012 = pkin(3) * t827;
t830 = t865 * pkin(2);
t809 = t830 + pkin(1);
t776 = t809 + t1012;
t756 = t776 * t876 + t870 * t850;
t886 = xDP(3);
t741 = t756 * t832 * t886;
t714 = t718 + t717 + t741;
t1036 = 0.2e1 * t714;
t877 = cos(qJ(3,2));
t861 = t877 ^ 2;
t872 = sin(qJ(1,2));
t999 = pkin(5) + qJ(2,2);
t851 = pkin(6) + t999;
t878 = cos(qJ(1,2));
t931 = pkin(1) * t872 - t878 * t851;
t871 = sin(qJ(3,2));
t971 = t864 * t871;
t745 = t931 * t971 + (t861 - 0.1e1) * t872 * pkin(3);
t1002 = t877 * pkin(2);
t1005 = t861 * pkin(3);
t760 = pkin(1) * t871 + (-pkin(3) + t1002 + 0.2e1 * t1005) * t864;
t780 = t1005 + t1002 / 0.2e1 + t897;
t867 = legFrame(2,2);
t836 = sin(t867);
t839 = cos(t867);
t943 = t872 * t971;
t913 = pkin(2) * t943 + (t943 * t1023 - t931) * t877;
t968 = t877 * (-t871 * pkin(3) + t831);
t977 = t836 * t872;
t984 = (t877 * pkin(3) + t898) * t871;
t725 = (-t780 * t977 + t839 * t984) * t1021 + (t839 * t760 + t913 * t836) * t865 + t745 * t836 + t839 * t968;
t773 = 0.1e1 / (t865 * t877 - t971);
t833 = 0.1e1 / t851;
t987 = t773 * t833;
t719 = t725 * t887 * t987;
t974 = t839 * t872;
t726 = (t780 * t974 + t836 * t984) * t1021 + (t836 * t760 - t913 * t839) * t865 - t745 * t839 + t836 * t968;
t720 = t726 * t888 * t987;
t858 = pkin(7) + qJ(3,2);
t828 = cos(t858);
t1011 = pkin(3) * t828;
t777 = t809 + t1011;
t757 = t777 * t878 + t872 * t851;
t742 = t757 * t833 * t886;
t715 = t720 + t719 + t742;
t1035 = 0.2e1 * t715;
t879 = cos(qJ(3,1));
t862 = t879 ^ 2;
t874 = sin(qJ(1,1));
t1000 = pkin(5) + qJ(2,1);
t852 = pkin(6) + t1000;
t880 = cos(qJ(1,1));
t930 = pkin(1) * t874 - t880 * t852;
t873 = sin(qJ(3,1));
t970 = t864 * t873;
t746 = t930 * t970 + (t862 - 0.1e1) * t874 * pkin(3);
t1001 = t879 * pkin(2);
t1004 = t862 * pkin(3);
t761 = pkin(1) * t873 + (-pkin(3) + t1001 + 0.2e1 * t1004) * t864;
t781 = t1004 + t1001 / 0.2e1 + t897;
t868 = legFrame(1,2);
t837 = sin(t868);
t840 = cos(t868);
t942 = t874 * t970;
t912 = pkin(2) * t942 + (t942 * t1023 - t930) * t879;
t967 = t879 * (-t873 * pkin(3) + t831);
t976 = t837 * t874;
t983 = (t879 * pkin(3) + t898) * t873;
t727 = (-t781 * t976 + t840 * t983) * t1021 + (t840 * t761 + t912 * t837) * t865 + t746 * t837 + t840 * t967;
t774 = 0.1e1 / (t865 * t879 - t970);
t834 = 0.1e1 / t852;
t986 = t774 * t834;
t721 = t727 * t887 * t986;
t973 = t840 * t874;
t728 = (t781 * t973 + t837 * t983) * t1021 + (t837 * t761 - t912 * t840) * t865 - t746 * t840 + t837 * t967;
t722 = t728 * t888 * t986;
t859 = pkin(7) + qJ(3,1);
t829 = cos(t859);
t1010 = pkin(3) * t829;
t778 = t809 + t1010;
t758 = t778 * t880 + t874 * t852;
t743 = t758 * t834 * t886;
t716 = t722 + t721 + t743;
t1034 = 0.2e1 * t716;
t821 = sin(t857);
t762 = t838 * t821 - t827 * t978;
t763 = t835 * t821 + t827 * t975;
t815 = 0.1e1 / t827;
t735 = (t876 * t886 + (t762 * t887 + t763 * t888) * t815) * t832;
t1033 = t735 / 0.2e1;
t822 = sin(t858);
t764 = t839 * t822 - t828 * t977;
t765 = t836 * t822 + t828 * t974;
t816 = 0.1e1 / t828;
t736 = (t878 * t886 + (t764 * t887 + t765 * t888) * t816) * t833;
t1032 = t736 / 0.2e1;
t823 = sin(t859);
t766 = t840 * t823 - t829 * t976;
t767 = t837 * t823 + t829 * t973;
t817 = 0.1e1 / t829;
t737 = (t880 * t886 + (t766 * t887 + t767 * t888) * t817) * t834;
t1031 = t737 / 0.2e1;
t1030 = -0.2e1 * t830;
t1029 = t1037 / 0.2e1;
t957 = 2 * m(3);
t1028 = (rSges(3,1) * t823 + rSges(3,2) * t829) * t957;
t1027 = (rSges(3,1) * t822 + rSges(3,2) * t828) * t957;
t1026 = (rSges(3,1) * t821 + rSges(3,2) * t827) * t957;
t1025 = -2 * pkin(1);
t1022 = 4 * rSges(2,3);
t1020 = -4 * pkin(5) - 4 * pkin(6);
t1019 = -m(2) / 0.2e1;
t1018 = pkin(2) * m(3);
t1017 = -rSges(3,1) / 0.2e1;
t1016 = -rSges(3,2) / 0.2e1;
t1015 = m(3) * rSges(3,2);
t1014 = rSges(2,2) * m(2);
t902 = rSges(3,2) ^ 2;
t904 = rSges(3,1) ^ 2;
t790 = (-t902 + t904) * m(3) - Icges(3,1) + Icges(3,2);
t1013 = -t790 / 0.2e1;
t846 = (rSges(3,3) + t998);
t1009 = t846 * m(3);
t847 = (rSges(3,3) + t999);
t1008 = t847 * m(3);
t848 = (rSges(3,3) + t1000);
t1007 = t848 * m(3);
t810 = 0.2e1 * t857;
t799 = sin(t810);
t997 = t735 * t799;
t811 = 0.2e1 * t858;
t800 = sin(t811);
t996 = t736 * t800;
t812 = 0.2e1 * t859;
t801 = sin(t812);
t995 = t737 * t801;
t907 = 0.1e1 / pkin(3);
t753 = (t835 * t888 + t838 * t887) * t907 * t815;
t994 = t753 ^ 2 * t815;
t754 = (t836 * t888 + t839 * t887) * t907 * t816;
t993 = t754 ^ 2 * t816;
t755 = (t837 * t888 + t840 * t887) * t907 * t817;
t992 = t755 ^ 2 * t817;
t991 = t753 * t821;
t990 = t754 * t822;
t989 = t755 * t823;
t982 = (m(2) * rSges(2,1) + t1018) * t865;
t802 = cos(t810);
t814 = rSges(3,1) * t1015 - Icges(3,4);
t981 = t814 * t802;
t803 = cos(t811);
t980 = t814 * t803;
t804 = cos(t812);
t979 = t814 * t804;
t732 = pkin(1) * t735;
t708 = t732 - t718 / 0.2e1 - t717 / 0.2e1 - t741 / 0.2e1;
t896 = 0.2e1 * pkin(7);
t854 = t896 + qJ(3,3);
t824 = cos(t854);
t899 = (qJ(2,3) ^ 2);
t906 = pkin(3) ^ 2;
t849 = cos(t896);
t908 = pkin(2) ^ 2;
t909 = pkin(1) ^ 2;
t911 = -(2 * pkin(6) ^ 2) - t908 * t849 - t906 - t908 - (2 * t909) + ((-4 * pkin(6) - 2 * pkin(5)) * pkin(5));
t693 = ((t708 * t1030 + ((qJ(2,3) * t1020) - t906 * t802 - (2 * t899) + t911) * t1033 + (t776 + t1029) * t714) * t735 + ((-0.2e1 * t708 * t827 + t850 * t991 + (-pkin(2) * t824 - t1003) * t735) * t735 - (-t850 * t997 / 0.2e1 + t753 * t776) * t815 * t753) * pkin(3)) * t832;
t705 = (-pkin(3) * t994 + (-t732 + t1036 + (-t830 - t1012) * t735) * t735) * t832;
t890 = m(2) + m(3);
t950 = t864 * t1014;
t915 = pkin(1) * t890 - t950 + t982;
t928 = rSges(3,1) * t827 - rSges(3,2) * t821;
t738 = t928 * m(3) + t915;
t792 = -rSges(3,2) * t1009 + Icges(3,6);
t795 = rSges(3,1) * t1009 - Icges(3,5);
t747 = -t792 * t827 + t795 * t821;
t783 = t1009 + m(2) * (rSges(2,3) + qJ(2,3));
t818 = sin(t854);
t903 = rSges(2,2) ^ 2;
t905 = rSges(2,1) ^ 2;
t910 = -(m(3) * t908 + (-t903 + t905) * m(2) - Icges(2,1) + Icges(2,2)) * t849 / 0.2e1 - (-rSges(2,1) * t1014 + Icges(2,4)) * sin(t896) + t982 * t1025 - (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + t950 * t1037 - Icges(3,2) / 0.2e1 - Icges(2,2) / 0.2e1 - Icges(3,1) / 0.2e1 - Icges(2,1) / 0.2e1 - Icges(1,3);
t894 = 2 * t909;
t916 = -t894 / 0.2e1 - t902 / 0.2e1 - t904 / 0.2e1 - t908 / 0.2e1;
t929 = (2 * rSges(2,3) ^ 2) + t894 + t903 + t905;
t947 = t821 * t994;
t948 = t1015 * t1025;
t949 = m(3) * rSges(3,1) * t1037;
t954 = t753 * t1018;
t960 = t824 + t875;
t963 = t818 + t869;
t966 = t738 * t693 - t747 * t947 + (t802 * t1013 + t814 * t799 + ((qJ(2,3) * t1022) + (2 * t899) + t929) * t1019 + t910 + (-(t846 ^ 2) - t928 * t1037 + (-t960 * rSges(3,1) + t963 * rSges(3,2)) * pkin(2) + t916) * m(3)) * t705 + (-t795 * t827 * t753 - t790 * t997 - t792 * t991) * t753 + (-t949 * t991 + t783 * t1036 + 0.2e1 * (t875 * t1016 + t869 * t1017) * t954 + (-rSges(3,1) * t818 - rSges(3,2) * t824) * t954 + (t827 * t948 - 0.2e1 * t981) * t753) * t735;
t733 = pkin(1) * t736;
t709 = t733 - t720 / 0.2e1 - t719 / 0.2e1 - t742 / 0.2e1;
t855 = t896 + qJ(3,2);
t825 = cos(t855);
t900 = (qJ(2,2) ^ 2);
t694 = ((t709 * t1030 + ((qJ(2,2) * t1020) - t906 * t803 - (2 * t900) + t911) * t1032 + (t777 + t1029) * t715) * t736 + ((-0.2e1 * t709 * t828 + t851 * t990 + (-pkin(2) * t825 - t1002) * t736) * t736 - (-t851 * t996 / 0.2e1 + t754 * t777) * t816 * t754) * pkin(3)) * t833;
t706 = (-pkin(3) * t993 + (-t733 + t1035 + (-t830 - t1011) * t736) * t736) * t833;
t926 = rSges(3,1) * t828 - rSges(3,2) * t822;
t739 = t926 * m(3) + t915;
t793 = -rSges(3,2) * t1008 + Icges(3,6);
t796 = rSges(3,1) * t1008 - Icges(3,5);
t748 = -t793 * t828 + t796 * t822;
t784 = t1008 + m(2) * (rSges(2,3) + qJ(2,2));
t819 = sin(t855);
t946 = t822 * t993;
t953 = t754 * t1018;
t959 = t825 + t877;
t962 = t819 + t871;
t965 = t739 * t694 - t748 * t946 + (t803 * t1013 + t814 * t800 + ((qJ(2,2) * t1022) + (2 * t900) + t929) * t1019 + t910 + (-(t847 ^ 2) - t926 * t1037 + (-t959 * rSges(3,1) + t962 * rSges(3,2)) * pkin(2) + t916) * m(3)) * t706 + (-t796 * t828 * t754 - t790 * t996 - t793 * t990) * t754 + (-t949 * t990 + t784 * t1035 + 0.2e1 * (t877 * t1016 + t871 * t1017) * t953 + (-rSges(3,1) * t819 - rSges(3,2) * t825) * t953 + (t828 * t948 - 0.2e1 * t980) * t754) * t736;
t734 = pkin(1) * t737;
t710 = t734 - t722 / 0.2e1 - t721 / 0.2e1 - t743 / 0.2e1;
t856 = t896 + qJ(3,1);
t826 = cos(t856);
t901 = (qJ(2,1) ^ 2);
t695 = ((t710 * t1030 + ((qJ(2,1) * t1020) - t906 * t804 - (2 * t901) + t911) * t1031 + (t778 + t1029) * t716) * t737 + ((-0.2e1 * t710 * t829 + t852 * t989 + (-pkin(2) * t826 - t1001) * t737) * t737 - (-t852 * t995 / 0.2e1 + t755 * t778) * t817 * t755) * pkin(3)) * t834;
t707 = (-pkin(3) * t992 + (-t734 + t1034 + (-t830 - t1010) * t737) * t737) * t834;
t924 = rSges(3,1) * t829 - rSges(3,2) * t823;
t740 = t924 * m(3) + t915;
t794 = -rSges(3,2) * t1007 + Icges(3,6);
t797 = rSges(3,1) * t1007 - Icges(3,5);
t749 = -t794 * t829 + t797 * t823;
t785 = t1007 + m(2) * (rSges(2,3) + qJ(2,1));
t820 = sin(t856);
t945 = t823 * t992;
t952 = t755 * t1018;
t958 = t826 + t879;
t961 = t820 + t873;
t964 = t740 * t695 - t749 * t945 + (t804 * t1013 + t814 * t801 + ((qJ(2,1) * t1022) + (2 * t901) + t929) * t1019 + t910 + (-(t848 ^ 2) - t924 * t1037 + (-t958 * rSges(3,1) + t961 * rSges(3,2)) * pkin(2) + t916) * m(3)) * t707 + (-t797 * t829 * t755 - t790 * t995 - t794 * t989) * t755 + (-t949 * t989 + t785 * t1034 + 0.2e1 * (t879 * t1016 + t873 * t1017) * t952 + (-rSges(3,1) * t820 - rSges(3,2) * t826) * t952 + (t829 * t948 - 0.2e1 * t979) * t755) * t737;
t938 = t966 * t815;
t937 = t965 * t816;
t936 = t964 * t817;
t935 = -t890 * t693 + t738 * t705 - (-t753 * t1026 + t735 * t783) * t735;
t934 = -t890 * t694 + t739 * t706 - (-t754 * t1027 + t736 * t784) * t736;
t933 = -t890 * t695 + t740 * t707 - (-t755 * t1028 + t737 * t785) * t737;
t922 = t935 * t772;
t921 = t934 * t773;
t920 = t933 * t774;
t798 = -(t902 + t904) * m(3) - Icges(3,3);
t919 = t815 * (t747 * t705 - t798 * t947 + ((t732 - 0.2e1 * t718 - 0.2e1 * t717 - 0.2e1 * t741) * t1026 + (t790 * t799 + 0.2e1 * t981 + (t963 * rSges(3,1) + t960 * rSges(3,2)) * t1018) * t735) * t1033);
t918 = t816 * (t748 * t706 - t798 * t946 + ((t733 - 0.2e1 * t720 - 0.2e1 * t719 - 0.2e1 * t742) * t1027 + (t790 * t800 + 0.2e1 * t980 + (t962 * rSges(3,1) + t959 * rSges(3,2)) * t1018) * t736) * t1032);
t917 = t817 * (t749 * t707 - t798 * t945 + ((t734 - 0.2e1 * t722 - 0.2e1 * t721 - 0.2e1 * t743) * t1028 + (t790 * t801 + 0.2e1 * t979 + (t961 * rSges(3,1) + t958 * rSges(3,2)) * t1018) * t737) * t1031);
t1 = [(t728 * t920 + t767 * t936) * t834 + (t726 * t921 + t765 * t937) * t833 + (t724 * t922 + t763 * t938) * t832 + (t835 * t919 + t836 * t918 + t837 * t917) * t907; (t727 * t920 + t766 * t936) * t834 + (t725 * t921 + t764 * t937) * t833 + (t723 * t922 + t762 * t938) * t832 + (t838 * t919 + t839 * t918 + t840 * t917) * t907; (t933 * t758 + t964 * t880) * t834 + (t934 * t757 + t965 * t878) * t833 + (t935 * t756 + t966 * t876) * t832;];
taucX  = t1;
