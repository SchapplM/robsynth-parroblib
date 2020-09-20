% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RPRRR6V1G1A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x12]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:32
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RPRRR6V1G1A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:32:20
% EndTime: 2020-08-06 18:32:24
% DurationCPUTime: 4.17s
% Computational Cost: add. (21987->351), mult. (19962->634), div. (1515->20), fcn. (9759->68), ass. (0->324)
t1085 = 2 * pkin(2);
t889 = 0.1e1 / pkin(3) ^ 2;
t1084 = t889 / 0.2e1;
t886 = 2 * qJ(3,1);
t861 = sin(t886);
t850 = pkin(3) * t861;
t876 = sin(qJ(3,1));
t780 = t876 * t1085 + t850 + (sin((pkin(7) + qJ(3,1))) + sin((-pkin(7) + qJ(3,1)))) * pkin(1);
t772 = 0.1e1 / t780;
t879 = cos(qJ(3,1));
t1027 = t772 * t879;
t855 = cos(pkin(7)) * pkin(1);
t844 = t855 + pkin(2);
t1001 = t876 * t844;
t847 = (qJ(1,1) + legFrame(1,3));
t830 = pkin(7) + t847;
t810 = sin(t830);
t813 = cos(t830);
t881 = xDP(2);
t882 = xDP(1);
t777 = -t810 * t882 + t813 * t881;
t1060 = t879 * pkin(3) + pkin(2);
t795 = t855 + t1060;
t790 = 0.1e1 / t795;
t1022 = t777 * t790;
t784 = 0.1e1 / (t850 + 0.2e1 * t1001);
t867 = t879 ^ 2;
t887 = pkin(3) ^ 2;
t888 = 0.1e1 / pkin(3);
t883 = (pkin(6) + pkin(5));
t1063 = 2 * t883;
t854 = sin(pkin(7)) * pkin(1);
t891 = pkin(1) ^ 2;
t906 = (pkin(2) ^ 2) + t854 * t1063 + (t883 ^ 2) + t891;
t1066 = 0.2e1 * t813;
t1074 = -0.2e1 * pkin(1);
t785 = -0.2e1 * pkin(2) * t881 - 0.2e1 * t882 * t883;
t792 = -pkin(2) * t882 + t881 * t883;
t833 = sin(t847);
t836 = cos(t847);
t823 = qJ(3,1) + t830;
t824 = -qJ(3,1) + t830;
t993 = cos(t823) + cos(t824);
t996 = sin(t823) + sin(t824);
t747 = t792 * t1066 + t785 * t810 + (t833 * t881 + t836 * t882) * t1074 + (-t881 * t996 - t882 * t993) * pkin(3);
t1037 = t747 * t772;
t954 = t888 * t1037;
t910 = t795 * t954;
t843 = t854 + pkin(5);
t827 = pkin(6) + t843;
t1028 = t772 * t876;
t955 = t747 * t1028;
t913 = t827 * t955;
t948 = t827 * t1022;
t990 = 0.2e1 * t855;
t999 = pkin(3) * t1085;
t726 = -(-t913 + (t1060 * t990 + t867 * t887 + t879 * t999 + t906) * t1022) / (t850 / 0.2e1 + t1001) * t888 * t1022 - 0.2e1 * (-t876 * t948 + t879 * t910) * t784 * t954;
t1083 = t726 * t1027 / 0.2e1;
t1082 = -t726 * t1028 / 0.2e1;
t885 = 2 * qJ(3,2);
t860 = sin(t885);
t849 = pkin(3) * t860;
t875 = sin(qJ(3,2));
t779 = t875 * t1085 + t849 + (sin((pkin(7) + qJ(3,2))) + sin((-pkin(7) + qJ(3,2)))) * pkin(1);
t769 = 0.1e1 / t779;
t878 = cos(qJ(3,2));
t1029 = t769 * t878;
t1003 = t875 * t844;
t846 = (qJ(1,2) + legFrame(2,3));
t829 = pkin(7) + t846;
t809 = sin(t829);
t812 = cos(t829);
t776 = -t809 * t882 + t812 * t881;
t1061 = t878 * pkin(3) + pkin(2);
t794 = t855 + t1061;
t788 = 0.1e1 / t794;
t1024 = t776 * t788;
t783 = 0.1e1 / (t849 + 0.2e1 * t1003);
t866 = t878 ^ 2;
t1067 = 0.2e1 * t812;
t832 = sin(t846);
t835 = cos(t846);
t819 = qJ(3,2) + t829;
t820 = -qJ(3,2) + t829;
t994 = cos(t819) + cos(t820);
t997 = sin(t819) + sin(t820);
t746 = t792 * t1067 + t785 * t809 + (t832 * t881 + t835 * t882) * t1074 + (-t881 * t997 - t882 * t994) * pkin(3);
t1039 = t746 * t769;
t956 = t888 * t1039;
t914 = t794 * t956;
t1030 = t769 * t875;
t957 = t746 * t1030;
t917 = t827 * t957;
t949 = t827 * t1024;
t725 = -(-t917 + (t1061 * t990 + t866 * t887 + t878 * t999 + t906) * t1024) / (t849 / 0.2e1 + t1003) * t888 * t1024 - 0.2e1 * (-t875 * t949 + t878 * t914) * t783 * t956;
t1081 = t725 * t1029 / 0.2e1;
t1080 = -t725 * t1030 / 0.2e1;
t884 = 2 * qJ(3,3);
t859 = sin(t884);
t848 = pkin(3) * t859;
t874 = sin(qJ(3,3));
t778 = t874 * t1085 + t848 + (sin((pkin(7) + qJ(3,3))) + sin((-pkin(7) + qJ(3,3)))) * pkin(1);
t766 = 0.1e1 / t778;
t877 = cos(qJ(3,3));
t1031 = t766 * t877;
t1005 = t874 * t844;
t845 = (qJ(1,3) + legFrame(3,3));
t828 = pkin(7) + t845;
t808 = sin(t828);
t811 = cos(t828);
t775 = -t808 * t882 + t811 * t881;
t1062 = t877 * pkin(3) + pkin(2);
t793 = t855 + t1062;
t786 = 0.1e1 / t793;
t1026 = t775 * t786;
t782 = 0.1e1 / (t848 + 0.2e1 * t1005);
t865 = t877 ^ 2;
t1068 = 0.2e1 * t811;
t831 = sin(t845);
t834 = cos(t845);
t815 = qJ(3,3) + t828;
t816 = -qJ(3,3) + t828;
t995 = cos(t815) + cos(t816);
t998 = sin(t815) + sin(t816);
t745 = t792 * t1068 + t785 * t808 + (t831 * t881 + t834 * t882) * t1074 + (-t881 * t998 - t882 * t995) * pkin(3);
t1041 = t745 * t766;
t958 = t888 * t1041;
t918 = t793 * t958;
t1032 = t766 * t874;
t959 = t745 * t1032;
t921 = t827 * t959;
t950 = t827 * t1026;
t724 = -(-t921 + (t1062 * t990 + t865 * t887 + t877 * t999 + t906) * t1026) / (t848 / 0.2e1 + t1005) * t888 * t1026 - 0.2e1 * (-t874 * t950 + t877 * t918) * t782 * t958;
t1079 = t724 * t1031 / 0.2e1;
t1078 = -t724 * t1032 / 0.2e1;
t1012 = t790 * t844;
t781 = pkin(2) * t990 + t887 / 0.2e1 + t906;
t723 = (-t1022 * (0.2e1 * (-t1022 * t781 + t913) * t879 + (-cos(t886) + (-0.2e1 * t867 - 0.1e1) * t1012) * t777 * pkin(3)) + 0.2e1 * (-t861 * t948 / 0.2e1 + t910) * t1037) * t784;
t1077 = t723 * t772 / 0.2e1;
t1015 = t788 * t844;
t722 = (-t1024 * (0.2e1 * (-t1024 * t781 + t917) * t878 + (-cos(t885) + (-0.2e1 * t866 - 0.1e1) * t1015) * t776 * pkin(3)) + 0.2e1 * (-t860 * t949 / 0.2e1 + t914) * t1039) * t783;
t1076 = t722 * t769 / 0.2e1;
t1018 = t786 * t844;
t721 = (-t1026 * (0.2e1 * (-t1026 * t781 + t921) * t877 + (-cos(t884) + (-0.2e1 * t865 - 0.1e1) * t1018) * t775 * pkin(3)) + 0.2e1 * (-t859 * t950 / 0.2e1 + t918) * t1041) * t782;
t1075 = t721 * t766 / 0.2e1;
t1073 = -2 * pkin(2);
t767 = 0.1e1 / t778 ^ 2;
t770 = 0.1e1 / t779 ^ 2;
t773 = 0.1e1 / t780 ^ 2;
t1071 = -0.2e1 * t808;
t1070 = -0.2e1 * t809;
t1069 = -0.2e1 * t810;
t1065 = 0.2e1 * t844;
t1064 = -2 * t883;
t851 = 0.2e1 * t865 - 0.1e1;
t852 = 0.2e1 * t866 - 0.1e1;
t853 = 0.2e1 * t867 - 0.1e1;
t1056 = t724 * t766;
t1055 = t724 * t874;
t1054 = t724 * t877;
t1053 = t725 * t769;
t1052 = t725 * t875;
t1051 = t725 * t878;
t1050 = t726 * t772;
t1049 = t726 * t876;
t1048 = t726 * t879;
t742 = t745 ^ 2;
t1047 = t742 * t874;
t1046 = t742 * t877;
t743 = t746 ^ 2;
t1045 = t743 * t875;
t1044 = t743 * t878;
t744 = t747 ^ 2;
t1043 = t744 * t876;
t1042 = t744 * t879;
t1040 = t745 * t775;
t1038 = t746 * t776;
t1036 = t747 * t777;
t763 = t775 ^ 2;
t787 = 0.1e1 / t793 ^ 2;
t1035 = t763 * t787;
t764 = t776 ^ 2;
t789 = 0.1e1 / t794 ^ 2;
t1034 = t764 * t789;
t765 = t777 ^ 2;
t791 = 0.1e1 / t795 ^ 2;
t1033 = t765 * t791;
t1025 = t775 * t787;
t1023 = t776 * t789;
t1021 = t777 * t791;
t1020 = t786 * t808;
t1019 = t786 * t811;
t1017 = t788 * t809;
t1016 = t788 * t812;
t1014 = t790 * t810;
t1013 = t790 * t813;
t1011 = t843 * t874;
t1010 = t843 * t875;
t1009 = t843 * t876;
t1008 = t843 * t877;
t1007 = t843 * t878;
t1006 = t843 * t879;
t1004 = t874 * t877;
t1002 = t875 * t878;
t1000 = t876 * t879;
t992 = 0.2e1 * pkin(1);
t987 = t786 * t1055;
t986 = t786 * t1054;
t983 = t788 * t1052;
t982 = t788 * t1051;
t979 = t790 * t1049;
t978 = t790 * t1048;
t938 = t883 + t854;
t727 = (-t1026 * t938 + 0.2e1 * t959) * t1025;
t977 = t727 * t1032;
t976 = t727 * t1031;
t728 = (-t1024 * t938 + 0.2e1 * t957) * t1023;
t975 = t728 * t1030;
t974 = t728 * t1029;
t729 = (-t1022 * t938 + 0.2e1 * t955) * t1021;
t973 = t729 * t1028;
t972 = t729 * t1027;
t971 = t767 * t1047;
t970 = t767 * t1046;
t768 = t766 * t767;
t969 = t768 * t1047;
t968 = t768 * t1046;
t967 = t770 * t1045;
t966 = t770 * t1044;
t771 = t769 * t770;
t965 = t771 * t1045;
t964 = t771 * t1044;
t963 = t773 * t1043;
t962 = t773 * t1042;
t774 = t772 * t773;
t961 = t774 * t1043;
t960 = t774 * t1042;
t953 = t844 * t1035;
t952 = t844 * t1034;
t951 = t844 * t1033;
t947 = t727 * t1020;
t946 = t728 * t1017;
t945 = t729 * t1014;
t944 = t727 * t1019;
t943 = t728 * t1016;
t942 = t729 * t1013;
t941 = t775 * t1018;
t940 = t776 * t1015;
t939 = t777 * t1012;
t937 = 0.2e1 * t745 * t1020;
t936 = -0.2e1 * t745 * t1019;
t935 = 0.2e1 * t746 * t1017;
t934 = -0.2e1 * t746 * t1016;
t933 = 0.2e1 * t747 * t1014;
t932 = -0.2e1 * t747 * t1013;
t930 = t727 * t786 * t1004;
t929 = t728 * t788 * t1002;
t928 = t729 * t790 * t1000;
t927 = t786 * t971;
t926 = t786 * t970;
t925 = t788 * t967;
t924 = t788 * t966;
t923 = t790 * t963;
t922 = t790 * t962;
t920 = t843 * t958;
t919 = t745 * t851 * t1025;
t916 = t843 * t956;
t915 = t746 * t852 * t1023;
t912 = t843 * t954;
t911 = t747 * t853 * t1021;
t909 = t766 * t787 * t1004;
t908 = t769 * t789 * t1002;
t907 = t772 * t791 * t1000;
t905 = -t1009 * t729 + t879 * t723 + t876 * t951;
t904 = -t1006 * t729 - t876 * t723 + t879 * t951;
t903 = -t1011 * t727 + t877 * t721 + t874 * t953;
t902 = -t1008 * t727 - t874 * t721 + t877 * t953;
t901 = -t1010 * t728 + t878 * t722 + t875 * t952;
t900 = -t1007 * t728 - t875 * t722 + t878 * t952;
t899 = -t945 - t946 - t947;
t898 = t942 + t943 + t944;
t864 = t876 ^ 2;
t863 = t875 ^ 2;
t862 = t874 ^ 2;
t842 = -qJ(3,1) + t847;
t841 = qJ(3,1) + t847;
t840 = -qJ(3,2) + t846;
t839 = qJ(3,2) + t846;
t838 = -qJ(3,3) + t845;
t837 = qJ(3,3) + t845;
t825 = -2 * qJ(3,1) + t830;
t822 = t886 + t830;
t821 = -2 * qJ(3,2) + t829;
t818 = t885 + t829;
t817 = -2 * qJ(3,3) + t828;
t814 = t884 + t828;
t762 = -pkin(3) * t993 + t810 * t1064 + t813 * t1073 + t836 * t1074;
t761 = -pkin(3) * t994 + t809 * t1064 + t812 * t1073 + t835 * t1074;
t760 = -pkin(3) * t995 + t808 * t1064 + t811 * t1073 + t834 * t1074;
t759 = pkin(2) * t1069 - pkin(3) * t996 + t813 * t1063 + t833 * t1074;
t758 = pkin(2) * t1070 - pkin(3) * t997 + t812 * t1063 + t832 * t1074;
t757 = pkin(2) * t1071 - pkin(3) * t998 + t811 * t1063 + t831 * t1074;
t756 = t853 * t1033;
t755 = t852 * t1034;
t754 = t851 * t1035;
t753 = t996 * t1085 + (sin(t842) + sin(t841)) * t992 + t993 * t1064 + (sin(t825) + sin(t822) + 0.2e1 * t810) * pkin(3);
t752 = t997 * t1085 + (sin(t840) + sin(t839)) * t992 + t994 * t1064 + (sin(t821) + sin(t818) + 0.2e1 * t809) * pkin(3);
t751 = t998 * t1085 + (sin(t838) + sin(t837)) * t992 + t995 * t1064 + (sin(t817) + sin(t814) + 0.2e1 * t808) * pkin(3);
t750 = t993 * t1085 + (cos(t842) + cos(t841)) * t992 + t996 * t1063 + (cos(t825) + cos(t822) + t1066) * pkin(3);
t749 = t994 * t1085 + (cos(t840) + cos(t839)) * t992 + t997 * t1063 + (cos(t821) + cos(t818) + t1067) * pkin(3);
t748 = t995 * t1085 + (cos(t838) + cos(t837)) * t992 + t998 * t1063 + (cos(t817) + cos(t814) + t1068) * pkin(3);
t741 = t879 * t939 - t876 * t912 / 0.2e1;
t740 = t878 * t940 - t875 * t916 / 0.2e1;
t739 = t877 * t941 - t874 * t920 / 0.2e1;
t738 = t879 * t912 / 0.2e1 + t876 * t939;
t737 = t878 * t916 / 0.2e1 + t875 * t940;
t736 = t877 * t920 / 0.2e1 + t874 * t941;
t720 = t879 * t729 * t1065 - t1009 * t726;
t719 = t877 * t727 * t1065 - t1011 * t724;
t718 = t728 * t878 * t1065 - t1010 * t725;
t717 = -0.2e1 * t1001 * t729 - t1006 * t726;
t716 = -0.2e1 * t1003 * t728 - t1007 * t725;
t715 = -0.2e1 * t1005 * t727 - t1008 * t724;
t1 = [t899, 0, 0, t748 * t1075 + t749 * t1076 + t750 * t1077 + t899 * t891, -t862 * t947 - t863 * t946 - t864 * t945 + ((t1036 * t1069 - t762 * t765) * t907 + (t1038 * t1070 - t761 * t764) * t908 + (t1040 * t1071 - t760 * t763) * t909) * t888, t930 * t1071 + t929 * t1070 + t928 * t1069 + ((t911 * t1069 - t756 * t762) * t772 + (t915 * t1070 - t755 * t761) * t769 + (t919 * t1071 - t754 * t760) * t766) * t888, -t808 * t987 - t809 * t983 - t810 * t979 + (-t808 * t926 - t809 * t924 - t810 * t922) * t889 + (t760 * t977 + t761 * t975 + t762 * t973) * t888, -t808 * t986 - t809 * t982 - t810 * t978 + (t808 * t927 + t809 * t925 + t810 * t923) * t889 + (t760 * t976 + t761 * t974 + t762 * t972) * t888, (t1050 * t762 + t1053 * t761 + t1056 * t760) * t888, -t718 * t1017 - t719 * t1020 - t720 * t1014 + t748 * t1079 + t749 * t1081 + t750 * t1083 + (-t748 * t969 - t749 * t965 - t750 * t961) * t1084 + ((t738 * t933 + t762 * t905) * t772 + (t737 * t935 + t761 * t901) * t769 + (t736 * t937 + t760 * t903) * t766) * t888, -t715 * t1020 - t716 * t1017 - t717 * t1014 + t748 * t1078 + t749 * t1080 + t750 * t1082 + (-t748 * t968 - t749 * t964 - t750 * t960) * t1084 + ((t741 * t933 + t762 * t904) * t772 + (t740 * t935 + t761 * t900) * t769 + (t739 * t937 + t760 * t902) * t766) * t888, 0; t898, 0, 0, t751 * t1075 + t752 * t1076 + t753 * t1077 + t898 * t891, t862 * t944 + t863 * t943 + t864 * t942 + ((t1036 * t1066 - t759 * t765) * t907 + (t1038 * t1067 - t758 * t764) * t908 + (t1040 * t1068 - t757 * t763) * t909) * t888, t930 * t1068 + t929 * t1067 + t928 * t1066 + ((t911 * t1066 - t756 * t759) * t772 + (t915 * t1067 - t755 * t758) * t769 + (t919 * t1068 - t754 * t757) * t766) * t888, t811 * t987 + t812 * t983 + t813 * t979 + (t811 * t926 + t812 * t924 + t813 * t922) * t889 + (t757 * t977 + t758 * t975 + t759 * t973) * t888, t811 * t986 + t812 * t982 + t813 * t978 + (-t811 * t927 - t812 * t925 - t813 * t923) * t889 + (t757 * t976 + t758 * t974 + t759 * t972) * t888, (t1050 * t759 + t1053 * t758 + t1056 * t757) * t888, t718 * t1016 + t719 * t1019 + t720 * t1013 + t751 * t1079 + t752 * t1081 + t753 * t1083 + (-t751 * t969 - t752 * t965 - t753 * t961) * t1084 + ((t738 * t932 + t759 * t905) * t772 + (t737 * t934 + t758 * t901) * t769 + (t736 * t936 + t757 * t903) * t766) * t888, t715 * t1019 + t716 * t1016 + t717 * t1013 + t751 * t1078 + t752 * t1080 + t753 * t1082 + (-t751 * t968 - t752 * t964 - t753 * t960) * t1084 + ((t741 * t932 + t759 * t904) * t772 + (t740 * t934 + t758 * t900) * t769 + (t739 * t936 + t757 * t902) * t766) * t888, 0; 0, 0, 0, t721 + t722 + t723, 0, 0, 0, 0, 0, t1054 + t1051 + t1048 + (-t963 - t967 - t971) * t889, -t1055 - t1052 - t1049 + (-t962 - t966 - t970) * t889, 0;];
tau_reg  = t1;
