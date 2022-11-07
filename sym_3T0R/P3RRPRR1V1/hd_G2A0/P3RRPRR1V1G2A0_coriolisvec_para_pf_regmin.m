% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RRPRR1V1G2A0
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x15]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2022-11-04 17:08
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRPRR1V1G2A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_regmin: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-04 17:08:17
% EndTime: 2022-11-04 17:08:20
% DurationCPUTime: 3.50s
% Computational Cost: add. (11618->353), mult. (18579->793), div. (2535->21), fcn. (14808->18), ass. (0->322)
t936 = (pkin(1) ^ 2);
t1153 = -2 * t936;
t914 = pkin(3) + qJ(3,3);
t884 = 0.1e1 / t914 ^ 2;
t926 = cos(qJ(2,3));
t896 = 0.1e1 / t926;
t1094 = t884 * t896;
t915 = pkin(3) + qJ(3,2);
t887 = 0.1e1 / t915 ^ 2;
t928 = cos(qJ(2,2));
t901 = 0.1e1 / t928;
t1092 = t887 * t901;
t916 = pkin(3) + qJ(3,1);
t890 = 0.1e1 / t916 ^ 2;
t930 = cos(qJ(2,1));
t906 = 0.1e1 / t930;
t1090 = t890 * t906;
t920 = sin(qJ(2,3));
t1089 = t896 * t920;
t922 = sin(qJ(2,2));
t1087 = t901 * t922;
t924 = sin(qJ(2,1));
t1085 = t906 * t924;
t1152 = 2 * pkin(1);
t883 = 0.1e1 / t914;
t886 = 0.1e1 / t915;
t889 = 0.1e1 / t916;
t935 = pkin(1) + pkin(2);
t911 = 1 / t935;
t917 = legFrame(3,2);
t877 = sin(t917);
t880 = cos(t917);
t933 = xDP(2);
t934 = xDP(1);
t868 = t877 * t934 + t880 * t933;
t1080 = t920 * t868;
t921 = sin(qJ(1,3));
t927 = cos(qJ(1,3));
t932 = xDP(3);
t850 = (t927 * t932 + (-t877 * t933 + t880 * t934) * t921) * t926 + t1080;
t1151 = 0.2e1 * t850;
t918 = legFrame(2,2);
t878 = sin(t918);
t881 = cos(t918);
t869 = t878 * t934 + t881 * t933;
t1077 = t922 * t869;
t923 = sin(qJ(1,2));
t929 = cos(qJ(1,2));
t851 = (t929 * t932 + (-t878 * t933 + t881 * t934) * t923) * t928 + t1077;
t1150 = 0.2e1 * t851;
t919 = legFrame(1,2);
t879 = sin(t919);
t882 = cos(t919);
t870 = t879 * t934 + t882 * t933;
t1074 = t924 * t870;
t925 = sin(qJ(1,1));
t931 = cos(qJ(1,1));
t852 = (t931 * t932 + (-t879 * t933 + t882 * t934) * t925) * t930 + t1074;
t1149 = 0.2e1 * t852;
t895 = t926 ^ 2;
t1148 = 0.2e1 * t895;
t900 = t928 ^ 2;
t1147 = 0.2e1 * t900;
t905 = t930 ^ 2;
t1146 = 0.2e1 * t905;
t1145 = t911 / 0.4e1;
t897 = 0.1e1 / t926 ^ 2;
t1144 = t897 - 0.2e1;
t902 = 0.1e1 / t928 ^ 2;
t1143 = t902 - 0.2e1;
t907 = 0.1e1 / t930 ^ 2;
t1142 = t907 - 0.2e1;
t1141 = pkin(1) * t911;
t1037 = t850 * t1094;
t1078 = t921 * t926;
t1083 = t914 * t927;
t1071 = t926 * t935;
t871 = t927 * t1071 + t914 * t921;
t838 = (-t934 * t1083 + (t934 * t1078 + t920 * t933) * t935) * t880 + (t933 * t1083 - (t933 * t1078 - t920 * t934) * t935) * t877 + t871 * t932;
t1125 = t838 * t850;
t865 = t868 ^ 2;
t885 = t883 * t884;
t898 = t896 / t895;
t805 = -t885 * t896 * t1125 + (-(-t850 * t935 + t838) * t1037 + t911 * t865 * t898) * t883;
t1140 = t805 * qJ(3,3);
t1139 = t805 * t920;
t1138 = t805 * t926;
t1034 = t851 * t1092;
t1075 = t923 * t928;
t1082 = t915 * t929;
t1070 = t928 * t935;
t872 = t929 * t1070 + t915 * t923;
t839 = (-t934 * t1082 + (t934 * t1075 + t922 * t933) * t935) * t881 + (t933 * t1082 - (t933 * t1075 - t922 * t934) * t935) * t878 + t872 * t932;
t1122 = t839 * t851;
t866 = t869 ^ 2;
t888 = t886 * t887;
t903 = t901 / t900;
t806 = -t888 * t901 * t1122 + (-(-t851 * t935 + t839) * t1034 + t911 * t866 * t903) * t886;
t1137 = t806 * qJ(3,2);
t1136 = t806 * t922;
t1135 = t806 * t928;
t1031 = t852 * t1090;
t1072 = t925 * t930;
t1081 = t916 * t931;
t1069 = t930 * t935;
t873 = t931 * t1069 + t916 * t925;
t840 = (-t934 * t1081 + (t934 * t1072 + t924 * t933) * t935) * t882 + (t933 * t1081 - (t933 * t1072 - t924 * t934) * t935) * t879 + t873 * t932;
t1119 = t840 * t852;
t867 = t870 ^ 2;
t891 = t889 * t890;
t908 = t906 / t905;
t807 = -t891 * t906 * t1119 + (-(-t852 * t935 + t840) * t1031 + t911 * t867 * t908) * t889;
t1134 = t807 * qJ(3,1);
t1133 = t807 * t924;
t1132 = t807 * t930;
t1124 = t838 * t883;
t1113 = t850 * t883;
t844 = pkin(1) * t1113;
t1131 = (t844 - t1124) * t883;
t1121 = t839 * t886;
t1112 = t851 * t886;
t845 = pkin(1) * t1112;
t1130 = (t845 - t1121) * t886;
t1118 = t840 * t889;
t1111 = t852 * t889;
t846 = pkin(1) * t1111;
t1129 = (t846 - t1118) * t889;
t847 = t850 ^ 2;
t1116 = t847 * t884;
t1128 = t1144 * t1116 * t896;
t848 = t851 ^ 2;
t1115 = t848 * t887;
t1127 = t1143 * t1115 * t901;
t849 = t852 ^ 2;
t1114 = t849 * t890;
t1126 = t1142 * t1114 * t906;
t1123 = t838 * t935;
t1120 = t839 * t935;
t1117 = t840 * t935;
t859 = -t877 * t1078 + t880 * t920;
t1110 = t859 * t896;
t860 = -t878 * t1075 + t881 * t922;
t1109 = t860 * t901;
t861 = -t879 * t1072 + t882 * t924;
t1108 = t861 * t906;
t862 = t880 * t1078 + t877 * t920;
t1107 = t862 * t896;
t863 = t881 * t1075 + t878 * t922;
t1106 = t863 * t901;
t864 = t882 * t1072 + t879 * t924;
t1105 = t864 * t906;
t912 = 1 / t935 ^ 2;
t1104 = t865 * t912;
t1103 = t866 * t912;
t1102 = t867 * t912;
t1101 = t868 * t883;
t1100 = t868 * t884;
t1099 = t869 * t886;
t1098 = t869 * t887;
t1097 = t870 * t889;
t1096 = t870 * t890;
t1095 = t883 * t927;
t1093 = t886 * t929;
t1091 = t889 * t931;
t1088 = t897 * t920;
t1086 = t902 * t922;
t1084 = t907 * t924;
t1079 = t920 * t935;
t1076 = t922 * t935;
t1073 = t924 * t935;
t910 = t935 ^ 2;
t793 = (-t850 * t920 + t868) * t896 ^ 2 * t1101 + (-(t914 * (-t850 + t1080) * t896 + (-t850 * t910 + t1123) * t926 * t883) * t1094 - t885 * t1123) * t850;
t1029 = t898 * t1104;
t892 = t920 ^ 2;
t984 = t892 * t1029;
t955 = -qJ(3,3) * t984 - t793 * t926;
t1061 = t884 * t1151;
t960 = -pkin(1) * t1104 + t838 * t1061;
t1023 = t911 * t1080;
t975 = t1023 * t1113;
t1068 = (qJ(3,3) ^ 2 + t895 * t936) * t805 + t955 * pkin(1) + (t960 * qJ(3,3) + t975 * t1153) * t896;
t794 = (-t851 * t922 + t869) * t901 ^ 2 * t1099 + (-(t915 * (-t851 + t1077) * t901 + (-t851 * t910 + t1120) * t928 * t886) * t1092 - t888 * t1120) * t851;
t1027 = t903 * t1103;
t893 = t922 ^ 2;
t983 = t893 * t1027;
t956 = -qJ(3,2) * t983 - t794 * t928;
t1059 = t887 * t1150;
t959 = -pkin(1) * t1103 + t839 * t1059;
t1022 = t911 * t1077;
t974 = t1022 * t1112;
t1067 = (qJ(3,2) ^ 2 + t900 * t936) * t806 + t956 * pkin(1) + (t959 * qJ(3,2) + t974 * t1153) * t901;
t795 = (-t852 * t924 + t870) * t906 ^ 2 * t1097 + (-(t916 * (-t852 + t1074) * t906 + (-t852 * t910 + t1117) * t930 * t889) * t1090 - t891 * t1117) * t852;
t1025 = t908 * t1102;
t894 = t924 ^ 2;
t982 = t894 * t1025;
t957 = -qJ(3,1) * t982 - t795 * t930;
t1057 = t890 * t1149;
t958 = -pkin(1) * t1102 + t840 * t1057;
t1021 = t911 * t1074;
t973 = t1021 * t1111;
t1066 = (qJ(3,1) ^ 2 + t905 * t936) * t807 + t957 * pkin(1) + (t958 * qJ(3,1) + t973 * t1153) * t906;
t1053 = pkin(1) * t1138;
t1065 = t793 - t1053 + (-qJ(3,3) * t1116 + t975 * t1152) * t897;
t1052 = pkin(1) * t1135;
t1064 = t794 - t1052 + (-qJ(3,2) * t1115 + t974 * t1152) * t902;
t1051 = pkin(1) * t1132;
t1063 = t795 - t1051 + (-qJ(3,1) * t1114 + t973 * t1152) * t907;
t1062 = t868 * t1151;
t1060 = t869 * t1150;
t1058 = t870 * t1149;
t1056 = -0.4e1 * t1095;
t1055 = -0.4e1 * t1093;
t1054 = -0.4e1 * t1091;
t1050 = t805 * t883 * t896;
t1049 = t805 * t1095;
t1048 = t805 * t1089;
t1047 = t806 * t886 * t901;
t1046 = t806 * t1093;
t1045 = t806 * t1087;
t1044 = t807 * t889 * t906;
t1043 = t807 * t1091;
t1042 = t807 * t1085;
t1041 = t847 * t885 * t897;
t1040 = t848 * t888 * t902;
t1039 = t849 * t891 * t907;
t1038 = t850 * t1100;
t1036 = t850 * t1088;
t1035 = t851 * t1098;
t1033 = t851 * t1086;
t1032 = t852 * t1096;
t1030 = t852 * t1084;
t899 = 0.1e1 / t895 ^ 2;
t1028 = t865 * t899 * t920;
t904 = 0.1e1 / t900 ^ 2;
t1026 = t866 * t904 * t922;
t909 = 0.1e1 / t905 ^ 2;
t1024 = t867 * t909 * t924;
t1020 = t883 * t1139;
t1019 = t884 * t1088;
t1018 = t886 * t1136;
t1017 = t887 * t1086;
t1016 = t889 * t1133;
t1015 = t890 * t1084;
t1014 = t883 * (-pkin(1) * t984 + t960 * t896 + 0.2e1 * t1140);
t1013 = t886 * (-pkin(1) * t983 + t959 * t901 + 0.2e1 * t1137);
t1012 = t889 * (-pkin(1) * t982 + t958 * t906 + 0.2e1 * t1134);
t1011 = t1068 * t896;
t1010 = t1067 * t901;
t1009 = t1066 * t906;
t835 = (-qJ(3,3) * t1023 / 0.4e1 + (t895 - 0.1e1 / 0.2e1) * t844) * t896;
t1008 = -0.4e1 * t835 * t883 * t897;
t836 = (-qJ(3,2) * t1022 / 0.4e1 + (t900 - 0.1e1 / 0.2e1) * t845) * t901;
t1007 = -0.4e1 * t836 * t886 * t902;
t837 = (-qJ(3,1) * t1021 / 0.4e1 + (t905 - 0.1e1 / 0.2e1) * t846) * t906;
t1006 = -0.4e1 * t837 * t889 * t907;
t841 = (qJ(3,3) * t868 * t1145 + t920 * t844) * t896;
t1005 = -0.4e1 * t841 * t1101;
t842 = (qJ(3,2) * t869 * t1145 + t922 * t845) * t901;
t1004 = -0.4e1 * t842 * t1099;
t843 = (qJ(3,1) * t870 * t1145 + t924 * t846) * t906;
t1003 = -0.4e1 * t843 * t1097;
t1002 = t871 * t1061;
t1001 = 0.2e1 * t1037;
t1000 = t872 * t1059;
t999 = 0.2e1 * t1034;
t998 = t873 * t1057;
t997 = 0.2e1 * t1031;
t996 = 0.2e1 * t1020;
t995 = 0.2e1 * t1018;
t994 = 0.2e1 * t1016;
t993 = pkin(1) * t1029;
t992 = pkin(1) * t1027;
t991 = pkin(1) * t1025;
t990 = t892 * t1050;
t989 = t893 * t1047;
t988 = t894 * t1044;
t987 = t927 * t1038;
t986 = t929 * t1035;
t985 = t931 * t1032;
t981 = t896 * t1014;
t980 = t901 * t1013;
t979 = t906 * t1012;
t978 = 0.2e1 * t1036;
t977 = 0.2e1 * t1033;
t976 = 0.2e1 * t1030;
t972 = (-t1142 * t849 * pkin(1) - 0.2e1 * t1119) * t1090 - t1134;
t971 = (-t1143 * t848 * pkin(1) - 0.2e1 * t1122) * t1092 - t1137;
t970 = (-t1144 * t847 * pkin(1) - 0.2e1 * t1125) * t1094 - t1140;
t874 = t1148 - 0.1e1;
t969 = 0.2e1 * t874 * t898 * t1038;
t875 = t1147 - 0.1e1;
t968 = 0.2e1 * t875 * t903 * t1035;
t876 = t1146 - 0.1e1;
t967 = 0.2e1 * t876 * t908 * t1032;
t966 = t865 * t883 * (t892 * t899 + t897);
t965 = t866 * t886 * (t893 * t904 + t902);
t964 = t867 * t889 * (t894 * t909 + t907);
t963 = t921 * t1071 - t1083;
t962 = t923 * t1070 - t1082;
t961 = t925 * t1069 - t1081;
t954 = (t844 - 0.2e1 * t1124) * t883 * t1036 + (t993 - t1140) * t1089;
t953 = (t845 - 0.2e1 * t1121) * t886 * t1033 + (t992 - t1137) * t1087;
t952 = (t846 - 0.2e1 * t1118) * t889 * t1030 + (t991 - t1134) * t1085;
t951 = t879 * t1042 + t878 * t1045 + t877 * t1048;
t950 = t882 * t1042 + t881 * t1045 + t880 * t1048;
t913 = t911 / t910;
t858 = t879 * t1073 + t961 * t882;
t857 = t878 * t1076 + t962 * t881;
t856 = t877 * t1079 + t963 * t880;
t855 = t881 * t1076 - t962 * t878;
t854 = t882 * t1073 - t961 * t879;
t853 = t880 * t1079 - t963 * t877;
t798 = (0.2e1 * t991 - t1134) * t924;
t797 = (0.2e1 * t992 - t1137) * t922;
t796 = (0.2e1 * t993 - t1140) * t920;
t789 = (-t907 * qJ(3,1) * t1102 - 0.2e1 * t1051 + t795) * t924;
t788 = (-t902 * qJ(3,2) * t1103 - 0.2e1 * t1052 + t794) * t922;
t787 = (-t897 * qJ(3,3) * t1104 - 0.2e1 * t1053 + t793) * t920;
t786 = pkin(1) * t807 * t1146 + t957;
t785 = pkin(1) * t806 * t1147 + t956;
t784 = pkin(1) * t805 * t1148 + t955;
t1 = [t864 * t1044 + t863 * t1047 + t862 * t1050, 0, 0, t862 * t990 + t863 * t989 + t864 * t988 + ((t864 * t1058 - t849 * t879) * t1015 + (t863 * t1060 - t848 * t878) * t1017 + (t862 * t1062 - t847 * t877) * t1019) * t911, t862 * t996 + t863 * t995 + t864 * t994 + (t879 * t1126 + t878 * t1127 + t877 * t1128 + t862 * t969 + t863 * t968 + t864 * t967) * t911, t951 * t911 + (t862 * t966 + t863 * t965 + t864 * t964) * t912, (t805 * t877 + t806 * t878 + t807 * t879) * t911, (t879 * t1024 + t878 * t1026 + t877 * t1028) * t913, 0, 0, (t786 * t1105 - t858 * t1132) * t889 + (t785 * t1106 - t857 * t1135) * t886 + (t784 * t1107 - t856 * t1138) * t883 + ((t864 * t1003 + t879 * t798) * t906 + (t863 * t1004 + t878 * t797) * t901 + (t862 * t1005 + t877 * t796) * t896 + (t858 * t1096 + t879 * t1129) * t976 + (t857 * t1098 + t878 * t1130) * t977 + (t856 * t1100 + t877 * t1131) * t978) * t911, (t789 * t1105 + t858 * t1133) * t889 + (t788 * t1106 + t857 * t1136) * t886 + (t787 * t1107 + t856 * t1139) * t883 + (t972 * t879 + t971 * t878 + t970 * t877 + (t1006 * t864 + t858 * t997) * t870 + (t1007 * t863 + t857 * t999) * t869 + (t1001 * t856 + t1008 * t862) * t868) * t911, -t1039 * t858 - t1040 * t857 - t1041 * t856 - t951 * t1141 + t862 * t981 + t863 * t980 + t864 * t979, (t1009 * t864 + t1063 * t858) * t889 + (t1010 * t863 + t1064 * t857) * t886 + (t1011 * t862 + t1065 * t856) * t883 + (t877 * t954 + t878 * t953 + t879 * t952) * t1141, 0; t861 * t1044 + t860 * t1047 + t859 * t1050, 0, 0, t859 * t990 + t860 * t989 + t861 * t988 + ((t861 * t1058 - t849 * t882) * t1015 + (t860 * t1060 - t848 * t881) * t1017 + (t859 * t1062 - t847 * t880) * t1019) * t911, t859 * t996 + t860 * t995 + t861 * t994 + (t882 * t1126 + t881 * t1127 + t880 * t1128 + t859 * t969 + t860 * t968 + t861 * t967) * t911, t950 * t911 + (t859 * t966 + t860 * t965 + t861 * t964) * t912, (t805 * t880 + t806 * t881 + t807 * t882) * t911, (t882 * t1024 + t881 * t1026 + t880 * t1028) * t913, 0, 0, (t786 * t1108 - t854 * t1132) * t889 + (t785 * t1109 - t855 * t1135) * t886 + (t784 * t1110 - t853 * t1138) * t883 + ((t1003 * t861 + t882 * t798) * t906 + (t1004 * t860 + t881 * t797) * t901 + (t1005 * t859 + t880 * t796) * t896 + (t854 * t1096 + t882 * t1129) * t976 + (t855 * t1098 + t881 * t1130) * t977 + (t853 * t1100 + t880 * t1131) * t978) * t911, (t789 * t1108 + t854 * t1133) * t889 + (t788 * t1109 + t855 * t1136) * t886 + (t787 * t1110 + t853 * t1139) * t883 + (t972 * t882 + t971 * t881 + t970 * t880 + (t1006 * t861 + t854 * t997) * t870 + (t1007 * t860 + t855 * t999) * t869 + (t1001 * t853 + t1008 * t859) * t868) * t911, -t1039 * t854 - t1040 * t855 - t1041 * t853 - t950 * t1141 + t859 * t981 + t860 * t980 + t861 * t979, (t1009 * t861 + t1063 * t854) * t889 + (t1010 * t860 + t1064 * t855) * t886 + (t1011 * t859 + t1065 * t853) * t883 + (t880 * t954 + t881 * t953 + t882 * t952) * t1141, 0; t1043 + t1046 + t1049, 0, 0, t892 * t1049 + t893 * t1046 + t894 * t1043 + 0.2e1 * (t985 * t1085 + t986 * t1087 + t987 * t1089) * t911, 0.2e1 * t927 * t926 * t1020 + 0.2e1 * t929 * t928 * t1018 + 0.2e1 * t931 * t930 * t1016 + 0.2e1 * (t874 * t897 * t987 + t875 * t902 * t986 + t876 * t907 * t985) * t911, ((t894 * t908 + t906) * t867 * t1091 + (t893 * t903 + t901) * t866 * t1093 + (t892 * t898 + t896) * t865 * t1095) * t912, 0, 0, 0, 0, (-t873 * t1132 + t786 * t931) * t889 + (-t872 * t1135 + t785 * t929) * t886 + (-t871 * t1138 + t784 * t927) * t883 + ((t1054 * t843 + t998 * t1084) * t870 + (t1000 * t1086 + t1055 * t842) * t869 + (t1002 * t1088 + t1056 * t841) * t868) * t911, (t873 * t1133 + t789 * t931) * t889 + (t872 * t1136 + t788 * t929) * t886 + (t871 * t1139 + t787 * t927) * t883 + ((t1054 * t837 + t998) * t906 * t870 + (t1055 * t836 + t1000) * t901 * t869 + (t1056 * t835 + t1002) * t896 * t868) * t911, t1012 * t931 + t1013 * t929 + t1014 * t927 - t1039 * t873 - t1040 * t872 - t1041 * t871, (t1063 * t873 + t1066 * t931) * t889 + (t1064 * t872 + t1067 * t929) * t886 + (t1065 * t871 + t1068 * t927) * t883, 0;];
tau_reg  = t1;
