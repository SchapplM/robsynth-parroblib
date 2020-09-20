% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PRRRR1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [12x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRRR1G1A0_convert_par2_MPV_fixb.m

% Output:
% tauX [3x1]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 20:34
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRRR1G1A0_invdyn_para_pf_mdp(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_mdp: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3PRRRR1G1A0_invdyn_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 20:34:14
% EndTime: 2020-03-09 20:34:16
% DurationCPUTime: 2.27s
% Computational Cost: add. (2639->274), mult. (6866->536), div. (2670->17), fcn. (7164->18), ass. (0->230)
t1028 = sin(qJ(2,3));
t1001 = 0.1e1 / t1028 ^ 2;
t1033 = cos(qJ(3,3));
t1039 = xDP(2);
t1040 = xDP(1);
t1027 = sin(qJ(3,3));
t1034 = cos(qJ(2,3));
t1111 = t1027 * t1034;
t1021 = legFrame(3,3);
t994 = sin(t1021);
t997 = cos(t1021);
t982 = -t1039 * t997 + t1040 * t994;
t964 = (-t1039 * t994 - t1040 * t997) * t1033 - t982 * t1111;
t961 = t964 ^ 2;
t1162 = t1001 * t961;
t1030 = sin(qJ(2,2));
t1004 = 0.1e1 / t1030 ^ 2;
t1035 = cos(qJ(3,2));
t1029 = sin(qJ(3,2));
t1036 = cos(qJ(2,2));
t1110 = t1029 * t1036;
t1022 = legFrame(2,3);
t995 = sin(t1022);
t998 = cos(t1022);
t983 = -t1039 * t998 + t1040 * t995;
t965 = (-t1039 * t995 - t1040 * t998) * t1035 - t983 * t1110;
t962 = t965 ^ 2;
t1161 = t1004 * t962;
t1032 = sin(qJ(2,1));
t1007 = 0.1e1 / t1032 ^ 2;
t1037 = cos(qJ(3,1));
t1031 = sin(qJ(3,1));
t1038 = cos(qJ(2,1));
t1109 = t1031 * t1038;
t1023 = legFrame(1,3);
t996 = sin(t1023);
t999 = cos(t1023);
t984 = -t1039 * t999 + t1040 * t996;
t966 = (-t1039 * t996 - t1040 * t999) * t1037 - t984 * t1109;
t963 = t966 ^ 2;
t1160 = t1007 * t963;
t1041 = 0.1e1 / pkin(2);
t1000 = 0.1e1 / t1028;
t1003 = 0.1e1 / t1030;
t1006 = 0.1e1 / t1032;
t1009 = 0.1e1 / t1033;
t1013 = 0.1e1 / t1035;
t1017 = 0.1e1 / t1037;
t1042 = 0.1e1 / pkin(2) ^ 2;
t1052 = t1033 ^ 2;
t1010 = 0.1e1 / t1052;
t1055 = t1035 ^ 2;
t1014 = 0.1e1 / t1055;
t1058 = t1037 ^ 2;
t1018 = 0.1e1 / t1058;
t1159 = xDDP(3) - g(3);
t1158 = MDP(2) * t1041;
t1157 = MDP(6) * t1041;
t1156 = MDP(7) * t1041;
t1155 = MDP(8) * t1041;
t1154 = MDP(9) * t1041;
t1025 = xDDP(2);
t1026 = xDDP(1);
t1011 = t1009 * t1010;
t1114 = t1011 * t1041;
t1151 = t1000 * t1162;
t1108 = t1033 * t1034;
t970 = t1027 * t994 + t997 * t1108;
t976 = -t1027 * t997 + t994 * t1108;
t979 = t982 ^ 2;
t934 = t1114 * t1151 + (t979 * t1114 + (t1025 * t976 + t1026 * t970) * t1009) * t1000 + t1159;
t991 = g(1) * t997 + g(2) * t994;
t1080 = t991 * t1028 + t1034 * t934;
t1153 = t1000 * t1080;
t1152 = t1000 * t982;
t1015 = t1013 * t1014;
t1113 = t1015 * t1041;
t1148 = t1003 * t1161;
t1107 = t1035 * t1036;
t972 = t1029 * t995 + t998 * t1107;
t977 = -t1029 * t998 + t995 * t1107;
t980 = t983 ^ 2;
t935 = t1113 * t1148 + (t980 * t1113 + (t1025 * t977 + t1026 * t972) * t1013) * t1003 + t1159;
t992 = g(1) * t998 + g(2) * t995;
t1078 = t992 * t1030 + t1036 * t935;
t1150 = t1003 * t1078;
t1149 = t1003 * t983;
t1019 = t1017 * t1018;
t1112 = t1019 * t1041;
t1145 = t1006 * t1160;
t1106 = t1037 * t1038;
t974 = t1031 * t996 + t999 * t1106;
t978 = -t1031 * t999 + t996 * t1106;
t981 = t984 ^ 2;
t936 = t1112 * t1145 + (t981 * t1112 + (t1025 * t978 + t1026 * t974) * t1017) * t1006 + t1159;
t993 = g(1) * t999 + g(2) * t996;
t1076 = t993 * t1032 + t1038 * t936;
t1147 = t1006 * t1076;
t1146 = t1006 * t984;
t1144 = t1009 * t994;
t1143 = t1009 * t997;
t1012 = 0.1e1 / t1052 ^ 2;
t1142 = t1012 * t961;
t1141 = t1013 * t995;
t1140 = t1013 * t998;
t1016 = 0.1e1 / t1055 ^ 2;
t1139 = t1016 * t962;
t1138 = t1017 * t996;
t1137 = t1017 * t999;
t1020 = 0.1e1 / t1058 ^ 2;
t1136 = t1020 * t963;
t1119 = t1000 * t1010;
t967 = -t1033 * t997 - t994 * t1111;
t971 = -t1033 * t994 + t997 * t1111;
t910 = ((t1025 * t971 + t1026 * t967) * t1119 + (-(t1000 * t1034 * t964 - t1027 * t1028 * t982) * t1001 * t964 - (-t1027 * t964 + t1034 * t982) * t1152) * t1009 * t1114) * t1041;
t1135 = t1027 * t910;
t1126 = t1042 * t979;
t1084 = t1027 * t1126;
t958 = t1011 * t1084 + (-t1025 * t997 + t1026 * t994) * t1041 * t1009;
t1134 = t1027 * t958;
t1117 = t1003 * t1014;
t968 = -t1035 * t998 - t995 * t1110;
t973 = -t1035 * t995 + t998 * t1110;
t911 = ((t1025 * t973 + t1026 * t968) * t1117 + (-(t1003 * t1036 * t965 - t1029 * t1030 * t983) * t1004 * t965 - (-t1029 * t965 + t1036 * t983) * t1149) * t1013 * t1113) * t1041;
t1133 = t1029 * t911;
t1125 = t1042 * t980;
t1083 = t1029 * t1125;
t959 = t1015 * t1083 + (-t1025 * t998 + t1026 * t995) * t1041 * t1013;
t1132 = t1029 * t959;
t1115 = t1006 * t1018;
t969 = -t1037 * t999 - t996 * t1109;
t975 = -t1037 * t996 + t999 * t1109;
t912 = ((t1025 * t975 + t1026 * t969) * t1115 + (-(t1006 * t1038 * t966 - t1031 * t1032 * t984) * t1007 * t966 - (-t1031 * t966 + t1038 * t984) * t1146) * t1017 * t1112) * t1041;
t1131 = t1031 * t912;
t1124 = t1042 * t981;
t1082 = t1031 * t1124;
t960 = t1019 * t1082 + (-t1025 * t999 + t1026 * t996) * t1041 * t1017;
t1130 = t1031 * t960;
t1129 = t1033 * t958;
t1128 = t1035 * t959;
t1127 = t1037 * t960;
t1123 = t910 * t1028;
t1122 = t911 * t1030;
t1121 = t912 * t1032;
t1120 = t1000 * t1009;
t1118 = t1003 * t1013;
t1116 = t1006 * t1017;
t1072 = t1042 * t964 * t1152;
t1100 = t1001 * t1142;
t955 = (t1010 * t979 + t1100) * t1042;
t883 = t910 * t1108 + (-t1033 * t955 - t1134) * t1028 - 0.2e1 * t1011 * t1072 * t1111;
t1071 = t1042 * t965 * t1149;
t1094 = t1004 * t1139;
t956 = (t1014 * t980 + t1094) * t1042;
t884 = t911 * t1107 + (-t1035 * t956 - t1132) * t1030 - 0.2e1 * t1015 * t1071 * t1110;
t1070 = t1042 * t966 * t1146;
t1088 = t1007 * t1136;
t957 = (t1018 * t981 + t1088) * t1042;
t885 = t912 * t1106 + (-t1037 * t957 - t1130) * t1032 - 0.2e1 * t1019 * t1070 * t1109;
t1105 = t970 * t1120;
t1104 = t976 * t1120;
t1103 = t967 * t1119;
t1102 = t971 * t1119;
t1101 = t1000 * t1142;
t1099 = t972 * t1118;
t1098 = t977 * t1118;
t1097 = t968 * t1117;
t1096 = t973 * t1117;
t1095 = t1003 * t1139;
t1093 = t974 * t1116;
t1092 = t978 * t1116;
t1091 = t969 * t1115;
t1090 = t975 * t1115;
t1089 = t1006 * t1136;
t1087 = t1009 * t1135;
t1086 = t1013 * t1133;
t1085 = t1017 * t1131;
t896 = 0.2e1 * t1010 * t1072 + t1135;
t886 = -t1034 * t896 + (t1027 * t955 - t1129) * t1028;
t1081 = -t1028 * t934 + t1034 * t991;
t897 = 0.2e1 * t1014 * t1071 + t1133;
t887 = -t1036 * t897 + (t1029 * t956 - t1128) * t1030;
t1079 = -t1030 * t935 + t1036 * t992;
t895 = 0.2e1 * t1018 * t1070 + t1131;
t888 = -t1038 * t895 + (t1031 * t957 - t1127) * t1032;
t1077 = -t1032 * t936 + t1038 * t993;
t1075 = 0.2e1 * (t1033 * t1135 + (0.2e1 * t1009 - t1011) * t1072) * t1119;
t1074 = 0.2e1 * (t1035 * t1133 + (0.2e1 * t1013 - t1015) * t1071) * t1117;
t1073 = 0.2e1 * (t1037 * t1131 + (0.2e1 * t1017 - t1019) * t1070) * t1115;
t1069 = t1027 * t1080 * t1119;
t1068 = t1027 * t1100;
t1067 = t1034 * t1100;
t1066 = t1029 * t1078 * t1117;
t1065 = t1029 * t1094;
t1064 = t1036 * t1094;
t1063 = t1031 * t1076 * t1115;
t1062 = t1031 * t1088;
t1061 = t1038 * t1088;
t1043 = t1041 * t1042;
t990 = g(1) * t996 - g(2) * t999;
t989 = g(1) * t995 - g(2) * t998;
t988 = g(1) * t994 - g(2) * t997;
t954 = t1017 * t1124 + t1130;
t953 = t1013 * t1125 + t1132;
t952 = t1009 * t1126 + t1134;
t948 = -t1018 * t1082 + t1127;
t947 = -t1014 * t1083 + t1128;
t946 = -t1010 * t1084 + t1129;
t945 = (-0.2e1 * t1018 + t1020) * t1042 * t1160;
t944 = (-0.2e1 * t1014 + t1016) * t1042 * t1161;
t943 = (-0.2e1 * t1010 + t1012) * t1042 * t1162;
t918 = t1031 * t990 + t1037 * t1077;
t917 = t1031 * t1077 - t1037 * t990;
t916 = t1029 * t989 + t1035 * t1079;
t915 = t1029 * t1079 - t1035 * t989;
t914 = t1027 * t988 + t1033 * t1081;
t913 = t1027 * t1081 - t1033 * t988;
t909 = t1038 * t912;
t908 = t1036 * t911;
t907 = t1034 * t910;
t903 = -t1042 * t1089 + t909;
t902 = -t1042 * t1095 + t908;
t901 = -t1042 * t1101 + t907;
t900 = -t1042 * t1061 - t1121;
t899 = -t1042 * t1064 - t1122;
t898 = -t1042 * t1067 - t1123;
t894 = t895 * t1031;
t893 = t897 * t1029;
t892 = t896 * t1027;
t1 = [(t936 * t1093 + t935 * t1099 + t934 * t1105) * MDP(1) + (t912 * t1091 + t911 * t1097 + t910 * t1103) * t1158 + (t901 * t1105 + t902 * t1099 + t903 * t1093 + (t1076 * t1091 + t1078 * t1097 + t1080 * t1103) * t1041) * MDP(3) + (t898 * t1105 + t899 * t1099 + t900 * t1093 + (t1077 * t1091 + t1079 * t1097 + t1081 * t1103) * t1041) * MDP(4) + ((-t996 * t1062 - t995 * t1065 - t994 * t1068) * t1043 + (t894 * t1091 + t893 * t1097 + t892 * t1103) * t1041) * MDP(5) + (t969 * t1073 + t968 * t1074 + t967 * t1075 + t945 * t1138 + t944 * t1141 + t943 * t1144) * t1157 + (t996 * t1085 + t995 * t1086 + t994 * t1087 + t954 * t1091 + t953 * t1097 + t952 * t1103) * t1156 + (t948 * t1091 + t947 * t1097 + t946 * t1103 + t910 * t994 + t911 * t995 + t912 * t996) * t1155 + (t960 * t1138 + t959 * t1141 + t958 * t1144) * t1154 + (t883 * t1105 + t884 * t1099 + t885 * t1093 + ((t969 * t1147 + t917 * t996) * t1017 + (t968 * t1150 + t915 * t995) * t1013 + (t967 * t1153 + t913 * t994) * t1009) * t1041) * MDP(10) + (t886 * t1105 + t887 * t1099 + t888 * t1093 + (-t969 * t1063 - t968 * t1066 - t967 * t1069 + t918 * t1138 + t916 * t1141 + t914 * t1144) * t1041) * MDP(11) + (t1026 - g(1)) * MDP(12); (t936 * t1092 + t935 * t1098 + t934 * t1104) * MDP(1) + (t912 * t1090 + t911 * t1096 + t910 * t1102) * t1158 + (t901 * t1104 + t902 * t1098 + t903 * t1092 + (t1076 * t1090 + t1078 * t1096 + t1080 * t1102) * t1041) * MDP(3) + (t898 * t1104 + t899 * t1098 + t900 * t1092 + (t1077 * t1090 + t1079 * t1096 + t1081 * t1102) * t1041) * MDP(4) + ((t999 * t1062 + t998 * t1065 + t997 * t1068) * t1043 + (t894 * t1090 + t893 * t1096 + t892 * t1102) * t1041) * MDP(5) + (t975 * t1073 + t973 * t1074 + t971 * t1075 - t945 * t1137 - t944 * t1140 - t943 * t1143) * t1157 + (-t999 * t1085 - t998 * t1086 - t997 * t1087 + t954 * t1090 + t953 * t1096 + t952 * t1102) * t1156 + (t948 * t1090 + t947 * t1096 + t946 * t1102 - t910 * t997 - t911 * t998 - t912 * t999) * t1155 + (-t960 * t1137 - t959 * t1140 - t958 * t1143) * t1154 + (t883 * t1104 + t884 * t1098 + t885 * t1092 + ((t975 * t1147 - t917 * t999) * t1017 + (t973 * t1150 - t915 * t998) * t1013 + (t971 * t1153 - t913 * t997) * t1009) * t1041) * MDP(10) + (t886 * t1104 + t887 * t1098 + t888 * t1092 + (-t975 * t1063 - t973 * t1066 - t971 * t1069 - t918 * t1137 - t916 * t1140 - t914 * t1143) * t1041) * MDP(11) + (t1025 - g(2)) * MDP(12); (t907 + t908 + t909) * MDP(3) + (-t1121 - t1122 - t1123) * MDP(4) + (t885 + t884 + t883) * MDP(10) + (t886 + t887 + t888) * MDP(11) + t1159 * MDP(12) + ((-t1089 - t1095 - t1101) * MDP(3) + (-t1061 - t1064 - t1067) * MDP(4)) * t1042 + (0.3e1 * t1159 + (t1093 + t1099 + t1105) * t1026 + (t1092 + t1098 + t1104) * t1025 + ((t1006 * t981 + t1145) * t1019 + (t1003 * t980 + t1148) * t1015 + (t1000 * t979 + t1151) * t1011) * t1041) * MDP(1);];
tauX  = t1;
