% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RRR1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [2x3]
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [10x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRR1G1A0_convert_par2_MPV_fixb.m

% Output:
% tauX [3x1]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 15:38
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRR1G1A0_invdyn_para_pf_mdp(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1),zeros(10,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRR1G1A0_invdyn_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRR1G1A0_invdyn_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRR1G1A0_invdyn_para_pf_mdp: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RRR1G1A0_invdyn_para_pf_mdp: qJ has to be [2x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRR1G1A0_invdyn_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRR1G1A0_invdyn_para_pf_mdp: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRR1G1A0_invdyn_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRR1G1A0_invdyn_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'P3RRR1G1A0_invdyn_para_pf_mdp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 15:38:38
% EndTime: 2019-05-03 15:38:41
% DurationCPUTime: 3.89s
% Computational Cost: add. (25095->290), mult. (39520->529), div. (2880->9), fcn. (29778->26), ass. (0->220)
t1117 = 0.1e1 / pkin(1);
t1086 = qJ(1,1) + qJ(2,1);
t1069 = sin(t1086);
t1072 = cos(t1086);
t1098 = sin(qJ(1,1));
t1104 = cos(qJ(1,1));
t1037 = t1069 * t1104 - t1072 * t1098;
t1199 = 0.1e1 / t1037;
t1209 = t1199 * t1117;
t1085 = qJ(1,2) + qJ(2,2);
t1068 = sin(t1085);
t1071 = cos(t1085);
t1096 = sin(qJ(1,2));
t1102 = cos(qJ(1,2));
t1036 = t1068 * t1102 - t1071 * t1096;
t1200 = 0.1e1 / t1036;
t1208 = t1200 * t1117;
t1084 = qJ(1,3) + qJ(2,3);
t1067 = sin(t1084);
t1070 = cos(t1084);
t1094 = sin(qJ(1,3));
t1100 = cos(qJ(1,3));
t1035 = t1067 * t1100 - t1070 * t1094;
t1201 = 0.1e1 / t1035;
t1207 = t1201 * t1117;
t1108 = xP(3);
t1081 = sin(t1108);
t1082 = cos(t1108);
t1109 = koppelP(3,2);
t1112 = koppelP(3,1);
t1043 = -t1081 * t1109 + t1082 * t1112;
t1105 = xDP(3);
t1106 = xDP(2);
t1000 = t1043 * t1105 + t1106;
t1087 = legFrame(3,3);
t1073 = sin(t1087);
t1076 = cos(t1087);
t1040 = t1081 * t1112 + t1082 * t1109;
t1107 = xDP(1);
t997 = t1040 * t1105 - t1107;
t958 = t1070 * (t1000 * t1073 - t1076 * t997) + (t1000 * t1076 + t1073 * t997) * t1067;
t1110 = koppelP(2,2);
t1113 = koppelP(2,1);
t1044 = -t1081 * t1110 + t1082 * t1113;
t1001 = t1044 * t1105 + t1106;
t1088 = legFrame(2,3);
t1074 = sin(t1088);
t1077 = cos(t1088);
t1041 = t1081 * t1113 + t1082 * t1110;
t998 = t1041 * t1105 - t1107;
t959 = t1071 * (t1001 * t1074 - t1077 * t998) + (t1001 * t1077 + t1074 * t998) * t1068;
t1111 = koppelP(1,2);
t1114 = koppelP(1,1);
t1045 = -t1081 * t1111 + t1082 * t1114;
t1002 = t1045 * t1105 + t1106;
t1089 = legFrame(1,3);
t1075 = sin(t1089);
t1078 = cos(t1089);
t1042 = t1081 * t1114 + t1082 * t1111;
t999 = t1042 * t1105 - t1107;
t960 = t1072 * (t1002 * t1075 - t1078 * t999) + (t1002 * t1078 + t1075 * t999) * t1069;
t1052 = -t1112 * t1094 + t1100 * t1109;
t1055 = t1094 * t1109 + t1100 * t1112;
t1206 = (t1052 * t1082 + t1055 * t1081) * t1076 - t1073 * (-t1052 * t1081 + t1055 * t1082);
t1053 = -t1113 * t1096 + t1102 * t1110;
t1056 = t1096 * t1110 + t1102 * t1113;
t1205 = (t1053 * t1082 + t1056 * t1081) * t1077 - t1074 * (-t1053 * t1081 + t1056 * t1082);
t1054 = -t1114 * t1098 + t1104 * t1111;
t1057 = t1098 * t1111 + t1104 * t1114;
t1204 = (t1054 * t1082 + t1057 * t1081) * t1078 - t1075 * (-t1054 * t1081 + t1057 * t1082);
t1203 = -2 * pkin(1);
t1202 = 2 * pkin(1);
t1115 = pkin(2) ^ 2;
t1124 = t1067 * t1094 + t1070 * t1100;
t1168 = 2 * pkin(2);
t1188 = t1201 * t958;
t1116 = 1 / pkin(2);
t1152 = t1116 * t1117;
t946 = pkin(1) * ((-t1094 * t1106 - t1100 * t1107) * t1076 - t1073 * (-t1094 * t1107 + t1100 * t1106) + t1206 * t1105) - t958 * pkin(2);
t1189 = t1201 * t946;
t1139 = t1152 * t1189;
t952 = t958 * t1207;
t940 = t952 + t1139 / 0.2e1;
t943 = t952 + t1139;
t1198 = (t943 * t1115 + (t1124 * t1168 * t940 + t1188) * pkin(1)) * t958;
t1123 = t1068 * t1096 + t1071 * t1102;
t1186 = t1200 * t959;
t947 = pkin(1) * ((-t1096 * t1106 - t1102 * t1107) * t1077 - t1074 * (-t1096 * t1107 + t1102 * t1106) + t1205 * t1105) - t959 * pkin(2);
t1187 = t1200 * t947;
t1138 = t1152 * t1187;
t953 = t959 * t1208;
t941 = t953 + t1138 / 0.2e1;
t944 = t953 + t1138;
t1197 = (t944 * t1115 + (t1123 * t1168 * t941 + t1186) * pkin(1)) * t959;
t1122 = t1069 * t1098 + t1072 * t1104;
t1184 = t1199 * t960;
t948 = pkin(1) * ((-t1098 * t1106 - t1104 * t1107) * t1078 - t1075 * (-t1098 * t1107 + t1104 * t1106) + t1204 * t1105) - t960 * pkin(2);
t1185 = t1199 * t948;
t1137 = t1152 * t1185;
t954 = t960 * t1209;
t942 = t954 + t1137 / 0.2e1;
t945 = t954 + t1137;
t1196 = (t945 * t1115 + (t1122 * t1168 * t942 + t1184) * pkin(1)) * t960;
t1003 = t1067 * t1076 + t1070 * t1073;
t973 = pkin(1) * (t1073 * t1100 + t1076 * t1094) + t1003 * pkin(2);
t1083 = t1105 ^ 2;
t1090 = xDDP(3);
t1091 = xDDP(2);
t985 = -t1040 * t1083 + t1043 * t1090 + t1091;
t1195 = t973 * t985;
t1005 = t1068 * t1077 + t1071 * t1074;
t974 = pkin(1) * (t1074 * t1102 + t1077 * t1096) + t1005 * pkin(2);
t986 = -t1041 * t1083 + t1044 * t1090 + t1091;
t1194 = t974 * t986;
t1007 = t1069 * t1078 + t1072 * t1075;
t975 = pkin(1) * (t1075 * t1104 + t1078 * t1098) + t1007 * pkin(2);
t987 = -t1042 * t1083 + t1045 * t1090 + t1091;
t1193 = t975 * t987;
t1004 = -t1067 * t1073 + t1070 * t1076;
t976 = pkin(1) * (-t1073 * t1094 + t1076 * t1100) + t1004 * pkin(2);
t1092 = xDDP(1);
t988 = -t1040 * t1090 - t1043 * t1083 + t1092;
t1192 = t976 * t988;
t1006 = -t1068 * t1074 + t1071 * t1077;
t977 = pkin(1) * (-t1074 * t1096 + t1077 * t1102) + t1006 * pkin(2);
t989 = -t1041 * t1090 - t1044 * t1083 + t1092;
t1191 = t977 * t989;
t1008 = -t1069 * t1075 + t1072 * t1078;
t978 = pkin(1) * (-t1075 * t1098 + t1078 * t1104) + t1008 * pkin(2);
t990 = -t1042 * t1090 - t1045 * t1083 + t1092;
t1190 = t978 * t990;
t1133 = t1067 * (t1040 * t1073 + t1043 * t1076) - t1070 * (t1040 * t1076 - t1043 * t1073);
t1183 = t1201 * (-pkin(1) * t1206 + t1133 * pkin(2));
t1182 = t1201 * t1133;
t1181 = t1201 * t973;
t1180 = t1201 * t976;
t1132 = t1068 * (t1041 * t1074 + t1044 * t1077) - t1071 * (t1041 * t1077 - t1044 * t1074);
t1179 = t1200 * (-pkin(1) * t1205 + t1132 * pkin(2));
t1178 = t1200 * t1132;
t1177 = t1200 * t974;
t1176 = t1200 * t977;
t1131 = t1069 * (t1042 * t1075 + t1045 * t1078) - t1072 * (t1042 * t1078 - t1045 * t1075);
t1175 = t1199 * (-pkin(1) * t1204 + t1131 * pkin(2));
t1174 = t1199 * t1131;
t1173 = t1199 * t975;
t1172 = t1199 * t978;
t1167 = t1003 * t1201;
t1166 = t1004 * t1201;
t1165 = t1005 * t1200;
t1164 = t1006 * t1200;
t1163 = t1007 * t1199;
t1162 = t1008 * t1199;
t1027 = 0.1e1 / t1035 ^ 2;
t1118 = 1 / pkin(1) ^ 2;
t1160 = t1027 * t1118;
t1029 = 0.1e1 / t1036 ^ 2;
t1158 = t1029 * t1118;
t1031 = 0.1e1 / t1037 ^ 2;
t1156 = t1031 * t1118;
t1046 = -t1073 * g(1) + g(2) * t1076;
t1049 = g(1) * t1076 + g(2) * t1073;
t1151 = t1046 * t1067 + t1049 * t1070;
t1047 = -t1074 * g(1) + g(2) * t1077;
t1050 = g(1) * t1077 + g(2) * t1074;
t1150 = t1047 * t1068 + t1050 * t1071;
t1048 = -t1075 * g(1) + g(2) * t1078;
t1051 = g(1) * t1078 + g(2) * t1075;
t1149 = t1048 * t1069 + t1051 * t1072;
t1148 = -t1046 * t1070 + t1049 * t1067;
t1147 = -t1047 * t1071 + t1050 * t1068;
t1146 = -t1048 * t1072 + t1051 * t1069;
t1145 = t943 * (pkin(1) * t1124 + pkin(2)) * t946;
t1144 = t944 * (pkin(1) * t1123 + pkin(2)) * t947;
t1143 = t945 * (pkin(1) * t1122 + pkin(2)) * t948;
t1142 = t958 ^ 2 * t1160;
t1141 = t959 ^ 2 * t1158;
t1140 = t960 ^ 2 * t1156;
t1136 = t940 * t1139;
t1135 = t941 * t1138;
t1134 = t942 * t1137;
t931 = (-(-pkin(2) * t943 - t1124 * t1188) * t1027 * t958 + t943 * t1201 * t1189) * t1118 + (t1003 * t985 + t1004 * t988) * t1207;
t932 = (-(-pkin(2) * t944 - t1123 * t1186) * t1029 * t959 + t944 * t1200 * t1187) * t1118 + (t1005 * t986 + t1006 * t989) * t1208;
t933 = (-t1031 * (-t945 * pkin(2) - t1122 * t1184) * t960 + t1199 * t945 * t1185) * t1118 + (t1007 * t987 + t1008 * t990) * t1209;
t1103 = cos(qJ(2,1));
t1101 = cos(qJ(2,2));
t1099 = cos(qJ(2,3));
t1097 = sin(qJ(2,1));
t1095 = sin(qJ(2,2));
t1093 = sin(qJ(2,3));
t1080 = t1092 - g(1);
t1079 = t1091 - g(2);
t1039 = -t1081 * t1090 - t1082 * t1083;
t1038 = -t1081 * t1083 + t1082 * t1090;
t1013 = t1079 * t1081 + t1080 * t1082;
t1012 = t1079 * t1082 - t1080 * t1081;
t996 = t1048 * t1098 + t1051 * t1104;
t995 = -t1048 * t1104 + t1051 * t1098;
t994 = t1047 * t1096 + t1050 * t1102;
t993 = -t1047 * t1102 + t1050 * t1096;
t992 = t1046 * t1094 + t1049 * t1100;
t991 = -t1046 * t1100 + t1049 * t1094;
t930 = pkin(1) * (-t933 * t1097 + t1103 * t1140) + t1149;
t929 = pkin(1) * (t1097 * t1140 + t933 * t1103) + t1146;
t928 = pkin(1) * (-t932 * t1095 + t1101 * t1141) + t1150;
t927 = pkin(1) * (t1095 * t1141 + t932 * t1101) + t1147;
t926 = pkin(1) * (-t931 * t1093 + t1099 * t1142) + t1151;
t925 = pkin(1) * (t1093 * t1142 + t931 * t1099) + t1148;
t924 = (-(t1190 + t1193) * t1209 + (-t1143 - t1196) * t1156) * t1116 + t933;
t923 = (-(t1191 + t1194) * t1208 + (-t1144 - t1197) * t1158) * t1116 + t932;
t922 = (-(t1192 + t1195) * t1207 + (-t1145 - t1198) * t1160) * t1116 + t931;
t921 = (-(t1190 / 0.2e1 + t1193 / 0.2e1) * t1209 + (-t1196 / 0.2e1 - t1143 / 0.2e1) * t1156) * t1116 + t933;
t920 = (-(t1191 / 0.2e1 + t1194 / 0.2e1) * t1208 + (-t1197 / 0.2e1 - t1144 / 0.2e1) * t1158) * t1116 + t932;
t919 = (-(t1192 / 0.2e1 + t1195 / 0.2e1) * t1207 + (-t1198 / 0.2e1 - t1145 / 0.2e1) * t1160) * t1116 + t931;
t918 = (t1097 * t921 + t1103 * t1134) * t1203 + t1149;
t917 = (-t1097 * t1134 + t1103 * t921) * t1202 + t1146;
t916 = (t1095 * t920 + t1101 * t1135) * t1203 + t1150;
t915 = (-t1095 * t1135 + t1101 * t920) * t1202 + t1147;
t914 = (t1093 * t919 + t1099 * t1136) * t1203 + t1151;
t913 = (-t1093 * t1136 + t1099 * t919) * t1202 + t1148;
t1 = [t1039 * MDP(8) - t1038 * MDP(9) + (-t1012 * t1081 + t1013 * t1082) * MDP(10) + ((t1162 * t933 + t1164 * t932 + t1166 * t931) * MDP(1) + (t1162 * t995 + t1164 * t993 + t1166 * t991) * MDP(2) + (t1162 * t996 + t1164 * t994 + t1166 * t992) * MDP(3) + (t924 * t1162 + t923 * t1164 + t922 * t1166) * MDP(4) + (t917 * t1162 + t915 * t1164 + t913 * t1166) * MDP(5) + (t918 * t1162 + t916 * t1164 + t914 * t1166) * MDP(6) + ((-t1172 * t924 - t1176 * t923 - t1180 * t922) * MDP(4) + (-t1172 * t929 - t1176 * t927 - t1180 * t925) * MDP(5) + (-t1172 * t930 - t1176 * t928 - t1180 * t926) * MDP(6)) * t1116) * t1117; t1038 * MDP(8) + t1039 * MDP(9) + (t1012 * t1082 + t1013 * t1081) * MDP(10) + ((t1163 * t933 + t1165 * t932 + t1167 * t931) * MDP(1) + (t1163 * t995 + t1165 * t993 + t1167 * t991) * MDP(2) + (t1163 * t996 + t1165 * t994 + t1167 * t992) * MDP(3) + (t924 * t1163 + t923 * t1165 + t922 * t1167) * MDP(4) + (t917 * t1163 + t915 * t1165 + t913 * t1167) * MDP(5) + (t918 * t1163 + t916 * t1165 + t914 * t1167) * MDP(6) + ((-t1173 * t924 - t1177 * t923 - t1181 * t922) * MDP(4) + (-t1173 * t929 - t1177 * t927 - t1181 * t925) * MDP(5) + (-t1173 * t930 - t1177 * t928 - t1181 * t926) * MDP(6)) * t1116) * t1117; t1090 * MDP(7) + t1012 * MDP(8) - t1013 * MDP(9) + ((t1174 * t933 + t1178 * t932 + t1182 * t931) * MDP(1) + (t1174 * t995 + t1178 * t993 + t1182 * t991) * MDP(2) + (t1174 * t996 + t1178 * t994 + t1182 * t992) * MDP(3) + (t924 * t1174 + t923 * t1178 + t922 * t1182) * MDP(4) + (t917 * t1174 + t915 * t1178 + t913 * t1182) * MDP(5) + (t918 * t1174 + t916 * t1178 + t914 * t1182) * MDP(6) + ((-t1175 * t924 - t1179 * t923 - t1183 * t922) * MDP(4) + (-t1175 * t929 - t1179 * t927 - t1183 * t925) * MDP(5) + (-t1175 * t930 - t1179 * t928 - t1183 * t926) * MDP(6)) * t1116) * t1117;];
tauX  = t1;
