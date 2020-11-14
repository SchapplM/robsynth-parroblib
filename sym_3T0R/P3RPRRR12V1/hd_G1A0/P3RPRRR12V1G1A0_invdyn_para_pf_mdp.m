% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RPRRR12V1G1A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [14x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPRRR12V1G1A0_convert_par2_MPV_fixb.m

% Output:
% tauX [3x1]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:21
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRRR12V1G1A0_invdyn_para_pf_mdp(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(14,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_mdp: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_mdp: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_mdp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:21:38
% EndTime: 2020-08-06 18:21:42
% DurationCPUTime: 3.40s
% Computational Cost: add. (12057->340), mult. (18295->579), div. (2190->14), fcn. (14151->18), ass. (0->234)
t1104 = sin(qJ(3,3));
t1089 = 0.1e1 / t1104 ^ 2;
t1116 = xDP(3);
t1192 = t1116 ^ 2 / pkin(3) ^ 2;
t1247 = t1089 * t1192;
t1106 = sin(qJ(3,2));
t1091 = 0.1e1 / t1106 ^ 2;
t1246 = t1091 * t1192;
t1108 = sin(qJ(3,1));
t1093 = 0.1e1 / t1108 ^ 2;
t1245 = t1093 * t1192;
t1110 = cos(qJ(3,3));
t1244 = t1110 * t1247;
t1112 = cos(qJ(3,2));
t1243 = t1112 * t1246;
t1114 = cos(qJ(3,1));
t1242 = t1114 * t1245;
t1100 = legFrame(1,3);
t1073 = sin(t1100);
t1076 = cos(t1100);
t1034 = t1073 * g(1) - t1076 * g(2);
t1037 = t1076 * g(1) + t1073 * g(2);
t1109 = sin(qJ(1,1));
t1115 = cos(qJ(1,1));
t1241 = -t1034 * t1109 + t1037 * t1115;
t1099 = legFrame(2,3);
t1072 = sin(t1099);
t1075 = cos(t1099);
t1033 = t1072 * g(1) - t1075 * g(2);
t1036 = t1075 * g(1) + t1072 * g(2);
t1107 = sin(qJ(1,2));
t1113 = cos(qJ(1,2));
t1240 = -t1033 * t1107 + t1036 * t1113;
t1098 = legFrame(3,3);
t1071 = sin(t1098);
t1074 = cos(t1098);
t1032 = t1071 * g(1) - t1074 * g(2);
t1035 = t1074 * g(1) + t1071 * g(2);
t1105 = sin(qJ(1,3));
t1111 = cos(qJ(1,3));
t1239 = -t1032 * t1105 + t1035 * t1111;
t1238 = 2 * MDP(8);
t1237 = 2 * qJ(2,1);
t1236 = 2 * qJ(2,2);
t1235 = 2 * qJ(2,3);
t1067 = t1104 * pkin(3) + qJ(2,3);
t1056 = 0.1e1 / t1067;
t1068 = t1106 * pkin(3) + qJ(2,2);
t1059 = 0.1e1 / t1068;
t1069 = -t1108 * pkin(3) - qJ(2,1);
t1062 = 0.1e1 / t1069;
t1088 = 0.1e1 / t1104;
t1090 = 0.1e1 / t1106;
t1092 = 0.1e1 / t1108;
t1057 = 0.1e1 / t1067 ^ 2;
t1060 = 0.1e1 / t1068 ^ 2;
t1063 = 0.1e1 / t1069 ^ 2;
t1119 = pkin(1) + pkin(5);
t1198 = t1074 * t1111;
t1023 = -t1071 * t1105 + t1198;
t1201 = t1071 * t1111;
t1026 = t1074 * t1105 + t1201;
t1058 = t1056 * t1057;
t1102 = xDDP(2);
t1103 = xDDP(1);
t1187 = t1116 * t1110;
t1087 = (pkin(6) + t1119);
t1118 = xDP(1);
t1065 = (t1087 * t1118);
t1117 = xDP(2);
t1038 = -qJ(2,3) * t1117 + t1065;
t1066 = t1117 * t1087;
t1041 = qJ(2,3) * t1118 + t1066;
t1158 = t1105 * t1117 + t1111 * t1118;
t1159 = -t1105 * t1118 + t1111 * t1117;
t963 = ((t1038 * t1111 + t1105 * t1041) * t1074 + (-t1105 * t1038 + t1041 * t1111) * t1071) * t1104 + qJ(2,3) * t1187 + (t1104 * t1187 + (-t1158 * t1071 + t1159 * t1074) * (t1110 - 0.1e1) * (t1110 + 0.1e1)) * pkin(3);
t1216 = t1088 * t963;
t996 = t1159 * t1071 + t1158 * t1074;
t1219 = t1087 * t996;
t948 = (0.2e1 * t1057 * t1187 - t1058 * t963) * t996 * t1088 + (t1023 * t1103 + t1026 * t1102 - (t1216 - t1219) * t1057 * t996) * t1056;
t1234 = pkin(1) * t948;
t1197 = t1075 * t1113;
t1024 = -t1072 * t1107 + t1197;
t1200 = t1072 * t1113;
t1027 = t1075 * t1107 + t1200;
t1061 = t1059 * t1060;
t1186 = t1116 * t1112;
t1039 = -qJ(2,2) * t1117 + t1065;
t1042 = qJ(2,2) * t1118 + t1066;
t1156 = t1107 * t1117 + t1113 * t1118;
t1157 = -t1107 * t1118 + t1113 * t1117;
t964 = ((t1039 * t1113 + t1107 * t1042) * t1075 + (-t1107 * t1039 + t1042 * t1113) * t1072) * t1106 + qJ(2,2) * t1186 + (t1106 * t1186 + (-t1156 * t1072 + t1157 * t1075) * (t1112 - 0.1e1) * (t1112 + 0.1e1)) * pkin(3);
t1215 = t1090 * t964;
t997 = t1157 * t1072 + t1156 * t1075;
t1218 = t1087 * t997;
t949 = (0.2e1 * t1060 * t1186 - t1061 * t964) * t997 * t1090 + (t1024 * t1103 + t1027 * t1102 - (t1215 - t1218) * t1060 * t997) * t1059;
t1233 = pkin(1) * t949;
t1196 = t1076 * t1115;
t1025 = -t1073 * t1109 + t1196;
t1199 = t1073 * t1115;
t1028 = t1076 * t1109 + t1199;
t1064 = t1062 * t1063;
t1185 = t1116 * t1114;
t1040 = -qJ(2,1) * t1117 + t1065;
t1043 = qJ(2,1) * t1118 + t1066;
t1154 = t1109 * t1117 + t1115 * t1118;
t1155 = -t1109 * t1118 + t1115 * t1117;
t965 = ((t1040 * t1115 + t1109 * t1043) * t1076 + (-t1109 * t1040 + t1043 * t1115) * t1073) * t1108 + qJ(2,1) * t1185 + (t1108 * t1185 + (-t1154 * t1073 + t1155 * t1076) * (t1114 - 0.1e1) * (t1114 + 0.1e1)) * pkin(3);
t1214 = t1092 * t965;
t998 = t1155 * t1073 + t1154 * t1076;
t1217 = t1087 * t998;
t950 = (0.2e1 * t1063 * t1185 + t1064 * t965) * t998 * t1092 - (t1025 * t1103 + t1028 * t1102 - (t1214 - t1217) * t1063 * t998) * t1062;
t1232 = pkin(1) * t950;
t1094 = t1110 ^ 2;
t1231 = 0.2e1 * t1094 - 0.1e1;
t1095 = t1112 ^ 2;
t1230 = 0.2e1 * t1095 - 0.1e1;
t1096 = t1114 ^ 2;
t1229 = 0.2e1 * t1096 - 0.1e1;
t1204 = t1057 * t1088;
t1179 = t996 * t1204;
t960 = 0.2e1 * t963 * t1179;
t1228 = t948 * t1235 + t960;
t1203 = t1060 * t1090;
t1178 = t997 * t1203;
t961 = 0.2e1 * t964 * t1178;
t1227 = t949 * t1236 + t961;
t1202 = t1063 * t1092;
t1177 = t998 * t1202;
t962 = 0.2e1 * t965 * t1177;
t1226 = t950 * t1237 + t962;
t993 = t996 ^ 2;
t1225 = t1057 * t993;
t1224 = t1058 * t993;
t994 = t997 ^ 2;
t1223 = t1060 * t994;
t1222 = t1061 * t994;
t995 = t998 ^ 2;
t1221 = t1063 * t995;
t1220 = t1064 * t995;
t1101 = xDDP(3);
t1124 = 0.1e1 / pkin(3);
t1191 = t1101 * t1124;
t1014 = (-t1191 - t1244) * t1088;
t1213 = t1014 * t1104;
t1015 = (-t1191 - t1243) * t1090;
t1212 = t1015 * t1106;
t1016 = (-t1191 - t1242) * t1092;
t1211 = t1016 * t1108;
t1029 = t1035 * t1105;
t1030 = t1036 * t1107;
t1031 = t1037 * t1109;
t1195 = t1088 * t1110;
t1194 = t1090 * t1112;
t1193 = t1092 * t1114;
t1190 = t1104 * t1110;
t1189 = t1106 * t1112;
t1188 = t1108 * t1114;
t1184 = t1116 * t1124;
t1123 = pkin(3) ^ 2;
t1183 = -(t1087 ^ 2) - t1123;
t1182 = qJ(2,3) * t1225;
t1181 = qJ(2,2) * t1223;
t1180 = qJ(2,1) * t1221;
t1176 = t1087 * t1216;
t1175 = t1087 * t1215;
t1174 = t1087 * t1214;
t1173 = t1110 * t1225;
t1172 = t1112 * t1223;
t1171 = t1114 * t1221;
t1167 = t1119 * t1192;
t1166 = t996 * t1056 * t1184;
t1165 = t997 * t1059 * t1184;
t1164 = t998 * t1062 * t1184;
t1008 = t1014 * t1110 - t1088 * t1192;
t1009 = t1015 * t1112 - t1090 * t1192;
t1010 = t1016 * t1114 - t1092 * t1192;
t1163 = t1088 * t1166;
t1162 = t1090 * t1165;
t1161 = t1092 * t1164;
t1160 = pkin(3) * t1087 * t1184;
t1120 = qJ(2,3) ^ 2;
t1017 = t1105 * t1067 + t1087 * t1111;
t1020 = -t1067 * t1111 + t1087 * t1105;
t987 = t1017 * t1074 - t1071 * t1020;
t990 = t1071 * t1017 + t1020 * t1074;
t1153 = t1101 * t1195 + ((t1056 * t1110 * t1219 + t1088 * t1116) * t1104 + qJ(2,3) * t1088 * t1184) * t1056 * t1089 * t1116 + t987 * t1056 * t1103 + t990 * t1056 * t1102 - (-t1088 * t1160 * t1190 + (t1104 * t1176 + ((t1094 * t1123 - t1120 + t1183) * t1104 + (t1094 - 0.1e1) * t1235 * pkin(3)) * t996) * t1056) * t1179 - t996 * t1058 * t1176;
t1121 = qJ(2,2) ^ 2;
t1018 = t1107 * t1068 + t1087 * t1113;
t1021 = -t1068 * t1113 + t1087 * t1107;
t988 = t1018 * t1075 - t1072 * t1021;
t991 = t1072 * t1018 + t1021 * t1075;
t1152 = t1101 * t1194 + ((t1059 * t1112 * t1218 + t1090 * t1116) * t1106 + qJ(2,2) * t1090 * t1184) * t1059 * t1091 * t1116 + t988 * t1059 * t1103 + t991 * t1059 * t1102 - (-t1090 * t1160 * t1189 + (t1106 * t1175 + ((t1095 * t1123 - t1121 + t1183) * t1106 + (t1095 - 0.1e1) * t1236 * pkin(3)) * t997) * t1059) * t1178 - t997 * t1061 * t1175;
t1122 = qJ(2,1) ^ 2;
t1019 = -t1109 * t1069 + t1087 * t1115;
t1022 = t1069 * t1115 + t1087 * t1109;
t989 = t1019 * t1076 - t1073 * t1022;
t992 = t1073 * t1019 + t1022 * t1076;
t1151 = t1101 * t1193 - ((-t1062 * t1114 * t1217 + t1092 * t1116) * t1108 + qJ(2,1) * t1092 * t1184) * t1062 * t1093 * t1116 - t989 * t1062 * t1103 - t992 * t1062 * t1102 - (-t1092 * t1160 * t1188 - (t1108 * t1174 + ((t1096 * t1123 - t1122 + t1183) * t1108 + (t1096 - 0.1e1) * t1237 * pkin(3)) * t998) * t1062) * t1177 + t998 * t1064 * t1174;
t1146 = -g(1) * t1201 + g(2) * t1198 - t1029 + t1153;
t921 = t1146 - t1182 - t1234;
t969 = -t1213 + (-t1225 - t1247) * t1110;
t972 = -t1104 * t1225 + t1008;
t1150 = t972 * MDP(12) + t969 * MDP(13) + t948 * MDP(4) + t921 * MDP(6);
t1145 = -g(1) * t1200 + g(2) * t1197 - t1030 + t1152;
t922 = t1145 - t1181 - t1233;
t970 = -t1212 + (-t1223 - t1246) * t1112;
t973 = -t1106 * t1223 + t1009;
t1149 = t973 * MDP(12) + t970 * MDP(13) + t949 * MDP(4) + t922 * MDP(6);
t1144 = -g(1) * t1199 + g(2) * t1196 - t1031 + t1151;
t923 = t1144 - t1180 - t1232;
t971 = -t1211 + (-t1221 - t1245) * t1114;
t974 = -t1108 * t1221 + t1010;
t1148 = t974 * MDP(12) + t971 * MDP(13) + t950 * MDP(4) + t923 * MDP(6);
t1147 = t950 * t1193 + t949 * t1194 + t948 * t1195;
t1048 = g(1) * t1109 - g(2) * t1115;
t1049 = g(1) * t1115 + g(2) * t1109;
t1143 = t1048 * t1076 + t1049 * t1073 + t1119 * t950 - t1151 + t1180;
t1046 = -t1107 * g(1) + t1113 * g(2);
t1047 = t1113 * g(1) + t1107 * g(2);
t1142 = -t1046 * t1075 + t1047 * t1072 + t1119 * t949 - t1152 + t1181;
t1044 = -t1105 * g(1) + t1111 * g(2);
t1045 = t1111 * g(1) + t1105 * g(2);
t1141 = -t1044 * t1074 + t1045 * t1071 + t1119 * t948 - t1153 + t1182;
t933 = -t1044 * t1071 - t1045 * t1074 + t1089 * t1167 + t1228;
t975 = t1014 * t1119 + t1163 * t1235;
t1140 = t1110 * (t1110 * t948 + 0.2e1 * t1166) * MDP(7) + t948 * MDP(1) + (-t1213 - t1244) * MDP(10) + (t933 * t1104 - t1110 * t975) * MDP(12) + (t1104 * t975 + t933 * t1110) * MDP(13) + (t1032 * t1111 + t1029) * MDP(2) + t1239 * MDP(3) + (t1146 - 0.2e1 * t1234) * MDP(4) + (t1228 - t1239) * MDP(5) + (t948 * t1120 + qJ(2,3) * t960 + (-t1153 + t1234) * pkin(1) + t1035 * (t1105 * pkin(1) - t1111 * qJ(2,3)) + t1032 * (t1111 * pkin(1) + t1105 * qJ(2,3))) * MDP(6) + t1008 * MDP(9) + (t1231 * t1163 - t948 * t1190) * t1238;
t934 = -t1046 * t1072 - t1047 * t1075 + t1091 * t1167 + t1227;
t976 = t1015 * t1119 + t1162 * t1236;
t1139 = t1112 * (t1112 * t949 + 0.2e1 * t1165) * MDP(7) + t949 * MDP(1) + (-t1212 - t1243) * MDP(10) + (t934 * t1106 - t1112 * t976) * MDP(12) + (t1106 * t976 + t934 * t1112) * MDP(13) + (t1033 * t1113 + t1030) * MDP(2) + t1240 * MDP(3) + (t1145 - 0.2e1 * t1233) * MDP(4) + (t1227 - t1240) * MDP(5) + (t949 * t1121 + qJ(2,2) * t961 + (-t1152 + t1233) * pkin(1) + t1036 * (t1107 * pkin(1) - t1113 * qJ(2,2)) + t1033 * (t1113 * pkin(1) + t1107 * qJ(2,2))) * MDP(6) + t1009 * MDP(9) + (t1230 * t1162 - t949 * t1189) * t1238;
t935 = t1048 * t1073 - t1049 * t1076 + t1093 * t1167 + t1226;
t977 = -0.2e1 * qJ(2,1) * t1161 + t1016 * t1119;
t1138 = t1114 * (t1114 * t950 - 0.2e1 * t1164) * MDP(7) + t950 * MDP(1) + (-t1211 - t1242) * MDP(10) + (t935 * t1108 - t1114 * t977) * MDP(12) + (t1108 * t977 + t935 * t1114) * MDP(13) + (t1034 * t1115 + t1031) * MDP(2) + t1241 * MDP(3) + (t1144 - 0.2e1 * t1232) * MDP(4) + (t1226 - t1241) * MDP(5) + (t950 * t1122 + qJ(2,1) * t962 + (-t1151 + t1232) * pkin(1) + t1037 * (t1109 * pkin(1) - t1115 * qJ(2,1)) + t1034 * (t1115 * pkin(1) + t1109 * qJ(2,1))) * MDP(6) + t1010 * MDP(9) + (-t1229 * t1161 - t950 * t1188) * t1238;
t1 = [(t989 * t1220 - t988 * t1222 - t987 * t1224) * MDP(5) + (t1103 - g(1)) * MDP(14) - (t1138 * t1025 + t1148 * t989) * t1062 + (t1139 * t1024 + t1149 * t988) * t1059 + (t1140 * t1023 + t1150 * t987) * t1056; (t992 * t1220 - t991 * t1222 - t990 * t1224) * MDP(5) + (t1102 - g(2)) * MDP(14) - (t1138 * t1028 + t1148 * t992) * t1062 + (t1139 * t1027 + t1149 * t991) * t1059 + (t1140 * t1026 + t1150 * t990) * t1056; t1147 * MDP(4) + (-t1088 * t1173 - t1090 * t1172 - t1092 * t1171) * MDP(5) + (t923 * t1193 + t922 * t1194 + t921 * t1195) * MDP(6) + (t974 * t1193 + t973 * t1194 + t972 * t1195) * MDP(12) + (t971 * t1193 + t970 * t1194 + t969 * t1195) * MDP(13) + (t1101 - g(3)) * MDP(14) + ((-t1171 - t1172 - t1173) * MDP(7) + (-t995 * t1229 * t1202 - t994 * t1230 * t1203 - t993 * t1231 * t1204) * MDP(8) - t1147 * MDP(9) + (t948 + t949 + t950) * MDP(10) + (-t1014 * t1088 - t1015 * t1090 - t1016 * t1092) * MDP(11) + (-t1092 * (g(3) * t1108 - t1143 * t1114) - t1090 * (g(3) * t1106 - t1142 * t1112) - t1088 * (g(3) * t1104 - t1141 * t1110)) * MDP(12) + (-t1092 * (g(3) * t1114 + t1143 * t1108) - t1090 * (g(3) * t1112 + t1142 * t1106) - t1088 * (g(3) * t1110 + t1141 * t1104)) * MDP(13)) * t1124;];
tauX  = t1;
