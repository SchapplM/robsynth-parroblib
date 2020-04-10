% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PRRRR1G2P3A0
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
%   see P3PRRRR1G2P3A0_convert_par2_MPV_fixb.m

% Output:
% tauX [3x1]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:16
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRRR1G2P3A0_invdyn_para_pf_mdp(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G2P3A0_invdyn_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR1G2P3A0_invdyn_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR1G2P3A0_invdyn_para_pf_mdp: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G2P3A0_invdyn_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR1G2P3A0_invdyn_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G2P3A0_invdyn_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G2P3A0_invdyn_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G2P3A0_invdyn_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3PRRRR1G2P3A0_invdyn_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:16:38
% EndTime: 2020-03-09 21:16:40
% DurationCPUTime: 2.85s
% Computational Cost: add. (2571->270), mult. (6998->546), div. (3429->17), fcn. (7395->18), ass. (0->234)
t1101 = legFrame(3,2);
t1074 = sin(t1101);
t1108 = sin(qJ(2,3));
t1113 = cos(qJ(3,3));
t1199 = t1108 * t1113;
t1077 = cos(t1101);
t1107 = sin(qJ(3,3));
t1214 = t1077 * t1107;
t1047 = t1074 * t1199 - t1214;
t1220 = t1074 * t1107;
t1050 = t1077 * t1199 + t1220;
t1120 = xDP(1);
t1245 = xDP(2);
t1056 = -t1074 * t1245 + t1077 * t1120;
t1053 = t1056 ^ 2;
t1081 = 0.1e1 / t1108;
t1089 = 0.1e1 / t1113;
t1104 = xDDP(3);
t1105 = xDDP(2);
t1106 = xDDP(1);
t1114 = cos(qJ(2,3));
t1132 = t1113 ^ 2;
t1090 = 0.1e1 / t1132;
t1091 = t1089 * t1090;
t1121 = 0.1e1 / pkin(2);
t1203 = t1091 * t1121;
t1119 = xDP(3);
t1200 = t1107 * t1114;
t1044 = t1056 * t1200 - t1113 * t1119;
t1041 = t1044 ^ 2;
t1082 = 0.1e1 / t1108 ^ 2;
t1251 = t1041 * t1082;
t1023 = -t1074 * g(1) - t1077 * g(2) + (t1114 * t1104 + (t1047 * t1106 + t1050 * t1105) * t1089 + (t1053 + t1251) * t1203) * t1081;
t1014 = g(3) * t1108 + t1023 * t1114;
t1209 = t1081 * t1089;
t1254 = t1014 * t1209;
t1102 = legFrame(2,2);
t1075 = sin(t1102);
t1110 = sin(qJ(2,2));
t1115 = cos(qJ(3,2));
t1197 = t1110 * t1115;
t1078 = cos(t1102);
t1109 = sin(qJ(3,2));
t1212 = t1078 * t1109;
t1048 = t1075 * t1197 - t1212;
t1218 = t1075 * t1109;
t1051 = t1078 * t1197 + t1218;
t1058 = -t1075 * t1245 + t1078 * t1120;
t1054 = t1058 ^ 2;
t1084 = 0.1e1 / t1110;
t1093 = 0.1e1 / t1115;
t1116 = cos(qJ(2,2));
t1135 = t1115 ^ 2;
t1094 = 0.1e1 / t1135;
t1095 = t1093 * t1094;
t1202 = t1095 * t1121;
t1198 = t1109 * t1116;
t1045 = t1058 * t1198 - t1115 * t1119;
t1042 = t1045 ^ 2;
t1085 = 0.1e1 / t1110 ^ 2;
t1250 = t1042 * t1085;
t1024 = -t1075 * g(1) - t1078 * g(2) + (t1116 * t1104 + (t1048 * t1106 + t1051 * t1105) * t1093 + (t1054 + t1250) * t1202) * t1084;
t1015 = g(3) * t1110 + t1024 * t1116;
t1207 = t1084 * t1093;
t1253 = t1015 * t1207;
t1103 = legFrame(1,2);
t1076 = sin(t1103);
t1112 = sin(qJ(2,1));
t1117 = cos(qJ(3,1));
t1195 = t1112 * t1117;
t1079 = cos(t1103);
t1111 = sin(qJ(3,1));
t1210 = t1079 * t1111;
t1049 = t1076 * t1195 - t1210;
t1216 = t1076 * t1111;
t1052 = t1079 * t1195 + t1216;
t1060 = -t1076 * t1245 + t1079 * t1120;
t1055 = t1060 ^ 2;
t1087 = 0.1e1 / t1112;
t1097 = 0.1e1 / t1117;
t1118 = cos(qJ(2,1));
t1138 = t1117 ^ 2;
t1098 = 0.1e1 / t1138;
t1099 = t1097 * t1098;
t1201 = t1099 * t1121;
t1196 = t1111 * t1118;
t1046 = t1060 * t1196 - t1117 * t1119;
t1043 = t1046 ^ 2;
t1088 = 0.1e1 / t1112 ^ 2;
t1249 = t1043 * t1088;
t1025 = -t1076 * g(1) - t1079 * g(2) + (t1118 * t1104 + (t1049 * t1106 + t1052 * t1105) * t1097 + (t1055 + t1249) * t1201) * t1087;
t1016 = g(3) * t1112 + t1025 * t1118;
t1205 = t1087 * t1097;
t1252 = t1016 * t1205;
t1122 = 0.1e1 / pkin(2) ^ 2;
t1248 = -0.2e1 * t1090;
t1247 = -0.2e1 * t1094;
t1246 = -0.2e1 * t1098;
t1244 = MDP(2) * t1121;
t1243 = MDP(6) * t1121;
t1242 = MDP(7) * t1121;
t1241 = MDP(8) * t1121;
t1240 = MDP(9) * t1121;
t1158 = t1074 * t1105 - t1077 * t1106;
t1178 = t1081 * t1200;
t1167 = t1090 * t1178;
t1208 = t1081 * t1114;
t1224 = t1056 * t1081;
t1005 = (-t1158 * t1167 + (-t1081 * t1104 + (-(t1056 * t1107 * t1108 + t1044 * t1208) * t1082 * t1044 + (-t1044 * t1107 - t1056 * t1114) * t1224) * t1203) * t1089) * t1121;
t1239 = t1005 * t1113;
t1157 = t1075 * t1105 - t1078 * t1106;
t1176 = t1084 * t1198;
t1166 = t1094 * t1176;
t1206 = t1084 * t1116;
t1223 = t1058 * t1084;
t1006 = (-t1157 * t1166 + (-t1084 * t1104 + (-(t1058 * t1109 * t1110 + t1045 * t1206) * t1085 * t1045 + (-t1045 * t1109 - t1058 * t1116) * t1223) * t1202) * t1093) * t1121;
t1238 = t1006 * t1115;
t1156 = t1076 * t1105 - t1079 * t1106;
t1174 = t1087 * t1196;
t1165 = t1098 * t1174;
t1204 = t1087 * t1118;
t1222 = t1060 * t1087;
t1007 = (-t1156 * t1165 + (-t1087 * t1104 + (-(t1060 * t1111 * t1112 + t1046 * t1204) * t1088 * t1046 + (-t1046 * t1111 - t1060 * t1118) * t1222) * t1201) * t1097) * t1121;
t1237 = t1007 * t1117;
t1227 = t1053 * t1122;
t1182 = t1107 * t1227;
t1038 = t1089 * t1121 * t1158 + t1091 * t1182;
t1236 = t1038 * t1107;
t1235 = t1038 * t1113;
t1226 = t1054 * t1122;
t1181 = t1109 * t1226;
t1039 = t1093 * t1121 * t1157 + t1095 * t1181;
t1234 = t1039 * t1109;
t1233 = t1039 * t1115;
t1225 = t1055 * t1122;
t1180 = t1111 * t1225;
t1040 = t1097 * t1121 * t1156 + t1099 * t1180;
t1232 = t1040 * t1111;
t1231 = t1040 * t1117;
t1230 = t1041 * t1122;
t1229 = t1042 * t1122;
t1228 = t1043 * t1122;
t1221 = t1074 * t1089;
t1219 = t1075 * t1093;
t1217 = t1076 * t1097;
t1215 = t1077 * t1089;
t1213 = t1078 * t1093;
t1211 = t1079 * t1097;
t1092 = 0.1e1 / t1132 ^ 2;
t1194 = t1092 * t1251;
t1193 = t1092 * t1230;
t1096 = 0.1e1 / t1135 ^ 2;
t1192 = t1096 * t1250;
t1191 = t1096 * t1229;
t1100 = 0.1e1 / t1138 ^ 2;
t1190 = t1100 * t1249;
t1189 = t1100 * t1228;
t1188 = t1047 * t1209;
t1187 = t1048 * t1207;
t1186 = t1049 * t1205;
t1185 = t1050 * t1209;
t1184 = t1051 * t1207;
t1183 = t1052 * t1205;
t1179 = t1090 * t1208;
t1177 = t1094 * t1206;
t1175 = t1098 * t1204;
t1173 = t1107 * t1194;
t1172 = t1109 * t1192;
t1171 = t1111 * t1190;
t1170 = t1044 * t1122 * t1224;
t1169 = t1045 * t1122 * t1223;
t1168 = t1046 * t1122 * t1222;
t1164 = t1074 * t1167;
t1163 = t1075 * t1166;
t1162 = t1076 * t1165;
t1161 = t1077 * t1167;
t1160 = t1078 * t1166;
t1159 = t1079 * t1165;
t1155 = t1170 * t1248;
t1154 = t1169 * t1247;
t1153 = t1168 * t1246;
t1017 = g(3) * t1114 - t1023 * t1108;
t1062 = g(1) * t1077 - g(2) * t1074;
t1152 = -t1014 * t1178 + t1017 * t1107 + t1062 * t1113;
t1019 = g(3) * t1116 - t1024 * t1110;
t1063 = g(1) * t1078 - g(2) * t1075;
t1151 = -t1015 * t1176 + t1019 * t1109 + t1063 * t1115;
t1021 = g(3) * t1118 - t1025 * t1112;
t1064 = g(1) * t1079 - g(2) * t1076;
t1150 = -t1016 * t1174 + t1021 * t1111 + t1064 * t1117;
t1026 = -t1090 * t1182 + t1235;
t1149 = t1026 * t1167 - t1005;
t1027 = -t1094 * t1181 + t1233;
t1148 = t1027 * t1166 - t1006;
t1028 = -t1098 * t1180 + t1231;
t1147 = t1028 * t1165 - t1007;
t1032 = t1089 * t1227 + t1236;
t1146 = -t1005 * t1089 + t1032 * t1179;
t1033 = t1093 * t1226 + t1234;
t1145 = -t1006 * t1093 + t1033 * t1177;
t1034 = t1097 * t1225 + t1232;
t1144 = -t1007 * t1097 + t1034 * t1175;
t1080 = t1107 ^ 2;
t1143 = -t1014 * t1080 * t1179 - (t1017 * t1113 - t1062 * t1107) * t1089;
t1083 = t1109 ^ 2;
t1142 = -t1015 * t1083 * t1177 - (t1019 * t1115 - t1063 * t1109) * t1093;
t1086 = t1111 ^ 2;
t1141 = -t1016 * t1086 * t1175 - (t1021 * t1117 - t1064 * t1111) * t1097;
t1123 = t1121 * t1122;
t1037 = (t1055 * t1098 + t1190) * t1122;
t1036 = (t1054 * t1094 + t1192) * t1122;
t1035 = (t1053 * t1090 + t1194) * t1122;
t1031 = (t1246 + t1100) * t1088 * t1228;
t1030 = (t1247 + t1096) * t1085 * t1229;
t1029 = (t1248 + t1092) * t1082 * t1230;
t1004 = t1007 * t1118 - t1087 * t1189;
t1003 = t1006 * t1116 - t1084 * t1191;
t1002 = t1005 * t1114 - t1081 * t1193;
t1001 = -t1088 * t1118 * t1189 - t1007 * t1112;
t1000 = -t1085 * t1116 * t1191 - t1006 * t1110;
t999 = -t1082 * t1114 * t1193 - t1005 * t1108;
t998 = t1007 * t1086 + t1111 * t1153;
t997 = t1006 * t1083 + t1109 * t1154;
t996 = t1005 * t1080 + t1107 * t1155;
t995 = t1111 * t1237 + (-0.2e1 * t1097 + t1099) * t1168;
t994 = t1109 * t1238 + (-0.2e1 * t1093 + t1095) * t1169;
t993 = t1107 * t1239 + (-0.2e1 * t1089 + t1091) * t1170;
t992 = (t1037 * t1111 - t1231) * t1112 - t1118 * (t1007 * t1111 + t1153);
t991 = (t1036 * t1109 - t1233) * t1110 - t1116 * (t1006 * t1109 + t1154);
t990 = (t1035 * t1107 - t1235) * t1108 - t1114 * (t1005 * t1107 + t1155);
t989 = (-t1037 * t1117 - t1232) * t1112 + (0.2e1 * t1099 * t1111 * t1168 + t1237) * t1118;
t988 = (-t1036 * t1115 - t1234) * t1110 + (0.2e1 * t1095 * t1109 * t1169 + t1238) * t1116;
t987 = (-t1035 * t1113 - t1236) * t1108 + (0.2e1 * t1091 * t1107 * t1170 + t1239) * t1114;
t1 = [(t1023 * t1188 + t1024 * t1187 + t1025 * t1186) * MDP(1) + (t1005 * t1161 + t1006 * t1160 + t1007 * t1159) * t1244 + (t1002 * t1188 + t1003 * t1187 + t1004 * t1186 + (t1014 * t1161 + t1015 * t1160 + t1016 * t1159) * t1121) * MDP(3) + (t1000 * t1187 + t1001 * t1186 + t999 * t1188 + (t1017 * t1161 + t1019 * t1160 + t1021 * t1159) * t1121) * MDP(4) + ((t1077 * t1173 + t1078 * t1172 + t1079 * t1171) * t1123 + (t1159 * t998 + t1160 * t997 + t1161 * t996) * t1121) * MDP(5) + (-t1029 * t1215 - t1030 * t1213 - t1031 * t1211 + 0.2e1 * t1159 * t995 + 0.2e1 * t1160 * t994 + 0.2e1 * t1161 * t993) * t1243 + (t1144 * t1210 + t1145 * t1212 + t1146 * t1214) * t1242 + (t1077 * t1149 + t1078 * t1148 + t1079 * t1147) * t1241 + (-t1038 * t1215 - t1039 * t1213 - t1040 * t1211) * t1240 + (t987 * t1188 + t988 * t1187 + t989 * t1186 + (-t1150 * t1211 - t1151 * t1213 - t1152 * t1215) * t1121) * MDP(10) + (t990 * t1188 + t991 * t1187 + t992 * t1186 + (t1077 * t1143 + t1078 * t1142 + t1079 * t1141) * t1121) * MDP(11) + (t1106 - g(1)) * MDP(12); (t1023 * t1185 + t1024 * t1184 + t1025 * t1183) * MDP(1) + (-t1005 * t1164 - t1006 * t1163 - t1007 * t1162) * t1244 + (t1002 * t1185 + t1003 * t1184 + t1004 * t1183 + (-t1014 * t1164 - t1015 * t1163 - t1016 * t1162) * t1121) * MDP(3) + (t1000 * t1184 + t1001 * t1183 + t999 * t1185 + (-t1017 * t1164 - t1019 * t1163 - t1021 * t1162) * t1121) * MDP(4) + ((-t1074 * t1173 - t1075 * t1172 - t1076 * t1171) * t1123 + (-t1162 * t998 - t1163 * t997 - t1164 * t996) * t1121) * MDP(5) + (t1029 * t1221 + t1030 * t1219 + t1031 * t1217 - 0.2e1 * t1162 * t995 - 0.2e1 * t1163 * t994 - 0.2e1 * t1164 * t993) * t1243 + (-t1144 * t1216 - t1145 * t1218 - t1146 * t1220) * t1242 + (-t1074 * t1149 - t1075 * t1148 - t1076 * t1147) * t1241 + (t1038 * t1221 + t1039 * t1219 + t1040 * t1217) * t1240 + (t987 * t1185 + t988 * t1184 + t989 * t1183 + (t1150 * t1217 + t1151 * t1219 + t1152 * t1221) * t1121) * MDP(10) + (t990 * t1185 + t991 * t1184 + t992 * t1183 + (-t1074 * t1143 - t1075 * t1142 - t1076 * t1141) * t1121) * MDP(11) + (t1105 - g(2)) * MDP(12); (t1023 * t1208 + t1024 * t1206 + t1025 * t1204) * MDP(1) + (t1002 * t1208 + t1003 * t1206 + t1004 * t1204) * MDP(3) + (t1000 * t1206 + t1001 * t1204 + t999 * t1208) * MDP(4) + (t989 * t1204 + t988 * t1206 + t987 * t1208) * MDP(10) + (t992 * t1204 + t991 * t1206 + t990 * t1208) * MDP(11) + (t1104 - g(3)) * MDP(12) + ((-t1005 * t1209 - t1006 * t1207 - t1007 * t1205) * MDP(2) + (-t1252 - t1253 - t1254) * MDP(3) + (-t1017 * t1209 - t1019 * t1207 - t1021 * t1205) * MDP(4) + (-t1205 * t998 - t1207 * t997 - t1209 * t996) * MDP(5) + (-t1032 * t1209 - t1033 * t1207 - t1034 * t1205) * MDP(7) + (-t1026 * t1209 - t1027 * t1207 - t1028 * t1205) * MDP(8) + (-t1014 * t1081 - t1015 * t1084 - t1016 * t1087) * MDP(10) + (t1107 * t1254 + t1109 * t1253 + t1111 * t1252) * MDP(11) + 0.2e1 * (-t1205 * t995 - t1207 * t994 - t1209 * t993) * MDP(6)) * t1121;];
tauX  = t1;
