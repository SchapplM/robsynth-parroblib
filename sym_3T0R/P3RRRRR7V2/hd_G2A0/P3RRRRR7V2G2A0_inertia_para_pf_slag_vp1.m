% Calculate inertia matrix for parallel robot
% P3RRRRR7V2G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
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
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 10:08
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRRRR7V2G2A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V2G2A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V2G2A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRRRR7V2G2A0_inertia_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR7V2G2A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRRRR7V2G2A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRRRR7V2G2A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V2G2A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V2G2A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 09:42:44
% EndTime: 2020-08-07 09:42:47
% DurationCPUTime: 3.93s
% Computational Cost: add. (9966->462), mult. (16596->660), div. (792->14), fcn. (8658->84), ass. (0->324)
t1329 = rSges(2,3) + pkin(5);
t1328 = 2 * pkin(1);
t1327 = m(2) / 0.2e1;
t1326 = m(3) / 0.2e1;
t1325 = Icges(2,2) / 0.2e1;
t1324 = Icges(3,2) / 0.2e1;
t1173 = pkin(2) ^ 2;
t1323 = t1173 / 0.2e1;
t1126 = qJ(2,1) + qJ(3,1);
t1082 = sin(t1126);
t1297 = qJ(3,1) + qJ(1,1);
t1090 = qJ(2,1) + t1297;
t1298 = -qJ(3,1) + qJ(1,1);
t1091 = -qJ(2,1) + t1298;
t1165 = 2 * qJ(2,1);
t1125 = t1165 + qJ(3,1);
t1139 = sin(qJ(2,1));
t1149 = cos(qJ(1,1));
t1164 = 2 * qJ(3,1);
t1292 = qJ(1,1) - (2 * qJ(2,1));
t1138 = sin(qJ(3,1));
t1303 = t1138 * pkin(1);
t1314 = pkin(6) + pkin(5);
t1111 = (pkin(7) + t1314);
t1319 = 2 * t1111;
t1206 = ((cos(t1091) + cos(t1090)) * t1328 + (sin(t1091) + sin(t1090)) * t1319 + (cos((2 * qJ(3,1)) - t1292) + cos(qJ(1,1) + t1165 + t1164) + 0.2e1 * t1149) * pkin(3) + (cos(qJ(3,1) - t1292) + cos(qJ(1,1) + t1125) + cos(t1298) + cos(t1297)) * pkin(2)) / (-t1173 * sin((qJ(2,1) - qJ(3,1))) + pkin(2) * (0.2e1 * t1303 + pkin(2) * t1082 + (sin((t1164 + qJ(2,1))) - t1139) * pkin(3))) / 0.2e1;
t1124 = qJ(2,2) + qJ(3,2);
t1081 = sin(t1124);
t1295 = qJ(3,2) + qJ(1,2);
t1088 = qJ(2,2) + t1295;
t1296 = -qJ(3,2) + qJ(1,2);
t1089 = -qJ(2,2) + t1296;
t1162 = 2 * qJ(2,2);
t1123 = qJ(3,2) + t1162;
t1136 = sin(qJ(2,2));
t1146 = cos(qJ(1,2));
t1161 = 2 * qJ(3,2);
t1291 = qJ(1,2) - (2 * qJ(2,2));
t1135 = sin(qJ(3,2));
t1305 = t1135 * pkin(1);
t1207 = ((cos(t1089) + cos(t1088)) * t1328 + (sin(t1089) + sin(t1088)) * t1319 + (cos((2 * qJ(3,2)) - t1291) + cos(qJ(1,2) + t1162 + t1161) + 0.2e1 * t1146) * pkin(3) + (cos(qJ(3,2) - t1291) + cos(qJ(1,2) + t1123) + cos(t1296) + cos(t1295)) * pkin(2)) / (-t1173 * sin((qJ(2,2) - qJ(3,2))) + pkin(2) * (0.2e1 * t1305 + pkin(2) * t1081 + (sin((t1161 + qJ(2,2))) - t1136) * pkin(3))) / 0.2e1;
t1122 = qJ(2,3) + qJ(3,3);
t1080 = sin(t1122);
t1293 = qJ(3,3) + qJ(1,3);
t1086 = qJ(2,3) + t1293;
t1294 = -qJ(3,3) + qJ(1,3);
t1087 = -qJ(2,3) + t1294;
t1159 = 2 * qJ(2,3);
t1121 = t1159 + qJ(3,3);
t1133 = sin(qJ(2,3));
t1143 = cos(qJ(1,3));
t1158 = 2 * qJ(3,3);
t1290 = qJ(1,3) - (2 * qJ(2,3));
t1132 = sin(qJ(3,3));
t1307 = t1132 * pkin(1);
t1208 = ((cos(t1087) + cos(t1086)) * t1328 + (sin(t1087) + sin(t1086)) * t1319 + (cos((2 * qJ(3,3)) - t1290) + cos(qJ(1,3) + t1159 + t1158) + 0.2e1 * t1143) * pkin(3) + (cos(qJ(3,3) - t1290) + cos(qJ(1,3) + t1121) + cos(t1294) + cos(t1293)) * pkin(2)) / (-t1173 * sin((qJ(2,3) - qJ(3,3))) + pkin(2) * (0.2e1 * t1307 + pkin(2) * t1080 + (sin((t1158 + qJ(2,3))) - t1133) * pkin(3))) / 0.2e1;
t1322 = -2 * pkin(1);
t1320 = 0.2e1 * pkin(3);
t1142 = cos(qJ(2,3));
t1116 = t1142 ^ 2;
t1318 = 0.2e1 * t1116;
t1145 = cos(qJ(2,2));
t1118 = t1145 ^ 2;
t1317 = 0.2e1 * t1118;
t1148 = cos(qJ(2,1));
t1120 = t1148 ^ 2;
t1316 = 0.2e1 * t1120;
t1315 = pkin(2) * m(3);
t1313 = m(2) * rSges(2,1);
t1312 = m(2) * rSges(2,2);
t1311 = m(3) * rSges(3,2);
t1168 = rSges(2,2) ^ 2;
t1170 = rSges(2,1) ^ 2;
t1310 = m(3) * t1323 + (-t1168 + t1170) * t1327 + t1325 - Icges(2,1) / 0.2e1;
t1167 = rSges(3,2) ^ 2;
t1169 = rSges(3,1) ^ 2;
t1309 = (-t1167 + t1169) * t1326 - Icges(3,1) / 0.2e1 + t1324;
t1109 = rSges(3,3) + t1314;
t1308 = m(3) * t1109;
t1141 = cos(qJ(3,3));
t1115 = t1141 ^ 2;
t1101 = pkin(3) * t1115;
t1144 = cos(qJ(3,2));
t1117 = t1144 ^ 2;
t1102 = pkin(3) * t1117;
t1147 = cos(qJ(3,1));
t1119 = t1147 ^ 2;
t1103 = pkin(3) * t1119;
t1306 = t1132 * pkin(3);
t1304 = t1135 * pkin(3);
t1302 = t1138 * pkin(3);
t1301 = t1141 * pkin(2);
t1098 = t1141 * pkin(3);
t1300 = t1144 * pkin(2);
t1099 = t1144 * pkin(3);
t1299 = t1147 * pkin(2);
t1100 = t1147 * pkin(3);
t1232 = t1167 + t1169;
t1061 = t1232 * m(3) + Icges(3,3);
t1134 = sin(qJ(1,3));
t1186 = pkin(1) * t1134 - t1143 * t1111;
t1239 = t1132 * t1133;
t1015 = t1186 * t1239 + (t1115 - 0.1e1) * t1134 * pkin(3);
t1129 = legFrame(3,2);
t1092 = sin(t1129);
t1285 = t1015 * t1092;
t1095 = cos(t1129);
t1284 = t1015 * t1095;
t1137 = sin(qJ(1,2));
t1185 = pkin(1) * t1137 - t1146 * t1111;
t1238 = t1135 * t1136;
t1016 = t1185 * t1238 + (t1117 - 0.1e1) * t1137 * pkin(3);
t1130 = legFrame(2,2);
t1093 = sin(t1130);
t1283 = t1016 * t1093;
t1096 = cos(t1130);
t1282 = t1016 * t1096;
t1140 = sin(qJ(1,1));
t1184 = pkin(1) * t1140 - t1149 * t1111;
t1237 = t1138 * t1139;
t1017 = t1184 * t1237 + (t1119 - 0.1e1) * t1140 * pkin(3);
t1131 = legFrame(1,2);
t1094 = sin(t1131);
t1281 = t1017 * t1094;
t1097 = cos(t1131);
t1280 = t1017 * t1097;
t1059 = rSges(3,2) * t1308 - Icges(3,6);
t1060 = rSges(3,1) * t1308 - Icges(3,5);
t1083 = cos(t1122);
t1024 = -t1059 * t1083 - t1060 * t1080;
t1172 = 0.1e1 / pkin(3);
t1279 = t1024 * t1172;
t1084 = cos(t1124);
t1025 = -t1059 * t1084 - t1060 * t1081;
t1278 = t1025 * t1172;
t1085 = cos(t1126);
t1026 = -t1059 * t1085 - t1060 * t1082;
t1277 = t1026 * t1172;
t1235 = rSges(3,1) * t1315;
t1066 = t1141 * t1235;
t1234 = pkin(2) * t1311;
t1190 = t1132 * t1234;
t1033 = t1066 - t1190 + t1061;
t1276 = t1033 * t1172;
t1067 = t1144 * t1235;
t1189 = t1135 * t1234;
t1034 = t1067 - t1189 + t1061;
t1275 = t1034 * t1172;
t1068 = t1147 * t1235;
t1188 = t1138 * t1234;
t1035 = t1068 - t1188 + t1061;
t1274 = t1035 * t1172;
t1224 = pkin(3) * t1239;
t1070 = t1098 + pkin(2);
t1255 = t1070 * t1142;
t1273 = 0.1e1 / (pkin(1) - t1224 + t1255) / t1132;
t1223 = pkin(3) * t1238;
t1071 = t1099 + pkin(2);
t1252 = t1071 * t1145;
t1272 = 0.1e1 / (pkin(1) - t1223 + t1252) / t1135;
t1222 = pkin(3) * t1237;
t1072 = t1100 + pkin(2);
t1249 = t1072 * t1148;
t1271 = 0.1e1 / (pkin(1) - t1222 + t1249) / t1138;
t1044 = 0.1e1 / (t1142 * pkin(2) + pkin(3) * t1083 + pkin(1));
t1270 = t1044 * t1134;
t1269 = t1044 * t1143;
t1045 = 0.1e1 / (t1145 * pkin(2) + pkin(3) * t1084 + pkin(1));
t1268 = t1045 * t1137;
t1267 = t1045 * t1146;
t1046 = 0.1e1 / (t1148 * pkin(2) + pkin(3) * t1085 + pkin(1));
t1266 = t1046 * t1140;
t1265 = t1046 * t1149;
t1054 = pkin(1) * t1133 - t1306;
t1264 = t1054 * t1141;
t1055 = pkin(1) * t1136 - t1304;
t1263 = t1055 * t1144;
t1056 = pkin(1) * t1139 - t1302;
t1262 = t1056 * t1147;
t1261 = t1061 * t1172;
t1157 = pkin(2) / 0.2e1;
t1260 = (t1098 + t1157) * t1132;
t1259 = (t1099 + t1157) * t1135;
t1258 = (t1100 + t1157) * t1138;
t1257 = t1070 * t1092;
t1256 = t1070 * t1095;
t1254 = t1071 * t1093;
t1253 = t1071 * t1096;
t1251 = t1072 * t1094;
t1250 = t1072 * t1097;
t1248 = t1092 * t1134;
t1247 = t1093 * t1137;
t1246 = t1094 * t1140;
t1245 = t1095 * t1134;
t1244 = t1096 * t1137;
t1243 = t1097 * t1140;
t1171 = pkin(3) ^ 2;
t1242 = t1115 * t1171;
t1241 = t1117 * t1171;
t1240 = t1119 * t1171;
t1174 = 0.1e1 / pkin(2);
t1236 = t1172 * t1174;
t1233 = -t1171 / 0.2e1 + t1323;
t1231 = t1168 + t1170;
t1230 = pkin(2) * t1098;
t1229 = pkin(2) * t1099;
t1228 = pkin(2) * t1100;
t1227 = t1070 * t1306;
t1226 = t1071 * t1304;
t1225 = t1072 * t1302;
t1107 = -t1171 + t1173;
t1021 = pkin(1) * t1306 + (t1107 + 0.2e1 * t1230 + 0.2e1 * t1242) * t1133;
t1193 = t1134 * t1239;
t1030 = t1193 * t1320 - t1186;
t1040 = t1230 + t1233 + t1242;
t991 = (t1040 * t1248 - t1095 * t1227) * t1318 + (-t1021 * t1095 - t1030 * t1257) * t1142 - pkin(3) * t1285 - t1054 * t1256;
t1221 = t991 * t1273;
t992 = (-t1040 * t1245 - t1092 * t1227) * t1318 + (-t1092 * t1021 + t1030 * t1256) * t1142 + pkin(3) * t1284 - t1054 * t1257;
t1220 = t992 * t1273;
t1027 = t1307 + (-pkin(3) + t1301 + 0.2e1 * t1101) * t1133;
t1156 = -pkin(3) / 0.2e1;
t1048 = t1101 + t1301 / 0.2e1 + t1156;
t1179 = pkin(2) * t1193 + t1030 * t1141;
t997 = (-t1048 * t1248 + t1095 * t1260) * t1318 + (t1095 * t1027 + t1179 * t1092) * t1142 + t1285 + t1095 * t1264;
t1219 = t997 * t1273;
t998 = (t1048 * t1245 + t1092 * t1260) * t1318 + (t1092 * t1027 - t1179 * t1095) * t1142 - t1284 + t1092 * t1264;
t1218 = t998 * t1273;
t1022 = pkin(1) * t1304 + (t1107 + 0.2e1 * t1229 + 0.2e1 * t1241) * t1136;
t1192 = t1137 * t1238;
t1031 = t1192 * t1320 - t1185;
t1041 = t1229 + t1233 + t1241;
t993 = (t1041 * t1247 - t1096 * t1226) * t1317 + (-t1022 * t1096 - t1031 * t1254) * t1145 - pkin(3) * t1283 - t1055 * t1253;
t1217 = t993 * t1272;
t994 = (-t1041 * t1244 - t1093 * t1226) * t1317 + (-t1093 * t1022 + t1031 * t1253) * t1145 + pkin(3) * t1282 - t1055 * t1254;
t1216 = t994 * t1272;
t1028 = t1305 + (-pkin(3) + t1300 + 0.2e1 * t1102) * t1136;
t1049 = t1102 + t1300 / 0.2e1 + t1156;
t1178 = pkin(2) * t1192 + t1031 * t1144;
t999 = (-t1049 * t1247 + t1096 * t1259) * t1317 + (t1096 * t1028 + t1178 * t1093) * t1145 + t1283 + t1096 * t1263;
t1215 = t999 * t1272;
t1023 = pkin(1) * t1302 + (t1107 + 0.2e1 * t1228 + 0.2e1 * t1240) * t1139;
t1191 = t1140 * t1237;
t1032 = t1191 * t1320 - t1184;
t1042 = t1228 + t1233 + t1240;
t995 = (t1042 * t1246 - t1097 * t1225) * t1316 + (-t1023 * t1097 - t1032 * t1251) * t1148 - pkin(3) * t1281 - t1056 * t1250;
t1214 = t995 * t1271;
t996 = (-t1042 * t1243 - t1094 * t1225) * t1316 + (-t1094 * t1023 + t1032 * t1250) * t1148 + pkin(3) * t1280 - t1056 * t1251;
t1213 = t996 * t1271;
t1212 = t1329 * t1312 - Icges(2,6);
t1000 = (t1049 * t1244 + t1093 * t1259) * t1317 + (t1093 * t1028 - t1178 * t1096) * t1145 - t1282 + t1093 * t1263;
t1211 = t1000 * t1272;
t1029 = t1303 + (-pkin(3) + t1299 + 0.2e1 * t1103) * t1139;
t1050 = t1103 + t1299 / 0.2e1 + t1156;
t1177 = pkin(2) * t1191 + t1032 * t1147;
t1001 = (-t1050 * t1246 + t1097 * t1258) * t1316 + (t1097 * t1029 + t1177 * t1094) * t1148 + t1281 + t1097 * t1262;
t1210 = t1001 * t1271;
t1002 = (t1050 * t1243 + t1094 * t1258) * t1316 + (t1094 * t1029 - t1177 * t1097) * t1148 - t1280 + t1094 * t1262;
t1209 = t1002 * t1271;
t1006 = -0.2e1 * t1040 * t1143 * t1116 - ((pkin(1) - 0.2e1 * t1224) * t1143 + t1134 * t1111) * t1255 + pkin(3) * ((pkin(1) * t1239 - pkin(3) + t1101) * t1143 + t1111 * t1193);
t1205 = t1006 * t1273;
t1007 = -0.2e1 * t1041 * t1146 * t1118 - ((pkin(1) - 0.2e1 * t1223) * t1146 + t1137 * t1111) * t1252 + pkin(3) * ((pkin(1) * t1238 - pkin(3) + t1102) * t1146 + t1111 * t1192);
t1204 = t1007 * t1272;
t1008 = -0.2e1 * t1042 * t1149 * t1120 - ((pkin(1) - 0.2e1 * t1222) * t1149 + t1140 * t1111) * t1249 + pkin(3) * ((pkin(1) * t1237 - pkin(3) + t1103) * t1149 + t1111 * t1191);
t1203 = t1008 * t1271;
t1202 = t1174 * t1273;
t1201 = t1174 * t1272;
t1200 = t1174 * t1271;
t1199 = t1092 * t1269;
t1198 = t1095 * t1269;
t1197 = t1093 * t1267;
t1196 = t1096 * t1267;
t1195 = t1094 * t1265;
t1194 = t1097 * t1265;
t1187 = t1173 + t1232;
t1183 = t1231 * m(2) + t1187 * m(3) + Icges(2,3) + Icges(3,3);
t1182 = t1006 * t1172 * t1202;
t1181 = t1007 * t1172 * t1201;
t1180 = t1008 * t1172 * t1200;
t1176 = -pkin(2) * t1308 - t1329 * t1313 + Icges(2,5);
t1155 = 2 * pkin(1) ^ 2;
t1175 = Icges(1,3) + ((2 * pkin(5) ^ 2) + t1155 + ((4 * pkin(5) + 2 * rSges(2,3)) * rSges(2,3)) + t1231) * t1327 + (0.2e1 * t1109 ^ 2 + t1155 + t1187) * t1326 + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + t1324 + t1325 + Icges(3,1) / 0.2e1 + Icges(2,1) / 0.2e1;
t1079 = -rSges(2,1) * t1312 + Icges(2,4);
t1078 = -rSges(3,1) * t1311 + Icges(3,4);
t1077 = 0.2e1 * t1126;
t1076 = 0.2e1 * t1124;
t1075 = 0.2e1 * t1122;
t1069 = t1313 + t1315;
t1020 = 0.2e1 * t1068 + t1183 - 0.2e1 * t1188;
t1019 = 0.2e1 * t1067 + t1183 - 0.2e1 * t1189;
t1018 = 0.2e1 * t1066 + t1183 - 0.2e1 * t1190;
t1011 = (t1059 * t1138 - t1060 * t1147 + t1176) * t1139 - t1148 * (t1059 * t1147 + t1060 * t1138 + t1212);
t1010 = (t1059 * t1135 - t1060 * t1144 + t1176) * t1136 - t1145 * (t1059 * t1144 + t1060 * t1135 + t1212);
t1009 = (t1059 * t1132 - t1060 * t1141 + t1176) * t1133 - t1142 * (t1059 * t1141 + t1060 * t1132 + t1212);
t990 = cos(t1077) * t1309 + cos(t1165) * t1310 + t1068 + (t1069 * t1148 - t1139 * t1312) * t1328 + t1078 * sin(t1077) + t1079 * sin(t1165) + ((cos(t1125) * pkin(2) + t1085 * t1328) * rSges(3,1) + (t1082 * t1322 + (-sin(t1125) - t1138) * pkin(2)) * rSges(3,2)) * m(3) + t1175;
t989 = cos(t1076) * t1309 + cos(t1162) * t1310 + t1067 + (t1069 * t1145 - t1136 * t1312) * t1328 + t1078 * sin(t1076) + t1079 * sin(t1162) + ((cos(t1123) * pkin(2) + t1084 * t1328) * rSges(3,1) + (t1081 * t1322 + (-sin(t1123) - t1135) * pkin(2)) * rSges(3,2)) * m(3) + t1175;
t988 = cos(t1075) * t1309 + cos(t1159) * t1310 + t1066 + (t1069 * t1142 - t1133 * t1312) * t1328 + t1078 * sin(t1075) + t1079 * sin(t1159) + ((cos(t1121) * pkin(2) + t1083 * t1328) * rSges(3,1) + (t1080 * t1322 + (-sin(t1121) - t1132) * pkin(2)) * rSges(3,2)) * m(3) + t1175;
t987 = -t1026 * t1266 + t1035 * t1206 + t1061 * t1180;
t986 = -t1025 * t1268 + t1034 * t1207 + t1061 * t1181;
t985 = -t1024 * t1270 + t1033 * t1208 + t1061 * t1182;
t984 = -t1011 * t1266 + t1020 * t1206 + t1035 * t1180;
t983 = -t1010 * t1268 + t1019 * t1207 + t1034 * t1181;
t982 = -t1009 * t1270 + t1018 * t1208 + t1033 * t1182;
t981 = t1026 * t1194 + (t1002 * t1035 + t996 * t1261) * t1200;
t980 = t1025 * t1196 + (t1000 * t1034 + t994 * t1261) * t1201;
t979 = t1024 * t1198 + (t1033 * t998 + t992 * t1261) * t1202;
t978 = -t1026 * t1195 + (t1001 * t1035 + t995 * t1261) * t1200;
t977 = -t1025 * t1197 + (t1034 * t999 + t993 * t1261) * t1201;
t976 = -t1024 * t1199 + (t1033 * t997 + t991 * t1261) * t1202;
t975 = t1011 * t1206 + t1026 * t1180 - t990 * t1266;
t974 = t1010 * t1207 + t1025 * t1181 - t989 * t1268;
t973 = t1009 * t1208 + t1024 * t1182 - t988 * t1270;
t972 = t1011 * t1194 + (t1002 * t1020 + t996 * t1274) * t1200;
t971 = t1010 * t1196 + (t1000 * t1019 + t994 * t1275) * t1201;
t970 = t1009 * t1198 + (t1018 * t998 + t992 * t1276) * t1202;
t969 = -t1011 * t1195 + (t1001 * t1020 + t995 * t1274) * t1200;
t968 = -t1010 * t1197 + (t1019 * t999 + t993 * t1275) * t1201;
t967 = -t1009 * t1199 + (t1018 * t997 + t991 * t1276) * t1202;
t966 = t990 * t1194 + (t1002 * t1011 + t996 * t1277) * t1200;
t965 = t989 * t1196 + (t1000 * t1010 + t994 * t1278) * t1201;
t964 = t988 * t1198 + (t1009 * t998 + t992 * t1279) * t1202;
t963 = -t990 * t1195 + (t1001 * t1011 + t995 * t1277) * t1200;
t962 = -t989 * t1197 + (t1010 * t999 + t993 * t1278) * t1201;
t961 = -t988 * t1199 + (t1009 * t997 + t991 * t1279) * t1202;
t1 = [t964 * t1198 + t965 * t1196 + t966 * t1194 + m(4) + (t971 * t1211 + t972 * t1209 + t970 * t1218 + (t981 * t1213 + t980 * t1216 + t979 * t1220) * t1172) * t1174, -t964 * t1199 - t965 * t1197 - t966 * t1195 + (t972 * t1210 + t970 * t1219 + t971 * t1215 + (t981 * t1214 + t980 * t1217 + t979 * t1221) * t1172) * t1174, -t964 * t1270 - t965 * t1268 - t966 * t1266 + (t1203 * t981 + t1204 * t980 + t1205 * t979) * t1236 + t972 * t1206 + t971 * t1207 + t970 * t1208; t961 * t1198 + t962 * t1196 + t963 * t1194 + (t968 * t1211 + t969 * t1209 + t967 * t1218 + (t978 * t1213 + t977 * t1216 + t976 * t1220) * t1172) * t1174, -t961 * t1199 - t962 * t1197 - t963 * t1195 + m(4) + (t969 * t1210 + t967 * t1219 + t968 * t1215 + (t978 * t1214 + t977 * t1217 + t976 * t1221) * t1172) * t1174, -t961 * t1270 - t962 * t1268 - t963 * t1266 + (t1203 * t978 + t1204 * t977 + t1205 * t976) * t1236 + t969 * t1206 + t968 * t1207 + t967 * t1208; t973 * t1198 + t974 * t1196 + t975 * t1194 + (t983 * t1211 + t984 * t1209 + t982 * t1218 + (t987 * t1213 + t986 * t1216 + t985 * t1220) * t1172) * t1174, -t973 * t1199 - t974 * t1197 - t975 * t1195 + (t984 * t1210 + t982 * t1219 + t983 * t1215 + (t987 * t1214 + t986 * t1217 + t985 * t1221) * t1172) * t1174, -t973 * t1270 - t974 * t1268 - t975 * t1266 + m(4) + (t1203 * t987 + t1204 * t986 + t1205 * t985) * t1236 + t984 * t1206 + t983 * t1207 + t982 * t1208;];
MX  = t1;
