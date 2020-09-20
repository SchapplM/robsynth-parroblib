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
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 10:08
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRRRR7V2G2A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V2G2A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V2G2A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRRRR7V2G2A0_inertia_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR7V2G2A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR7V2G2A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRRRR7V2G2A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V2G2A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V2G2A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 09:44:55
% EndTime: 2020-08-07 09:45:00
% DurationCPUTime: 4.78s
% Computational Cost: add. (9210->462), mult. (15273->661), div. (792->14), fcn. (9090->66), ass. (0->309)
t1192 = pkin(6) + pkin(5);
t1117 = t1192 * mrSges(3,1) - Ifges(3,5);
t1173 = sin(qJ(3,3));
t1182 = cos(qJ(3,3));
t1361 = t1192 * mrSges(3,2) - Ifges(3,6);
t1373 = -t1173 * t1117 - t1182 * t1361;
t1176 = sin(qJ(3,2));
t1185 = cos(qJ(3,2));
t1372 = -t1176 * t1117 - t1185 * t1361;
t1179 = sin(qJ(3,1));
t1188 = cos(qJ(3,1));
t1371 = -t1179 * t1117 - t1188 * t1361;
t1370 = t1117 * t1182 - t1173 * t1361;
t1369 = t1117 * t1185 - t1176 * t1361;
t1368 = t1117 * t1188 - t1179 * t1361;
t1320 = 2 * pkin(1);
t1331 = qJ(3,1) + qJ(1,1);
t1130 = qJ(2,1) + t1331;
t1332 = -qJ(3,1) + qJ(1,1);
t1131 = -qJ(2,1) + t1332;
t1166 = qJ(2,1) + qJ(3,1);
t1180 = sin(qJ(2,1));
t1190 = cos(qJ(1,1));
t1203 = 0.2e1 * qJ(3,1);
t1208 = pkin(2) ^ 2;
t1325 = qJ(1,1) - 0.2e1 * qJ(2,1);
t1326 = qJ(1,1) + 0.2e1 * qJ(2,1);
t1344 = t1179 * pkin(1);
t1150 = (pkin(7) + t1192);
t1357 = 2 * t1150;
t1237 = ((cos(t1131) + cos(t1130)) * t1320 + (sin(t1131) + sin(t1130)) * t1357 + (cos(0.2e1 * qJ(3,1) - t1325) + cos(t1203 + t1326) + 0.2e1 * t1190) * pkin(3) + (cos(qJ(3,1) - t1325) + cos(qJ(3,1) + t1326) + cos(t1332) + cos(t1331)) * pkin(2)) / (-t1208 * sin(qJ(2,1) - qJ(3,1)) + pkin(2) * (pkin(2) * sin(t1166) + 0.2e1 * t1344 + (sin(t1203 + qJ(2,1)) - t1180) * pkin(3))) / 0.2e1;
t1329 = qJ(3,2) + qJ(1,2);
t1128 = qJ(2,2) + t1329;
t1330 = -qJ(3,2) + qJ(1,2);
t1129 = -qJ(2,2) + t1330;
t1165 = qJ(2,2) + qJ(3,2);
t1177 = sin(qJ(2,2));
t1187 = cos(qJ(1,2));
t1200 = 0.2e1 * qJ(3,2);
t1323 = qJ(1,2) - 0.2e1 * qJ(2,2);
t1324 = qJ(1,2) + 0.2e1 * qJ(2,2);
t1346 = t1176 * pkin(1);
t1238 = ((cos(t1129) + cos(t1128)) * t1320 + (sin(t1129) + sin(t1128)) * t1357 + (cos(0.2e1 * qJ(3,2) - t1323) + cos(t1200 + t1324) + 0.2e1 * t1187) * pkin(3) + (cos(qJ(3,2) - t1323) + cos(qJ(3,2) + t1324) + cos(t1330) + cos(t1329)) * pkin(2)) / (-t1208 * sin(qJ(2,2) - qJ(3,2)) + pkin(2) * (pkin(2) * sin(t1165) + 0.2e1 * t1346 + (sin(t1200 + qJ(2,2)) - t1177) * pkin(3))) / 0.2e1;
t1327 = qJ(3,3) + qJ(1,3);
t1126 = qJ(2,3) + t1327;
t1328 = -qJ(3,3) + qJ(1,3);
t1127 = -qJ(2,3) + t1328;
t1164 = qJ(2,3) + qJ(3,3);
t1174 = sin(qJ(2,3));
t1184 = cos(qJ(1,3));
t1197 = 0.2e1 * qJ(3,3);
t1321 = qJ(1,3) - 0.2e1 * qJ(2,3);
t1322 = 0.2e1 * qJ(2,3) + qJ(1,3);
t1348 = t1173 * pkin(1);
t1239 = ((cos(t1127) + cos(t1126)) * t1320 + (sin(t1127) + sin(t1126)) * t1357 + (cos(0.2e1 * qJ(3,3) - t1321) + cos(t1197 + t1322) + 0.2e1 * t1184) * pkin(3) + (cos(qJ(3,3) - t1321) + cos(qJ(3,3) + t1322) + cos(t1328) + cos(t1327)) * pkin(2)) / (-t1208 * sin(qJ(2,3) - qJ(3,3)) + pkin(2) * (pkin(2) * sin(t1164) + 0.2e1 * t1348 + (sin(t1197 + qJ(2,3)) - t1174) * pkin(3))) / 0.2e1;
t1360 = -2 * pkin(1);
t1359 = 0.2e1 * pkin(3);
t1358 = 4 * Ifges(3,4);
t1183 = cos(qJ(2,3));
t1159 = t1183 ^ 2;
t1356 = 0.2e1 * t1159;
t1186 = cos(qJ(2,2));
t1161 = t1186 ^ 2;
t1355 = 0.2e1 * t1161;
t1189 = cos(qJ(2,1));
t1163 = t1189 ^ 2;
t1354 = 0.2e1 * t1163;
t1353 = pkin(2) * mrSges(3,1);
t1352 = 2 * Ifges(2,4) - 2 * Ifges(3,4);
t1169 = Ifges(3,2) - Ifges(3,1);
t1351 = pkin(1) * t1174;
t1350 = pkin(1) * t1177;
t1349 = pkin(1) * t1180;
t1158 = t1182 ^ 2;
t1141 = pkin(3) * t1158;
t1160 = t1185 ^ 2;
t1142 = pkin(3) * t1160;
t1162 = t1188 ^ 2;
t1143 = pkin(3) * t1162;
t1347 = t1173 * pkin(3);
t1345 = t1176 * pkin(3);
t1343 = t1179 * pkin(3);
t1342 = t1182 * pkin(2);
t1138 = t1182 * pkin(3);
t1341 = t1185 * pkin(2);
t1139 = t1185 * pkin(3);
t1340 = t1188 * pkin(2);
t1140 = t1188 * pkin(3);
t1339 = mrSges(3,2) * t1173;
t1338 = mrSges(3,2) * t1176;
t1337 = mrSges(3,2) * t1179;
t1336 = Ifges(3,4) * t1173;
t1335 = Ifges(3,4) * t1176;
t1334 = Ifges(3,4) * t1179;
t1207 = 0.1e1 / pkin(3);
t1333 = Ifges(3,3) * t1207;
t1152 = 0.2e1 * t1353;
t1151 = -0.2e1 * pkin(2) * mrSges(3,2);
t1068 = -t1174 * t1370 + t1373 * t1183;
t1316 = t1068 * t1207;
t1069 = -t1177 * t1369 + t1372 * t1186;
t1315 = t1069 * t1207;
t1070 = -t1180 * t1368 + t1371 * t1189;
t1314 = t1070 * t1207;
t1175 = sin(qJ(1,3));
t1221 = pkin(1) * t1175 - t1184 * t1150;
t1268 = t1173 * t1174;
t1074 = t1221 * t1268 + (t1158 - 0.1e1) * t1175 * pkin(3);
t1170 = legFrame(3,2);
t1132 = sin(t1170);
t1313 = t1074 * t1132;
t1135 = cos(t1170);
t1312 = t1074 * t1135;
t1178 = sin(qJ(1,2));
t1220 = pkin(1) * t1178 - t1187 * t1150;
t1267 = t1176 * t1177;
t1075 = t1220 * t1267 + (t1160 - 0.1e1) * t1178 * pkin(3);
t1171 = legFrame(2,2);
t1133 = sin(t1171);
t1311 = t1075 * t1133;
t1136 = cos(t1171);
t1310 = t1075 * t1136;
t1181 = sin(qJ(1,1));
t1219 = pkin(1) * t1181 - t1190 * t1150;
t1266 = t1179 * t1180;
t1076 = t1219 * t1266 + (t1162 - 0.1e1) * t1181 * pkin(3);
t1172 = legFrame(1,2);
t1134 = sin(t1172);
t1309 = t1076 * t1134;
t1137 = cos(t1172);
t1308 = t1076 * t1137;
t1255 = pkin(3) * t1268;
t1122 = t1138 + pkin(2);
t1287 = t1122 * t1183;
t1307 = 0.1e1 / (pkin(1) - t1255 + t1287) / t1173;
t1254 = pkin(3) * t1267;
t1123 = t1139 + pkin(2);
t1284 = t1123 * t1186;
t1306 = 0.1e1 / (pkin(1) - t1254 + t1284) / t1176;
t1253 = pkin(3) * t1266;
t1124 = t1140 + pkin(2);
t1281 = t1124 * t1189;
t1305 = 0.1e1 / (pkin(1) - t1253 + t1281) / t1179;
t1095 = 0.1e1 / (pkin(3) * cos(t1164) + t1183 * pkin(2) + pkin(1));
t1304 = t1095 * t1175;
t1303 = t1095 * t1184;
t1096 = 0.1e1 / (pkin(3) * cos(t1165) + t1186 * pkin(2) + pkin(1));
t1302 = t1096 * t1178;
t1301 = t1096 * t1187;
t1097 = 0.1e1 / (pkin(3) * cos(t1166) + t1189 * pkin(2) + pkin(1));
t1300 = t1097 * t1181;
t1299 = t1097 * t1190;
t1098 = Ifges(3,3) + (mrSges(3,1) * t1182 - t1339) * pkin(2);
t1298 = t1098 * t1207;
t1099 = Ifges(3,3) + (mrSges(3,1) * t1185 - t1338) * pkin(2);
t1297 = t1099 * t1207;
t1100 = Ifges(3,3) + (mrSges(3,1) * t1188 - t1337) * pkin(2);
t1296 = t1100 * t1207;
t1107 = -t1347 + t1351;
t1295 = t1107 * t1182;
t1108 = -t1345 + t1350;
t1294 = t1108 * t1185;
t1109 = -t1343 + t1349;
t1293 = t1109 * t1188;
t1196 = pkin(2) / 0.2e1;
t1292 = (t1138 + t1196) * t1173;
t1291 = (t1139 + t1196) * t1176;
t1290 = (t1140 + t1196) * t1179;
t1289 = t1122 * t1132;
t1288 = t1122 * t1135;
t1286 = t1123 * t1133;
t1285 = t1123 * t1136;
t1283 = t1124 * t1134;
t1282 = t1124 * t1137;
t1280 = t1132 * t1175;
t1279 = t1133 * t1178;
t1278 = t1134 * t1181;
t1277 = t1135 * t1175;
t1276 = t1136 * t1178;
t1275 = t1137 * t1181;
t1206 = pkin(3) ^ 2;
t1274 = t1158 * t1206;
t1273 = t1160 * t1206;
t1272 = t1162 * t1206;
t1271 = t1169 * t1158;
t1270 = t1169 * t1160;
t1269 = t1169 * t1162;
t1209 = 0.1e1 / pkin(2);
t1265 = t1207 * t1209;
t1264 = -t1206 / 0.2e1 + t1208 / 0.2e1;
t1263 = pkin(2) * t1138;
t1262 = pkin(2) * t1139;
t1261 = pkin(2) * t1140;
t1193 = m(3) * t1208;
t1260 = Ifges(2,3) + Ifges(3,3) + t1193;
t1259 = -m(3) * pkin(2) - mrSges(2,1);
t1258 = t1122 * t1347;
t1257 = t1123 * t1345;
t1256 = t1124 * t1343;
t1252 = -pkin(5) * mrSges(2,2) + Ifges(2,6);
t1149 = -t1206 + t1208;
t1077 = pkin(1) * t1347 + (t1149 + 0.2e1 * t1263 + 0.2e1 * t1274) * t1174;
t1224 = t1175 * t1268;
t1083 = t1224 * t1359 - t1221;
t1092 = t1263 + t1264 + t1274;
t1044 = (t1092 * t1280 - t1135 * t1258) * t1356 + (-t1077 * t1135 - t1083 * t1289) * t1183 - pkin(3) * t1313 - t1107 * t1288;
t1251 = t1044 * t1307;
t1045 = (-t1092 * t1277 - t1132 * t1258) * t1356 + (-t1132 * t1077 + t1083 * t1288) * t1183 + pkin(3) * t1312 - t1107 * t1289;
t1250 = t1045 * t1307;
t1078 = pkin(1) * t1345 + (t1149 + 0.2e1 * t1262 + 0.2e1 * t1273) * t1177;
t1223 = t1178 * t1267;
t1084 = t1223 * t1359 - t1220;
t1093 = t1262 + t1264 + t1273;
t1046 = (t1093 * t1279 - t1136 * t1257) * t1355 + (-t1078 * t1136 - t1084 * t1286) * t1186 - pkin(3) * t1311 - t1108 * t1285;
t1249 = t1046 * t1306;
t1047 = (-t1093 * t1276 - t1133 * t1257) * t1355 + (-t1133 * t1078 + t1084 * t1285) * t1186 + pkin(3) * t1310 - t1108 * t1286;
t1248 = t1047 * t1306;
t1079 = pkin(1) * t1343 + (t1149 + 0.2e1 * t1261 + 0.2e1 * t1272) * t1180;
t1222 = t1181 * t1266;
t1085 = t1222 * t1359 - t1219;
t1094 = t1261 + t1264 + t1272;
t1048 = (t1094 * t1278 - t1137 * t1256) * t1354 + (-t1079 * t1137 - t1085 * t1283) * t1189 - pkin(3) * t1309 - t1109 * t1282;
t1247 = t1048 * t1305;
t1049 = (-t1094 * t1275 - t1134 * t1256) * t1354 + (-t1134 * t1079 + t1085 * t1282) * t1189 + pkin(3) * t1308 - t1109 * t1283;
t1246 = t1049 * t1305;
t1080 = t1348 + (-pkin(3) + t1342 + 0.2e1 * t1141) * t1174;
t1195 = -pkin(3) / 0.2e1;
t1101 = t1141 + t1342 / 0.2e1 + t1195;
t1213 = pkin(2) * t1224 + t1083 * t1182;
t1050 = (-t1101 * t1280 + t1135 * t1292) * t1356 + (t1135 * t1080 + t1213 * t1132) * t1183 + t1313 + t1135 * t1295;
t1245 = t1050 * t1307;
t1051 = (t1101 * t1277 + t1132 * t1292) * t1356 + (t1132 * t1080 - t1213 * t1135) * t1183 - t1312 + t1132 * t1295;
t1244 = t1051 * t1307;
t1081 = t1346 + (-pkin(3) + t1341 + 0.2e1 * t1142) * t1177;
t1102 = t1142 + t1341 / 0.2e1 + t1195;
t1212 = pkin(2) * t1223 + t1084 * t1185;
t1052 = (-t1102 * t1279 + t1136 * t1291) * t1355 + (t1136 * t1081 + t1212 * t1133) * t1186 + t1311 + t1136 * t1294;
t1243 = t1052 * t1306;
t1053 = (t1102 * t1276 + t1133 * t1291) * t1355 + (t1133 * t1081 - t1212 * t1136) * t1186 - t1310 + t1133 * t1294;
t1242 = t1053 * t1306;
t1082 = t1344 + (-pkin(3) + t1340 + 0.2e1 * t1143) * t1180;
t1103 = t1143 + t1340 / 0.2e1 + t1195;
t1211 = pkin(2) * t1222 + t1085 * t1188;
t1054 = (-t1103 * t1278 + t1137 * t1290) * t1354 + (t1137 * t1082 + t1211 * t1134) * t1189 + t1309 + t1137 * t1293;
t1241 = t1054 * t1305;
t1055 = (t1103 * t1275 + t1134 * t1290) * t1354 + (t1134 * t1082 - t1211 * t1137) * t1189 - t1308 + t1134 * t1293;
t1240 = t1055 * t1305;
t1062 = -0.2e1 * t1092 * t1184 * t1159 - ((pkin(1) - 0.2e1 * t1255) * t1184 + t1175 * t1150) * t1287 + pkin(3) * ((pkin(1) * t1268 - pkin(3) + t1141) * t1184 + t1150 * t1224);
t1236 = t1062 * t1307;
t1063 = -0.2e1 * t1093 * t1187 * t1161 - ((pkin(1) - 0.2e1 * t1254) * t1187 + t1178 * t1150) * t1284 + pkin(3) * ((pkin(1) * t1267 - pkin(3) + t1142) * t1187 + t1150 * t1223);
t1235 = t1063 * t1306;
t1064 = -0.2e1 * t1094 * t1190 * t1163 - ((pkin(1) - 0.2e1 * t1253) * t1190 + t1181 * t1150) * t1281 + pkin(3) * ((pkin(1) * t1266 - pkin(3) + t1143) * t1190 + t1150 * t1222);
t1234 = t1064 * t1305;
t1233 = t1209 * t1307;
t1232 = t1209 * t1306;
t1231 = t1209 * t1305;
t1230 = t1132 * t1303;
t1229 = t1135 * t1303;
t1228 = t1133 * t1301;
t1227 = t1136 * t1301;
t1226 = t1134 * t1299;
t1225 = t1137 * t1299;
t1218 = Ifges(2,2) + t1193 - Ifges(2,1) - t1169;
t1217 = t1062 * t1207 * t1233;
t1216 = t1063 * t1207 * t1232;
t1215 = t1064 * t1207 * t1231;
t1214 = Ifges(2,5) + t1259 * pkin(5) + (-m(3) * pkin(6) - mrSges(3,3)) * pkin(2);
t1210 = Ifges(2,1) + Ifges(3,2) + Ifges(1,3) + 0.2e1 * (mrSges(2,3) + mrSges(3,3)) * pkin(5) + t1192 ^ 2 * m(3) + (m(2) + m(3)) * (pkin(1) ^ 2) + 0.2e1 * mrSges(3,3) * pkin(6) + m(2) * pkin(5) ^ 2;
t1153 = mrSges(3,1) * t1320;
t1121 = t1179 * t1151;
t1120 = t1176 * t1151;
t1119 = t1173 * t1151;
t1091 = t1188 * t1152 + t1121 + t1260;
t1090 = t1185 * t1152 + t1120 + t1260;
t1089 = t1182 * t1152 + t1119 + t1260;
t1067 = (t1214 - t1368) * t1180 + (t1252 + t1371) * t1189;
t1066 = (t1214 - t1369) * t1177 + (t1252 + t1372) * t1186;
t1065 = (t1214 - t1370) * t1174 + (t1252 + t1373) * t1183;
t1058 = (0.2e1 * t1269 + (t1152 + 0.4e1 * t1334) * t1188 + t1121 + t1218) * t1163 + (t1153 * t1188 + (t1259 + t1337) * t1360 + (t1162 * t1358 + t1151 * t1188 + 0.2e1 * (-t1169 * t1188 - t1353) * t1179 + t1352) * t1180) * t1189 - t1269 + 0.2e1 * (-mrSges(3,2) * t1349 - t1334) * t1188 - 0.2e1 * (t1179 * mrSges(3,1) + mrSges(2,2)) * t1349 + t1210;
t1057 = (0.2e1 * t1270 + (t1152 + 0.4e1 * t1335) * t1185 + t1120 + t1218) * t1161 + (t1153 * t1185 + (t1259 + t1338) * t1360 + (t1160 * t1358 + t1151 * t1185 + 0.2e1 * (-t1169 * t1185 - t1353) * t1176 + t1352) * t1177) * t1186 - t1270 + 0.2e1 * (-mrSges(3,2) * t1350 - t1335) * t1185 - 0.2e1 * (t1176 * mrSges(3,1) + mrSges(2,2)) * t1350 + t1210;
t1056 = (0.2e1 * t1271 + (t1152 + 0.4e1 * t1336) * t1182 + t1119 + t1218) * t1159 + (t1153 * t1182 + (t1259 + t1339) * t1360 + (t1158 * t1358 + t1151 * t1182 + 0.2e1 * (-t1169 * t1182 - t1353) * t1173 + t1352) * t1174) * t1183 - t1271 + 0.2e1 * (-mrSges(3,2) * t1351 - t1336) * t1182 - 0.2e1 * (t1173 * mrSges(3,1) + mrSges(2,2)) * t1351 + t1210;
t1043 = Ifges(3,3) * t1215 - t1070 * t1300 + t1100 * t1237;
t1042 = Ifges(3,3) * t1216 - t1069 * t1302 + t1099 * t1238;
t1041 = Ifges(3,3) * t1217 - t1068 * t1304 + t1098 * t1239;
t1040 = -t1067 * t1300 + t1091 * t1237 + t1100 * t1215;
t1039 = -t1066 * t1302 + t1090 * t1238 + t1099 * t1216;
t1038 = -t1065 * t1304 + t1089 * t1239 + t1098 * t1217;
t1037 = t1070 * t1225 + (t1049 * t1333 + t1055 * t1100) * t1231;
t1036 = t1069 * t1227 + (t1047 * t1333 + t1053 * t1099) * t1232;
t1035 = t1068 * t1229 + (t1045 * t1333 + t1051 * t1098) * t1233;
t1034 = -t1070 * t1226 + (t1048 * t1333 + t1054 * t1100) * t1231;
t1033 = -t1069 * t1228 + (t1046 * t1333 + t1052 * t1099) * t1232;
t1032 = -t1068 * t1230 + (t1044 * t1333 + t1050 * t1098) * t1233;
t1031 = t1067 * t1225 + (t1049 * t1296 + t1055 * t1091) * t1231;
t1030 = t1066 * t1227 + (t1047 * t1297 + t1053 * t1090) * t1232;
t1029 = t1065 * t1229 + (t1045 * t1298 + t1051 * t1089) * t1233;
t1028 = -t1067 * t1226 + (t1048 * t1296 + t1054 * t1091) * t1231;
t1027 = -t1066 * t1228 + (t1046 * t1297 + t1052 * t1090) * t1232;
t1026 = -t1065 * t1230 + (t1044 * t1298 + t1050 * t1089) * t1233;
t1025 = -t1058 * t1300 + t1067 * t1237 + t1070 * t1215;
t1024 = -t1057 * t1302 + t1066 * t1238 + t1069 * t1216;
t1023 = -t1056 * t1304 + t1065 * t1239 + t1068 * t1217;
t1022 = t1058 * t1225 + (t1049 * t1314 + t1055 * t1067) * t1231;
t1021 = t1057 * t1227 + (t1047 * t1315 + t1053 * t1066) * t1232;
t1020 = t1056 * t1229 + (t1045 * t1316 + t1051 * t1065) * t1233;
t1019 = -t1058 * t1226 + (t1048 * t1314 + t1054 * t1067) * t1231;
t1018 = -t1057 * t1228 + (t1046 * t1315 + t1052 * t1066) * t1232;
t1017 = -t1056 * t1230 + (t1044 * t1316 + t1050 * t1065) * t1233;
t1 = [t1020 * t1229 + t1021 * t1227 + t1022 * t1225 + m(4) + (t1029 * t1244 + t1030 * t1242 + t1031 * t1240 + (t1035 * t1250 + t1036 * t1248 + t1037 * t1246) * t1207) * t1209, -t1020 * t1230 - t1021 * t1228 - t1022 * t1226 + (t1029 * t1245 + t1030 * t1243 + t1031 * t1241 + (t1035 * t1251 + t1036 * t1249 + t1037 * t1247) * t1207) * t1209, -t1020 * t1304 - t1021 * t1302 - t1022 * t1300 + (t1035 * t1236 + t1036 * t1235 + t1037 * t1234) * t1265 + t1029 * t1239 + t1030 * t1238 + t1031 * t1237; t1017 * t1229 + t1018 * t1227 + t1019 * t1225 + (t1026 * t1244 + t1027 * t1242 + t1028 * t1240 + (t1032 * t1250 + t1033 * t1248 + t1034 * t1246) * t1207) * t1209, -t1017 * t1230 - t1018 * t1228 - t1019 * t1226 + m(4) + (t1026 * t1245 + t1027 * t1243 + t1028 * t1241 + (t1032 * t1251 + t1033 * t1249 + t1034 * t1247) * t1207) * t1209, -t1017 * t1304 - t1018 * t1302 - t1019 * t1300 + (t1032 * t1236 + t1033 * t1235 + t1034 * t1234) * t1265 + t1026 * t1239 + t1027 * t1238 + t1028 * t1237; t1023 * t1229 + t1024 * t1227 + t1025 * t1225 + (t1038 * t1244 + t1039 * t1242 + t1040 * t1240 + (t1041 * t1250 + t1042 * t1248 + t1043 * t1246) * t1207) * t1209, -t1023 * t1230 - t1024 * t1228 - t1025 * t1226 + (t1038 * t1245 + t1039 * t1243 + t1040 * t1241 + (t1041 * t1251 + t1042 * t1249 + t1043 * t1247) * t1207) * t1209, -t1023 * t1304 - t1024 * t1302 - t1025 * t1300 + m(4) + (t1041 * t1236 + t1042 * t1235 + t1043 * t1234) * t1265 + t1038 * t1239 + t1039 * t1238 + t1040 * t1237;];
MX  = t1;
