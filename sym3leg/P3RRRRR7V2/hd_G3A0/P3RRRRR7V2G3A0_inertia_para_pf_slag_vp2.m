% Calculate inertia matrix for parallel robot
% P3RRRRR7V2G3A0
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
% Datum: 2020-08-07 10:47
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRRRR7V2G3A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V2G3A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V2G3A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRRRR7V2G3A0_inertia_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR7V2G3A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR7V2G3A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRRRR7V2G3A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V2G3A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V2G3A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 10:22:49
% EndTime: 2020-08-07 10:22:54
% DurationCPUTime: 4.73s
% Computational Cost: add. (9210->462), mult. (15273->658), div. (792->14), fcn. (9090->66), ass. (0->309)
t1198 = pkin(6) + pkin(5);
t1123 = t1198 * mrSges(3,1) - Ifges(3,5);
t1179 = sin(qJ(3,3));
t1188 = cos(qJ(3,3));
t1367 = t1198 * mrSges(3,2) - Ifges(3,6);
t1379 = -t1123 * t1179 - t1188 * t1367;
t1182 = sin(qJ(3,2));
t1191 = cos(qJ(3,2));
t1378 = -t1123 * t1182 - t1191 * t1367;
t1185 = sin(qJ(3,1));
t1194 = cos(qJ(3,1));
t1377 = -t1123 * t1185 - t1194 * t1367;
t1376 = t1123 * t1188 - t1179 * t1367;
t1375 = t1123 * t1191 - t1182 * t1367;
t1374 = t1123 * t1194 - t1185 * t1367;
t1326 = 2 * pkin(1);
t1337 = qJ(3,1) + qJ(1,1);
t1136 = qJ(2,1) + t1337;
t1338 = -qJ(3,1) + qJ(1,1);
t1137 = -qJ(2,1) + t1338;
t1172 = qJ(2,1) + qJ(3,1);
t1186 = sin(qJ(2,1));
t1187 = sin(qJ(1,1));
t1209 = 0.2e1 * qJ(3,1);
t1214 = pkin(2) ^ 2;
t1331 = qJ(1,1) - 0.2e1 * qJ(2,1);
t1332 = qJ(1,1) + 0.2e1 * qJ(2,1);
t1350 = t1185 * pkin(1);
t1156 = (pkin(7) + t1198);
t1363 = 2 * t1156;
t1240 = ((-sin(t1137) - sin(t1136)) * t1326 + (cos(t1137) + cos(t1136)) * t1363 + (sin(0.2e1 * qJ(3,1) - t1331) - sin(t1209 + t1332) - 0.2e1 * t1187) * pkin(3) + (sin(qJ(3,1) - t1331) - sin(qJ(3,1) + t1332) - sin(t1338) - sin(t1337)) * pkin(2)) / (-t1214 * sin(qJ(2,1) - qJ(3,1)) + pkin(2) * (pkin(2) * sin(t1172) + 0.2e1 * t1350 + (sin(t1209 + qJ(2,1)) - t1186) * pkin(3))) / 0.2e1;
t1335 = qJ(3,2) + qJ(1,2);
t1134 = qJ(2,2) + t1335;
t1336 = -qJ(3,2) + qJ(1,2);
t1135 = -qJ(2,2) + t1336;
t1171 = qJ(2,2) + qJ(3,2);
t1183 = sin(qJ(2,2));
t1184 = sin(qJ(1,2));
t1206 = 0.2e1 * qJ(3,2);
t1329 = qJ(1,2) - 0.2e1 * qJ(2,2);
t1330 = qJ(1,2) + 0.2e1 * qJ(2,2);
t1352 = t1182 * pkin(1);
t1241 = ((-sin(t1135) - sin(t1134)) * t1326 + (cos(t1135) + cos(t1134)) * t1363 + (sin(0.2e1 * qJ(3,2) - t1329) - sin(t1206 + t1330) - 0.2e1 * t1184) * pkin(3) + (sin(qJ(3,2) - t1329) - sin(qJ(3,2) + t1330) - sin(t1336) - sin(t1335)) * pkin(2)) / (-t1214 * sin(qJ(2,2) - qJ(3,2)) + pkin(2) * (pkin(2) * sin(t1171) + 0.2e1 * t1352 + (sin(t1206 + qJ(2,2)) - t1183) * pkin(3))) / 0.2e1;
t1333 = qJ(3,3) + qJ(1,3);
t1132 = qJ(2,3) + t1333;
t1334 = -qJ(3,3) + qJ(1,3);
t1133 = -qJ(2,3) + t1334;
t1170 = qJ(2,3) + qJ(3,3);
t1180 = sin(qJ(2,3));
t1181 = sin(qJ(1,3));
t1203 = 0.2e1 * qJ(3,3);
t1327 = -0.2e1 * qJ(2,3) + qJ(1,3);
t1328 = 0.2e1 * qJ(2,3) + qJ(1,3);
t1354 = t1179 * pkin(1);
t1242 = ((-sin(t1133) - sin(t1132)) * t1326 + (cos(t1133) + cos(t1132)) * t1363 + (sin(0.2e1 * qJ(3,3) - t1327) - sin(t1203 + t1328) - 0.2e1 * t1181) * pkin(3) + (sin(qJ(3,3) - t1327) - sin(qJ(3,3) + t1328) - sin(t1334) - sin(t1333)) * pkin(2)) / (-t1214 * sin(qJ(2,3) - qJ(3,3)) + pkin(2) * (pkin(2) * sin(t1170) + 0.2e1 * t1354 + (sin(t1203 + qJ(2,3)) - t1180) * pkin(3))) / 0.2e1;
t1366 = -2 * pkin(1);
t1365 = 0.2e1 * pkin(3);
t1364 = 4 * Ifges(3,4);
t1189 = cos(qJ(2,3));
t1165 = t1189 ^ 2;
t1362 = 0.2e1 * t1165;
t1192 = cos(qJ(2,2));
t1167 = t1192 ^ 2;
t1361 = 0.2e1 * t1167;
t1195 = cos(qJ(2,1));
t1169 = t1195 ^ 2;
t1360 = 0.2e1 * t1169;
t1359 = pkin(2) * mrSges(3,1);
t1358 = 2 * Ifges(2,4) - 2 * Ifges(3,4);
t1175 = Ifges(3,2) - Ifges(3,1);
t1357 = pkin(1) * t1180;
t1356 = pkin(1) * t1183;
t1355 = pkin(1) * t1186;
t1164 = t1188 ^ 2;
t1147 = pkin(3) * t1164;
t1166 = t1191 ^ 2;
t1148 = pkin(3) * t1166;
t1168 = t1194 ^ 2;
t1149 = pkin(3) * t1168;
t1353 = t1179 * pkin(3);
t1351 = t1182 * pkin(3);
t1349 = t1185 * pkin(3);
t1348 = t1188 * pkin(2);
t1144 = t1188 * pkin(3);
t1347 = t1191 * pkin(2);
t1145 = t1191 * pkin(3);
t1346 = t1194 * pkin(2);
t1146 = t1194 * pkin(3);
t1345 = mrSges(3,2) * t1179;
t1344 = mrSges(3,2) * t1182;
t1343 = mrSges(3,2) * t1185;
t1342 = Ifges(3,4) * t1179;
t1341 = Ifges(3,4) * t1182;
t1340 = Ifges(3,4) * t1185;
t1213 = 0.1e1 / pkin(3);
t1339 = Ifges(3,3) * t1213;
t1158 = 0.2e1 * t1359;
t1157 = -0.2e1 * pkin(2) * mrSges(3,2);
t1071 = -t1180 * t1376 + t1379 * t1189;
t1322 = t1071 * t1213;
t1072 = -t1183 * t1375 + t1378 * t1192;
t1321 = t1072 * t1213;
t1073 = -t1186 * t1374 + t1377 * t1195;
t1320 = t1073 * t1213;
t1190 = cos(qJ(1,3));
t1270 = pkin(1) * t1190 + t1181 * t1156;
t1277 = t1179 * t1180;
t1077 = t1270 * t1277 + (t1164 - 0.1e1) * t1190 * pkin(3);
t1176 = legFrame(3,2);
t1138 = sin(t1176);
t1319 = t1077 * t1138;
t1141 = cos(t1176);
t1318 = t1077 * t1141;
t1193 = cos(qJ(1,2));
t1269 = pkin(1) * t1193 + t1184 * t1156;
t1276 = t1182 * t1183;
t1078 = t1269 * t1276 + (t1166 - 0.1e1) * t1193 * pkin(3);
t1177 = legFrame(2,2);
t1139 = sin(t1177);
t1317 = t1078 * t1139;
t1142 = cos(t1177);
t1316 = t1078 * t1142;
t1196 = cos(qJ(1,1));
t1268 = pkin(1) * t1196 + t1187 * t1156;
t1275 = t1185 * t1186;
t1079 = t1268 * t1275 + (t1168 - 0.1e1) * t1196 * pkin(3);
t1178 = legFrame(1,2);
t1140 = sin(t1178);
t1315 = t1079 * t1140;
t1143 = cos(t1178);
t1314 = t1079 * t1143;
t1258 = pkin(3) * t1277;
t1128 = t1144 + pkin(2);
t1296 = t1128 * t1189;
t1313 = 0.1e1 / (pkin(1) - t1258 + t1296) / t1179;
t1257 = pkin(3) * t1276;
t1129 = t1145 + pkin(2);
t1293 = t1129 * t1192;
t1312 = 0.1e1 / (pkin(1) - t1257 + t1293) / t1182;
t1256 = pkin(3) * t1275;
t1130 = t1146 + pkin(2);
t1290 = t1130 * t1195;
t1311 = 0.1e1 / (pkin(1) - t1256 + t1290) / t1185;
t1098 = 0.1e1 / (pkin(3) * cos(t1170) + t1189 * pkin(2) + pkin(1));
t1310 = t1098 * t1181;
t1309 = t1098 * t1190;
t1099 = 0.1e1 / (pkin(3) * cos(t1171) + t1192 * pkin(2) + pkin(1));
t1308 = t1099 * t1184;
t1307 = t1099 * t1193;
t1100 = 0.1e1 / (pkin(3) * cos(t1172) + t1195 * pkin(2) + pkin(1));
t1306 = t1100 * t1187;
t1305 = t1100 * t1196;
t1101 = Ifges(3,3) + (mrSges(3,1) * t1188 - t1345) * pkin(2);
t1304 = t1101 * t1213;
t1102 = Ifges(3,3) + (mrSges(3,1) * t1191 - t1344) * pkin(2);
t1303 = t1102 * t1213;
t1103 = Ifges(3,3) + (mrSges(3,1) * t1194 - t1343) * pkin(2);
t1302 = t1103 * t1213;
t1202 = pkin(2) / 0.2e1;
t1301 = (t1144 + t1202) * t1179;
t1300 = (t1145 + t1202) * t1182;
t1299 = (t1146 + t1202) * t1185;
t1298 = t1128 * t1138;
t1297 = t1128 * t1141;
t1295 = t1129 * t1139;
t1294 = t1129 * t1142;
t1292 = t1130 * t1140;
t1291 = t1130 * t1143;
t1289 = t1138 * t1190;
t1288 = t1139 * t1193;
t1287 = t1140 * t1196;
t1286 = t1141 * t1190;
t1285 = t1142 * t1193;
t1284 = t1143 * t1196;
t1212 = pkin(3) ^ 2;
t1283 = t1164 * t1212;
t1282 = t1166 * t1212;
t1281 = t1168 * t1212;
t1280 = t1175 * t1164;
t1279 = t1175 * t1166;
t1278 = t1175 * t1168;
t1110 = -t1353 + t1357;
t1274 = t1188 * t1110;
t1111 = -t1351 + t1356;
t1273 = t1191 * t1111;
t1112 = -t1349 + t1355;
t1272 = t1194 * t1112;
t1215 = 0.1e1 / pkin(2);
t1271 = t1213 * t1215;
t1267 = -t1212 / 0.2e1 + t1214 / 0.2e1;
t1266 = pkin(2) * t1144;
t1265 = pkin(2) * t1145;
t1264 = pkin(2) * t1146;
t1199 = m(3) * t1214;
t1263 = Ifges(2,3) + Ifges(3,3) + t1199;
t1262 = -m(3) * pkin(2) - mrSges(2,1);
t1261 = t1128 * t1353;
t1260 = t1129 * t1351;
t1259 = t1130 * t1349;
t1255 = -pkin(5) * mrSges(2,2) + Ifges(2,6);
t1155 = -t1212 + t1214;
t1080 = pkin(1) * t1353 + (t1155 + 0.2e1 * t1266 + 0.2e1 * t1283) * t1180;
t1227 = t1190 * t1277;
t1086 = t1227 * t1365 - t1270;
t1095 = t1266 + t1267 + t1283;
t1047 = (t1095 * t1289 - t1141 * t1261) * t1362 + (-t1141 * t1080 - t1086 * t1298) * t1189 - pkin(3) * t1319 - t1110 * t1297;
t1254 = t1047 * t1313;
t1048 = (-t1095 * t1286 - t1138 * t1261) * t1362 + (-t1138 * t1080 + t1086 * t1297) * t1189 + pkin(3) * t1318 - t1110 * t1298;
t1253 = t1048 * t1313;
t1081 = pkin(1) * t1351 + (t1155 + 0.2e1 * t1265 + 0.2e1 * t1282) * t1183;
t1226 = t1193 * t1276;
t1087 = t1226 * t1365 - t1269;
t1096 = t1265 + t1267 + t1282;
t1049 = (t1096 * t1288 - t1142 * t1260) * t1361 + (-t1142 * t1081 - t1087 * t1295) * t1192 - pkin(3) * t1317 - t1111 * t1294;
t1252 = t1049 * t1312;
t1050 = (-t1096 * t1285 - t1139 * t1260) * t1361 + (-t1139 * t1081 + t1087 * t1294) * t1192 + pkin(3) * t1316 - t1111 * t1295;
t1251 = t1050 * t1312;
t1082 = pkin(1) * t1349 + (t1155 + 0.2e1 * t1264 + 0.2e1 * t1281) * t1186;
t1225 = t1196 * t1275;
t1088 = t1225 * t1365 - t1268;
t1097 = t1264 + t1267 + t1281;
t1051 = (t1097 * t1287 - t1143 * t1259) * t1360 + (-t1143 * t1082 - t1088 * t1292) * t1195 - pkin(3) * t1315 - t1112 * t1291;
t1250 = t1051 * t1311;
t1052 = (-t1097 * t1284 - t1140 * t1259) * t1360 + (-t1140 * t1082 + t1088 * t1291) * t1195 + pkin(3) * t1314 - t1112 * t1292;
t1249 = t1052 * t1311;
t1083 = t1354 + (-pkin(3) + t1348 + 0.2e1 * t1147) * t1180;
t1201 = -pkin(3) / 0.2e1;
t1104 = t1147 + t1348 / 0.2e1 + t1201;
t1219 = pkin(2) * t1227 + t1086 * t1188;
t1053 = (-t1104 * t1289 + t1141 * t1301) * t1362 + (t1141 * t1083 + t1219 * t1138) * t1189 + t1319 + t1141 * t1274;
t1248 = t1053 * t1313;
t1054 = (t1104 * t1286 + t1138 * t1301) * t1362 + (t1138 * t1083 - t1219 * t1141) * t1189 - t1318 + t1138 * t1274;
t1247 = t1054 * t1313;
t1084 = t1352 + (-pkin(3) + t1347 + 0.2e1 * t1148) * t1183;
t1105 = t1148 + t1347 / 0.2e1 + t1201;
t1218 = pkin(2) * t1226 + t1087 * t1191;
t1055 = (-t1105 * t1288 + t1142 * t1300) * t1361 + (t1142 * t1084 + t1218 * t1139) * t1192 + t1317 + t1142 * t1273;
t1246 = t1055 * t1312;
t1056 = (t1105 * t1285 + t1139 * t1300) * t1361 + (t1139 * t1084 - t1218 * t1142) * t1192 - t1316 + t1139 * t1273;
t1245 = t1056 * t1312;
t1085 = t1350 + (-pkin(3) + t1346 + 0.2e1 * t1149) * t1186;
t1106 = t1149 + t1346 / 0.2e1 + t1201;
t1217 = pkin(2) * t1225 + t1088 * t1194;
t1057 = (-t1106 * t1287 + t1143 * t1299) * t1360 + (t1143 * t1085 + t1217 * t1140) * t1195 + t1315 + t1143 * t1272;
t1244 = t1057 * t1311;
t1058 = (t1106 * t1284 + t1140 * t1299) * t1360 + (t1140 * t1085 - t1217 * t1143) * t1195 - t1314 + t1140 * t1272;
t1243 = t1058 * t1311;
t1065 = t1095 * t1181 * t1362 + ((pkin(1) - 0.2e1 * t1258) * t1181 - t1190 * t1156) * t1296 - pkin(3) * ((pkin(1) * t1277 - pkin(3) + t1147) * t1181 - t1156 * t1227);
t1239 = t1065 * t1313;
t1066 = t1096 * t1184 * t1361 + ((pkin(1) - 0.2e1 * t1257) * t1184 - t1193 * t1156) * t1293 - pkin(3) * ((pkin(1) * t1276 - pkin(3) + t1148) * t1184 - t1156 * t1226);
t1238 = t1066 * t1312;
t1067 = t1097 * t1187 * t1360 + ((pkin(1) - 0.2e1 * t1256) * t1187 - t1196 * t1156) * t1290 - pkin(3) * ((pkin(1) * t1275 - pkin(3) + t1149) * t1187 - t1156 * t1225);
t1237 = t1067 * t1311;
t1236 = t1215 * t1313;
t1235 = t1215 * t1312;
t1234 = t1215 * t1311;
t1233 = t1138 * t1310;
t1232 = t1141 * t1310;
t1231 = t1139 * t1308;
t1230 = t1142 * t1308;
t1229 = t1140 * t1306;
t1228 = t1143 * t1306;
t1224 = Ifges(2,2) + t1199 - Ifges(2,1) - t1175;
t1223 = t1065 * t1213 * t1236;
t1222 = t1066 * t1213 * t1235;
t1221 = t1067 * t1213 * t1234;
t1220 = Ifges(2,5) + t1262 * pkin(5) + (-m(3) * pkin(6) - mrSges(3,3)) * pkin(2);
t1216 = Ifges(2,1) + Ifges(3,2) + Ifges(1,3) + 0.2e1 * (mrSges(2,3) + mrSges(3,3)) * pkin(5) + t1198 ^ 2 * m(3) + (m(2) + m(3)) * (pkin(1) ^ 2) + 0.2e1 * mrSges(3,3) * pkin(6) + m(2) * pkin(5) ^ 2;
t1159 = mrSges(3,1) * t1326;
t1127 = t1185 * t1157;
t1126 = t1182 * t1157;
t1125 = t1179 * t1157;
t1094 = t1194 * t1158 + t1127 + t1263;
t1093 = t1191 * t1158 + t1126 + t1263;
t1092 = t1188 * t1158 + t1125 + t1263;
t1070 = (t1220 - t1374) * t1186 + (t1255 + t1377) * t1195;
t1069 = (t1220 - t1375) * t1183 + (t1255 + t1378) * t1192;
t1068 = (t1220 - t1376) * t1180 + (t1255 + t1379) * t1189;
t1061 = (0.2e1 * t1278 + (t1158 + 0.4e1 * t1340) * t1194 + t1127 + t1224) * t1169 + (t1159 * t1194 + (t1262 + t1343) * t1366 + (t1168 * t1364 + t1157 * t1194 + 0.2e1 * (-t1175 * t1194 - t1359) * t1185 + t1358) * t1186) * t1195 - t1278 + 0.2e1 * (-mrSges(3,2) * t1355 - t1340) * t1194 - 0.2e1 * (t1185 * mrSges(3,1) + mrSges(2,2)) * t1355 + t1216;
t1060 = (0.2e1 * t1279 + (t1158 + 0.4e1 * t1341) * t1191 + t1126 + t1224) * t1167 + (t1159 * t1191 + (t1262 + t1344) * t1366 + (t1166 * t1364 + t1157 * t1191 + 0.2e1 * (-t1175 * t1191 - t1359) * t1182 + t1358) * t1183) * t1192 - t1279 + 0.2e1 * (-mrSges(3,2) * t1356 - t1341) * t1191 - 0.2e1 * (t1182 * mrSges(3,1) + mrSges(2,2)) * t1356 + t1216;
t1059 = (0.2e1 * t1280 + (t1158 + 0.4e1 * t1342) * t1188 + t1125 + t1224) * t1165 + (t1159 * t1188 + (t1262 + t1345) * t1366 + (t1164 * t1364 + t1157 * t1188 + 0.2e1 * (-t1175 * t1188 - t1359) * t1179 + t1358) * t1180) * t1189 - t1280 + 0.2e1 * (-mrSges(3,2) * t1357 - t1342) * t1188 - 0.2e1 * (t1179 * mrSges(3,1) + mrSges(2,2)) * t1357 + t1216;
t1046 = Ifges(3,3) * t1221 - t1073 * t1305 + t1103 * t1240;
t1045 = Ifges(3,3) * t1222 - t1072 * t1307 + t1102 * t1241;
t1044 = Ifges(3,3) * t1223 - t1071 * t1309 + t1101 * t1242;
t1043 = -t1070 * t1305 + t1094 * t1240 + t1103 * t1221;
t1042 = -t1069 * t1307 + t1093 * t1241 + t1102 * t1222;
t1041 = -t1068 * t1309 + t1092 * t1242 + t1101 * t1223;
t1040 = -t1073 * t1228 + (t1052 * t1339 + t1058 * t1103) * t1234;
t1039 = -t1072 * t1230 + (t1050 * t1339 + t1056 * t1102) * t1235;
t1038 = -t1071 * t1232 + (t1048 * t1339 + t1054 * t1101) * t1236;
t1037 = t1073 * t1229 + (t1051 * t1339 + t1057 * t1103) * t1234;
t1036 = t1072 * t1231 + (t1049 * t1339 + t1055 * t1102) * t1235;
t1035 = t1071 * t1233 + (t1047 * t1339 + t1053 * t1101) * t1236;
t1034 = -t1070 * t1228 + (t1052 * t1302 + t1058 * t1094) * t1234;
t1033 = -t1069 * t1230 + (t1050 * t1303 + t1056 * t1093) * t1235;
t1032 = -t1068 * t1232 + (t1048 * t1304 + t1054 * t1092) * t1236;
t1031 = t1070 * t1229 + (t1051 * t1302 + t1057 * t1094) * t1234;
t1030 = t1069 * t1231 + (t1049 * t1303 + t1055 * t1093) * t1235;
t1029 = t1068 * t1233 + (t1047 * t1304 + t1053 * t1092) * t1236;
t1028 = -t1061 * t1305 + t1070 * t1240 + t1073 * t1221;
t1027 = -t1060 * t1307 + t1069 * t1241 + t1072 * t1222;
t1026 = -t1059 * t1309 + t1068 * t1242 + t1071 * t1223;
t1025 = -t1061 * t1228 + (t1052 * t1320 + t1058 * t1070) * t1234;
t1024 = -t1060 * t1230 + (t1050 * t1321 + t1056 * t1069) * t1235;
t1023 = -t1059 * t1232 + (t1048 * t1322 + t1054 * t1068) * t1236;
t1022 = t1061 * t1229 + (t1051 * t1320 + t1057 * t1070) * t1234;
t1021 = t1060 * t1231 + (t1049 * t1321 + t1055 * t1069) * t1235;
t1020 = t1059 * t1233 + (t1047 * t1322 + t1053 * t1068) * t1236;
t1 = [-t1023 * t1232 - t1024 * t1230 - t1025 * t1228 + m(4) + (t1032 * t1247 + t1033 * t1245 + t1034 * t1243 + (t1038 * t1253 + t1039 * t1251 + t1040 * t1249) * t1213) * t1215, t1023 * t1233 + t1024 * t1231 + t1025 * t1229 + (t1032 * t1248 + t1033 * t1246 + t1034 * t1244 + (t1038 * t1254 + t1039 * t1252 + t1040 * t1250) * t1213) * t1215, -t1023 * t1309 - t1024 * t1307 - t1025 * t1305 + (t1038 * t1239 + t1039 * t1238 + t1040 * t1237) * t1271 + t1032 * t1242 + t1033 * t1241 + t1034 * t1240; -t1020 * t1232 - t1021 * t1230 - t1022 * t1228 + (t1029 * t1247 + t1030 * t1245 + t1031 * t1243 + (t1035 * t1253 + t1036 * t1251 + t1037 * t1249) * t1213) * t1215, t1020 * t1233 + t1021 * t1231 + t1022 * t1229 + m(4) + (t1029 * t1248 + t1030 * t1246 + t1031 * t1244 + (t1035 * t1254 + t1036 * t1252 + t1037 * t1250) * t1213) * t1215, -t1020 * t1309 - t1021 * t1307 - t1022 * t1305 + (t1035 * t1239 + t1036 * t1238 + t1037 * t1237) * t1271 + t1029 * t1242 + t1030 * t1241 + t1031 * t1240; -t1026 * t1232 - t1027 * t1230 - t1028 * t1228 + (t1041 * t1247 + t1042 * t1245 + t1043 * t1243 + (t1044 * t1253 + t1045 * t1251 + t1046 * t1249) * t1213) * t1215, t1026 * t1233 + t1027 * t1231 + t1028 * t1229 + (t1041 * t1248 + t1042 * t1246 + t1043 * t1244 + (t1044 * t1254 + t1045 * t1252 + t1046 * t1250) * t1213) * t1215, -t1026 * t1309 - t1027 * t1307 - t1028 * t1305 + m(4) + (t1044 * t1239 + t1045 * t1238 + t1046 * t1237) * t1271 + t1041 * t1242 + t1042 * t1241 + t1043 * t1240;];
MX  = t1;
