% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RPRR1G3P3A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [8x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPRR1G3P3A0_convert_par2_MPV_fixb.m

% Output:
% tauX [3x1]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:27
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRR1G3P3A0_invdyn_para_pf_mdp(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_mdp: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:27:05
% EndTime: 2020-03-09 21:27:09
% DurationCPUTime: 4.34s
% Computational Cost: add. (57783->347), mult. (43497->516), div. (2799->10), fcn. (22812->71), ass. (0->245)
t1305 = legFrame(3,2);
t1361 = qJ(1,3) + pkin(7);
t1253 = t1305 + t1361;
t1254 = -t1305 + t1361;
t1284 = qJ(1,3) + t1305;
t1285 = qJ(1,3) - t1305;
t1247 = qJ(3,3) + t1253;
t1220 = sin(t1247);
t1248 = qJ(3,3) + t1254;
t1221 = sin(t1248);
t1360 = t1220 - t1221;
t1138 = (-sin(t1285) + sin(t1284)) * pkin(1) + (-sin(t1254) + sin(t1253)) * pkin(2) + t1360 * pkin(3);
t1306 = legFrame(2,2);
t1362 = qJ(1,2) + pkin(7);
t1255 = t1306 + t1362;
t1256 = -t1306 + t1362;
t1286 = qJ(1,2) + t1306;
t1287 = qJ(1,2) - t1306;
t1249 = qJ(3,2) + t1255;
t1222 = sin(t1249);
t1250 = qJ(3,2) + t1256;
t1223 = sin(t1250);
t1359 = t1222 - t1223;
t1139 = (-sin(t1287) + sin(t1286)) * pkin(1) + (-sin(t1256) + sin(t1255)) * pkin(2) + t1359 * pkin(3);
t1307 = legFrame(1,2);
t1363 = qJ(1,1) + pkin(7);
t1257 = t1307 + t1363;
t1258 = -t1307 + t1363;
t1288 = qJ(1,1) + t1307;
t1289 = qJ(1,1) - t1307;
t1251 = qJ(3,1) + t1257;
t1224 = sin(t1251);
t1252 = qJ(3,1) + t1258;
t1225 = sin(t1252);
t1358 = t1224 - t1225;
t1140 = (-sin(t1289) + sin(t1288)) * pkin(1) + (-sin(t1258) + sin(t1257)) * pkin(2) + t1358 * pkin(3);
t1445 = MDP(1) / 0.2e1;
t1444 = MDP(2) / 0.2e1;
t1443 = MDP(3) / 0.2e1;
t1442 = MDP(5) / 0.2e1;
t1441 = MDP(6) / 0.2e1;
t1440 = MDP(7) / 0.2e1;
t1329 = 0.1e1 / pkin(3);
t1439 = t1329 / 0.2e1;
t1422 = MDP(4) * pkin(1);
t1438 = t1422 / 0.2e1;
t1304 = cos(pkin(7));
t1437 = pkin(2) * t1304;
t1230 = cos(t1251);
t1231 = cos(t1252);
t1167 = t1231 + t1230;
t1302 = pkin(7) + qJ(3,1);
t1283 = qJ(1,1) + t1302;
t1261 = sin(t1283);
t1326 = xDP(2);
t1327 = xDP(1);
t1325 = xDP(3);
t1425 = -2 * t1325;
t1136 = t1167 * t1327 + t1261 * t1425 - t1358 * t1326;
t1274 = sin(t1302);
t1315 = sin(qJ(3,1));
t1213 = pkin(1) * t1274 + pkin(2) * t1315;
t1185 = 0.1e1 / t1213;
t1130 = t1136 * t1185 / 0.2e1;
t1280 = cos(t1302);
t1321 = cos(qJ(3,1));
t1143 = -t1167 * pkin(3) + (-cos(t1257) - cos(t1258)) * pkin(2) + (-cos(t1288) - cos(t1289)) * pkin(1);
t1316 = sin(qJ(1,1));
t1161 = pkin(1) * t1316 + pkin(2) * sin(t1363) + pkin(3) * t1261;
t1357 = 2 * t1325;
t1400 = (t1140 * t1326 + t1143 * t1327 + t1161 * t1357) * t1185 * t1329;
t1349 = t1400 / 0.2e1;
t1310 = xDDP(1);
t1367 = t1310 / 0.2e1;
t1309 = xDDP(2);
t1368 = t1309 / 0.2e1;
t1375 = t1167 * t1185;
t1378 = t1358 * t1185;
t1186 = 0.1e1 / t1213 ^ 2;
t1397 = t1136 * t1186;
t1124 = t1130 + t1349;
t1407 = pkin(3) * t1124;
t1433 = t1185 * t1349 * t1407 + t1367 * t1375 - t1368 * t1378 - (-t1407 + (-pkin(1) * t1280 - pkin(2) * t1321) * t1130) * t1397 / 0.2e1;
t1228 = cos(t1249);
t1229 = cos(t1250);
t1166 = t1229 + t1228;
t1301 = pkin(7) + qJ(3,2);
t1282 = qJ(1,2) + t1301;
t1260 = sin(t1282);
t1137 = t1166 * t1327 + t1260 * t1425 - t1359 * t1326;
t1273 = sin(t1301);
t1313 = sin(qJ(3,2));
t1212 = pkin(1) * t1273 + pkin(2) * t1313;
t1183 = 0.1e1 / t1212;
t1131 = t1137 * t1183 / 0.2e1;
t1279 = cos(t1301);
t1319 = cos(qJ(3,2));
t1142 = -t1166 * pkin(3) + (-cos(t1255) - cos(t1256)) * pkin(2) + (-cos(t1286) - cos(t1287)) * pkin(1);
t1314 = sin(qJ(1,2));
t1160 = pkin(1) * t1314 + pkin(2) * sin(t1362) + pkin(3) * t1260;
t1399 = (t1139 * t1326 + t1142 * t1327 + t1160 * t1357) * t1183 * t1329;
t1348 = t1399 / 0.2e1;
t1376 = t1166 * t1183;
t1379 = t1359 * t1183;
t1184 = 0.1e1 / t1212 ^ 2;
t1396 = t1137 * t1184;
t1125 = t1131 + t1348;
t1406 = pkin(3) * t1125;
t1432 = t1183 * t1348 * t1406 + t1367 * t1376 - t1368 * t1379 - (-t1406 + (-pkin(1) * t1279 - pkin(2) * t1319) * t1131) * t1396 / 0.2e1;
t1226 = cos(t1247);
t1227 = cos(t1248);
t1165 = t1227 + t1226;
t1300 = pkin(7) + qJ(3,3);
t1281 = qJ(1,3) + t1300;
t1259 = sin(t1281);
t1135 = t1165 * t1327 + t1259 * t1425 - t1360 * t1326;
t1272 = sin(t1300);
t1311 = sin(qJ(3,3));
t1211 = pkin(1) * t1272 + pkin(2) * t1311;
t1181 = 0.1e1 / t1211;
t1129 = t1135 * t1181 / 0.2e1;
t1278 = cos(t1300);
t1317 = cos(qJ(3,3));
t1141 = -t1165 * pkin(3) + (-cos(t1253) - cos(t1254)) * pkin(2) + (-cos(t1284) - cos(t1285)) * pkin(1);
t1312 = sin(qJ(1,3));
t1159 = pkin(1) * t1312 + pkin(2) * sin(t1361) + pkin(3) * t1259;
t1401 = (t1138 * t1326 + t1141 * t1327 + t1159 * t1357) * t1181 * t1329;
t1350 = t1401 / 0.2e1;
t1377 = t1165 * t1181;
t1380 = t1360 * t1181;
t1182 = 0.1e1 / t1211 ^ 2;
t1398 = t1135 * t1182;
t1123 = t1129 + t1350;
t1408 = pkin(3) * t1123;
t1431 = t1181 * t1350 * t1408 + t1367 * t1377 - t1368 * t1380 - (-t1408 + (-pkin(1) * t1278 - pkin(2) * t1317) * t1129) * t1398 / 0.2e1;
t1427 = -0.2e1 * pkin(2);
t1426 = 0.2e1 * pkin(2);
t1424 = g(1) / 0.2e1;
t1423 = -g(2) / 0.2e1;
t1421 = t1221 / 0.2e1;
t1420 = t1223 / 0.2e1;
t1419 = t1225 / 0.2e1;
t1418 = t1226 / 0.2e1;
t1417 = t1228 / 0.2e1;
t1416 = t1230 / 0.2e1;
t1409 = pkin(1) * sin(pkin(7));
t1405 = t1309 - g(2);
t1404 = t1310 - g(1);
t1403 = pkin(3) * t1426;
t1402 = 0.2e1 * pkin(1);
t1395 = t1138 * t1181;
t1394 = t1138 * t1309;
t1393 = t1139 * t1183;
t1392 = t1139 * t1309;
t1391 = t1140 * t1185;
t1390 = t1140 * t1309;
t1389 = t1141 * t1181;
t1388 = t1141 * t1310;
t1387 = t1142 * t1183;
t1386 = t1142 * t1310;
t1385 = t1143 * t1185;
t1384 = t1143 * t1310;
t1383 = t1159 * t1181;
t1382 = t1160 * t1183;
t1381 = t1161 * t1185;
t1374 = t1181 * t1259;
t1308 = xDDP(3);
t1373 = t1181 * t1308;
t1372 = t1183 * t1260;
t1371 = t1183 * t1308;
t1370 = t1185 * t1308;
t1369 = t1261 * t1185;
t1366 = t1329 * t1159;
t1365 = t1329 * t1160;
t1364 = t1329 * t1161;
t1290 = sin(t1305);
t1293 = cos(t1305);
t1178 = g(1) * t1293 - g(2) * t1290;
t1318 = cos(qJ(1,3));
t1153 = g(3) * t1318 + t1178 * t1312;
t1120 = t1129 + t1350 / 0.2e1;
t1299 = pkin(1) ^ 2 + pkin(2) ^ 2;
t1328 = pkin(3) ^ 2;
t1356 = (t1120 * t1317 * t1403 + t1299 * t1129 + t1123 * t1328 + (pkin(3) * t1120 * t1278 + t1129 * t1437) * t1402) * t1398;
t1121 = t1130 + t1349 / 0.2e1;
t1355 = (t1121 * t1321 * t1403 + t1299 * t1130 + t1124 * t1328 + (pkin(3) * t1121 * t1280 + t1130 * t1437) * t1402) * t1397;
t1122 = t1131 + t1348 / 0.2e1;
t1354 = (t1122 * t1319 * t1403 + t1299 * t1131 + t1125 * t1328 + (pkin(3) * t1122 * t1279 + t1131 * t1437) * t1402) * t1396;
t1347 = t1135 ^ 2 * t1182 / 0.4e1;
t1346 = t1136 ^ 2 * t1186 / 0.4e1;
t1345 = t1137 ^ 2 * t1184 / 0.4e1;
t1291 = sin(t1306);
t1294 = cos(t1306);
t1179 = g(1) * t1294 - g(2) * t1291;
t1320 = cos(qJ(1,2));
t1155 = g(3) * t1320 + t1179 * t1314;
t1292 = sin(t1307);
t1295 = cos(t1307);
t1180 = g(1) * t1295 - g(2) * t1292;
t1322 = cos(qJ(1,1));
t1157 = g(3) * t1322 + t1180 * t1316;
t1341 = t1120 * t1350;
t1340 = t1121 * t1349;
t1339 = t1122 * t1348;
t1271 = pkin(1) * t1304 + pkin(2);
t1338 = (t1271 * t1317 - t1311 * t1409 + pkin(3)) * t1123 / (t1271 * t1311 + t1317 * t1409) * t1401;
t1337 = (t1271 * t1319 - t1313 * t1409 + pkin(3)) * t1125 / (t1271 * t1313 + t1319 * t1409) * t1399;
t1336 = (t1271 * t1321 - t1315 * t1409 + pkin(3)) * t1124 / (t1271 * t1315 + t1321 * t1409) * t1400;
t1335 = g(1) * t1421 + g(2) * t1418 + t1220 * t1424 + t1227 * t1423 + g(3) * cos(t1281);
t1334 = g(1) * t1420 + g(2) * t1417 + t1222 * t1424 + t1229 * t1423 + g(3) * cos(t1282);
t1333 = g(1) * t1419 + g(2) * t1416 + t1224 * t1424 + t1231 * t1423 + g(3) * cos(t1283);
t1332 = g(1) * t1418 + g(2) * t1421 - g(3) * t1259 + t1220 * t1423 + t1227 * t1424;
t1331 = g(1) * t1417 + g(2) * t1420 - g(3) * t1260 + t1222 * t1423 + t1229 * t1424;
t1330 = g(1) * t1416 + g(2) * t1419 - g(3) * t1261 + t1224 * t1423 + t1231 * t1424;
t1158 = -g(3) * t1316 + t1180 * t1322;
t1156 = -g(3) * t1314 + t1179 * t1320;
t1154 = -g(3) * t1312 + t1178 * t1318;
t1152 = t1404 * t1292 + t1405 * t1295;
t1151 = t1404 * t1291 + t1405 * t1294;
t1150 = t1404 * t1290 + t1405 * t1293;
t1110 = -t1308 * t1369 + t1433;
t1109 = -t1260 * t1371 + t1432;
t1108 = -t1259 * t1373 + t1431;
t1107 = pkin(1) * t1108 + t1153;
t1106 = pkin(1) * t1110 + t1157;
t1105 = pkin(1) * t1109 + t1155;
t1104 = (-t1110 * t1315 + t1321 * t1346) * pkin(2) + (-t1110 * t1274 + t1280 * t1346) * pkin(1) + t1330;
t1103 = (-t1109 * t1313 + t1319 * t1345) * pkin(2) + (-t1109 * t1273 + t1279 * t1345) * pkin(1) + t1331;
t1102 = (-t1108 * t1311 + t1317 * t1347) * pkin(2) + (-t1108 * t1272 + t1278 * t1347) * pkin(1) + t1332;
t1101 = (t1110 * t1321 + t1315 * t1346) * pkin(2) + (t1110 * t1280 + t1274 * t1346) * pkin(1) + t1333;
t1100 = (t1109 * t1319 + t1313 * t1345) * pkin(2) + (t1109 * t1279 + t1273 * t1345) * pkin(1) + t1334;
t1099 = (t1108 * t1317 + t1311 * t1347) * pkin(2) + (t1108 * t1278 + t1272 * t1347) * pkin(1) + t1335;
t1098 = (-t1261 + t1364) * t1370 - t1336 / 0.2e1 + (-t1355 + (t1384 + t1390) * t1185) * t1439 + t1433;
t1097 = (-t1260 + t1365) * t1371 - t1337 / 0.2e1 + (-t1354 + (t1386 + t1392) * t1183) * t1439 + t1432;
t1096 = (-t1259 + t1366) * t1373 - t1338 / 0.2e1 + (-t1356 + (t1388 + t1394) * t1181) * t1439 + t1431;
t1095 = (-t1261 + t1364 / 0.2e1) * t1370 - t1336 / 0.4e1 + (-t1355 / 0.2e1 + (t1384 / 0.2e1 + t1390 / 0.2e1) * t1185) * t1439 + t1433;
t1094 = (-t1260 + t1365 / 0.2e1) * t1371 - t1337 / 0.4e1 + (-t1354 / 0.2e1 + (t1386 / 0.2e1 + t1392 / 0.2e1) * t1183) * t1439 + t1432;
t1093 = (-t1259 + t1366 / 0.2e1) * t1373 - t1338 / 0.4e1 + (-t1356 / 0.2e1 + (t1388 / 0.2e1 + t1394 / 0.2e1) * t1181) * t1439 + t1431;
t1092 = (t1095 * t1315 + t1321 * t1340) * t1427 + (-t1095 * t1274 - t1280 * t1340) * t1402 + t1330;
t1091 = (t1094 * t1313 + t1319 * t1339) * t1427 + (-t1094 * t1273 - t1279 * t1339) * t1402 + t1331;
t1090 = (t1093 * t1311 + t1317 * t1341) * t1427 + (-t1093 * t1272 - t1278 * t1341) * t1402 + t1332;
t1089 = (t1095 * t1321 - t1315 * t1340) * t1426 + (t1095 * t1280 - t1274 * t1340) * t1402 + t1333;
t1088 = (t1094 * t1319 - t1313 * t1339) * t1426 + (t1094 * t1279 - t1273 * t1339) * t1402 + t1334;
t1087 = (t1093 * t1317 - t1311 * t1341) * t1426 + (t1093 * t1278 - t1272 * t1341) * t1402 + t1335;
t1 = [(t1150 * t1290 + t1151 * t1291 + t1152 * t1292) * MDP(4) + t1404 * MDP(8) + (t1108 * t1377 + t1109 * t1376 + t1110 * t1375) * t1445 + (t1153 * t1377 + t1155 * t1376 + t1157 * t1375) * t1444 + (t1154 * t1377 + t1156 * t1376 + t1158 * t1375) * t1443 + (t1096 * t1377 + t1097 * t1376 + t1098 * t1375) * t1442 + (t1087 * t1377 + t1088 * t1376 + t1089 * t1375) * t1441 + (t1090 * t1377 + t1091 * t1376 + t1092 * t1375) * t1440 + (t1105 * t1376 + t1106 * t1375 + t1107 * t1377) * t1438 + ((t1096 * t1389 + t1097 * t1387 + t1098 * t1385) * MDP(5) + (t1099 * t1389 + t1100 * t1387 + t1101 * t1385) * MDP(6) + (t1102 * t1389 + t1103 * t1387 + t1104 * t1385) * MDP(7)) * t1439; (t1150 * t1293 + t1151 * t1294 + t1152 * t1295) * MDP(4) + t1405 * MDP(8) + (-t1108 * t1380 - t1109 * t1379 - t1110 * t1378) * t1445 + (-t1153 * t1380 - t1155 * t1379 - t1157 * t1378) * t1444 + (-t1154 * t1380 - t1156 * t1379 - t1158 * t1378) * t1443 + (-t1096 * t1380 - t1097 * t1379 - t1098 * t1378) * t1442 + (-t1087 * t1380 - t1088 * t1379 - t1089 * t1378) * t1441 + (-t1090 * t1380 - t1091 * t1379 - t1092 * t1378) * t1440 + (-t1105 * t1379 - t1106 * t1378 - t1107 * t1380) * t1438 + ((t1096 * t1395 + t1097 * t1393 + t1098 * t1391) * MDP(5) + (t1099 * t1395 + t1100 * t1393 + t1101 * t1391) * MDP(6) + (t1102 * t1395 + t1103 * t1393 + t1104 * t1391) * MDP(7)) * t1439; (-t1108 * t1374 - t1109 * t1372 - t1110 * t1369) * MDP(1) + (-t1153 * t1374 - t1155 * t1372 - t1157 * t1369) * MDP(2) + (-t1154 * t1374 - t1156 * t1372 - t1158 * t1369) * MDP(3) + (-t1096 * t1374 - t1097 * t1372 - t1098 * t1369) * MDP(5) + (-t1087 * t1374 - t1088 * t1372 - t1089 * t1369) * MDP(6) + (-t1090 * t1374 - t1091 * t1372 - t1092 * t1369) * MDP(7) + (t1308 - g(3)) * MDP(8) + (-t1105 * t1372 - t1106 * t1369 - t1107 * t1374) * t1422 + ((t1096 * t1383 + t1097 * t1382 + t1098 * t1381) * MDP(5) + (t1099 * t1383 + t1100 * t1382 + t1101 * t1381) * MDP(6) + (t1102 * t1383 + t1103 * t1382 + t1104 * t1381) * MDP(7)) * t1329;];
tauX  = t1;
