% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RPRR1G2P2A0
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
%   see P3RPRR1G2P2A0_convert_par2_MPV_fixb.m

% Output:
% tauX [3x1]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:25
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRR1G2P2A0_invdyn_para_pf_mdp(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G2P2A0_invdyn_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRR1G2P2A0_invdyn_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRR1G2P2A0_invdyn_para_pf_mdp: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G2P2A0_invdyn_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRR1G2P2A0_invdyn_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G2P2A0_invdyn_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G2P2A0_invdyn_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G2P2A0_invdyn_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3RPRR1G2P2A0_invdyn_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:25:08
% EndTime: 2020-03-09 21:25:12
% DurationCPUTime: 3.90s
% Computational Cost: add. (57783->341), mult. (43500->509), div. (2799->10), fcn. (22815->71), ass. (0->240)
t1426 = 2 * xDP(3);
t1425 = MDP(1) / 0.2e1;
t1424 = MDP(2) / 0.2e1;
t1423 = MDP(3) / 0.2e1;
t1422 = MDP(5) / 0.2e1;
t1421 = MDP(6) / 0.2e1;
t1420 = MDP(7) / 0.2e1;
t1321 = 0.1e1 / pkin(3);
t1419 = t1321 / 0.2e1;
t1409 = MDP(4) * pkin(1);
t1418 = t1409 / 0.2e1;
t1296 = cos(pkin(7));
t1417 = pkin(2) * t1296;
t1298 = legFrame(2,2);
t1360 = qJ(1,2) + pkin(7);
t1244 = t1298 + t1360;
t1238 = qJ(3,2) + t1244;
t1208 = sin(t1238);
t1245 = -t1298 + t1360;
t1239 = qJ(3,2) + t1245;
t1209 = sin(t1239);
t1156 = t1208 + t1209;
t1214 = cos(t1238);
t1215 = cos(t1239);
t1159 = -t1215 + t1214;
t1293 = pkin(7) + qJ(3,2);
t1274 = qJ(1,2) + t1293;
t1255 = cos(t1274);
t1318 = xDP(2);
t1319 = xDP(1);
t1127 = t1156 * t1319 + t1159 * t1318 + t1255 * t1426;
t1265 = sin(t1293);
t1305 = sin(qJ(3,2));
t1204 = pkin(1) * t1265 + pkin(2) * t1305;
t1175 = 0.1e1 / t1204;
t1121 = t1127 * t1175 / 0.2e1;
t1268 = cos(t1293);
t1300 = xDDP(3);
t1311 = cos(qJ(3,2));
t1278 = qJ(1,2) + t1298;
t1279 = qJ(1,2) - t1298;
t1129 = (sin(t1278) + sin(t1279)) * pkin(1) + (sin(t1244) + sin(t1245)) * pkin(2) + t1156 * pkin(3);
t1132 = -t1159 * pkin(3) + (-cos(t1244) + cos(t1245)) * pkin(2) + (-cos(t1278) + cos(t1279)) * pkin(1);
t1312 = cos(qJ(1,2));
t1153 = -pkin(1) * t1312 - pkin(2) * cos(t1360) - pkin(3) * t1255;
t1394 = (-t1129 * t1319 + t1132 * t1318 + t1153 * t1426) * t1175 * t1321;
t1345 = t1394 / 0.2e1;
t1302 = xDDP(1);
t1362 = t1302 / 0.2e1;
t1301 = xDDP(2);
t1363 = t1301 / 0.2e1;
t1366 = t1175 * t1255;
t1369 = t1159 * t1175;
t1372 = t1156 * t1175;
t1176 = 0.1e1 / t1204 ^ 2;
t1389 = t1127 * t1176;
t1115 = t1121 + t1345;
t1400 = pkin(3) * t1115;
t1099 = t1175 * t1345 * t1400 + t1300 * t1366 + t1362 * t1372 + t1363 * t1369 - (-t1400 + (-pkin(1) * t1268 - pkin(2) * t1311) * t1121) * t1389 / 0.2e1;
t1297 = legFrame(3,2);
t1359 = qJ(1,3) + pkin(7);
t1242 = t1297 + t1359;
t1236 = qJ(3,3) + t1242;
t1206 = sin(t1236);
t1243 = -t1297 + t1359;
t1237 = qJ(3,3) + t1243;
t1207 = sin(t1237);
t1155 = t1206 + t1207;
t1212 = cos(t1236);
t1213 = cos(t1237);
t1158 = -t1213 + t1212;
t1292 = pkin(7) + qJ(3,3);
t1273 = qJ(1,3) + t1292;
t1254 = cos(t1273);
t1125 = t1155 * t1319 + t1158 * t1318 + t1254 * t1426;
t1264 = sin(t1292);
t1303 = sin(qJ(3,3));
t1203 = pkin(1) * t1264 + pkin(2) * t1303;
t1173 = 0.1e1 / t1203;
t1119 = t1125 * t1173 / 0.2e1;
t1267 = cos(t1292);
t1309 = cos(qJ(3,3));
t1276 = qJ(1,3) + t1297;
t1277 = qJ(1,3) - t1297;
t1128 = (sin(t1276) + sin(t1277)) * pkin(1) + (sin(t1242) + sin(t1243)) * pkin(2) + t1155 * pkin(3);
t1131 = -t1158 * pkin(3) + (-cos(t1242) + cos(t1243)) * pkin(2) + (-cos(t1276) + cos(t1277)) * pkin(1);
t1310 = cos(qJ(1,3));
t1152 = -pkin(1) * t1310 - pkin(2) * cos(t1359) - pkin(3) * t1254;
t1393 = (-t1128 * t1319 + t1131 * t1318 + t1152 * t1426) * t1173 * t1321;
t1344 = t1393 / 0.2e1;
t1367 = t1173 * t1254;
t1370 = t1158 * t1173;
t1373 = t1155 * t1173;
t1174 = 0.1e1 / t1203 ^ 2;
t1391 = t1125 * t1174;
t1113 = t1119 + t1344;
t1401 = pkin(3) * t1113;
t1098 = t1173 * t1344 * t1401 + t1300 * t1367 + t1362 * t1373 + t1363 * t1370 - (-t1401 + (-pkin(1) * t1267 - pkin(2) * t1309) * t1119) * t1391 / 0.2e1;
t1299 = legFrame(1,2);
t1361 = qJ(1,1) + pkin(7);
t1246 = t1299 + t1361;
t1240 = qJ(3,1) + t1246;
t1210 = sin(t1240);
t1247 = -t1299 + t1361;
t1241 = qJ(3,1) + t1247;
t1211 = sin(t1241);
t1157 = t1210 + t1211;
t1216 = cos(t1240);
t1217 = cos(t1241);
t1160 = -t1217 + t1216;
t1294 = pkin(7) + qJ(3,1);
t1275 = qJ(1,1) + t1294;
t1256 = cos(t1275);
t1126 = t1157 * t1319 + t1160 * t1318 + t1256 * t1426;
t1266 = sin(t1294);
t1307 = sin(qJ(3,1));
t1205 = pkin(1) * t1266 + pkin(2) * t1307;
t1177 = 0.1e1 / t1205;
t1120 = t1126 * t1177 / 0.2e1;
t1269 = cos(t1294);
t1313 = cos(qJ(3,1));
t1280 = qJ(1,1) + t1299;
t1281 = qJ(1,1) - t1299;
t1130 = (sin(t1280) + sin(t1281)) * pkin(1) + (sin(t1246) + sin(t1247)) * pkin(2) + t1157 * pkin(3);
t1133 = -t1160 * pkin(3) + (-cos(t1246) + cos(t1247)) * pkin(2) + (-cos(t1280) + cos(t1281)) * pkin(1);
t1314 = cos(qJ(1,1));
t1154 = -pkin(1) * t1314 - pkin(2) * cos(t1361) - pkin(3) * t1256;
t1392 = (-t1130 * t1319 + t1133 * t1318 + t1154 * t1426) * t1177 * t1321;
t1343 = t1392 / 0.2e1;
t1365 = t1177 * t1256;
t1368 = t1160 * t1177;
t1371 = t1157 * t1177;
t1178 = 0.1e1 / t1205 ^ 2;
t1390 = t1126 * t1178;
t1114 = t1120 + t1343;
t1399 = t1114 * pkin(3);
t1100 = t1177 * t1343 * t1399 + t1300 * t1365 + t1362 * t1371 + t1363 * t1368 - (-t1399 + (-pkin(1) * t1269 - pkin(2) * t1313) * t1120) * t1390 / 0.2e1;
t1416 = -0.2e1 * pkin(2);
t1415 = 0.2e1 * pkin(2);
t1413 = -g(1) / 0.2e1;
t1412 = g(1) / 0.2e1;
t1411 = -g(2) / 0.2e1;
t1410 = g(2) / 0.2e1;
t1408 = t1206 / 0.2e1;
t1407 = t1208 / 0.2e1;
t1406 = t1210 / 0.2e1;
t1405 = -t1213 / 0.2e1;
t1404 = -t1215 / 0.2e1;
t1403 = -t1217 / 0.2e1;
t1402 = pkin(1) * sin(pkin(7));
t1398 = t1301 - g(2);
t1397 = t1302 - g(1);
t1396 = pkin(3) * t1415;
t1395 = 0.2e1 * pkin(1);
t1388 = t1128 * t1173;
t1387 = t1128 * t1302;
t1386 = t1129 * t1175;
t1385 = t1129 * t1302;
t1384 = t1130 * t1177;
t1383 = t1130 * t1302;
t1382 = t1131 * t1173;
t1381 = t1131 * t1301;
t1380 = t1132 * t1175;
t1379 = t1132 * t1301;
t1378 = t1133 * t1177;
t1377 = t1133 * t1301;
t1376 = t1152 * t1173;
t1375 = t1153 * t1175;
t1374 = t1154 * t1177;
t1364 = t1300 * t1321;
t1110 = t1119 + t1344 / 0.2e1;
t1291 = pkin(1) ^ 2 + pkin(2) ^ 2;
t1320 = pkin(3) ^ 2;
t1351 = (t1110 * t1309 * t1396 + t1291 * t1119 + t1113 * t1320 + (pkin(3) * t1110 * t1267 + t1119 * t1417) * t1395) * t1391;
t1111 = t1120 + t1343 / 0.2e1;
t1350 = (t1111 * t1313 * t1396 + t1291 * t1120 + t1114 * t1320 + (pkin(3) * t1111 * t1269 + t1120 * t1417) * t1395) * t1390;
t1112 = t1121 + t1345 / 0.2e1;
t1349 = (t1112 * t1311 * t1396 + t1291 * t1121 + t1115 * t1320 + (pkin(3) * t1112 * t1268 + t1121 * t1417) * t1395) * t1389;
t1342 = t1125 ^ 2 * t1174 / 0.4e1;
t1341 = t1126 ^ 2 * t1178 / 0.4e1;
t1340 = t1127 ^ 2 * t1176 / 0.4e1;
t1282 = sin(t1297);
t1285 = cos(t1297);
t1170 = g(1) * t1285 - g(2) * t1282;
t1304 = sin(qJ(1,3));
t1146 = g(3) * t1304 - t1170 * t1310;
t1283 = sin(t1298);
t1286 = cos(t1298);
t1171 = g(1) * t1286 - g(2) * t1283;
t1306 = sin(qJ(1,2));
t1147 = g(3) * t1306 - t1171 * t1312;
t1284 = sin(t1299);
t1287 = cos(t1299);
t1172 = g(1) * t1287 - g(2) * t1284;
t1308 = sin(qJ(1,1));
t1148 = g(3) * t1308 - t1172 * t1314;
t1339 = t1110 * t1344;
t1338 = t1111 * t1343;
t1337 = t1112 * t1345;
t1336 = t1364 * t1376;
t1335 = t1364 * t1375;
t1334 = t1364 * t1374;
t1263 = pkin(1) * t1296 + pkin(2);
t1330 = t1113 * (t1263 * t1309 - t1303 * t1402 + pkin(3)) / (t1263 * t1303 + t1309 * t1402) * t1393;
t1329 = t1114 * (t1263 * t1313 - t1307 * t1402 + pkin(3)) / (t1263 * t1307 + t1313 * t1402) * t1392;
t1328 = t1115 * (t1263 * t1311 - t1305 * t1402 + pkin(3)) / (t1263 * t1305 + t1311 * t1402) * t1394;
t1327 = g(1) * t1405 + g(2) * t1408 + t1207 * t1411 + t1212 * t1413 + g(3) * sin(t1273);
t1326 = g(1) * t1404 + g(2) * t1407 + t1209 * t1411 + t1214 * t1413 + g(3) * sin(t1274);
t1325 = g(1) * t1403 + g(2) * t1406 + t1211 * t1411 + t1216 * t1413 + g(3) * sin(t1275);
t1324 = g(1) * t1408 + g(2) * t1405 + g(3) * t1254 + t1207 * t1412 + t1212 * t1410;
t1323 = g(1) * t1407 + g(2) * t1404 + g(3) * t1255 + t1209 * t1412 + t1214 * t1410;
t1322 = g(1) * t1406 + g(2) * t1403 + g(3) * t1256 + t1211 * t1412 + t1216 * t1410;
t1151 = g(3) * t1314 + t1172 * t1308;
t1150 = g(3) * t1312 + t1171 * t1306;
t1149 = g(3) * t1310 + t1170 * t1304;
t1142 = t1397 * t1284 + t1398 * t1287;
t1141 = t1397 * t1283 + t1398 * t1286;
t1140 = t1397 * t1282 + t1398 * t1285;
t1097 = t1100 * pkin(1) + t1148;
t1096 = pkin(1) * t1099 + t1147;
t1095 = pkin(1) * t1098 + t1146;
t1094 = (t1099 * t1311 + t1305 * t1340) * pkin(2) + (t1099 * t1268 + t1265 * t1340) * pkin(1) + t1326;
t1093 = (t1098 * t1309 + t1303 * t1342) * pkin(2) + (t1098 * t1267 + t1264 * t1342) * pkin(1) + t1327;
t1092 = (-t1100 * t1307 + t1313 * t1341) * pkin(2) + (-t1100 * t1266 + t1269 * t1341) * pkin(1) + t1322;
t1091 = (-t1099 * t1305 + t1311 * t1340) * pkin(2) + (-t1099 * t1265 + t1268 * t1340) * pkin(1) + t1323;
t1090 = (-t1098 * t1303 + t1309 * t1342) * pkin(2) + (-t1098 * t1264 + t1267 * t1342) * pkin(1) + t1324;
t1089 = (t1100 * t1313 + t1307 * t1341) * pkin(2) + (t1100 * t1269 + t1266 * t1341) * pkin(1) + t1325;
t1088 = t1334 - t1329 / 0.2e1 + (-t1350 + (t1377 - t1383) * t1177) * t1419 + t1100;
t1087 = t1335 - t1328 / 0.2e1 + (-t1349 + (t1379 - t1385) * t1175) * t1419 + t1099;
t1086 = t1336 - t1330 / 0.2e1 + (-t1351 + (t1381 - t1387) * t1173) * t1419 + t1098;
t1085 = t1334 / 0.2e1 - t1329 / 0.4e1 + (-t1350 / 0.2e1 + (-t1383 / 0.2e1 + t1377 / 0.2e1) * t1177) * t1419 + t1100;
t1084 = t1335 / 0.2e1 - t1328 / 0.4e1 + (-t1349 / 0.2e1 + (-t1385 / 0.2e1 + t1379 / 0.2e1) * t1175) * t1419 + t1099;
t1083 = t1336 / 0.2e1 - t1330 / 0.4e1 + (-t1351 / 0.2e1 + (-t1387 / 0.2e1 + t1381 / 0.2e1) * t1173) * t1419 + t1098;
t1082 = (t1085 * t1313 - t1307 * t1338) * t1415 + (t1085 * t1269 - t1266 * t1338) * t1395 + t1325;
t1081 = (t1084 * t1311 - t1305 * t1337) * t1415 + (t1084 * t1268 - t1265 * t1337) * t1395 + t1326;
t1080 = (t1083 * t1309 - t1303 * t1339) * t1415 + (t1083 * t1267 - t1264 * t1339) * t1395 + t1327;
t1079 = (t1085 * t1307 + t1313 * t1338) * t1416 + (-t1085 * t1266 - t1269 * t1338) * t1395 + t1322;
t1078 = (t1084 * t1305 + t1311 * t1337) * t1416 + (-t1084 * t1265 - t1268 * t1337) * t1395 + t1323;
t1077 = (t1083 * t1303 + t1309 * t1339) * t1416 + (-t1083 * t1264 - t1267 * t1339) * t1395 + t1324;
t1 = [(t1140 * t1282 + t1141 * t1283 + t1142 * t1284) * MDP(4) + t1397 * MDP(8) + (t1098 * t1373 + t1099 * t1372 + t1100 * t1371) * t1425 + (t1146 * t1373 + t1147 * t1372 + t1148 * t1371) * t1424 + (t1149 * t1373 + t1150 * t1372 + t1151 * t1371) * t1423 + (t1086 * t1373 + t1087 * t1372 + t1088 * t1371) * t1422 + (t1080 * t1373 + t1081 * t1372 + t1082 * t1371) * t1421 + (t1077 * t1373 + t1078 * t1372 + t1079 * t1371) * t1420 + (t1095 * t1373 + t1096 * t1372 + t1097 * t1371) * t1418 + ((-t1086 * t1388 - t1087 * t1386 - t1088 * t1384) * MDP(5) + (-t1089 * t1384 - t1093 * t1388 - t1094 * t1386) * MDP(6) + (-t1090 * t1388 - t1091 * t1386 - t1092 * t1384) * MDP(7)) * t1419; (t1140 * t1285 + t1141 * t1286 + t1142 * t1287) * MDP(4) + t1398 * MDP(8) + (t1098 * t1370 + t1099 * t1369 + t1100 * t1368) * t1425 + (t1146 * t1370 + t1147 * t1369 + t1148 * t1368) * t1424 + (t1149 * t1370 + t1150 * t1369 + t1151 * t1368) * t1423 + (t1086 * t1370 + t1087 * t1369 + t1088 * t1368) * t1422 + (t1080 * t1370 + t1081 * t1369 + t1082 * t1368) * t1421 + (t1077 * t1370 + t1078 * t1369 + t1079 * t1368) * t1420 + (t1095 * t1370 + t1096 * t1369 + t1097 * t1368) * t1418 + ((t1086 * t1382 + t1087 * t1380 + t1088 * t1378) * MDP(5) + (t1089 * t1378 + t1093 * t1382 + t1094 * t1380) * MDP(6) + (t1090 * t1382 + t1091 * t1380 + t1092 * t1378) * MDP(7)) * t1419; (t1098 * t1367 + t1099 * t1366 + t1100 * t1365) * MDP(1) + (t1146 * t1367 + t1147 * t1366 + t1148 * t1365) * MDP(2) + (t1149 * t1367 + t1150 * t1366 + t1151 * t1365) * MDP(3) + (t1086 * t1367 + t1087 * t1366 + t1088 * t1365) * MDP(5) + (t1080 * t1367 + t1081 * t1366 + t1082 * t1365) * MDP(6) + (t1077 * t1367 + t1078 * t1366 + t1079 * t1365) * MDP(7) + (t1300 - g(3)) * MDP(8) + (t1095 * t1367 + t1096 * t1366 + t1097 * t1365) * t1409 + ((t1086 * t1376 + t1087 * t1375 + t1088 * t1374) * MDP(5) + (t1089 * t1374 + t1093 * t1376 + t1094 * t1375) * MDP(6) + (t1090 * t1376 + t1091 * t1375 + t1092 * t1374) * MDP(7)) * t1321;];
tauX  = t1;
