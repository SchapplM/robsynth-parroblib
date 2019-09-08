% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RRR1G1P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [2x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
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
%   see P3RRR1G1P1A0_convert_par2_MPV_fixb.m

% Output:
% MMX [3x3]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 15:38
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P3RRR1G1P1A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(5,1),zeros(10,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRR1G1P1A0_inertia_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RRR1G1P1A0_inertia_para_pf_mdp: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRR1G1P1A0_inertia_para_pf_mdp: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRR1G1P1A0_inertia_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRR1G1P1A0_inertia_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'P3RRR1G1P1A0_inertia_para_pf_mdp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 15:38:41
% EndTime: 2019-05-03 15:38:43
% DurationCPUTime: 1.90s
% Computational Cost: add. (3231->212), mult. (3673->417), div. (759->9), fcn. (4484->26), ass. (0->198)
t1315 = qJ(1,3) + qJ(2,3);
t1301 = sin(t1315);
t1304 = cos(t1315);
t1322 = sin(qJ(1,3));
t1328 = cos(qJ(1,3));
t1282 = t1301 * t1328 - t1304 * t1322;
t1434 = 0.1e1 / t1282;
t1316 = qJ(1,2) + qJ(2,2);
t1302 = sin(t1316);
t1305 = cos(t1316);
t1324 = sin(qJ(1,2));
t1330 = cos(qJ(1,2));
t1283 = t1302 * t1330 - t1305 * t1324;
t1433 = 0.1e1 / t1283;
t1317 = qJ(1,1) + qJ(2,1);
t1303 = sin(t1317);
t1306 = cos(t1317);
t1326 = sin(qJ(1,1));
t1332 = cos(qJ(1,1));
t1284 = t1303 * t1332 - t1306 * t1326;
t1432 = 0.1e1 / t1284;
t1431 = MDP(1) / pkin(1) ^ 2;
t1333 = xP(3);
t1313 = sin(t1333);
t1314 = cos(t1333);
t1334 = koppelP(3,2);
t1337 = koppelP(3,1);
t1285 = t1313 * t1337 + t1314 * t1334;
t1288 = -t1313 * t1334 + t1314 * t1337;
t1318 = legFrame(3,3);
t1307 = sin(t1318);
t1310 = cos(t1318);
t1234 = (t1285 * t1310 - t1288 * t1307) * t1304 - (t1285 * t1307 + t1288 * t1310) * t1301;
t1292 = t1322 * t1337 - t1328 * t1334;
t1293 = t1322 * t1334 + t1328 * t1337;
t1216 = pkin(1) * ((t1292 * t1314 - t1293 * t1313) * t1310 + t1307 * (t1292 * t1313 + t1293 * t1314)) - t1234 * pkin(2);
t1430 = t1216 * t1434;
t1335 = koppelP(2,2);
t1338 = koppelP(2,1);
t1286 = t1313 * t1338 + t1314 * t1335;
t1289 = -t1313 * t1335 + t1314 * t1338;
t1319 = legFrame(2,3);
t1308 = sin(t1319);
t1311 = cos(t1319);
t1235 = (t1286 * t1311 - t1289 * t1308) * t1305 - (t1286 * t1308 + t1289 * t1311) * t1302;
t1294 = t1324 * t1338 - t1330 * t1335;
t1295 = t1324 * t1335 + t1330 * t1338;
t1217 = pkin(1) * ((t1294 * t1314 - t1295 * t1313) * t1311 + t1308 * (t1294 * t1313 + t1295 * t1314)) - t1235 * pkin(2);
t1429 = t1217 * t1433;
t1336 = koppelP(1,2);
t1339 = koppelP(1,1);
t1287 = t1313 * t1339 + t1314 * t1336;
t1290 = -t1313 * t1336 + t1314 * t1339;
t1320 = legFrame(1,3);
t1309 = sin(t1320);
t1312 = cos(t1320);
t1236 = (t1287 * t1312 - t1290 * t1309) * t1306 - (t1287 * t1309 + t1290 * t1312) * t1303;
t1296 = t1326 * t1339 - t1332 * t1336;
t1297 = t1326 * t1336 + t1332 * t1339;
t1218 = pkin(1) * ((t1296 * t1314 - t1297 * t1313) * t1312 + t1309 * (t1296 * t1313 + t1297 * t1314)) - t1236 * pkin(2);
t1428 = t1218 * t1432;
t1268 = -t1301 * t1307 + t1304 * t1310;
t1258 = pkin(1) * (-t1307 * t1322 + t1310 * t1328) + t1268 * pkin(2);
t1340 = 0.1e1 / pkin(2);
t1341 = 0.1e1 / pkin(1);
t1400 = t1340 * t1341;
t1363 = t1434 * t1400;
t1357 = t1258 * t1363;
t1411 = t1434 * t1341;
t1378 = t1268 * t1411;
t1229 = -t1357 + t1378;
t1427 = t1229 * t1434;
t1270 = -t1302 * t1308 + t1305 * t1311;
t1259 = pkin(1) * (-t1308 * t1324 + t1311 * t1330) + t1270 * pkin(2);
t1362 = t1433 * t1400;
t1356 = t1259 * t1362;
t1409 = t1433 * t1341;
t1374 = t1270 * t1409;
t1231 = -t1356 + t1374;
t1426 = t1231 * t1433;
t1272 = -t1303 * t1309 + t1306 * t1312;
t1260 = pkin(1) * (-t1309 * t1326 + t1312 * t1332) + t1272 * pkin(2);
t1361 = t1432 * t1400;
t1355 = t1260 * t1361;
t1407 = t1432 * t1341;
t1370 = t1272 * t1407;
t1233 = -t1355 + t1370;
t1425 = t1233 * t1432;
t1424 = t1234 * t1434;
t1423 = t1235 * t1433;
t1422 = t1236 * t1432;
t1267 = t1301 * t1310 + t1304 * t1307;
t1255 = pkin(1) * (t1307 * t1328 + t1310 * t1322) + t1267 * pkin(2);
t1421 = t1255 * t1434;
t1269 = t1302 * t1311 + t1305 * t1308;
t1256 = pkin(1) * (t1308 * t1330 + t1311 * t1324) + t1269 * pkin(2);
t1420 = t1256 * t1433;
t1271 = t1303 * t1312 + t1306 * t1309;
t1257 = pkin(1) * (t1309 * t1332 + t1312 * t1326) + t1271 * pkin(2);
t1419 = t1257 * t1432;
t1418 = t1267 * t1434;
t1417 = t1268 * t1434;
t1416 = t1269 * t1433;
t1415 = t1270 * t1433;
t1414 = t1271 * t1432;
t1413 = t1272 * t1432;
t1412 = t1434 ^ 2;
t1410 = t1433 ^ 2;
t1408 = t1432 ^ 2;
t1321 = sin(qJ(2,3));
t1406 = t1434 * t1321;
t1327 = cos(qJ(2,3));
t1405 = t1434 * t1327;
t1323 = sin(qJ(2,2));
t1404 = t1433 * t1323;
t1329 = cos(qJ(2,2));
t1403 = t1433 * t1329;
t1325 = sin(qJ(2,1));
t1402 = t1432 * t1325;
t1331 = cos(qJ(2,1));
t1401 = t1432 * t1331;
t1360 = t1255 * t1363;
t1399 = t1285 * t1357 - t1288 * t1360;
t1359 = t1256 * t1362;
t1398 = t1286 * t1356 - t1289 * t1359;
t1358 = t1257 * t1361;
t1397 = t1287 * t1355 - t1290 * t1358;
t1249 = t1285 * t1378;
t1379 = t1267 * t1411;
t1252 = t1288 * t1379;
t1219 = -t1249 + t1252;
t1250 = t1286 * t1374;
t1375 = t1269 * t1409;
t1253 = t1289 * t1375;
t1220 = -t1250 + t1253;
t1251 = t1287 * t1370;
t1371 = t1271 * t1407;
t1254 = t1290 * t1371;
t1221 = -t1251 + t1254;
t1396 = t1219 * t1430;
t1395 = t1220 * t1429;
t1394 = t1221 * t1428;
t1225 = -t1357 + 0.2e1 * t1378;
t1393 = t1225 * t1417;
t1226 = -t1356 + 0.2e1 * t1374;
t1392 = t1226 * t1415;
t1227 = -t1355 + 0.2e1 * t1370;
t1391 = t1227 * t1413;
t1390 = t1234 * t1412;
t1389 = t1234 * t1406;
t1388 = t1234 * t1405;
t1387 = t1235 * t1410;
t1386 = t1235 * t1404;
t1385 = t1235 * t1403;
t1384 = t1236 * t1408;
t1383 = t1236 * t1402;
t1382 = t1236 * t1401;
t1381 = t1321 * t1418;
t1380 = t1327 * t1418;
t1377 = t1323 * t1416;
t1376 = t1329 * t1416;
t1373 = t1325 * t1414;
t1372 = t1331 * t1414;
t1369 = t1434 * t1406;
t1368 = t1434 * t1405;
t1367 = t1433 * t1404;
t1366 = t1433 * t1403;
t1365 = t1432 * t1402;
t1364 = t1432 * t1401;
t1354 = t1267 * t1369;
t1353 = t1267 * t1368;
t1352 = t1268 * t1369;
t1351 = t1268 * t1368;
t1350 = t1269 * t1367;
t1349 = t1269 * t1366;
t1348 = t1270 * t1367;
t1347 = t1270 * t1366;
t1346 = t1271 * t1365;
t1345 = t1271 * t1364;
t1344 = t1272 * t1365;
t1343 = t1272 * t1364;
t1291 = (t1313 ^ 2 + t1314 ^ 2) * MDP(10);
t1278 = 0.1e1 / t1284 ^ 2;
t1276 = 0.1e1 / t1283 ^ 2;
t1274 = 0.1e1 / t1282 ^ 2;
t1232 = -t1358 + t1371;
t1230 = -t1359 + t1375;
t1228 = -t1360 + t1379;
t1224 = -t1358 + 0.2e1 * t1371;
t1223 = -t1359 + 0.2e1 * t1375;
t1222 = -t1360 + 0.2e1 * t1379;
t1215 = t1221 + t1397;
t1214 = t1220 + t1398;
t1213 = t1219 + t1399;
t1212 = -0.2e1 * t1251 + 0.2e1 * t1254 + t1397;
t1211 = -0.2e1 * t1250 + 0.2e1 * t1253 + t1398;
t1210 = -0.2e1 * t1249 + 0.2e1 * t1252 + t1399;
t1 = [(t1327 * t1393 + t1329 * t1392 + t1331 * t1391) * MDP(5) + (-t1321 * t1393 - t1323 * t1392 - t1325 * t1391) * MDP(6) + t1291 + (t1268 ^ 2 * t1412 + t1270 ^ 2 * t1410 + t1272 ^ 2 * t1408) * t1431 + ((t1229 * t1417 + t1231 * t1415 + t1233 * t1413) * MDP(4) + ((-t1258 * t1427 - t1259 * t1426 - t1260 * t1425) * MDP(4) + (-t1258 * t1351 - t1259 * t1347 - t1260 * t1343) * MDP(5) + (t1258 * t1352 + t1259 * t1348 + t1260 * t1344) * MDP(6)) * t1340) * t1341; (t1225 * t1380 + t1226 * t1376 + t1227 * t1372) * MDP(5) + (-t1225 * t1381 - t1226 * t1377 - t1227 * t1373) * MDP(6) + (t1267 * t1268 * t1274 + t1269 * t1270 * t1276 + t1271 * t1272 * t1278) * t1431 + ((t1229 * t1418 + t1231 * t1416 + t1233 * t1414) * MDP(4) + ((-t1229 * t1421 - t1231 * t1420 - t1233 * t1419) * MDP(4) + (-t1255 * t1351 - t1256 * t1347 - t1257 * t1343) * MDP(5) + (t1255 * t1352 + t1256 * t1348 + t1257 * t1344) * MDP(6)) * t1340) * t1341; (t1222 * t1380 + t1223 * t1376 + t1224 * t1372) * MDP(5) + (-t1222 * t1381 - t1223 * t1377 - t1224 * t1373) * MDP(6) + t1291 + (t1267 ^ 2 * t1274 + t1269 ^ 2 * t1276 + t1271 ^ 2 * t1278) * t1431 + ((t1228 * t1418 + t1230 * t1416 + t1232 * t1414) * MDP(4) + ((-t1228 * t1421 - t1230 * t1420 - t1232 * t1419) * MDP(4) + (-t1255 * t1353 - t1256 * t1349 - t1257 * t1345) * MDP(5) + (t1255 * t1354 + t1256 * t1350 + t1257 * t1346) * MDP(6)) * t1340) * t1341; (-t1225 * t1388 - t1226 * t1385 - t1227 * t1382) * MDP(5) + (t1225 * t1389 + t1226 * t1386 + t1227 * t1383) * MDP(6) - t1313 * MDP(8) - t1314 * MDP(9) + (-t1268 * t1390 - t1270 * t1387 - t1272 * t1384) * t1431 + ((-t1229 * t1424 - t1231 * t1423 - t1233 * t1422) * MDP(4) + ((-t1216 * t1427 - t1217 * t1426 - t1218 * t1425) * MDP(4) + (-t1216 * t1351 - t1217 * t1347 - t1218 * t1343) * MDP(5) + (t1216 * t1352 + t1217 * t1348 + t1218 * t1344) * MDP(6)) * t1340) * t1341; (-t1222 * t1388 - t1223 * t1385 - t1224 * t1382) * MDP(5) + (t1222 * t1389 + t1223 * t1386 + t1224 * t1383) * MDP(6) + t1314 * MDP(8) - t1313 * MDP(9) + (-t1267 * t1390 - t1269 * t1387 - t1271 * t1384) * t1431 + ((-t1228 * t1424 - t1230 * t1423 - t1232 * t1422) * MDP(4) + ((-t1228 * t1430 - t1230 * t1429 - t1232 * t1428) * MDP(4) + (-t1216 * t1353 - t1217 * t1349 - t1218 * t1345) * MDP(5) + (t1216 * t1354 + t1217 * t1350 + t1218 * t1346) * MDP(6)) * t1340) * t1341; (-t1210 * t1388 - t1211 * t1385 - t1212 * t1382) * MDP(5) + (t1210 * t1389 + t1211 * t1386 + t1212 * t1383) * MDP(6) + MDP(7) + ((-t1327 * t1396 - t1329 * t1395 - t1331 * t1394) * MDP(5) + (t1321 * t1396 + t1323 * t1395 + t1325 * t1394) * MDP(6)) * t1340 + ((-t1219 * t1424 - t1220 * t1423 - t1221 * t1422) * MDP(1) + (-t1213 * t1424 - t1214 * t1423 - t1215 * t1422 + (-t1213 * t1430 - t1214 * t1429 - t1215 * t1428) * t1340) * MDP(4)) * t1341;];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t1(1), t1(2), t1(4); t1(2), t1(3), t1(5); t1(4), t1(5), t1(6);];
MMX  = res;
