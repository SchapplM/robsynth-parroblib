% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRRR8V1G3A0
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
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x12]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 17:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRRR8V1G3A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_regmin: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:16:46
% EndTime: 2020-08-06 17:16:47
% DurationCPUTime: 1.07s
% Computational Cost: add. (762->161), mult. (2088->334), div. (126->7), fcn. (2336->22), ass. (0->169)
t1371 = cos(qJ(2,3));
t1365 = sin(qJ(2,3));
t1370 = cos(qJ(3,3));
t1428 = t1365 * t1370;
t1335 = pkin(2) * t1428 - t1371 * pkin(5);
t1358 = sin(pkin(3));
t1360 = cos(pkin(3));
t1364 = sin(qJ(3,3));
t1440 = t1360 * t1364;
t1309 = pkin(2) * t1440 + t1335 * t1358;
t1461 = 0.1e1 / t1309;
t1373 = cos(qJ(2,2));
t1367 = sin(qJ(2,2));
t1372 = cos(qJ(3,2));
t1425 = t1367 * t1372;
t1336 = pkin(2) * t1425 - t1373 * pkin(5);
t1366 = sin(qJ(3,2));
t1438 = t1360 * t1366;
t1310 = pkin(2) * t1438 + t1336 * t1358;
t1460 = 0.1e1 / t1310;
t1375 = cos(qJ(2,1));
t1369 = sin(qJ(2,1));
t1374 = cos(qJ(3,1));
t1421 = t1369 * t1374;
t1337 = pkin(2) * t1421 - t1375 * pkin(5);
t1368 = sin(qJ(3,1));
t1436 = t1360 * t1368;
t1311 = pkin(2) * t1436 + t1337 * t1358;
t1459 = 0.1e1 / t1311;
t1458 = pkin(2) * t1370;
t1457 = pkin(2) * t1372;
t1456 = pkin(2) * t1374;
t1357 = sin(pkin(6));
t1455 = t1357 * g(3);
t1442 = t1357 * t1360;
t1327 = -t1358 * g(1) - g(2) * t1442;
t1328 = g(1) * t1442 - t1358 * g(2);
t1361 = legFrame(3,2);
t1348 = sin(t1361);
t1351 = cos(t1361);
t1332 = t1351 * g(1) - t1348 * g(2);
t1359 = cos(pkin(6));
t1344 = t1360 * t1359 * g(3);
t1278 = (t1327 * t1348 + t1328 * t1351 + t1344) * t1371 + t1365 * (t1332 * t1359 - t1455);
t1454 = t1278 * t1461;
t1362 = legFrame(2,2);
t1349 = sin(t1362);
t1352 = cos(t1362);
t1333 = t1352 * g(1) - t1349 * g(2);
t1279 = (t1327 * t1349 + t1328 * t1352 + t1344) * t1373 + t1367 * (t1333 * t1359 - t1455);
t1453 = t1279 * t1460;
t1363 = legFrame(1,2);
t1350 = sin(t1363);
t1353 = cos(t1363);
t1334 = t1353 * g(1) - t1350 * g(2);
t1280 = (t1327 * t1350 + t1328 * t1353 + t1344) * t1375 + t1369 * (t1334 * t1359 - t1455);
t1452 = t1280 * t1459;
t1451 = t1461 / t1370;
t1450 = t1460 / t1372;
t1449 = t1459 / t1374;
t1329 = t1348 * g(1) + t1351 * g(2);
t1448 = t1461 * t1329;
t1330 = t1349 * g(1) + t1352 * g(2);
t1447 = t1460 * t1330;
t1331 = t1350 * g(1) + t1353 * g(2);
t1446 = t1459 * t1331;
t1445 = t1329 * t1358;
t1444 = t1330 * t1358;
t1443 = t1331 * t1358;
t1441 = t1358 * t1366;
t1439 = t1360 * t1365;
t1437 = t1360 * t1367;
t1435 = t1360 * t1369;
t1434 = t1360 * t1371;
t1433 = t1360 * t1373;
t1432 = t1360 * t1375;
t1431 = t1364 * t1358;
t1430 = t1364 * t1365;
t1429 = t1364 * t1371;
t1427 = t1366 * t1367;
t1426 = t1366 * t1373;
t1424 = t1368 * t1358;
t1423 = t1368 * t1369;
t1422 = t1368 * t1375;
t1420 = t1370 * t1358;
t1419 = t1370 * t1371;
t1418 = t1372 * t1373;
t1417 = t1374 * t1358;
t1416 = t1374 * t1375;
t1324 = t1360 * t1430 + t1420;
t1300 = t1324 * t1359 + t1357 * t1429;
t1415 = t1300 * t1454;
t1325 = t1358 * t1372 + t1360 * t1427;
t1301 = t1325 * t1359 + t1357 * t1426;
t1414 = t1301 * t1453;
t1326 = t1360 * t1423 + t1417;
t1302 = t1326 * t1359 + t1357 * t1422;
t1413 = t1302 * t1452;
t1319 = t1357 * t1434 + t1359 * t1365;
t1322 = -t1357 * t1439 + t1359 * t1371;
t1412 = (-t1322 * pkin(5) + t1319 * t1458) * t1451;
t1320 = t1357 * t1433 + t1359 * t1367;
t1323 = -t1357 * t1437 + t1359 * t1373;
t1411 = (-t1323 * pkin(5) + t1320 * t1457) * t1450;
t1312 = t1357 * t1435 - t1359 * t1375;
t1321 = t1357 * t1432 + t1359 * t1369;
t1410 = (t1312 * pkin(5) + t1321 * t1456) * t1449;
t1379 = -t1324 * t1357 + t1359 * t1429;
t1409 = t1379 * t1451;
t1378 = -t1325 * t1357 + t1359 * t1426;
t1408 = t1378 * t1450;
t1377 = -t1326 * t1357 + t1359 * t1422;
t1407 = t1377 * t1449;
t1406 = t1348 * t1451;
t1405 = t1351 * t1451;
t1404 = t1349 * t1450;
t1403 = t1352 * t1450;
t1402 = t1350 * t1449;
t1401 = t1353 * t1449;
t1400 = t1278 * t1364 * t1451;
t1399 = t1279 * t1366 * t1450;
t1398 = t1280 * t1368 * t1449;
t1313 = -t1357 * t1365 + t1359 * t1434;
t1316 = t1357 * t1371 + t1359 * t1439;
t1291 = pkin(5) * t1316 + t1313 * t1458;
t1397 = t1291 * t1405;
t1396 = t1291 * t1406;
t1314 = -t1357 * t1367 + t1359 * t1433;
t1317 = t1357 * t1373 + t1359 * t1437;
t1292 = pkin(5) * t1317 + t1314 * t1457;
t1395 = t1292 * t1404;
t1394 = t1292 * t1403;
t1315 = -t1357 * t1369 + t1359 * t1432;
t1318 = t1357 * t1375 + t1359 * t1435;
t1293 = pkin(5) * t1318 + t1315 * t1456;
t1393 = t1293 * t1402;
t1392 = t1293 * t1401;
t1391 = t1300 * t1406;
t1390 = t1300 * t1405;
t1389 = t1301 * t1404;
t1388 = t1301 * t1403;
t1387 = t1302 * t1402;
t1386 = t1302 * t1401;
t1385 = t1300 * t1400;
t1384 = t1301 * t1399;
t1383 = t1302 * t1398;
t1382 = pkin(2) * t1431 - t1335 * t1360;
t1381 = pkin(2) * t1441 - t1336 * t1360;
t1380 = pkin(2) * t1424 - t1337 * t1360;
t1376 = 0.1e1 / pkin(2);
t1340 = pkin(2) * t1416 + pkin(5) * t1369;
t1339 = pkin(2) * t1418 + pkin(5) * t1367;
t1338 = pkin(2) * t1419 + pkin(5) * t1365;
t1289 = t1359 * t1340 + t1380 * t1357;
t1288 = t1359 * t1339 + t1381 * t1357;
t1287 = t1359 * t1338 + t1382 * t1357;
t1286 = -g(3) * t1318 - t1334 * t1312 + t1369 * t1443;
t1285 = g(3) * t1315 + t1334 * t1321 - t1375 * t1443;
t1284 = -g(3) * t1317 + t1333 * t1323 + t1367 * t1444;
t1283 = g(3) * t1314 + t1333 * t1320 - t1373 * t1444;
t1282 = -g(3) * t1316 + t1332 * t1322 + t1365 * t1445;
t1281 = g(3) * t1313 + t1332 * t1319 - t1371 * t1445;
t1277 = t1377 * t1334 - t1302 * g(3) - t1331 * (-t1358 * t1423 + t1360 * t1374);
t1276 = t1378 * t1333 - t1301 * g(3) - t1330 * (-t1358 * t1427 + t1360 * t1372);
t1275 = t1379 * t1332 - t1300 * g(3) - t1329 * (-t1358 * t1430 + t1360 * t1370);
t1274 = ((-t1360 * t1421 + t1424) * t1357 + t1359 * t1416) * t1334 + g(3) * (-t1318 * t1374 + t1359 * t1424) + t1331 * (t1369 * t1417 + t1436);
t1273 = ((-t1360 * t1425 + t1441) * t1357 + t1359 * t1418) * t1333 + g(3) * (-t1317 * t1372 + t1359 * t1441) + t1330 * (t1358 * t1425 + t1438);
t1272 = ((-t1360 * t1428 + t1431) * t1357 + t1359 * t1419) * t1332 + g(3) * (-t1316 * t1370 + t1359 * t1431) + t1329 * (t1365 * t1420 + t1440);
t1 = [-(t1289 * t1353 + t1350 * t1311) * t1446 - (t1288 * t1352 + t1349 * t1310) * t1447 - (t1287 * t1351 + t1348 * t1309) * t1448, 0, -t1281 * t1390 - t1283 * t1388 - t1285 * t1386, -t1282 * t1390 - t1284 * t1388 - t1286 * t1386, 0, 0, 0, 0, 0, -t1351 * t1415 - t1352 * t1414 - t1353 * t1413 + (-t1275 * t1397 - t1276 * t1394 - t1277 * t1392) * t1376, t1351 * t1385 + t1352 * t1384 + t1353 * t1383 + (-t1272 * t1397 - t1273 * t1394 - t1274 * t1392) * t1376, -g(1); -(-t1289 * t1350 + t1353 * t1311) * t1446 - (-t1288 * t1349 + t1352 * t1310) * t1447 - (-t1287 * t1348 + t1351 * t1309) * t1448, 0, t1281 * t1391 + t1283 * t1389 + t1285 * t1387, t1282 * t1391 + t1284 * t1389 + t1286 * t1387, 0, 0, 0, 0, 0, t1348 * t1415 + t1349 * t1414 + t1350 * t1413 + (t1275 * t1396 + t1276 * t1395 + t1277 * t1393) * t1376, -t1348 * t1385 - t1349 * t1384 - t1350 * t1383 + (t1272 * t1396 + t1273 * t1395 + t1274 * t1393) * t1376, -g(2); -(-t1357 * t1340 + t1380 * t1359) * t1446 - (-t1357 * t1339 + t1381 * t1359) * t1447 - (-t1357 * t1338 + t1382 * t1359) * t1448, 0, -t1281 * t1409 - t1283 * t1408 - t1285 * t1407, -t1282 * t1409 - t1284 * t1408 - t1286 * t1407, 0, 0, 0, 0, 0, -t1379 * t1454 - t1378 * t1453 - t1377 * t1452 + (t1275 * t1412 + t1276 * t1411 + t1277 * t1410) * t1376, t1379 * t1400 + t1378 * t1399 + t1377 * t1398 + (t1272 * t1412 + t1273 * t1411 + t1274 * t1410) * t1376, -g(3);];
tau_reg  = t1;
