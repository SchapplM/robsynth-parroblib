% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RPRRR12V1G3A0
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
%   see P3RPRRR12V1G3A0_convert_par2_MPV_fixb.m

% Output:
% tauX [3x1]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:28
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRRR12V1G3A0_invdyn_para_pf_mdp(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(14,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_mdp: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_mdp: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'P3RPRRR12V1G3A0_invdyn_para_pf_mdp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:28:28
% EndTime: 2020-08-06 18:28:33
% DurationCPUTime: 4.41s
% Computational Cost: add. (15597->344), mult. (25556->618), div. (3228->14), fcn. (19026->18), ass. (0->255)
t1359 = sin(qJ(3,3));
t1344 = 0.1e1 / t1359 ^ 2;
t1353 = legFrame(3,2);
t1323 = sin(t1353);
t1326 = cos(t1353);
t1372 = xDP(2);
t1373 = xDP(1);
t1295 = t1323 * t1373 + t1326 * t1372;
t1383 = 0.1e1 / pkin(3) ^ 2;
t1497 = t1295 ^ 2 * t1383;
t1445 = t1344 * t1497;
t1361 = sin(qJ(3,2));
t1346 = 0.1e1 / t1361 ^ 2;
t1354 = legFrame(2,2);
t1324 = sin(t1354);
t1327 = cos(t1354);
t1297 = t1324 * t1373 + t1327 * t1372;
t1496 = t1297 ^ 2 * t1383;
t1443 = t1346 * t1496;
t1363 = sin(qJ(3,1));
t1348 = 0.1e1 / t1363 ^ 2;
t1355 = legFrame(1,2);
t1325 = sin(t1355);
t1328 = cos(t1355);
t1299 = t1325 * t1373 + t1328 * t1372;
t1495 = t1299 ^ 2 * t1383;
t1441 = t1348 * t1495;
t1371 = xDP(3);
t1335 = pkin(3) * t1371;
t1522 = (-pkin(5) - pkin(6));
t1342 = pkin(1) - t1522;
t1360 = sin(qJ(1,3));
t1366 = cos(qJ(1,3));
t1424 = -t1323 * t1372 + t1326 * t1373;
t1365 = cos(qJ(3,3));
t1439 = pkin(3) * (t1365 + 0.1e1) * (t1365 - 0.1e1) * t1360;
t1467 = t1342 * t1371;
t1493 = t1295 * t1365;
t1349 = t1365 ^ 2;
t1541 = -t1349 + 0.1e1;
t1232 = ((qJ(2,3) * t1371 + t1424 * t1342) * t1366 + (t1424 * qJ(2,3) - t1467) * t1360 + pkin(3) * t1493) * t1359 - t1424 * t1439 + qJ(2,3) * t1493 + t1541 * t1366 * t1335;
t1277 = -t1360 * t1371 + t1424 * t1366;
t1304 = g(1) * t1326 - g(2) * t1323;
t1284 = -g(3) * t1360 + t1304 * t1366;
t1518 = pkin(3) * t1359;
t1316 = qJ(2,3) + t1518;
t1308 = 0.1e1 / t1316 ^ 2;
t1343 = 0.1e1 / t1359;
t1484 = t1308 * t1343;
t1544 = 0.2e1 * t1277 * t1232 * t1484 - t1284;
t1362 = sin(qJ(1,2));
t1368 = cos(qJ(1,2));
t1423 = -t1324 * t1372 + t1327 * t1373;
t1367 = cos(qJ(3,2));
t1438 = pkin(3) * (t1367 + 0.1e1) * (t1367 - 0.1e1) * t1362;
t1491 = t1297 * t1367;
t1350 = t1367 ^ 2;
t1540 = -t1350 + 0.1e1;
t1233 = ((qJ(2,2) * t1371 + t1423 * t1342) * t1368 + (t1423 * qJ(2,2) - t1467) * t1362 + pkin(3) * t1491) * t1361 - t1423 * t1438 + qJ(2,2) * t1491 + t1540 * t1368 * t1335;
t1278 = -t1362 * t1371 + t1423 * t1368;
t1305 = g(1) * t1327 - g(2) * t1324;
t1286 = -g(3) * t1362 + t1305 * t1368;
t1517 = pkin(3) * t1361;
t1317 = qJ(2,2) + t1517;
t1311 = 0.1e1 / t1317 ^ 2;
t1345 = 0.1e1 / t1361;
t1481 = t1311 * t1345;
t1543 = 0.2e1 * t1278 * t1233 * t1481 - t1286;
t1364 = sin(qJ(1,1));
t1370 = cos(qJ(1,1));
t1422 = -t1325 * t1372 + t1328 * t1373;
t1369 = cos(qJ(3,1));
t1437 = pkin(3) * (t1369 + 0.1e1) * (t1369 - 0.1e1) * t1364;
t1489 = t1299 * t1369;
t1351 = t1369 ^ 2;
t1539 = -t1351 + 0.1e1;
t1234 = ((qJ(2,1) * t1371 + t1422 * t1342) * t1370 + (t1422 * qJ(2,1) - t1467) * t1364 + pkin(3) * t1489) * t1363 - t1422 * t1437 + qJ(2,1) * t1489 + t1539 * t1370 * t1335;
t1279 = -t1364 * t1371 + t1422 * t1370;
t1306 = g(1) * t1328 - g(2) * t1325;
t1288 = -g(3) * t1364 + t1306 * t1370;
t1516 = pkin(3) * t1363;
t1318 = qJ(2,1) + t1516;
t1314 = 0.1e1 / t1318 ^ 2;
t1347 = 0.1e1 / t1363;
t1478 = t1314 * t1347;
t1542 = 0.2e1 * t1279 * t1234 * t1478 - t1288;
t1538 = t1445 * t1365;
t1537 = t1443 * t1367;
t1536 = t1441 * t1369;
t1307 = 0.1e1 / t1316;
t1309 = t1307 * t1308;
t1356 = xDDP(3);
t1357 = xDDP(2);
t1358 = xDDP(1);
t1483 = t1308 * t1365;
t1502 = t1277 * t1342;
t1503 = t1277 * t1308;
t1515 = t1232 * t1343;
t1217 = (-t1232 * t1309 + 0.2e1 * t1295 * t1483) * t1343 * t1277 + (-t1360 * t1356 - (-t1502 + t1515) * t1503 + (-t1323 * t1357 + t1326 * t1358) * t1366) * t1307;
t1523 = 0.2e1 * qJ(2,3);
t1205 = t1217 * t1523 + t1544;
t1374 = pkin(1) + pkin(5);
t1202 = t1374 * t1445 + t1205;
t1382 = 0.1e1 / pkin(3);
t1465 = t1343 * t1382;
t1271 = -t1343 * t1538 + (-t1323 * t1358 - t1326 * t1357) * t1465;
t1494 = t1295 * t1307;
t1433 = t1277 * t1382 * t1494;
t1427 = t1343 * t1433;
t1244 = t1271 * t1374 + t1427 * t1523;
t1253 = t1271 * t1365 - t1343 * t1497;
t1320 = 0.2e1 * t1349 - 0.1e1;
t1332 = g(3) * t1366;
t1375 = pkin(1) * g(3);
t1376 = qJ(2,3) ^ 2;
t1428 = qJ(2,3) * t1360 + t1342 * t1366;
t1460 = t1359 * t1365;
t1262 = t1428 * t1326 * t1359 + t1323 * t1365 * qJ(2,3) + (t1541 * t1326 * t1360 + t1323 * t1460) * pkin(3);
t1472 = t1326 * t1365;
t1265 = (pkin(3) * t1472 - t1428 * t1323) * t1359 + t1323 * t1439 + qJ(2,3) * t1472;
t1289 = t1316 * t1366 - t1342 * t1360;
t1381 = pkin(3) ^ 2;
t1415 = -(pkin(6) ^ 2) - t1381 + ((-2 * pkin(6) - pkin(5)) * pkin(5)) + ((2 * t1522) - pkin(1)) * pkin(1);
t1454 = t1342 * t1515;
t1466 = t1343 * t1365;
t1418 = (-t1295 * t1342 * t1466 + (t1454 + (-0.2e1 * qJ(2,3) * t1518 + t1349 * t1381 - t1376 + t1415) * t1277) * t1307) * t1503 + t1277 * t1309 * t1454 - ((t1307 * t1365 * t1502 + t1295 * t1343) * t1359 + qJ(2,3) * t1295 * t1465) * t1344 * t1494 - t1289 * t1307 * t1356 + (-t1262 * t1358 - t1265 * t1357) * t1307 * t1343;
t1414 = -t1332 - t1418;
t1421 = t1304 * t1360;
t1512 = t1271 * t1359;
t1521 = pkin(1) * t1217;
t1526 = 2 * MDP(8);
t1528 = -t1421 - t1521;
t1399 = MDP(7) * (t1217 * t1365 + 0.2e1 * t1433) * t1365 + MDP(1) * t1217 + MDP(10) * (-t1512 - t1538) + MDP(12) * (t1202 * t1359 - t1244 * t1365) + MDP(13) * (t1202 * t1365 + t1244 * t1359) + MDP(2) * (t1332 + t1421) + MDP(3) * t1284 + MDP(4) * (t1414 - t1421 - 0.2e1 * t1521) + MDP(5) * t1205 + MDP(6) * (t1217 * t1376 + t1375 * t1366 + t1544 * qJ(2,3) + (t1418 - t1528) * pkin(1)) + MDP(9) * t1253 + (-t1217 * t1460 + t1320 * t1427) * t1526;
t1535 = t1366 * t1399;
t1310 = 0.1e1 / t1317;
t1312 = t1310 * t1311;
t1480 = t1311 * t1367;
t1500 = t1278 * t1342;
t1501 = t1278 * t1311;
t1514 = t1233 * t1345;
t1218 = (-t1233 * t1312 + 0.2e1 * t1297 * t1480) * t1345 * t1278 + (-t1362 * t1356 - (-t1500 + t1514) * t1501 + (-t1324 * t1357 + t1327 * t1358) * t1368) * t1310;
t1524 = 0.2e1 * qJ(2,2);
t1206 = t1218 * t1524 + t1543;
t1203 = t1374 * t1443 + t1206;
t1463 = t1345 * t1382;
t1272 = -t1345 * t1537 + (-t1324 * t1358 - t1327 * t1357) * t1463;
t1492 = t1297 * t1310;
t1432 = t1278 * t1382 * t1492;
t1426 = t1345 * t1432;
t1245 = t1272 * t1374 + t1426 * t1524;
t1254 = t1272 * t1367 - t1345 * t1496;
t1321 = 0.2e1 * t1350 - 0.1e1;
t1333 = g(3) * t1368;
t1377 = qJ(2,2) ^ 2;
t1429 = qJ(2,2) * t1362 + t1342 * t1368;
t1459 = t1361 * t1367;
t1263 = t1429 * t1327 * t1361 + t1324 * t1367 * qJ(2,2) + (t1540 * t1327 * t1362 + t1324 * t1459) * pkin(3);
t1470 = t1327 * t1367;
t1266 = (pkin(3) * t1470 - t1429 * t1324) * t1361 + t1324 * t1438 + qJ(2,2) * t1470;
t1290 = t1317 * t1368 - t1342 * t1362;
t1453 = t1342 * t1514;
t1464 = t1345 * t1367;
t1417 = (-t1297 * t1342 * t1464 + (t1453 + (-0.2e1 * qJ(2,2) * t1517 + t1350 * t1381 - t1377 + t1415) * t1278) * t1310) * t1501 + t1278 * t1312 * t1453 - ((t1310 * t1367 * t1500 + t1297 * t1345) * t1361 + qJ(2,2) * t1297 * t1463) * t1346 * t1492 - t1290 * t1310 * t1356 + (-t1263 * t1358 - t1266 * t1357) * t1310 * t1345;
t1413 = -t1333 - t1417;
t1420 = t1305 * t1362;
t1511 = t1272 * t1361;
t1520 = pkin(1) * t1218;
t1530 = -t1420 - t1520;
t1398 = MDP(7) * (t1218 * t1367 + 0.2e1 * t1432) * t1367 + MDP(1) * t1218 + MDP(10) * (-t1511 - t1537) + MDP(12) * (t1203 * t1361 - t1245 * t1367) + MDP(13) * (t1203 * t1367 + t1245 * t1361) + MDP(2) * (t1333 + t1420) + MDP(3) * t1286 + MDP(4) * (t1413 - t1420 - 0.2e1 * t1520) + MDP(5) * t1206 + MDP(6) * (t1218 * t1377 + t1375 * t1368 + t1543 * qJ(2,2) + (t1417 - t1530) * pkin(1)) + MDP(9) * t1254 + (-t1218 * t1459 + t1321 * t1426) * t1526;
t1534 = t1368 * t1398;
t1313 = 0.1e1 / t1318;
t1315 = t1313 * t1314;
t1477 = t1314 * t1369;
t1498 = t1279 * t1342;
t1499 = t1279 * t1314;
t1513 = t1234 * t1347;
t1219 = (-t1234 * t1315 + 0.2e1 * t1299 * t1477) * t1347 * t1279 + (-t1364 * t1356 - (-t1498 + t1513) * t1499 + (-t1325 * t1357 + t1328 * t1358) * t1370) * t1313;
t1525 = 0.2e1 * qJ(2,1);
t1207 = t1219 * t1525 + t1542;
t1204 = t1374 * t1441 + t1207;
t1461 = t1347 * t1382;
t1273 = -t1347 * t1536 + (-t1325 * t1358 - t1328 * t1357) * t1461;
t1490 = t1299 * t1313;
t1431 = t1279 * t1382 * t1490;
t1425 = t1347 * t1431;
t1246 = t1273 * t1374 + t1425 * t1525;
t1255 = t1273 * t1369 - t1347 * t1495;
t1322 = 0.2e1 * t1351 - 0.1e1;
t1334 = g(3) * t1370;
t1378 = qJ(2,1) ^ 2;
t1430 = qJ(2,1) * t1364 + t1342 * t1370;
t1458 = t1363 * t1369;
t1264 = t1430 * t1328 * t1363 + t1325 * t1369 * qJ(2,1) + (t1539 * t1328 * t1364 + t1325 * t1458) * pkin(3);
t1468 = t1328 * t1369;
t1267 = (pkin(3) * t1468 - t1430 * t1325) * t1363 + t1325 * t1437 + qJ(2,1) * t1468;
t1291 = t1318 * t1370 - t1342 * t1364;
t1452 = t1342 * t1513;
t1462 = t1347 * t1369;
t1416 = (-t1299 * t1342 * t1462 + (t1452 + (-0.2e1 * qJ(2,1) * t1516 + t1351 * t1381 - t1378 + t1415) * t1279) * t1313) * t1499 + t1279 * t1315 * t1452 - ((t1313 * t1369 * t1498 + t1299 * t1347) * t1363 + qJ(2,1) * t1299 * t1461) * t1348 * t1490 - t1291 * t1313 * t1356 + (-t1264 * t1358 - t1267 * t1357) * t1313 * t1347;
t1412 = -t1334 - t1416;
t1419 = t1306 * t1364;
t1510 = t1273 * t1363;
t1519 = pkin(1) * t1219;
t1532 = -t1419 - t1519;
t1397 = MDP(7) * (t1219 * t1369 + 0.2e1 * t1431) * t1369 + MDP(1) * t1219 + MDP(10) * (-t1510 - t1536) + MDP(12) * (t1204 * t1363 - t1246 * t1369) + MDP(13) * (t1204 * t1369 + t1246 * t1363) + MDP(2) * (t1334 + t1419) + MDP(3) * t1288 + MDP(4) * (t1412 - t1419 - 0.2e1 * t1519) + MDP(5) * t1207 + MDP(6) * (t1219 * t1378 + t1375 * t1370 + t1542 * qJ(2,1) + (t1416 - t1532) * pkin(1)) + MDP(9) * t1255 + (-t1219 * t1458 + t1322 * t1425) * t1526;
t1533 = t1370 * t1397;
t1276 = t1279 ^ 2;
t1505 = t1276 * t1314;
t1408 = qJ(2,1) * t1505 - t1412;
t1531 = t1374 * t1219 + t1408 + t1419;
t1275 = t1278 ^ 2;
t1507 = t1275 * t1311;
t1407 = qJ(2,2) * t1507 - t1413;
t1529 = t1374 * t1218 + t1407 + t1420;
t1274 = t1277 ^ 2;
t1509 = t1274 * t1308;
t1406 = qJ(2,3) * t1509 - t1414;
t1527 = t1374 * t1217 + t1406 + t1421;
t1508 = t1274 * t1309;
t1506 = t1275 * t1312;
t1504 = t1276 * t1315;
t1476 = t1323 * t1343;
t1475 = t1324 * t1345;
t1474 = t1325 * t1347;
t1473 = t1326 * t1343;
t1471 = t1327 * t1345;
t1469 = t1328 * t1347;
t1457 = t1217 * t1466;
t1456 = t1218 * t1464;
t1455 = t1219 * t1462;
t1451 = t1274 * t1483;
t1450 = t1343 * t1508;
t1449 = t1275 * t1480;
t1448 = t1345 * t1506;
t1447 = t1276 * t1477;
t1446 = t1347 * t1504;
t1436 = t1274 * t1320 * t1484;
t1435 = t1275 * t1321 * t1481;
t1434 = t1276 * t1322 * t1478;
t1411 = MDP(12) * (-t1359 * t1509 + t1253) + MDP(13) * (-t1512 + (-t1445 - t1509) * t1365) + MDP(4) * t1217 + MDP(6) * (-t1406 + t1528);
t1410 = MDP(12) * (-t1361 * t1507 + t1254) + MDP(13) * (-t1511 + (-t1443 - t1507) * t1367) + MDP(4) * t1218 + MDP(6) * (-t1407 + t1530);
t1409 = MDP(12) * (-t1363 * t1505 + t1255) + MDP(13) * (-t1510 + (-t1441 - t1505) * t1369) + MDP(4) * t1219 + MDP(6) * (-t1408 + t1532);
t1405 = t1411 * t1343;
t1404 = t1410 * t1345;
t1403 = t1409 * t1347;
t1303 = g(1) * t1325 + g(2) * t1328;
t1302 = g(1) * t1324 + g(2) * t1327;
t1301 = g(1) * t1323 + g(2) * t1326;
t1189 = t1303 * t1363 - t1531 * t1369;
t1188 = t1303 * t1369 + t1531 * t1363;
t1187 = t1302 * t1361 - t1529 * t1367;
t1186 = t1302 * t1367 + t1529 * t1361;
t1185 = t1301 * t1359 - t1527 * t1365;
t1184 = t1301 * t1365 + t1527 * t1359;
t1 = [(-t1262 * t1450 - t1263 * t1448 - t1264 * t1446) * MDP(5) + (t1358 - g(1)) * MDP(14) + (t1264 * t1403 + t1328 * t1533) * t1313 + (t1263 * t1404 + t1327 * t1534) * t1310 + (t1262 * t1405 + t1326 * t1535) * t1307 + ((-t1323 * t1451 - t1324 * t1449 - t1325 * t1447) * MDP(7) + (-t1323 * t1436 - t1324 * t1435 - t1325 * t1434) * MDP(8) + (-t1323 * t1457 - t1324 * t1456 - t1325 * t1455) * MDP(9) + (t1217 * t1323 + t1218 * t1324 + t1219 * t1325) * MDP(10) + (-t1271 * t1476 - t1272 * t1475 - t1273 * t1474) * MDP(11) + (-t1185 * t1476 - t1187 * t1475 - t1189 * t1474) * MDP(12) + (-t1184 * t1476 - t1186 * t1475 - t1188 * t1474) * MDP(13)) * t1382; (-t1265 * t1450 - t1266 * t1448 - t1267 * t1446) * MDP(5) + (t1357 - g(2)) * MDP(14) + (t1267 * t1403 - t1325 * t1533) * t1313 + (t1266 * t1404 - t1324 * t1534) * t1310 + (t1265 * t1405 - t1323 * t1535) * t1307 + ((-t1326 * t1451 - t1327 * t1449 - t1328 * t1447) * MDP(7) + (-t1326 * t1436 - t1327 * t1435 - t1328 * t1434) * MDP(8) + (-t1326 * t1457 - t1327 * t1456 - t1328 * t1455) * MDP(9) + (t1217 * t1326 + t1218 * t1327 + t1219 * t1328) * MDP(10) + (-t1271 * t1473 - t1272 * t1471 - t1273 * t1469) * MDP(11) + (-t1185 * t1473 - t1187 * t1471 - t1189 * t1469) * MDP(12) + (-t1184 * t1473 - t1186 * t1471 - t1188 * t1469) * MDP(13)) * t1382; (-t1289 * t1508 - t1290 * t1506 - t1291 * t1504) * MDP(5) + (t1356 - g(3)) * MDP(14) + (t1409 * t1291 - t1397 * t1364) * t1313 + (t1410 * t1290 - t1398 * t1362) * t1310 + (t1411 * t1289 - t1399 * t1360) * t1307;];
tauX  = t1;
