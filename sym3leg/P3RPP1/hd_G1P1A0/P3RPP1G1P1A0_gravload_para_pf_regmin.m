% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RPP1G1P1A0
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x13]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:53
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RPP1G1P1A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPP1G1P1A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPP1G1P1A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPP1G1P1A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RPP1G1P1A0_gravload_para_pf_regmin: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPP1G1P1A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPP1G1P1A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:52:42
% EndTime: 2019-05-03 14:52:44
% DurationCPUTime: 1.73s
% Computational Cost: add. (2153->224), mult. (3130->363), div. (99->3), fcn. (1468->14), ass. (0->184)
t1457 = pkin(1) + qJ(3,3);
t1477 = koppelP(3,1);
t1501 = t1457 * t1477;
t1474 = koppelP(3,2);
t1521 = qJ(2,3) * t1474;
t1413 = t1501 - t1521;
t1502 = t1457 * t1474;
t1520 = (qJ(2,3) * t1477);
t1416 = t1502 + t1520;
t1460 = sin(qJ(1,3));
t1463 = cos(qJ(1,3));
t1369 = t1413 * t1463 + t1416 * t1460;
t1372 = t1413 * t1460 - t1416 * t1463;
t1454 = legFrame(3,3);
t1440 = sin(t1454);
t1443 = cos(t1454);
t1467 = xP(3);
t1446 = sin(t1467);
t1447 = cos(t1467);
t1324 = (t1369 * t1447 + t1372 * t1446) * t1443 - (-t1369 * t1446 + t1372 * t1447) * t1440;
t1458 = pkin(1) + qJ(3,2);
t1478 = koppelP(2,1);
t1497 = t1458 * t1478;
t1475 = koppelP(2,2);
t1523 = qJ(2,2) * t1475;
t1414 = t1497 - t1523;
t1498 = t1458 * t1475;
t1522 = (qJ(2,2) * t1478);
t1417 = t1498 + t1522;
t1461 = sin(qJ(1,2));
t1464 = cos(qJ(1,2));
t1370 = t1414 * t1464 + t1417 * t1461;
t1373 = t1414 * t1461 - t1417 * t1464;
t1455 = legFrame(2,3);
t1441 = sin(t1455);
t1444 = cos(t1455);
t1325 = (t1370 * t1447 + t1373 * t1446) * t1444 - (-t1370 * t1446 + t1373 * t1447) * t1441;
t1459 = pkin(1) + qJ(3,1);
t1479 = koppelP(1,1);
t1493 = t1459 * t1479;
t1476 = koppelP(1,2);
t1525 = qJ(2,1) * t1476;
t1415 = t1493 - t1525;
t1494 = t1459 * t1476;
t1524 = (qJ(2,1) * t1479);
t1418 = t1494 + t1524;
t1462 = sin(qJ(1,1));
t1465 = cos(qJ(1,1));
t1371 = t1415 * t1465 + t1418 * t1462;
t1374 = t1415 * t1462 - t1418 * t1465;
t1456 = legFrame(1,3);
t1442 = sin(t1456);
t1445 = cos(t1456);
t1326 = (t1371 * t1447 + t1374 * t1446) * t1445 - (-t1371 * t1446 + t1374 * t1447) * t1442;
t1398 = g(1) * t1442 - g(2) * t1445;
t1401 = g(1) * t1445 + g(2) * t1442;
t1526 = pkin(1) ^ 2 + 1;
t1490 = qJ(3,1) ^ 2 + t1526;
t1534 = 2 * pkin(1);
t1430 = qJ(3,1) * t1534 + t1490;
t1473 = (qJ(2,1) ^ 2);
t1412 = 1 / (t1473 + t1430);
t1538 = (t1398 * t1465 + t1401 * t1462) * t1412;
t1397 = g(1) * t1441 - g(2) * t1444;
t1400 = g(1) * t1444 + g(2) * t1441;
t1491 = qJ(3,2) ^ 2 + t1526;
t1429 = qJ(3,2) * t1534 + t1491;
t1471 = (qJ(2,2) ^ 2);
t1411 = 1 / (t1471 + t1429);
t1540 = (t1397 * t1464 + t1400 * t1461) * t1411;
t1396 = g(1) * t1440 - g(2) * t1443;
t1399 = g(1) * t1443 + g(2) * t1440;
t1492 = qJ(3,3) ^ 2 + t1526;
t1428 = qJ(3,3) * t1534 + t1492;
t1469 = (qJ(2,3) ^ 2);
t1410 = 1 / (t1469 + t1428);
t1542 = (t1396 * t1463 + t1399 * t1460) * t1410;
t1316 = t1324 * t1542 + t1325 * t1540 + t1326 * t1538;
t1361 = t1398 * t1462 - t1401 * t1465;
t1539 = t1361 * t1412;
t1359 = t1397 * t1461 - t1400 * t1464;
t1541 = t1359 * t1411;
t1357 = t1396 * t1460 - t1399 * t1463;
t1543 = t1357 * t1410;
t1315 = t1324 * t1543 + t1325 * t1541 + t1326 * t1539;
t1504 = t1457 * t1460;
t1404 = qJ(2,3) * t1463 - t1504;
t1503 = t1457 * t1463;
t1407 = qJ(2,3) * t1460 + t1503;
t1351 = t1404 * t1443 - t1407 * t1440;
t1500 = t1458 * t1461;
t1405 = qJ(2,2) * t1464 - t1500;
t1499 = t1458 * t1464;
t1408 = qJ(2,2) * t1461 + t1499;
t1352 = t1405 * t1444 - t1408 * t1441;
t1496 = t1459 * t1462;
t1406 = qJ(2,1) * t1465 - t1496;
t1495 = t1459 * t1465;
t1409 = qJ(2,1) * t1462 + t1495;
t1353 = t1406 * t1445 - t1409 * t1442;
t1318 = t1351 * t1542 + t1352 * t1540 + t1353 * t1538;
t1317 = t1351 * t1543 + t1352 * t1541 + t1353 * t1539;
t1354 = t1404 * t1440 + t1407 * t1443;
t1355 = t1405 * t1441 + t1408 * t1444;
t1356 = t1406 * t1442 + t1409 * t1445;
t1320 = t1354 * t1542 + t1355 * t1540 + t1356 * t1538;
t1319 = t1354 * t1543 + t1355 * t1541 + t1356 * t1539;
t1537 = t1474 * t1492 + pkin(1) * t1520 + (t1474 * t1534 + t1520) * qJ(3,3);
t1536 = t1475 * t1491 + pkin(1) * t1522 + (t1475 * t1534 + t1522) * qJ(3,2);
t1535 = t1476 * t1490 + pkin(1) * t1524 + (t1476 * t1534 + t1524) * qJ(3,1);
t1533 = pkin(1) * g(2);
t1532 = g(1) * qJ(2,1);
t1531 = g(1) * qJ(2,2);
t1530 = g(1) * qJ(2,3);
t1529 = g(2) * qJ(2,1);
t1528 = g(2) * qJ(2,2);
t1527 = g(2) * qJ(2,3);
t1466 = pkin(1) * g(1);
t1434 = t1466 - t1527;
t1435 = t1530 + t1533;
t1519 = ((t1434 * t1460 - t1435 * t1463) * t1443 + (t1434 * t1463 + t1435 * t1460) * t1440) * t1410;
t1436 = t1466 - t1528;
t1437 = t1531 + t1533;
t1518 = ((t1436 * t1461 - t1437 * t1464) * t1444 + (t1436 * t1464 + t1437 * t1461) * t1441) * t1411;
t1438 = t1466 - t1529;
t1439 = t1532 + t1533;
t1517 = ((t1438 * t1462 - t1439 * t1465) * t1445 + (t1438 * t1465 + t1439 * t1462) * t1442) * t1412;
t1489 = qJ(2,1) * t1495;
t1488 = qJ(2,2) * t1499;
t1487 = qJ(2,3) * t1503;
t1390 = qJ(2,3) * t1501 - (t1469 * t1474) - t1474;
t1391 = qJ(2,3) * t1502 + (t1469 * t1477) + t1477;
t1345 = t1390 * t1463 + t1391 * t1460;
t1346 = t1390 * t1460 - t1391 * t1463;
t1392 = qJ(2,2) * t1497 - (t1471 * t1475) - t1475;
t1393 = qJ(2,2) * t1498 + (t1471 * t1478) + t1478;
t1347 = t1392 * t1464 + t1393 * t1461;
t1348 = t1392 * t1461 - t1393 * t1464;
t1394 = qJ(2,1) * t1493 - (t1473 * t1476) - t1476;
t1395 = qJ(2,1) * t1494 + (t1473 * t1479) + t1479;
t1349 = t1394 * t1465 + t1395 * t1462;
t1350 = t1394 * t1462 - t1395 * t1465;
t1486 = -((-t1345 * t1446 + t1346 * t1447) * t1443 + (t1345 * t1447 + t1346 * t1446) * t1440) * t1542 - ((-t1347 * t1446 + t1348 * t1447) * t1444 + (t1347 * t1447 + t1348 * t1446) * t1441) * t1540 - ((-t1349 * t1446 + t1350 * t1447) * t1445 + (t1349 * t1447 + t1350 * t1446) * t1442) * t1538;
t1451 = 1 + t1469;
t1384 = t1451 * t1460 + t1487;
t1452 = 1 + t1471;
t1385 = t1452 * t1461 + t1488;
t1453 = 1 + t1473;
t1386 = t1453 * t1462 + t1489;
t1431 = qJ(2,3) * t1504;
t1387 = -t1451 * t1463 + t1431;
t1432 = qJ(2,2) * t1500;
t1388 = -t1452 * t1464 + t1432;
t1433 = qJ(2,1) * t1496;
t1389 = -t1453 * t1465 + t1433;
t1485 = -(t1384 * t1443 - t1387 * t1440) * t1542 - (t1385 * t1444 - t1388 * t1441) * t1540 - (t1386 * t1445 - t1389 * t1442) * t1538;
t1484 = -(t1384 * t1440 + t1387 * t1443) * t1542 - (t1385 * t1441 + t1388 * t1444) * t1540 - (t1386 * t1442 + t1389 * t1445) * t1538;
t1427 = g(2) * t1459 + t1532;
t1426 = g(2) * t1458 + t1531;
t1425 = g(2) * t1457 + t1530;
t1424 = g(1) * t1459 - t1529;
t1423 = g(1) * t1458 - t1528;
t1422 = g(1) * t1457 - t1527;
t1403 = g(1) * t1447 + g(2) * t1446;
t1402 = g(1) * t1446 - g(2) * t1447;
t1383 = t1430 * t1465 + t1433;
t1382 = t1429 * t1464 + t1432;
t1381 = t1428 * t1463 + t1431;
t1380 = t1430 * t1462 - t1489;
t1379 = t1429 * t1461 - t1488;
t1378 = t1428 * t1460 - t1487;
t1377 = (t1430 * t1479) - t1459 * t1525;
t1376 = (t1429 * t1478) - t1458 * t1523;
t1375 = (t1428 * t1477) - t1457 * t1521;
t1341 = t1377 * t1465 + t1462 * t1535;
t1340 = t1377 * t1462 - t1465 * t1535;
t1339 = t1376 * t1464 + t1461 * t1536;
t1338 = t1376 * t1461 - t1464 * t1536;
t1337 = t1375 * t1463 + t1460 * t1537;
t1336 = t1375 * t1460 - t1463 * t1537;
t1335 = (t1424 * t1462 - t1427 * t1465) * t1445 + t1442 * (t1424 * t1465 + t1427 * t1462);
t1334 = (t1423 * t1461 - t1426 * t1464) * t1444 + t1441 * (t1423 * t1464 + t1426 * t1461);
t1333 = (t1422 * t1460 - t1425 * t1463) * t1443 + t1440 * (t1422 * t1463 + t1425 * t1460);
t1 = [0, t1318, -t1317, -t1318, t1317, t1351 * t1519 + t1352 * t1518 + t1353 * t1517 + t1485, t1317, t1318, (t1353 * t1335 + (-t1380 * t1442 + t1383 * t1445) * t1361) * t1412 + (t1352 * t1334 + (-t1379 * t1441 + t1382 * t1444) * t1359) * t1411 + (t1351 * t1333 + (-t1378 * t1440 + t1381 * t1443) * t1357) * t1410 + t1485, 0, 0, 0, -t1402 * t1446 - t1403 * t1447; 0, t1320, -t1319, -t1320, t1319, t1354 * t1519 + t1355 * t1518 + t1356 * t1517 + t1484, t1319, t1320, (t1356 * t1335 + (t1380 * t1445 + t1383 * t1442) * t1361) * t1412 + (t1355 * t1334 + (t1379 * t1444 + t1382 * t1441) * t1359) * t1411 + (t1354 * t1333 + (t1378 * t1443 + t1381 * t1440) * t1357) * t1410 + t1484, 0, 0, 0, t1402 * t1447 - t1403 * t1446; 0, t1316, -t1315, -t1316, t1315, t1324 * t1519 + t1325 * t1518 + t1326 * t1517 + t1486, t1315, t1316, (t1326 * t1335 + ((t1340 * t1447 - t1341 * t1446) * t1445 + (t1340 * t1446 + t1341 * t1447) * t1442) * t1361) * t1412 + (t1325 * t1334 + ((t1338 * t1447 - t1339 * t1446) * t1444 + (t1338 * t1446 + t1339 * t1447) * t1441) * t1359) * t1411 + (t1324 * t1333 + ((t1336 * t1447 - t1337 * t1446) * t1443 + (t1336 * t1446 + t1337 * t1447) * t1440) * t1357) * t1410 + t1486, 0, t1402, t1403, 0;];
tau_reg  = t1;
