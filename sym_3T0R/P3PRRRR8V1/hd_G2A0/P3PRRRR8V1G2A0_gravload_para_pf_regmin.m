% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRRR8V1G2A0
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
% Datum: 2020-08-06 17:03
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRRR8V1G2A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_regmin: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G2A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:03:27
% EndTime: 2020-08-06 17:03:28
% DurationCPUTime: 1.20s
% Computational Cost: add. (762->162), mult. (2074->337), div. (126->7), fcn. (2326->22), ass. (0->172)
t1439 = sin(pkin(3));
t1440 = cos(pkin(6));
t1441 = cos(pkin(3));
t1528 = t1440 * t1441;
t1408 = -t1439 * g(1) + g(2) * t1528;
t1409 = g(1) * t1528 + t1439 * g(2);
t1444 = legFrame(1,2);
t1431 = sin(t1444);
t1434 = cos(t1444);
t1415 = t1434 * g(1) - t1431 * g(2);
t1438 = sin(pkin(6));
t1423 = t1441 * t1438 * g(3);
t1428 = t1440 * g(3);
t1450 = sin(qJ(2,1));
t1456 = cos(qJ(2,1));
t1355 = (t1408 * t1431 - t1409 * t1434 + t1423) * t1456 + t1450 * (t1415 * t1438 + t1428);
t1449 = sin(qJ(3,1));
t1550 = t1355 * t1449;
t1442 = legFrame(3,2);
t1429 = sin(t1442);
t1432 = cos(t1442);
t1413 = t1432 * g(1) - t1429 * g(2);
t1446 = sin(qJ(2,3));
t1452 = cos(qJ(2,3));
t1549 = (t1408 * t1429 - t1409 * t1432 + t1423) * t1452 + t1446 * (t1413 * t1438 + t1428);
t1443 = legFrame(2,2);
t1430 = sin(t1443);
t1433 = cos(t1443);
t1414 = t1433 * g(1) - t1430 * g(2);
t1448 = sin(qJ(2,2));
t1454 = cos(qJ(2,2));
t1548 = (t1408 * t1430 - t1409 * t1433 + t1423) * t1454 + t1448 * (t1414 * t1438 + t1428);
t1451 = cos(qJ(3,3));
t1514 = t1446 * t1451;
t1416 = pkin(2) * t1514 - t1452 * pkin(5);
t1445 = sin(qJ(3,3));
t1527 = t1441 * t1445;
t1384 = pkin(2) * t1527 + t1416 * t1439;
t1547 = 0.1e1 / t1384;
t1453 = cos(qJ(3,2));
t1509 = t1448 * t1453;
t1417 = pkin(2) * t1509 - t1454 * pkin(5);
t1447 = sin(qJ(3,2));
t1525 = t1441 * t1447;
t1385 = pkin(2) * t1525 + t1417 * t1439;
t1546 = 0.1e1 / t1385;
t1455 = cos(qJ(3,1));
t1505 = t1450 * t1455;
t1418 = pkin(2) * t1505 - t1456 * pkin(5);
t1523 = t1441 * t1449;
t1386 = pkin(2) * t1523 + t1418 * t1439;
t1545 = 0.1e1 / t1386;
t1544 = pkin(2) * t1451;
t1543 = pkin(2) * t1453;
t1542 = pkin(2) * t1455;
t1541 = t1549 * t1547;
t1540 = t1548 * t1546;
t1539 = t1355 * t1545;
t1538 = t1547 / t1451;
t1537 = t1546 / t1453;
t1437 = 0.1e1 / t1455;
t1536 = t1545 * t1437;
t1410 = t1429 * g(1) + t1432 * g(2);
t1535 = t1547 * t1410;
t1411 = t1430 * g(1) + t1433 * g(2);
t1534 = t1546 * t1411;
t1412 = t1431 * g(1) + t1434 * g(2);
t1533 = t1545 * t1412;
t1532 = t1410 * t1439;
t1531 = t1411 * t1439;
t1530 = t1412 * t1439;
t1529 = t1439 * t1449;
t1526 = t1441 * t1446;
t1524 = t1441 * t1448;
t1522 = t1441 * t1450;
t1521 = t1441 * t1452;
t1520 = t1441 * t1454;
t1519 = t1441 * t1456;
t1518 = t1445 * t1439;
t1517 = t1445 * t1446;
t1516 = t1445 * t1452;
t1513 = t1447 * t1439;
t1512 = t1447 * t1448;
t1511 = t1447 * t1454;
t1508 = t1449 * t1450;
t1507 = t1449 * t1456;
t1504 = t1451 * t1439;
t1503 = t1451 * t1452;
t1502 = t1453 * t1439;
t1501 = t1453 * t1454;
t1500 = t1455 * t1439;
t1499 = t1455 * t1456;
t1402 = t1441 * t1517 + t1504;
t1465 = -t1402 * t1438 + t1440 * t1516;
t1498 = t1465 * t1541;
t1403 = t1441 * t1512 + t1502;
t1463 = -t1403 * t1438 + t1440 * t1511;
t1497 = t1463 * t1540;
t1391 = -t1438 * t1446 + t1440 * t1521;
t1394 = t1438 * t1452 + t1440 * t1526;
t1496 = (-pkin(5) * t1394 - t1391 * t1544) * t1538;
t1392 = -t1438 * t1448 + t1440 * t1520;
t1395 = t1438 * t1454 + t1440 * t1524;
t1495 = (-pkin(5) * t1395 - t1392 * t1543) * t1537;
t1393 = -t1438 * t1450 + t1440 * t1519;
t1396 = t1438 * t1456 + t1440 * t1522;
t1494 = (-pkin(5) * t1396 - t1393 * t1542) * t1536;
t1464 = t1402 * t1440 + t1438 * t1516;
t1493 = t1464 * t1538;
t1462 = t1403 * t1440 + t1438 * t1511;
t1492 = t1462 * t1537;
t1404 = t1441 * t1508 + t1500;
t1461 = t1404 * t1440 + t1438 * t1507;
t1491 = t1461 * t1536;
t1390 = t1438 * t1522 - t1440 * t1456;
t1490 = (t1390 * t1449 + t1438 * t1500) * t1545 * t1431;
t1489 = t1429 * t1538;
t1488 = t1432 * t1538;
t1487 = t1430 * t1537;
t1486 = t1433 * t1537;
t1485 = t1434 * t1536;
t1373 = -t1404 * t1438 + t1440 * t1507;
t1484 = t1549 * t1445 * t1538;
t1483 = t1548 * t1447 * t1537;
t1397 = t1438 * t1521 + t1440 * t1446;
t1400 = -t1438 * t1526 + t1440 * t1452;
t1368 = -t1400 * pkin(5) + t1397 * t1544;
t1482 = t1368 * t1489;
t1481 = t1368 * t1488;
t1398 = t1438 * t1520 + t1440 * t1448;
t1401 = -t1438 * t1524 + t1440 * t1454;
t1369 = -t1401 * pkin(5) + t1398 * t1543;
t1480 = t1369 * t1487;
t1479 = t1369 * t1486;
t1399 = t1438 * t1519 + t1440 * t1450;
t1370 = t1390 * pkin(5) + t1399 * t1542;
t1478 = t1370 * t1431 * t1536;
t1477 = t1370 * t1485;
t1476 = t1465 * t1489;
t1475 = t1465 * t1488;
t1474 = t1463 * t1487;
t1473 = t1463 * t1486;
t1472 = t1373 * t1485;
t1471 = t1437 * t1490;
t1470 = t1465 * t1484;
t1469 = t1463 * t1483;
t1468 = pkin(2) * t1518 - t1416 * t1441;
t1467 = pkin(2) * t1513 - t1417 * t1441;
t1466 = pkin(2) * t1529 - t1418 * t1441;
t1457 = 0.1e1 / pkin(2);
t1421 = pkin(2) * t1499 + pkin(5) * t1450;
t1420 = pkin(2) * t1501 + pkin(5) * t1448;
t1419 = pkin(2) * t1503 + pkin(5) * t1446;
t1407 = t1441 * t1509 - t1513;
t1406 = t1441 * t1514 - t1518;
t1405 = t1441 * t1505 - t1529;
t1364 = -t1438 * t1421 + t1466 * t1440;
t1363 = -t1438 * t1420 + t1467 * t1440;
t1362 = -t1438 * t1419 + t1468 * t1440;
t1361 = -g(3) * t1390 + t1415 * t1396 + t1450 * t1530;
t1360 = g(3) * t1401 + t1414 * t1395 + t1448 * t1531;
t1359 = g(3) * t1400 + t1413 * t1394 + t1446 * t1532;
t1358 = g(3) * t1399 - t1415 * t1393 - t1456 * t1530;
t1357 = g(3) * t1398 - t1414 * t1392 - t1454 * t1531;
t1356 = g(3) * t1397 - t1413 * t1391 - t1452 * t1532;
t1349 = g(3) * (-t1406 * t1438 + t1440 * t1503) + t1413 * (t1406 * t1440 + t1438 * t1503) + t1410 * (t1446 * t1504 + t1527);
t1348 = g(3) * t1373 + t1415 * t1461 - t1412 * (-t1439 * t1508 + t1441 * t1455);
t1347 = g(3) * t1463 + t1414 * t1462 - t1411 * (-t1439 * t1512 + t1441 * t1453);
t1346 = g(3) * t1465 + t1413 * t1464 - t1410 * (-t1439 * t1517 + t1441 * t1451);
t1345 = (-t1407 * t1438 + t1440 * t1501) * g(3) + t1414 * (t1407 * t1440 + t1438 * t1501) + t1411 * (t1448 * t1502 + t1525);
t1344 = (-t1405 * t1438 + t1440 * t1499) * g(3) + t1415 * (t1405 * t1440 + t1438 * t1499) + t1412 * (t1450 * t1500 + t1523);
t1 = [-(-t1364 * t1434 + t1386 * t1431) * t1533 - (-t1363 * t1433 + t1385 * t1430) * t1534 - (-t1362 * t1432 + t1384 * t1429) * t1535, 0, t1356 * t1475 + t1357 * t1473 + t1358 * t1472, t1359 * t1475 + t1360 * t1473 + t1361 * t1472, 0, 0, 0, 0, 0, t1432 * t1498 + t1433 * t1497 + t1434 * t1373 * t1539 + (-t1346 * t1481 - t1347 * t1479 - t1348 * t1477) * t1457, -t1432 * t1470 - t1433 * t1469 - t1472 * t1550 + (-t1344 * t1477 - t1345 * t1479 - t1349 * t1481) * t1457, -g(1); -(t1364 * t1431 + t1386 * t1434) * t1533 - (t1363 * t1430 + t1385 * t1433) * t1534 - (t1362 * t1429 + t1384 * t1432) * t1535, 0, -t1356 * t1476 - t1357 * t1474 + t1358 * t1471, -t1359 * t1476 - t1360 * t1474 + t1361 * t1471, 0, 0, 0, 0, 0, -t1429 * t1498 - t1430 * t1497 + t1355 * t1490 + (t1346 * t1482 + t1347 * t1480 + t1348 * t1478) * t1457, t1429 * t1470 + t1430 * t1469 - t1471 * t1550 + (t1344 * t1478 + t1345 * t1480 + t1349 * t1482) * t1457, -g(2); -(t1440 * t1421 + t1466 * t1438) * t1533 - (t1440 * t1420 + t1467 * t1438) * t1534 - (t1440 * t1419 + t1468 * t1438) * t1535, 0, -t1356 * t1493 - t1357 * t1492 - t1358 * t1491, -t1359 * t1493 - t1360 * t1492 - t1361 * t1491, 0, 0, 0, 0, 0, -t1464 * t1541 - t1462 * t1540 - t1461 * t1539 + (t1346 * t1496 + t1347 * t1495 + t1348 * t1494) * t1457, t1464 * t1484 + t1462 * t1483 + t1491 * t1550 + (t1344 * t1494 + t1345 * t1495 + t1349 * t1496) * t1457, -g(3);];
tau_reg  = t1;
