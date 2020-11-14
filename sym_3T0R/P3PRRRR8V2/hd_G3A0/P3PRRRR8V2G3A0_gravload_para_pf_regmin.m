% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRRR8V2G3A0
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2020-08-06 18:05
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRRR8V2G3A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_regmin: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:05:21
% EndTime: 2020-08-06 18:05:22
% DurationCPUTime: 1.21s
% Computational Cost: add. (1086->210), mult. (2529->423), div. (81->4), fcn. (2667->22), ass. (0->176)
t1480 = sin(qJ(3,3));
t1581 = pkin(2) * t1480;
t1482 = sin(qJ(3,2));
t1580 = pkin(2) * t1482;
t1484 = sin(qJ(3,1));
t1579 = pkin(2) * t1484;
t1486 = cos(qJ(3,3));
t1578 = pkin(3) * t1486 ^ 2;
t1488 = cos(qJ(3,2));
t1577 = pkin(3) * t1488 ^ 2;
t1490 = cos(qJ(3,1));
t1576 = pkin(3) * t1490 ^ 2;
t1575 = pkin(3) * t1486;
t1574 = pkin(3) * t1488;
t1573 = pkin(3) * t1490;
t1473 = sin(pkin(8));
t1572 = t1473 * g(3);
t1481 = sin(qJ(2,3));
t1487 = cos(qJ(2,3));
t1492 = pkin(7) + pkin(6);
t1457 = pkin(2) * t1481 - t1492 * t1487;
t1474 = sin(pkin(4));
t1476 = cos(pkin(4));
t1543 = t1476 * t1480;
t1431 = pkin(3) * t1543 + t1474 * t1457;
t1551 = t1474 * t1481;
t1419 = 0.1e1 / (pkin(2) * t1543 + t1431 * t1486 + t1551 * t1578);
t1553 = t1473 * t1476;
t1449 = -t1474 * g(1) - g(2) * t1553;
t1450 = g(1) * t1553 - t1474 * g(2);
t1477 = legFrame(3,2);
t1464 = sin(t1477);
t1467 = cos(t1477);
t1454 = t1467 * g(1) - t1464 * g(2);
t1475 = cos(pkin(8));
t1547 = t1475 * t1476;
t1463 = g(3) * t1547;
t1571 = ((t1449 * t1464 + t1450 * t1467 + t1463) * t1487 + t1481 * (t1454 * t1475 - t1572)) * t1419;
t1483 = sin(qJ(2,2));
t1489 = cos(qJ(2,2));
t1458 = pkin(2) * t1483 - t1492 * t1489;
t1541 = t1476 * t1482;
t1432 = pkin(3) * t1541 + t1474 * t1458;
t1550 = t1474 * t1483;
t1420 = 0.1e1 / (pkin(2) * t1541 + t1432 * t1488 + t1550 * t1577);
t1478 = legFrame(2,2);
t1465 = sin(t1478);
t1468 = cos(t1478);
t1455 = t1468 * g(1) - t1465 * g(2);
t1570 = ((t1449 * t1465 + t1450 * t1468 + t1463) * t1489 + t1483 * (t1455 * t1475 - t1572)) * t1420;
t1485 = sin(qJ(2,1));
t1491 = cos(qJ(2,1));
t1459 = pkin(2) * t1485 - t1492 * t1491;
t1539 = t1476 * t1484;
t1433 = pkin(3) * t1539 + t1474 * t1459;
t1548 = t1474 * t1485;
t1421 = 0.1e1 / (pkin(2) * t1539 + t1433 * t1490 + t1548 * t1576);
t1479 = legFrame(1,2);
t1466 = sin(t1479);
t1469 = cos(t1479);
t1456 = t1469 * g(1) - t1466 * g(2);
t1569 = ((t1449 * t1466 + t1450 * t1469 + t1463) * t1491 + t1485 * (t1456 * t1475 - t1572)) * t1421;
t1536 = t1476 * t1487;
t1443 = t1473 * t1536 + t1475 * t1481;
t1460 = pkin(2) * t1487 + t1481 * t1492;
t1568 = (t1443 * t1575 + t1457 * t1475 + t1460 * t1553) * t1419;
t1534 = t1476 * t1489;
t1444 = t1473 * t1534 + t1475 * t1483;
t1461 = pkin(2) * t1489 + t1483 * t1492;
t1567 = (t1444 * t1574 + t1458 * t1475 + t1461 * t1553) * t1420;
t1532 = t1476 * t1491;
t1445 = t1473 * t1532 + t1475 * t1485;
t1462 = pkin(2) * t1491 + t1485 * t1492;
t1566 = (t1445 * t1573 + t1459 * t1475 + t1462 * t1553) * t1421;
t1542 = t1476 * t1481;
t1546 = t1475 * t1487;
t1434 = t1473 * t1542 - t1546;
t1523 = t1486 * t1474;
t1426 = t1480 * t1434 + t1473 * t1523;
t1565 = t1419 * t1426;
t1451 = t1464 * g(1) + t1467 * g(2);
t1564 = t1419 * t1451;
t1563 = t1419 * t1464;
t1562 = t1419 * t1467;
t1540 = t1476 * t1483;
t1545 = t1475 * t1489;
t1435 = t1473 * t1540 - t1545;
t1522 = t1488 * t1474;
t1427 = t1482 * t1435 + t1473 * t1522;
t1561 = t1420 * t1427;
t1452 = t1465 * g(1) + t1468 * g(2);
t1560 = t1420 * t1452;
t1559 = t1420 * t1465;
t1558 = t1420 * t1468;
t1538 = t1476 * t1485;
t1544 = t1475 * t1491;
t1436 = t1473 * t1538 - t1544;
t1521 = t1490 * t1474;
t1425 = t1484 * t1436 + t1473 * t1521;
t1557 = t1421 * t1425;
t1453 = t1466 * g(1) + t1469 * g(2);
t1556 = t1421 * t1453;
t1555 = t1421 * t1466;
t1554 = t1421 * t1469;
t1552 = t1474 * t1473;
t1549 = t1474 * t1484;
t1537 = t1476 * t1486;
t1535 = t1476 * t1488;
t1533 = t1476 * t1490;
t1531 = t1480 * t1474;
t1530 = t1480 * t1481;
t1529 = t1480 * t1487;
t1528 = t1482 * t1474;
t1527 = t1482 * t1483;
t1526 = t1482 * t1489;
t1525 = t1484 * t1485;
t1524 = t1484 * t1491;
t1520 = t1480 * t1571;
t1519 = t1486 * t1571;
t1518 = t1482 * t1570;
t1517 = t1488 * t1570;
t1516 = t1484 * t1569;
t1515 = t1490 * t1569;
t1437 = t1473 * t1481 - t1475 * t1536;
t1413 = t1437 * t1575 + t1457 * t1473 - t1460 * t1547;
t1514 = t1413 * t1563;
t1513 = t1413 * t1562;
t1438 = t1473 * t1483 - t1475 * t1534;
t1414 = t1438 * t1574 + t1458 * t1473 - t1461 * t1547;
t1512 = t1414 * t1559;
t1511 = t1414 * t1558;
t1439 = t1473 * t1485 - t1475 * t1532;
t1415 = t1439 * t1573 + t1459 * t1473 - t1462 * t1547;
t1510 = t1415 * t1555;
t1509 = t1415 * t1554;
t1440 = t1473 * t1487 + t1475 * t1542;
t1428 = t1480 * t1440 + t1475 * t1523;
t1508 = t1428 * t1563;
t1507 = t1428 * t1562;
t1441 = t1473 * t1489 + t1475 * t1540;
t1429 = t1482 * t1441 + t1475 * t1522;
t1506 = t1429 * t1559;
t1505 = t1429 * t1558;
t1442 = t1473 * t1491 + t1475 * t1538;
t1430 = t1442 * t1484 + t1475 * t1521;
t1504 = t1430 * t1555;
t1503 = t1430 * t1554;
t1502 = t1428 * t1520;
t1501 = t1428 * t1519;
t1500 = t1429 * t1518;
t1499 = t1429 * t1517;
t1498 = t1430 * t1516;
t1497 = t1430 * t1515;
t1496 = pkin(3) * t1531 - t1457 * t1476;
t1495 = pkin(3) * t1528 - t1458 * t1476;
t1494 = pkin(3) * t1549 - t1459 * t1476;
t1493 = 0.1e1 / pkin(3);
t1448 = t1476 * t1525 + t1521;
t1447 = t1476 * t1527 + t1522;
t1446 = t1476 * t1530 + t1523;
t1424 = t1475 * t1462 + t1494 * t1473;
t1423 = t1475 * t1461 + t1495 * t1473;
t1422 = t1475 * t1460 + t1496 * t1473;
t1412 = -g(3) * t1442 - t1456 * t1436 + t1453 * t1548;
t1411 = -t1453 * t1474 * t1491 - g(3) * t1439 + t1456 * t1445;
t1410 = -g(3) * t1441 - t1455 * t1435 + t1452 * t1550;
t1409 = -t1452 * t1474 * t1489 - g(3) * t1438 + t1455 * t1444;
t1408 = -g(3) * t1440 - t1454 * t1434 + t1451 * t1551;
t1407 = -t1451 * t1474 * t1487 - g(3) * t1437 + t1454 * t1443;
t1403 = (-t1448 * t1473 + t1475 * t1524) * t1456 - g(3) * (t1448 * t1475 + t1473 * t1524) - t1453 * (-t1474 * t1525 + t1533);
t1402 = (-t1447 * t1473 + t1475 * t1526) * t1455 - g(3) * (t1447 * t1475 + t1473 * t1526) - t1452 * (-t1474 * t1527 + t1535);
t1401 = (-t1446 * t1473 + t1475 * t1529) * t1454 - g(3) * (t1446 * t1475 + t1473 * t1529) - t1451 * (-t1474 * t1530 + t1537);
t1400 = ((-t1483 * t1535 + t1528) * t1473 + t1488 * t1545) * t1455 + g(3) * (-t1441 * t1488 + t1475 * t1528) + t1452 * (t1483 * t1522 + t1541);
t1399 = ((-t1485 * t1533 + t1549) * t1473 + t1490 * t1544) * t1456 + g(3) * (-t1442 * t1490 + t1475 * t1549) + t1453 * (t1485 * t1521 + t1539);
t1398 = t1454 * ((-t1481 * t1537 + t1531) * t1473 + t1486 * t1546) + g(3) * (-t1440 * t1486 + t1475 * t1531) + t1451 * (t1481 * t1523 + t1543);
t1 = [-(-(t1436 * t1469 - t1466 * t1548) * t1576 + (t1424 * t1469 + t1466 * t1433) * t1490 + (t1476 * t1466 + t1469 * t1552) * t1579) * t1556 - (-(t1435 * t1468 - t1465 * t1550) * t1577 + (t1423 * t1468 + t1465 * t1432) * t1488 + (t1476 * t1465 + t1468 * t1552) * t1580) * t1560 - (-(t1434 * t1467 - t1464 * t1551) * t1578 + (t1422 * t1467 + t1464 * t1431) * t1486 + (t1476 * t1464 + t1467 * t1552) * t1581) * t1564, 0, -t1407 * t1507 - t1409 * t1505 - t1411 * t1503, -t1408 * t1507 - t1410 * t1505 - t1412 * t1503, 0, 0, 0, 0, 0, -t1467 * t1501 - t1468 * t1499 - t1469 * t1497 + (t1401 * t1513 + t1402 * t1511 + t1403 * t1509) * t1493, t1467 * t1502 + t1468 * t1500 + t1469 * t1498 + (t1398 * t1513 + t1399 * t1509 + t1400 * t1511) * t1493, -g(1); -((t1436 * t1466 + t1469 * t1548) * t1576 + (-t1424 * t1466 + t1469 * t1433) * t1490 + (-t1466 * t1552 + t1469 * t1476) * t1579) * t1556 - ((t1435 * t1465 + t1468 * t1550) * t1577 + (-t1423 * t1465 + t1468 * t1432) * t1488 + (-t1465 * t1552 + t1468 * t1476) * t1580) * t1560 - ((t1434 * t1464 + t1467 * t1551) * t1578 + (-t1422 * t1464 + t1467 * t1431) * t1486 + (-t1464 * t1552 + t1467 * t1476) * t1581) * t1564, 0, t1407 * t1508 + t1409 * t1506 + t1411 * t1504, t1408 * t1508 + t1410 * t1506 + t1412 * t1504, 0, 0, 0, 0, 0, t1464 * t1501 + t1465 * t1499 + t1466 * t1497 + (-t1401 * t1514 - t1402 * t1512 - t1403 * t1510) * t1493, -t1464 * t1502 - t1465 * t1500 - t1466 * t1498 + (-t1398 * t1514 - t1399 * t1510 - t1400 * t1512) * t1493, -g(2); -(-t1442 * t1576 - t1462 * t1473 * t1490 + (pkin(2) * t1549 + t1494 * t1490) * t1475) * t1556 - (-t1441 * t1577 - t1461 * t1473 * t1488 + (pkin(2) * t1528 + t1495 * t1488) * t1475) * t1560 - (-t1440 * t1578 - t1460 * t1473 * t1486 + (pkin(2) * t1531 + t1496 * t1486) * t1475) * t1564, 0, t1407 * t1565 + t1409 * t1561 + t1411 * t1557, t1408 * t1565 + t1410 * t1561 + t1412 * t1557, 0, 0, 0, 0, 0, t1426 * t1519 + t1427 * t1517 + t1425 * t1515 + (t1401 * t1568 + t1402 * t1567 + t1403 * t1566) * t1493, -t1426 * t1520 - t1427 * t1518 - t1425 * t1516 + (t1398 * t1568 + t1399 * t1566 + t1400 * t1567) * t1493, -g(3);];
tau_reg  = t1;
