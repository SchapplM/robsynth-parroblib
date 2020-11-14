% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRRR8V2G2A0
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
% Datum: 2020-08-06 17:50
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRRR8V2G2A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_regmin: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:50:16
% EndTime: 2020-08-06 17:50:17
% DurationCPUTime: 1.63s
% Computational Cost: add. (1086->210), mult. (2529->426), div. (81->4), fcn. (2667->22), ass. (0->179)
t1501 = sin(qJ(2,3));
t1507 = cos(qJ(2,3));
t1512 = pkin(7) + pkin(6);
t1475 = pkin(2) * t1501 - t1512 * t1507;
t1494 = sin(pkin(4));
t1496 = cos(pkin(4));
t1500 = sin(qJ(3,3));
t1573 = t1496 * t1500;
t1446 = pkin(3) * t1573 + t1475 * t1494;
t1506 = cos(qJ(3,3));
t1577 = t1494 * t1501;
t1600 = pkin(3) * t1506 ^ 2;
t1431 = 0.1e1 / (pkin(2) * t1573 + t1446 * t1506 + t1577 * t1600);
t1495 = cos(pkin(8));
t1574 = t1495 * t1496;
t1467 = -t1494 * g(1) + g(2) * t1574;
t1468 = g(1) * t1574 + t1494 * g(2);
t1497 = legFrame(3,2);
t1484 = sin(t1497);
t1487 = cos(t1497);
t1472 = t1487 * g(1) - t1484 * g(2);
t1493 = sin(pkin(8));
t1579 = t1493 * t1496;
t1481 = g(3) * t1579;
t1483 = g(3) * t1495;
t1609 = t1431 * ((t1467 * t1484 - t1468 * t1487 + t1481) * t1507 + t1501 * (t1472 * t1493 + t1483));
t1503 = sin(qJ(2,2));
t1509 = cos(qJ(2,2));
t1476 = pkin(2) * t1503 - t1512 * t1509;
t1502 = sin(qJ(3,2));
t1571 = t1496 * t1502;
t1447 = pkin(3) * t1571 + t1476 * t1494;
t1508 = cos(qJ(3,2));
t1576 = t1494 * t1503;
t1599 = pkin(3) * t1508 ^ 2;
t1432 = 0.1e1 / (pkin(2) * t1571 + t1447 * t1508 + t1576 * t1599);
t1498 = legFrame(2,2);
t1485 = sin(t1498);
t1488 = cos(t1498);
t1473 = t1488 * g(1) - t1485 * g(2);
t1608 = t1432 * ((t1467 * t1485 - t1468 * t1488 + t1481) * t1509 + t1503 * (t1473 * t1493 + t1483));
t1505 = sin(qJ(2,1));
t1511 = cos(qJ(2,1));
t1477 = pkin(2) * t1505 - t1512 * t1511;
t1504 = sin(qJ(3,1));
t1569 = t1496 * t1504;
t1448 = pkin(3) * t1569 + t1477 * t1494;
t1510 = cos(qJ(3,1));
t1575 = t1494 * t1505;
t1598 = pkin(3) * t1510 ^ 2;
t1433 = 0.1e1 / (pkin(2) * t1569 + t1448 * t1510 + t1575 * t1598);
t1499 = legFrame(1,2);
t1486 = sin(t1499);
t1489 = cos(t1499);
t1474 = t1489 * g(1) - t1486 * g(2);
t1607 = t1433 * ((t1467 * t1486 - t1468 * t1489 + t1481) * t1511 + t1505 * (t1474 * t1493 + t1483));
t1603 = pkin(2) * t1500;
t1602 = pkin(2) * t1502;
t1601 = pkin(2) * t1504;
t1597 = pkin(3) * t1506;
t1596 = pkin(3) * t1508;
t1595 = pkin(3) * t1510;
t1566 = t1496 * t1507;
t1452 = t1493 * t1501 - t1495 * t1566;
t1478 = pkin(2) * t1507 + t1501 * t1512;
t1594 = (t1452 * t1597 + t1493 * t1475 - t1478 * t1574) * t1431;
t1564 = t1496 * t1509;
t1453 = t1493 * t1503 - t1495 * t1564;
t1479 = pkin(2) * t1509 + t1503 * t1512;
t1593 = (t1453 * t1596 + t1493 * t1476 - t1479 * t1574) * t1432;
t1562 = t1496 * t1511;
t1454 = t1493 * t1505 - t1495 * t1562;
t1480 = pkin(2) * t1511 + t1505 * t1512;
t1592 = (t1454 * t1595 + t1493 * t1477 - t1480 * t1574) * t1433;
t1572 = t1496 * t1501;
t1455 = t1493 * t1507 + t1495 * t1572;
t1549 = t1506 * t1494;
t1440 = -t1500 * t1455 - t1495 * t1549;
t1591 = t1431 * t1440;
t1469 = t1484 * g(1) + t1487 * g(2);
t1590 = t1431 * t1469;
t1589 = t1431 * t1484;
t1588 = t1431 * t1487;
t1570 = t1496 * t1503;
t1456 = t1493 * t1509 + t1495 * t1570;
t1547 = t1508 * t1494;
t1441 = -t1502 * t1456 - t1495 * t1547;
t1587 = t1432 * t1441;
t1470 = t1485 * g(1) + t1488 * g(2);
t1586 = t1432 * t1470;
t1585 = t1432 * t1485;
t1584 = t1432 * t1488;
t1568 = t1496 * t1505;
t1457 = t1493 * t1511 + t1495 * t1568;
t1545 = t1510 * t1494;
t1442 = -t1457 * t1504 - t1495 * t1545;
t1583 = t1433 * t1442;
t1471 = t1486 * g(1) + t1489 * g(2);
t1582 = t1433 * t1471;
t1581 = t1433 * t1486;
t1580 = t1433 * t1489;
t1578 = t1494 * t1495;
t1567 = t1496 * t1506;
t1565 = t1496 * t1508;
t1563 = t1496 * t1510;
t1561 = t1500 * t1494;
t1560 = t1500 * t1501;
t1559 = t1500 * t1507;
t1557 = t1502 * t1494;
t1556 = t1502 * t1503;
t1555 = t1502 * t1509;
t1553 = t1504 * t1494;
t1552 = t1504 * t1505;
t1551 = t1504 * t1511;
t1548 = t1506 * t1507;
t1546 = t1508 * t1509;
t1544 = t1510 * t1511;
t1543 = t1500 * t1609;
t1542 = t1502 * t1608;
t1541 = t1504 * t1607;
t1540 = t1506 * t1609;
t1539 = t1508 * t1608;
t1538 = t1510 * t1607;
t1458 = t1493 * t1566 + t1495 * t1501;
t1428 = t1458 * t1597 + t1475 * t1495 + t1478 * t1579;
t1537 = t1428 * t1589;
t1536 = t1428 * t1588;
t1459 = t1493 * t1564 + t1495 * t1503;
t1429 = t1459 * t1596 + t1476 * t1495 + t1479 * t1579;
t1535 = t1429 * t1585;
t1534 = t1429 * t1584;
t1460 = t1493 * t1562 + t1495 * t1505;
t1430 = t1460 * t1595 + t1477 * t1495 + t1480 * t1579;
t1533 = t1430 * t1581;
t1532 = t1430 * t1580;
t1449 = t1493 * t1572 - t1495 * t1507;
t1437 = t1500 * t1449 + t1493 * t1549;
t1531 = t1437 * t1589;
t1530 = t1437 * t1588;
t1450 = t1493 * t1570 - t1495 * t1509;
t1438 = t1502 * t1450 + t1493 * t1547;
t1529 = t1438 * t1585;
t1528 = t1438 * t1584;
t1451 = t1493 * t1568 - t1495 * t1511;
t1439 = t1504 * t1451 + t1493 * t1545;
t1527 = t1439 * t1581;
t1526 = t1439 * t1580;
t1525 = t1437 * t1543;
t1524 = t1438 * t1542;
t1523 = t1439 * t1541;
t1522 = t1437 * t1540;
t1521 = t1438 * t1539;
t1520 = t1439 * t1538;
t1519 = pkin(3) * t1561 - t1475 * t1496;
t1518 = pkin(3) * t1557 - t1476 * t1496;
t1517 = pkin(3) * t1553 - t1477 * t1496;
t1513 = 0.1e1 / pkin(3);
t1466 = t1505 * t1563 - t1553;
t1465 = t1503 * t1565 - t1557;
t1464 = t1501 * t1567 - t1561;
t1463 = t1496 * t1552 + t1545;
t1462 = t1496 * t1556 + t1547;
t1461 = t1496 * t1560 + t1549;
t1436 = -t1493 * t1480 + t1495 * t1517;
t1435 = -t1493 * t1479 + t1495 * t1518;
t1434 = -t1493 * t1478 + t1495 * t1519;
t1424 = -g(3) * t1451 + t1474 * t1457 + t1471 * t1575;
t1423 = -g(3) * t1450 + t1473 * t1456 + t1470 * t1576;
t1422 = -g(3) * t1449 + t1472 * t1455 + t1469 * t1577;
t1421 = -t1471 * t1494 * t1511 + g(3) * t1460 + t1474 * t1454;
t1420 = -t1470 * t1494 * t1509 + g(3) * t1459 + t1473 * t1453;
t1419 = -t1469 * t1494 * t1507 + g(3) * t1458 + t1472 * t1452;
t1412 = g(3) * (-t1461 * t1493 + t1495 * t1559) + t1472 * (t1461 * t1495 + t1493 * t1559) - t1469 * (-t1494 * t1560 + t1567);
t1411 = (-t1466 * t1493 + t1495 * t1544) * g(3) + (t1466 * t1495 + t1493 * t1544) * t1474 + t1471 * (t1505 * t1545 + t1569);
t1410 = (-t1465 * t1493 + t1495 * t1546) * g(3) + (t1465 * t1495 + t1493 * t1546) * t1473 + t1470 * (t1503 * t1547 + t1571);
t1409 = (-t1464 * t1493 + t1495 * t1548) * g(3) + (t1464 * t1495 + t1493 * t1548) * t1472 + t1469 * (t1501 * t1549 + t1573);
t1408 = (-t1463 * t1493 + t1495 * t1551) * g(3) + t1474 * (t1463 * t1495 + t1493 * t1551) - t1471 * (-t1494 * t1552 + t1563);
t1407 = (-t1462 * t1493 + t1495 * t1555) * g(3) + t1473 * (t1462 * t1495 + t1493 * t1555) - t1470 * (-t1494 * t1556 + t1565);
t1 = [-((t1457 * t1489 + t1486 * t1575) * t1598 + (-t1436 * t1489 + t1448 * t1486) * t1510 + (t1496 * t1486 - t1489 * t1578) * t1601) * t1582 - ((t1456 * t1488 + t1485 * t1576) * t1599 + (-t1435 * t1488 + t1447 * t1485) * t1508 + (t1496 * t1485 - t1488 * t1578) * t1602) * t1586 - ((t1455 * t1487 + t1484 * t1577) * t1600 + (-t1434 * t1487 + t1446 * t1484) * t1506 + (t1496 * t1484 - t1487 * t1578) * t1603) * t1590, 0, -t1419 * t1530 - t1420 * t1528 - t1421 * t1526, -t1422 * t1530 - t1423 * t1528 - t1424 * t1526, 0, 0, 0, 0, 0, -t1487 * t1522 - t1488 * t1521 - t1489 * t1520 + (-t1407 * t1534 - t1408 * t1532 - t1412 * t1536) * t1513, t1487 * t1525 + t1488 * t1524 + t1489 * t1523 + (-t1409 * t1536 - t1410 * t1534 - t1411 * t1532) * t1513, -g(1); -(-(t1457 * t1486 - t1489 * t1575) * t1598 + (t1436 * t1486 + t1489 * t1448) * t1510 + (t1486 * t1578 + t1489 * t1496) * t1601) * t1582 - (-(t1456 * t1485 - t1488 * t1576) * t1599 + (t1435 * t1485 + t1488 * t1447) * t1508 + (t1485 * t1578 + t1488 * t1496) * t1602) * t1586 - (-(t1455 * t1484 - t1487 * t1577) * t1600 + (t1434 * t1484 + t1487 * t1446) * t1506 + (t1484 * t1578 + t1487 * t1496) * t1603) * t1590, 0, t1419 * t1531 + t1420 * t1529 + t1421 * t1527, t1422 * t1531 + t1423 * t1529 + t1424 * t1527, 0, 0, 0, 0, 0, t1484 * t1522 + t1485 * t1521 + t1486 * t1520 + (t1407 * t1535 + t1408 * t1533 + t1412 * t1537) * t1513, -t1484 * t1525 - t1485 * t1524 - t1486 * t1523 + (t1409 * t1537 + t1410 * t1535 + t1411 * t1533) * t1513, -g(2); -(-t1451 * t1598 + t1495 * t1480 * t1510 + (pkin(2) * t1553 + t1510 * t1517) * t1493) * t1582 - (-t1450 * t1599 + t1495 * t1479 * t1508 + (pkin(2) * t1557 + t1508 * t1518) * t1493) * t1586 - (-t1449 * t1600 + t1495 * t1478 * t1506 + (pkin(2) * t1561 + t1506 * t1519) * t1493) * t1590, 0, t1419 * t1591 + t1420 * t1587 + t1421 * t1583, t1422 * t1591 + t1423 * t1587 + t1424 * t1583, 0, 0, 0, 0, 0, t1440 * t1540 + t1441 * t1539 + t1442 * t1538 + (t1407 * t1593 + t1408 * t1592 + t1412 * t1594) * t1513, -t1440 * t1543 - t1441 * t1542 - t1442 * t1541 + (t1409 * t1594 + t1410 * t1593 + t1411 * t1592) * t1513, -g(3);];
tau_reg  = t1;
