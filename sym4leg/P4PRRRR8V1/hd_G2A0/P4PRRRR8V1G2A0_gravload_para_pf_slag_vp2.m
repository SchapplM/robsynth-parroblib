% Calculate Gravitation load for parallel robot
% P4PRRRR8V1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [4x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:06
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P4PRRRR8V1G2A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V1G2A0_gravload_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V1G2A0_gravload_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P4PRRRR8V1G2A0_gravload_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V1G2A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR8V1G2A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRRR8V1G2A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V1G2A0_gravload_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V1G2A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:03:18
% EndTime: 2020-08-07 11:03:19
% DurationCPUTime: 1.54s
% Computational Cost: add. (1315->209), mult. (2993->397), div. (120->9), fcn. (2654->30), ass. (0->201)
t1552 = sin(qJ(2,1));
t1546 = legFrame(1,2);
t1522 = sin(t1546);
t1526 = cos(t1546);
t1502 = t1526 * g(1) - t1522 * g(2);
t1536 = cos(pkin(6));
t1623 = t1536 * t1502;
t1534 = sin(pkin(6));
t1652 = t1534 * g(3);
t1484 = t1623 - t1652;
t1498 = t1522 * g(1) + t1526 * g(2);
t1535 = sin(pkin(3));
t1537 = cos(pkin(3));
t1587 = t1484 * t1537 + t1498 * t1535;
t1669 = t1587 * t1552;
t1550 = sin(qJ(2,2));
t1545 = legFrame(2,2);
t1521 = sin(t1545);
t1525 = cos(t1545);
t1501 = t1525 * g(1) - t1521 * g(2);
t1624 = t1536 * t1501;
t1483 = t1624 - t1652;
t1497 = t1521 * g(1) + t1525 * g(2);
t1588 = t1483 * t1537 + t1497 * t1535;
t1668 = t1588 * t1550;
t1548 = sin(qJ(2,3));
t1544 = legFrame(3,2);
t1520 = sin(t1544);
t1524 = cos(t1544);
t1500 = t1524 * g(1) - t1520 * g(2);
t1625 = t1536 * t1500;
t1482 = t1625 - t1652;
t1496 = t1520 * g(1) + t1524 * g(2);
t1589 = t1482 * t1537 + t1496 * t1535;
t1667 = t1589 * t1548;
t1539 = sin(qJ(2,4));
t1543 = legFrame(4,2);
t1519 = sin(t1543);
t1523 = cos(t1543);
t1499 = t1523 * g(1) - t1519 * g(2);
t1626 = t1536 * t1499;
t1481 = t1626 - t1652;
t1495 = t1519 * g(1) + t1523 * g(2);
t1590 = t1481 * t1537 + t1495 * t1535;
t1666 = t1590 * t1539;
t1541 = cos(qJ(2,4));
t1540 = cos(qJ(3,4));
t1660 = pkin(2) * t1540;
t1493 = -t1541 * pkin(5) + t1539 * t1660;
t1538 = sin(qJ(3,4));
t1661 = pkin(2) * t1538;
t1474 = t1493 * t1535 + t1537 * t1661;
t1665 = 0.1e1 / t1474;
t1554 = cos(qJ(2,3));
t1553 = cos(qJ(3,3));
t1656 = pkin(2) * t1553;
t1503 = -t1554 * pkin(5) + t1548 * t1656;
t1547 = sin(qJ(3,3));
t1659 = pkin(2) * t1547;
t1478 = t1503 * t1535 + t1537 * t1659;
t1664 = 0.1e1 / t1478;
t1556 = cos(qJ(2,2));
t1555 = cos(qJ(3,2));
t1655 = pkin(2) * t1555;
t1504 = -t1556 * pkin(5) + t1550 * t1655;
t1549 = sin(qJ(3,2));
t1658 = pkin(2) * t1549;
t1479 = t1504 * t1535 + t1537 * t1658;
t1663 = 0.1e1 / t1479;
t1558 = cos(qJ(2,1));
t1557 = cos(qJ(3,1));
t1654 = pkin(2) * t1557;
t1505 = -t1558 * pkin(5) + t1552 * t1654;
t1551 = sin(qJ(3,1));
t1657 = pkin(2) * t1551;
t1480 = t1505 * t1535 + t1537 * t1657;
t1662 = 0.1e1 / t1480;
t1653 = g(3) * t1536;
t1651 = t1665 / t1540;
t1650 = t1664 / t1553;
t1649 = t1663 / t1555;
t1648 = t1662 / t1557;
t1647 = t1665 * t1495;
t1646 = t1664 * t1496;
t1645 = t1663 * t1497;
t1644 = t1662 * t1498;
t1638 = t1495 * t1537;
t1636 = t1496 * t1537;
t1634 = t1497 * t1537;
t1632 = t1498 * t1537;
t1631 = mrSges(3,2) * t1652 * t1535;
t1630 = t1534 * t1541;
t1629 = t1534 * t1554;
t1628 = t1534 * t1556;
t1627 = t1534 * t1558;
t1622 = t1537 * t1539;
t1621 = t1537 * t1541;
t1620 = t1537 * t1548;
t1619 = t1537 * t1550;
t1618 = t1537 * t1552;
t1617 = t1537 * t1554;
t1616 = t1537 * t1556;
t1615 = t1537 * t1558;
t1614 = t1538 * t1541;
t1613 = t1547 * t1554;
t1612 = t1549 * t1556;
t1611 = t1551 * t1558;
t1594 = t1499 * t1534 + t1653;
t1574 = t1594 * t1541 + t1666;
t1441 = ((t1481 * t1535 - t1638) * mrSges(3,1) + t1574 * mrSges(3,2)) * t1540 + (t1631 + (-t1535 * t1626 + t1638) * mrSges(3,2) + t1574 * mrSges(3,1)) * t1538;
t1610 = t1441 * t1651;
t1593 = t1500 * t1534 + t1653;
t1573 = t1593 * t1554 + t1667;
t1442 = ((t1482 * t1535 - t1636) * mrSges(3,1) + t1573 * mrSges(3,2)) * t1553 + (t1631 + (-t1535 * t1625 + t1636) * mrSges(3,2) + t1573 * mrSges(3,1)) * t1547;
t1609 = t1442 * t1650;
t1592 = t1501 * t1534 + t1653;
t1572 = t1592 * t1556 + t1668;
t1443 = ((t1483 * t1535 - t1634) * mrSges(3,1) + t1572 * mrSges(3,2)) * t1555 + (t1631 + (-t1535 * t1624 + t1634) * mrSges(3,2) + t1572 * mrSges(3,1)) * t1549;
t1608 = t1443 * t1649;
t1591 = t1502 * t1534 + t1653;
t1571 = t1591 * t1558 + t1669;
t1444 = ((t1484 * t1535 - t1632) * mrSges(3,1) + t1571 * mrSges(3,2)) * t1557 + (t1631 + (-t1535 * t1623 + t1632) * mrSges(3,2) + t1571 * mrSges(3,1)) * t1551;
t1607 = t1444 * t1648;
t1542 = mrSges(2,2) - mrSges(3,3);
t1513 = t1542 * t1653;
t1445 = t1513 * t1541 + (t1499 * t1630 + t1666) * t1542 + (t1594 * t1539 - t1590 * t1541) * (t1540 * mrSges(3,1) - t1538 * mrSges(3,2) + mrSges(2,1));
t1606 = t1445 * t1651;
t1446 = t1513 * t1554 + (t1500 * t1629 + t1667) * t1542 + (t1593 * t1548 - t1589 * t1554) * (t1553 * mrSges(3,1) - t1547 * mrSges(3,2) + mrSges(2,1));
t1605 = t1446 * t1650;
t1447 = t1513 * t1556 + (t1501 * t1628 + t1668) * t1542 + (t1592 * t1550 - t1588 * t1556) * (t1555 * mrSges(3,1) - t1549 * mrSges(3,2) + mrSges(2,1));
t1604 = t1447 * t1649;
t1448 = t1513 * t1558 + (t1502 * t1627 + t1669) * t1542 + (t1591 * t1552 - t1587 * t1558) * (t1557 * mrSges(3,1) - t1551 * mrSges(3,2) + mrSges(2,1));
t1603 = t1448 * t1648;
t1461 = (t1534 * t1621 + t1536 * t1539) * t1660 + (t1534 * t1622 - t1536 * t1541) * pkin(5);
t1602 = t1461 * t1610;
t1462 = (t1534 * t1617 + t1536 * t1548) * t1656 + (t1534 * t1620 - t1536 * t1554) * pkin(5);
t1601 = t1462 * t1609;
t1463 = (t1534 * t1616 + t1536 * t1550) * t1655 + (t1534 * t1619 - t1536 * t1556) * pkin(5);
t1600 = t1463 * t1608;
t1464 = (t1534 * t1615 + t1536 * t1552) * t1654 + (t1534 * t1618 - t1536 * t1558) * pkin(5);
t1599 = t1464 * t1607;
t1582 = t1535 * t1540 + t1538 * t1622;
t1465 = t1582 * t1534 - t1536 * t1614;
t1598 = t1465 * t1606;
t1581 = t1535 * t1553 + t1547 * t1620;
t1466 = t1581 * t1534 - t1536 * t1613;
t1597 = t1466 * t1605;
t1580 = t1535 * t1555 + t1549 * t1619;
t1467 = t1580 * t1534 - t1536 * t1612;
t1596 = t1467 * t1604;
t1579 = t1535 * t1557 + t1551 * t1618;
t1468 = t1579 * t1534 - t1536 * t1611;
t1595 = t1468 * t1603;
t1586 = -t1493 * t1537 + t1535 * t1661;
t1585 = -t1503 * t1537 + t1535 * t1659;
t1584 = -t1504 * t1537 + t1535 * t1658;
t1583 = -t1505 * t1537 + t1535 * t1657;
t1559 = xP(4);
t1527 = sin(t1559);
t1528 = cos(t1559);
t1562 = koppelP(4,2);
t1566 = koppelP(4,1);
t1485 = -t1527 * t1566 - t1528 * t1562;
t1489 = -t1527 * t1562 + t1528 * t1566;
t1578 = (-t1485 * t1523 + t1489 * t1519) * t1651;
t1563 = koppelP(3,2);
t1567 = koppelP(3,1);
t1486 = -t1527 * t1567 - t1528 * t1563;
t1490 = -t1527 * t1563 + t1528 * t1567;
t1577 = (-t1486 * t1524 + t1490 * t1520) * t1650;
t1564 = koppelP(2,2);
t1568 = koppelP(2,1);
t1487 = -t1527 * t1568 - t1528 * t1564;
t1491 = -t1527 * t1564 + t1528 * t1568;
t1576 = (-t1487 * t1525 + t1491 * t1521) * t1649;
t1565 = koppelP(1,2);
t1569 = koppelP(1,1);
t1488 = -t1527 * t1569 - t1528 * t1565;
t1492 = -t1527 * t1565 + t1528 * t1569;
t1575 = (-t1488 * t1526 + t1492 * t1522) * t1648;
t1570 = 0.1e1 / pkin(2);
t1561 = mrSges(4,1);
t1560 = mrSges(4,2);
t1530 = m(1) + m(2) + m(3);
t1508 = pkin(5) * t1552 + t1558 * t1654;
t1507 = pkin(5) * t1550 + t1556 * t1655;
t1506 = pkin(5) * t1548 + t1554 * t1656;
t1494 = pkin(5) * t1539 + t1541 * t1660;
t1460 = -t1534 * t1508 + t1583 * t1536;
t1459 = -t1534 * t1507 + t1584 * t1536;
t1458 = -t1534 * t1506 + t1585 * t1536;
t1457 = -t1534 * t1494 + t1586 * t1536;
t1456 = -t1460 * t1526 + t1522 * t1480;
t1455 = t1460 * t1522 + t1526 * t1480;
t1454 = -t1459 * t1525 + t1521 * t1479;
t1453 = t1459 * t1521 + t1525 * t1479;
t1452 = -t1458 * t1524 + t1520 * t1478;
t1451 = t1458 * t1520 + t1524 * t1478;
t1450 = -t1457 * t1523 + t1519 * t1474;
t1449 = t1457 * t1519 + t1523 * t1474;
t1 = [-t1523 * t1598 - t1524 * t1597 - t1525 * t1596 - t1526 * t1595 - g(1) * m(4) + (-t1523 * t1602 - t1524 * t1601 - t1525 * t1600 - t1526 * t1599) * t1570 + (-t1450 * t1647 - t1452 * t1646 - t1454 * t1645 - t1456 * t1644) * t1530; t1519 * t1598 + t1520 * t1597 + t1521 * t1596 + t1522 * t1595 - g(2) * m(4) + (t1519 * t1602 + t1520 * t1601 + t1521 * t1600 + t1522 * t1599) * t1570 + (-t1449 * t1647 - t1451 * t1646 - t1453 * t1645 - t1455 * t1644) * t1530; (-t1534 * t1611 - t1579 * t1536) * t1603 + (-t1534 * t1612 - t1580 * t1536) * t1604 + (-t1534 * t1613 - t1581 * t1536) * t1605 + (-t1534 * t1614 - t1582 * t1536) * t1606 - g(3) * m(4) + ((-(-t1534 * t1552 + t1536 * t1615) * t1654 - pkin(5) * (t1536 * t1618 + t1627)) * t1607 + (-(-t1534 * t1550 + t1536 * t1616) * t1655 - pkin(5) * (t1536 * t1619 + t1628)) * t1608 + (-(-t1534 * t1548 + t1536 * t1617) * t1656 - pkin(5) * (t1536 * t1620 + t1629)) * t1609 + (-(-t1534 * t1539 + t1536 * t1621) * t1660 - pkin(5) * (t1536 * t1622 + t1630)) * t1610) * t1570 + (-(t1536 * t1508 + t1583 * t1534) * t1644 - (t1536 * t1507 + t1584 * t1534) * t1645 - (t1536 * t1506 + t1585 * t1534) * t1646 - (t1536 * t1494 + t1586 * t1534) * t1647) * t1530; -(-g(1) * t1561 - g(2) * t1560) * t1527 + t1528 * (g(1) * t1560 - g(2) * t1561) + t1448 * t1468 * t1575 + t1447 * t1467 * t1576 + t1446 * t1466 * t1577 + t1445 * t1465 * t1578 + (-(t1455 * t1492 + t1456 * t1488) * t1644 - (t1453 * t1491 + t1454 * t1487) * t1645 - (t1451 * t1490 + t1452 * t1486) * t1646 - (t1449 * t1489 + t1450 * t1485) * t1647) * t1530 + (t1441 * t1461 * t1578 + t1442 * t1462 * t1577 + t1443 * t1463 * t1576 + t1464 * t1444 * t1575) * t1570;];
taugX  = t1;
