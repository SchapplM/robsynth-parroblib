% Calculate Gravitation load for parallel robot
% P4PRRRR8V2G3A0
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2020-08-07 11:30
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P4PRRRR8V2G3A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G3A0_gravload_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G3A0_gravload_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G3A0_gravload_para_pf_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V2G3A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR8V2G3A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRRR8V2G3A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G3A0_gravload_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G3A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:26:11
% EndTime: 2020-08-07 11:26:14
% DurationCPUTime: 2.73s
% Computational Cost: add. (1851->264), mult. (3705->500), div. (80->5), fcn. (3150->30), ass. (0->210)
t1558 = sin(qJ(2,4));
t1561 = legFrame(4,2);
t1538 = sin(t1561);
t1542 = cos(t1561);
t1518 = t1538 * g(1) + t1542 * g(2);
t1554 = sin(pkin(4));
t1556 = cos(pkin(4));
t1522 = t1542 * g(1) - t1538 * g(2);
t1555 = cos(pkin(8));
t1537 = g(3) * t1555;
t1553 = sin(pkin(8));
t1614 = -t1522 * t1553 - t1537;
t1701 = t1518 * t1554 + t1614 * t1556;
t1708 = t1701 * t1558;
t1566 = sin(qJ(2,3));
t1562 = legFrame(3,2);
t1539 = sin(t1562);
t1543 = cos(t1562);
t1519 = t1539 * g(1) + t1543 * g(2);
t1523 = t1543 * g(1) - t1539 * g(2);
t1612 = -t1523 * t1553 - t1537;
t1702 = t1519 * t1554 + t1612 * t1556;
t1707 = t1702 * t1566;
t1568 = sin(qJ(2,2));
t1563 = legFrame(2,2);
t1540 = sin(t1563);
t1544 = cos(t1563);
t1520 = t1540 * g(1) + t1544 * g(2);
t1524 = t1544 * g(1) - t1540 * g(2);
t1610 = -t1524 * t1553 - t1537;
t1703 = t1520 * t1554 + t1610 * t1556;
t1706 = t1703 * t1568;
t1570 = sin(qJ(2,1));
t1564 = legFrame(1,2);
t1541 = sin(t1564);
t1545 = cos(t1564);
t1521 = t1541 * g(1) + t1545 * g(2);
t1525 = t1545 * g(1) - t1541 * g(2);
t1608 = -t1525 * t1553 - t1537;
t1704 = t1521 * t1554 + t1608 * t1556;
t1705 = t1704 * t1570;
t1557 = sin(qJ(3,4));
t1696 = pkin(2) * t1557;
t1565 = sin(qJ(3,3));
t1695 = pkin(2) * t1565;
t1567 = sin(qJ(3,2));
t1694 = pkin(2) * t1567;
t1569 = sin(qJ(3,1));
t1693 = pkin(2) * t1569;
t1559 = cos(qJ(3,4));
t1692 = pkin(3) * t1559 ^ 2;
t1571 = cos(qJ(3,3));
t1691 = pkin(3) * t1571 ^ 2;
t1573 = cos(qJ(3,2));
t1690 = pkin(3) * t1573 ^ 2;
t1575 = cos(qJ(3,1));
t1689 = pkin(3) * t1575 ^ 2;
t1688 = pkin(3) * t1559;
t1687 = pkin(3) * t1571;
t1686 = pkin(3) * t1573;
t1685 = pkin(3) * t1575;
t1684 = g(3) * t1553;
t1683 = m(3) * pkin(2) + mrSges(2,1);
t1560 = cos(qJ(2,4));
t1577 = pkin(7) + pkin(6);
t1526 = pkin(2) * t1558 - t1577 * t1560;
t1634 = t1556 * t1557;
t1490 = pkin(3) * t1634 + t1526 * t1554;
t1650 = t1554 * t1558;
t1478 = 0.1e1 / (pkin(2) * t1634 + t1490 * t1559 + t1650 * t1692);
t1613 = t1522 * t1555 - t1684;
t1594 = t1613 * t1560 + t1708;
t1653 = t1553 * t1554;
t1654 = mrSges(3,2) * t1554 * t1537;
t1668 = t1518 * t1556;
t1682 = (((t1614 * t1554 - t1668) * mrSges(3,1) + t1594 * mrSges(3,2)) * t1559 + t1557 * (t1654 + (t1522 * t1653 + t1668) * mrSges(3,2) + t1594 * mrSges(3,1))) * t1478;
t1572 = cos(qJ(2,3));
t1528 = pkin(2) * t1566 - t1577 * t1572;
t1631 = t1556 * t1565;
t1491 = pkin(3) * t1631 + t1528 * t1554;
t1647 = t1554 * t1566;
t1479 = 0.1e1 / (pkin(2) * t1631 + t1491 * t1571 + t1647 * t1691);
t1611 = t1523 * t1555 - t1684;
t1593 = t1611 * t1572 + t1707;
t1665 = t1519 * t1556;
t1681 = (((t1612 * t1554 - t1665) * mrSges(3,1) + t1593 * mrSges(3,2)) * t1571 + t1565 * (t1654 + (t1523 * t1653 + t1665) * mrSges(3,2) + t1593 * mrSges(3,1))) * t1479;
t1574 = cos(qJ(2,2));
t1529 = pkin(2) * t1568 - t1577 * t1574;
t1629 = t1556 * t1567;
t1492 = pkin(3) * t1629 + t1529 * t1554;
t1645 = t1554 * t1568;
t1480 = 0.1e1 / (pkin(2) * t1629 + t1492 * t1573 + t1645 * t1690);
t1609 = t1524 * t1555 - t1684;
t1592 = t1609 * t1574 + t1706;
t1662 = t1520 * t1556;
t1680 = (((t1610 * t1554 - t1662) * mrSges(3,1) + t1592 * mrSges(3,2)) * t1573 + t1567 * (t1654 + (t1524 * t1653 + t1662) * mrSges(3,2) + t1592 * mrSges(3,1))) * t1480;
t1576 = cos(qJ(2,1));
t1530 = pkin(2) * t1570 - t1577 * t1576;
t1627 = t1556 * t1569;
t1493 = pkin(3) * t1627 + t1530 * t1554;
t1643 = t1554 * t1570;
t1481 = 0.1e1 / (pkin(2) * t1627 + t1493 * t1575 + t1643 * t1689);
t1607 = t1525 * t1555 - t1684;
t1591 = t1607 * t1576 + t1705;
t1659 = t1521 * t1556;
t1679 = (((t1608 * t1554 - t1659) * mrSges(3,1) + t1591 * mrSges(3,2)) * t1575 + t1569 * (t1654 + (t1525 * t1653 + t1659) * mrSges(3,2) + t1591 * mrSges(3,1))) * t1481;
t1536 = m(3) * pkin(6) - mrSges(2,2) + mrSges(3,3);
t1534 = t1536 * t1684;
t1638 = t1555 * t1560;
t1462 = t1534 * t1560 + (-t1522 * t1638 - t1708) * t1536 + (t1558 * t1613 - t1560 * t1701) * (t1559 * mrSges(3,1) - mrSges(3,2) * t1557 + t1683);
t1633 = t1556 * t1558;
t1499 = t1553 * t1560 + t1555 * t1633;
t1649 = t1554 * t1559;
t1678 = t1462 * (t1557 * t1499 + t1555 * t1649);
t1637 = t1555 * t1572;
t1463 = t1534 * t1572 + (-t1523 * t1637 - t1707) * t1536 + (t1566 * t1611 - t1572 * t1702) * (t1571 * mrSges(3,1) - mrSges(3,2) * t1565 + t1683);
t1630 = t1556 * t1566;
t1503 = t1553 * t1572 + t1555 * t1630;
t1642 = t1554 * t1571;
t1677 = t1463 * (t1565 * t1503 + t1555 * t1642);
t1636 = t1555 * t1574;
t1464 = t1534 * t1574 + (-t1524 * t1636 - t1706) * t1536 + (t1568 * t1609 - t1574 * t1703) * (t1573 * mrSges(3,1) - mrSges(3,2) * t1567 + t1683);
t1628 = t1556 * t1568;
t1504 = t1553 * t1574 + t1555 * t1628;
t1641 = t1554 * t1573;
t1676 = t1464 * (t1567 * t1504 + t1555 * t1641);
t1635 = t1555 * t1576;
t1465 = t1534 * t1576 + (-t1525 * t1635 - t1705) * t1536 + (t1570 * t1607 - t1576 * t1704) * (t1575 * mrSges(3,1) - mrSges(3,2) * t1569 + t1683);
t1626 = t1556 * t1570;
t1505 = t1553 * t1576 + t1555 * t1626;
t1640 = t1554 * t1575;
t1675 = t1465 * (t1569 * t1505 + t1555 * t1640);
t1674 = t1478 * t1518;
t1673 = t1479 * t1519;
t1672 = t1480 * t1520;
t1671 = t1481 * t1521;
t1549 = m(1) + m(2) + m(3);
t1670 = t1518 * t1549;
t1667 = t1519 * t1549;
t1664 = t1520 * t1549;
t1661 = t1521 * t1549;
t1652 = t1553 * t1556;
t1651 = t1554 * t1557;
t1648 = t1554 * t1565;
t1646 = t1554 * t1567;
t1644 = t1554 * t1569;
t1639 = t1555 * t1556;
t1632 = t1556 * t1560;
t1625 = t1556 * t1572;
t1624 = t1556 * t1574;
t1623 = t1556 * t1576;
t1527 = pkin(2) * t1560 + t1558 * t1577;
t1622 = ((t1553 * t1558 - t1555 * t1632) * t1688 - t1527 * t1639 + t1526 * t1553) * t1682;
t1531 = pkin(2) * t1572 + t1566 * t1577;
t1621 = ((t1553 * t1566 - t1555 * t1625) * t1687 - t1531 * t1639 + t1528 * t1553) * t1681;
t1532 = pkin(2) * t1574 + t1568 * t1577;
t1620 = ((t1553 * t1568 - t1555 * t1624) * t1686 - t1532 * t1639 + t1529 * t1553) * t1680;
t1533 = pkin(2) * t1576 + t1570 * t1577;
t1619 = ((t1553 * t1570 - t1555 * t1623) * t1685 - t1533 * t1639 + t1530 * t1553) * t1679;
t1618 = t1478 * t1678;
t1617 = t1479 * t1677;
t1616 = t1480 * t1676;
t1615 = t1481 * t1675;
t1579 = xP(4);
t1546 = sin(t1579);
t1547 = cos(t1579);
t1582 = koppelP(4,2);
t1586 = koppelP(4,1);
t1506 = -t1546 * t1586 - t1547 * t1582;
t1510 = -t1546 * t1582 + t1547 * t1586;
t1602 = t1506 * t1542 - t1510 * t1538;
t1583 = koppelP(3,2);
t1587 = koppelP(3,1);
t1507 = -t1546 * t1587 - t1547 * t1583;
t1511 = -t1546 * t1583 + t1547 * t1587;
t1601 = t1507 * t1543 - t1511 * t1539;
t1584 = koppelP(2,2);
t1588 = koppelP(2,1);
t1508 = -t1546 * t1588 - t1547 * t1584;
t1512 = -t1546 * t1584 + t1547 * t1588;
t1600 = t1508 * t1544 - t1512 * t1540;
t1585 = koppelP(1,2);
t1589 = koppelP(1,1);
t1509 = -t1546 * t1589 - t1547 * t1585;
t1513 = -t1546 * t1585 + t1547 * t1589;
t1599 = t1509 * t1545 - t1513 * t1541;
t1598 = pkin(3) * t1651 - t1526 * t1556;
t1597 = pkin(3) * t1648 - t1528 * t1556;
t1596 = pkin(3) * t1646 - t1529 * t1556;
t1595 = pkin(3) * t1644 - t1530 * t1556;
t1590 = 0.1e1 / pkin(3);
t1581 = mrSges(4,1);
t1580 = mrSges(4,2);
t1502 = t1553 * t1626 - t1635;
t1501 = t1553 * t1628 - t1636;
t1500 = t1553 * t1630 - t1637;
t1498 = t1553 * t1633 - t1638;
t1485 = t1533 * t1555 + t1595 * t1553;
t1484 = t1532 * t1555 + t1596 * t1553;
t1483 = t1531 * t1555 + t1597 * t1553;
t1482 = t1527 * t1555 + t1598 * t1553;
t1473 = -(t1502 * t1545 - t1541 * t1643) * t1689 + (t1485 * t1545 + t1493 * t1541) * t1575 + (t1556 * t1541 + t1545 * t1653) * t1693;
t1472 = (t1502 * t1541 + t1545 * t1643) * t1689 + (-t1485 * t1541 + t1493 * t1545) * t1575 + (-t1541 * t1653 + t1556 * t1545) * t1693;
t1471 = -(t1501 * t1544 - t1540 * t1645) * t1690 + (t1484 * t1544 + t1492 * t1540) * t1573 + (t1556 * t1540 + t1544 * t1653) * t1694;
t1470 = (t1501 * t1540 + t1544 * t1645) * t1690 + (-t1484 * t1540 + t1492 * t1544) * t1573 + (-t1540 * t1653 + t1556 * t1544) * t1694;
t1469 = -(t1500 * t1543 - t1539 * t1647) * t1691 + (t1483 * t1543 + t1491 * t1539) * t1571 + (t1556 * t1539 + t1543 * t1653) * t1695;
t1468 = (t1500 * t1539 + t1543 * t1647) * t1691 + (-t1483 * t1539 + t1491 * t1543) * t1571 + (-t1539 * t1653 + t1556 * t1543) * t1695;
t1467 = -(t1498 * t1542 - t1538 * t1650) * t1692 + (t1482 * t1542 + t1490 * t1538) * t1559 + (t1556 * t1538 + t1542 * t1653) * t1696;
t1466 = (t1498 * t1538 + t1542 * t1650) * t1692 + (-t1482 * t1538 + t1490 * t1542) * t1559 + (-t1538 * t1653 + t1556 * t1542) * t1696;
t1 = [-t1542 * t1618 - t1543 * t1617 - t1544 * t1616 - t1545 * t1615 - g(1) * m(4) + (t1542 * t1622 + t1543 * t1621 + t1544 * t1620 + t1545 * t1619) * t1590 + (-t1467 * t1674 - t1469 * t1673 - t1471 * t1672 - t1473 * t1671) * t1549; t1538 * t1618 + t1539 * t1617 + t1540 * t1616 + t1541 * t1615 - g(2) * m(4) + (-t1538 * t1622 - t1539 * t1621 - t1540 * t1620 - t1541 * t1619) * t1590 + (-t1466 * t1674 - t1468 * t1673 - t1470 * t1672 - t1472 * t1671) * t1549; -g(3) * m(4) + (-(-t1505 * t1689 - t1533 * t1553 * t1575 + (pkin(2) * t1644 + t1595 * t1575) * t1555) * t1661 + (t1569 * t1502 + t1553 * t1640) * t1465) * t1481 + (-(-t1504 * t1690 - t1532 * t1553 * t1573 + (pkin(2) * t1646 + t1596 * t1573) * t1555) * t1664 + (t1567 * t1501 + t1553 * t1641) * t1464) * t1480 + (-(-t1503 * t1691 - t1531 * t1553 * t1571 + (pkin(2) * t1648 + t1597 * t1571) * t1555) * t1667 + (t1565 * t1500 + t1553 * t1642) * t1463) * t1479 + (-(-t1499 * t1692 - t1527 * t1553 * t1559 + (pkin(2) * t1651 + t1598 * t1559) * t1555) * t1670 + (t1557 * t1498 + t1553 * t1649) * t1462) * t1478 + (((t1553 * t1623 + t1555 * t1570) * t1685 + t1533 * t1652 + t1555 * t1530) * t1679 + ((t1553 * t1624 + t1555 * t1568) * t1686 + t1532 * t1652 + t1555 * t1529) * t1680 + ((t1553 * t1625 + t1555 * t1566) * t1687 + t1531 * t1652 + t1555 * t1528) * t1681 + ((t1553 * t1632 + t1555 * t1558) * t1688 + t1527 * t1652 + t1555 * t1526) * t1682) * t1590; -(-g(1) * t1581 - g(2) * t1580) * t1546 + t1547 * (g(1) * t1580 - g(2) * t1581) + (-(t1472 * t1513 + t1473 * t1509) * t1661 - t1599 * t1675) * t1481 + (-(t1470 * t1512 + t1471 * t1508) * t1664 - t1600 * t1676) * t1480 + (-(t1468 * t1511 + t1469 * t1507) * t1667 - t1601 * t1677) * t1479 + (-(t1466 * t1510 + t1467 * t1506) * t1670 - t1602 * t1678) * t1478 + (t1599 * t1619 + t1600 * t1620 + t1601 * t1621 + t1602 * t1622) * t1590;];
taugX  = t1;
