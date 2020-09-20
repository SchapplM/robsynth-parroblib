% Calculate Gravitation load for parallel robot
% P4RRRRR2G1A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [4x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 17:26
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P4RRRRR2G1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4RRRRR2G1A0_gravload_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4RRRRR2G1A0_gravload_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4RRRRR2G1A0_gravload_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4RRRRR2G1A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4RRRRR2G1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4RRRRR2G1A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4RRRRR2G1A0_gravload_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4RRRRR2G1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 17:24:24
% EndTime: 2020-08-07 17:24:26
% DurationCPUTime: 1.33s
% Computational Cost: add. (1479->229), mult. (2506->324), div. (116->18), fcn. (1414->98), ass. (0->181)
t1615 = sin(qJ(3,1));
t1618 = cos(qJ(3,1));
t1685 = t1618 * rSges(3,1) - t1615 * rSges(3,2);
t1614 = sin(qJ(3,2));
t1617 = cos(qJ(3,2));
t1684 = t1617 * rSges(3,1) - t1614 * rSges(3,2);
t1613 = sin(qJ(3,3));
t1616 = cos(qJ(3,3));
t1683 = t1616 * rSges(3,1) - t1613 * rSges(3,2);
t1611 = sin(qJ(3,4));
t1612 = cos(qJ(3,4));
t1682 = t1612 * rSges(3,1) - t1611 * rSges(3,2);
t1681 = -2 * pkin(1);
t1680 = -m(3) / 0.2e1;
t1679 = m(3) / 0.2e1;
t1678 = m(1) * rSges(1,2);
t1607 = legFrame(4,3);
t1579 = sin(t1607);
t1583 = cos(t1607);
t1528 = -t1579 * g(1) + t1583 * g(2);
t1676 = rSges(3,1) * t1528;
t1608 = legFrame(3,3);
t1580 = sin(t1608);
t1584 = cos(t1608);
t1529 = -t1580 * g(1) + t1584 * g(2);
t1675 = rSges(3,1) * t1529;
t1609 = legFrame(2,3);
t1581 = sin(t1609);
t1585 = cos(t1609);
t1530 = -t1581 * g(1) + t1585 * g(2);
t1674 = rSges(3,1) * t1530;
t1610 = legFrame(1,3);
t1582 = sin(t1610);
t1586 = cos(t1610);
t1531 = -t1582 * g(1) + t1586 * g(2);
t1673 = rSges(3,1) * t1531;
t1532 = t1583 * g(1) + t1579 * g(2);
t1672 = rSges(3,1) * t1532;
t1533 = t1584 * g(1) + t1580 * g(2);
t1671 = rSges(3,1) * t1533;
t1534 = t1585 * g(1) + t1581 * g(2);
t1670 = rSges(3,1) * t1534;
t1535 = t1586 * g(1) + t1582 * g(2);
t1669 = rSges(3,1) * t1535;
t1668 = rSges(3,3) * t1528;
t1667 = rSges(3,3) * t1529;
t1666 = rSges(3,3) * t1530;
t1665 = rSges(3,3) * t1531;
t1664 = rSges(3,3) * t1532;
t1663 = rSges(3,3) * t1533;
t1662 = rSges(3,3) * t1534;
t1661 = rSges(3,3) * t1535;
t1652 = qJ(2,1) - qJ(3,1);
t1651 = qJ(2,1) + qJ(3,1);
t1650 = qJ(2,2) - qJ(3,2);
t1649 = qJ(2,2) + qJ(3,2);
t1648 = qJ(2,3) - qJ(3,3);
t1647 = qJ(2,3) + qJ(3,3);
t1572 = qJ(1,1) + t1610;
t1571 = qJ(1,2) + t1609;
t1570 = qJ(1,3) + t1608;
t1496 = m(2) * (-rSges(2,1) * t1528 + rSges(2,2) * t1532);
t1500 = m(2) * (rSges(2,1) * t1532 + rSges(2,2) * t1528);
t1512 = rSges(3,2) * t1528;
t1516 = rSges(3,2) * t1532;
t1595 = qJ(1,4) + qJ(2,4);
t1566 = sin(t1595);
t1567 = cos(t1595);
t1637 = qJ(2,4) + qJ(3,4);
t1568 = qJ(1,4) + t1637;
t1638 = qJ(2,4) - qJ(3,4);
t1569 = qJ(1,4) + t1638;
t1636 = m(1) * rSges(1,1) + (m(3) + m(2)) * pkin(1);
t1488 = sin(qJ(1,4)) * (t1528 * t1678 + t1636 * t1532) + (-m(3) * t1664 + t1496) * t1567 + t1566 * (-m(3) * t1668 + t1500) + (-t1636 * t1528 + t1532 * t1678) * cos(qJ(1,4)) + ((t1516 + t1676) * cos(t1569) + (t1512 - t1672) * sin(t1569)) * t1680 + ((t1516 - t1676) * cos(t1568) + (t1512 + t1672) * sin(t1568)) * t1679;
t1596 = 0.1e1 / sin(qJ(2,4));
t1646 = t1488 * t1596;
t1497 = m(2) * (-rSges(2,1) * t1529 + rSges(2,2) * t1533);
t1501 = m(2) * (rSges(2,1) * t1533 + rSges(2,2) * t1529);
t1513 = rSges(3,2) * t1529;
t1517 = rSges(3,2) * t1533;
t1604 = qJ(1,3) + qJ(2,3);
t1573 = sin(t1604);
t1576 = cos(t1604);
t1587 = qJ(1,3) + t1647;
t1588 = qJ(1,3) + t1648;
t1489 = sin(qJ(1,3)) * (t1529 * t1678 + t1636 * t1533) + (-m(3) * t1663 + t1497) * t1576 + t1573 * (-m(3) * t1667 + t1501) + (-t1636 * t1529 + t1533 * t1678) * cos(qJ(1,3)) + ((t1517 + t1675) * cos(t1588) + (t1513 - t1671) * sin(t1588)) * t1680 + ((t1517 - t1675) * cos(t1587) + (t1513 + t1671) * sin(t1587)) * t1679;
t1598 = 0.1e1 / sin(qJ(2,3));
t1645 = t1489 * t1598;
t1498 = m(2) * (-rSges(2,1) * t1530 + rSges(2,2) * t1534);
t1502 = m(2) * (rSges(2,1) * t1534 + rSges(2,2) * t1530);
t1514 = rSges(3,2) * t1530;
t1518 = rSges(3,2) * t1534;
t1605 = qJ(2,2) + qJ(1,2);
t1574 = sin(t1605);
t1577 = cos(t1605);
t1589 = qJ(1,2) + t1649;
t1590 = qJ(1,2) + t1650;
t1490 = sin(qJ(1,2)) * (t1530 * t1678 + t1636 * t1534) + (-m(3) * t1662 + t1498) * t1577 + t1574 * (-m(3) * t1666 + t1502) + (-t1636 * t1530 + t1534 * t1678) * cos(qJ(1,2)) + ((t1518 + t1674) * cos(t1590) + (t1514 - t1670) * sin(t1590)) * t1680 + ((t1518 - t1674) * cos(t1589) + (t1514 + t1670) * sin(t1589)) * t1679;
t1599 = 0.1e1 / sin(qJ(2,2));
t1644 = t1490 * t1599;
t1499 = m(2) * (-rSges(2,1) * t1531 + rSges(2,2) * t1535);
t1503 = m(2) * (rSges(2,1) * t1535 + rSges(2,2) * t1531);
t1515 = rSges(3,2) * t1531;
t1519 = rSges(3,2) * t1535;
t1606 = qJ(1,1) + qJ(2,1);
t1575 = sin(t1606);
t1578 = cos(t1606);
t1591 = qJ(1,1) + t1651;
t1592 = qJ(1,1) + t1652;
t1491 = sin(qJ(1,1)) * (t1531 * t1678 + t1636 * t1535) + (-m(3) * t1661 + t1499) * t1578 + t1575 * (-m(3) * t1665 + t1503) + (-t1636 * t1531 + t1535 * t1678) * cos(qJ(1,1)) + ((t1519 + t1673) * cos(t1592) + (t1515 - t1669) * sin(t1592)) * t1680 + ((t1519 - t1673) * cos(t1591) + (t1515 + t1669) * sin(t1591)) * t1679;
t1600 = 0.1e1 / sin(qJ(2,1));
t1643 = t1491 * t1600;
t1492 = t1496 * t1567 + t1566 * t1500 + ((-t1682 * t1528 - t1664) * t1567 + t1566 * (t1682 * t1532 - t1668)) * m(3);
t1642 = t1492 / (sin(t1637) + sin(t1638));
t1493 = t1497 * t1576 + t1573 * t1501 + ((-t1683 * t1529 - t1663) * t1576 + t1573 * (t1683 * t1533 - t1667)) * m(3);
t1641 = t1493 / (sin(t1647) + sin(t1648));
t1494 = t1498 * t1577 + t1574 * t1502 + ((-t1684 * t1530 - t1662) * t1577 + t1574 * (t1684 * t1534 - t1666)) * m(3);
t1640 = t1494 / (sin(t1649) + sin(t1650));
t1495 = t1499 * t1578 + t1575 * t1503 + ((-t1685 * t1531 - t1661) * t1578 + t1575 * (t1685 * t1535 - t1665)) * m(3);
t1639 = t1495 / (sin(t1651) + sin(t1652));
t1565 = qJ(1,4) + t1607;
t1563 = qJ(2,1) + t1572;
t1562 = qJ(2,2) + t1571;
t1561 = qJ(2,3) + t1570;
t1631 = 1 / pkin(1);
t1635 = t1596 * t1611 * t1631;
t1634 = t1598 * t1613 * t1631;
t1633 = t1599 * t1614 * t1631;
t1632 = t1600 * t1615 * t1631;
t1560 = qJ(2,4) + t1565;
t1630 = 0.1e1 / pkin(2);
t1629 = koppelP(1,1);
t1628 = koppelP(2,1);
t1627 = koppelP(3,1);
t1626 = koppelP(4,1);
t1625 = koppelP(1,2);
t1624 = koppelP(2,2);
t1623 = koppelP(3,2);
t1622 = koppelP(4,2);
t1621 = rSges(4,1);
t1620 = rSges(4,2);
t1619 = xP(4);
t1603 = 0.1e1 / t1618;
t1602 = 0.1e1 / t1617;
t1601 = 0.1e1 / t1616;
t1597 = 0.1e1 / t1612;
t1594 = cos(t1619);
t1593 = sin(t1619);
t1559 = -qJ(3,1) + t1563;
t1558 = qJ(3,1) + t1563;
t1557 = -qJ(3,2) + t1562;
t1556 = qJ(3,2) + t1562;
t1555 = -qJ(3,3) + t1561;
t1554 = qJ(3,3) + t1561;
t1553 = cos(t1563);
t1552 = cos(t1562);
t1551 = cos(t1561);
t1550 = sin(t1563);
t1549 = sin(t1562);
t1548 = sin(t1561);
t1547 = -qJ(3,4) + t1560;
t1546 = qJ(3,4) + t1560;
t1545 = cos(t1560);
t1544 = sin(t1560);
t1527 = -t1593 * t1625 + t1594 * t1629;
t1526 = -t1593 * t1624 + t1594 * t1628;
t1525 = -t1593 * t1623 + t1594 * t1627;
t1524 = -t1593 * t1622 + t1594 * t1626;
t1523 = -t1593 * t1629 - t1594 * t1625;
t1522 = -t1593 * t1628 - t1594 * t1624;
t1521 = -t1593 * t1627 - t1594 * t1623;
t1520 = -t1593 * t1626 - t1594 * t1622;
t1511 = cos(t1572) * t1681 + (-cos(t1558) - cos(t1559)) * pkin(2);
t1510 = cos(t1571) * t1681 + (-cos(t1556) - cos(t1557)) * pkin(2);
t1509 = cos(t1570) * t1681 + (-cos(t1554) - cos(t1555)) * pkin(2);
t1508 = sin(t1572) * t1681 + (-sin(t1559) - sin(t1558)) * pkin(2);
t1507 = sin(t1571) * t1681 + (-sin(t1557) - sin(t1556)) * pkin(2);
t1506 = sin(t1570) * t1681 + (-sin(t1555) - sin(t1554)) * pkin(2);
t1505 = cos(t1565) * t1681 + (-cos(t1546) - cos(t1547)) * pkin(2);
t1504 = sin(t1565) * t1681 + (-sin(t1547) - sin(t1546)) * pkin(2);
t1 = [-m(4) * g(1) + (t1545 * t1646 + t1551 * t1645 + t1552 * t1644 + t1553 * t1643 + (t1505 * t1642 + t1509 * t1641 + t1510 * t1640 + t1511 * t1639) * t1630) * t1631; -m(4) * g(2) + (t1544 * t1646 + t1548 * t1645 + t1549 * t1644 + t1550 * t1643 + (t1504 * t1642 + t1506 * t1641 + t1507 * t1640 + t1508 * t1639) * t1630) * t1631; t1597 * t1488 * t1635 + t1601 * t1489 * t1634 + t1602 * t1490 * t1633 + t1603 * t1491 * t1632 - m(4) * g(3) + ((t1603 * (-g(3) * t1685 + (t1575 * t1531 + t1535 * t1578) * (t1615 * rSges(3,1) + t1618 * rSges(3,2))) + t1602 * (-g(3) * t1684 + (t1574 * t1530 + t1534 * t1577) * (t1614 * rSges(3,1) + t1617 * rSges(3,2))) + t1601 * (-g(3) * t1683 + (t1573 * t1529 + t1533 * t1576) * (t1613 * rSges(3,1) + t1616 * rSges(3,2))) + t1597 * (-g(3) * t1682 + (t1566 * t1528 + t1532 * t1567) * (t1611 * rSges(3,1) + t1612 * rSges(3,2)))) * m(3) - (cos(qJ(2,1)) * pkin(1) + t1618 * pkin(2)) / t1618 ^ 2 * t1495 * t1632 - (cos(qJ(2,2)) * pkin(1) + t1617 * pkin(2)) / t1617 ^ 2 * t1494 * t1633 - (cos(qJ(2,3)) * pkin(1) + t1616 * pkin(2)) / t1616 ^ 2 * t1493 * t1634 - (cos(qJ(2,4)) * pkin(1) + t1612 * pkin(2)) / t1612 ^ 2 * t1492 * t1635) * t1630; m(4) * ((g(1) * t1621 + g(2) * t1620) * t1593 + (g(1) * t1620 - g(2) * t1621) * t1594) + ((t1523 * t1553 + t1527 * t1550) * t1643 + (t1522 * t1552 + t1526 * t1549) * t1644 + (t1521 * t1551 + t1525 * t1548) * t1645 + (t1520 * t1545 + t1524 * t1544) * t1646 + ((t1508 * t1527 + t1511 * t1523) * t1639 + (t1507 * t1526 + t1510 * t1522) * t1640 + (t1506 * t1525 + t1509 * t1521) * t1641 + (t1504 * t1524 + t1505 * t1520) * t1642) * t1630) * t1631;];
taugX  = t1;
