% Calculate Gravitation load for parallel robot
% P4PRRRR8V1G3A0
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
% Datum: 2020-09-20 23:03
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P4PRRRR8V1G3A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V1G3A0_gravload_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V1G3A0_gravload_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P4PRRRR8V1G3A0_gravload_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V1G3A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR8V1G3A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRRR8V1G3A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V1G3A0_gravload_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V1G3A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-09-20 22:59:14
% EndTime: 2020-09-20 22:59:16
% DurationCPUTime: 2.22s
% Computational Cost: add. (1315->201), mult. (3330->397), div. (120->9), fcn. (2654->30), ass. (0->190)
t1610 = sin(qJ(2,4));
t1612 = cos(qJ(2,4));
t1613 = legFrame(4,2);
t1590 = sin(t1613);
t1594 = cos(t1613);
t1573 = t1594 * g(1) - t1590 * g(2);
t1607 = cos(pkin(6));
t1605 = sin(pkin(6));
t1724 = g(3) * t1605;
t1672 = -t1573 * t1607 + t1724;
t1569 = t1590 * g(1) + t1594 * g(2);
t1606 = sin(pkin(3));
t1608 = cos(pkin(3));
t1589 = g(3) * t1607;
t1673 = -t1573 * t1605 - t1589;
t1742 = t1569 * t1606 + t1673 * t1608;
t1645 = t1742 * t1610 - t1672 * t1612;
t1618 = sin(qJ(2,3));
t1624 = cos(qJ(2,3));
t1614 = legFrame(3,2);
t1591 = sin(t1614);
t1595 = cos(t1614);
t1574 = t1595 * g(1) - t1591 * g(2);
t1670 = -t1574 * t1607 + t1724;
t1570 = t1591 * g(1) + t1595 * g(2);
t1671 = -t1574 * t1605 - t1589;
t1743 = t1570 * t1606 + t1671 * t1608;
t1644 = t1743 * t1618 - t1670 * t1624;
t1620 = sin(qJ(2,2));
t1626 = cos(qJ(2,2));
t1615 = legFrame(2,2);
t1592 = sin(t1615);
t1596 = cos(t1615);
t1575 = t1596 * g(1) - t1592 * g(2);
t1668 = -t1575 * t1607 + t1724;
t1571 = t1592 * g(1) + t1596 * g(2);
t1669 = -t1575 * t1605 - t1589;
t1744 = t1571 * t1606 + t1669 * t1608;
t1643 = t1744 * t1620 - t1668 * t1626;
t1622 = sin(qJ(2,1));
t1628 = cos(qJ(2,1));
t1616 = legFrame(1,2);
t1593 = sin(t1616);
t1597 = cos(t1616);
t1576 = t1597 * g(1) - t1593 * g(2);
t1666 = -t1576 * t1607 + t1724;
t1572 = t1593 * g(1) + t1597 * g(2);
t1667 = -t1576 * t1605 - t1589;
t1745 = t1572 * t1606 + t1667 * t1608;
t1642 = t1745 * t1622 - t1666 * t1628;
t1611 = cos(qJ(3,4));
t1731 = pkin(2) * t1611;
t1567 = -t1612 * pkin(5) + t1610 * t1731;
t1609 = sin(qJ(3,4));
t1732 = pkin(2) * t1609;
t1544 = t1567 * t1606 + t1608 * t1732;
t1737 = 0.1e1 / t1544;
t1623 = cos(qJ(3,3));
t1727 = pkin(2) * t1623;
t1577 = -t1624 * pkin(5) + t1618 * t1727;
t1617 = sin(qJ(3,3));
t1730 = pkin(2) * t1617;
t1548 = t1577 * t1606 + t1608 * t1730;
t1736 = 0.1e1 / t1548;
t1625 = cos(qJ(3,2));
t1726 = pkin(2) * t1625;
t1578 = -t1626 * pkin(5) + t1620 * t1726;
t1619 = sin(qJ(3,2));
t1729 = pkin(2) * t1619;
t1549 = t1578 * t1606 + t1608 * t1729;
t1735 = 0.1e1 / t1549;
t1627 = cos(qJ(3,1));
t1725 = pkin(2) * t1627;
t1579 = -t1628 * pkin(5) + t1622 * t1725;
t1621 = sin(qJ(3,1));
t1728 = pkin(2) * t1621;
t1550 = t1579 * t1606 + t1608 * t1728;
t1734 = 0.1e1 / t1550;
t1733 = m(3) / pkin(2);
t1723 = t1737 / t1611;
t1722 = t1736 / t1623;
t1721 = t1735 / t1625;
t1720 = t1734 / t1627;
t1719 = t1737 * t1569;
t1718 = t1736 * t1570;
t1717 = t1735 * t1571;
t1716 = t1734 * t1572;
t1714 = t1569 * t1608;
t1712 = t1570 * t1608;
t1710 = t1571 * t1608;
t1708 = t1572 * t1608;
t1703 = rSges(3,2) * t1606 * t1589;
t1702 = t1605 * t1606;
t1701 = t1608 * t1610;
t1700 = t1608 * t1612;
t1699 = t1608 * t1618;
t1698 = t1608 * t1620;
t1697 = t1608 * t1622;
t1696 = t1608 * t1624;
t1695 = t1608 * t1626;
t1694 = t1608 * t1628;
t1693 = t1609 * t1612;
t1692 = t1617 * t1624;
t1691 = t1619 * t1626;
t1690 = t1621 * t1628;
t1511 = ((t1673 * t1606 - t1714) * rSges(3,1) + t1645 * rSges(3,2)) * t1611 + (t1703 + (t1573 * t1702 + t1714) * rSges(3,2) + t1645 * rSges(3,1)) * t1609;
t1689 = t1511 * t1723;
t1512 = ((t1671 * t1606 - t1712) * rSges(3,1) + t1644 * rSges(3,2)) * t1623 + (t1703 + (t1574 * t1702 + t1712) * rSges(3,2) + t1644 * rSges(3,1)) * t1617;
t1688 = t1512 * t1722;
t1513 = ((t1669 * t1606 - t1710) * rSges(3,1) + t1643 * rSges(3,2)) * t1625 + (t1703 + (t1575 * t1702 + t1710) * rSges(3,2) + t1643 * rSges(3,1)) * t1619;
t1687 = t1513 * t1721;
t1514 = ((t1667 * t1606 - t1708) * rSges(3,1) + t1642 * rSges(3,2)) * t1627 + (t1703 + (t1576 * t1702 + t1708) * rSges(3,2) + t1642 * rSges(3,1)) * t1621;
t1686 = t1514 * t1720;
t1587 = m(2) * rSges(2,2) - m(3) * rSges(3,3);
t1629 = m(2) * rSges(2,1);
t1515 = t1645 * t1587 + (-t1672 * t1610 - t1612 * t1742) * (t1629 + (rSges(3,1) * t1611 - rSges(3,2) * t1609) * m(3));
t1685 = t1515 * t1723;
t1516 = t1644 * t1587 + (-t1670 * t1618 - t1624 * t1743) * (t1629 + (rSges(3,1) * t1623 - rSges(3,2) * t1617) * m(3));
t1684 = t1516 * t1722;
t1517 = t1643 * t1587 + (-t1668 * t1620 - t1626 * t1744) * (t1629 + (rSges(3,1) * t1625 - rSges(3,2) * t1619) * m(3));
t1683 = t1517 * t1721;
t1518 = t1642 * t1587 + (-t1666 * t1622 - t1628 * t1745) * (t1629 + (rSges(3,1) * t1627 - rSges(3,2) * t1621) * m(3));
t1682 = t1518 * t1720;
t1531 = (-t1605 * t1610 + t1607 * t1700) * t1731 + pkin(5) * (t1605 * t1612 + t1607 * t1701);
t1681 = t1531 * t1689;
t1532 = (-t1605 * t1618 + t1607 * t1696) * t1727 + pkin(5) * (t1605 * t1624 + t1607 * t1699);
t1680 = t1532 * t1688;
t1533 = (-t1605 * t1620 + t1607 * t1695) * t1726 + pkin(5) * (t1605 * t1626 + t1607 * t1698);
t1679 = t1533 * t1687;
t1534 = (-t1605 * t1622 + t1607 * t1694) * t1725 + pkin(5) * (t1605 * t1628 + t1607 * t1697);
t1678 = t1534 * t1686;
t1657 = t1606 * t1611 + t1609 * t1701;
t1535 = t1605 * t1693 + t1657 * t1607;
t1677 = t1535 * t1685;
t1656 = t1606 * t1623 + t1617 * t1699;
t1536 = t1605 * t1692 + t1656 * t1607;
t1676 = t1536 * t1684;
t1655 = t1606 * t1625 + t1619 * t1698;
t1537 = t1605 * t1691 + t1655 * t1607;
t1675 = t1537 * t1683;
t1654 = t1606 * t1627 + t1621 * t1697;
t1538 = t1605 * t1690 + t1654 * t1607;
t1674 = t1538 * t1682;
t1661 = -t1567 * t1608 + t1606 * t1732;
t1660 = -t1577 * t1608 + t1606 * t1730;
t1659 = -t1578 * t1608 + t1606 * t1729;
t1658 = -t1579 * t1608 + t1606 * t1728;
t1630 = xP(4);
t1598 = sin(t1630);
t1599 = cos(t1630);
t1633 = koppelP(4,2);
t1637 = koppelP(4,1);
t1559 = -t1598 * t1637 - t1599 * t1633;
t1563 = -t1598 * t1633 + t1599 * t1637;
t1649 = (-t1559 * t1594 + t1563 * t1590) * t1723;
t1634 = koppelP(3,2);
t1638 = koppelP(3,1);
t1560 = -t1598 * t1638 - t1599 * t1634;
t1564 = -t1598 * t1634 + t1599 * t1638;
t1648 = (-t1560 * t1595 + t1564 * t1591) * t1722;
t1635 = koppelP(2,2);
t1639 = koppelP(2,1);
t1561 = -t1598 * t1639 - t1599 * t1635;
t1565 = -t1598 * t1635 + t1599 * t1639;
t1647 = (-t1561 * t1596 + t1565 * t1592) * t1721;
t1636 = koppelP(1,2);
t1640 = koppelP(1,1);
t1562 = -t1598 * t1640 - t1599 * t1636;
t1566 = -t1598 * t1636 + t1599 * t1640;
t1646 = (-t1562 * t1597 + t1566 * t1593) * t1720;
t1632 = rSges(4,1);
t1631 = rSges(4,2);
t1601 = m(1) + m(2) + m(3);
t1582 = pkin(5) * t1622 + t1628 * t1725;
t1581 = pkin(5) * t1620 + t1626 * t1726;
t1580 = pkin(5) * t1618 + t1624 * t1727;
t1568 = pkin(5) * t1610 + t1612 * t1731;
t1530 = t1607 * t1582 + t1658 * t1605;
t1529 = t1607 * t1581 + t1659 * t1605;
t1528 = t1607 * t1580 + t1660 * t1605;
t1527 = t1607 * t1568 + t1661 * t1605;
t1526 = -t1530 * t1593 + t1550 * t1597;
t1525 = t1530 * t1597 + t1550 * t1593;
t1524 = -t1529 * t1592 + t1549 * t1596;
t1523 = t1529 * t1596 + t1549 * t1592;
t1522 = -t1528 * t1591 + t1548 * t1595;
t1521 = t1528 * t1595 + t1548 * t1591;
t1520 = -t1527 * t1590 + t1544 * t1594;
t1519 = t1527 * t1594 + t1544 * t1590;
t1 = [-t1594 * t1677 - t1595 * t1676 - t1596 * t1675 - t1597 * t1674 - m(4) * g(1) + (-t1519 * t1719 - t1521 * t1718 - t1523 * t1717 - t1525 * t1716) * t1601 + (-t1594 * t1681 - t1595 * t1680 - t1596 * t1679 - t1597 * t1678) * t1733; t1590 * t1677 + t1591 * t1676 + t1592 * t1675 + t1593 * t1674 - m(4) * g(2) + (-t1520 * t1719 - t1522 * t1718 - t1524 * t1717 - t1526 * t1716) * t1601 + (t1590 * t1681 + t1591 * t1680 + t1592 * t1679 + t1593 * t1678) * t1733; (t1654 * t1605 - t1607 * t1690) * t1682 + (t1655 * t1605 - t1607 * t1691) * t1683 + (t1656 * t1605 - t1607 * t1692) * t1684 + (t1657 * t1605 - t1607 * t1693) * t1685 - m(4) * g(3) + (-(-t1605 * t1582 + t1658 * t1607) * t1716 - (-t1605 * t1581 + t1659 * t1607) * t1717 - (-t1605 * t1580 + t1660 * t1607) * t1718 - (-t1605 * t1568 + t1661 * t1607) * t1719) * t1601 + (((t1605 * t1694 + t1607 * t1622) * t1725 + (t1605 * t1697 - t1607 * t1628) * pkin(5)) * t1686 + ((t1605 * t1695 + t1607 * t1620) * t1726 + (t1605 * t1698 - t1607 * t1626) * pkin(5)) * t1687 + ((t1605 * t1696 + t1607 * t1618) * t1727 + (t1605 * t1699 - t1607 * t1624) * pkin(5)) * t1688 + ((t1605 * t1700 + t1607 * t1610) * t1731 + (t1605 * t1701 - t1607 * t1612) * pkin(5)) * t1689) * t1733; m(4) * ((g(1) * t1632 + g(2) * t1631) * t1598 + (g(1) * t1631 - g(2) * t1632) * t1599) + t1518 * t1538 * t1646 + t1517 * t1537 * t1647 + t1516 * t1536 * t1648 + t1515 * t1535 * t1649 + (-(t1525 * t1562 + t1526 * t1566) * t1716 - (t1523 * t1561 + t1524 * t1565) * t1717 - (t1521 * t1560 + t1522 * t1564) * t1718 - (t1519 * t1559 + t1520 * t1563) * t1719) * t1601 + (t1511 * t1531 * t1649 + t1512 * t1532 * t1648 + t1513 * t1533 * t1647 + t1534 * t1514 * t1646) * t1733;];
taugX  = t1;
