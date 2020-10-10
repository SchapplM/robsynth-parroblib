% Calculate Gravitation load for parallel robot
% P4PRRRR8V2G2A0
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
% Datum: 2020-08-07 11:22
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P4PRRRR8V2G2A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G2A0_gravload_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G2A0_gravload_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G2A0_gravload_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V2G2A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR8V2G2A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRRR8V2G2A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G2A0_gravload_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G2A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:16:33
% EndTime: 2020-08-07 11:16:36
% DurationCPUTime: 2.83s
% Computational Cost: add. (1851->264), mult. (3978->504), div. (80->5), fcn. (3150->30), ass. (0->211)
t1613 = sin(qJ(2,4));
t1616 = legFrame(4,2);
t1593 = sin(t1616);
t1597 = cos(t1616);
t1573 = t1593 * g(1) + t1597 * g(2);
t1609 = sin(pkin(4));
t1611 = cos(pkin(4));
t1577 = t1597 * g(1) - t1593 * g(2);
t1608 = sin(pkin(8));
t1592 = g(3) * t1608;
t1610 = cos(pkin(8));
t1669 = t1577 * t1610 - t1592;
t1759 = t1573 * t1609 + t1669 * t1611;
t1766 = t1759 * t1613;
t1621 = sin(qJ(2,3));
t1617 = legFrame(3,2);
t1594 = sin(t1617);
t1598 = cos(t1617);
t1574 = t1594 * g(1) + t1598 * g(2);
t1578 = t1598 * g(1) - t1594 * g(2);
t1667 = t1578 * t1610 - t1592;
t1760 = t1574 * t1609 + t1667 * t1611;
t1765 = t1760 * t1621;
t1623 = sin(qJ(2,2));
t1618 = legFrame(2,2);
t1595 = sin(t1618);
t1599 = cos(t1618);
t1575 = t1595 * g(1) + t1599 * g(2);
t1579 = t1599 * g(1) - t1595 * g(2);
t1665 = t1579 * t1610 - t1592;
t1761 = t1575 * t1609 + t1665 * t1611;
t1764 = t1761 * t1623;
t1625 = sin(qJ(2,1));
t1619 = legFrame(1,2);
t1596 = sin(t1619);
t1600 = cos(t1619);
t1576 = t1596 * g(1) + t1600 * g(2);
t1580 = t1600 * g(1) - t1596 * g(2);
t1663 = t1580 * t1610 - t1592;
t1762 = t1576 * t1609 + t1663 * t1611;
t1763 = t1762 * t1625;
t1754 = m(3) / pkin(3);
t1612 = sin(qJ(3,4));
t1753 = pkin(2) * t1612;
t1620 = sin(qJ(3,3));
t1752 = pkin(2) * t1620;
t1622 = sin(qJ(3,2));
t1751 = pkin(2) * t1622;
t1624 = sin(qJ(3,1));
t1750 = pkin(2) * t1624;
t1614 = cos(qJ(3,4));
t1749 = pkin(3) * t1614 ^ 2;
t1626 = cos(qJ(3,3));
t1748 = pkin(3) * t1626 ^ 2;
t1628 = cos(qJ(3,2));
t1747 = pkin(3) * t1628 ^ 2;
t1630 = cos(qJ(3,1));
t1746 = pkin(3) * t1630 ^ 2;
t1745 = pkin(3) * t1614;
t1744 = pkin(3) * t1626;
t1743 = pkin(3) * t1628;
t1742 = pkin(3) * t1630;
t1741 = g(3) * t1610;
t1615 = cos(qJ(2,4));
t1633 = pkin(7) + pkin(6);
t1582 = pkin(2) * t1613 - t1633 * t1615;
t1691 = t1611 * t1612;
t1545 = pkin(3) * t1691 + t1582 * t1609;
t1703 = t1609 * t1613;
t1533 = 0.1e1 / (pkin(2) * t1691 + t1545 * t1614 + t1703 * t1749);
t1670 = t1577 * t1608 + t1741;
t1650 = t1670 * t1615 + t1766;
t1705 = t1609 * t1610;
t1712 = rSges(3,2) * t1592 * t1609;
t1726 = t1573 * t1611;
t1740 = (((t1669 * t1609 - t1726) * rSges(3,1) + t1650 * rSges(3,2)) * t1614 + t1612 * (t1712 + (-t1577 * t1705 + t1726) * rSges(3,2) + t1650 * rSges(3,1))) * t1533;
t1627 = cos(qJ(2,3));
t1584 = pkin(2) * t1621 - t1633 * t1627;
t1688 = t1611 * t1620;
t1546 = pkin(3) * t1688 + t1584 * t1609;
t1701 = t1609 * t1621;
t1534 = 0.1e1 / (pkin(2) * t1688 + t1546 * t1626 + t1701 * t1748);
t1668 = t1578 * t1608 + t1741;
t1649 = t1668 * t1627 + t1765;
t1723 = t1574 * t1611;
t1739 = (((t1667 * t1609 - t1723) * rSges(3,1) + t1649 * rSges(3,2)) * t1626 + t1620 * (t1712 + (-t1578 * t1705 + t1723) * rSges(3,2) + t1649 * rSges(3,1))) * t1534;
t1629 = cos(qJ(2,2));
t1585 = pkin(2) * t1623 - t1633 * t1629;
t1686 = t1611 * t1622;
t1547 = pkin(3) * t1686 + t1585 * t1609;
t1699 = t1609 * t1623;
t1535 = 0.1e1 / (pkin(2) * t1686 + t1547 * t1628 + t1699 * t1747);
t1666 = t1579 * t1608 + t1741;
t1648 = t1666 * t1629 + t1764;
t1720 = t1575 * t1611;
t1738 = (((t1665 * t1609 - t1720) * rSges(3,1) + t1648 * rSges(3,2)) * t1628 + t1622 * (t1712 + (-t1579 * t1705 + t1720) * rSges(3,2) + t1648 * rSges(3,1))) * t1535;
t1631 = cos(qJ(2,1));
t1586 = pkin(2) * t1625 - t1633 * t1631;
t1684 = t1611 * t1624;
t1548 = pkin(3) * t1684 + t1586 * t1609;
t1697 = t1609 * t1625;
t1536 = 0.1e1 / (pkin(2) * t1684 + t1548 * t1630 + t1697 * t1746);
t1664 = t1580 * t1608 + t1741;
t1647 = t1664 * t1631 + t1763;
t1717 = t1576 * t1611;
t1737 = (((t1663 * t1609 - t1717) * rSges(3,1) + t1647 * rSges(3,2)) * t1630 + t1624 * (t1712 + (-t1580 * t1705 + t1717) * rSges(3,2) + t1647 * rSges(3,1))) * t1536;
t1590 = (-pkin(6) - rSges(3,3)) * m(3) + m(2) * rSges(2,2);
t1581 = t1590 * t1741;
t1679 = m(2) * rSges(2,1) + pkin(2) * m(3);
t1709 = t1608 * t1615;
t1517 = t1581 * t1615 + (t1577 * t1709 + t1766) * t1590 + (t1613 * t1670 - t1615 * t1759) * ((rSges(3,1) * t1614 - rSges(3,2) * t1612) * m(3) + t1679);
t1690 = t1611 * t1613;
t1554 = t1608 * t1690 - t1610 * t1615;
t1711 = t1608 * t1609;
t1736 = t1517 * (t1612 * t1554 + t1614 * t1711);
t1708 = t1608 * t1627;
t1518 = t1581 * t1627 + (t1578 * t1708 + t1765) * t1590 + (t1621 * t1668 - t1627 * t1760) * ((rSges(3,1) * t1626 - rSges(3,2) * t1620) * m(3) + t1679);
t1687 = t1611 * t1621;
t1559 = t1608 * t1687 - t1610 * t1627;
t1735 = t1518 * (t1620 * t1559 + t1626 * t1711);
t1707 = t1608 * t1629;
t1519 = t1581 * t1629 + (t1579 * t1707 + t1764) * t1590 + (t1623 * t1666 - t1629 * t1761) * ((rSges(3,1) * t1628 - rSges(3,2) * t1622) * m(3) + t1679);
t1685 = t1611 * t1623;
t1560 = t1608 * t1685 - t1610 * t1629;
t1734 = t1519 * (t1622 * t1560 + t1628 * t1711);
t1706 = t1608 * t1631;
t1520 = t1581 * t1631 + (t1580 * t1706 + t1763) * t1590 + (t1625 * t1664 - t1631 * t1762) * ((rSges(3,1) * t1630 - rSges(3,2) * t1624) * m(3) + t1679);
t1683 = t1611 * t1625;
t1561 = t1608 * t1683 - t1610 * t1631;
t1733 = t1520 * (t1624 * t1561 + t1630 * t1711);
t1732 = t1533 * t1573;
t1731 = t1534 * t1574;
t1730 = t1535 * t1575;
t1729 = t1536 * t1576;
t1604 = m(1) + m(2) + m(3);
t1728 = t1573 * t1604;
t1725 = t1574 * t1604;
t1722 = t1575 * t1604;
t1719 = t1576 * t1604;
t1710 = t1608 * t1611;
t1704 = t1609 * t1612;
t1702 = t1609 * t1620;
t1700 = t1609 * t1622;
t1698 = t1609 * t1624;
t1696 = t1610 * t1611;
t1695 = t1610 * t1614;
t1694 = t1610 * t1626;
t1693 = t1610 * t1628;
t1692 = t1610 * t1630;
t1689 = t1611 * t1615;
t1682 = t1611 * t1627;
t1681 = t1611 * t1629;
t1680 = t1611 * t1631;
t1583 = pkin(2) * t1615 + t1613 * t1633;
t1678 = ((t1608 * t1689 + t1610 * t1613) * t1745 + t1583 * t1710 + t1582 * t1610) * t1740;
t1587 = pkin(2) * t1627 + t1621 * t1633;
t1677 = ((t1608 * t1682 + t1610 * t1621) * t1744 + t1587 * t1710 + t1584 * t1610) * t1739;
t1588 = pkin(2) * t1629 + t1623 * t1633;
t1676 = ((t1608 * t1681 + t1610 * t1623) * t1743 + t1588 * t1710 + t1585 * t1610) * t1738;
t1589 = pkin(2) * t1631 + t1625 * t1633;
t1675 = ((t1608 * t1680 + t1610 * t1625) * t1742 + t1589 * t1710 + t1586 * t1610) * t1737;
t1674 = t1533 * t1736;
t1673 = t1534 * t1735;
t1672 = t1535 * t1734;
t1671 = t1536 * t1733;
t1635 = xP(4);
t1601 = sin(t1635);
t1602 = cos(t1635);
t1638 = koppelP(4,2);
t1642 = koppelP(4,1);
t1565 = -t1601 * t1642 - t1602 * t1638;
t1569 = -t1601 * t1638 + t1602 * t1642;
t1658 = -t1565 * t1597 + t1569 * t1593;
t1639 = koppelP(3,2);
t1643 = koppelP(3,1);
t1566 = -t1601 * t1643 - t1602 * t1639;
t1570 = -t1601 * t1639 + t1602 * t1643;
t1657 = -t1566 * t1598 + t1570 * t1594;
t1640 = koppelP(2,2);
t1644 = koppelP(2,1);
t1567 = -t1601 * t1644 - t1602 * t1640;
t1571 = -t1601 * t1640 + t1602 * t1644;
t1656 = -t1567 * t1599 + t1571 * t1595;
t1641 = koppelP(1,2);
t1645 = koppelP(1,1);
t1568 = -t1601 * t1645 - t1602 * t1641;
t1572 = -t1601 * t1641 + t1602 * t1645;
t1655 = -t1568 * t1600 + t1572 * t1596;
t1654 = pkin(3) * t1704 - t1582 * t1611;
t1653 = pkin(3) * t1702 - t1584 * t1611;
t1652 = pkin(3) * t1700 - t1585 * t1611;
t1651 = pkin(3) * t1698 - t1586 * t1611;
t1637 = rSges(4,1);
t1636 = rSges(4,2);
t1564 = t1610 * t1683 + t1706;
t1563 = t1610 * t1685 + t1707;
t1562 = t1610 * t1687 + t1708;
t1555 = t1610 * t1690 + t1709;
t1540 = -t1589 * t1608 + t1651 * t1610;
t1539 = -t1588 * t1608 + t1652 * t1610;
t1538 = -t1587 * t1608 + t1653 * t1610;
t1537 = -t1583 * t1608 + t1654 * t1610;
t1528 = -(t1564 * t1596 - t1600 * t1697) * t1746 + (t1540 * t1596 + t1600 * t1548) * t1630 + (t1596 * t1705 + t1611 * t1600) * t1750;
t1527 = -(t1563 * t1595 - t1599 * t1699) * t1747 + (t1539 * t1595 + t1599 * t1547) * t1628 + (t1595 * t1705 + t1611 * t1599) * t1751;
t1526 = -(t1562 * t1594 - t1598 * t1701) * t1748 + (t1538 * t1594 + t1598 * t1546) * t1626 + (t1594 * t1705 + t1611 * t1598) * t1752;
t1525 = (t1564 * t1600 + t1596 * t1697) * t1746 + (-t1540 * t1600 + t1596 * t1548) * t1630 + (t1611 * t1596 - t1600 * t1705) * t1750;
t1524 = (t1563 * t1599 + t1595 * t1699) * t1747 + (-t1539 * t1599 + t1595 * t1547) * t1628 + (t1611 * t1595 - t1599 * t1705) * t1751;
t1523 = (t1562 * t1598 + t1594 * t1701) * t1748 + (-t1538 * t1598 + t1594 * t1546) * t1626 + (t1611 * t1594 - t1598 * t1705) * t1752;
t1522 = -(t1555 * t1593 - t1597 * t1703) * t1749 + (t1537 * t1593 + t1597 * t1545) * t1614 + (t1593 * t1705 + t1611 * t1597) * t1753;
t1521 = (t1555 * t1597 + t1593 * t1703) * t1749 + (-t1537 * t1597 + t1593 * t1545) * t1614 + (t1611 * t1593 - t1597 * t1705) * t1753;
t1 = [-t1597 * t1674 - t1598 * t1673 - t1599 * t1672 - t1600 * t1671 - m(4) * g(1) + (-t1521 * t1732 - t1523 * t1731 - t1524 * t1730 - t1525 * t1729) * t1604 + (-t1597 * t1678 - t1598 * t1677 - t1599 * t1676 - t1600 * t1675) * t1754; t1593 * t1674 + t1594 * t1673 + t1595 * t1672 + t1596 * t1671 - m(4) * g(2) + (-t1522 * t1732 - t1526 * t1731 - t1527 * t1730 - t1528 * t1729) * t1604 + (t1593 * t1678 + t1594 * t1677 + t1595 * t1676 + t1596 * t1675) * t1754; -m(4) * g(3) + (-(-t1561 * t1746 + t1589 * t1692 + (pkin(2) * t1698 + t1651 * t1630) * t1608) * t1719 + (-t1624 * t1564 - t1609 * t1692) * t1520) * t1536 + (-(-t1560 * t1747 + t1588 * t1693 + (pkin(2) * t1700 + t1652 * t1628) * t1608) * t1722 + (-t1622 * t1563 - t1609 * t1693) * t1519) * t1535 + (-(-t1559 * t1748 + t1587 * t1694 + (pkin(2) * t1702 + t1653 * t1626) * t1608) * t1725 + (-t1620 * t1562 - t1609 * t1694) * t1518) * t1534 + (-(-t1554 * t1749 + t1583 * t1695 + (pkin(2) * t1704 + t1654 * t1614) * t1608) * t1728 + (-t1612 * t1555 - t1609 * t1695) * t1517) * t1533 + (((t1608 * t1625 - t1610 * t1680) * t1742 - t1589 * t1696 + t1586 * t1608) * t1737 + ((t1608 * t1623 - t1610 * t1681) * t1743 - t1588 * t1696 + t1585 * t1608) * t1738 + ((t1608 * t1621 - t1610 * t1682) * t1744 - t1587 * t1696 + t1584 * t1608) * t1739 + ((t1608 * t1613 - t1610 * t1689) * t1745 - t1583 * t1696 + t1582 * t1608) * t1740) * t1754; ((g(1) * t1636 - g(2) * t1637) * t1602 + t1601 * (g(1) * t1637 + g(2) * t1636)) * m(4) + (-(t1525 * t1568 + t1528 * t1572) * t1719 + t1655 * t1733) * t1536 + (-(t1524 * t1567 + t1527 * t1571) * t1722 + t1656 * t1734) * t1535 + (-(t1523 * t1566 + t1526 * t1570) * t1725 + t1657 * t1735) * t1534 + (-(t1521 * t1565 + t1522 * t1569) * t1728 + t1658 * t1736) * t1533 + (t1655 * t1675 + t1656 * t1676 + t1657 * t1677 + t1658 * t1678) * t1754;];
taugX  = t1;
