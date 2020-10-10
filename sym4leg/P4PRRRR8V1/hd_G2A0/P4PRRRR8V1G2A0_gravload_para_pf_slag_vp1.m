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
% Datum: 2020-08-07 11:06
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P4PRRRR8V1G2A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V1G2A0_gravload_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V1G2A0_gravload_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P4PRRRR8V1G2A0_gravload_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V1G2A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR8V1G2A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRRR8V1G2A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V1G2A0_gravload_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V1G2A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:02:44
% EndTime: 2020-08-07 11:02:46
% DurationCPUTime: 2.21s
% Computational Cost: add. (1315->209), mult. (3330->406), div. (120->9), fcn. (2654->30), ass. (0->199)
t1607 = sin(qJ(2,4));
t1610 = legFrame(4,2);
t1587 = sin(t1610);
t1591 = cos(t1610);
t1565 = t1587 * g(1) + t1591 * g(2);
t1603 = sin(pkin(3));
t1605 = cos(pkin(3));
t1569 = t1591 * g(1) - t1587 * g(2);
t1602 = sin(pkin(6));
t1586 = g(3) * t1602;
t1604 = cos(pkin(6));
t1665 = t1569 * t1604 - t1586;
t1739 = t1565 * t1603 + t1665 * t1605;
t1746 = t1739 * t1607;
t1615 = sin(qJ(2,3));
t1611 = legFrame(3,2);
t1588 = sin(t1611);
t1592 = cos(t1611);
t1566 = t1588 * g(1) + t1592 * g(2);
t1570 = t1592 * g(1) - t1588 * g(2);
t1663 = t1570 * t1604 - t1586;
t1740 = t1566 * t1603 + t1663 * t1605;
t1745 = t1740 * t1615;
t1617 = sin(qJ(2,2));
t1612 = legFrame(2,2);
t1589 = sin(t1612);
t1593 = cos(t1612);
t1567 = t1589 * g(1) + t1593 * g(2);
t1571 = t1593 * g(1) - t1589 * g(2);
t1661 = t1571 * t1604 - t1586;
t1741 = t1567 * t1603 + t1661 * t1605;
t1744 = t1741 * t1617;
t1619 = sin(qJ(2,1));
t1613 = legFrame(1,2);
t1590 = sin(t1613);
t1594 = cos(t1613);
t1568 = t1590 * g(1) + t1594 * g(2);
t1572 = t1594 * g(1) - t1590 * g(2);
t1659 = t1572 * t1604 - t1586;
t1742 = t1568 * t1603 + t1659 * t1605;
t1743 = t1742 * t1619;
t1609 = cos(qJ(2,4));
t1608 = cos(qJ(3,4));
t1728 = pkin(2) * t1608;
t1563 = -t1609 * pkin(5) + t1607 * t1728;
t1606 = sin(qJ(3,4));
t1729 = pkin(2) * t1606;
t1540 = t1563 * t1603 + t1605 * t1729;
t1734 = 0.1e1 / t1540;
t1621 = cos(qJ(2,3));
t1620 = cos(qJ(3,3));
t1724 = pkin(2) * t1620;
t1573 = -t1621 * pkin(5) + t1615 * t1724;
t1614 = sin(qJ(3,3));
t1727 = pkin(2) * t1614;
t1544 = t1573 * t1603 + t1605 * t1727;
t1733 = 0.1e1 / t1544;
t1623 = cos(qJ(2,2));
t1622 = cos(qJ(3,2));
t1723 = pkin(2) * t1622;
t1574 = -t1623 * pkin(5) + t1617 * t1723;
t1616 = sin(qJ(3,2));
t1726 = pkin(2) * t1616;
t1545 = t1574 * t1603 + t1605 * t1726;
t1732 = 0.1e1 / t1545;
t1625 = cos(qJ(2,1));
t1624 = cos(qJ(3,1));
t1722 = pkin(2) * t1624;
t1575 = -t1625 * pkin(5) + t1619 * t1722;
t1618 = sin(qJ(3,1));
t1725 = pkin(2) * t1618;
t1546 = t1575 * t1603 + t1605 * t1725;
t1731 = 0.1e1 / t1546;
t1730 = m(3) / pkin(2);
t1721 = g(3) * t1604;
t1720 = t1734 / t1608;
t1719 = t1733 / t1620;
t1718 = t1732 / t1622;
t1717 = t1731 / t1624;
t1716 = t1734 * t1565;
t1715 = t1733 * t1566;
t1714 = t1732 * t1567;
t1713 = t1731 * t1568;
t1711 = t1565 * t1605;
t1709 = t1566 * t1605;
t1707 = t1567 * t1605;
t1705 = t1568 * t1605;
t1700 = rSges(3,2) * t1586 * t1603;
t1699 = t1602 * t1609;
t1698 = t1602 * t1621;
t1697 = t1602 * t1623;
t1696 = t1602 * t1625;
t1695 = t1603 * t1604;
t1694 = t1605 * t1607;
t1693 = t1605 * t1609;
t1692 = t1605 * t1615;
t1691 = t1605 * t1617;
t1690 = t1605 * t1619;
t1689 = t1605 * t1621;
t1688 = t1605 * t1623;
t1687 = t1605 * t1625;
t1686 = t1606 * t1609;
t1685 = t1614 * t1621;
t1684 = t1616 * t1623;
t1683 = t1618 * t1625;
t1666 = t1569 * t1602 + t1721;
t1642 = t1666 * t1609 + t1746;
t1507 = ((t1665 * t1603 - t1711) * rSges(3,1) + t1642 * rSges(3,2)) * t1608 + t1606 * (t1700 + (-t1569 * t1695 + t1711) * rSges(3,2) + t1642 * rSges(3,1));
t1682 = t1507 * t1720;
t1664 = t1570 * t1602 + t1721;
t1641 = t1664 * t1621 + t1745;
t1508 = ((t1663 * t1603 - t1709) * rSges(3,1) + t1641 * rSges(3,2)) * t1620 + t1614 * (t1700 + (-t1570 * t1695 + t1709) * rSges(3,2) + t1641 * rSges(3,1));
t1681 = t1508 * t1719;
t1662 = t1571 * t1602 + t1721;
t1640 = t1662 * t1623 + t1744;
t1509 = ((t1661 * t1603 - t1707) * rSges(3,1) + t1640 * rSges(3,2)) * t1622 + t1616 * (t1700 + (-t1571 * t1695 + t1707) * rSges(3,2) + t1640 * rSges(3,1));
t1680 = t1509 * t1718;
t1660 = t1572 * t1602 + t1721;
t1639 = t1660 * t1625 + t1743;
t1510 = ((t1659 * t1603 - t1705) * rSges(3,1) + t1639 * rSges(3,2)) * t1624 + t1618 * (t1700 + (-t1572 * t1695 + t1705) * rSges(3,2) + t1639 * rSges(3,1));
t1679 = t1510 * t1717;
t1584 = m(2) * rSges(2,2) - m(3) * rSges(3,3);
t1579 = t1584 * t1721;
t1626 = m(2) * rSges(2,1);
t1511 = t1579 * t1609 + (t1569 * t1699 + t1746) * t1584 + (t1607 * t1666 - t1609 * t1739) * (t1626 + (rSges(3,1) * t1608 - rSges(3,2) * t1606) * m(3));
t1678 = t1511 * t1720;
t1512 = t1579 * t1621 + (t1570 * t1698 + t1745) * t1584 + (t1615 * t1664 - t1621 * t1740) * (t1626 + (rSges(3,1) * t1620 - rSges(3,2) * t1614) * m(3));
t1677 = t1512 * t1719;
t1513 = t1579 * t1623 + (t1571 * t1697 + t1744) * t1584 + (t1617 * t1662 - t1623 * t1741) * (t1626 + (rSges(3,1) * t1622 - rSges(3,2) * t1616) * m(3));
t1676 = t1513 * t1718;
t1514 = t1579 * t1625 + (t1572 * t1696 + t1743) * t1584 + (t1619 * t1660 - t1625 * t1742) * (t1626 + (rSges(3,1) * t1624 - rSges(3,2) * t1618) * m(3));
t1675 = t1514 * t1717;
t1527 = (t1602 * t1693 + t1604 * t1607) * t1728 + (t1602 * t1694 - t1604 * t1609) * pkin(5);
t1674 = t1527 * t1682;
t1528 = (t1602 * t1689 + t1604 * t1615) * t1724 + (t1602 * t1692 - t1604 * t1621) * pkin(5);
t1673 = t1528 * t1681;
t1529 = (t1602 * t1688 + t1604 * t1617) * t1723 + (t1602 * t1691 - t1604 * t1623) * pkin(5);
t1672 = t1529 * t1680;
t1530 = (t1602 * t1687 + t1604 * t1619) * t1722 + (t1602 * t1690 - t1604 * t1625) * pkin(5);
t1671 = t1530 * t1679;
t1650 = t1603 * t1608 + t1606 * t1694;
t1531 = t1650 * t1602 - t1604 * t1686;
t1670 = t1531 * t1678;
t1649 = t1603 * t1620 + t1614 * t1692;
t1532 = t1649 * t1602 - t1604 * t1685;
t1669 = t1532 * t1677;
t1648 = t1603 * t1622 + t1616 * t1691;
t1533 = t1648 * t1602 - t1604 * t1684;
t1668 = t1533 * t1676;
t1647 = t1603 * t1624 + t1618 * t1690;
t1534 = t1647 * t1602 - t1604 * t1683;
t1667 = t1534 * t1675;
t1654 = -t1563 * t1605 + t1603 * t1729;
t1653 = -t1573 * t1605 + t1603 * t1727;
t1652 = -t1574 * t1605 + t1603 * t1726;
t1651 = -t1575 * t1605 + t1603 * t1725;
t1627 = xP(4);
t1595 = sin(t1627);
t1596 = cos(t1627);
t1630 = koppelP(4,2);
t1634 = koppelP(4,1);
t1555 = -t1595 * t1634 - t1596 * t1630;
t1559 = -t1595 * t1630 + t1596 * t1634;
t1646 = (-t1555 * t1591 + t1559 * t1587) * t1720;
t1631 = koppelP(3,2);
t1635 = koppelP(3,1);
t1556 = -t1595 * t1635 - t1596 * t1631;
t1560 = -t1595 * t1631 + t1596 * t1635;
t1645 = (-t1556 * t1592 + t1560 * t1588) * t1719;
t1632 = koppelP(2,2);
t1636 = koppelP(2,1);
t1557 = -t1595 * t1636 - t1596 * t1632;
t1561 = -t1595 * t1632 + t1596 * t1636;
t1644 = (-t1557 * t1593 + t1561 * t1589) * t1718;
t1633 = koppelP(1,2);
t1637 = koppelP(1,1);
t1558 = -t1595 * t1637 - t1596 * t1633;
t1562 = -t1595 * t1633 + t1596 * t1637;
t1643 = (-t1558 * t1594 + t1562 * t1590) * t1717;
t1629 = rSges(4,1);
t1628 = rSges(4,2);
t1598 = m(1) + m(2) + m(3);
t1578 = pkin(5) * t1619 + t1625 * t1722;
t1577 = pkin(5) * t1617 + t1623 * t1723;
t1576 = pkin(5) * t1615 + t1621 * t1724;
t1564 = pkin(5) * t1607 + t1609 * t1728;
t1526 = -t1602 * t1578 + t1651 * t1604;
t1525 = -t1602 * t1577 + t1652 * t1604;
t1524 = -t1602 * t1576 + t1653 * t1604;
t1523 = -t1602 * t1564 + t1654 * t1604;
t1522 = -t1526 * t1594 + t1590 * t1546;
t1521 = t1526 * t1590 + t1594 * t1546;
t1520 = -t1525 * t1593 + t1589 * t1545;
t1519 = t1525 * t1589 + t1593 * t1545;
t1518 = -t1524 * t1592 + t1588 * t1544;
t1517 = t1524 * t1588 + t1592 * t1544;
t1516 = -t1523 * t1591 + t1587 * t1540;
t1515 = t1523 * t1587 + t1591 * t1540;
t1 = [-t1591 * t1670 - t1592 * t1669 - t1593 * t1668 - t1594 * t1667 - m(4) * g(1) + (-t1516 * t1716 - t1518 * t1715 - t1520 * t1714 - t1522 * t1713) * t1598 + (-t1591 * t1674 - t1592 * t1673 - t1593 * t1672 - t1594 * t1671) * t1730; t1587 * t1670 + t1588 * t1669 + t1589 * t1668 + t1590 * t1667 - m(4) * g(2) + (-t1515 * t1716 - t1517 * t1715 - t1519 * t1714 - t1521 * t1713) * t1598 + (t1587 * t1674 + t1588 * t1673 + t1589 * t1672 + t1590 * t1671) * t1730; (-t1602 * t1683 - t1647 * t1604) * t1675 + (-t1602 * t1684 - t1648 * t1604) * t1676 + (-t1602 * t1685 - t1649 * t1604) * t1677 + (-t1602 * t1686 - t1650 * t1604) * t1678 - m(4) * g(3) + (-(t1604 * t1578 + t1651 * t1602) * t1713 - (t1604 * t1577 + t1652 * t1602) * t1714 - (t1604 * t1576 + t1653 * t1602) * t1715 - (t1604 * t1564 + t1654 * t1602) * t1716) * t1598 + ((-(-t1602 * t1619 + t1604 * t1687) * t1722 - pkin(5) * (t1604 * t1690 + t1696)) * t1679 + (-(-t1602 * t1617 + t1604 * t1688) * t1723 - pkin(5) * (t1604 * t1691 + t1697)) * t1680 + (-(-t1602 * t1615 + t1604 * t1689) * t1724 - pkin(5) * (t1604 * t1692 + t1698)) * t1681 + (-(-t1602 * t1607 + t1604 * t1693) * t1728 - pkin(5) * (t1604 * t1694 + t1699)) * t1682) * t1730; m(4) * ((g(1) * t1629 + g(2) * t1628) * t1595 + (g(1) * t1628 - g(2) * t1629) * t1596) + t1514 * t1534 * t1643 + t1513 * t1533 * t1644 + t1512 * t1532 * t1645 + t1511 * t1531 * t1646 + (-(t1521 * t1562 + t1522 * t1558) * t1713 - (t1519 * t1561 + t1520 * t1557) * t1714 - (t1517 * t1560 + t1518 * t1556) * t1715 - (t1515 * t1559 + t1516 * t1555) * t1716) * t1598 + (t1507 * t1527 * t1646 + t1508 * t1528 * t1645 + t1509 * t1529 * t1644 + t1530 * t1510 * t1643) * t1730;];
taugX  = t1;
