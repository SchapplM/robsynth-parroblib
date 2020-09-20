% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRRRR2G2P2A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x14]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:10
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRRRR2G2P2A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G2P2A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G2P2A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR2G2P2A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G2P2A0_gravload_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G2P2A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G2P2A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:09:25
% EndTime: 2020-03-09 21:09:26
% DurationCPUTime: 1.06s
% Computational Cost: add. (795->162), mult. (1194->305), div. (294->14), fcn. (1368->42), ass. (0->169)
t1659 = cos(qJ(3,1));
t1754 = t1659 ^ 2;
t1656 = cos(qJ(3,2));
t1753 = t1656 ^ 2;
t1653 = cos(qJ(3,3));
t1752 = t1653 ^ 2;
t1643 = legFrame(1,2);
t1622 = sin(t1643);
t1625 = cos(t1643);
t1603 = -t1625 * g(1) + g(2) * t1622;
t1640 = qJ(1,1) + qJ(2,1);
t1616 = sin(t1640);
t1619 = cos(t1640);
t1670 = g(3) * t1616 + t1603 * t1619;
t1651 = sin(qJ(2,1));
t1628 = 0.1e1 / t1651;
t1636 = 0.1e1 / t1659;
t1704 = t1628 * t1636;
t1751 = t1670 * t1704;
t1642 = legFrame(2,2);
t1621 = sin(t1642);
t1624 = cos(t1642);
t1601 = -t1624 * g(1) + g(2) * t1621;
t1639 = qJ(1,2) + qJ(2,2);
t1615 = sin(t1639);
t1618 = cos(t1639);
t1671 = g(3) * t1615 + t1601 * t1618;
t1648 = sin(qJ(2,2));
t1627 = 0.1e1 / t1648;
t1633 = 0.1e1 / t1656;
t1706 = t1627 * t1633;
t1750 = t1671 * t1706;
t1641 = legFrame(3,2);
t1620 = sin(t1641);
t1623 = cos(t1641);
t1599 = -t1623 * g(1) + g(2) * t1620;
t1638 = qJ(1,3) + qJ(2,3);
t1614 = sin(t1638);
t1617 = cos(t1638);
t1672 = g(3) * t1614 + t1599 * t1617;
t1645 = sin(qJ(2,3));
t1626 = 0.1e1 / t1645;
t1630 = 0.1e1 / t1653;
t1708 = t1626 * t1630;
t1749 = t1672 * t1708;
t1748 = -2 * pkin(1);
t1646 = sin(qJ(1,3));
t1747 = pkin(1) * t1646;
t1649 = sin(qJ(1,2));
t1746 = pkin(1) * t1649;
t1652 = sin(qJ(1,1));
t1745 = pkin(1) * t1652;
t1744 = qJ(2,1) - qJ(3,1);
t1743 = qJ(2,1) + qJ(3,1);
t1742 = qJ(2,2) - qJ(3,2);
t1741 = qJ(2,2) + qJ(3,2);
t1740 = qJ(2,3) - qJ(3,3);
t1739 = qJ(2,3) + qJ(3,3);
t1738 = t1672 * t1626;
t1737 = t1672 * t1653;
t1736 = t1671 * t1627;
t1735 = t1671 * t1656;
t1734 = t1670 * t1628;
t1733 = t1670 * t1659;
t1644 = sin(qJ(3,3));
t1732 = t1672 * t1644;
t1647 = sin(qJ(3,2));
t1731 = t1671 * t1647;
t1650 = sin(qJ(3,1));
t1730 = t1670 * t1650;
t1655 = cos(qJ(1,3));
t1729 = (t1655 * t1748 + (-cos(qJ(1,3) + t1740) - cos(qJ(1,3) + t1739)) * pkin(2)) / (sin(t1739) + sin(t1740));
t1658 = cos(qJ(1,2));
t1728 = (t1658 * t1748 + (-cos(qJ(1,2) + t1742) - cos(qJ(1,2) + t1741)) * pkin(2)) / (sin(t1741) + sin(t1742));
t1661 = cos(qJ(1,1));
t1727 = (t1661 * t1748 + (-cos(qJ(1,1) + t1744) - cos(qJ(1,1) + t1743)) * pkin(2)) / (sin(t1743) + sin(t1744));
t1654 = cos(qJ(2,3));
t1596 = t1645 * t1655 + t1646 * t1654;
t1726 = t1596 * t1653;
t1657 = cos(qJ(2,2));
t1597 = t1648 * t1658 + t1649 * t1657;
t1725 = t1597 * t1656;
t1660 = cos(qJ(2,1));
t1598 = t1651 * t1661 + t1652 * t1660;
t1724 = t1598 * t1659;
t1723 = t1617 * t1626;
t1722 = t1618 * t1627;
t1721 = t1619 * t1628;
t1720 = t1620 * t1630;
t1719 = t1620 * t1644;
t1718 = t1621 * t1633;
t1717 = t1621 * t1647;
t1716 = t1622 * t1636;
t1715 = t1622 * t1650;
t1714 = t1623 * t1630;
t1713 = t1623 * t1644;
t1712 = t1624 * t1633;
t1711 = t1624 * t1647;
t1710 = t1625 * t1636;
t1709 = t1625 * t1650;
t1631 = 0.1e1 / t1752;
t1707 = t1626 * t1631;
t1634 = 0.1e1 / t1753;
t1705 = t1627 * t1634;
t1637 = 0.1e1 / t1754;
t1703 = t1628 * t1637;
t1702 = pkin(1) * t1644 * t1654;
t1701 = pkin(1) * t1647 * t1657;
t1700 = pkin(1) * t1650 * t1660;
t1699 = pkin(2) * t1596 * t1752;
t1698 = pkin(2) * t1597 * t1753;
t1697 = pkin(2) * t1598 * t1754;
t1693 = t1626 * t1732;
t1692 = t1627 * t1731;
t1691 = t1628 * t1730;
t1689 = t1672 * t1707;
t1687 = t1671 * t1705;
t1685 = t1670 * t1703;
t1584 = g(3) * t1617 - t1599 * t1614;
t1684 = t1584 * t1708;
t1683 = t1584 * t1707;
t1585 = g(3) * t1618 - t1601 * t1615;
t1682 = t1585 * t1706;
t1681 = t1585 * t1705;
t1586 = g(3) * t1619 - t1603 * t1616;
t1680 = t1586 * t1704;
t1679 = t1586 * t1703;
t1587 = g(3) * t1646 + t1599 * t1655;
t1678 = t1587 * t1708;
t1588 = g(3) * t1649 + t1601 * t1658;
t1677 = t1588 * t1706;
t1589 = g(3) * t1652 + t1603 * t1661;
t1676 = t1589 * t1704;
t1590 = g(3) * t1655 - t1599 * t1646;
t1675 = t1590 * t1708;
t1591 = g(3) * t1658 - t1601 * t1649;
t1674 = t1591 * t1706;
t1592 = g(3) * t1661 - t1603 * t1652;
t1673 = t1592 * t1704;
t1669 = t1630 * t1693;
t1668 = t1631 * t1693;
t1667 = t1633 * t1692;
t1666 = t1634 * t1692;
t1665 = t1636 * t1691;
t1664 = t1637 * t1691;
t1663 = 1 / pkin(1);
t1662 = 0.1e1 / pkin(2);
t1604 = g(1) * t1622 + g(2) * t1625;
t1602 = g(1) * t1621 + g(2) * t1624;
t1600 = g(1) * t1620 + g(2) * t1623;
t1574 = t1625 * t1724 + t1715;
t1573 = -t1622 * t1724 + t1709;
t1572 = t1624 * t1725 + t1717;
t1571 = -t1621 * t1725 + t1711;
t1570 = t1623 * t1726 + t1719;
t1569 = -t1620 * t1726 + t1713;
t1568 = t1586 * t1659 + t1604 * t1650;
t1567 = t1586 * t1650 - t1604 * t1659;
t1566 = t1585 * t1656 + t1602 * t1647;
t1565 = t1585 * t1647 - t1602 * t1656;
t1564 = t1584 * t1653 + t1600 * t1644;
t1563 = t1584 * t1644 - t1600 * t1653;
t1562 = -t1625 * t1697 + (-pkin(2) * t1715 - t1625 * t1745) * t1659 - t1622 * t1700;
t1561 = t1622 * t1697 + (-pkin(2) * t1709 + t1622 * t1745) * t1659 - t1625 * t1700;
t1560 = -t1624 * t1698 + (-pkin(2) * t1717 - t1624 * t1746) * t1656 - t1621 * t1701;
t1559 = t1621 * t1698 + (-pkin(2) * t1711 + t1621 * t1746) * t1656 - t1624 * t1701;
t1558 = -t1623 * t1699 + (-pkin(2) * t1719 - t1623 * t1747) * t1653 - t1620 * t1702;
t1557 = t1620 * t1699 + (-pkin(2) * t1713 + t1620 * t1747) * t1653 - t1623 * t1702;
t1 = [0, (t1570 * t1678 + t1572 * t1677 + t1574 * t1676) * t1663, (t1570 * t1675 + t1572 * t1674 + t1574 * t1673) * t1663, 0, (t1570 * t1749 + t1572 * t1750 + t1574 * t1751 + (t1558 * t1689 + t1560 * t1687 + t1562 * t1685) * t1662) * t1663, (t1570 * t1684 + t1572 * t1682 + t1574 * t1680 + (t1558 * t1683 + t1560 * t1681 + t1562 * t1679) * t1662) * t1663, 0, 0, 0, 0, 0, (t1570 * t1738 + t1572 * t1736 + t1574 * t1734) * t1663 + (t1563 * t1720 + t1565 * t1718 + t1567 * t1716 + (t1558 * t1749 + t1560 * t1750 + t1562 * t1751) * t1663) * t1662, (-t1570 * t1669 - t1572 * t1667 - t1574 * t1665) * t1663 + (t1564 * t1720 + t1566 * t1718 + t1568 * t1716 + (-t1558 * t1668 - t1560 * t1666 - t1562 * t1664) * t1663) * t1662, -g(1); 0, (t1569 * t1678 + t1571 * t1677 + t1573 * t1676) * t1663, (t1569 * t1675 + t1571 * t1674 + t1573 * t1673) * t1663, 0, (t1569 * t1749 + t1571 * t1750 + t1573 * t1751 + (t1557 * t1689 + t1559 * t1687 + t1561 * t1685) * t1662) * t1663, (t1569 * t1684 + t1571 * t1682 + t1573 * t1680 + (t1557 * t1683 + t1559 * t1681 + t1561 * t1679) * t1662) * t1663, 0, 0, 0, 0, 0, (t1569 * t1738 + t1571 * t1736 + t1573 * t1734) * t1663 + (t1563 * t1714 + t1565 * t1712 + t1567 * t1710 + (t1557 * t1749 + t1559 * t1750 + t1561 * t1751) * t1663) * t1662, (-t1569 * t1669 - t1571 * t1667 - t1573 * t1665) * t1663 + (t1564 * t1714 + t1566 * t1712 + t1568 * t1710 + (-t1557 * t1668 - t1559 * t1666 - t1561 * t1664) * t1663) * t1662, -g(2); 0, (t1587 * t1723 + t1588 * t1722 + t1589 * t1721) * t1663, (t1590 * t1723 + t1591 * t1722 + t1592 * t1721) * t1663, 0, (t1672 * t1723 + t1671 * t1722 + t1670 * t1721 + (t1670 * t1727 + t1671 * t1728 + t1672 * t1729) * t1662) * t1663, (t1584 * t1723 + t1585 * t1722 + t1586 * t1721 + (t1584 * t1729 + t1585 * t1728 + t1586 * t1727) * t1662) * t1663, 0, 0, 0, 0, 0, (t1723 * t1737 + t1722 * t1735 + t1721 * t1733 + (t1727 * t1733 + t1728 * t1735 + t1729 * t1737) * t1662) * t1663, (-t1617 * t1693 - t1618 * t1692 - t1619 * t1691 + (-t1727 * t1730 - t1728 * t1731 - t1729 * t1732) * t1662) * t1663, -g(3);];
tau_reg  = t1;
