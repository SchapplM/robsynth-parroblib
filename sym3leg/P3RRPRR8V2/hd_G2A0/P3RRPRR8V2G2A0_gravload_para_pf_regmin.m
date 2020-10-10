% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRPRR8V2G2A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x13]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 21:13
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRPRR8V2G2A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:12:54
% EndTime: 2020-08-06 21:12:55
% DurationCPUTime: 0.98s
% Computational Cost: add. (912->181), mult. (1587->312), div. (123->9), fcn. (1374->23), ass. (0->166)
t1709 = sin(qJ(1,3));
t1715 = cos(qJ(1,3));
t1705 = legFrame(3,2);
t1680 = sin(t1705);
t1683 = cos(t1705);
t1733 = t1683 * g(1) - t1680 * g(2);
t1644 = g(3) * t1709 - t1733 * t1715;
t1711 = sin(qJ(1,2));
t1717 = cos(qJ(1,2));
t1706 = legFrame(2,2);
t1681 = sin(t1706);
t1684 = cos(t1706);
t1732 = t1684 * g(1) - t1681 * g(2);
t1645 = g(3) * t1711 - t1732 * t1717;
t1713 = sin(qJ(1,1));
t1719 = cos(qJ(1,1));
t1707 = legFrame(1,2);
t1682 = sin(t1707);
t1685 = cos(t1707);
t1731 = t1685 * g(1) - t1682 * g(2);
t1646 = g(3) * t1713 - t1731 * t1719;
t1714 = cos(qJ(2,3));
t1806 = 0.2e1 * t1714 ^ 2;
t1716 = cos(qJ(2,2));
t1805 = 0.2e1 * t1716 ^ 2;
t1718 = cos(qJ(2,1));
t1804 = 0.2e1 * t1718 ^ 2;
t1704 = -qJ(3,1) - pkin(5);
t1703 = -qJ(3,2) - pkin(5);
t1702 = -qJ(3,3) - pkin(5);
t1693 = pkin(6) - t1702;
t1742 = pkin(1) * t1709 - t1715 * t1693;
t1701 = cos(pkin(7));
t1696 = t1701 ^ 2;
t1764 = pkin(3) * (t1696 - 0.1e1);
t1700 = sin(pkin(7));
t1708 = sin(qJ(2,3));
t1768 = t1700 * t1708;
t1803 = pkin(3) * (t1709 * t1764 + t1742 * t1768);
t1694 = pkin(6) - t1703;
t1741 = pkin(1) * t1711 - t1717 * t1694;
t1710 = sin(qJ(2,2));
t1767 = t1700 * t1710;
t1802 = pkin(3) * (t1711 * t1764 + t1741 * t1767);
t1695 = pkin(6) - t1704;
t1740 = pkin(1) * t1713 - t1719 * t1695;
t1712 = sin(qJ(2,1));
t1766 = t1700 * t1712;
t1801 = pkin(3) * (t1713 * t1764 + t1740 * t1766);
t1800 = pkin(3) * cos(qJ(2,3) + pkin(7));
t1799 = pkin(3) * cos(qJ(2,2) + pkin(7));
t1798 = pkin(3) * cos(qJ(2,1) + pkin(7));
t1794 = t1700 * pkin(3);
t1793 = t1701 * pkin(3);
t1669 = pkin(2) + t1793;
t1757 = pkin(3) * t1768;
t1650 = 0.1e1 / (t1669 * t1714 - t1757);
t1677 = 0.1e1 / t1693;
t1792 = t1650 * t1677;
t1756 = pkin(3) * t1767;
t1651 = 0.1e1 / (t1669 * t1716 - t1756);
t1678 = 0.1e1 / t1694;
t1791 = t1651 * t1678;
t1755 = pkin(3) * t1766;
t1652 = 0.1e1 / (t1669 * t1718 - t1755);
t1679 = 0.1e1 / t1695;
t1790 = t1652 * t1679;
t1686 = t1714 * pkin(2);
t1660 = 0.1e1 / (t1686 + t1800);
t1789 = t1660 * t1680;
t1788 = t1660 * t1683;
t1687 = t1716 * pkin(2);
t1661 = 0.1e1 / (t1687 + t1799);
t1787 = t1661 * t1681;
t1786 = t1661 * t1684;
t1688 = t1718 * pkin(2);
t1662 = 0.1e1 / (t1688 + t1798);
t1785 = t1662 * t1682;
t1784 = t1662 * t1685;
t1783 = t1669 * t1683;
t1782 = t1669 * t1684;
t1781 = t1669 * t1685;
t1780 = t1677 * t1715;
t1779 = t1678 * t1717;
t1778 = t1679 * t1719;
t1777 = t1680 * t1669;
t1776 = t1680 * t1709;
t1775 = t1681 * t1669;
t1774 = t1681 * t1711;
t1773 = t1682 * t1669;
t1772 = t1682 * t1713;
t1771 = t1683 * t1709;
t1770 = t1684 * t1711;
t1769 = t1685 * t1713;
t1765 = pkin(2) * t1793;
t1763 = t1683 * t1794;
t1762 = t1684 * t1794;
t1761 = t1685 * t1794;
t1760 = t1680 * t1794;
t1759 = t1681 * t1794;
t1758 = t1682 * t1794;
t1754 = t1644 * t1677 * t1714;
t1753 = t1645 * t1678 * t1716;
t1752 = t1646 * t1679 * t1718;
t1751 = t1644 * t1792;
t1750 = t1645 * t1791;
t1749 = t1646 * t1790;
t1729 = g(3) * t1715 + t1733 * t1709;
t1748 = t1729 * t1792;
t1728 = g(3) * t1717 + t1732 * t1711;
t1747 = t1728 * t1791;
t1727 = g(3) * t1719 + t1731 * t1713;
t1746 = t1727 * t1790;
t1745 = t1644 * t1780;
t1744 = t1645 * t1779;
t1743 = t1646 * t1778;
t1739 = t1650 * t1754;
t1738 = t1651 * t1753;
t1737 = t1652 * t1752;
t1736 = t1708 * t1751;
t1735 = t1710 * t1750;
t1734 = t1712 * t1749;
t1720 = pkin(3) ^ 2;
t1721 = pkin(2) ^ 2;
t1730 = 0.2e1 * t1696 * t1720 - t1720 + t1721 + 0.2e1 * t1765;
t1657 = t1680 * g(1) + t1683 * g(2);
t1623 = -t1657 * t1714 + t1708 * t1729;
t1658 = t1681 * g(1) + t1684 * g(2);
t1625 = -t1658 * t1716 + t1710 * t1728;
t1659 = t1682 * g(1) + t1685 * g(2);
t1627 = -t1659 * t1718 + t1712 * t1727;
t1726 = t1623 * t1789 + t1625 * t1787 + t1627 * t1785;
t1725 = t1623 * t1788 + t1625 * t1786 + t1627 * t1784;
t1724 = t1727 * t1778 + t1728 * t1779 + t1729 * t1780;
t1614 = (-t1669 * t1776 + t1763) * t1714 + (t1709 * t1760 + t1783) * t1708;
t1615 = (-t1669 * t1774 + t1762) * t1716 + (t1711 * t1759 + t1782) * t1710;
t1616 = (-t1669 * t1772 + t1761) * t1718 + (t1713 * t1758 + t1781) * t1712;
t1723 = t1614 * t1748 + t1615 * t1747 + t1616 * t1746;
t1617 = (t1669 * t1771 + t1760) * t1714 + t1708 * (-t1709 * t1763 + t1777);
t1618 = (t1669 * t1770 + t1759) * t1716 + t1710 * (-t1711 * t1762 + t1775);
t1619 = (t1669 * t1769 + t1758) * t1718 + t1712 * (-t1713 * t1761 + t1773);
t1722 = t1617 * t1748 + t1618 * t1747 + t1619 * t1746;
t1673 = pkin(1) * t1794;
t1672 = t1688 + pkin(1);
t1671 = t1687 + pkin(1);
t1670 = t1686 + pkin(1);
t1668 = t1671 * t1717;
t1667 = t1670 * t1715;
t1666 = t1719 * t1672;
t1665 = pkin(1) * t1712 - t1794;
t1664 = pkin(1) * t1710 - t1794;
t1663 = pkin(1) * t1708 - t1794;
t1656 = t1765 + t1721 / 0.2e1 + (t1696 - 0.1e1 / 0.2e1) * t1720;
t1637 = -0.2e1 * t1713 * t1755 + t1740;
t1636 = -0.2e1 * t1711 * t1756 + t1741;
t1635 = -0.2e1 * t1709 * t1757 + t1742;
t1634 = t1730 * t1712 + t1673;
t1633 = t1730 * t1710 + t1673;
t1632 = t1730 * t1708 + t1673;
t1628 = t1659 * t1712 + t1718 * t1727;
t1626 = t1658 * t1710 + t1716 * t1728;
t1624 = t1657 * t1708 + t1714 * t1729;
t1622 = g(3) * (t1709 * t1670 + t1702 * t1715) - t1733 * (-t1709 * t1702 + t1667);
t1621 = (t1713 * t1672 + t1704 * t1719) * g(3) - t1731 * (-t1713 * t1704 + t1666);
t1620 = (t1711 * t1671 + t1703 * t1717) * g(3) - t1732 * (-t1711 * t1703 + t1668);
t1 = [0, t1617 * t1751 + t1618 * t1750 + t1619 * t1749, t1722, 0, 0, 0, 0, 0, t1617 * t1739 + t1618 * t1738 + t1619 * t1737 + t1726, -t1617 * t1736 - t1618 * t1735 - t1619 * t1734 + t1624 * t1789 + t1626 * t1787 + t1628 * t1785, -t1722, (t1619 * t1621 - ((t1656 * t1769 + t1669 * t1758) * t1804 + (t1682 * t1634 + t1637 * t1781) * t1718 - t1685 * t1801 + t1665 * t1773) * t1646) * t1790 + (t1618 * t1620 - ((t1656 * t1770 + t1669 * t1759) * t1805 + (t1681 * t1633 + t1636 * t1782) * t1716 - t1684 * t1802 + t1664 * t1775) * t1645) * t1791 + (t1617 * t1622 - ((t1656 * t1771 + t1669 * t1760) * t1806 + (t1680 * t1632 + t1635 * t1783) * t1714 - t1683 * t1803 + t1663 * t1777) * t1644) * t1792 + t1726 * pkin(2), -g(1); 0, t1614 * t1751 + t1615 * t1750 + t1616 * t1749, t1723, 0, 0, 0, 0, 0, t1614 * t1739 + t1615 * t1738 + t1616 * t1737 + t1725, -t1614 * t1736 - t1615 * t1735 - t1616 * t1734 + t1624 * t1788 + t1626 * t1786 + t1628 * t1784, -t1723, (t1616 * t1621 - ((-t1656 * t1772 + t1669 * t1761) * t1804 + (t1634 * t1685 - t1637 * t1773) * t1718 + t1682 * t1801 + t1665 * t1781) * t1646) * t1790 + (t1615 * t1620 - ((-t1656 * t1774 + t1669 * t1762) * t1805 + (t1633 * t1684 - t1636 * t1775) * t1716 + t1681 * t1802 + t1664 * t1782) * t1645) * t1791 + (t1614 * t1622 - ((-t1656 * t1776 + t1669 * t1763) * t1806 + (t1632 * t1683 - t1635 * t1777) * t1714 + t1680 * t1803 + t1663 * t1783) * t1644) * t1792 + t1725 * pkin(2), -g(2); 0, t1743 + t1744 + t1745, t1724, 0, 0, 0, 0, 0, t1715 * t1754 + t1717 * t1753 + t1719 * t1752, -t1708 * t1745 - t1710 * t1744 - t1712 * t1743, -t1724, (t1719 * t1621 - (t1713 * t1695 + t1719 * t1798 + t1666) * t1646) * t1679 + (t1717 * t1620 - (t1711 * t1694 + t1717 * t1799 + t1668) * t1645) * t1678 + (t1715 * t1622 - (t1709 * t1693 + t1715 * t1800 + t1667) * t1644) * t1677, -g(3);];
tau_reg  = t1;
