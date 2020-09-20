% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P4PRRRR8V2G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% xDP [4x1]
%   Generalized platform velocities
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [4x15]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:13
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P4PRRRR8V2G1A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G1A0_coriolisvec_para_pf_regmin: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR8V2G1A0_coriolisvec_para_pf_regmin: xDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G1A0_coriolisvec_para_pf_regmin: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G1A0_coriolisvec_para_pf_regmin: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G1A0_coriolisvec_para_pf_regmin: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G1A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:11:17
% EndTime: 2020-08-07 11:11:43
% DurationCPUTime: 28.13s
% Computational Cost: add. (172801->819), mult. (352068->1624), div. (11204->22), fcn. (325172->30), ass. (0->642)
t1792 = cos(qJ(3,4));
t1746 = pkin(3) * t1792 + pkin(2);
t1793 = cos(qJ(2,4));
t1809 = pkin(7) + pkin(6);
t1747 = t1793 * t1809;
t1791 = sin(qJ(2,4));
t1687 = t1746 * t1791 - t1747;
t1785 = cos(pkin(4));
t1790 = sin(qJ(3,4));
t2076 = t1785 * t1790;
t1710 = t1746 * t2076;
t1783 = sin(pkin(4));
t2091 = t1783 * t1792;
t1641 = t1687 * t2091 + t1710;
t1638 = 0.1e1 / t1641;
t1800 = cos(qJ(3,3));
t1748 = pkin(3) * t1800 + pkin(2);
t1801 = cos(qJ(2,3));
t1751 = t1801 * t1809;
t1795 = sin(qJ(2,3));
t1691 = t1748 * t1795 - t1751;
t1794 = sin(qJ(3,3));
t2074 = t1785 * t1794;
t1711 = t1748 * t2074;
t2083 = t1783 * t1800;
t1651 = t1691 * t2083 + t1711;
t1642 = 0.1e1 / t1651;
t1802 = cos(qJ(3,2));
t1749 = pkin(3) * t1802 + pkin(2);
t1803 = cos(qJ(2,2));
t1752 = t1803 * t1809;
t1797 = sin(qJ(2,2));
t1692 = t1749 * t1797 - t1752;
t1796 = sin(qJ(3,2));
t2072 = t1785 * t1796;
t1712 = t1749 * t2072;
t2081 = t1783 * t1802;
t1652 = t1692 * t2081 + t1712;
t1644 = 0.1e1 / t1652;
t1804 = cos(qJ(3,1));
t1750 = pkin(3) * t1804 + pkin(2);
t1805 = cos(qJ(2,1));
t1753 = t1805 * t1809;
t1799 = sin(qJ(2,1));
t1693 = t1750 * t1799 - t1753;
t1798 = sin(qJ(3,1));
t2070 = t1785 * t1798;
t1713 = t1750 * t2070;
t2079 = t1783 * t1804;
t1653 = t1693 * t2079 + t1713;
t1646 = 0.1e1 / t1653;
t1722 = pkin(2) * t1791 - t1747;
t2092 = t1783 * t1791;
t1774 = t1792 ^ 2;
t2177 = pkin(3) * t1774;
t1614 = 0.1e1 / ((pkin(3) * t2076 + t1722 * t1783) * t1792 + pkin(2) * t2076 + t2092 * t2177);
t1724 = pkin(2) * t1795 - t1751;
t2088 = t1783 * t1795;
t1778 = t1800 ^ 2;
t2176 = pkin(3) * t1778;
t1615 = 0.1e1 / ((pkin(3) * t2074 + t1724 * t1783) * t1800 + pkin(2) * t2074 + t2088 * t2176);
t1725 = pkin(2) * t1797 - t1752;
t2086 = t1783 * t1797;
t1779 = t1802 ^ 2;
t2175 = pkin(3) * t1779;
t1616 = 0.1e1 / ((pkin(3) * t2072 + t1725 * t1783) * t1802 + pkin(2) * t2072 + t2086 * t2175);
t1726 = pkin(2) * t1799 - t1753;
t2084 = t1783 * t1799;
t1780 = t1804 ^ 2;
t2174 = pkin(3) * t1780;
t1617 = 0.1e1 / ((pkin(3) * t2070 + t1726 * t1783) * t1804 + pkin(2) * t2070 + t2084 * t2174);
t1810 = xP(4);
t1771 = sin(t1810);
t1772 = cos(t1810);
t1811 = koppelP(4,2);
t1815 = koppelP(4,1);
t1706 = -t1771 * t1811 + t1772 * t1815;
t1806 = xDP(4);
t1807 = xDP(2);
t1658 = -t1706 * t1806 - t1807;
t1702 = t1771 * t1815 + t1772 * t1811;
t1808 = xDP(1);
t1661 = t1702 * t1806 - t1808;
t1786 = legFrame(4,3);
t1759 = sin(t1786);
t1763 = cos(t1786);
t1782 = sin(pkin(8));
t1784 = cos(pkin(8));
t1670 = -t1759 * t1782 + t1763 * t1784;
t1674 = t1759 * t1784 + t1763 * t1782;
t2075 = t1785 * t1791;
t1832 = t1670 * t2075 + t1674 * t1793;
t1833 = -t1670 * t1793 + t1674 * t2075;
t1559 = ((t1658 * t1674 + t1661 * t1670) * t2091 + (t1658 * t1833 + t1661 * t1832) * t1790) * t1638;
t2214 = 0.2e1 * t1559;
t1813 = koppelP(2,2);
t1817 = koppelP(2,1);
t1708 = -t1771 * t1813 + t1772 * t1817;
t1659 = -t1708 * t1806 - t1807;
t1704 = t1771 * t1817 + t1772 * t1813;
t1662 = t1704 * t1806 - t1808;
t1788 = legFrame(2,3);
t1761 = sin(t1788);
t1765 = cos(t1788);
t1672 = -t1761 * t1782 + t1765 * t1784;
t1676 = t1761 * t1784 + t1765 * t1782;
t2071 = t1785 * t1797;
t1828 = t1672 * t2071 + t1676 * t1803;
t1829 = -t1672 * t1803 + t1676 * t2071;
t1564 = ((t1659 * t1676 + t1662 * t1672) * t2081 + (t1659 * t1829 + t1662 * t1828) * t1796) * t1644;
t2213 = 0.2e1 * t1564;
t1814 = koppelP(1,2);
t1818 = koppelP(1,1);
t1709 = -t1771 * t1814 + t1772 * t1818;
t1660 = -t1709 * t1806 - t1807;
t1705 = t1771 * t1818 + t1772 * t1814;
t1663 = t1705 * t1806 - t1808;
t1789 = legFrame(1,3);
t1762 = sin(t1789);
t1766 = cos(t1789);
t1673 = -t1762 * t1782 + t1766 * t1784;
t1677 = t1762 * t1784 + t1766 * t1782;
t2069 = t1785 * t1799;
t1826 = t1673 * t2069 + t1677 * t1805;
t1827 = -t1673 * t1805 + t1677 * t2069;
t1565 = ((t1660 * t1677 + t1663 * t1673) * t2079 + (t1660 * t1827 + t1663 * t1826) * t1798) * t1646;
t2212 = 0.2e1 * t1565;
t1812 = koppelP(3,2);
t1816 = koppelP(3,1);
t1707 = -t1771 * t1812 + t1772 * t1816;
t1664 = -t1707 * t1806 - t1807;
t1703 = t1771 * t1816 + t1772 * t1812;
t1665 = t1703 * t1806 - t1808;
t1787 = legFrame(3,3);
t1760 = sin(t1787);
t1764 = cos(t1787);
t1671 = -t1760 * t1782 + t1764 * t1784;
t1675 = t1760 * t1784 + t1764 * t1782;
t2073 = t1785 * t1795;
t1830 = t1671 * t2073 + t1675 * t1801;
t1831 = -t1671 * t1801 + t1675 * t2073;
t1566 = ((t1664 * t1675 + t1665 * t1671) * t2083 + (t1664 * t1831 + t1665 * t1830) * t1794) * t1642;
t2211 = 0.2e1 * t1566;
t2144 = t1559 * t1809;
t1925 = t1790 * t2144;
t1610 = t1658 * t1782 + t1661 * t1784;
t1853 = t1658 * t1784 - t1661 * t1782;
t2062 = t1791 * t1809;
t2116 = (t1746 * t1793 + t2062) * t1785;
t1571 = (t1610 * t1763 + t1759 * t1853) * t2116 - t1687 * (t1759 * t1610 - t1763 * t1853);
t2134 = t1571 * t1638;
t1536 = t1925 - t2134;
t1594 = -t1670 * t2091 - t1790 * t1832;
t1595 = -t1674 * t2091 - t1790 * t1833;
t1781 = t1806 ^ 2;
t1820 = 0.1e1 / pkin(3);
t1981 = t1820 * t2134;
t1890 = t1793 * t1981;
t2063 = t1791 * t1792;
t1929 = t1783 * t2063;
t2124 = t1638 * t1790;
t1982 = t1571 * t2124;
t2090 = t1783 * t1793;
t2182 = pkin(2) * t1793;
t1494 = (-((t1559 * t2090 + t1785 * t1981) * t2177 + ((-t1982 + t2144) * t1791 + t1559 * t2182) * t2091 + t1785 * t1536) * t1559 - (t1783 * t1890 + (t1774 * t1785 - t1790 * t1929 - t1785) * t1559) * t2134 + (-t1594 * t1706 - t1595 * t1702) * t1781) * t1614;
t1602 = -t1670 * t2116 + t1674 * t1687;
t1603 = -t1670 * t1687 - t1674 * t2116;
t1666 = pkin(3) * t2063 + t1722;
t1626 = 0.1e1 / (t1666 * t2091 + t1710);
t1754 = pkin(2) ^ 2 + t1809 ^ 2;
t1819 = pkin(3) ^ 2;
t2173 = 0.2e1 * pkin(2) * pkin(3);
t1993 = (-t1809 * t1982 + (t1774 * t1819 + t1792 * t2173 + t1754) * t1559) * t1559 * t1614;
t2067 = t1785 * t1820;
t2077 = t1783 * t1820;
t2094 = t1781 * t1820;
t2178 = pkin(2) * t1820;
t1498 = -t1993 * t2067 - (-t1785 * t1925 + (-t1666 * t1790 * t2077 + t1785 * (t1792 * t2178 + t1774)) * t2134) * t1626 * t1981 + (-t1602 * t1706 - t1603 * t1702) * t1638 * t2094;
t1480 = t1494 * t2090 + t1498 * t1785;
t1558 = t1559 ^ 2;
t1821 = 0.1e1 / pkin(3) ^ 2;
t2136 = t1571 ^ 2 / t1641 ^ 2;
t1983 = t1821 * t2136;
t1543 = t1558 + t1983;
t1845 = -0.2e1 * t1559 * t1890;
t2202 = t1480 * t1792 + (-t1543 * t2063 + t1790 * t1845) * t1783;
t2140 = t1566 * t1809;
t1924 = t1794 * t2140;
t1613 = t1664 * t1782 + t1665 * t1784;
t1850 = t1664 * t1784 - t1665 * t1782;
t2057 = t1795 * t1809;
t2115 = (t1748 * t1801 + t2057) * t1785;
t1575 = (t1613 * t1764 + t1760 * t1850) * t2115 - t1691 * (t1760 * t1613 - t1764 * t1850);
t2129 = t1575 * t1642;
t1544 = t1924 - t2129;
t1596 = -t1671 * t2083 - t1794 * t1830;
t1599 = -t1675 * t2083 - t1794 * t1831;
t1970 = t1820 * t2129;
t1880 = t1801 * t1970;
t2058 = t1795 * t1800;
t1928 = t1783 * t2058;
t2122 = t1642 * t1794;
t1971 = t1575 * t2122;
t2082 = t1783 * t1801;
t2181 = pkin(2) * t1801;
t1495 = (-((t1566 * t2082 + t1785 * t1970) * t2176 + ((-t1971 + t2140) * t1795 + t1566 * t2181) * t2083 + t1785 * t1544) * t1566 - (t1783 * t1880 + (t1778 * t1785 - t1794 * t1928 - t1785) * t1566) * t2129 + (-t1596 * t1707 - t1599 * t1703) * t1781) * t1615;
t1604 = -t1671 * t2115 + t1675 * t1691;
t1607 = -t1671 * t1691 - t1675 * t2115;
t1667 = pkin(3) * t2058 + t1724;
t1627 = 0.1e1 / (t1667 * t2083 + t1711);
t1992 = (-t1809 * t1971 + (t1778 * t1819 + t1800 * t2173 + t1754) * t1566) * t1566 * t1615;
t1499 = -t1992 * t2067 - (-t1785 * t1924 + (-t1667 * t1794 * t2077 + t1785 * (t1800 * t2178 + t1778)) * t2129) * t1627 * t1970 + (-t1604 * t1707 - t1607 * t1703) * t1642 * t2094;
t1487 = t1495 * t2082 + t1499 * t1785;
t1562 = t1566 ^ 2;
t2133 = t1575 ^ 2 / t1651 ^ 2;
t1978 = t1821 * t2133;
t1549 = t1562 + t1978;
t1842 = -0.2e1 * t1566 * t1880;
t2201 = t1487 * t1800 + (-t1549 * t2058 + t1794 * t1842) * t1783;
t2142 = t1564 * t1809;
t1923 = t1796 * t2142;
t1611 = t1659 * t1782 + t1662 * t1784;
t1852 = t1659 * t1784 - t1662 * t1782;
t2053 = t1797 * t1809;
t2114 = (t1749 * t1803 + t2053) * t1785;
t1576 = (t1611 * t1765 + t1761 * t1852) * t2114 - t1692 * (t1761 * t1611 - t1765 * t1852);
t2127 = t1576 * t1644;
t1545 = t1923 - t2127;
t1597 = -t1672 * t2081 - t1796 * t1828;
t1600 = -t1676 * t2081 - t1796 * t1829;
t1968 = t1820 * t2127;
t1879 = t1803 * t1968;
t2054 = t1797 * t1802;
t1927 = t1783 * t2054;
t2120 = t1644 * t1796;
t1969 = t1576 * t2120;
t2080 = t1783 * t1803;
t2180 = pkin(2) * t1803;
t1496 = (-((t1564 * t2080 + t1785 * t1968) * t2175 + ((-t1969 + t2142) * t1797 + t1564 * t2180) * t2081 + t1785 * t1545) * t1564 - (t1783 * t1879 + (t1779 * t1785 - t1796 * t1927 - t1785) * t1564) * t2127 + (-t1597 * t1708 - t1600 * t1704) * t1781) * t1616;
t1605 = -t1672 * t2114 + t1676 * t1692;
t1608 = -t1672 * t1692 - t1676 * t2114;
t1668 = pkin(3) * t2054 + t1725;
t1628 = 0.1e1 / (t1668 * t2081 + t1712);
t1991 = (-t1809 * t1969 + (t1779 * t1819 + t1802 * t2173 + t1754) * t1564) * t1564 * t1616;
t1500 = -t1991 * t2067 - (-t1785 * t1923 + (-t1668 * t1796 * t2077 + t1785 * (t1802 * t2178 + t1779)) * t2127) * t1628 * t1968 + (-t1605 * t1708 - t1608 * t1704) * t1644 * t2094;
t1488 = t1496 * t2080 + t1500 * t1785;
t1560 = t1564 ^ 2;
t2132 = t1576 ^ 2 / t1652 ^ 2;
t1975 = t1821 * t2132;
t1547 = t1560 + t1975;
t1844 = -0.2e1 * t1564 * t1879;
t2200 = t1488 * t1802 + (-t1547 * t2054 + t1796 * t1844) * t1783;
t2141 = t1565 * t1809;
t1922 = t1798 * t2141;
t1612 = t1660 * t1782 + t1663 * t1784;
t1851 = t1660 * t1784 - t1663 * t1782;
t2049 = t1799 * t1809;
t2113 = (t1750 * t1805 + t2049) * t1785;
t1577 = (t1612 * t1766 + t1762 * t1851) * t2113 - t1693 * (t1762 * t1612 - t1766 * t1851);
t2125 = t1577 * t1646;
t1546 = t1922 - t2125;
t1598 = -t1673 * t2079 - t1798 * t1826;
t1601 = -t1677 * t2079 - t1798 * t1827;
t1966 = t1820 * t2125;
t1878 = t1805 * t1966;
t2050 = t1799 * t1804;
t1926 = t1783 * t2050;
t2118 = t1646 * t1798;
t1967 = t1577 * t2118;
t2078 = t1783 * t1805;
t2179 = pkin(2) * t1805;
t1497 = (-((t1565 * t2078 + t1785 * t1966) * t2174 + ((-t1967 + t2141) * t1799 + t1565 * t2179) * t2079 + t1785 * t1546) * t1565 - (t1783 * t1878 + (t1780 * t1785 - t1798 * t1926 - t1785) * t1565) * t2125 + (-t1598 * t1709 - t1601 * t1705) * t1781) * t1617;
t1606 = -t1673 * t2113 + t1677 * t1693;
t1609 = -t1673 * t1693 - t1677 * t2113;
t1669 = pkin(3) * t2050 + t1726;
t1629 = 0.1e1 / (t1669 * t2079 + t1713);
t1990 = (-t1809 * t1967 + (t1780 * t1819 + t1804 * t2173 + t1754) * t1565) * t1565 * t1617;
t1501 = -t1990 * t2067 - (-t1785 * t1922 + (-t1669 * t1798 * t2077 + t1785 * (t1804 * t2178 + t1780)) * t2125) * t1629 * t1966 + (-t1606 * t1709 - t1609 * t1705) * t1646 * t2094;
t1489 = t1497 * t2078 + t1501 * t1785;
t1561 = t1565 ^ 2;
t2131 = t1577 ^ 2 / t1653 ^ 2;
t1972 = t1821 * t2131;
t1548 = t1561 + t1972;
t1843 = -0.2e1 * t1565 * t1878;
t2199 = t1489 * t1804 + (-t1548 * t2050 + t1798 * t1843) * t1783;
t2065 = t1790 * t1791;
t2198 = -t1480 * t1790 + (t1543 * t2065 + t1792 * t1845) * t1783;
t2056 = t1796 * t1797;
t2197 = -t1488 * t1796 + (t1547 * t2056 + t1802 * t1844) * t1783;
t2052 = t1798 * t1799;
t2196 = -t1489 * t1798 + (t1548 * t2052 + t1804 * t1843) * t1783;
t2060 = t1794 * t1795;
t2195 = -t1487 * t1794 + (t1549 * t2060 + t1800 * t1842) * t1783;
t1729 = t2049 + t2179;
t2085 = t1783 * t1798;
t1834 = pkin(3) * t2085 - t1726 * t1785;
t2194 = t1729 * t1784 + t1782 * t1834;
t1728 = t2053 + t2180;
t2087 = t1783 * t1796;
t1835 = pkin(3) * t2087 - t1725 * t1785;
t2193 = t1728 * t1784 + t1782 * t1835;
t1727 = t2057 + t2181;
t2089 = t1783 * t1794;
t1836 = pkin(3) * t2089 - t1724 * t1785;
t2192 = t1727 * t1784 + t1782 * t1836;
t1723 = t2062 + t2182;
t2093 = t1783 * t1790;
t1837 = pkin(3) * t2093 - t1722 * t1785;
t2191 = t1723 * t1784 + t1782 * t1837;
t2190 = 0.2e1 * t1494;
t2189 = 0.2e1 * t1495;
t2188 = 0.2e1 * t1496;
t2187 = 0.2e1 * t1497;
t1755 = 0.2e1 * t1774 - 0.1e1;
t1756 = 0.2e1 * t1778 - 0.1e1;
t1757 = 0.2e1 * t1779 - 0.1e1;
t1758 = 0.2e1 * t1780 - 0.1e1;
t2186 = pkin(2) * t1559;
t2185 = pkin(2) * t1564;
t2184 = pkin(2) * t1565;
t2183 = pkin(2) * t1566;
t2172 = t1494 * t1614;
t2171 = t1495 * t1615;
t2170 = t1496 * t1616;
t2169 = t1497 * t1617;
t2168 = t1498 * t1638;
t2167 = t1498 * t1790;
t2166 = t1498 * t1792;
t2165 = t1499 * t1642;
t2164 = t1499 * t1794;
t2163 = t1499 * t1800;
t2162 = t1500 * t1644;
t2161 = t1500 * t1796;
t2160 = t1500 * t1802;
t2159 = t1501 * t1646;
t2158 = t1501 * t1798;
t2157 = t1501 * t1804;
t1618 = t1782 * t1723 - t1784 * t1837;
t1678 = t1782 * t2075 - t1784 * t1793;
t1679 = t1782 * t1793 + t1784 * t2075;
t2025 = pkin(2) * t2093;
t1582 = -(t1678 * t1763 + t1679 * t1759) * t2177 + (-t1759 * t1618 + t2191 * t1763) * t1792 + t1674 * t2025;
t1583 = (-t1678 * t1759 + t1679 * t1763) * t2177 + (t1618 * t1763 + t2191 * t1759) * t1792 - t1670 * t2025;
t1825 = (-t1582 * t1706 - t1583 * t1702) * t1614;
t2043 = t1792 * t1993 + (pkin(2) * t1981 - t1536 * t1792) * t1614 * t2134;
t1502 = t1781 * t1825 + t2043;
t2156 = t1502 * t1614;
t2155 = t1502 * t1791;
t2154 = t1502 * t1793;
t1619 = t1782 * t1727 - t1784 * t1836;
t1680 = t1782 * t2073 - t1784 * t1801;
t1683 = t1782 * t1801 + t1784 * t2073;
t2024 = pkin(2) * t2089;
t1584 = -(t1680 * t1764 + t1683 * t1760) * t2176 + (-t1760 * t1619 + t2192 * t1764) * t1800 + t1675 * t2024;
t1587 = (-t1680 * t1760 + t1683 * t1764) * t2176 + (t1619 * t1764 + t2192 * t1760) * t1800 - t1671 * t2024;
t1824 = (-t1584 * t1707 - t1587 * t1703) * t1615;
t2036 = t1800 * t1992 + (pkin(2) * t1970 - t1544 * t1800) * t1615 * t2129;
t1503 = t1781 * t1824 + t2036;
t2153 = t1503 * t1615;
t2152 = t1503 * t1795;
t2151 = t1503 * t1801;
t1620 = t1782 * t1728 - t1784 * t1835;
t1681 = t1782 * t2071 - t1784 * t1803;
t1684 = t1782 * t1803 + t1784 * t2071;
t2023 = pkin(2) * t2087;
t1585 = -(t1681 * t1765 + t1684 * t1761) * t2175 + (-t1761 * t1620 + t2193 * t1765) * t1802 + t1676 * t2023;
t1588 = (-t1681 * t1761 + t1684 * t1765) * t2175 + (t1620 * t1765 + t2193 * t1761) * t1802 - t1672 * t2023;
t1823 = (-t1585 * t1708 - t1588 * t1704) * t1616;
t2035 = t1802 * t1991 + (pkin(2) * t1968 - t1545 * t1802) * t1616 * t2127;
t1504 = t1781 * t1823 + t2035;
t2150 = t1504 * t1616;
t2149 = t1504 * t1797;
t2148 = t1504 * t1803;
t1621 = t1782 * t1729 - t1784 * t1834;
t1682 = t1782 * t2069 - t1784 * t1805;
t1685 = t1782 * t1805 + t1784 * t2069;
t2022 = pkin(2) * t2085;
t1586 = -(t1682 * t1766 + t1685 * t1762) * t2174 + (-t1762 * t1621 + t2194 * t1766) * t1804 + t1677 * t2022;
t1589 = (-t1682 * t1762 + t1685 * t1766) * t2174 + (t1621 * t1766 + t2194 * t1762) * t1804 - t1673 * t2022;
t1822 = (-t1586 * t1709 - t1589 * t1705) * t1617;
t2034 = t1804 * t1990 + (pkin(2) * t1966 - t1546 * t1804) * t1617 * t2125;
t1505 = t1781 * t1822 + t2034;
t2147 = t1505 * t1617;
t2146 = t1505 * t1799;
t2145 = t1505 * t1805;
t1714 = t1782 * t1811 + t1784 * t1815;
t1715 = -t1782 * t1815 + t1784 * t1811;
t1654 = t1714 * t1771 + t1715 * t1772;
t1849 = t1714 * t1772 - t1715 * t1771;
t1590 = t1654 * t1763 - t1759 * t1849;
t1694 = t1793 * t1811 - t1815 * t2075;
t1695 = t1793 * t1815 + t1811 * t2075;
t1630 = t1694 * t1784 - t1695 * t1782;
t1631 = t1694 * t1782 + t1695 * t1784;
t2143 = (((-t1630 * t1771 + t1631 * t1772) * t1763 + (t1630 * t1772 + t1631 * t1771) * t1759) * t1790 + t1590 * t2091) * t1626;
t1716 = t1782 * t1812 + t1784 * t1816;
t1717 = -t1782 * t1816 + t1784 * t1812;
t1655 = t1716 * t1771 + t1717 * t1772;
t1848 = t1716 * t1772 - t1717 * t1771;
t1591 = t1655 * t1764 - t1760 * t1848;
t1696 = t1801 * t1812 - t1816 * t2073;
t1699 = t1801 * t1816 + t1812 * t2073;
t1632 = t1696 * t1784 - t1699 * t1782;
t1635 = t1696 * t1782 + t1699 * t1784;
t2139 = (((-t1632 * t1771 + t1635 * t1772) * t1764 + (t1632 * t1772 + t1635 * t1771) * t1760) * t1794 + t1591 * t2083) * t1627;
t1718 = t1782 * t1813 + t1784 * t1817;
t1719 = -t1782 * t1817 + t1784 * t1813;
t1656 = t1718 * t1771 + t1719 * t1772;
t1847 = t1718 * t1772 - t1719 * t1771;
t1592 = t1656 * t1765 - t1761 * t1847;
t1697 = t1803 * t1813 - t1817 * t2071;
t1700 = t1803 * t1817 + t1813 * t2071;
t1633 = t1697 * t1784 - t1700 * t1782;
t1636 = t1697 * t1782 + t1700 * t1784;
t2138 = (((-t1633 * t1771 + t1636 * t1772) * t1765 + t1761 * (t1633 * t1772 + t1636 * t1771)) * t1796 + t1592 * t2081) * t1628;
t1720 = t1782 * t1814 + t1784 * t1818;
t1721 = -t1782 * t1818 + t1784 * t1814;
t1657 = t1720 * t1771 + t1721 * t1772;
t1846 = t1720 * t1772 - t1721 * t1771;
t1593 = t1657 * t1766 - t1762 * t1846;
t1698 = t1805 * t1814 - t1818 * t2069;
t1701 = t1805 * t1818 + t1814 * t2069;
t1634 = t1698 * t1784 - t1701 * t1782;
t1637 = t1698 * t1782 + t1701 * t1784;
t2137 = (((-t1634 * t1771 + t1637 * t1772) * t1766 + (t1634 * t1772 + t1637 * t1771) * t1762) * t1798 + t1593 * t2079) * t1629;
t2135 = t1571 * t1614;
t2130 = t1575 * t1615;
t2128 = t1576 * t1616;
t2126 = t1577 * t1617;
t2123 = t1638 * t1792;
t2121 = t1642 * t1800;
t2119 = t1644 * t1802;
t2117 = t1646 * t1804;
t2104 = t1771 * t1781;
t2099 = t1772 * t1781;
t2068 = t1785 * t1809;
t2066 = t1785 * t1821;
t2064 = t1790 * t1792;
t2061 = t1792 * t1793;
t2059 = t1794 * t1800;
t2055 = t1796 * t1802;
t2051 = t1798 * t1804;
t2048 = t1800 * t1801;
t2047 = t1802 * t1803;
t2046 = t1804 * t1805;
t1891 = t1785 * t1983;
t2045 = -t1498 * t1929 - t1792 * t1891 + t2198;
t2010 = t1498 * t2065;
t2044 = -t1783 * t2010 - t1790 * t1891 + t2202;
t1884 = t1785 * t1975;
t2042 = -t1500 * t1927 - t1802 * t1884 + t2197;
t1881 = t1785 * t1972;
t2041 = -t1501 * t1926 - t1804 * t1881 + t2196;
t1887 = t1785 * t1978;
t2040 = -t1499 * t1928 - t1800 * t1887 + t2195;
t2006 = t1499 * t2060;
t2039 = -t1783 * t2006 - t1794 * t1887 + t2201;
t2002 = t1500 * t2056;
t2038 = -t1783 * t2002 - t1796 * t1884 + t2200;
t1998 = t1501 * t2052;
t2037 = -t1783 * t1998 - t1798 * t1881 + t2199;
t2029 = -0.2e1 * t2135;
t2028 = -0.2e1 * t2130;
t2027 = -0.2e1 * t2128;
t2026 = -0.2e1 * t2126;
t2021 = t1494 * t2143;
t1773 = t1790 ^ 2;
t2020 = t1773 * t2172;
t2019 = t1495 * t2139;
t1775 = t1794 ^ 2;
t2018 = t1775 * t2171;
t2017 = t1496 * t2138;
t1776 = t1796 ^ 2;
t2016 = t1776 * t2170;
t2015 = t1497 * t2137;
t1777 = t1798 ^ 2;
t2014 = t1777 * t2169;
t2013 = t1498 * t2143;
t2012 = t1614 * t2167;
t2011 = t1614 * t2166;
t2009 = t1499 * t2139;
t2008 = t1615 * t2164;
t2007 = t1615 * t2163;
t2005 = t1500 * t2138;
t2004 = t1616 * t2161;
t2003 = t1616 * t2160;
t2001 = t1501 * t2137;
t2000 = t1617 * t2158;
t1999 = t1617 * t2157;
t1997 = t1502 * t2143;
t1996 = t1503 * t2139;
t1995 = t1504 * t2138;
t1994 = t1505 * t2137;
t1989 = t1571 * t2143;
t1988 = t1575 * t2139;
t1987 = t1576 * t2138;
t1986 = t1577 * t2137;
t1985 = t1790 * t2136;
t1984 = t1792 * t2136;
t1980 = t1794 * t2133;
t1979 = t1800 * t2133;
t1977 = t1796 * t2132;
t1976 = t1802 * t2132;
t1974 = t1798 * t2131;
t1973 = t1804 * t2131;
t1578 = t1590 * t2116 - t1687 * (t1654 * t1759 + t1763 * t1849);
t1965 = t1578 * t2124;
t1964 = t1578 * t2123;
t1579 = t1591 * t2115 - t1691 * (t1655 * t1760 + t1764 * t1848);
t1963 = t1579 * t2122;
t1962 = t1579 * t2121;
t1580 = t1592 * t2114 - (t1656 * t1761 + t1765 * t1847) * t1692;
t1961 = t1580 * t2120;
t1960 = t1580 * t2119;
t1581 = t1593 * t2113 - (t1762 * t1657 + t1766 * t1846) * t1693;
t1959 = t1581 * t2118;
t1958 = t1581 * t2117;
t1957 = t1602 * t2124;
t1956 = t1602 * t2123;
t1955 = t1603 * t2124;
t1954 = t1603 * t2123;
t1953 = t1604 * t2122;
t1952 = t1604 * t2121;
t1951 = t1605 * t2120;
t1950 = t1605 * t2119;
t1949 = t1606 * t2118;
t1948 = t1606 * t2117;
t1947 = t1607 * t2122;
t1946 = t1607 * t2121;
t1945 = t1608 * t2120;
t1944 = t1608 * t2119;
t1943 = t1609 * t2118;
t1942 = t1609 * t2117;
t1940 = t1638 * t2064;
t1938 = t1642 * t2059;
t1936 = t1644 * t2055;
t1934 = t1646 * t2051;
t1933 = t1746 * t2093;
t1932 = t1748 * t2089;
t1931 = t1749 * t2087;
t1930 = t1750 * t2085;
t1921 = t2064 * t2190;
t1920 = t2059 * t2189;
t1919 = t2055 * t2188;
t1918 = t2051 * t2187;
t1917 = t2135 * t2214;
t1916 = -0.2e1 * t1989;
t1915 = t2128 * t2213;
t1914 = t2126 * t2212;
t1913 = t2130 * t2211;
t1912 = -0.2e1 * t1988;
t1911 = -0.2e1 * t1987;
t1910 = -0.2e1 * t1986;
t1909 = t1594 * t2029;
t1908 = t1595 * t2029;
t1907 = t1596 * t2028;
t1906 = t1599 * t2028;
t1905 = t1597 * t2027;
t1904 = t1600 * t2027;
t1903 = t1598 * t2026;
t1902 = t1601 * t2026;
t1901 = pkin(6) * t1981;
t1900 = pkin(6) * t1970;
t1899 = pkin(6) * t1968;
t1898 = pkin(6) * t1966;
t1897 = t2136 * t2143;
t1896 = t2133 * t2139;
t1895 = t2132 * t2138;
t1894 = t2131 * t2137;
t1893 = t1614 * t1985;
t1892 = t1614 * t1984;
t1889 = t1615 * t1980;
t1888 = t1615 * t1979;
t1886 = t1616 * t1977;
t1885 = t1616 * t1976;
t1883 = t1617 * t1974;
t1882 = t1617 * t1973;
t1877 = t1746 * t1785;
t1876 = t1748 * t1785;
t1875 = t1749 * t1785;
t1874 = t1750 * t1785;
t1873 = t1614 * t1921;
t1872 = t1615 * t1920;
t1871 = t1616 * t1919;
t1870 = t1617 * t1918;
t1869 = t1989 * t2214;
t1868 = t1755 * t1917;
t1867 = t1987 * t2213;
t1866 = t1757 * t1915;
t1865 = t1986 * t2212;
t1864 = t1758 * t1914;
t1863 = t1988 * t2211;
t1862 = t1756 * t1913;
t1861 = t1494 * t1793 - t1558 * t1791;
t1860 = -t1494 * t1791 - t1558 * t1793;
t1859 = t1495 * t1801 - t1562 * t1795;
t1858 = -t1495 * t1795 - t1562 * t1801;
t1857 = t1496 * t1803 - t1560 * t1797;
t1856 = -t1496 * t1797 - t1560 * t1803;
t1855 = t1497 * t1805 - t1561 * t1799;
t1854 = -t1497 * t1799 - t1561 * t1805;
t1841 = pkin(2) * t2190 + t1502 * t2090;
t1840 = pkin(2) * t2189 + t1503 * t2082;
t1839 = pkin(2) * t2188 + t1504 * t2080;
t1838 = pkin(2) * t2187 + t1505 * t2078;
t1557 = (-t1705 * ((t1750 * t1673 + t1677 * t2068) * t2046 - (-t1673 * t1809 + t1677 * t1874) * t2050 + t1677 * t1930) + t1709 * ((-t1673 * t2068 + t1750 * t1677) * t2046 + (t1673 * t1874 + t1677 * t1809) * t2050 - t1673 * t1930)) / (-t1753 * t2079 + (t1926 + t2070) * t1750);
t1556 = (-t1704 * ((t1749 * t1672 + t1676 * t2068) * t2047 - (-t1672 * t1809 + t1676 * t1875) * t2054 + t1676 * t1931) + t1708 * ((-t1672 * t2068 + t1749 * t1676) * t2047 + (t1672 * t1875 + t1676 * t1809) * t2054 - t1672 * t1931)) / (-t1752 * t2081 + (t1927 + t2072) * t1749);
t1555 = (-t1703 * ((t1748 * t1671 + t1675 * t2068) * t2048 - (-t1671 * t1809 + t1675 * t1876) * t2058 + t1675 * t1932) + t1707 * ((-t1671 * t2068 + t1748 * t1675) * t2048 + (t1671 * t1876 + t1675 * t1809) * t2058 - t1671 * t1932)) / (-t1751 * t2083 + (t1928 + t2074) * t1748);
t1554 = (-t1702 * ((t1746 * t1670 + t1674 * t2068) * t2061 - (-t1670 * t1809 + t1674 * t1877) * t2063 + t1674 * t1933) + t1706 * ((-t1670 * t2068 + t1746 * t1674) * t2061 + (t1670 * t1877 + t1674 * t1809) * t2063 - t1670 * t1933)) / (-t1747 * t2091 + (t1929 + t2076) * t1746);
t1553 = t1758 * t1561;
t1552 = t1757 * t1560;
t1551 = t1756 * t1562;
t1550 = t1755 * t1558;
t1542 = t1804 * t2184 - t1798 * t1898 / 0.2e1;
t1541 = t1802 * t2185 - t1796 * t1899 / 0.2e1;
t1540 = t1800 * t2183 - t1794 * t1900 / 0.2e1;
t1539 = t1798 * t2184 + t1804 * t1898 / 0.2e1;
t1538 = t1796 * t2185 + t1802 * t1899 / 0.2e1;
t1537 = t1794 * t2183 + t1800 * t1900 / 0.2e1;
t1535 = t1792 * t2186 - t1790 * t1901 / 0.2e1;
t1534 = t1790 * t2186 + t1792 * t1901 / 0.2e1;
t1493 = pkin(6) * t1497 + t1505 * t2084;
t1492 = pkin(6) * t1496 + t1504 * t2086;
t1491 = pkin(6) * t1495 + t1503 * t2088;
t1490 = pkin(6) * t1494 + t1502 * t2092;
t1477 = -t1493 * t1804 - t1505 * t2070;
t1476 = -t1492 * t1802 - t1504 * t2072;
t1475 = -t1491 * t1800 - t1503 * t2074;
t1474 = t1505 * t1785 * t1804 - t1493 * t1798;
t1473 = t1504 * t1785 * t1802 - t1492 * t1796;
t1472 = t1503 * t1785 * t1800 - t1491 * t1794;
t1471 = -pkin(6) * t2157 - t1798 * t1838;
t1470 = -pkin(6) * t2158 + t1804 * t1838;
t1469 = -pkin(6) * t2163 - t1794 * t1840;
t1468 = -pkin(6) * t2164 + t1800 * t1840;
t1467 = -pkin(6) * t2160 - t1796 * t1839;
t1466 = -pkin(6) * t2161 + t1802 * t1839;
t1465 = -t1490 * t1792 - t1502 * t2076;
t1464 = t1502 * t1785 * t1792 - t1490 * t1790;
t1457 = -pkin(6) * t2166 - t1790 * t1841;
t1456 = -pkin(6) * t2167 + t1792 * t1841;
t1 = [t1582 * t2156 + t1584 * t2153 + t1585 * t2150 + t1586 * t2147, t1594 * t2172 + t1596 * t2171 + t1597 * t2170 + t1598 * t2169, ((t1586 * t1855 + t1598 * t2145) * t1617 + (t1585 * t1857 + t1597 * t2148) * t1616 + (t1584 * t1859 + t1596 * t2151) * t1615 + (t1582 * t1861 + t1594 * t2154) * t1614) * t1783, ((t1586 * t1854 - t1598 * t2146) * t1617 + (t1585 * t1856 - t1597 * t2149) * t1616 + (t1584 * t1858 - t1596 * t2152) * t1615 + (t1582 * t1860 - t1594 * t2155) * t1614) * t1783, t1594 * t2020 + t1596 * t2018 + t1597 * t2016 + t1598 * t2014 + ((-t1561 * t1606 + t1598 * t1914) * t1934 + (-t1560 * t1605 + t1597 * t1915) * t1936 + (-t1562 * t1604 + t1596 * t1913) * t1938 + (-t1558 * t1602 + t1594 * t1917) * t1940) * t1820, t1594 * t1873 + t1596 * t1872 + t1597 * t1871 + t1598 * t1870 + ((-t1553 * t1606 + t1598 * t1864) * t1646 + (-t1552 * t1605 + t1597 * t1866) * t1644 + (-t1551 * t1604 + t1596 * t1862) * t1642 + (-t1550 * t1602 + t1594 * t1868) * t1638) * t1820, t1594 * t2012 + t1596 * t2008 + t1597 * t2004 + t1598 * t2000 + (t1594 * t1892 + t1596 * t1888 + t1597 * t1885 + t1598 * t1882) * t1821 + (t1494 * t1957 + t1495 * t1953 + t1496 * t1951 + t1497 * t1949) * t1820, t1594 * t2011 + t1596 * t2007 + t1597 * t2003 + t1598 * t1999 + (-t1594 * t1893 - t1596 * t1889 - t1597 * t1886 - t1598 * t1883) * t1821 + (t1494 * t1956 + t1495 * t1952 + t1496 * t1950 + t1497 * t1948) * t1820, (t1602 * t2168 + t1604 * t2165 + t1605 * t2162 + t1606 * t2159) * t1820, (t1598 * t1470 + t1586 * t2037) * t1617 + (t1597 * t1466 + t1585 * t2038) * t1616 + (t1596 * t1468 + t1584 * t2039) * t1615 + (t1594 * t1456 + t1582 * t2044) * t1614 + ((t1474 * t1606 + t1539 * t1903) * t1646 + (t1473 * t1605 + t1538 * t1905) * t1644 + (t1472 * t1604 + t1537 * t1907) * t1642 + (t1464 * t1602 + t1534 * t1909) * t1638 + (t1558 * t1957 + t1560 * t1951 + t1561 * t1949 + t1562 * t1953) * pkin(2)) * t1820, (t1598 * t1471 + t1586 * t2041) * t1617 + (t1597 * t1467 + t1585 * t2042) * t1616 + (t1596 * t1469 + t1584 * t2040) * t1615 + (t1594 * t1457 + t1582 * t2045) * t1614 + ((t1477 * t1606 + t1542 * t1903) * t1646 + (t1476 * t1605 + t1541 * t1905) * t1644 + (t1475 * t1604 + t1540 * t1907) * t1642 + (t1465 * t1602 + t1535 * t1909) * t1638 + (t1558 * t1956 + t1560 * t1950 + t1561 * t1948 + t1562 * t1952) * pkin(2)) * t1820, 0, -t2099, t2104, 0; t1583 * t2156 + t1587 * t2153 + t1588 * t2150 + t1589 * t2147, t1595 * t2172 + t1599 * t2171 + t1600 * t2170 + t1601 * t2169, ((t1589 * t1855 + t1601 * t2145) * t1617 + (t1588 * t1857 + t1600 * t2148) * t1616 + (t1587 * t1859 + t1599 * t2151) * t1615 + (t1583 * t1861 + t1595 * t2154) * t1614) * t1783, ((t1589 * t1854 - t1601 * t2146) * t1617 + (t1588 * t1856 - t1600 * t2149) * t1616 + (t1587 * t1858 - t1599 * t2152) * t1615 + (t1583 * t1860 - t1595 * t2155) * t1614) * t1783, t1595 * t2020 + t1599 * t2018 + t1600 * t2016 + t1601 * t2014 + ((-t1561 * t1609 + t1601 * t1914) * t1934 + (-t1560 * t1608 + t1600 * t1915) * t1936 + (-t1562 * t1607 + t1599 * t1913) * t1938 + (-t1558 * t1603 + t1595 * t1917) * t1940) * t1820, t1595 * t1873 + t1599 * t1872 + t1600 * t1871 + t1601 * t1870 + ((-t1553 * t1609 + t1601 * t1864) * t1646 + (-t1552 * t1608 + t1600 * t1866) * t1644 + (-t1551 * t1607 + t1599 * t1862) * t1642 + (-t1550 * t1603 + t1595 * t1868) * t1638) * t1820, t1595 * t2012 + t1599 * t2008 + t1600 * t2004 + t1601 * t2000 + (t1595 * t1892 + t1599 * t1888 + t1600 * t1885 + t1601 * t1882) * t1821 + (t1494 * t1955 + t1495 * t1947 + t1496 * t1945 + t1497 * t1943) * t1820, t1595 * t2011 + t1599 * t2007 + t1600 * t2003 + t1601 * t1999 + (-t1595 * t1893 - t1599 * t1889 - t1600 * t1886 - t1601 * t1883) * t1821 + (t1494 * t1954 + t1495 * t1946 + t1496 * t1944 + t1497 * t1942) * t1820, (t1603 * t2168 + t1607 * t2165 + t1608 * t2162 + t1609 * t2159) * t1820, (t1601 * t1470 + t1589 * t2037) * t1617 + (t1600 * t1466 + t1588 * t2038) * t1616 + (t1599 * t1468 + t1587 * t2039) * t1615 + (t1595 * t1456 + t1583 * t2044) * t1614 + ((t1474 * t1609 + t1539 * t1902) * t1646 + (t1473 * t1608 + t1538 * t1904) * t1644 + (t1472 * t1607 + t1537 * t1906) * t1642 + (t1464 * t1603 + t1534 * t1908) * t1638 + (t1558 * t1955 + t1560 * t1945 + t1561 * t1943 + t1562 * t1947) * pkin(2)) * t1820, (t1601 * t1471 + t1589 * t2041) * t1617 + (t1600 * t1467 + t1588 * t2042) * t1616 + (t1599 * t1469 + t1587 * t2040) * t1615 + (t1595 * t1457 + t1583 * t2045) * t1614 + ((t1477 * t1609 + t1542 * t1902) * t1646 + (t1476 * t1608 + t1541 * t1904) * t1644 + (t1475 * t1607 + t1540 * t1906) * t1642 + (t1465 * t1603 + t1535 * t1908) * t1638 + (t1558 * t1954 + t1560 * t1944 + t1561 * t1942 + t1562 * t1946) * pkin(2)) * t1820, 0, -t2104, -t2099, 0; (t1822 + t1823 + t1824 + t1825) * t1781 + t2034 + t2035 + t2036 + t2043, 0, (t1855 + t1857 + t1859 + t1861) * t1783, (t1854 + t1856 + t1858 + t1860) * t1783, 0, 0, 0, 0, 0, (-t1974 - t1977 - t1980 - t1985) * t2066 + (-t1998 - t2002 - t2006 - t2010) * t1783 + t2199 + t2200 + t2201 + t2202, (-t1973 - t1976 - t1979 - t1984) * t2066 + (-t1498 * t2063 - t1499 * t2058 - t1500 * t2054 - t1501 * t2050) * t1783 + t2195 + t2196 + t2197 + t2198, 0, 0, 0, 0; t1502 * t1554 + t1503 * t1555 + t1504 * t1556 + t1505 * t1557, t2015 + t2017 + t2019 + t2021, (t1554 * t1861 + t1555 * t1859 + t1556 * t1857 + t1557 * t1855 + t1793 * t1997 + t1801 * t1996 + t1803 * t1995 + t1805 * t1994) * t1783, (t1554 * t1860 + t1555 * t1858 + t1556 * t1856 + t1557 * t1854 - t1791 * t1997 - t1795 * t1996 - t1797 * t1995 - t1799 * t1994) * t1783, t1773 * t2021 + t1775 * t2019 + t1776 * t2017 + t1777 * t2015 + ((-t1561 * t1581 + t1865) * t1934 + (-t1560 * t1580 + t1867) * t1936 + (-t1562 * t1579 + t1863) * t1938 + (-t1558 * t1578 + t1869) * t1940) * t1820, t1921 * t2143 + t1920 * t2139 + t1919 * t2138 + t1918 * t2137 + ((-t1553 * t1581 + t1758 * t1865) * t1646 + (-t1552 * t1580 + t1757 * t1867) * t1644 + (-t1551 * t1579 + t1756 * t1863) * t1642 + (-t1550 * t1578 + t1755 * t1869) * t1638) * t1820, t1790 * t2013 + t1794 * t2009 + t1796 * t2005 + t1798 * t2001 + (t1792 * t1897 + t1800 * t1896 + t1802 * t1895 + t1804 * t1894) * t1821 + (t1494 * t1965 + t1495 * t1963 + t1496 * t1961 + t1497 * t1959) * t1820, t1792 * t2013 + t1800 * t2009 + t1802 * t2005 + t1804 * t2001 + (-t1790 * t1897 - t1794 * t1896 - t1796 * t1895 - t1798 * t1894) * t1821 + (t1494 * t1964 + t1495 * t1962 + t1496 * t1960 + t1497 * t1958) * t1820, (t1578 * t2168 + t1579 * t2165 + t1580 * t2162 + t1581 * t2159) * t1820, t1456 * t2143 + t1466 * t2138 + t1468 * t2139 + t1470 * t2137 + t2037 * t1557 + t2038 * t1556 + t2039 * t1555 + t2044 * t1554 + ((t1474 * t1581 + t1539 * t1910) * t1646 + (t1473 * t1580 + t1538 * t1911) * t1644 + (t1472 * t1579 + t1537 * t1912) * t1642 + (t1464 * t1578 + t1534 * t1916) * t1638 + (t1558 * t1965 + t1560 * t1961 + t1561 * t1959 + t1562 * t1963) * pkin(2)) * t1820, t1457 * t2143 + t1467 * t2138 + t1469 * t2139 + t1471 * t2137 + t2041 * t1557 + t2042 * t1556 + t2040 * t1555 + t2045 * t1554 + ((t1477 * t1581 + t1542 * t1910) * t1646 + (t1476 * t1580 + t1541 * t1911) * t1644 + (t1475 * t1579 + t1540 * t1912) * t1642 + (t1465 * t1578 + t1535 * t1916) * t1638 + (t1558 * t1964 + t1560 * t1960 + t1561 * t1958 + t1562 * t1962) * pkin(2)) * t1820, 0, 0, 0, 0;];
tau_reg  = t1;
