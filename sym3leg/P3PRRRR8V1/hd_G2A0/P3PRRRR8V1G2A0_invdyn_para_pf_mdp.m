% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PRRRR8V1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [12x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRRR8V1G2A0_convert_par2_MPV_fixb.m

% Output:
% tauX [3x1]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 17:03
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRRR8V1G2A0_invdyn_para_pf_mdp(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_mdp: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_mdp: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:03:28
% EndTime: 2020-08-06 17:03:37
% DurationCPUTime: 8.76s
% Computational Cost: add. (33969->435), mult. (88302->875), div. (6471->14), fcn. (86574->22), ass. (0->364)
t1755 = cos(qJ(2,1));
t1749 = sin(qJ(2,1));
t1754 = cos(qJ(3,1));
t1873 = t1749 * t1754;
t1704 = pkin(2) * t1873 - pkin(5) * t1755;
t1735 = sin(pkin(3));
t1737 = cos(pkin(3));
t1748 = sin(qJ(3,1));
t1944 = pkin(2) * t1748;
t1831 = t1704 * t1735 + t1737 * t1944;
t1677 = 0.1e1 / t1831 ^ 2;
t1955 = 0.1e1 / t1831;
t1753 = cos(qJ(2,2));
t1747 = sin(qJ(2,2));
t1752 = cos(qJ(3,2));
t1876 = t1747 * t1752;
t1703 = pkin(2) * t1876 - pkin(5) * t1753;
t1746 = sin(qJ(3,2));
t1945 = pkin(2) * t1746;
t1832 = t1703 * t1735 + t1737 * t1945;
t1675 = 0.1e1 / t1832 ^ 2;
t1954 = 0.1e1 / t1832;
t1751 = cos(qJ(2,3));
t1745 = sin(qJ(2,3));
t1750 = cos(qJ(3,3));
t1879 = t1745 * t1750;
t1702 = pkin(2) * t1879 - pkin(5) * t1751;
t1744 = sin(qJ(3,3));
t1946 = pkin(2) * t1744;
t1833 = t1702 * t1735 + t1737 * t1946;
t1673 = 0.1e1 / t1833 ^ 2;
t1953 = 0.1e1 / t1833;
t1738 = legFrame(3,2);
t1722 = sin(t1738);
t1725 = cos(t1738);
t1757 = xDP(2);
t1758 = xDP(1);
t1697 = -t1722 * t1757 + t1725 * t1758;
t1736 = cos(pkin(6));
t1734 = sin(pkin(6));
t1756 = xDP(3);
t1894 = t1734 * t1756;
t1663 = t1697 * t1736 - t1894;
t1713 = t1756 * t1736;
t1666 = t1697 * t1734 + t1713;
t1886 = t1737 * t1751;
t1889 = t1737 * t1745;
t1943 = pkin(2) * t1750;
t1629 = -(t1663 * t1745 + t1666 * t1886) * t1943 + (t1663 * t1751 - t1666 * t1889) * pkin(5);
t1883 = t1737 * t1756;
t1893 = t1735 * t1750;
t1635 = ((-t1697 * t1889 - t1751 * t1756) * t1734 + (t1697 * t1751 - t1745 * t1883) * t1736) * t1744 - t1666 * t1893;
t1761 = 0.1e1 / pkin(2);
t1903 = t1953 * t1761;
t1848 = t1953 * t1903;
t1821 = t1635 * t1848;
t1792 = t1629 * t1821;
t1728 = 0.1e1 / t1750;
t1917 = t1635 * t1953;
t1872 = pkin(5) * t1917;
t1830 = t1744 * t1872;
t1926 = t1629 * t1953;
t1617 = (t1830 - t1926) * t1728;
t1763 = t1750 ^ 2;
t1729 = 0.1e1 / t1763;
t1741 = xDDP(3);
t1742 = xDDP(2);
t1743 = xDDP(1);
t1844 = t1728 * t1903;
t1791 = t1629 * t1735 * t1844;
t1881 = t1744 * t1745;
t1771 = t1737 * t1881 + t1893;
t1880 = t1744 * t1751;
t1656 = t1771 * t1734 - t1736 * t1880;
t1905 = t1953 * t1728;
t1846 = t1725 * t1905;
t1811 = t1656 * t1846;
t1847 = t1722 * t1905;
t1812 = t1656 * t1847;
t1659 = -t1734 * t1880 - t1771 * t1736;
t1852 = t1659 * t1905;
t1858 = t1953 * t1917;
t1882 = t1737 * t1761;
t1904 = t1953 * t1744;
t1925 = t1629 * t1673;
t1575 = t1741 * t1852 + t1742 * t1812 - t1743 * t1811 + (-(t1737 * t1617 + (pkin(2) * (t1735 * t1751 * t1917 + t1882 * t1926) * t1763 - t1735 * (t1629 * t1904 - t1872) * t1879) * t1728) * t1858 + (-t1751 * t1791 + (t1735 * t1881 + (t1728 - t1750) * t1737) * t1917) * t1925) * t1729;
t1937 = t1575 * t1744;
t1949 = t1729 - 0.2e1;
t1563 = t1750 * t1937 - t1949 * t1792;
t1952 = 0.2e1 * t1563;
t1739 = legFrame(2,2);
t1723 = sin(t1739);
t1726 = cos(t1739);
t1698 = -t1723 * t1757 + t1726 * t1758;
t1664 = t1698 * t1736 - t1894;
t1667 = t1698 * t1734 + t1713;
t1885 = t1737 * t1753;
t1888 = t1737 * t1747;
t1942 = pkin(2) * t1752;
t1630 = -(t1664 * t1747 + t1667 * t1885) * t1942 + (t1664 * t1753 - t1667 * t1888) * pkin(5);
t1892 = t1735 * t1752;
t1636 = ((-t1698 * t1888 - t1753 * t1756) * t1734 + (t1698 * t1753 - t1747 * t1883) * t1736) * t1746 - t1667 * t1892;
t1900 = t1954 * t1761;
t1843 = t1954 * t1900;
t1820 = t1636 * t1843;
t1790 = t1630 * t1820;
t1730 = 0.1e1 / t1752;
t1916 = t1636 * t1954;
t1871 = pkin(5) * t1916;
t1829 = t1746 * t1871;
t1924 = t1630 * t1954;
t1618 = (t1829 - t1924) * t1730;
t1764 = t1752 ^ 2;
t1731 = 0.1e1 / t1764;
t1839 = t1730 * t1900;
t1789 = t1630 * t1735 * t1839;
t1878 = t1746 * t1747;
t1770 = t1737 * t1878 + t1892;
t1877 = t1746 * t1753;
t1657 = t1770 * t1734 - t1736 * t1877;
t1902 = t1954 * t1730;
t1841 = t1726 * t1902;
t1809 = t1657 * t1841;
t1842 = t1723 * t1902;
t1810 = t1657 * t1842;
t1660 = -t1734 * t1877 - t1770 * t1736;
t1851 = t1660 * t1902;
t1857 = t1954 * t1916;
t1901 = t1954 * t1746;
t1923 = t1630 * t1675;
t1576 = t1741 * t1851 + t1742 * t1810 - t1743 * t1809 + (-(t1737 * t1618 + (pkin(2) * (t1735 * t1753 * t1916 + t1882 * t1924) * t1764 - t1735 * (t1630 * t1901 - t1871) * t1876) * t1730) * t1857 + (-t1753 * t1789 + (t1735 * t1878 + (t1730 - t1752) * t1737) * t1916) * t1923) * t1731;
t1934 = t1576 * t1746;
t1948 = t1731 - 0.2e1;
t1564 = t1752 * t1934 - t1948 * t1790;
t1951 = 0.2e1 * t1564;
t1740 = legFrame(1,2);
t1724 = sin(t1740);
t1727 = cos(t1740);
t1699 = -t1724 * t1757 + t1727 * t1758;
t1665 = t1699 * t1736 - t1894;
t1668 = t1699 * t1734 + t1713;
t1884 = t1737 * t1755;
t1887 = t1737 * t1749;
t1941 = pkin(2) * t1754;
t1631 = -(t1665 * t1749 + t1668 * t1884) * t1941 + (t1665 * t1755 - t1668 * t1887) * pkin(5);
t1891 = t1735 * t1754;
t1637 = ((-t1699 * t1887 - t1755 * t1756) * t1734 + (t1699 * t1755 - t1749 * t1883) * t1736) * t1748 - t1668 * t1891;
t1897 = t1955 * t1761;
t1838 = t1955 * t1897;
t1819 = t1637 * t1838;
t1788 = t1631 * t1819;
t1732 = 0.1e1 / t1754;
t1915 = t1637 * t1955;
t1870 = pkin(5) * t1915;
t1828 = t1748 * t1870;
t1922 = t1631 * t1955;
t1619 = (t1828 - t1922) * t1732;
t1875 = t1748 * t1749;
t1769 = t1737 * t1875 + t1891;
t1874 = t1748 * t1755;
t1658 = t1769 * t1734 - t1736 * t1874;
t1765 = t1754 ^ 2;
t1733 = 0.1e1 / t1765;
t1834 = t1732 * t1897;
t1787 = t1631 * t1735 * t1834;
t1899 = t1955 * t1732;
t1836 = t1727 * t1899;
t1808 = t1658 * t1836;
t1895 = t1724 * t1732;
t1837 = t1955 * t1895;
t1661 = -t1734 * t1874 - t1769 * t1736;
t1850 = t1661 * t1899;
t1856 = t1955 * t1915;
t1898 = t1955 * t1748;
t1921 = t1631 * t1677;
t1577 = t1658 * t1742 * t1837 + t1741 * t1850 - t1743 * t1808 + (-(t1737 * t1619 + (pkin(2) * (t1735 * t1755 * t1915 + t1882 * t1922) * t1765 - t1735 * (t1631 * t1898 - t1870) * t1873) * t1732) * t1856 + (-t1755 * t1787 + (t1735 * t1875 + (t1732 - t1754) * t1737) * t1915) * t1921) * t1733;
t1931 = t1577 * t1748;
t1947 = t1733 - 0.2e1;
t1565 = t1754 * t1931 - t1947 * t1788;
t1950 = 0.2e1 * t1565;
t1940 = g(3) * t1734;
t1939 = MDP(9) * t1761;
t1938 = t1575 * t1953;
t1936 = t1575 * t1751;
t1935 = t1576 * t1954;
t1933 = t1576 * t1753;
t1932 = t1577 * t1955;
t1930 = t1577 * t1755;
t1721 = t1736 * g(3);
t1801 = g(1) * t1725 - g(2) * t1722;
t1693 = t1801 * t1734 + t1721;
t1759 = pkin(5) ^ 2;
t1760 = pkin(2) ^ 2;
t1845 = t1728 * t1904;
t1614 = -pkin(5) * t1629 * t1845 + (t1728 * t1759 + t1750 * t1760) * t1917;
t1705 = pkin(5) * t1745 + t1751 * t1943;
t1777 = -t1702 * t1737 + t1735 * t1946;
t1644 = -t1705 * t1734 + t1777 * t1736;
t1638 = t1644 * t1722 + t1725 * t1833;
t1639 = -t1644 * t1725 + t1722 * t1833;
t1647 = t1705 * t1736 + t1777 * t1734;
t1605 = (t1638 * t1742 + t1639 * t1743 + t1647 * t1741) * t1953 + (t1614 * t1858 - t1617 * t1925) * t1728;
t1890 = t1736 * t1737;
t1700 = -t1735 * g(1) + g(2) * t1890;
t1701 = g(1) * t1890 + g(2) * t1735;
t1709 = t1737 * t1940;
t1768 = t1605 * t1735 + t1700 * t1722 - t1701 * t1725 + t1709;
t1593 = t1745 * t1693 + t1768 * t1751;
t1929 = t1593 * t1953;
t1800 = g(1) * t1726 - g(2) * t1723;
t1694 = t1800 * t1734 + t1721;
t1840 = t1730 * t1901;
t1615 = -pkin(5) * t1630 * t1840 + (t1730 * t1759 + t1752 * t1760) * t1916;
t1706 = pkin(5) * t1747 + t1753 * t1942;
t1776 = -t1703 * t1737 + t1735 * t1945;
t1645 = -t1706 * t1734 + t1776 * t1736;
t1640 = t1645 * t1723 + t1726 * t1832;
t1641 = -t1645 * t1726 + t1723 * t1832;
t1648 = t1706 * t1736 + t1776 * t1734;
t1606 = (t1640 * t1742 + t1641 * t1743 + t1648 * t1741) * t1954 + (t1615 * t1857 - t1618 * t1923) * t1730;
t1767 = t1606 * t1735 + t1700 * t1723 - t1701 * t1726 + t1709;
t1594 = t1747 * t1694 + t1767 * t1753;
t1928 = t1594 * t1954;
t1799 = g(1) * t1727 - g(2) * t1724;
t1695 = t1799 * t1734 + t1721;
t1835 = t1732 * t1898;
t1616 = -pkin(5) * t1631 * t1835 + (t1732 * t1759 + t1754 * t1760) * t1915;
t1707 = pkin(5) * t1749 + t1755 * t1941;
t1775 = -t1704 * t1737 + t1735 * t1944;
t1646 = -t1707 * t1734 + t1775 * t1736;
t1642 = t1646 * t1724 + t1727 * t1831;
t1643 = -t1646 * t1727 + t1724 * t1831;
t1649 = t1707 * t1736 + t1775 * t1734;
t1607 = (t1642 * t1742 + t1643 * t1743 + t1649 * t1741) * t1955 + (t1616 * t1856 - t1619 * t1921) * t1732;
t1766 = t1607 * t1735 + t1700 * t1724 - t1701 * t1727 + t1709;
t1595 = t1749 * t1695 + t1766 * t1755;
t1927 = t1595 * t1955;
t1920 = t1635 ^ 2 * t1673;
t1919 = t1636 ^ 2 * t1675;
t1918 = t1637 ^ 2 * t1677;
t1914 = t1638 * t1953;
t1913 = t1639 * t1953;
t1912 = t1640 * t1954;
t1911 = t1641 * t1954;
t1910 = t1642 * t1955;
t1909 = t1643 * t1955;
t1908 = t1647 * t1953;
t1907 = t1648 * t1954;
t1906 = t1649 * t1955;
t1651 = -(-t1734 * t1747 + t1736 * t1885) * t1942 - pkin(5) * (t1734 * t1753 + t1736 * t1888);
t1654 = (t1734 * t1885 + t1736 * t1747) * t1942 + (t1734 * t1888 - t1736 * t1753) * pkin(5);
t1588 = (-t1737 * t1615 * t1820 - (-t1746 * t1703 * t1789 + t1737 * (-t1730 * t1829 + t1752 * t1924)) * t1630 * t1843) * t1731 + (t1651 * t1741 + (t1723 * t1742 - t1726 * t1743) * t1654) * t1839;
t1779 = 0.2e1 * t1790;
t1579 = t1731 * t1753 * t1779 + t1588 * t1747;
t1762 = 0.1e1 / pkin(2) ^ 2;
t1863 = t1630 ^ 2 * t1675 * t1762;
t1585 = t1588 * t1746 + t1730 * t1863;
t1782 = -(t1863 + t1919) * t1731 * t1747 + t1933;
t1896 = t1954 * ((-t1752 * t1579 - t1782 * t1746) * t1735 - t1737 * t1585);
t1653 = (t1734 * t1886 + t1736 * t1745) * t1943 + (t1734 * t1889 - t1736 * t1751) * pkin(5);
t1869 = t1653 * t1938;
t1868 = t1654 * t1935;
t1696 = t1734 * t1887 - t1736 * t1755;
t1655 = (t1734 * t1884 + t1736 * t1749) * t1941 + t1696 * pkin(5);
t1867 = t1655 * t1932;
t1866 = t1656 * t1929;
t1865 = t1657 * t1928;
t1864 = t1629 ^ 2 * t1673 * t1762;
t1862 = t1631 ^ 2 * t1677 * t1762;
t1861 = t1729 * t1920;
t1860 = t1731 * t1919;
t1859 = t1733 * t1918;
t1650 = -(-t1734 * t1745 + t1736 * t1886) * t1943 - pkin(5) * (t1734 * t1751 + t1736 * t1889);
t1855 = t1650 * t1905;
t1854 = t1651 * t1902;
t1652 = -(-t1734 * t1749 + t1736 * t1884) * t1941 - pkin(5) * (t1734 * t1755 + t1736 * t1887);
t1853 = t1652 * t1899;
t1662 = t1696 * t1748 + t1734 * t1891;
t1849 = t1662 * t1895;
t1827 = t1575 * t1845;
t1826 = t1576 * t1840;
t1825 = t1577 * t1835;
t1824 = t1593 * t1845;
t1823 = t1594 * t1840;
t1822 = t1595 * t1850;
t1818 = t1653 * t1847;
t1817 = t1653 * t1846;
t1816 = t1654 * t1842;
t1815 = t1654 * t1841;
t1814 = t1655 * t1837;
t1813 = t1655 * t1836;
t1807 = t1955 * t1849;
t1806 = t1722 * t1866;
t1805 = t1723 * t1865;
t1804 = t1725 * t1866;
t1803 = t1726 * t1865;
t1802 = t1727 * t1658 * t1927;
t1798 = t1653 * t1827;
t1797 = t1654 * t1826;
t1796 = t1655 * t1825;
t1795 = t1656 * t1824;
t1794 = t1657 * t1823;
t1793 = t1732 * t1802;
t1786 = t1861 * t1904;
t1785 = t1860 * t1901;
t1784 = t1859 * t1898;
t1783 = -(t1864 + t1920) * t1729 * t1745 + t1936;
t1781 = -(t1862 + t1918) * t1733 * t1749 + t1930;
t1780 = 0.2e1 * t1792;
t1778 = 0.2e1 * t1788;
t1774 = t1653 * t1786;
t1773 = t1654 * t1785;
t1772 = t1655 * t1784;
t1692 = t1799 * t1736 - t1940;
t1691 = t1800 * t1736 - t1940;
t1690 = t1801 * t1736 - t1940;
t1671 = t1755 * t1695;
t1670 = t1753 * t1694;
t1669 = t1751 * t1693;
t1625 = t1947 * t1918;
t1624 = t1948 * t1919;
t1623 = t1949 * t1920;
t1604 = -t1724 * g(1) - t1727 * g(2) + t1607;
t1603 = -t1723 * g(1) - t1726 * g(2) + t1606;
t1602 = -t1722 * g(1) - t1725 * g(2) + t1605;
t1601 = -t1604 * t1737 - t1692 * t1735;
t1600 = -t1603 * t1737 - t1691 * t1735;
t1599 = -t1602 * t1737 - t1690 * t1735;
t1598 = t1671 + (-t1604 * t1735 + t1692 * t1737) * t1749;
t1597 = t1670 + (-t1603 * t1735 + t1691 * t1737) * t1747;
t1596 = t1669 + (-t1602 * t1735 + t1690 * t1737) * t1745;
t1592 = -t1766 * t1749 + t1671;
t1591 = -t1767 * t1747 + t1670;
t1590 = -t1768 * t1745 + t1669;
t1589 = (-t1737 * t1616 * t1819 - (-t1748 * t1704 * t1787 + t1737 * (-t1732 * t1828 + t1754 * t1922)) * t1631 * t1838) * t1733 + (t1652 * t1741 + (t1724 * t1742 - t1727 * t1743) * t1655) * t1834;
t1587 = (-t1737 * t1614 * t1821 - (-t1744 * t1702 * t1791 + t1737 * (-t1728 * t1830 + t1750 * t1926)) * t1629 * t1848) * t1729 + (t1650 * t1741 + (t1722 * t1742 - t1725 * t1743) * t1653) * t1844;
t1586 = t1589 * t1748 + t1732 * t1862;
t1584 = t1587 * t1744 + t1728 * t1864;
t1583 = -t1733 * t1748 * t1862 + t1589 * t1754;
t1582 = -t1731 * t1746 * t1863 + t1588 * t1752;
t1581 = -t1729 * t1744 * t1864 + t1587 * t1750;
t1580 = t1733 * t1755 * t1778 + t1589 * t1749;
t1578 = t1729 * t1751 * t1780 + t1587 * t1745;
t1574 = t1576 * t1747 + t1753 * t1860;
t1573 = t1747 * t1860 - t1933;
t1572 = t1577 * t1749 + t1755 * t1859;
t1571 = t1575 * t1745 + t1751 * t1861;
t1570 = t1749 * t1859 - t1930;
t1569 = t1745 * t1861 - t1936;
t1568 = (t1732 * t1778 + t1931) * t1748;
t1567 = (t1730 * t1779 + t1934) * t1746;
t1566 = (t1728 * t1780 + t1937) * t1744;
t1562 = t1598 * t1754 + t1601 * t1748;
t1561 = t1598 * t1748 - t1601 * t1754;
t1560 = t1597 * t1752 + t1600 * t1746;
t1559 = t1597 * t1746 - t1600 * t1752;
t1558 = t1596 * t1750 + t1599 * t1744;
t1557 = t1596 * t1744 - t1599 * t1750;
t1556 = (-t1754 * t1580 - t1781 * t1748) * t1735 - t1737 * t1586;
t1554 = (-t1750 * t1578 - t1783 * t1744) * t1735 - t1737 * t1584;
t1553 = (-t1748 * t1580 + t1781 * t1754) * t1735 + t1737 * t1583;
t1552 = (-t1746 * t1579 + t1782 * t1752) * t1735 + t1737 * t1582;
t1551 = (-t1744 * t1578 + t1783 * t1750) * t1735 + t1737 * t1581;
t1 = [(t1602 * t1913 + t1603 * t1911 + t1604 * t1909) * MDP(1) + (-t1575 * t1811 - t1576 * t1809 - t1577 * t1808) * MDP(2) + (-t1728 * t1804 - t1730 * t1803 - t1793 + (-t1569 * t1913 - t1570 * t1909 - t1573 * t1911) * t1735) * MDP(3) + (-t1590 * t1811 - t1591 * t1809 - t1592 * t1808 + (-t1571 * t1913 - t1572 * t1909 - t1574 * t1911) * t1735) * MDP(4) + (-t1566 * t1811 - t1567 * t1809 - t1568 * t1808 + (t1725 * t1774 + t1726 * t1773 + t1727 * t1772) * t1761) * MDP(5) + ((-t1623 * t1817 - t1624 * t1815 - t1625 * t1813) * t1761 - 0.2e1 * t1563 * t1811 - 0.2e1 * t1564 * t1809 - 0.2e1 * t1565 * t1808) * MDP(6) + (-t1584 * t1811 - t1585 * t1809 - t1586 * t1808 + (-t1725 * t1798 - t1726 * t1797 - t1727 * t1796) * t1761) * MDP(7) + (-t1581 * t1811 - t1582 * t1809 - t1583 * t1808 + (-t1725 * t1869 - t1726 * t1868 - t1727 * t1867) * t1761) * MDP(8) + (-t1587 * t1817 - t1588 * t1815 - t1589 * t1813) * t1939 + (-t1804 - t1803 - t1802 + t1551 * t1913 + t1552 * t1911 + t1553 * t1909 + (-t1557 * t1817 - t1559 * t1815 - t1561 * t1813) * t1761) * MDP(10) + (t1725 * t1795 + t1726 * t1794 + t1748 * t1793 + t1554 * t1913 + t1641 * t1896 + t1556 * t1909 + (-t1558 * t1817 - t1560 * t1815 - t1562 * t1813) * t1761) * MDP(11) + (t1743 - g(1)) * MDP(12); (t1602 * t1914 + t1603 * t1912 + t1604 * t1910) * MDP(1) + (t1575 * t1812 + t1576 * t1810 + t1577 * t1807) * MDP(2) + (t1728 * t1806 + t1730 * t1805 + t1595 * t1807 + (-t1569 * t1914 - t1570 * t1910 - t1573 * t1912) * t1735) * MDP(3) + (t1590 * t1812 + t1591 * t1810 + t1592 * t1807 + (-t1571 * t1914 - t1572 * t1910 - t1574 * t1912) * t1735) * MDP(4) + (t1566 * t1812 + t1567 * t1810 + t1568 * t1807 + (-t1722 * t1774 - t1723 * t1773 - t1724 * t1772) * t1761) * MDP(5) + ((t1623 * t1818 + t1624 * t1816 + t1625 * t1814) * t1761 + t1812 * t1952 + t1810 * t1951 + t1807 * t1950) * MDP(6) + (t1584 * t1812 + t1585 * t1810 + t1586 * t1807 + (t1722 * t1798 + t1723 * t1797 + t1724 * t1796) * t1761) * MDP(7) + (t1581 * t1812 + t1582 * t1810 + t1583 * t1807 + (t1722 * t1869 + t1723 * t1868 + t1724 * t1867) * t1761) * MDP(8) + (t1587 * t1818 + t1588 * t1816 + t1589 * t1814) * t1939 + (t1806 + t1805 + t1551 * t1914 + t1552 * t1912 + (t1595 * t1662 * t1724 + t1553 * t1642) * t1955 + (t1557 * t1818 + t1559 * t1816 + t1561 * t1814) * t1761) * MDP(10) + (-t1722 * t1795 - t1723 * t1794 + t1554 * t1914 + t1640 * t1896 + (-t1595 * t1748 * t1849 + t1556 * t1642) * t1955 + (t1558 * t1818 + t1560 * t1816 + t1562 * t1814) * t1761) * MDP(11) + (t1742 - g(2)) * MDP(12); (t1602 * t1908 + t1603 * t1907 + t1604 * t1906) * MDP(1) + (t1575 * t1852 + t1576 * t1851 + t1577 * t1850) * MDP(2) + (t1593 * t1852 + t1594 * t1851 + t1822 + (-t1569 * t1908 - t1570 * t1906 - t1573 * t1907) * t1735) * MDP(3) + (t1590 * t1852 + t1591 * t1851 + t1592 * t1850 + (-t1571 * t1908 - t1572 * t1906 - t1574 * t1907) * t1735) * MDP(4) + (t1566 * t1852 + t1567 * t1851 + t1568 * t1850 + (-t1650 * t1786 - t1651 * t1785 - t1652 * t1784) * t1761) * MDP(5) + ((t1623 * t1855 + t1624 * t1854 + t1625 * t1853) * t1761 + t1852 * t1952 + t1851 * t1951 + t1850 * t1950) * MDP(6) + (t1584 * t1852 + t1585 * t1851 + t1586 * t1850 + (t1650 * t1827 + t1651 * t1826 + t1652 * t1825) * t1761) * MDP(7) + (t1581 * t1852 + t1582 * t1851 + t1583 * t1850 + (t1650 * t1938 + t1651 * t1935 + t1652 * t1932) * t1761) * MDP(8) + (t1587 * t1855 + t1588 * t1854 + t1589 * t1853) * t1939 + (t1551 * t1908 + t1552 * t1907 + t1553 * t1906 + t1659 * t1929 + t1660 * t1928 + t1661 * t1927 + (t1557 * t1855 + t1559 * t1854 + t1561 * t1853) * t1761) * MDP(10) + (-t1659 * t1824 - t1660 * t1823 - t1748 * t1822 + t1554 * t1908 + t1648 * t1896 + t1556 * t1906 + (t1558 * t1855 + t1560 * t1854 + t1562 * t1853) * t1761) * MDP(11) + (t1741 - g(3)) * MDP(12);];
tauX  = t1;
