% Calculate minimal parameter regressor of inverse dynamics forces for
% P4PRRRR1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% xDP [4x1]
%   Generalized platform velocities
% xDDP [4x1]
%   Generalized platform accelerations
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [15x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P4PRRRR1G2A0_convert_par2_MPV_fixb.m

% Output:
% tauX [4x1]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:00
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P4PRRRR1G2A0_invdyn_para_pf_mdp(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(2,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_mdp: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_mdp: xDP has to be [4x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [4 1]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_mdp: xDDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_mdp: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_mdp: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_mdp: Koppelpunkt has to be [4x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P4PRRRR1G2A0_invdyn_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 10:59:24
% EndTime: 2020-08-07 10:59:31
% DurationCPUTime: 7.34s
% Computational Cost: add. (13187->473), mult. (29289->881), div. (6272->22), fcn. (22400->26), ass. (0->361)
t1872 = xP(4);
t1813 = sin(t1872);
t1814 = cos(t1872);
t1873 = koppelP(4,2);
t1877 = koppelP(4,1);
t1785 = t1813 * t1877 + t1814 * t1873;
t1789 = -t1813 * t1873 + t1814 * t1877;
t1848 = legFrame(4,2);
t1805 = sin(t1848);
t1809 = cos(t1848);
t1764 = t1785 * t1809 + t1805 * t1789;
t1868 = xDP(4);
t1870 = xDP(2);
t1871 = xDP(1);
t1759 = t1764 * t1868 + t1805 * t1870 - t1871 * t1809;
t1755 = t1759 ^ 2;
t1843 = t1868 ^ 2;
t1852 = xDDP(4);
t1854 = xDDP(2);
t1767 = -t1785 * t1843 + t1789 * t1852 + t1854;
t1855 = xDDP(1);
t1771 = -t1785 * t1852 - t1789 * t1843 + t1855;
t1845 = sin(qJ(2,4));
t1846 = cos(qJ(3,4));
t2000 = t1845 * t1846;
t1844 = sin(qJ(3,4));
t2020 = t1809 * t1844;
t1775 = t1805 * t2000 - t2020;
t1776 = t1805 * t1844 + t1809 * t2000;
t1816 = 0.1e1 / t1845;
t1818 = 0.1e1 / t1846;
t1847 = cos(qJ(2,4));
t1853 = xDDP(3);
t1888 = t1846 ^ 2;
t1819 = 0.1e1 / t1888;
t1820 = t1818 * t1819;
t1881 = 0.1e1 / pkin(2);
t2011 = t1820 * t1881;
t1869 = xDP(3);
t2001 = t1844 * t1847;
t1751 = -t1759 * t2001 - t1846 * t1869;
t1747 = t1751 ^ 2;
t1817 = 0.1e1 / t1845 ^ 2;
t2064 = t1747 * t1817;
t1719 = -t1805 * g(1) - t1809 * g(2) + (t1847 * t1853 + (t1767 * t1776 + t1771 * t1775) * t1818 + (t1755 + t2064) * t2011) * t1816;
t1707 = g(3) * t1845 + t1719 * t1847;
t2013 = t1816 * t1818;
t2068 = t1707 * t2013;
t1874 = koppelP(3,2);
t1878 = koppelP(3,1);
t1786 = t1813 * t1878 + t1814 * t1874;
t1790 = -t1813 * t1874 + t1814 * t1878;
t1849 = legFrame(3,2);
t1806 = sin(t1849);
t1810 = cos(t1849);
t1765 = t1786 * t1810 + t1806 * t1790;
t1760 = t1765 * t1868 + t1806 * t1870 - t1871 * t1810;
t1756 = t1760 ^ 2;
t1768 = -t1786 * t1843 + t1790 * t1852 + t1854;
t1772 = -t1786 * t1852 - t1790 * t1843 + t1855;
t1857 = sin(qJ(2,3));
t1862 = cos(qJ(3,3));
t1998 = t1857 * t1862;
t1856 = sin(qJ(3,3));
t2018 = t1810 * t1856;
t1777 = t1806 * t1998 - t2018;
t1780 = t1806 * t1856 + t1810 * t1998;
t1823 = 0.1e1 / t1857;
t1831 = 0.1e1 / t1862;
t1863 = cos(qJ(2,3));
t1897 = t1862 ^ 2;
t1832 = 0.1e1 / t1897;
t1833 = t1831 * t1832;
t2004 = t1833 * t1881;
t1999 = t1856 * t1863;
t1752 = -t1760 * t1999 - t1862 * t1869;
t1748 = t1752 ^ 2;
t1824 = 0.1e1 / t1857 ^ 2;
t2063 = t1748 * t1824;
t1720 = -t1806 * g(1) - t1810 * g(2) + (t1863 * t1853 + (t1768 * t1780 + t1772 * t1777) * t1831 + (t1756 + t2063) * t2004) * t1823;
t1710 = g(3) * t1857 + t1720 * t1863;
t2010 = t1823 * t1831;
t2067 = t1710 * t2010;
t1875 = koppelP(2,2);
t1879 = koppelP(2,1);
t1787 = t1813 * t1879 + t1814 * t1875;
t1791 = -t1813 * t1875 + t1814 * t1879;
t1850 = legFrame(2,2);
t1807 = sin(t1850);
t1811 = cos(t1850);
t1766 = t1787 * t1811 + t1807 * t1791;
t1761 = t1766 * t1868 + t1807 * t1870 - t1871 * t1811;
t1757 = t1761 ^ 2;
t1769 = -t1787 * t1843 + t1791 * t1852 + t1854;
t1773 = -t1787 * t1852 - t1791 * t1843 + t1855;
t1859 = sin(qJ(2,2));
t1864 = cos(qJ(3,2));
t1996 = t1859 * t1864;
t1858 = sin(qJ(3,2));
t2016 = t1811 * t1858;
t1778 = t1807 * t1996 - t2016;
t1781 = t1807 * t1858 + t1811 * t1996;
t1826 = 0.1e1 / t1859;
t1835 = 0.1e1 / t1864;
t1865 = cos(qJ(2,2));
t1900 = t1864 ^ 2;
t1836 = 0.1e1 / t1900;
t1837 = t1835 * t1836;
t2003 = t1837 * t1881;
t1997 = t1858 * t1865;
t1753 = -t1761 * t1997 - t1864 * t1869;
t1749 = t1753 ^ 2;
t1827 = 0.1e1 / t1859 ^ 2;
t2062 = t1749 * t1827;
t1721 = -t1807 * g(1) - t1811 * g(2) + (t1865 * t1853 + (t1769 * t1781 + t1773 * t1778) * t1835 + (t1757 + t2062) * t2003) * t1826;
t1711 = g(3) * t1859 + t1721 * t1865;
t2008 = t1826 * t1835;
t2066 = t1711 * t2008;
t1876 = koppelP(1,2);
t1880 = koppelP(1,1);
t1788 = t1813 * t1880 + t1814 * t1876;
t1792 = -t1813 * t1876 + t1814 * t1880;
t1851 = legFrame(1,2);
t1808 = sin(t1851);
t1812 = cos(t1851);
t1763 = t1788 * t1812 + t1808 * t1792;
t1762 = t1763 * t1868 + t1808 * t1870 - t1871 * t1812;
t1758 = t1762 ^ 2;
t1770 = -t1788 * t1843 + t1792 * t1852 + t1854;
t1774 = -t1788 * t1852 - t1792 * t1843 + t1855;
t1861 = sin(qJ(2,1));
t1866 = cos(qJ(3,1));
t1994 = t1861 * t1866;
t1860 = sin(qJ(3,1));
t2014 = t1812 * t1860;
t1779 = t1808 * t1994 - t2014;
t1782 = t1808 * t1860 + t1812 * t1994;
t1829 = 0.1e1 / t1861;
t1839 = 0.1e1 / t1866;
t1867 = cos(qJ(2,1));
t1903 = t1866 ^ 2;
t1840 = 0.1e1 / t1903;
t1841 = t1839 * t1840;
t2002 = t1841 * t1881;
t1995 = t1860 * t1867;
t1754 = -t1762 * t1995 - t1866 * t1869;
t1750 = t1754 ^ 2;
t1830 = 0.1e1 / t1861 ^ 2;
t2061 = t1750 * t1830;
t1722 = -t1808 * g(1) - t1812 * g(2) + (t1867 * t1853 + (t1770 * t1782 + t1774 * t1779) * t1839 + (t1758 + t2061) * t2002) * t1829;
t1712 = g(3) * t1861 + t1722 * t1867;
t2006 = t1829 * t1839;
t2065 = t1712 * t2006;
t1882 = 0.1e1 / pkin(2) ^ 2;
t2060 = t1854 - g(2);
t2059 = t1855 - g(1);
t2058 = MDP(2) * t1881;
t2057 = MDP(6) * t1881;
t2056 = MDP(7) * t1881;
t2055 = MDP(8) * t1881;
t2054 = MDP(9) * t1881;
t1941 = t1767 * t1805 - t1771 * t1809;
t1972 = t1816 * t2001;
t1957 = t1819 * t1972;
t2012 = t1816 * t1847;
t2033 = t1759 * t1816;
t1695 = (-t1941 * t1957 + (-t1816 * t1853 + (-(-t1759 * t1844 * t1845 + t1751 * t2012) * t1817 * t1751 - (-t1751 * t1844 + t1759 * t1847) * t2033) * t2011) * t1818) * t1881;
t2053 = t1695 * t1846;
t1940 = t1768 * t1806 - t1772 * t1810;
t1970 = t1823 * t1999;
t1956 = t1832 * t1970;
t2009 = t1823 * t1863;
t2032 = t1760 * t1823;
t1696 = (-t1940 * t1956 + (-t1823 * t1853 + (-(-t1760 * t1856 * t1857 + t1752 * t2009) * t1824 * t1752 - (-t1752 * t1856 + t1760 * t1863) * t2032) * t2004) * t1831) * t1881;
t2052 = t1696 * t1862;
t1939 = t1769 * t1807 - t1773 * t1811;
t1968 = t1826 * t1997;
t1955 = t1836 * t1968;
t2007 = t1826 * t1865;
t2031 = t1761 * t1826;
t1697 = (-t1939 * t1955 + (-t1826 * t1853 + (-(-t1761 * t1858 * t1859 + t1753 * t2007) * t1827 * t1753 - (-t1753 * t1858 + t1761 * t1865) * t2031) * t2003) * t1835) * t1881;
t2051 = t1697 * t1864;
t1938 = t1770 * t1808 - t1774 * t1812;
t1966 = t1829 * t1995;
t1954 = t1840 * t1966;
t2005 = t1829 * t1867;
t2030 = t1762 * t1829;
t1698 = (-t1938 * t1954 + (-t1829 * t1853 + (-(-t1762 * t1860 * t1861 + t1754 * t2005) * t1830 * t1754 - (-t1754 * t1860 + t1762 * t1867) * t2030) * t2002) * t1839) * t1881;
t2050 = t1698 * t1866;
t2037 = t1755 * t1882;
t1985 = t1844 * t2037;
t1731 = t1941 * t1881 * t1818 + t1820 * t1985;
t2049 = t1731 * t1844;
t2048 = t1731 * t1846;
t2036 = t1756 * t1882;
t1984 = t1856 * t2036;
t1732 = t1940 * t1881 * t1831 + t1833 * t1984;
t2047 = t1732 * t1856;
t2046 = t1732 * t1862;
t2035 = t1757 * t1882;
t1983 = t1858 * t2035;
t1733 = t1939 * t1881 * t1835 + t1837 * t1983;
t2045 = t1733 * t1858;
t2044 = t1733 * t1864;
t2034 = t1758 * t1882;
t1982 = t1860 * t2034;
t1734 = t1938 * t1881 * t1839 + t1841 * t1982;
t2043 = t1734 * t1860;
t2042 = t1734 * t1866;
t2041 = t1747 * t1882;
t2040 = t1748 * t1882;
t2039 = t1749 * t1882;
t2038 = t1750 * t1882;
t2029 = t1763 * t1839;
t2028 = t1764 * t1818;
t2027 = t1765 * t1831;
t2026 = t1766 * t1835;
t2025 = t1805 * t1818;
t2024 = t1806 * t1831;
t2023 = t1807 * t1835;
t2022 = t1808 * t1839;
t2021 = t1809 * t1818;
t2019 = t1810 * t1831;
t2017 = t1811 * t1835;
t2015 = t1812 * t1839;
t1821 = 0.1e1 / t1888 ^ 2;
t1993 = t1821 * t2064;
t1992 = t1821 * t2041;
t1834 = 0.1e1 / t1897 ^ 2;
t1991 = t1834 * t2063;
t1990 = t1834 * t2040;
t1838 = 0.1e1 / t1900 ^ 2;
t1989 = t1838 * t2062;
t1988 = t1838 * t2039;
t1842 = 0.1e1 / t1903 ^ 2;
t1987 = t1842 * t2061;
t1986 = t1842 * t2038;
t1981 = t1775 * t2013;
t1980 = t1776 * t2013;
t1979 = t1777 * t2010;
t1978 = t1778 * t2008;
t1977 = t1779 * t2006;
t1976 = t1780 * t2010;
t1975 = t1781 * t2008;
t1974 = t1782 * t2006;
t1973 = t1819 * t2012;
t1971 = t1832 * t2009;
t1969 = t1836 * t2007;
t1967 = t1840 * t2005;
t1965 = t1844 * t1993;
t1964 = t1856 * t1991;
t1963 = t1858 * t1989;
t1962 = t1860 * t1987;
t1961 = t1751 * t1882 * t2033;
t1960 = t1752 * t1882 * t2032;
t1959 = t1753 * t1882 * t2031;
t1958 = t1754 * t1882 * t2030;
t1953 = t1763 * t1954;
t1952 = t1764 * t1957;
t1951 = t1765 * t1956;
t1950 = t1766 * t1955;
t1949 = t1805 * t1957;
t1948 = t1806 * t1956;
t1947 = t1807 * t1955;
t1946 = t1808 * t1954;
t1945 = t1809 * t1957;
t1944 = t1810 * t1956;
t1943 = t1811 * t1955;
t1942 = t1812 * t1954;
t1679 = t1844 * t2053 + (0.2e1 * t1818 - t1820) * t1961;
t1937 = -0.2e1 * t1679 * t1957;
t1680 = t1856 * t2052 + (0.2e1 * t1831 - t1833) * t1960;
t1936 = -0.2e1 * t1680 * t1956;
t1681 = t1858 * t2051 + (0.2e1 * t1835 - t1837) * t1959;
t1935 = -0.2e1 * t1681 * t1955;
t1682 = t1860 * t2050 + (0.2e1 * t1839 - t1841) * t1958;
t1934 = -0.2e1 * t1682 * t1954;
t1933 = 0.2e1 * t1819 * t1961;
t1932 = 0.2e1 * t1832 * t1960;
t1931 = 0.2e1 * t1836 * t1959;
t1930 = 0.2e1 * t1840 * t1958;
t1708 = g(3) * t1847 - t1719 * t1845;
t1793 = g(1) * t1809 - g(2) * t1805;
t1929 = -t1707 * t1972 + t1708 * t1844 + t1793 * t1846;
t1713 = g(3) * t1863 - t1720 * t1857;
t1794 = g(1) * t1810 - g(2) * t1806;
t1928 = -t1710 * t1970 + t1713 * t1856 + t1794 * t1862;
t1715 = g(3) * t1865 - t1721 * t1859;
t1795 = g(1) * t1811 - g(2) * t1807;
t1927 = -t1711 * t1968 + t1715 * t1858 + t1795 * t1864;
t1717 = g(3) * t1867 - t1722 * t1861;
t1796 = g(1) * t1812 - g(2) * t1808;
t1926 = -t1712 * t1966 + t1717 * t1860 + t1796 * t1866;
t1723 = -t1819 * t1985 + t2048;
t1925 = t1723 * t1957 - t1695;
t1724 = -t1832 * t1984 + t2046;
t1924 = t1724 * t1956 - t1696;
t1725 = -t1836 * t1983 + t2044;
t1923 = t1725 * t1955 - t1697;
t1726 = -t1840 * t1982 + t2042;
t1922 = t1726 * t1954 - t1698;
t1921 = t1818 * t1929;
t1920 = t1831 * t1928;
t1919 = t1835 * t1927;
t1918 = t1839 * t1926;
t1727 = t1818 * t2037 + t2049;
t1917 = -t1695 * t1818 + t1727 * t1973;
t1728 = t1831 * t2036 + t2047;
t1916 = -t1696 * t1831 + t1728 * t1971;
t1729 = t1835 * t2035 + t2045;
t1915 = -t1697 * t1835 + t1729 * t1969;
t1730 = t1839 * t2034 + t2043;
t1914 = -t1698 * t1839 + t1730 * t1967;
t1815 = t1844 ^ 2;
t1913 = -t1815 * t1707 * t1973 - t1818 * (t1708 * t1846 - t1793 * t1844);
t1822 = t1856 ^ 2;
t1912 = -t1822 * t1710 * t1971 - t1831 * (t1713 * t1862 - t1794 * t1856);
t1825 = t1858 ^ 2;
t1911 = -t1825 * t1711 * t1969 - t1835 * (t1715 * t1864 - t1795 * t1858);
t1828 = t1860 ^ 2;
t1910 = -t1828 * t1712 * t1967 - t1839 * (t1717 * t1866 - t1796 * t1860);
t1909 = t1844 * t1917;
t1908 = t1856 * t1916;
t1907 = t1858 * t1915;
t1906 = t1860 * t1914;
t1883 = t1881 * t1882;
t1784 = -t1813 * t1852 - t1814 * t1843;
t1783 = -t1813 * t1843 + t1814 * t1852;
t1746 = (-t1779 * t1788 + t1782 * t1792) * t2006;
t1745 = (-t1778 * t1787 + t1781 * t1791) * t2008;
t1744 = (-t1777 * t1786 + t1780 * t1790) * t2010;
t1743 = (-t1775 * t1785 + t1776 * t1789) * t2013;
t1742 = (t1758 * t1840 + t1987) * t1882;
t1741 = (t1757 * t1836 + t1989) * t1882;
t1740 = (t1756 * t1832 + t1991) * t1882;
t1739 = (t1755 * t1819 + t1993) * t1882;
t1738 = (-0.2e1 * t1840 + t1842) * t1830 * t2038;
t1737 = (-0.2e1 * t1836 + t1838) * t1827 * t2039;
t1736 = (-0.2e1 * t1832 + t1834) * t1824 * t2040;
t1735 = (-0.2e1 * t1819 + t1821) * t1817 * t2041;
t1694 = t1698 * t1867 - t1829 * t1986;
t1693 = t1697 * t1865 - t1826 * t1988;
t1692 = t1696 * t1863 - t1823 * t1990;
t1691 = -t1830 * t1867 * t1986 - t1698 * t1861;
t1690 = -t1827 * t1865 * t1988 - t1697 * t1859;
t1689 = -t1824 * t1863 * t1990 - t1696 * t1857;
t1688 = t1695 * t1847 - t1816 * t1992;
t1687 = -t1817 * t1847 * t1992 - t1695 * t1845;
t1686 = t1698 * t1828 + t1860 * t1930;
t1685 = t1697 * t1825 + t1858 * t1931;
t1684 = t1696 * t1822 + t1856 * t1932;
t1683 = t1695 * t1815 + t1844 * t1933;
t1678 = (t1742 * t1860 - t2042) * t1861 - t1867 * (t1698 * t1860 + t1930);
t1677 = (t1741 * t1858 - t2044) * t1859 - t1865 * (t1697 * t1858 + t1931);
t1676 = (t1740 * t1856 - t2046) * t1857 - t1863 * (t1696 * t1856 + t1932);
t1675 = (-t1742 * t1866 - t2043) * t1861 + (-0.2e1 * t1841 * t1860 * t1958 + t2050) * t1867;
t1674 = (-t1741 * t1864 - t2045) * t1859 + (-0.2e1 * t1837 * t1858 * t1959 + t2051) * t1865;
t1673 = (-t1740 * t1862 - t2047) * t1857 + (-0.2e1 * t1833 * t1856 * t1960 + t2052) * t1863;
t1672 = (t1739 * t1844 - t2048) * t1845 - t1847 * (t1695 * t1844 + t1933);
t1671 = (-t1739 * t1846 - t2049) * t1845 + (-0.2e1 * t1820 * t1844 * t1961 + t2053) * t1847;
t1 = [(t1719 * t1981 + t1720 * t1979 + t1721 * t1978 + t1722 * t1977) * MDP(1) + (t1695 * t1945 + t1696 * t1944 + t1697 * t1943 + t1698 * t1942) * t2058 + (t1688 * t1981 + t1692 * t1979 + t1693 * t1978 + t1694 * t1977 + (t1707 * t1945 + t1710 * t1944 + t1711 * t1943 + t1712 * t1942) * t1881) * MDP(3) + (t1687 * t1981 + t1689 * t1979 + t1690 * t1978 + t1691 * t1977 + (t1708 * t1945 + t1713 * t1944 + t1715 * t1943 + t1717 * t1942) * t1881) * MDP(4) + ((t1809 * t1965 + t1810 * t1964 + t1811 * t1963 + t1812 * t1962) * t1883 + (t1683 * t1945 + t1684 * t1944 + t1685 * t1943 + t1686 * t1942) * t1881) * MDP(5) + (0.2e1 * t1679 * t1945 + 0.2e1 * t1680 * t1944 + 0.2e1 * t1681 * t1943 + 0.2e1 * t1682 * t1942 - t1735 * t2021 - t1736 * t2019 - t1737 * t2017 - t1738 * t2015) * t2057 + (t1914 * t2014 + t1915 * t2016 + t1916 * t2018 + t1917 * t2020) * t2056 + (t1925 * t1809 + t1924 * t1810 + t1923 * t1811 + t1922 * t1812) * t2055 + (-t1731 * t2021 - t1732 * t2019 - t1733 * t2017 - t1734 * t2015) * t2054 + (t1671 * t1981 + t1673 * t1979 + t1674 * t1978 + t1675 * t1977 + (-t1926 * t2015 - t1927 * t2017 - t1928 * t2019 - t1929 * t2021) * t1881) * MDP(10) + (t1672 * t1981 + t1676 * t1979 + t1677 * t1978 + t1678 * t1977 + (t1913 * t1809 + t1912 * t1810 + t1911 * t1811 + t1910 * t1812) * t1881) * MDP(11) + t1784 * MDP(13) - t1783 * MDP(14) + t2059 * MDP(15); (t1719 * t1980 + t1720 * t1976 + t1721 * t1975 + t1722 * t1974) * MDP(1) + (-t1695 * t1949 - t1696 * t1948 - t1697 * t1947 - t1698 * t1946) * t2058 + (t1688 * t1980 + t1692 * t1976 + t1693 * t1975 + t1694 * t1974 + (-t1707 * t1949 - t1710 * t1948 - t1711 * t1947 - t1712 * t1946) * t1881) * MDP(3) + (t1687 * t1980 + t1689 * t1976 + t1690 * t1975 + t1691 * t1974 + (-t1708 * t1949 - t1713 * t1948 - t1715 * t1947 - t1717 * t1946) * t1881) * MDP(4) + ((-t1805 * t1965 - t1806 * t1964 - t1807 * t1963 - t1808 * t1962) * t1883 + (-t1683 * t1949 - t1684 * t1948 - t1685 * t1947 - t1686 * t1946) * t1881) * MDP(5) + (t1735 * t2025 + t1736 * t2024 + t1737 * t2023 + t1738 * t2022 + t1805 * t1937 + t1806 * t1936 + t1807 * t1935 + t1808 * t1934) * t2057 + (-t1805 * t1909 - t1806 * t1908 - t1807 * t1907 - t1808 * t1906) * t2056 + (-t1925 * t1805 - t1924 * t1806 - t1923 * t1807 - t1922 * t1808) * t2055 + (t1731 * t2025 + t1732 * t2024 + t1733 * t2023 + t1734 * t2022) * t2054 + (t1671 * t1980 + t1673 * t1976 + t1674 * t1975 + t1675 * t1974 + (t1805 * t1921 + t1806 * t1920 + t1807 * t1919 + t1808 * t1918) * t1881) * MDP(10) + (t1672 * t1980 + t1676 * t1976 + t1677 * t1975 + t1678 * t1974 + (-t1913 * t1805 - t1912 * t1806 - t1911 * t1807 - t1910 * t1808) * t1881) * MDP(11) + t1783 * MDP(13) + t1784 * MDP(14) + t2060 * MDP(15); (t1719 * t2012 + t1720 * t2009 + t1721 * t2007 + t1722 * t2005) * MDP(1) + (t1688 * t2012 + t1692 * t2009 + t1693 * t2007 + t1694 * t2005) * MDP(3) + (t1687 * t2012 + t1689 * t2009 + t1690 * t2007 + t1691 * t2005) * MDP(4) + (t1671 * t2012 + t1673 * t2009 + t1674 * t2007 + t1675 * t2005) * MDP(10) + (t1672 * t2012 + t1676 * t2009 + t1677 * t2007 + t1678 * t2005) * MDP(11) + (t1853 - g(3)) * MDP(15) + ((-t1695 * t2013 - t1696 * t2010 - t1697 * t2008 - t1698 * t2006) * MDP(2) + (-t2065 - t2066 - t2067 - t2068) * MDP(3) + (-t1708 * t2013 - t1713 * t2010 - t1715 * t2008 - t1717 * t2006) * MDP(4) + (-t1683 * t2013 - t1684 * t2010 - t1685 * t2008 - t1686 * t2006) * MDP(5) + (-t1727 * t2013 - t1728 * t2010 - t1729 * t2008 - t1730 * t2006) * MDP(7) + (-t1723 * t2013 - t1724 * t2010 - t1725 * t2008 - t1726 * t2006) * MDP(8) + (-t1707 * t1816 - t1710 * t1823 - t1711 * t1826 - t1712 * t1829) * MDP(10) + (t1844 * t2068 + t1856 * t2067 + t1858 * t2066 + t1860 * t2065) * MDP(11) + 0.2e1 * (-t1679 * t2013 - t1680 * t2010 - t1681 * t2008 - t1682 * t2006) * MDP(6)) * t1881; (t1719 * t1743 + t1720 * t1744 + t1721 * t1745 + t1722 * t1746) * MDP(1) + (-t1695 * t1952 - t1696 * t1951 - t1697 * t1950 - t1698 * t1953) * t2058 + (t1743 * t1688 + t1744 * t1692 + t1745 * t1693 + t1746 * t1694 + (-t1707 * t1952 - t1710 * t1951 - t1711 * t1950 - t1712 * t1953) * t1881) * MDP(3) + (t1743 * t1687 + t1744 * t1689 + t1745 * t1690 + t1746 * t1691 + (-t1708 * t1952 - t1713 * t1951 - t1715 * t1950 - t1717 * t1953) * t1881) * MDP(4) + ((-t1763 * t1962 - t1764 * t1965 - t1765 * t1964 - t1766 * t1963) * t1883 + (-t1683 * t1952 - t1684 * t1951 - t1685 * t1950 - t1686 * t1953) * t1881) * MDP(5) + (t1735 * t2028 + t1736 * t2027 + t1737 * t2026 + t1738 * t2029 + t1763 * t1934 + t1764 * t1937 + t1765 * t1936 + t1766 * t1935) * t2057 + (-t1763 * t1906 - t1764 * t1909 - t1765 * t1908 - t1766 * t1907) * t2056 + (-t1922 * t1763 - t1925 * t1764 - t1924 * t1765 - t1923 * t1766) * t2055 + (t1731 * t2028 + t1732 * t2027 + t1733 * t2026 + t1734 * t2029) * t2054 + (t1743 * t1671 + t1744 * t1673 + t1745 * t1674 + t1746 * t1675 + (t1763 * t1918 + t1764 * t1921 + t1765 * t1920 + t1766 * t1919) * t1881) * MDP(10) + (t1743 * t1672 + t1744 * t1676 + t1745 * t1677 + t1746 * t1678 + (-t1910 * t1763 - t1913 * t1764 - t1912 * t1765 - t1911 * t1766) * t1881) * MDP(11) + t1852 * MDP(12) + (-t2059 * t1813 + t2060 * t1814) * MDP(13) + (-t2060 * t1813 - t2059 * t1814) * MDP(14);];
tauX  = t1;
