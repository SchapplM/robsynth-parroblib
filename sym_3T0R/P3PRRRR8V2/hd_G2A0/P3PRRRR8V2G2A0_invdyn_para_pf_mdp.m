% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PRRRR8V2G2A0
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [12x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRRR8V2G2A0_convert_par2_MPV_fixb.m

% Output:
% tauX [3x1]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 17:50
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRRR8V2G2A0_invdyn_para_pf_mdp(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_mdp: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_mdp: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3PRRRR8V2G2A0_invdyn_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:50:17
% EndTime: 2020-08-06 17:50:30
% DurationCPUTime: 12.29s
% Computational Cost: add. (66792->445), mult. (135357->831), div. (5832->17), fcn. (125931->22), ass. (0->343)
t1859 = cos(qJ(2,3));
t1867 = pkin(7) + pkin(6);
t1824 = t1859 * t1867;
t1853 = sin(qJ(2,3));
t1801 = pkin(2) * t1853 - t1824;
t1843 = sin(pkin(4));
t1845 = cos(pkin(4));
t1852 = sin(qJ(3,3));
t1980 = t1845 * t1852;
t1778 = pkin(3) * t1980 + t1801 * t1843;
t1858 = cos(qJ(3,3));
t1994 = t1843 * t1853;
t1839 = t1858 ^ 2;
t2035 = pkin(3) * t1839;
t1745 = pkin(2) * t1980 + t1778 * t1858 + t1994 * t2035;
t1736 = 0.1e1 / t1745;
t1861 = cos(qJ(2,2));
t1825 = t1861 * t1867;
t1855 = sin(qJ(2,2));
t1802 = pkin(2) * t1855 - t1825;
t1854 = sin(qJ(3,2));
t1978 = t1845 * t1854;
t1779 = pkin(3) * t1978 + t1802 * t1843;
t1860 = cos(qJ(3,2));
t1992 = t1843 * t1855;
t1840 = t1860 ^ 2;
t2034 = pkin(3) * t1840;
t1746 = pkin(2) * t1978 + t1779 * t1860 + t1992 * t2034;
t1739 = 0.1e1 / t1746;
t1863 = cos(qJ(2,1));
t1826 = t1863 * t1867;
t1857 = sin(qJ(2,1));
t1803 = pkin(2) * t1857 - t1826;
t1856 = sin(qJ(3,1));
t1976 = t1845 * t1856;
t1780 = pkin(3) * t1976 + t1803 * t1843;
t1862 = cos(qJ(3,1));
t1990 = t1843 * t1857;
t1841 = t1862 ^ 2;
t2033 = pkin(3) * t1841;
t1747 = pkin(2) * t1976 + t1780 * t1862 + t1990 * t2033;
t1742 = 0.1e1 / t1747;
t1846 = legFrame(3,2);
t1830 = sin(t1846);
t1833 = cos(t1846);
t1865 = xDP(2);
t1866 = xDP(1);
t1790 = t1830 * t1865 - t1833 * t1866;
t1842 = sin(pkin(8));
t1864 = xDP(3);
t1817 = t1864 * t1842;
t1844 = cos(pkin(8));
t1766 = t1790 * t1844 + t1817;
t1981 = t1844 * t1864;
t1769 = t1790 * t1842 - t1981;
t1979 = t1845 * t1853;
t1989 = t1843 * t1858;
t1718 = t1769 * t1989 - (t1766 * t1859 - t1769 * t1979) * t1852;
t2018 = t1718 * t1736;
t1950 = t1867 * t2018;
t1915 = t1852 * t1950;
t2032 = pkin(3) * t1858;
t1818 = pkin(2) + t2032;
t1793 = t1818 * t1853 - t1824;
t1821 = t1853 * t1867;
t1724 = t1769 * (t1818 * t1859 + t1821) * t1845 + t1766 * t1793;
t1798 = t1818 * t1980;
t1757 = t1793 * t1989 + t1798;
t2015 = t1724 / t1757;
t1694 = t1915 - t2015;
t1804 = pkin(2) * t1859 + t1821;
t1995 = t1843 * t1852;
t1897 = pkin(3) * t1995 - t1801 * t1845;
t1748 = -t1804 * t1842 + t1897 * t1844;
t1787 = t1842 * t1859 + t1844 * t1979;
t1996 = t1843 * t1844;
t2039 = pkin(2) * t1852;
t1706 = -(t1787 * t1830 - t1833 * t1994) * t2035 + (t1748 * t1830 + t1778 * t1833) * t1858 + (t1830 * t1996 + t1833 * t1845) * t2039;
t1707 = (t1787 * t1833 + t1830 * t1994) * t2035 + (-t1748 * t1833 + t1778 * t1830) * t1858 + (t1830 * t1845 - t1833 * t1996) * t2039;
t1784 = t1842 * t1979 - t1844 * t1859;
t1984 = t1844 * t1858;
t1727 = -t1784 * t2035 + t1804 * t1984 + (pkin(2) * t1995 + t1897 * t1858) * t1842;
t1850 = xDDP(2);
t1851 = xDDP(1);
t1869 = 0.1e1 / pkin(3);
t1940 = t1869 * t2015;
t1942 = t1736 * t2015;
t1737 = 0.1e1 / t1745 ^ 2;
t1827 = pkin(2) ^ 2 + t1867 ^ 2;
t1868 = pkin(3) ^ 2;
t1941 = t1852 * t2015;
t2048 = 0.2e1 * pkin(2);
t2028 = pkin(3) * t2048;
t1954 = (-t1867 * t1941 + (t1839 * t1868 + t1858 * t2028 + t1827) * t2018) * t1718 * t1737;
t1849 = xDDP(3);
t2006 = t1736 * t1849;
t1676 = t1858 * t1954 + (pkin(2) * t1940 - t1694 * t1858) * t1942 + t1727 * t2006 + (t1706 * t1850 + t1707 * t1851) * t1736;
t1985 = t1844 * t1845;
t1796 = -t1843 * g(1) + g(2) * t1985;
t1797 = g(1) * t1985 + g(2) * t1843;
t1997 = t1842 * t1845;
t1810 = g(3) * t1997;
t1885 = t1676 * t1843 + t1796 * t1830 - t1797 * t1833 + t1810;
t1664 = t1885 * t1859;
t1829 = g(3) * t1844;
t1912 = g(1) * t1833 - g(2) * t1830;
t2052 = t1912 * t1842 + t1829;
t2072 = t1853 * t2052 + t1664;
t1847 = legFrame(2,2);
t1831 = sin(t1847);
t1834 = cos(t1847);
t1791 = t1831 * t1865 - t1834 * t1866;
t1767 = t1791 * t1844 + t1817;
t1770 = t1791 * t1842 - t1981;
t1977 = t1845 * t1855;
t1988 = t1843 * t1860;
t1719 = t1770 * t1988 - (t1767 * t1861 - t1770 * t1977) * t1854;
t2017 = t1719 * t1739;
t1948 = t1867 * t2017;
t1914 = t1854 * t1948;
t2031 = pkin(3) * t1860;
t1819 = pkin(2) + t2031;
t1794 = t1819 * t1855 - t1825;
t1822 = t1855 * t1867;
t1725 = t1770 * (t1819 * t1861 + t1822) * t1845 + t1767 * t1794;
t1799 = t1819 * t1978;
t1758 = t1794 * t1988 + t1799;
t2014 = t1725 / t1758;
t1695 = t1914 - t2014;
t1805 = pkin(2) * t1861 + t1822;
t1993 = t1843 * t1854;
t1896 = pkin(3) * t1993 - t1802 * t1845;
t1749 = -t1805 * t1842 + t1896 * t1844;
t1788 = t1842 * t1861 + t1844 * t1977;
t2038 = pkin(2) * t1854;
t1708 = -(t1788 * t1831 - t1834 * t1992) * t2034 + (t1749 * t1831 + t1779 * t1834) * t1860 + (t1831 * t1996 + t1834 * t1845) * t2038;
t1709 = (t1788 * t1834 + t1831 * t1992) * t2034 + (-t1749 * t1834 + t1779 * t1831) * t1860 + (t1831 * t1845 - t1834 * t1996) * t2038;
t1785 = t1842 * t1977 - t1844 * t1861;
t1983 = t1844 * t1860;
t1728 = -t1785 * t2034 + t1805 * t1983 + (pkin(2) * t1993 + t1896 * t1860) * t1842;
t1937 = t1869 * t2014;
t1939 = t1739 * t2014;
t1740 = 0.1e1 / t1746 ^ 2;
t1938 = t1854 * t2014;
t1953 = (-t1867 * t1938 + (t1840 * t1868 + t1860 * t2028 + t1827) * t2017) * t1719 * t1740;
t2005 = t1739 * t1849;
t1677 = t1860 * t1953 + (pkin(2) * t1937 - t1695 * t1860) * t1939 + t1728 * t2005 + (t1708 * t1850 + t1709 * t1851) * t1739;
t1884 = t1677 * t1843 + t1796 * t1831 - t1797 * t1834 + t1810;
t1665 = t1884 * t1861;
t1911 = g(1) * t1834 - g(2) * t1831;
t2053 = t1911 * t1842 + t1829;
t2071 = t1855 * t2053 + t1665;
t1848 = legFrame(1,2);
t1832 = sin(t1848);
t1835 = cos(t1848);
t1792 = t1832 * t1865 - t1835 * t1866;
t1768 = t1792 * t1844 + t1817;
t1771 = t1792 * t1842 - t1981;
t1975 = t1845 * t1857;
t1987 = t1843 * t1862;
t1720 = t1771 * t1987 - (t1768 * t1863 - t1771 * t1975) * t1856;
t2016 = t1720 * t1742;
t1946 = t1867 * t2016;
t1913 = t1856 * t1946;
t2030 = pkin(3) * t1862;
t1820 = pkin(2) + t2030;
t1795 = t1820 * t1857 - t1826;
t1823 = t1857 * t1867;
t1726 = t1771 * (t1820 * t1863 + t1823) * t1845 + t1768 * t1795;
t1800 = t1820 * t1976;
t1759 = t1795 * t1987 + t1800;
t2013 = t1726 / t1759;
t1696 = t1913 - t2013;
t1806 = pkin(2) * t1863 + t1823;
t1991 = t1843 * t1856;
t1895 = pkin(3) * t1991 - t1803 * t1845;
t1750 = -t1806 * t1842 + t1895 * t1844;
t1789 = t1842 * t1863 + t1844 * t1975;
t2037 = pkin(2) * t1856;
t1710 = -(t1789 * t1832 - t1835 * t1990) * t2033 + (t1750 * t1832 + t1780 * t1835) * t1862 + (t1832 * t1996 + t1835 * t1845) * t2037;
t1711 = (t1789 * t1835 + t1832 * t1990) * t2033 + (-t1750 * t1835 + t1780 * t1832) * t1862 + (t1832 * t1845 - t1835 * t1996) * t2037;
t1786 = t1842 * t1975 - t1844 * t1863;
t1982 = t1844 * t1862;
t1729 = -t1786 * t2033 + t1806 * t1982 + (pkin(2) * t1991 + t1895 * t1862) * t1842;
t1934 = t1869 * t2013;
t1936 = t1742 * t2013;
t1743 = 0.1e1 / t1747 ^ 2;
t1935 = t1856 * t2013;
t1952 = (-t1867 * t1935 + (t1841 * t1868 + t1862 * t2028 + t1827) * t2016) * t1720 * t1743;
t2004 = t1742 * t1849;
t1678 = t1862 * t1952 + (pkin(2) * t1934 - t1696 * t1862) * t1936 + t1729 * t2004 + (t1710 * t1850 + t1711 * t1851) * t1742;
t1883 = t1678 * t1843 + t1796 * t1832 - t1797 * t1835 + t1810;
t1666 = t1883 * t1863;
t1910 = g(1) * t1835 - g(2) * t1832;
t2054 = t1910 * t1842 + t1829;
t2070 = t1857 * t2054 + t1666;
t2063 = t2052 * t1859;
t2062 = t2053 * t1861;
t2061 = t2054 * t1863;
t2021 = t1718 ^ 2 * t1737;
t2020 = t1719 ^ 2 * t1740;
t2019 = t1720 ^ 2 * t1743;
t1760 = t1784 * t1852 + t1842 * t1989;
t1763 = -t1787 * t1852 - t1843 * t1984;
t1951 = t1859 * t2018;
t1968 = t1853 * t1858;
t2049 = t1830 * t1850 - t1833 * t1851;
t1643 = t1763 * t2006 - (t1843 * t1859 * t1940 + (t1845 * t1839 - t1968 * t1995 - t1845) * t2018) * t1942 + (-((t1843 * t1951 + t1845 * t1940) * t2035 + ((-t1941 + t1950) * t1853 + pkin(2) * t1951) * t1989 + t1845 * t1694) * t2018 + t2049 * t1760) * t1736;
t1781 = pkin(3) * t1968 + t1801;
t1970 = t1849 * t1869;
t1971 = t1845 * t1869;
t1986 = t1843 * t1869;
t1974 = t1845 * t1859;
t1733 = (t1842 * t1974 + t1844 * t1853) * t2032 + t1804 * t1997 + t1801 * t1844;
t2009 = t1733 * t1736;
t1730 = (t1842 * t1853 - t1844 * t1974) * t2032 - t1804 * t1985 + t1842 * t1801;
t2012 = t1730 * t1736;
t2036 = pkin(2) * t1869;
t1679 = t1970 * t2012 - t1954 * t1971 - (-t1845 * t1915 + (-t1781 * t1852 * t1986 + t1845 * (t1858 * t2036 + t1839)) * t2015) / (t1781 * t1989 + t1798) * t1940 + t2049 * t1869 * t2009;
t1870 = 0.1e1 / pkin(3) ^ 2;
t1945 = t1724 ^ 2 / t1757 ^ 2 * t1870;
t1652 = t1679 * t1852 + t1858 * t1945;
t1653 = t1679 * t1858 - t1852 * t1945;
t1891 = -pkin(6) * t1945 + t1643 * t2048;
t1906 = t1940 * t2018;
t1900 = 0.2e1 * t1906;
t1969 = t1852 * t1858;
t2027 = t1643 * t1852;
t2042 = 0.2e1 * t1839 - 0.1e1;
t2043 = pkin(6) / 0.2e1;
t2046 = -0.2e1 * pkin(2) * t1906 - 0.2e1 * t1679 * t2043;
t2047 = 2 * MDP(6);
t1879 = MDP(10) * ((t1664 + t1891) * t1858 + t1852 * t2046 + t2052 * t1968) + MDP(11) * (t1858 * t2046 + (-t1891 - t2072) * t1852) + MDP(2) * t1643 + MDP(3) * t2072 + MDP(4) * (-t1853 * t1885 + t2063) + MDP(5) * (t1858 * t1900 + t2027) * t1852 + MDP(7) * t1652 + MDP(8) * t1653 + (t1643 * t1969 + t2042 * t1906) * t2047;
t2060 = t1760 * t1879;
t1761 = t1785 * t1854 + t1842 * t1988;
t1764 = -t1788 * t1854 - t1843 * t1983;
t1949 = t1861 * t2017;
t1966 = t1855 * t1860;
t2050 = t1831 * t1850 - t1834 * t1851;
t1644 = t1764 * t2005 - (t1843 * t1861 * t1937 + (t1840 * t1845 - t1966 * t1993 - t1845) * t2017) * t1939 + (-((t1843 * t1949 + t1845 * t1937) * t2034 + ((-t1938 + t1948) * t1855 + pkin(2) * t1949) * t1988 + t1845 * t1695) * t2017 + t2050 * t1761) * t1739;
t1782 = pkin(3) * t1966 + t1802;
t1973 = t1845 * t1861;
t1734 = (t1842 * t1973 + t1844 * t1855) * t2031 + t1805 * t1997 + t1802 * t1844;
t2008 = t1734 * t1739;
t1731 = (t1842 * t1855 - t1844 * t1973) * t2031 - t1805 * t1985 + t1842 * t1802;
t2011 = t1731 * t1739;
t1680 = t1970 * t2011 - t1953 * t1971 - (-t1845 * t1914 + (-t1782 * t1854 * t1986 + t1845 * (t1860 * t2036 + t1840)) * t2014) / (t1782 * t1988 + t1799) * t1937 + t2050 * t1869 * t2008;
t1944 = t1725 ^ 2 / t1758 ^ 2 * t1870;
t1654 = t1680 * t1854 + t1860 * t1944;
t1655 = t1680 * t1860 - t1854 * t1944;
t1890 = -pkin(6) * t1944 + t1644 * t2048;
t1905 = t1937 * t2017;
t1899 = 0.2e1 * t1905;
t1967 = t1854 * t1860;
t2025 = t1644 * t1854;
t2041 = 0.2e1 * t1840 - 0.1e1;
t2045 = -0.2e1 * pkin(2) * t1905 - 0.2e1 * t1680 * t2043;
t1878 = MDP(10) * ((t1665 + t1890) * t1860 + t1854 * t2045 + t2053 * t1966) + MDP(11) * (t1860 * t2045 + (-t1890 - t2071) * t1854) + MDP(2) * t1644 + MDP(3) * t2071 + MDP(4) * (-t1855 * t1884 + t2062) + MDP(5) * (t1860 * t1899 + t2025) * t1854 + MDP(7) * t1654 + MDP(8) * t1655 + (t1644 * t1967 + t2041 * t1905) * t2047;
t2059 = t1761 * t1878;
t1762 = t1786 * t1856 + t1842 * t1987;
t1765 = -t1789 * t1856 - t1843 * t1982;
t1947 = t1863 * t2016;
t1964 = t1857 * t1862;
t2051 = t1832 * t1850 - t1835 * t1851;
t1645 = t1765 * t2004 - (t1843 * t1863 * t1934 + (t1841 * t1845 - t1964 * t1991 - t1845) * t2016) * t1936 + (-((t1843 * t1947 + t1845 * t1934) * t2033 + ((-t1935 + t1946) * t1857 + pkin(2) * t1947) * t1987 + t1845 * t1696) * t2016 + t2051 * t1762) * t1742;
t1783 = pkin(3) * t1964 + t1803;
t1972 = t1845 * t1863;
t1735 = (t1842 * t1972 + t1844 * t1857) * t2030 + t1806 * t1997 + t1803 * t1844;
t2007 = t1735 * t1742;
t1732 = (t1842 * t1857 - t1844 * t1972) * t2030 - t1806 * t1985 + t1842 * t1803;
t2010 = t1732 * t1742;
t1681 = t1970 * t2010 - t1952 * t1971 - (-t1845 * t1913 + (-t1783 * t1856 * t1986 + t1845 * (t1862 * t2036 + t1841)) * t2013) / (t1783 * t1987 + t1800) * t1934 + t2051 * t1869 * t2007;
t1943 = t1726 ^ 2 / t1759 ^ 2 * t1870;
t1656 = t1681 * t1856 + t1862 * t1943;
t1657 = t1681 * t1862 - t1856 * t1943;
t1889 = -pkin(6) * t1943 + t1645 * t2048;
t1904 = t1934 * t2016;
t1898 = 0.2e1 * t1904;
t1965 = t1856 * t1862;
t2023 = t1645 * t1856;
t2040 = 0.2e1 * t1841 - 0.1e1;
t2044 = -0.2e1 * pkin(2) * t1904 - 0.2e1 * t1681 * t2043;
t1877 = MDP(10) * ((t1666 + t1889) * t1862 + t1856 * t2044 + t2054 * t1964) + MDP(11) * (t1862 * t2044 + (-t1889 - t2070) * t1856) + MDP(2) * t1645 + MDP(3) * t2070 + MDP(4) * (-t1857 * t1883 + t2061) + MDP(5) * (t1862 * t1898 + t2023) * t1856 + MDP(7) * t1656 + MDP(8) * t1657 + (t1645 * t1965 + t2040 * t1904) * t2047;
t2058 = t1762 * t1877;
t2029 = g(3) * t1842;
t2026 = t1643 * t1859;
t2024 = t1644 * t1861;
t2022 = t1645 * t1863;
t1960 = t1736 * t2027;
t1959 = t1643 * t1736 * t1858;
t1958 = t1739 * t2025;
t1957 = t1644 * t1739 * t1860;
t1956 = t1742 * t2023;
t1955 = t1645 * t1742 * t1862;
t1933 = t1830 * t2009;
t1932 = t1833 * t2009;
t1930 = t1831 * t2008;
t1929 = t1834 * t2008;
t1927 = t1832 * t2007;
t1926 = t1835 * t2007;
t1924 = t1733 * t1960;
t1923 = t1733 * t1959;
t1922 = t1734 * t1958;
t1921 = t1734 * t1957;
t1920 = t1735 * t1956;
t1919 = t1735 * t1955;
t1918 = t1736 * t1969 * t2021;
t1917 = t1739 * t1967 * t2020;
t1916 = t1742 * t1965 * t2019;
t1909 = t1733 * t1918;
t1908 = t1734 * t1917;
t1907 = t1735 * t1916;
t1903 = -(t1945 + t2021) * t1853 + t2026;
t1902 = -(t1944 + t2020) * t1855 + t2024;
t1901 = -(t1943 + t2019) * t1857 + t2022;
t1646 = t1679 * t1853 + t1859 * t1900;
t1673 = -g(1) * t1830 - g(2) * t1833 + t1676;
t1882 = t1673 * MDP(1) + ((-t1852 * t1646 + t1903 * t1858) * t1843 + t1845 * t1653) * MDP(10) + ((-t1858 * t1646 - t1903 * t1852) * t1843 - t1845 * t1652) * MDP(11) + (-MDP(3) * (t1853 * t2021 - t2026) - MDP(4) * (t1643 * t1853 + t1859 * t2021)) * t1843;
t1647 = t1680 * t1855 + t1861 * t1899;
t1674 = -g(1) * t1831 - g(2) * t1834 + t1677;
t1881 = t1674 * MDP(1) + ((-t1854 * t1647 + t1902 * t1860) * t1843 + t1845 * t1655) * MDP(10) + ((-t1860 * t1647 - t1902 * t1854) * t1843 - t1845 * t1654) * MDP(11) + (-MDP(3) * (t1855 * t2020 - t2024) - MDP(4) * (t1644 * t1855 + t1861 * t2020)) * t1843;
t1648 = t1681 * t1857 + t1863 * t1898;
t1675 = -g(1) * t1832 - g(2) * t1835 + t1678;
t1880 = t1675 * MDP(1) + ((-t1856 * t1648 + t1901 * t1862) * t1843 + t1845 * t1657) * MDP(10) + ((-t1862 * t1648 - t1901 * t1856) * t1843 - t1845 * t1656) * MDP(11) + (-MDP(3) * (t1857 * t2019 - t2022) - MDP(4) * (t1645 * t1857 + t1863 * t2019)) * t1843;
t1774 = t1910 * t1844 - t2029;
t1773 = t1911 * t1844 - t2029;
t1772 = t1912 * t1844 - t2029;
t1693 = t2040 * t2019;
t1692 = t2041 * t2020;
t1691 = t2042 * t2021;
t1672 = t1675 * t1845 + t1774 * t1843;
t1671 = t1674 * t1845 + t1773 * t1843;
t1670 = t1673 * t1845 + t1772 * t1843;
t1630 = pkin(2) * t2019 - pkin(6) * t1645 + (-t1675 * t1843 + t1774 * t1845) * t1857 + t2061;
t1629 = pkin(2) * t2020 - pkin(6) * t1644 + (-t1674 * t1843 + t1773 * t1845) * t1855 + t2062;
t1628 = pkin(2) * t2021 - pkin(6) * t1643 + (-t1673 * t1843 + t1772 * t1845) * t1853 + t2063;
t1627 = t1630 * t1862 - t1672 * t1856;
t1626 = t1630 * t1856 + t1672 * t1862;
t1625 = t1629 * t1860 - t1671 * t1854;
t1624 = t1629 * t1854 + t1671 * t1860;
t1623 = t1628 * t1858 - t1670 * t1852;
t1622 = t1628 * t1852 + t1670 * t1858;
t1 = [(t1851 - g(1)) * MDP(12) + (t1880 * t1711 - t1835 * t2058) * t1742 + (t1881 * t1709 - t1834 * t2059) * t1739 + (t1882 * t1707 - t1833 * t2060) * t1736 + ((t1833 * t1909 + t1834 * t1908 + t1835 * t1907) * MDP(5) + (t1691 * t1932 + t1692 * t1929 + t1693 * t1926) * MDP(6) + (-t1833 * t1924 - t1834 * t1922 - t1835 * t1920) * MDP(7) + (-t1833 * t1923 - t1834 * t1921 - t1835 * t1919) * MDP(8) + (-t1679 * t1932 - t1680 * t1929 - t1681 * t1926) * MDP(9) + (-t1622 * t1932 - t1624 * t1929 - t1626 * t1926) * MDP(10) + (-t1623 * t1932 - t1625 * t1929 - t1627 * t1926) * MDP(11)) * t1869; (t1850 - g(2)) * MDP(12) + (t1880 * t1710 + t1832 * t2058) * t1742 + (t1881 * t1708 + t1831 * t2059) * t1739 + (t1882 * t1706 + t1830 * t2060) * t1736 + ((-t1830 * t1909 - t1831 * t1908 - t1832 * t1907) * MDP(5) + (-t1691 * t1933 - t1692 * t1930 - t1693 * t1927) * MDP(6) + (t1830 * t1924 + t1831 * t1922 + t1832 * t1920) * MDP(7) + (t1830 * t1923 + t1831 * t1921 + t1832 * t1919) * MDP(8) + (t1679 * t1933 + t1680 * t1930 + t1681 * t1927) * MDP(9) + (t1622 * t1933 + t1624 * t1930 + t1626 * t1927) * MDP(10) + (t1623 * t1933 + t1625 * t1930 + t1627 * t1927) * MDP(11)) * t1869; (t1849 - g(3)) * MDP(12) + (t1880 * t1729 + t1877 * t1765) * t1742 + (t1881 * t1728 + t1878 * t1764) * t1739 + (t1882 * t1727 + t1879 * t1763) * t1736 + ((-t1730 * t1918 - t1731 * t1917 - t1732 * t1916) * MDP(5) + (-t1691 * t2012 - t1692 * t2011 - t1693 * t2010) * MDP(6) + (t1730 * t1960 + t1731 * t1958 + t1732 * t1956) * MDP(7) + (t1730 * t1959 + t1731 * t1957 + t1732 * t1955) * MDP(8) + (t1679 * t2012 + t1680 * t2011 + t1681 * t2010) * MDP(9) + (t1622 * t2012 + t1624 * t2011 + t1626 * t2010) * MDP(10) + (t1623 * t2012 + t1625 * t2011 + t1627 * t2010) * MDP(11)) * t1869;];
tauX  = t1;
