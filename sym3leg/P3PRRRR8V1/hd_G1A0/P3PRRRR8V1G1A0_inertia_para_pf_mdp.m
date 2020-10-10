% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3PRRRR8V1G1A0
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
%   see P3PRRRR8V1G1A0_convert_par2_MPV_fixb.m

% Output:
% MMX [3x3]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 16:50
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P3PRRRR8V1G1A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G1A0_inertia_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G1A0_inertia_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G1A0_inertia_para_pf_mdp: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G1A0_inertia_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G1A0_inertia_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3PRRRR8V1G1A0_inertia_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 16:50:37
% EndTime: 2020-08-06 16:50:41
% DurationCPUTime: 4.35s
% Computational Cost: add. (4729->313), mult. (13109->685), div. (876->14), fcn. (15234->22), ass. (0->310)
t1892 = cos(qJ(2,3));
t1886 = sin(qJ(2,3));
t1891 = cos(qJ(3,3));
t2006 = t1886 * t1891;
t1854 = pkin(2) * t2006 - pkin(5) * t1892;
t1879 = sin(pkin(3));
t1881 = cos(pkin(3));
t1885 = sin(qJ(3,3));
t2017 = t1881 * t1885;
t1959 = pkin(2) * t2017 + t1854 * t1879;
t1828 = 0.1e1 / t1959 ^ 2;
t1873 = 0.1e1 / t1891 ^ 2;
t2094 = t1828 * t1873;
t1894 = cos(qJ(2,2));
t1888 = sin(qJ(2,2));
t1893 = cos(qJ(3,2));
t2004 = t1888 * t1893;
t1855 = pkin(2) * t2004 - pkin(5) * t1894;
t1887 = sin(qJ(3,2));
t2015 = t1881 * t1887;
t1958 = pkin(2) * t2015 + t1855 * t1879;
t1830 = 0.1e1 / t1958 ^ 2;
t1875 = 0.1e1 / t1893 ^ 2;
t2093 = t1830 * t1875;
t1896 = cos(qJ(2,1));
t1890 = sin(qJ(2,1));
t1895 = cos(qJ(3,1));
t2002 = t1890 * t1895;
t1856 = pkin(2) * t2002 - pkin(5) * t1896;
t1889 = sin(qJ(3,1));
t2013 = t1881 * t1889;
t1957 = pkin(2) * t2013 + t1856 * t1879;
t1832 = 0.1e1 / t1957 ^ 2;
t1877 = 0.1e1 / t1895 ^ 2;
t2092 = t1832 * t1877;
t2091 = 0.1e1 / t1957;
t2090 = 0.1e1 / t1958;
t2089 = 0.1e1 / t1959;
t1882 = legFrame(3,3);
t1863 = sin(t1882);
t1866 = cos(t1882);
t1878 = sin(pkin(6));
t1880 = cos(pkin(6));
t1836 = -t1863 * t1878 + t1866 * t1880;
t2016 = t1881 * t1886;
t1843 = t1878 * t1892 + t1880 * t2016;
t1846 = -t1878 * t2016 + t1880 * t1892;
t2024 = t1879 * t1891;
t1806 = (-t1843 * t1866 - t1846 * t1863) * t1885 - t1836 * t2024;
t1839 = t1863 * t1880 + t1866 * t1878;
t1809 = (-t1843 * t1863 + t1846 * t1866) * t1885 - t1839 * t2024;
t2088 = t1806 * t1809;
t1883 = legFrame(2,3);
t1864 = sin(t1883);
t1867 = cos(t1883);
t1837 = -t1864 * t1878 + t1867 * t1880;
t2014 = t1881 * t1888;
t1844 = t1878 * t1894 + t1880 * t2014;
t1847 = -t1878 * t2014 + t1880 * t1894;
t2022 = t1879 * t1893;
t1807 = (-t1844 * t1867 - t1847 * t1864) * t1887 - t1837 * t2022;
t1840 = t1864 * t1880 + t1867 * t1878;
t1810 = (-t1844 * t1864 + t1847 * t1867) * t1887 - t1840 * t2022;
t2087 = t1807 * t1810;
t1884 = legFrame(1,3);
t1865 = sin(t1884);
t1868 = cos(t1884);
t1838 = -t1865 * t1878 + t1868 * t1880;
t2012 = t1881 * t1890;
t1842 = t1878 * t2012 - t1880 * t1896;
t1845 = t1878 * t1896 + t1880 * t2012;
t2020 = t1879 * t1895;
t1808 = (t1842 * t1865 - t1845 * t1868) * t1889 - t1838 * t2020;
t1841 = t1865 * t1880 + t1868 * t1878;
t1811 = (-t1842 * t1868 - t1845 * t1865) * t1889 - t1841 * t2020;
t2086 = t1808 * t1811;
t2077 = pkin(2) * t1891;
t1857 = pkin(5) * t1886 + t1892 * t2077;
t2078 = pkin(2) * t1879;
t1908 = -t1854 * t1881 + t1885 * t2078;
t1812 = t1857 * t1880 + t1908 * t1878;
t1815 = t1857 * t1878 - t1908 * t1880;
t1791 = t1812 * t1863 + t1815 * t1866;
t2063 = t1791 * t1806;
t1788 = t1812 * t1866 - t1815 * t1863;
t2069 = t1788 * t1809;
t2085 = t1828 * (t2063 + t2069);
t2076 = pkin(2) * t1893;
t1858 = pkin(5) * t1888 + t1894 * t2076;
t1907 = -t1855 * t1881 + t1887 * t2078;
t1813 = t1858 * t1880 + t1907 * t1878;
t1816 = t1858 * t1878 - t1907 * t1880;
t1792 = t1813 * t1864 + t1816 * t1867;
t2061 = t1792 * t1807;
t1789 = t1813 * t1867 - t1816 * t1864;
t2067 = t1789 * t1810;
t2084 = t1830 * (t2061 + t2067);
t2075 = pkin(2) * t1895;
t1859 = pkin(5) * t1890 + t1896 * t2075;
t1906 = -t1856 * t1881 + t1889 * t2078;
t1814 = t1859 * t1880 + t1906 * t1878;
t1817 = t1859 * t1878 - t1906 * t1880;
t1793 = t1814 * t1865 + t1817 * t1868;
t2059 = t1793 * t1808;
t1790 = t1814 * t1868 - t1817 * t1865;
t2065 = t1790 * t1811;
t2083 = t1832 * (t2059 + t2065);
t2082 = t1873 * t1885;
t2081 = t1875 * t1887;
t2080 = t1877 * t1889;
t2079 = 2 * MDP(6);
t2074 = MDP(3) * t1879;
t2073 = MDP(4) * t1879;
t1897 = 0.1e1 / pkin(2);
t2072 = MDP(7) * t1897;
t2071 = MDP(8) * t1897;
t2070 = MDP(9) / pkin(2) ^ 2;
t2068 = t1788 * t2089;
t2066 = t1789 * t2090;
t2064 = t1790 * t2091;
t2062 = t1791 * t2089;
t2060 = t1792 * t2090;
t2058 = t1793 * t2091;
t2057 = t2089 ^ 2;
t1872 = 0.1e1 / t1891;
t2056 = t2089 * t1872;
t2055 = t2089 * t1886;
t2053 = t2090 ^ 2;
t1874 = 0.1e1 / t1893;
t2052 = t2090 * t1874;
t2051 = t2090 * t1888;
t2049 = t2091 ^ 2;
t1876 = 0.1e1 / t1895;
t2048 = t2091 * t1876;
t2047 = t2091 * t1890;
t2045 = t2089 * t1892;
t2044 = t1828 * t1872;
t2042 = t1828 * t1892;
t2041 = t2090 * t1894;
t2040 = t1830 * t1874;
t2038 = t1830 * t1894;
t2037 = t2091 * t1896;
t2036 = t1832 * t1876;
t2034 = t1832 * t1896;
t2033 = t1872 * t1885;
t2032 = t1872 * t1886;
t2031 = t1872 * t1892;
t2030 = t1874 * t1887;
t2029 = t1874 * t1888;
t2028 = t1874 * t1894;
t2027 = t1876 * t1889;
t2026 = t1876 * t1890;
t2025 = t1876 * t1896;
t2023 = t1879 * t1892;
t2021 = t1879 * t1894;
t2019 = t1879 * t1896;
t2018 = t1879 * t1897;
t2011 = t1881 * t1892;
t2010 = t1881 * t1894;
t2009 = t1881 * t1896;
t2008 = t1881 * t1897;
t2007 = t1885 * t1886;
t2005 = t1887 * t1888;
t2003 = t1889 * t1890;
t2001 = 0.2e1 * t2074;
t2000 = 0.2e1 * t2073;
t1999 = 0.2e1 * t2072;
t1998 = 0.2e1 * t2071;
t1997 = t2089 * t2056;
t1848 = -t1879 * t2007 + t1881 * t1891;
t1996 = t1848 * t2056;
t1850 = t1879 * t2006 + t2017;
t1995 = t1850 * t2056;
t1994 = t1897 * t2056;
t1993 = t2089 * t2008;
t1992 = t2090 * t2052;
t1849 = -t1879 * t2005 + t1881 * t1893;
t1991 = t1849 * t2052;
t1851 = t1879 * t2004 + t2015;
t1990 = t1851 * t2052;
t1989 = t1897 * t2052;
t1988 = t2090 * t2008;
t1987 = t2091 * t2048;
t1852 = -t1879 * t2003 + t1881 * t1895;
t1986 = t1852 * t2048;
t1853 = t1879 * t2002 + t2013;
t1985 = t1853 * t2048;
t1984 = t1897 * t2048;
t1983 = t2091 * t2008;
t1982 = t2089 * t2032;
t1981 = t2089 * t2031;
t1980 = t2089 * t2023;
t1869 = t1885 ^ 2;
t1979 = t1869 * t2094;
t1978 = t1828 * t2033;
t1977 = t1828 * t2023;
t1976 = t2090 * t2029;
t1975 = t2090 * t2028;
t1974 = t2090 * t2021;
t1870 = t1887 ^ 2;
t1973 = t1870 * t2093;
t1972 = t1830 * t2030;
t1971 = t1830 * t2021;
t1970 = t2091 * t2026;
t1969 = t2091 * t2025;
t1968 = t2091 * t2019;
t1871 = t1889 ^ 2;
t1967 = t1871 * t2092;
t1966 = t1832 * t2027;
t1965 = t1832 * t2019;
t1964 = t1885 * t2031;
t1963 = t1887 * t2028;
t1962 = t1889 * t2025;
t1961 = (t2064 + t2066 + t2068) * MDP(1) + (t1806 * t1981 + t1807 * t1975 + t1808 * t1969) * t2074 + (-t1806 * t1982 - t1807 * t1976 - t1808 * t1970) * t2073;
t1960 = (t2058 + t2060 + t2062) * MDP(1) + (t1809 * t1981 + t1810 * t1975 + t1811 * t1969) * t2074 + (-t1809 * t1982 - t1810 * t1976 - t1811 * t1970) * t2073;
t1956 = t1788 * t1806 * t2044;
t1955 = t1789 * t1807 * t2040;
t1954 = t1790 * t1808 * t2036;
t1953 = t1791 * t1809 * t2044;
t1952 = t1792 * t1810 * t2040;
t1951 = t1793 * t1811 * t2036;
t1794 = -(t1836 * t2011 - t1839 * t1886) * t2077 - pkin(5) * (t1836 * t2016 + t1839 * t1892);
t1950 = t1794 * t1997;
t1795 = -(t1836 * t1886 + t1839 * t2011) * t2077 - (-t1836 * t1892 + t1839 * t2016) * pkin(5);
t1949 = t1795 * t1997;
t1796 = -(t1837 * t2010 - t1840 * t1888) * t2076 - pkin(5) * (t1837 * t2014 + t1840 * t1894);
t1948 = t1796 * t1992;
t1797 = -(t1837 * t1888 + t1840 * t2010) * t2076 - (-t1837 * t1894 + t1840 * t2014) * pkin(5);
t1947 = t1797 * t1992;
t1798 = -(t1838 * t2009 - t1841 * t1890) * t2075 - pkin(5) * (t1838 * t2012 + t1841 * t1896);
t1946 = t1798 * t1987;
t1799 = -(t1838 * t1890 + t1841 * t2009) * t2075 - (-t1838 * t1896 + t1841 * t2012) * pkin(5);
t1945 = t1799 * t1987;
t1944 = t2094 * t2088;
t1943 = t2093 * t2087;
t1942 = t2092 * t2086;
t1941 = t2057 * t2082;
t1940 = t1848 * t1994;
t1939 = t1850 * t1994;
t1938 = t2007 * t2056;
t1937 = t2018 * t2055;
t1936 = t2053 * t2081;
t1935 = t1849 * t1989;
t1934 = t1851 * t1989;
t1933 = t2005 * t2052;
t1932 = t2018 * t2051;
t1931 = t2049 * t2080;
t1930 = t1852 * t1984;
t1929 = t1853 * t1984;
t1928 = t2003 * t2048;
t1927 = t2018 * t2047;
t1926 = t2089 * t1964;
t1925 = t1828 * t1964;
t1924 = t2090 * t1963;
t1923 = t1830 * t1963;
t1922 = t2091 * t1962;
t1921 = t1832 * t1962;
t1920 = t1788 * t1949;
t1919 = t1789 * t1947;
t1918 = t1790 * t1945;
t1917 = t1791 * t1950;
t1916 = t1792 * t1948;
t1915 = t1793 * t1946;
t1914 = t1879 * t1925;
t1913 = t1879 * t1923;
t1912 = t1879 * t1921;
t1905 = t1937 * t2033;
t1904 = t1932 * t2030;
t1903 = t1927 * t2027;
t1902 = (t1794 * t1809 + t1795 * t1806) * t2057;
t1901 = (t1796 * t1810 + t1797 * t1807) * t2053;
t1900 = (t1798 * t1811 + t1799 * t1808) * t2049;
t1899 = (t1900 * t2080 + t1901 * t2081 + t1902 * t2082) * t2072 + (t1872 * t1902 + t1874 * t1901 + t1876 * t1900) * t2071 + (-t2026 * t2083 - t2029 * t2084 - t2032 * t2085) * t2073 + (t2025 * t2083 + t2028 * t2084 + t2031 * t2085) * t2074 + (t1788 * t1791 * t1828 + t1789 * t1792 * t1830 + t1790 * t1793 * t1832) * MDP(1) + (t1794 * t1795 * t2094 + t1796 * t1797 * t2093 + t1798 * t1799 * t2092) * t2070 + (t1966 * t2086 + t1972 * t2087 + t1978 * t2088) * t2079 + (t1869 * t1944 + t1870 * t1943 + t1871 * t1942) * MDP(5) + (t1942 + t1943 + t1944) * MDP(2);
t1805 = t1811 ^ 2;
t1804 = t1810 ^ 2;
t1803 = t1809 ^ 2;
t1802 = t1808 ^ 2;
t1801 = t1807 ^ 2;
t1800 = t1806 ^ 2;
t1787 = (t1799 * t1983 + t1811 * t1968) * t1876;
t1786 = (t1798 * t1983 + t1808 * t1968) * t1876;
t1785 = (t1797 * t1988 + t1810 * t1974) * t1874;
t1784 = (t1796 * t1988 + t1807 * t1974) * t1874;
t1783 = (t1795 * t1993 + t1809 * t1980) * t1872;
t1782 = (t1794 * t1993 + t1806 * t1980) * t1872;
t1781 = t1787 * t1889;
t1780 = t1787 * t1895;
t1779 = t1786 * t1889;
t1778 = t1786 * t1895;
t1777 = t1785 * t1887;
t1776 = t1785 * t1893;
t1775 = t1784 * t1887;
t1774 = t1784 * t1893;
t1773 = t1783 * t1885;
t1772 = t1783 * t1891;
t1771 = t1782 * t1885;
t1770 = t1782 * t1891;
t1765 = -t1799 * t1927 - t1781;
t1764 = -t1798 * t1927 - t1779;
t1763 = -t1797 * t1932 - t1777;
t1762 = -t1796 * t1932 - t1775;
t1761 = -t1795 * t1937 - t1773;
t1760 = -t1794 * t1937 - t1771;
t1759 = -t1799 * t1903 + t1780;
t1758 = -t1798 * t1903 + t1778;
t1757 = -t1797 * t1904 + t1776;
t1756 = -t1796 * t1904 + t1774;
t1755 = -t1795 * t1905 + t1772;
t1754 = -t1794 * t1905 + t1770;
t1 = [(t1788 ^ 2 * t1828 + t1789 ^ 2 * t1830 + t1790 ^ 2 * t1832) * MDP(1) + (t1800 * t2094 + t1801 * t2093 + t1802 * t2092) * MDP(2) + (t1892 * t1956 + t1894 * t1955 + t1896 * t1954) * t2001 + (-t1886 * t1956 - t1888 * t1955 - t1890 * t1954) * t2000 + (t1800 * t1979 + t1801 * t1973 + t1802 * t1967) * MDP(5) + (t1800 * t1978 + t1801 * t1972 + t1802 * t1966) * t2079 + (t1794 * t1806 * t1941 + t1796 * t1807 * t1936 + t1798 * t1808 * t1931) * t1999 + (t1806 * t1950 + t1807 * t1948 + t1808 * t1946) * t1998 + (t1794 ^ 2 * t2094 + t1796 ^ 2 * t2093 + t1798 ^ 2 * t2092) * t2070 + ((t1808 * t1965 + (t1798 * t1930 + t1758) * t2091) * t1790 + (t1807 * t1971 + (t1796 * t1935 + t1756) * t2090) * t1789 + (t1806 * t1977 + (t1794 * t1940 + t1754) * t2089) * t1788) * MDP(10) + ((-t1808 * t1912 + (-t1798 * t1929 + t1764) * t2091) * t1790 + (-t1807 * t1913 + (-t1796 * t1934 + t1762) * t2090) * t1789 + (-t1806 * t1914 + (-t1794 * t1939 + t1760) * t2089) * t1788) * MDP(11) + MDP(12); (t1755 * t2068 + t1757 * t2066 + t1759 * t2064 + (t1848 * t1917 + t1849 * t1916 + t1852 * t1915) * t1897 + (t2034 * t2059 + t2038 * t2061 + t2042 * t2063) * t1879) * MDP(10) + (t1761 * t2068 + t1763 * t2066 + t1765 * t2064 + (-t1850 * t1917 - t1851 * t1916 - t1853 * t1915) * t1897 + (-t1921 * t2059 - t1923 * t2061 - t1925 * t2063) * t1879) * MDP(11) + t1899; ((t1794 * t1996 + t1796 * t1991 + t1798 * t1986) * MDP(10) + (-t1794 * t1995 - t1796 * t1990 - t1798 * t1985) * MDP(11)) * t1897 + ((t1806 * t2045 + t1807 * t2041 + t1808 * t2037) * MDP(10) + (-t1806 * t1926 - t1807 * t1924 - t1808 * t1922) * MDP(11)) * t1879 + t1961; (t1754 * t2062 + t1756 * t2060 + t1758 * t2058 + (t1848 * t1920 + t1849 * t1919 + t1852 * t1918) * t1897 + (t2034 * t2065 + t2038 * t2067 + t2042 * t2069) * t1879) * MDP(10) + (t1760 * t2062 + t1762 * t2060 + t1764 * t2058 + (-t1850 * t1920 - t1851 * t1919 - t1853 * t1918) * t1897 + (-t1921 * t2065 - t1923 * t2067 - t1925 * t2069) * t1879) * MDP(11) + t1899; (t1791 ^ 2 * t1828 + t1792 ^ 2 * t1830 + t1793 ^ 2 * t1832) * MDP(1) + (t1803 * t2094 + t1804 * t2093 + t1805 * t2092) * MDP(2) + (t1892 * t1953 + t1894 * t1952 + t1896 * t1951) * t2001 + (-t1886 * t1953 - t1888 * t1952 - t1890 * t1951) * t2000 + (t1803 * t1979 + t1804 * t1973 + t1805 * t1967) * MDP(5) + (t1803 * t1978 + t1804 * t1972 + t1805 * t1966) * t2079 + (t1795 * t1809 * t1941 + t1797 * t1810 * t1936 + t1799 * t1811 * t1931) * t1999 + (t1809 * t1949 + t1810 * t1947 + t1811 * t1945) * t1998 + (t1795 ^ 2 * t2094 + t1797 ^ 2 * t2093 + t1799 ^ 2 * t2092) * t2070 + ((t1811 * t1965 + (t1799 * t1930 + t1759) * t2091) * t1793 + (t1810 * t1971 + (t1797 * t1935 + t1757) * t2090) * t1792 + (t1809 * t1977 + (t1795 * t1940 + t1755) * t2089) * t1791) * MDP(10) + ((-t1811 * t1912 + (-t1799 * t1929 + t1765) * t2091) * t1793 + (-t1810 * t1913 + (-t1797 * t1934 + t1763) * t2090) * t1792 + (-t1809 * t1914 + (-t1795 * t1939 + t1761) * t2089) * t1791) * MDP(11) + MDP(12); ((t1795 * t1996 + t1797 * t1991 + t1799 * t1986) * MDP(10) + (-t1795 * t1995 - t1797 * t1990 - t1799 * t1985) * MDP(11)) * t1897 + ((t1809 * t2045 + t1810 * t2041 + t1811 * t2037) * MDP(10) + (-t1809 * t1926 - t1810 * t1924 - t1811 * t1922) * MDP(11)) * t1879 + t1960; (t1770 + t1774 + t1778) * MDP(10) + (-t1771 - t1775 - t1779) * MDP(11) + ((-t1794 * t1938 - t1796 * t1933 - t1798 * t1928) * MDP(10) + (-t1794 * t2055 - t1796 * t2051 - t1798 * t2047) * MDP(11)) * t2018 + t1961; (t1772 + t1776 + t1780) * MDP(10) + (-t1773 - t1777 - t1781) * MDP(11) + ((-t1795 * t1938 - t1797 * t1933 - t1799 * t1928) * MDP(10) + (-t1795 * t2055 - t1797 * t2051 - t1799 * t2047) * MDP(11)) * t2018 + t1960; 0.3e1 * MDP(1) + MDP(12);];
%% Postprocessing: Reshape Output
% From vec2mat_3_matlab.m
res = [t1(1), t1(2), t1(3); t1(4), t1(5), t1(6); t1(7), t1(8), t1(9);];
MMX  = res;
