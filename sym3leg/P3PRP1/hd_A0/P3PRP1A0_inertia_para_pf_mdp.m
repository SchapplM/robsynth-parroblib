% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3PRP1A0
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [11x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRP1A0_convert_par2_MPV_fixb.m

% Output:
% MMX [3x3]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:42
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P3PRP1A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,1),zeros(11,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRP1A0_inertia_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRP1A0_inertia_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PRP1A0_inertia_para_pf_mdp: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRP1A0_inertia_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRP1A0_inertia_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [11 1]), ...
  'P3PRP1A0_inertia_para_pf_mdp: MDP has to be [11x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:42:35
% EndTime: 2019-05-03 14:42:41
% DurationCPUTime: 6.00s
% Computational Cost: add. (11133->482), mult. (22327->799), div. (651->6), fcn. (9146->14), ass. (0->307)
t1993 = (pkin(2) ^ 2);
t2150 = t1993 + 1;
t2159 = MDP(5) + 2 * MDP(3);
t1990 = koppelP(3,1);
t1987 = koppelP(3,2);
t2121 = qJ(3,3) * t1987;
t1948 = pkin(2) * t1990 - t2121;
t2120 = qJ(3,3) * t1990;
t1949 = pkin(2) * t1987 + t2120;
t1974 = legFrame(3,3);
t1960 = sin(t1974);
t1963 = cos(t1974);
t1983 = xP(3);
t1966 = sin(t1983);
t1967 = cos(t1983);
t2158 = t1963 * (t1948 * t1966 + t1949 * t1967) - (t1948 * t1967 - t1949 * t1966) * t1960;
t1991 = koppelP(2,1);
t1988 = koppelP(2,2);
t2127 = qJ(3,2) * t1988;
t1950 = pkin(2) * t1991 - t2127;
t2126 = qJ(3,2) * t1991;
t1951 = pkin(2) * t1988 + t2126;
t1975 = legFrame(2,3);
t1961 = sin(t1975);
t1964 = cos(t1975);
t2157 = t1964 * (t1950 * t1966 + t1951 * t1967) - (t1950 * t1967 - t1951 * t1966) * t1961;
t1992 = koppelP(1,1);
t1989 = koppelP(1,2);
t2133 = qJ(3,1) * t1989;
t1952 = pkin(2) * t1992 - t2133;
t2132 = qJ(3,1) * t1992;
t1953 = pkin(2) * t1989 + t2132;
t1976 = legFrame(1,3);
t1962 = sin(t1976);
t1965 = cos(t1976);
t2156 = (t1952 * t1966 + t1953 * t1967) * t1965 - (t1952 * t1967 - t1953 * t1966) * t1962;
t2154 = 2 * MDP(4);
t2153 = 2 * qJ(3,1);
t2152 = 2 * qJ(3,2);
t2151 = 2 * qJ(3,3);
t1977 = sin(qJ(2,3));
t1980 = cos(qJ(2,3));
t2122 = qJ(3,3) * t1980;
t2049 = -0.2e1 * t2122;
t2124 = qJ(3,3) * t1963;
t1902 = t1960 * t2049 + t1977 * (pkin(2) * t1960 - t2124);
t2149 = pkin(2) * t1902;
t2125 = qJ(3,3) * t1960;
t1903 = t1963 * t2049 + t1977 * (pkin(2) * t1963 + t2125);
t2148 = pkin(2) * t1903;
t1978 = sin(qJ(2,2));
t1981 = cos(qJ(2,2));
t2128 = qJ(3,2) * t1981;
t2051 = -0.2e1 * t2128;
t2130 = qJ(3,2) * t1964;
t1904 = t1961 * t2051 + t1978 * (pkin(2) * t1961 - t2130);
t2147 = pkin(2) * t1904;
t2131 = qJ(3,2) * t1961;
t1905 = t1964 * t2051 + t1978 * (pkin(2) * t1964 + t2131);
t2146 = pkin(2) * t1905;
t1979 = sin(qJ(2,1));
t1982 = cos(qJ(2,1));
t2134 = qJ(3,1) * t1982;
t2053 = -0.2e1 * t2134;
t2136 = qJ(3,1) * t1965;
t1906 = t1962 * t2053 + t1979 * (pkin(2) * t1962 - t2136);
t2145 = pkin(2) * t1906;
t2137 = qJ(3,1) * t1962;
t1907 = t1965 * t2053 + t1979 * (pkin(2) * t1965 + t2137);
t2144 = pkin(2) * t1907;
t1986 = qJ(3,1) ^ 2;
t1956 = -t1986 + t2150;
t1973 = t1982 ^ 2;
t2052 = 0.2e1 * t2134;
t2054 = t1993 + t1986;
t1919 = pkin(2) * t1979 * t2052 + t1956 * t1973 - t2054 - 0.1e1;
t1915 = 0.1e1 / t1919;
t1931 = t1966 * t1992 + t1967 * t1989;
t1934 = -t1966 * t1989 + t1967 * t1992;
t1857 = (t1906 * t1934 - t1907 * t1931) * t1915;
t2143 = t1857 * pkin(2);
t1984 = qJ(3,3) ^ 2;
t1954 = -t1984 + t2150;
t1971 = t1980 ^ 2;
t2048 = 0.2e1 * t2122;
t2056 = t1993 + t1984;
t1917 = pkin(2) * t1977 * t2048 + t1954 * t1971 - t2056 - 0.1e1;
t1911 = 0.1e1 / t1917;
t1929 = t1966 * t1990 + t1967 * t1987;
t1932 = -t1966 * t1987 + t1967 * t1990;
t1858 = (t1902 * t1932 - t1903 * t1929) * t1911;
t2142 = t1858 * pkin(2);
t1985 = qJ(3,2) ^ 2;
t1955 = -t1985 + t2150;
t1972 = t1981 ^ 2;
t2050 = 0.2e1 * t2128;
t2055 = t1993 + t1985;
t1918 = pkin(2) * t1978 * t2050 + t1955 * t1972 - t2055 - 0.1e1;
t1913 = 0.1e1 / t1918;
t1930 = t1966 * t1991 + t1967 * t1988;
t1933 = -t1966 * t1988 + t1967 * t1991;
t1859 = (t1904 * t1933 - t1905 * t1930) * t1913;
t2141 = t1859 * pkin(2);
t2140 = MDP(7) * t1911;
t2139 = MDP(7) * t1913;
t2138 = MDP(7) * t1915;
t2135 = qJ(3,1) * t1979;
t2129 = qJ(3,2) * t1978;
t2123 = qJ(3,3) * t1977;
t1957 = pkin(2) * t2121;
t1923 = t1954 * t1990 - 0.2e1 * t1957;
t2045 = pkin(2) * t2120;
t2059 = t2150 * t1987;
t2065 = t1984 * t1987;
t1924 = 0.2e1 * t2045 + t2059 - t2065;
t1896 = t1923 * t1967 - t1924 * t1966;
t1897 = t1923 * t1966 + t1924 * t1967;
t1942 = t1990 * t1993 - t1957 + t1990;
t1943 = t2045 + t2059;
t2068 = t1977 * t1980;
t1824 = (t1896 * t1963 + t1897 * t1960) * t1971 + (-t1896 * t1960 + t1897 * t1963) * t2068 + (-t1942 * t1967 + t1943 * t1966) * t1963 - (t1942 * t1966 + t1943 * t1967) * t1960;
t2119 = t1824 * t1858;
t2118 = t1824 * t1911;
t1958 = pkin(2) * t2127;
t1925 = t1955 * t1991 - 0.2e1 * t1958;
t2046 = pkin(2) * t2126;
t2058 = t2150 * t1988;
t2064 = t1985 * t1988;
t1926 = 0.2e1 * t2046 + t2058 - t2064;
t1898 = t1925 * t1967 - t1926 * t1966;
t1899 = t1925 * t1966 + t1926 * t1967;
t1944 = t1991 * t1993 - t1958 + t1991;
t1945 = t2046 + t2058;
t2067 = t1978 * t1981;
t1825 = (t1898 * t1964 + t1899 * t1961) * t1972 + (-t1898 * t1961 + t1899 * t1964) * t2067 + (-t1944 * t1967 + t1945 * t1966) * t1964 - (t1944 * t1966 + t1945 * t1967) * t1961;
t2117 = t1825 * t1859;
t2116 = t1825 * t1913;
t1959 = pkin(2) * t2133;
t1927 = t1956 * t1992 - 0.2e1 * t1959;
t2047 = pkin(2) * t2132;
t2057 = t2150 * t1989;
t2063 = t1986 * t1989;
t1928 = 0.2e1 * t2047 + t2057 - t2063;
t1900 = t1927 * t1967 - t1928 * t1966;
t1901 = t1927 * t1966 + t1928 * t1967;
t1946 = t1992 * t1993 - t1959 + t1992;
t1947 = t2047 + t2057;
t2066 = t1979 * t1982;
t1826 = (t1900 * t1965 + t1901 * t1962) * t1973 + (-t1900 * t1962 + t1901 * t1965) * t2066 + (-t1946 * t1967 + t1947 * t1966) * t1965 - (t1946 * t1966 + t1947 * t1967) * t1962;
t2115 = t1826 * t1857;
t2114 = t1826 * t1915;
t2044 = pkin(2) * t2137;
t1922 = t1956 * t1965 + 0.2e1 * t2044;
t2039 = pkin(2) * t2136;
t2005 = -t1956 * t1962 + 0.2e1 * t2039;
t1888 = t1922 * t1973 - t1965 * t1993 + t2005 * t2066 - t1965 - t2044;
t1889 = -t1922 * t2066 + t1993 * t1962 + t2005 * t1973 + t1962 - t2039;
t2072 = t1915 * t1934;
t2073 = t1915 * t1931;
t1845 = t1888 * t2072 - t1889 * t2073;
t1871 = -t2156 * t1979 + (t1931 * t1965 - t1934 * t1962) * t2052;
t2113 = t1845 * t1871;
t2112 = t1845 * t1979;
t2042 = pkin(2) * t2125;
t1920 = t1954 * t1963 + 0.2e1 * t2042;
t2041 = pkin(2) * t2124;
t1999 = -t1954 * t1960 + 0.2e1 * t2041;
t1884 = t1920 * t1971 - t1963 * t1993 + t1999 * t2068 - t1963 - t2042;
t1885 = -t1920 * t2068 + t1993 * t1960 + t1999 * t1971 + t1960 - t2041;
t2082 = t1911 * t1932;
t2083 = t1911 * t1929;
t1846 = t1884 * t2082 - t1885 * t2083;
t1869 = -t2158 * t1977 + (t1929 * t1963 - t1932 * t1960) * t2048;
t2111 = t1846 * t1869;
t2110 = t1846 * t1977;
t2043 = pkin(2) * t2131;
t1921 = t1955 * t1964 + 0.2e1 * t2043;
t2040 = pkin(2) * t2130;
t2002 = -t1955 * t1961 + 0.2e1 * t2040;
t1886 = t1921 * t1972 - t1964 * t1993 + t2002 * t2067 - t1964 - t2043;
t1887 = -t1921 * t2067 + t1993 * t1961 + t2002 * t1972 + t1961 - t2040;
t2077 = t1913 * t1933;
t2078 = t1913 * t1930;
t1847 = t1886 * t2077 - t1887 * t2078;
t1870 = -t2157 * t1978 + (t1930 * t1964 - t1933 * t1961) * t2050;
t2109 = t1847 * t1870;
t2108 = t1847 * t1978;
t2107 = t1857 * t1979;
t2106 = t1858 * t1977;
t2105 = t1859 * t1978;
t2104 = t1869 * t1911;
t2103 = t1870 * t1913;
t2102 = t1871 * t1915;
t2101 = t1884 * t1977;
t1912 = 0.1e1 / t1917 ^ 2;
t2100 = t1885 * t1912;
t2099 = t1885 * t1977;
t2098 = t1886 * t1978;
t1914 = 0.1e1 / t1918 ^ 2;
t2097 = t1887 * t1914;
t2096 = t1887 * t1978;
t2095 = t1888 * t1979;
t1916 = 0.1e1 / t1919 ^ 2;
t2094 = t1889 * t1916;
t2093 = t1889 * t1979;
t2092 = t1902 * t1912;
t2091 = t1903 * t1912;
t2090 = t1904 * t1914;
t2089 = t1905 * t1914;
t2088 = t1906 * t1916;
t2087 = t1907 * t1916;
t2081 = t1911 * t1980;
t2080 = t1912 * t1977;
t2079 = t1912 * t1980;
t2076 = t1913 * t1981;
t2075 = t1914 * t1978;
t2074 = t1914 * t1981;
t2071 = t1915 * t1982;
t2070 = t1916 * t1979;
t2069 = t1916 * t1982;
t1997 = t1963 * t1984 - t2042;
t1998 = -t1960 * t1984 - t2041;
t1890 = t1998 * t1980 + t1977 * (-t1963 - t1997);
t1893 = t1997 * t1980 - t1977 * (t1960 - t1998);
t2062 = t1890 * t2083 - t1893 * t2082;
t2000 = t1964 * t1985 - t2043;
t2001 = -t1961 * t1985 - t2040;
t1891 = t2001 * t1981 + t1978 * (-t1964 - t2000);
t1894 = t2000 * t1981 - t1978 * (t1961 - t2001);
t2061 = t1891 * t2078 - t1894 * t2077;
t2003 = t1965 * t1986 - t2044;
t2004 = -t1962 * t1986 - t2039;
t1892 = t2004 * t1982 + t1979 * (-t1965 - t2003);
t1895 = t2003 * t1982 - t1979 * (t1962 - t2004);
t2060 = t1892 * t2073 - t1895 * t2072;
t2038 = t1869 * t2079;
t2037 = t1870 * t2074;
t2036 = t1871 * t2069;
t2035 = t1902 * t2080;
t2034 = t1902 * t2079;
t2033 = t1903 * t2080;
t2032 = t1903 * t2079;
t2031 = t1904 * t2075;
t2030 = t1904 * t2074;
t2029 = t1905 * t2075;
t2028 = t1905 * t2074;
t2027 = t1906 * t2070;
t2026 = t1906 * t2069;
t2025 = t1907 * t2070;
t2024 = t1907 * t2069;
t2023 = t2071 * t2115 + t2076 * t2117 + t2081 * t2119;
t2022 = t1824 * t2034 + t1825 * t2030 + t1826 * t2026;
t2021 = t1824 * t2032 + t1825 * t2028 + t1826 * t2024;
t2020 = t1884 * t2032 + t1886 * t2028 + t1888 * t2024;
t2017 = t1884 * t1980 - t1893;
t2016 = t1885 * t1980 - t1890;
t2015 = t1886 * t1981 - t1894;
t2014 = t1887 * t1981 - t1891;
t2013 = t1888 * t1982 - t1895;
t2012 = t1889 * t1982 - t1892;
t2011 = t1845 * t1982 + t2060;
t2010 = t1846 * t1980 + t2062;
t2009 = t1847 * t1981 + t2061;
t1813 = (pkin(2) * t2016 + qJ(3,3) * t2099 + t1903 * t2056) * t1911;
t1840 = (t2016 + 0.2e1 * t2148) * t1911;
t1852 = (t1903 * t2151 + t2099) * t1911;
t1996 = MDP(5) * t1840 + MDP(6) * t1852 + MDP(7) * t1813;
t1815 = (pkin(2) * t2014 + qJ(3,2) * t2096 + t1905 * t2055) * t1913;
t1842 = (t2014 + 0.2e1 * t2146) * t1913;
t1854 = (t1905 * t2152 + t2096) * t1913;
t1995 = MDP(5) * t1842 + MDP(6) * t1854 + MDP(7) * t1815;
t1817 = (pkin(2) * t2012 + qJ(3,1) * t2093 + t1907 * t2054) * t1915;
t1844 = (t2012 + 0.2e1 * t2144) * t1915;
t1856 = (t1907 * t2153 + t2093) * t1915;
t1994 = MDP(5) * t1844 + MDP(6) * t1856 + MDP(7) * t1817;
t1941 = t1986 * t1992 + t1959 + t1992;
t1940 = -t1989 + t2047 - t2063;
t1939 = t1985 * t1991 + t1958 + t1991;
t1938 = -t1988 + t2046 - t2064;
t1937 = t1984 * t1990 + t1957 + t1990;
t1936 = -t1987 + t2045 - t2065;
t1935 = (t1966 ^ 2 + t1967 ^ 2) * MDP(11);
t1855 = (t1906 * t2153 + t2095) * t1915;
t1853 = (t1904 * t2152 + t2098) * t1913;
t1851 = (t1902 * t2151 + t2101) * t1911;
t1850 = ((-t1940 * t1967 + t1941 * t1966) * t1965 - (t1940 * t1966 + t1941 * t1967) * t1962) * t1979 + t2156 * t2134;
t1849 = ((-t1938 * t1967 + t1939 * t1966) * t1964 - (t1938 * t1966 + t1939 * t1967) * t1961) * t1978 + t2157 * t2128;
t1848 = ((-t1936 * t1967 + t1937 * t1966) * t1963 - (t1936 * t1966 + t1937 * t1967) * t1960) * t1977 + t2158 * t2122;
t1843 = (t2013 + 0.2e1 * t2145) * t1915;
t1841 = (t2015 + 0.2e1 * t2147) * t1913;
t1839 = (t2017 + 0.2e1 * t2149) * t1911;
t1838 = (-t2012 - t2144) * t1915;
t1837 = (-t2013 - t2145) * t1915;
t1836 = (-t2014 - t2146) * t1913;
t1835 = (-t2015 - t2147) * t1913;
t1834 = (-t2016 - t2148) * t1911;
t1833 = (-t2017 - t2149) * t1911;
t1832 = ((-t1892 + t2144) * t1982 + t1907 * t2135 + t1889) * t1915;
t1831 = ((-t1895 + t2145) * t1982 + t1906 * t2135 + t1888) * t1915;
t1830 = ((-t1891 + t2146) * t1981 + t1905 * t2129 + t1887) * t1913;
t1829 = ((-t1894 + t2147) * t1981 + t1904 * t2129 + t1886) * t1913;
t1828 = ((-t1890 + t2148) * t1980 + t1903 * t2123 + t1885) * t1911;
t1827 = ((-t1893 + t2149) * t1980 + t1902 * t2123 + t1884) * t1911;
t1816 = (pkin(2) * t2013 + qJ(3,1) * t2095 + t1906 * t2054) * t1915;
t1814 = (pkin(2) * t2015 + qJ(3,2) * t2098 + t1904 * t2055) * t1913;
t1812 = (pkin(2) * t2017 + qJ(3,3) * t2101 + t1902 * t2056) * t1911;
t1 = [(t1885 ^ 2 * t1912 + t1887 ^ 2 * t1914 + t1889 ^ 2 * t1916) * MDP(1) + ((t1832 * t1889 + t1838 * t1892) * t1915 + (t1830 * t1887 + t1836 * t1891) * t1913 + (t1828 * t1885 + t1834 * t1890) * t1911) * MDP(7) + t1935 + ((t1844 * t1915 - t1892 * t1916) * MDP(5) + (t1856 * t1915 + t1889 * t2070) * MDP(6) + t1817 * t2138 + MDP(2) * t2087) * t1907 + ((t1842 * t1913 - t1891 * t1914) * MDP(5) + (t1854 * t1913 + t1887 * t2075) * MDP(6) + t1815 * t2139 + MDP(2) * t2089) * t1905 + ((t1840 * t1911 - t1890 * t1912) * MDP(5) + (t1852 * t1911 + t1885 * t2080) * MDP(6) + t1813 * t2140 + MDP(2) * t2091) * t1903 + (-t1885 * t2033 - t1887 * t2029 - t1889 * t2025) * t2154 + t2159 * (t1885 * t2032 + t1887 * t2028 + t1889 * t2024); (t1884 * t2100 + t1886 * t2097 + t1888 * t2094) * MDP(1) + (t1902 * t2091 + t1904 * t2089 + t1906 * t2087) * MDP(2) + (t1885 * t2034 + t1887 * t2030 + t1889 * t2026 + t2020) * MDP(3) + ((-t1888 * t1907 - t1889 * t1906) * t2070 + (-t1886 * t1905 - t1887 * t1904) * t2075 + (-t1884 * t1903 - t1885 * t1902) * t2080) * MDP(4) + (-t1893 * t2091 - t1894 * t2089 - t1895 * t2087 + t2020) * MDP(5) + (t1884 * t2033 + t1886 * t2029 + t1888 * t2025) * MDP(6) + ((t1832 * t1888 + t1838 * t1895) * MDP(7) + t1994 * t1906) * t1915 + ((t1830 * t1886 + t1836 * t1894) * MDP(7) + t1995 * t1904) * t1913 + ((t1828 * t1884 + t1834 * t1893) * MDP(7) + t1996 * t1902) * t1911; (t1884 ^ 2 * t1912 + t1886 ^ 2 * t1914 + t1888 ^ 2 * t1916) * MDP(1) + ((t1831 * t1888 + t1837 * t1895) * t1915 + (t1829 * t1886 + t1835 * t1894) * t1913 + (t1827 * t1884 + t1833 * t1893) * t1911) * MDP(7) + t1935 + ((t1843 * t1915 - t1895 * t1916) * MDP(5) + (t1855 * t1915 + t1888 * t2070) * MDP(6) + t1816 * t2138 + MDP(2) * t2088) * t1906 + ((t1841 * t1913 - t1894 * t1914) * MDP(5) + (t1853 * t1913 + t1886 * t2075) * MDP(6) + t1814 * t2139 + MDP(2) * t2090) * t1904 + ((t1839 * t1911 - t1893 * t1912) * MDP(5) + (t1851 * t1911 + t1884 * t2080) * MDP(6) + t1812 * t2140 + MDP(2) * t2092) * t1902 + (-t1884 * t2035 - t1886 * t2031 - t1888 * t2027) * t2154 + t2159 * (t1884 * t2034 + t1886 * t2030 + t1888 * t2026); (t1824 * t2100 + t1825 * t2097 + t1826 * t2094) * MDP(1) + (t1869 * t2091 + t1870 * t2089 + t1871 * t2087) * MDP(2) + (t1885 * t2038 + t1887 * t2037 + t1889 * t2036 + t2021) * MDP(3) + ((-t1826 * t1907 - t1871 * t1889) * t2070 + (-t1825 * t1905 - t1870 * t1887) * t2075 + (-t1824 * t1903 - t1869 * t1885) * t2080) * MDP(4) + (-t1848 * t2091 - t1849 * t2089 - t1850 * t2087 + t2021) * MDP(5) + (t1824 * t2033 + t1825 * t2029 + t1826 * t2025) * MDP(6) - t1966 * MDP(9) - t1967 * MDP(10) + ((t1826 * t1832 + t1838 * t1850) * MDP(7) + t1994 * t1871) * t1915 + ((t1825 * t1830 + t1836 * t1849) * MDP(7) + t1995 * t1870) * t1913 + ((t1824 * t1828 + t1834 * t1848) * MDP(7) + t1996 * t1869) * t1911; (t1824 * t1884 * t1912 + t1825 * t1886 * t1914 + t1826 * t1888 * t1916) * MDP(1) + (t1869 * t2092 + t1870 * t2090 + t1871 * t2088) * MDP(2) + (t1884 * t2038 + t1886 * t2037 + t1888 * t2036 + t2022) * MDP(3) + ((-t1826 * t1906 - t1871 * t1888) * t2070 + (-t1825 * t1904 - t1870 * t1886) * t2075 + (-t1824 * t1902 - t1869 * t1884) * t2080) * MDP(4) + (-t1848 * t2092 - t1849 * t2090 - t1850 * t2088 + t2022) * MDP(5) + (t1824 * t2035 + t1825 * t2031 + t1826 * t2027) * MDP(6) + t1967 * MDP(9) - t1966 * MDP(10) + ((t1826 * t1831 + t1837 * t1850) * MDP(7) + (MDP(5) * t1843 + MDP(6) * t1855 + MDP(7) * t1816) * t1871) * t1915 + ((t1825 * t1829 + t1835 * t1849) * MDP(7) + (MDP(5) * t1841 + MDP(6) * t1853 + MDP(7) * t1814) * t1870) * t1913 + ((t1824 * t1827 + t1833 * t1848) * MDP(7) + (MDP(5) * t1839 + MDP(6) * t1851 + MDP(7) * t1812) * t1869) * t1911; (t1845 * t2114 + t1846 * t2118 + t1847 * t2116) * MDP(1) + (t1857 * t2102 + t1858 * t2104 + t1859 * t2103) * MDP(2) + (t2071 * t2113 + t2076 * t2109 + t2081 * t2111 + t2023) * MDP(3) + ((-t2113 - t2115) * t1979 * t1915 + (-t2109 - t2117) * t1978 * t1913 + (-t2111 - t2119) * t1977 * t1911) * MDP(4) + ((t1871 * (t2011 + 0.2e1 * t2143) - t1850 * t1857) * t1915 + (t1870 * (t2009 + 0.2e1 * t2141) - t1849 * t1859) * t1913 + (t1869 * (t2010 + 0.2e1 * t2142) - t1848 * t1858) * t1911 + t2023) * MDP(5) + ((t1826 * t2107 + t1871 * (t1857 * t2153 + t2112)) * t1915 + (t1825 * t2105 + t1870 * (t1859 * t2152 + t2108)) * t1913 + (t1824 * t2106 + t1869 * (t1858 * t2151 + t2110)) * t1911) * MDP(6) + (((t2060 + t2143) * t1982 + qJ(3,1) * t2107 + t1845) * t2114 + (pkin(2) * t2011 + qJ(3,1) * t2112 + t1857 * t2054) * t2102 + t1850 * t1915 * (-t2011 - t2143) + ((t2061 + t2141) * t1981 + qJ(3,2) * t2105 + t1847) * t2116 + (pkin(2) * t2009 + qJ(3,2) * t2108 + t1859 * t2055) * t2103 + t1849 * t1913 * (-t2009 - t2141) + ((t2062 + t2142) * t1980 + qJ(3,3) * t2106 + t1846) * t2118 + (pkin(2) * t2010 + qJ(3,3) * t2110 + t1858 * t2056) * t2104 + t1848 * t1911 * (-t2010 - t2142)) * MDP(7) + MDP(8);];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t1(1), t1(2), t1(4); t1(2), t1(3), t1(5); t1(4), t1(5), t1(6);];
MMX  = res;