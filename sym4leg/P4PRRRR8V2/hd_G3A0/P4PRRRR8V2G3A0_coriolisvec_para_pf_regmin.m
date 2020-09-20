% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P4PRRRR8V2G3A0
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
% Datum: 2020-08-07 11:30
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P4PRRRR8V2G3A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G3A0_coriolisvec_para_pf_regmin: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR8V2G3A0_coriolisvec_para_pf_regmin: xDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G3A0_coriolisvec_para_pf_regmin: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G3A0_coriolisvec_para_pf_regmin: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G3A0_coriolisvec_para_pf_regmin: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G3A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:27:53
% EndTime: 2020-08-07 11:28:28
% DurationCPUTime: 37.40s
% Computational Cost: add. (178091->811), mult. (365768->1676), div. (10244->22), fcn. (289976->30), ass. (0->644)
t2039 = xP(4);
t1999 = sin(t2039);
t2000 = cos(t2039);
t2043 = koppelP(1,2);
t2047 = koppelP(1,1);
t1932 = t1999 * t2047 + t2000 * t2043;
t2021 = legFrame(1,2);
t1991 = sin(t2021);
t1995 = cos(t2021);
t2034 = xDP(4);
t2036 = xDP(2);
t2037 = xDP(1);
t2115 = t1999 * t2043 - t2000 * t2047;
t1835 = -(t2034 * t2115 - t2036) * t1991 + (t1932 * t2034 - t2037) * t1995;
t2010 = sin(pkin(8));
t2035 = xDP(3);
t1969 = t2035 * t2010;
t2012 = cos(pkin(8));
t1828 = t1835 * t2012 + t1969;
t2032 = cos(qJ(3,1));
t2458 = pkin(3) * t2032;
t1975 = pkin(2) + t2458;
t2033 = cos(qJ(2,1));
t2038 = pkin(7) + pkin(6);
t1981 = t2033 * t2038;
t2027 = sin(qJ(2,1));
t1926 = t1975 * t2027 - t1981;
t1978 = t2027 * t2038;
t2013 = cos(pkin(4));
t2363 = t2012 * t2035;
t1810 = (t1975 * t2033 + t1978) * t1828 * t2013 - (t1835 * t2010 - t2363) * t1926;
t2026 = sin(qJ(3,1));
t2342 = t2026 * t2013;
t1940 = t1975 * t2342;
t2011 = sin(pkin(4));
t2366 = t2011 * t2032;
t1884 = t1926 * t2366 + t1940;
t2424 = t1810 / t1884;
t2512 = 0.2e1 * t2424;
t2042 = koppelP(2,2);
t2046 = koppelP(2,1);
t1931 = t1999 * t2046 + t2000 * t2042;
t2020 = legFrame(2,2);
t1990 = sin(t2020);
t1994 = cos(t2020);
t2116 = t1999 * t2042 - t2000 * t2046;
t1834 = -(t2116 * t2034 - t2036) * t1990 + (t1931 * t2034 - t2037) * t1994;
t1827 = t1834 * t2012 + t1969;
t2030 = cos(qJ(3,2));
t2459 = pkin(3) * t2030;
t1974 = pkin(2) + t2459;
t2031 = cos(qJ(2,2));
t1980 = t2031 * t2038;
t2025 = sin(qJ(2,2));
t1925 = t1974 * t2025 - t1980;
t1977 = t2025 * t2038;
t1809 = (t1974 * t2031 + t1977) * t1827 * t2013 - (t1834 * t2010 - t2363) * t1925;
t2024 = sin(qJ(3,2));
t2345 = t2024 * t2013;
t1939 = t1974 * t2345;
t2367 = t2011 * t2030;
t1883 = t1925 * t2367 + t1939;
t2425 = t1809 / t1883;
t2511 = 0.2e1 * t2425;
t2041 = koppelP(3,2);
t2045 = koppelP(3,1);
t1930 = t1999 * t2045 + t2000 * t2041;
t2019 = legFrame(3,2);
t1989 = sin(t2019);
t1993 = cos(t2019);
t2117 = t1999 * t2041 - t2000 * t2045;
t1833 = -(t2117 * t2034 - t2036) * t1989 + (t1930 * t2034 - t2037) * t1993;
t1826 = t1833 * t2012 + t1969;
t2028 = cos(qJ(3,3));
t2460 = pkin(3) * t2028;
t1973 = pkin(2) + t2460;
t2029 = cos(qJ(2,3));
t1979 = t2029 * t2038;
t2023 = sin(qJ(2,3));
t1924 = t1973 * t2023 - t1979;
t1976 = t2023 * t2038;
t1808 = (t1973 * t2029 + t1976) * t1826 * t2013 - (t1833 * t2010 - t2363) * t1924;
t2022 = sin(qJ(3,3));
t2348 = t2022 * t2013;
t1938 = t1973 * t2348;
t2368 = t2011 * t2028;
t1882 = t1924 * t2368 + t1938;
t2426 = t1808 / t1882;
t2510 = 0.2e1 * t2426;
t2015 = sin(qJ(2,4));
t2016 = cos(qJ(3,4));
t2350 = t2015 * t2016;
t2017 = cos(qJ(2,4));
t1972 = t2017 * t2038;
t2483 = pkin(2) * t2015 - t1972;
t1907 = pkin(3) * t2350 + t2483;
t2461 = pkin(3) * t2016;
t1970 = pkin(2) + t2461;
t2014 = sin(qJ(3,4));
t2352 = t2014 * t2013;
t1937 = t1970 * t2352;
t2375 = t2011 * t2016;
t1856 = t1907 * t2375 + t1937;
t1853 = 0.1e1 / t1856;
t2346 = t2023 * t2028;
t2482 = pkin(2) * t2023 - t1979;
t1908 = pkin(3) * t2346 + t2482;
t1870 = t1908 * t2368 + t1938;
t1861 = 0.1e1 / t1870;
t2343 = t2025 * t2030;
t2481 = pkin(2) * t2025 - t1980;
t1909 = pkin(3) * t2343 + t2481;
t1871 = t1909 * t2367 + t1939;
t1863 = 0.1e1 / t1871;
t2340 = t2027 * t2032;
t2480 = pkin(2) * t2027 - t1981;
t1910 = pkin(3) * t2340 + t2480;
t1872 = t1910 * t2366 + t1940;
t1865 = 0.1e1 / t1872;
t1898 = pkin(3) * t2342 + t2011 * t2480;
t2008 = t2032 ^ 2;
t2462 = pkin(3) * t2008;
t2318 = t2027 * t2462;
t1848 = 0.1e1 / (pkin(2) * t2342 + t1898 * t2032 + t2011 * t2318);
t1897 = pkin(3) * t2345 + t2011 * t2481;
t2007 = t2030 ^ 2;
t2463 = pkin(3) * t2007;
t2319 = t2025 * t2463;
t1847 = 0.1e1 / (pkin(2) * t2345 + t1897 * t2030 + t2011 * t2319);
t1896 = pkin(3) * t2348 + t2011 * t2482;
t2006 = t2028 ^ 2;
t2464 = pkin(3) * t2006;
t2320 = t2023 * t2464;
t1846 = 0.1e1 / (pkin(2) * t2348 + t1896 * t2028 + t2011 * t2320);
t1895 = pkin(3) * t2352 + t2011 * t2483;
t2002 = t2016 ^ 2;
t2465 = pkin(3) * t2002;
t2321 = t2015 * t2465;
t1845 = 0.1e1 / (pkin(2) * t2352 + t1895 * t2016 + t2011 * t2321);
t2358 = t2013 * t2027;
t1922 = t2010 * t2033 + t2012 * t2358;
t1892 = t1922 * t2026 + t2012 * t2366;
t2509 = t1892 * t2512;
t2359 = t2013 * t2025;
t1921 = t2010 * t2031 + t2012 * t2359;
t1891 = t1921 * t2024 + t2012 * t2367;
t2508 = t1891 * t2511;
t2360 = t2013 * t2023;
t1920 = t2010 * t2029 + t2012 * t2360;
t1890 = t1920 * t2022 + t2012 * t2368;
t2507 = t1890 * t2510;
t2502 = t2011 * (t2032 * t2480 + t2318);
t2501 = t2011 * (t2030 * t2481 + t2319);
t2500 = t2011 * (t2028 * t2482 + t2320);
t2499 = t2011 * (t2016 * t2483 + t2321);
t2354 = t2013 * t2035;
t1800 = ((t1833 * t2360 - t2029 * t2035) * t2012 + (t1833 * t2029 + t2023 * t2354) * t2010) * t2022 + t1826 * t2368;
t2434 = t1800 * t1861;
t2272 = t2038 * t2434;
t2172 = t2022 * t2272;
t1781 = t2172 - t2426;
t2009 = t2034 ^ 2;
t2049 = 0.1e1 / pkin(3);
t2253 = t2049 * t2426;
t2156 = t2029 * t2253;
t2234 = t2011 * t2346;
t2254 = t2022 * t2426;
t2273 = t2029 * t2434;
t2477 = t1930 * t1989 + t1993 * t2117;
t1744 = (-((t2011 * t2273 + t2013 * t2253) * t2464 + ((-t2254 + t2272) * t2023 + pkin(2) * t2273) * t2368 + t1781 * t2013) * t2434 - (t2011 * t2156 + (t2013 * t2006 - t2022 * t2234 - t2013) * t2434) * t2426 - t2477 * t2009 * t1890) * t1846;
t1982 = pkin(2) ^ 2 + t2038 ^ 2;
t2048 = pkin(3) ^ 2;
t2471 = 0.2e1 * pkin(2);
t2453 = pkin(3) * t2471;
t1764 = -t2038 * t2254 + (t2006 * t2048 + t2028 * t2453 + t1982) * t2434;
t1946 = pkin(2) * t2029 + t1976;
t2374 = t2011 * t2022;
t2094 = pkin(3) * t2374 - t2013 * t2482;
t1850 = t1946 * t2012 + t2010 * t2094;
t1914 = t2010 * t2360 - t2012 * t2029;
t2373 = t2011 * t2023;
t2379 = t2010 * t2011;
t2469 = pkin(2) * t2022;
t1813 = (t1914 * t1989 + t1993 * t2373) * t2464 + (-t1850 * t1989 + t1896 * t1993) * t2028 + (-t1989 * t2379 + t1993 * t2013) * t2469;
t1814 = -(t1914 * t1993 - t1989 * t2373) * t2464 + (t1850 * t1993 + t1896 * t1989) * t2028 + (t1989 * t2013 + t1993 * t2379) * t2469;
t2274 = t2028 * t2434;
t1752 = (t1764 * t2274 + (pkin(2) * t2253 - t1781 * t2028) * t2426 + (-t1813 * t1930 + t1814 * t2117) * t2009) * t1846;
t1740 = pkin(6) * t1744 + t1752 * t2373;
t1729 = t1752 * t2013 * t2028 - t1740 * t2022;
t2221 = pkin(6) * t2253;
t1774 = t2434 * t2469 + t2028 * t2221 / 0.2e1;
t2357 = t2013 * t2029;
t2364 = t2012 * t2013;
t1836 = (t2010 * t2023 - t2012 * t2357) * t2460 + t2482 * t2010 - t1946 * t2364;
t2072 = t1729 * t1836 + t1774 * t2507;
t1801 = ((t1834 * t2359 - t2031 * t2035) * t2012 + (t1834 * t2031 + t2025 * t2354) * t2010) * t2024 + t1827 * t2367;
t2433 = t1801 * t1863;
t2269 = t2038 * t2433;
t2171 = t2024 * t2269;
t1782 = t2171 - t2425;
t2250 = t2049 * t2425;
t2155 = t2031 * t2250;
t2233 = t2011 * t2343;
t2251 = t2024 * t2425;
t2270 = t2031 * t2433;
t2478 = t1931 * t1990 + t1994 * t2116;
t1745 = (-((t2011 * t2270 + t2013 * t2250) * t2463 + ((-t2251 + t2269) * t2025 + pkin(2) * t2270) * t2367 + t1782 * t2013) * t2433 - (t2011 * t2155 + (t2013 * t2007 - t2024 * t2233 - t2013) * t2433) * t2425 - t2478 * t2009 * t1891) * t1847;
t1765 = -t2038 * t2251 + (t2007 * t2048 + t2030 * t2453 + t1982) * t2433;
t1947 = pkin(2) * t2031 + t1977;
t2372 = t2011 * t2024;
t2093 = pkin(3) * t2372 - t2013 * t2481;
t1851 = t1947 * t2012 + t2010 * t2093;
t1915 = t2010 * t2359 - t2012 * t2031;
t2371 = t2011 * t2025;
t2468 = pkin(2) * t2024;
t1815 = (t1915 * t1990 + t1994 * t2371) * t2463 + (-t1851 * t1990 + t1897 * t1994) * t2030 + (-t1990 * t2379 + t1994 * t2013) * t2468;
t1816 = -(t1915 * t1994 - t1990 * t2371) * t2463 + (t1851 * t1994 + t1897 * t1990) * t2030 + (t1990 * t2013 + t1994 * t2379) * t2468;
t2271 = t2030 * t2433;
t1753 = (t1765 * t2271 + (pkin(2) * t2250 - t1782 * t2030) * t2425 + (-t1815 * t1931 + t1816 * t2116) * t2009) * t1847;
t1741 = pkin(6) * t1745 + t1753 * t2371;
t1730 = t1753 * t2013 * t2030 - t1741 * t2024;
t2220 = pkin(6) * t2250;
t1775 = t2433 * t2468 + t2030 * t2220 / 0.2e1;
t2356 = t2013 * t2031;
t1838 = (t2010 * t2025 - t2012 * t2356) * t2459 + t2481 * t2010 - t1947 * t2364;
t2071 = t1730 * t1838 + t1775 * t2508;
t1802 = ((t1835 * t2358 - t2033 * t2035) * t2012 + (t1835 * t2033 + t2027 * t2354) * t2010) * t2026 + t1828 * t2366;
t2432 = t1802 * t1865;
t2266 = t2038 * t2432;
t2170 = t2026 * t2266;
t1783 = t2170 - t2424;
t2247 = t2049 * t2424;
t2154 = t2033 * t2247;
t2232 = t2011 * t2340;
t2248 = t2026 * t2424;
t2267 = t2033 * t2432;
t2479 = t1932 * t1991 + t1995 * t2115;
t1746 = (-((t2011 * t2267 + t2013 * t2247) * t2462 + ((-t2248 + t2266) * t2027 + pkin(2) * t2267) * t2366 + t1783 * t2013) * t2432 - (t2011 * t2154 + (t2013 * t2008 - t2026 * t2232 - t2013) * t2432) * t2424 - t2479 * t2009 * t1892) * t1848;
t1766 = -t2038 * t2248 + (t2008 * t2048 + t2032 * t2453 + t1982) * t2432;
t1948 = pkin(2) * t2033 + t1978;
t2370 = t2011 * t2026;
t2092 = pkin(3) * t2370 - t2013 * t2480;
t1852 = t1948 * t2012 + t2010 * t2092;
t1916 = t2010 * t2358 - t2012 * t2033;
t2369 = t2011 * t2027;
t2467 = pkin(2) * t2026;
t1817 = (t1916 * t1991 + t1995 * t2369) * t2462 + (-t1852 * t1991 + t1898 * t1995) * t2032 + (-t1991 * t2379 + t1995 * t2013) * t2467;
t1818 = -(t1916 * t1995 - t1991 * t2369) * t2462 + (t1852 * t1995 + t1898 * t1991) * t2032 + (t1991 * t2013 + t1995 * t2379) * t2467;
t2268 = t2032 * t2432;
t1754 = (t1766 * t2268 + (pkin(2) * t2247 - t1783 * t2032) * t2424 + (-t1817 * t1932 + t1818 * t2115) * t2009) * t1848;
t1742 = pkin(6) * t1746 + t1754 * t2369;
t1731 = t1754 * t2013 * t2032 - t1742 * t2026;
t2219 = pkin(6) * t2247;
t1776 = t2432 * t2467 + t2032 * t2219 / 0.2e1;
t2355 = t2013 * t2033;
t1840 = (t2010 * t2027 - t2012 * t2355) * t2458 + t2480 * t2010 - t1948 * t2364;
t2070 = t1731 * t1840 + t1776 * t2509;
t1732 = -t1740 * t2028 - t1752 * t2348;
t1777 = pkin(2) * t2274 - t2022 * t2221 / 0.2e1;
t2069 = t1732 * t1836 + t1777 * t2507;
t1733 = -t1741 * t2030 - t1753 * t2345;
t1778 = pkin(2) * t2271 - t2024 * t2220 / 0.2e1;
t2068 = t1733 * t1838 + t1778 * t2508;
t1734 = -t1742 * t2032 - t1754 * t2342;
t1779 = pkin(2) * t2268 - t2026 * t2219 / 0.2e1;
t2067 = t1734 * t1840 + t1779 * t2509;
t1792 = t1800 ^ 2 / t1870 ^ 2;
t2493 = 0.2e1 * t2006 - 0.1e1;
t1768 = t2493 * t1792;
t2148 = t2434 * t2510;
t2104 = t2493 * t2148;
t2063 = t1768 * t1836 + t1890 * t2104;
t1793 = t1801 ^ 2 / t1871 ^ 2;
t2492 = 0.2e1 * t2007 - 0.1e1;
t1769 = t2492 * t1793;
t2147 = t2433 * t2511;
t2102 = t2492 * t2147;
t2062 = t1769 * t1838 + t1891 * t2102;
t1794 = t1802 ^ 2 / t1872 ^ 2;
t2491 = 0.2e1 * t2008 - 0.1e1;
t1770 = t2491 * t1794;
t2146 = t2432 * t2512;
t2100 = t2491 * t2146;
t2061 = t1770 * t1840 + t1892 * t2100;
t2411 = t1846 * t2028;
t2240 = t2022 * t2411;
t2484 = (t1836 * t1792 + t1890 * t2148) * t2240;
t2405 = t1847 * t2030;
t2238 = t2024 * t2405;
t2485 = (t1838 * t1793 + t1891 * t2147) * t2238;
t2399 = t1848 * t2032;
t2236 = t2026 * t2399;
t2486 = (t1840 * t1794 + t1892 * t2146) * t2236;
t2494 = 0.2e1 * t2002 - 0.1e1;
t2490 = t1970 * t2013;
t2489 = t1973 * t2013;
t2488 = t1974 * t2013;
t2487 = t1975 * t2013;
t2040 = koppelP(4,2);
t2044 = koppelP(4,1);
t1929 = t1999 * t2044 + t2000 * t2040;
t2018 = legFrame(4,2);
t1988 = sin(t2018);
t1992 = cos(t2018);
t2118 = t1999 * t2040 - t2000 * t2044;
t1832 = -(t2118 * t2034 - t2036) * t1988 + (t1929 * t2034 - t2037) * t1992;
t1825 = t1832 * t2012 + t1969;
t2362 = t2013 * t2015;
t1796 = ((t1832 * t2362 - t2017 * t2035) * t2012 + (t1832 * t2017 + t2015 * t2354) * t2010) * t2014 + t1825 * t2375;
t1791 = t1796 ^ 2 / t1856 ^ 2;
t1971 = t2015 * t2038;
t1942 = pkin(2) * t2017 + t1971;
t2361 = t2013 * t2017;
t1829 = (t2010 * t2015 - t2012 * t2361) * t2461 + t2483 * t2010 - t1942 * t2364;
t1913 = t2010 * t2017 + t2012 * t2362;
t1886 = t1913 * t2014 + t2012 * t2375;
t1923 = t1970 * t2015 - t1972;
t1804 = t1825 * (t1970 * t2017 + t1971) * t2013 - (t1832 * t2010 - t2363) * t1923;
t1875 = t1923 * t2375 + t1937;
t2430 = t1804 / t1875;
t2322 = 0.2e1 * t2430;
t2435 = t1796 * t1853;
t2149 = t2322 * t2435;
t2417 = t1845 * t2016;
t2242 = t2014 * t2417;
t2051 = (t1791 * t1829 + t1886 * t2149) * t2242;
t2476 = t1929 * t1988 + t1992 * t2118;
t2377 = t2011 * t2014;
t2095 = pkin(3) * t2377 - t2013 * t2483;
t1849 = t1942 * t2012 + t2010 * t2095;
t1911 = t2010 * t2362 - t2012 * t2017;
t2376 = t2011 * t2015;
t2470 = pkin(2) * t2014;
t1811 = (-t1911 * t1992 + t1988 * t2376) * t2465 + (t1849 * t1992 + t1895 * t1988) * t2016 + (t1988 * t2013 + t1992 * t2379) * t2470;
t1812 = (t1911 * t1988 + t1992 * t2376) * t2465 + (-t1849 * t1988 + t1895 * t1992) * t2016 + (-t1988 * t2379 + t1992 * t2013) * t2470;
t2466 = pkin(2) * t2049;
t2263 = t2014 * t2430;
t1763 = -t2038 * t2263 + (t2002 * t2048 + t2016 * t2453 + t1982) * t2435;
t2278 = t2038 * t2435;
t2179 = t2014 * t2278;
t2262 = t2049 * t2430;
t2353 = t2013 * t2049;
t2365 = t2011 * t2049;
t2380 = t2009 * t2049;
t1747 = -(-t2013 * t2179 + (-t1907 * t2014 * t2365 + t2013 * (t2016 * t2466 + t2002)) * t2430) * t1853 * t2262 + (-t1763 * t2353 * t2435 + t2476 * t1829 * t2380) * t1845;
t2457 = pkin(6) * t1747;
t2423 = t1836 * t1846;
t1748 = -t1764 * t1846 * t2353 * t2434 - (-t2013 * t2172 + (-t1908 * t2022 * t2365 + t2013 * (t2028 * t2466 + t2006)) * t2426) * t1861 * t2253 + t2477 * t2380 * t2423;
t2456 = pkin(6) * t1748;
t2422 = t1838 * t1847;
t1749 = -t1765 * t1847 * t2353 * t2433 - (-t2013 * t2171 + (-t1909 * t2024 * t2365 + t2013 * (t2030 * t2466 + t2007)) * t2425) * t1863 * t2250 + t2478 * t2380 * t2422;
t2455 = pkin(6) * t1749;
t2421 = t1840 * t1848;
t1750 = -t1766 * t1848 * t2353 * t2432 - (-t2013 * t2170 + (-t1910 * t2026 * t2365 + t2013 * (t2032 * t2466 + t2008)) * t2424) * t1865 * t2247 + t2479 * t2380 * t2421;
t2454 = pkin(6) * t1750;
t1773 = t2179 - t2430;
t2166 = t2017 * t2262;
t2235 = t2011 * t2350;
t2279 = t2017 * t2435;
t1743 = (-((t2011 * t2279 + t2013 * t2262) * t2465 + ((-t2263 + t2278) * t2015 + pkin(2) * t2279) * t2375 + t1773 * t2013) * t2435 - (t2011 * t2166 + (t2013 * t2002 - t2014 * t2235 - t2013) * t2435) * t2430 - t2476 * t2009 * t1886) * t1845;
t2452 = t1743 * t1845;
t2451 = t1744 * t1846;
t2450 = t1745 * t1847;
t2449 = t1746 * t1848;
t2448 = t1747 * t1845;
t2280 = t2016 * t2435;
t1751 = (t1763 * t2280 + (pkin(2) * t2262 - t1773 * t2016) * t2430 + (t1811 * t2118 - t1812 * t1929) * t2009) * t1845;
t2447 = t1751 * t1845;
t2446 = t1751 * t2015;
t2445 = t1751 * t2017;
t2444 = t1752 * t1846;
t2443 = t1752 * t2023;
t2442 = t1752 * t2029;
t2441 = t1753 * t1847;
t2440 = t1753 * t2025;
t2439 = t1753 * t2031;
t2438 = t1754 * t1848;
t2437 = t1754 * t2027;
t2436 = t1754 * t2033;
t2431 = t1804 ^ 2 / t1875 ^ 2;
t2429 = t1808 ^ 2 / t1882 ^ 2;
t2428 = t1809 ^ 2 / t1883 ^ 2;
t2427 = t1810 ^ 2 / t1884 ^ 2;
t2420 = t1845 * t1886;
t2419 = t1845 * t1992;
t2418 = t1845 * t2014;
t1858 = t1930 * t1993 - t2117 * t1989;
t2416 = t1846 * t1858;
t2415 = t1846 * t1890;
t2414 = t1846 * t1989;
t2413 = t1846 * t1993;
t2412 = t1846 * t2022;
t1859 = t1931 * t1994 - t2116 * t1990;
t2410 = t1847 * t1859;
t2409 = t1847 * t1891;
t2408 = t1847 * t1990;
t2407 = t1847 * t1994;
t2406 = t1847 * t2024;
t1860 = t1932 * t1995 - t2115 * t1991;
t2404 = t1848 * t1860;
t2403 = t1848 * t1892;
t2402 = t1848 * t1991;
t2401 = t1848 * t1995;
t2400 = t1848 * t2026;
t2398 = t1886 * t1988;
t2397 = t1886 * t1992;
t2396 = t1890 * t1989;
t2395 = t1890 * t1993;
t2394 = t1891 * t1990;
t2393 = t1891 * t1994;
t2392 = t1892 * t1991;
t2391 = t1892 * t1995;
t2382 = t1999 * t2009;
t2381 = t2000 * t2009;
t2378 = t2010 * t2013;
t2351 = t2014 * t2015;
t2349 = t2017 * t1743;
t2347 = t2022 * t2023;
t2344 = t2024 * t2025;
t2341 = t2026 * t2027;
t2339 = t2029 * t1744;
t2338 = t2031 * t1745;
t2337 = t2033 * t1746;
t1735 = t1747 * t2013 + t2011 * t2349;
t2050 = 0.1e1 / pkin(3) ^ 2;
t2264 = t2050 * t2431;
t1780 = t1791 + t2264;
t2081 = -0.2e1 * t2166 * t2435;
t2167 = t2013 * t2264;
t2336 = -t1735 * t2014 - t1747 * t2235 + (t1780 * t2351 + t2016 * t2081) * t2011 - t2016 * t2167;
t2335 = t1735 * t2016 - t2014 * t2167 + (-t1747 * t2351 - t1780 * t2350 + t2014 * t2081) * t2011;
t1737 = t1748 * t2013 + t2011 * t2339;
t2260 = t2050 * t2429;
t1784 = t1792 + t2260;
t2079 = -0.2e1 * t2156 * t2434;
t2163 = t2013 * t2260;
t2334 = -t1737 * t2022 - t1748 * t2234 + (t1784 * t2347 + t2028 * t2079) * t2011 - t2028 * t2163;
t1738 = t1749 * t2013 + t2011 * t2338;
t2258 = t2050 * t2428;
t1785 = t1793 + t2258;
t2077 = -0.2e1 * t2155 * t2433;
t2160 = t2013 * t2258;
t2333 = -t1738 * t2024 - t1749 * t2233 + (t1785 * t2344 + t2030 * t2077) * t2011 - t2030 * t2160;
t1739 = t1750 * t2013 + t2011 * t2337;
t2256 = t2050 * t2427;
t1786 = t1794 + t2256;
t2075 = -0.2e1 * t2154 * t2432;
t2157 = t2013 * t2256;
t2332 = -t1739 * t2026 - t1750 * t2232 + (t1786 * t2341 + t2032 * t2075) * t2011 - t2032 * t2157;
t2331 = t1737 * t2028 - t2022 * t2163 + (-t1748 * t2347 - t1784 * t2346 + t2022 * t2079) * t2011;
t2330 = t1738 * t2030 - t2024 * t2160 + (-t1749 * t2344 - t1785 * t2343 + t2024 * t2077) * t2011;
t2329 = t1739 * t2032 - t2026 * t2157 + (-t1750 * t2341 - t1786 * t2340 + t2026 * t2075) * t2011;
t2317 = t1743 * t2420;
t2316 = t2014 ^ 2 * t2452;
t2315 = t1743 * t2418;
t2314 = t1743 * t2417;
t2313 = t1744 * t2415;
t2312 = t2022 ^ 2 * t2451;
t2311 = t1744 * t2412;
t2310 = t1744 * t2411;
t2309 = t1745 * t2409;
t2308 = t2024 ^ 2 * t2450;
t2307 = t1745 * t2406;
t2306 = t1745 * t2405;
t2305 = t1746 * t2403;
t2304 = t2026 ^ 2 * t2449;
t2303 = t1746 * t2400;
t2302 = t1746 * t2399;
t2301 = t1829 * t2448;
t2300 = t1747 * t2418;
t2299 = t1747 * t2417;
t2298 = t1748 * t2423;
t2297 = t1748 * t2412;
t2296 = t1748 * t2411;
t2295 = t1749 * t2422;
t2294 = t1749 * t2406;
t2293 = t1749 * t2405;
t2292 = t1750 * t2421;
t2291 = t1750 * t2400;
t2290 = t1750 * t2399;
t2289 = t1886 * t2446;
t2288 = t1886 * t2445;
t2287 = t1890 * t2443;
t2286 = t1890 * t2442;
t2285 = t1891 * t2440;
t2284 = t1891 * t2439;
t2283 = t1892 * t2437;
t2282 = t1892 * t2436;
t2281 = t1845 * t1791;
t2277 = t1846 * t1792;
t2276 = t1847 * t1793;
t2275 = t1848 * t1794;
t2265 = t1845 * t2431;
t2261 = t1846 * t2429;
t2259 = t1847 * t2428;
t2257 = t1848 * t2427;
t2246 = t1836 * t2416;
t2245 = t1838 * t2410;
t2244 = t1840 * t2404;
t1857 = t1929 * t1992 - t2118 * t1988;
t2243 = t1857 * t2420;
t2241 = t1858 * t2415;
t2239 = t1859 * t2409;
t2237 = t1860 * t2403;
t1885 = t1911 * t2014 + t2010 * t2375;
t2231 = -0.2e1 * t1885 * t2430;
t2230 = t1886 * t2322;
t1888 = t1914 * t2022 + t2010 * t2368;
t2229 = -0.2e1 * t1888 * t2426;
t1889 = t1915 * t2024 + t2010 * t2367;
t2228 = -0.2e1 * t1889 * t2425;
t1887 = t1916 * t2026 + t2010 * t2366;
t2227 = -0.2e1 * t1887 * t2424;
t2222 = pkin(6) * t2262;
t2218 = t1829 * t2315;
t2217 = t1829 * t2314;
t2216 = t1886 * t2316;
t2215 = t1743 * t2242;
t2214 = t1836 * t2311;
t2213 = t1836 * t2310;
t2212 = t1744 * t2246;
t2211 = t1890 * t2312;
t2210 = t1744 * t2240;
t2209 = t1838 * t2307;
t2208 = t1838 * t2306;
t2207 = t1745 * t2245;
t2206 = t1891 * t2308;
t2205 = t1745 * t2238;
t2204 = t1840 * t2303;
t2203 = t1840 * t2302;
t2202 = t1746 * t2244;
t2201 = t1892 * t2304;
t2200 = t1746 * t2236;
t2199 = t1886 * t2300;
t2198 = t1886 * t2299;
t2197 = t1890 * t2297;
t2196 = t1890 * t2296;
t2195 = t1891 * t2294;
t2194 = t1891 * t2293;
t2193 = t1892 * t2291;
t2192 = t1892 * t2290;
t2191 = t1751 * t2243;
t2190 = t1752 * t2241;
t2189 = t1753 * t2239;
t2188 = t1754 * t2237;
t2181 = t2014 * t2281;
t2180 = t2016 * t2281;
t2178 = t2022 * t2277;
t2177 = t2028 * t2277;
t2176 = t2024 * t2276;
t2175 = t2030 * t2276;
t2174 = t2026 * t2275;
t2173 = t2032 * t2275;
t2169 = t2014 * t2265;
t2168 = t2016 * t2265;
t2165 = t2022 * t2261;
t2164 = t2028 * t2261;
t2162 = t2024 * t2259;
t2161 = t2030 * t2259;
t2159 = t2026 * t2257;
t2158 = t2032 * t2257;
t2141 = t1886 * t2215;
t2140 = t1890 * t2210;
t2139 = t1891 * t2205;
t2138 = t1892 * t2200;
t2137 = t1829 * t2181;
t2136 = t1829 * t2180;
t2135 = t1836 * t2178;
t2134 = t1836 * t2177;
t2133 = t2246 * t1792;
t2132 = t1838 * t2176;
t2131 = t1838 * t2175;
t2130 = t2245 * t1793;
t2129 = t1840 * t2174;
t2128 = t1840 * t2173;
t2127 = t2244 * t1794;
t2126 = t1886 * t2169;
t2125 = t1886 * t2168;
t2124 = t1890 * t2165;
t2123 = t1890 * t2164;
t2122 = t1891 * t2162;
t2121 = t1891 * t2161;
t2120 = t1892 * t2159;
t2119 = t1892 * t2158;
t2114 = t1970 * t2379;
t2113 = t1973 * t2379;
t2112 = t1974 * t2379;
t2111 = t1975 * t2379;
t2110 = 0.2e1 * t2141;
t2109 = 0.2e1 * t2140;
t2108 = 0.2e1 * t2139;
t2107 = 0.2e1 * t2138;
t2106 = t2494 * t2149;
t2099 = t1743 * t2471 + t2011 * t2445;
t2098 = t1744 * t2471 + t2011 * t2442;
t2097 = t1745 * t2471 + t2011 * t2439;
t2096 = t1746 * t2471 + t2011 * t2436;
t2091 = -t2015 * t1743 - t2017 * t1791;
t2090 = -t2015 * t1791 + t2349;
t2089 = -t2023 * t1744 - t2029 * t1792;
t2088 = -t2023 * t1792 + t2339;
t2087 = -t2025 * t1745 - t2031 * t1793;
t2086 = -t2025 * t1793 + t2338;
t2085 = -t2027 * t1746 - t2033 * t1794;
t2084 = -t2027 * t1794 + t2337;
t1736 = pkin(6) * t1743 + t1751 * t2376;
t1721 = t1751 * t2013 * t2016 - t1736 * t2014;
t1771 = t2435 * t2470 + t2016 * t2222 / 0.2e1;
t2074 = t1721 * t1829 + t1771 * t2230;
t1722 = -t1736 * t2016 - t1751 * t2352;
t1772 = pkin(2) * t2280 - t2014 * t2222 / 0.2e1;
t2073 = t1722 * t1829 + t1772 * t2230;
t2066 = t1845 * t2074;
t2065 = t1845 * t2073;
t1767 = t2494 * t1791;
t2064 = t1767 * t1829 + t1886 * t2106;
t1927 = pkin(2) * t2012 + t2038 * t2378;
t1928 = pkin(2) * t2378 - t2012 * t2038;
t2056 = t1911 * t2465 - (t1927 * t2017 - t1928 * t2015) * t2016;
t2055 = t1914 * t2464 - (t1927 * t2029 - t1928 * t2023) * t2028;
t2054 = t1915 * t2463 - (t1927 * t2031 - t1928 * t2025) * t2030;
t2053 = t1916 * t2462 - (t1927 * t2033 - t1928 * t2027) * t2032;
t2052 = t1845 * t2064;
t1844 = (t2010 * t2355 + t2012 * t2027) * t2458 + t1948 * t2378 + t2012 * t2480;
t1843 = (t2010 * t2356 + t2012 * t2025) * t2459 + t1947 * t2378 + t2012 * t2481;
t1842 = (t2010 * t2357 + t2012 * t2023) * t2460 + t1946 * t2378 + t2012 * t2482;
t1831 = (t2010 * t2361 + t2012 * t2015) * t2461 + t1942 * t2378 + t2012 * t2483;
t1824 = -t1922 * t2462 - t1948 * t2010 * t2032 + (pkin(2) * t2370 + t2032 * t2092) * t2012;
t1823 = -t1921 * t2463 - t1947 * t2010 * t2030 + (pkin(2) * t2372 + t2030 * t2093) * t2012;
t1822 = -t1920 * t2464 - t1946 * t2010 * t2028 + (pkin(2) * t2374 + t2028 * t2094) * t2012;
t1821 = -t1913 * t2465 - t1942 * t2010 * t2016 + (pkin(2) * t2377 + t2016 * t2095) * t2012;
t1790 = ((t1932 * t2053 - t2115 * t2502 + (-t1932 * t2111 - t2115 * t2487) * t2026) * t1995 + (-t1932 * t2502 - t2053 * t2115 + (-t1932 * t2487 + t2111 * t2115) * t2026) * t1991) * t1865;
t1789 = ((t1931 * t2054 - t2116 * t2501 + (-t1931 * t2112 - t2116 * t2488) * t2024) * t1994 + (-t1931 * t2501 - t2054 * t2116 + (-t1931 * t2488 + t2112 * t2116) * t2024) * t1990) * t1863;
t1788 = ((t1930 * t2055 - t2117 * t2500 + (-t1930 * t2113 - t2117 * t2489) * t2022) * t1993 + (-t1930 * t2500 - t2055 * t2117 + (-t1930 * t2489 + t2113 * t2117) * t2022) * t1989) * t1861;
t1787 = ((t2056 * t1929 - t2118 * t2499 + (-t1929 * t2114 - t2118 * t2490) * t2014) * t1992 + (-t1929 * t2499 - t2118 * t2056 + (-t1929 * t2490 + t2114 * t2118) * t2014) * t1988) * t1853;
t1728 = -t2026 * t2454 + t2032 * t2096;
t1727 = -t2026 * t2096 - t2032 * t2454;
t1726 = -t2024 * t2455 + t2030 * t2097;
t1725 = -t2024 * t2097 - t2030 * t2455;
t1724 = -t2022 * t2456 + t2028 * t2098;
t1723 = -t2022 * t2098 - t2028 * t2456;
t1714 = -t2014 * t2457 + t2016 * t2099;
t1713 = -t2014 * t2099 - t2016 * t2457;
t1 = [t1811 * t2447 + t1814 * t2444 + t1816 * t2441 + t1818 * t2438, -t1992 * t2317 - t1993 * t2313 - t1994 * t2309 - t1995 * t2305, ((t1818 * t2084 - t1995 * t2282) * t1848 + (t1816 * t2086 - t1994 * t2284) * t1847 + (t1814 * t2088 - t1993 * t2286) * t1846 + (t1811 * t2090 - t1992 * t2288) * t1845) * t2011, ((t1818 * t2085 + t1995 * t2283) * t1848 + (t1816 * t2087 + t1994 * t2285) * t1847 + (t1814 * t2089 + t1993 * t2287) * t1846 + (t1811 * t2091 + t1992 * t2289) * t1845) * t2011, -t1992 * t2216 - t1993 * t2211 - t1994 * t2206 - t1995 * t2201 + (-t1992 * t2051 - t1993 * t2484 - t1994 * t2485 - t1995 * t2486) * t2049, -0.2e1 * t1992 * t2141 - 0.2e1 * t1993 * t2140 - 0.2e1 * t1994 * t2139 - 0.2e1 * t1995 * t2138 + (-t2061 * t2401 - t2062 * t2407 - t2063 * t2413 - t2064 * t2419) * t2049, -t1992 * t2199 - t1993 * t2197 - t1994 * t2195 - t1995 * t2193 + (-t1992 * t2125 - t1993 * t2123 - t1994 * t2121 - t1995 * t2119) * t2050 + (t1992 * t2218 + t1993 * t2214 + t1994 * t2209 + t1995 * t2204) * t2049, -t1992 * t2198 - t1993 * t2196 - t1994 * t2194 - t1995 * t2192 + (t1992 * t2126 + t1993 * t2124 + t1994 * t2122 + t1995 * t2120) * t2050 + (t1992 * t2217 + t1993 * t2213 + t1994 * t2208 + t1995 * t2203) * t2049, (t1992 * t2301 + t1993 * t2298 + t1994 * t2295 + t1995 * t2292) * t2049, (-t1728 * t2391 + t1818 * t2329) * t1848 + (-t1726 * t2393 + t1816 * t2330) * t1847 + (-t1724 * t2395 + t1814 * t2331) * t1846 + (-t1714 * t2397 + t1811 * t2335) * t1845 + (t2070 * t2401 + t2071 * t2407 + t2072 * t2413 + t2074 * t2419 + (t1992 * t2137 + t1993 * t2135 + t1994 * t2132 + t1995 * t2129) * pkin(2)) * t2049, (-t1727 * t2391 + t1818 * t2332) * t1848 + (-t1725 * t2393 + t1816 * t2333) * t1847 + (-t1723 * t2395 + t1814 * t2334) * t1846 + (-t1713 * t2397 + t1811 * t2336) * t1845 + (t2067 * t2401 + t2068 * t2407 + t2069 * t2413 + t2073 * t2419 + (t1992 * t2136 + t1993 * t2134 + t1994 * t2131 + t1995 * t2128) * pkin(2)) * t2049, 0, -t2381, t2382, 0; t1812 * t2447 + t1813 * t2444 + t1815 * t2441 + t1817 * t2438, t1988 * t2317 + t1989 * t2313 + t1990 * t2309 + t1991 * t2305, ((t1817 * t2084 + t1991 * t2282) * t1848 + (t1815 * t2086 + t1990 * t2284) * t1847 + (t1813 * t2088 + t1989 * t2286) * t1846 + (t1812 * t2090 + t1988 * t2288) * t1845) * t2011, ((t1817 * t2085 - t1991 * t2283) * t1848 + (t1815 * t2087 - t1990 * t2285) * t1847 + (t1813 * t2089 - t1989 * t2287) * t1846 + (t1812 * t2091 - t1988 * t2289) * t1845) * t2011, t1988 * t2216 + t1989 * t2211 + t1990 * t2206 + t1991 * t2201 + (t1988 * t2051 + t1989 * t2484 + t1990 * t2485 + t1991 * t2486) * t2049, t1988 * t2110 + t1989 * t2109 + t1990 * t2108 + t1991 * t2107 + (t1988 * t2052 + t2061 * t2402 + t2062 * t2408 + t2063 * t2414) * t2049, t1988 * t2199 + t1989 * t2197 + t1990 * t2195 + t1991 * t2193 + (t1988 * t2125 + t1989 * t2123 + t1990 * t2121 + t1991 * t2119) * t2050 + (-t1988 * t2218 - t1989 * t2214 - t1990 * t2209 - t1991 * t2204) * t2049, t1988 * t2198 + t1989 * t2196 + t1990 * t2194 + t1991 * t2192 + (-t1988 * t2126 - t1989 * t2124 - t1990 * t2122 - t1991 * t2120) * t2050 + (-t1988 * t2217 - t1989 * t2213 - t1990 * t2208 - t1991 * t2203) * t2049, (-t1988 * t2301 - t1989 * t2298 - t1990 * t2295 - t1991 * t2292) * t2049, (t1728 * t2392 + t1817 * t2329) * t1848 + (t1726 * t2394 + t1815 * t2330) * t1847 + (t1724 * t2396 + t1813 * t2331) * t1846 + (t1714 * t2398 + t1812 * t2335) * t1845 + (-t2070 * t2402 - t2071 * t2408 - t2072 * t2414 - t1988 * t2066 + (-t1988 * t2137 - t1989 * t2135 - t1990 * t2132 - t1991 * t2129) * pkin(2)) * t2049, (t1727 * t2392 + t1817 * t2332) * t1848 + (t1725 * t2394 + t1815 * t2333) * t1847 + (t1723 * t2396 + t1813 * t2334) * t1846 + (t1713 * t2398 + t1812 * t2336) * t1845 + (-t2067 * t2402 - t2068 * t2408 - t2069 * t2414 - t1988 * t2065 + (-t1988 * t2136 - t1989 * t2134 - t1990 * t2131 - t1991 * t2128) * pkin(2)) * t2049, 0, -t2382, -t2381, 0; t1821 * t2447 + t1822 * t2444 + t1823 * t2441 + t1824 * t2438, t1885 * t2452 + t1887 * t2449 + t1888 * t2451 + t1889 * t2450, ((t1824 * t2084 + t1887 * t2436) * t1848 + (t1823 * t2086 + t1889 * t2439) * t1847 + (t1822 * t2088 + t1888 * t2442) * t1846 + (t1821 * t2090 + t1885 * t2445) * t1845) * t2011, ((t1824 * t2085 - t1887 * t2437) * t1848 + (t1823 * t2087 - t1889 * t2440) * t1847 + (t1822 * t2089 - t1888 * t2443) * t1846 + (t1821 * t2091 - t1885 * t2446) * t1845) * t2011, t1885 * t2316 + t1888 * t2312 + t1889 * t2308 + t1887 * t2304 + ((-t1794 * t1844 + t1887 * t2146) * t2236 + (-t1793 * t1843 + t1889 * t2147) * t2238 + (-t1792 * t1842 + t1888 * t2148) * t2240 + (-t1791 * t1831 + t1885 * t2149) * t2242) * t2049, 0.2e1 * t1885 * t2215 + 0.2e1 * t1888 * t2210 + 0.2e1 * t1889 * t2205 + 0.2e1 * t1887 * t2200 + ((-t1844 * t1770 + t1887 * t2100) * t1848 + (-t1843 * t1769 + t1889 * t2102) * t1847 + (-t1842 * t1768 + t1888 * t2104) * t1846 + (-t1831 * t1767 + t1885 * t2106) * t1845) * t2049, t1885 * t2300 + t1888 * t2297 + t1889 * t2294 + t1887 * t2291 + (t1885 * t2168 + t1887 * t2158 + t1888 * t2164 + t1889 * t2161) * t2050 + (t1831 * t2315 + t1842 * t2311 + t1843 * t2307 + t1844 * t2303) * t2049, t1885 * t2299 + t1888 * t2296 + t1889 * t2293 + t1887 * t2290 + (-t1885 * t2169 - t1887 * t2159 - t1888 * t2165 - t1889 * t2162) * t2050 + (t1831 * t2314 + t1842 * t2310 + t1843 * t2306 + t1844 * t2302) * t2049, (t1748 * t1842 * t1846 + t1749 * t1843 * t1847 + t1750 * t1844 * t1848 + t1831 * t2448) * t2049, (t1887 * t1728 + t1824 * t2329) * t1848 + (t1889 * t1726 + t1823 * t2330) * t1847 + (t1888 * t1724 + t1822 * t2331) * t1846 + (t1885 * t1714 + t1821 * t2335) * t1845 + ((t1844 * t1731 + t1776 * t2227) * t1848 + (t1843 * t1730 + t1775 * t2228) * t1847 + (t1842 * t1729 + t1774 * t2229) * t1846 + (t1831 * t1721 + t1771 * t2231) * t1845 + (t1831 * t2181 + t1842 * t2178 + t1843 * t2176 + t1844 * t2174) * pkin(2)) * t2049, (t1887 * t1727 + t1824 * t2332) * t1848 + (t1889 * t1725 + t1823 * t2333) * t1847 + (t1888 * t1723 + t1822 * t2334) * t1846 + (t1885 * t1713 + t1821 * t2336) * t1845 + ((t1844 * t1734 + t1779 * t2227) * t1848 + (t1843 * t1733 + t1778 * t2228) * t1847 + (t1842 * t1732 + t1777 * t2229) * t1846 + (t1831 * t1722 + t1772 * t2231) * t1845 + (t1831 * t2180 + t1842 * t2177 + t1843 * t2175 + t1844 * t2173) * pkin(2)) * t2049, 0, 0, 0, 0; t1751 * t1787 + t1752 * t1788 + t1753 * t1789 + t1754 * t1790, t1743 * t2243 + t1744 * t2241 + t1745 * t2239 + t1746 * t2237, (t1787 * t2090 + t1788 * t2088 + t1789 * t2086 + t1790 * t2084 + t2017 * t2191 + t2029 * t2190 + t2031 * t2189 + t2033 * t2188) * t2011, (t1787 * t2091 + t1788 * t2089 + t1789 * t2087 + t1790 * t2085 - t2015 * t2191 - t2023 * t2190 - t2025 * t2189 - t2027 * t2188) * t2011, t1857 * t2216 + t1858 * t2211 + t1859 * t2206 + t1860 * t2201 + (t1857 * t2051 + t1858 * t2484 + t1859 * t2485 + t1860 * t2486) * t2049, t1857 * t2110 + t1858 * t2109 + t1859 * t2108 + t1860 * t2107 + (t1857 * t2052 + t2061 * t2404 + t2062 * t2410 + t2063 * t2416) * t2049, t1857 * t2199 + t1858 * t2197 + t1859 * t2195 + t1860 * t2193 + (t1857 * t2125 + t1858 * t2123 + t1859 * t2121 + t1860 * t2119) * t2050 + (-t1857 * t2218 - t2022 * t2212 - t2024 * t2207 - t2026 * t2202) * t2049, t1857 * t2198 + t1858 * t2196 + t1859 * t2194 + t1860 * t2192 + (-t1857 * t2126 - t1858 * t2124 - t1859 * t2122 - t1860 * t2120) * t2050 + (-t1857 * t2217 - t2028 * t2212 - t2030 * t2207 - t2032 * t2202) * t2049, (-t1748 * t2246 - t1749 * t2245 - t1750 * t2244 - t1857 * t2301) * t2049, t1714 * t2243 + t1724 * t2241 + t1726 * t2239 + t1728 * t2237 + t2329 * t1790 + t2330 * t1789 + t2331 * t1788 + t2335 * t1787 + (-t2070 * t2404 - t2071 * t2410 - t2072 * t2416 - t1857 * t2066 + (-t1857 * t2137 - t2022 * t2133 - t2024 * t2130 - t2026 * t2127) * pkin(2)) * t2049, t1713 * t2243 + t1723 * t2241 + t1725 * t2239 + t1727 * t2237 + t2332 * t1790 + t2333 * t1789 + t2334 * t1788 + t2336 * t1787 + (-t2067 * t2404 - t2068 * t2410 - t2069 * t2416 - t1857 * t2065 + (-t1857 * t2136 - t2028 * t2133 - t2030 * t2130 - t2032 * t2127) * pkin(2)) * t2049, 0, 0, 0, 0;];
tau_reg  = t1;
