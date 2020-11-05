% Calculate inertia matrix for parallel robot
% P4RRRRR10V1G1A0
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
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d1,d2,d4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% MX [4x4]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 13:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P4RRRRR10V1G1A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4RRRRR10V1G1A0_inertia_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4RRRRR10V1G1A0_inertia_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P4RRRRR10V1G1A0_inertia_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4RRRRR10V1G1A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4RRRRR10V1G1A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P4RRRRR10V1G1A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4RRRRR10V1G1A0_inertia_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4RRRRR10V1G1A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:37:40
% EndTime: 2020-08-07 11:38:02
% DurationCPUTime: 22.75s
% Computational Cost: add. (38142->848), mult. (87951->1347), div. (1920->17), fcn. (63432->44), ass. (0->544)
t2369 = m(3) * rSges(3,1);
t2240 = (pkin(5) * t2369);
t2430 = -2 * t2240;
t2040 = cos(qJ(3,1));
t2368 = m(3) * rSges(3,2);
t1954 = rSges(3,1) * t2368 - Icges(3,4);
t2031 = sin(qJ(3,1));
t2271 = t1954 * t2031;
t2065 = rSges(3,2) ^ 2;
t2067 = rSges(3,1) ^ 2;
t1975 = -t2065 + t2067;
t1907 = m(3) * t1975 - Icges(3,1) + Icges(3,2);
t2010 = t2040 ^ 2;
t2414 = t1907 * t2010;
t2429 = -0.2e1 * t2040 * t2271 + t2414;
t2037 = cos(qJ(3,2));
t2028 = sin(qJ(3,2));
t2272 = t1954 * t2028;
t2007 = t2037 ^ 2;
t2415 = t1907 * t2007;
t2428 = -0.2e1 * t2037 * t2272 + t2415;
t2034 = cos(qJ(3,3));
t2025 = sin(qJ(3,3));
t2273 = t1954 * t2025;
t2004 = t2034 ^ 2;
t2416 = t1907 * t2004;
t2427 = -0.2e1 * t2034 * t2273 + t2416;
t2022 = cos(qJ(3,4));
t2019 = sin(qJ(3,4));
t2274 = t1954 * t2019;
t2001 = t2022 ^ 2;
t2417 = t1907 * t2001;
t2426 = -0.2e1 * t2022 * t2274 + t2417;
t2023 = cos(qJ(2,4));
t2003 = t2023 ^ 2;
t2064 = rSges(3,3) ^ 2;
t2066 = rSges(2,2) ^ 2;
t2068 = rSges(2,1) ^ 2;
t2334 = Icges(2,1) + Icges(3,2);
t2413 = Icges(2,2) + Icges(3,3) - t2334 + (-t2066 + t2068) * m(2) + (-t2064 + t2065) * m(3);
t2425 = -(t2413 + t2426) * t2003 + t2019 * t2430;
t2035 = cos(qJ(2,3));
t2006 = t2035 ^ 2;
t2424 = -(t2413 + t2427) * t2006 + t2025 * t2430;
t2038 = cos(qJ(2,2));
t2009 = t2038 ^ 2;
t2423 = -(t2413 + t2428) * t2009 + t2028 * t2430;
t2041 = cos(qJ(2,1));
t2012 = t2041 ^ 2;
t2422 = -(t2413 + t2429) * t2012 + t2031 * t2430;
t2014 = cos(pkin(3));
t2000 = t2014 ^ 2;
t2342 = pkin(6) * t2000;
t2013 = sin(pkin(3));
t2266 = t2013 * t2014;
t2367 = rSges(3,3) * m(3);
t1952 = rSges(3,2) * t2367 - Icges(3,6);
t1953 = rSges(3,1) * t2367 - Icges(3,5);
t1853 = -t1952 * t2022 - t1953 * t2019;
t1857 = -t1952 * t2034 - t1953 * t2025;
t1858 = -t1952 * t2037 - t1953 * t2028;
t1859 = -t1952 * t2040 - t1953 * t2031;
t2032 = sin(qJ(2,1));
t2244 = t2032 * t2041;
t2161 = t2031 * t2244;
t2347 = pkin(2) * t2040;
t2226 = t2013 * t2347;
t2261 = t2013 * t2031;
t2355 = pkin(2) * t2010;
t1985 = t2032 * pkin(6);
t2399 = t1985 + pkin(1);
t2412 = t2014 * (t2040 * t2399 + t2041 * t2355) - t2161 * t2226 + pkin(6) * (t2041 - 0.1e1) * (t2041 + 0.1e1) * t2261;
t2029 = sin(qJ(2,2));
t2246 = t2029 * t2038;
t2162 = t2028 * t2246;
t2348 = pkin(2) * t2037;
t2227 = t2013 * t2348;
t2262 = t2013 * t2028;
t2356 = pkin(2) * t2007;
t1984 = t2029 * pkin(6);
t2400 = t1984 + pkin(1);
t2411 = t2014 * (t2037 * t2400 + t2038 * t2356) - t2162 * t2227 + pkin(6) * (t2038 - 0.1e1) * (t2038 + 0.1e1) * t2262;
t2026 = sin(qJ(2,3));
t2248 = t2026 * t2035;
t2163 = t2025 * t2248;
t2349 = pkin(2) * t2034;
t2228 = t2013 * t2349;
t2263 = t2013 * t2025;
t2357 = pkin(2) * t2004;
t1983 = t2026 * pkin(6);
t2401 = t1983 + pkin(1);
t2410 = (t2034 * t2401 + t2035 * t2357) * t2014 - t2163 * t2228 + pkin(6) * (t2035 - 0.1e1) * (t2035 + 0.1e1) * t2263;
t2020 = sin(qJ(2,4));
t2251 = t2020 * t2023;
t2164 = t2019 * t2251;
t2353 = pkin(2) * t2022;
t2229 = t2013 * t2353;
t2265 = t2013 * t2019;
t2358 = pkin(2) * t2001;
t1972 = t2020 * pkin(6);
t2402 = t1972 + pkin(1);
t2409 = t2014 * (t2022 * t2402 + t2023 * t2358) - t2164 * t2229 + pkin(6) * (t2023 - 0.1e1) * (t2023 + 0.1e1) * t2265;
t2043 = rSges(2,3) + pkin(5);
t2364 = m(2) * t2043;
t2083 = -rSges(2,1) * t2364 + Icges(2,5) + t1954;
t2130 = -rSges(2,2) * t2364 + pkin(5) * t2367 + Icges(2,6);
t2234 = t2019 * t2368;
t2371 = -0.2e1 * t1954;
t2408 = t2023 * (t2130 - t1853) + ((-t1907 * t2019 - t2240) * t2022 + t2001 * t2371 + pkin(5) * t2234 + t2083) * t2020;
t2233 = t2025 * t2368;
t2407 = t2035 * (t2130 - t1857) + ((-t1907 * t2025 - t2240) * t2034 + t2004 * t2371 + pkin(5) * t2233 + t2083) * t2026;
t2232 = t2028 * t2368;
t2406 = t2038 * (t2130 - t1858) + ((-t1907 * t2028 - t2240) * t2037 + t2007 * t2371 + pkin(5) * t2232 + t2083) * t2029;
t2231 = t2031 * t2368;
t2405 = t2041 * (t2130 - t1859) + ((-t1907 * t2031 - t2240) * t2040 + t2010 * t2371 + pkin(5) * t2231 + t2083) * t2032;
t2404 = -0.2e1 * pkin(6);
t2403 = 0.2e1 * pkin(6);
t2350 = pkin(2) * t2031;
t1960 = pkin(1) * t2350;
t2241 = t2040 * t2041;
t2217 = t2032 * t2347;
t1893 = -t2041 * pkin(6) + t2217;
t2288 = t1893 * t2014;
t1832 = 0.1e1 / (pkin(1) * t2288 + (-t1960 + (pkin(2) * t2241 + t1985) * pkin(5)) * t2013);
t2351 = pkin(2) * t2028;
t1959 = pkin(1) * t2351;
t2242 = t2037 * t2038;
t2218 = t2029 * t2348;
t1892 = -t2038 * pkin(6) + t2218;
t2289 = t1892 * t2014;
t1831 = 0.1e1 / (pkin(1) * t2289 + (-t1959 + (pkin(2) * t2242 + t1984) * pkin(5)) * t2013);
t2352 = pkin(2) * t2025;
t1958 = pkin(1) * t2352;
t2243 = t2034 * t2035;
t2219 = t2026 * t2349;
t1891 = -t2035 * pkin(6) + t2219;
t2290 = t1891 * t2014;
t1830 = 0.1e1 / (pkin(1) * t2290 + (-t1958 + (pkin(2) * t2243 + t1983) * pkin(5)) * t2013);
t2354 = pkin(2) * t2019;
t1956 = pkin(1) * t2354;
t2250 = t2022 * t2023;
t2221 = t2020 * t2353;
t1884 = -t2023 * pkin(6) + t2221;
t2291 = t1884 * t2014;
t1828 = 0.1e1 / (pkin(1) * t2291 + (-t1956 + (pkin(2) * t2250 + t1972) * pkin(5)) * t2013);
t2362 = pkin(1) * t2020;
t2361 = pkin(1) * t2026;
t2360 = pkin(1) * t2029;
t2359 = pkin(1) * t2032;
t2370 = m(2) * rSges(2,2);
t2332 = rSges(2,1) * t2370 - Icges(2,4);
t2018 = legFrame(1,3);
t1964 = sin(t2018);
t1968 = cos(t2018);
t2033 = sin(qJ(1,1));
t2042 = cos(qJ(1,1));
t1868 = t1964 * t2033 - t1968 * t2042;
t2069 = pkin(6) ^ 2;
t2071 = pkin(2) ^ 2;
t2144 = t2010 * t2071 - t2069;
t2398 = t2144 * t1868;
t2017 = legFrame(2,3);
t1963 = sin(t2017);
t1967 = cos(t2017);
t2030 = sin(qJ(1,2));
t2039 = cos(qJ(1,2));
t1867 = t1963 * t2030 - t1967 * t2039;
t2145 = t2007 * t2071 - t2069;
t2397 = t2145 * t1867;
t2016 = legFrame(3,3);
t1962 = sin(t2016);
t1966 = cos(t2016);
t2027 = sin(qJ(1,3));
t2036 = cos(qJ(1,3));
t1866 = t1962 * t2027 - t1966 * t2036;
t2146 = t2004 * t2071 - t2069;
t2396 = t2146 * t1866;
t2015 = legFrame(4,3);
t1961 = sin(t2015);
t1965 = cos(t2015);
t2021 = sin(qJ(1,4));
t2024 = cos(qJ(1,4));
t1861 = t1961 * t2021 - t1965 * t2024;
t2147 = t2001 * t2071 - t2069;
t2395 = t2147 * t1861;
t1900 = (t2003 - 0.2e1) * t2354 - pkin(5);
t1929 = pkin(5) * t2019 + pkin(2);
t2156 = t1929 - 0.2e1 * t2358;
t2377 = (-pkin(6) * t2164 - t1900 * t2022) * t2266 + t2250 * t2342 + (t2156 * t2000 - t1929 + t2358) * t2020;
t1901 = (t2006 - 0.2e1) * t2352 - pkin(5);
t1933 = pkin(5) * t2025 + pkin(2);
t2155 = t1933 - 0.2e1 * t2357;
t2376 = (-pkin(6) * t2163 - t1901 * t2034) * t2266 + t2243 * t2342 + (t2155 * t2000 - t1933 + t2357) * t2026;
t1902 = (t2009 - 0.2e1) * t2351 - pkin(5);
t1936 = pkin(5) * t2028 + pkin(2);
t2154 = t1936 - 0.2e1 * t2356;
t2375 = (-pkin(6) * t2162 - t1902 * t2037) * t2266 + t2242 * t2342 + (t2154 * t2000 - t1936 + t2356) * t2029;
t1903 = (t2012 - 0.2e1) * t2350 - pkin(5);
t1939 = pkin(5) * t2031 + pkin(2);
t2153 = t1939 - 0.2e1 * t2355;
t2374 = (-pkin(6) * t2161 - t1903 * t2040) * t2266 + t2241 * t2342 + (t2153 * t2000 - t1939 + t2355) * t2032;
t2373 = -0.2e1 * pkin(1);
t2372 = (t2000 - 0.1e1) * t2403;
t2046 = pkin(1) * rSges(3,1);
t2366 = pkin(1) * rSges(3,2);
t2365 = t1907 / 0.2e1;
t2363 = pkin(1) * t2014;
t2346 = pkin(5) * t2023;
t2345 = pkin(5) * t2035;
t2344 = pkin(5) * t2038;
t2343 = pkin(5) * t2041;
t2340 = (t2066 + t2068) * m(2);
t2333 = Icges(3,1) + Icges(2,3);
t1862 = t1961 * t2024 + t1965 * t2021;
t2184 = t1862 * t2265;
t1772 = t2377 * t1861 + t2409 * t1862 - t2184 * t2362;
t2002 = 0.1e1 / t2022;
t2331 = t1772 * t2002;
t2185 = t1861 * t2265;
t1773 = -t2409 * t1861 + t2377 * t1862 + t2185 * t2362;
t2330 = t1773 * t2002;
t1869 = t1962 * t2036 + t1966 * t2027;
t2180 = t1869 * t2263;
t1775 = t2376 * t1866 + t2410 * t1869 - t2180 * t2361;
t2005 = 0.1e1 / t2034;
t2329 = t1775 * t2005;
t1870 = t1963 * t2039 + t1967 * t2030;
t2179 = t1870 * t2262;
t1776 = t2375 * t1867 + t2411 * t1870 - t2179 * t2360;
t2008 = 0.1e1 / t2037;
t2328 = t1776 * t2008;
t1871 = t1964 * t2042 + t1968 * t2033;
t2178 = t1871 * t2261;
t1777 = t2374 * t1868 + t2412 * t1871 - t2178 * t2359;
t2011 = 0.1e1 / t2040;
t2327 = t1777 * t2011;
t2183 = t1866 * t2263;
t1778 = -t2410 * t1866 + t2376 * t1869 + t2183 * t2361;
t2326 = t1778 * t2005;
t2182 = t1867 * t2262;
t1779 = -t2411 * t1867 + t2375 * t1870 + t2182 * t2360;
t2325 = t1779 * t2008;
t2181 = t1868 * t2261;
t1780 = -t2412 * t1868 + t2374 * t1871 + t2181 * t2359;
t2324 = t1780 * t2011;
t2047 = m(2) * rSges(2,1);
t1917 = -t2047 + t2234;
t1928 = -t2367 + t2370;
t2135 = (t2064 + t2065) * m(3) + t2340 + t2333;
t1792 = t2408 * t2013 + ((-t1928 * t2020 + (t2022 * t2369 - t1917) * t2023) * pkin(1) + t2135 + t2426) * t2014;
t2323 = t1792 * t2002;
t1918 = -t2047 + t2233;
t1793 = t2407 * t2013 + ((-t1928 * t2026 + (t2034 * t2369 - t1918) * t2035) * pkin(1) + t2135 + t2427) * t2014;
t2322 = t1793 * t2005;
t1919 = -t2047 + t2232;
t1794 = t2406 * t2013 + ((-t1928 * t2029 + (t2037 * t2369 - t1919) * t2038) * pkin(1) + t2135 + t2428) * t2014;
t2321 = t1794 * t2008;
t1920 = -t2047 + t2231;
t1795 = t2405 * t2013 + ((-t1928 * t2032 + (t2040 * t2369 - t1920) * t2041) * pkin(1) + t2135 + t2429) * t2014;
t2320 = t1795 * t2011;
t1930 = pkin(5) + t2354;
t1932 = pkin(6) + t2362;
t2264 = t2013 * t2023;
t2168 = t2014 * t2264;
t2169 = (t2014 + 0.1e1) * (t2014 - 0.1e1) * t2071;
t2256 = t2014 * t2020;
t2286 = t1930 * t2013;
t1808 = -t2001 * t2169 * t2251 + (t1930 * t2168 + t2003 * t2372 + t1932 - t2342) * t2353 - ((-t2000 * t1972 + t2402) * t2023 - t2256 * t2286) * pkin(6);
t2072 = 0.1e1 / pkin(2);
t2319 = t1808 * t2072;
t1934 = pkin(5) + t2352;
t1942 = pkin(6) + t2361;
t2260 = t2013 * t2035;
t2167 = t2014 * t2260;
t2255 = t2014 * t2026;
t2285 = t1934 * t2013;
t1810 = -t2004 * t2169 * t2248 + (t1934 * t2167 + t2006 * t2372 + t1942 - t2342) * t2349 - ((-t2000 * t1983 + t2401) * t2035 - t2255 * t2285) * pkin(6);
t2318 = t1810 * t2072;
t1937 = pkin(5) + t2351;
t1943 = pkin(6) + t2360;
t2259 = t2013 * t2038;
t2166 = t2014 * t2259;
t2254 = t2014 * t2029;
t2284 = t1937 * t2013;
t1811 = -t2007 * t2169 * t2246 + (t1937 * t2166 + t2009 * t2372 + t1943 - t2342) * t2348 - ((-t2000 * t1984 + t2400) * t2038 - t2254 * t2284) * pkin(6);
t2317 = t1811 * t2072;
t1940 = pkin(5) + t2350;
t1944 = pkin(6) + t2359;
t2258 = t2013 * t2041;
t2165 = t2014 * t2258;
t2253 = t2014 * t2032;
t2283 = t1940 * t2013;
t1812 = -t2010 * t2169 * t2244 + (t1940 * t2165 + t2012 * t2372 + t1944 - t2342) * t2347 - ((-t2000 * t1985 + t2399) * t2041 - t2253 * t2283) * pkin(6);
t2316 = t1812 * t2072;
t2051 = 0.2e1 * qJ(3,4);
t2237 = t2065 + t2067;
t2125 = Icges(2,3) + t2340 + (0.2e1 * t2064 + t2237) * m(3) / 0.2e1 + Icges(3,2) / 0.2e1 + Icges(3,1) / 0.2e1;
t1816 = cos(t2051) * t2365 - t1954 * sin(t2051) + t2125;
t2315 = t1816 * t2002;
t2052 = 0.2e1 * qJ(3,3);
t1817 = cos(t2052) * t2365 - t1954 * sin(t2052) + t2125;
t2314 = t1817 * t2005;
t2053 = 0.2e1 * qJ(3,2);
t1818 = cos(t2053) * t2365 - t1954 * sin(t2053) + t2125;
t2313 = t1818 * t2008;
t2054 = 0.2e1 * qJ(3,1);
t1819 = cos(t2054) * t2365 - t1954 * sin(t2054) + t2125;
t2312 = t1819 * t2011;
t1955 = pkin(6) * t2363;
t2270 = t2013 * (-pkin(5) * t1972 + t1956);
t2311 = 0.1e1 / ((-pkin(5) * t2229 + t1955) * t2023 - t2221 * t2363 + t2270) * t2002;
t2310 = 0.1e1 / ((pkin(1) * t2256 + pkin(5) * t2264) * t2353 - t2023 * t1955 - t2270) * t2002;
t2269 = t2013 * (-pkin(5) * t1983 + t1958);
t2309 = 0.1e1 / ((-pkin(5) * t2228 + t1955) * t2035 - t2219 * t2363 + t2269) * t2005;
t2268 = t2013 * (-pkin(5) * t1984 + t1959);
t2308 = 0.1e1 / ((-pkin(5) * t2227 + t1955) * t2038 - t2218 * t2363 + t2268) * t2008;
t2267 = t2013 * (-pkin(5) * t1985 + t1960);
t2307 = 0.1e1 / ((-pkin(5) * t2226 + t1955) * t2041 - t2217 * t2363 + t2267) * t2011;
t2306 = 0.1e1 / ((pkin(1) * t2255 + pkin(5) * t2260) * t2349 - t2035 * t1955 - t2269) * t2005;
t2305 = 0.1e1 / ((pkin(1) * t2254 + pkin(5) * t2259) * t2348 - t2038 * t1955 - t2268) * t2008;
t2304 = 0.1e1 / ((pkin(1) * t2253 + pkin(5) * t2258) * t2347 - t2041 * t1955 - t2267) * t2011;
t2303 = t1828 * (t1884 * t2013 + t2014 * t2354);
t2302 = t1830 * (t1891 * t2013 + t2014 * t2352);
t2301 = t1831 * (t1892 * t2013 + t2014 * t2351);
t2300 = t1832 * (t1893 * t2013 + t2014 * t2350);
t2295 = t1853 * t2002;
t2294 = t1857 * t2005;
t2293 = t1858 * t2008;
t2292 = t1859 * t2011;
t1921 = t2237 * m(3) + Icges(3,3);
t2287 = t1921 * t2072;
t2282 = t1952 * t2019;
t2281 = t1952 * t2025;
t2280 = t1952 * t2028;
t2279 = t1952 * t2031;
t2278 = t1953 * t2020;
t2277 = t1953 * t2026;
t2276 = t1953 * t2029;
t2275 = t1953 * t2032;
t2257 = t2013 * t2072;
t2252 = t2019 * t2023;
t2249 = t2025 * t2035;
t2247 = t2028 * t2038;
t2245 = t2031 * t2041;
t2239 = pkin(5) * t2368;
t2238 = 0.2e1 * t2332;
t2070 = pkin(5) ^ 2;
t2236 = pkin(1) ^ 2 + t2070;
t2235 = t1928 * t2373;
t2230 = -0.2e1 * t2266;
t2225 = t2014 * t2353;
t2224 = t2014 * t2349;
t2223 = t2014 * t2348;
t2222 = t2014 * t2347;
t2176 = t2024 * t2286;
t1836 = t1932 * t2021 - t2020 * t2176;
t1931 = 0.2e1 * t1972 + pkin(1);
t1848 = t1931 * t2021 - t2176;
t2177 = t2021 * t2286;
t2091 = t1932 * t2024 + t2020 * t2177;
t2095 = t1931 * t2024 + t2177;
t2124 = t2147 * t1862;
t2142 = t1862 * t2225;
t1784 = (t2142 * t2404 + t2395) * t2003 + ((t1961 * t1848 - t2095 * t1965) * t2353 + t2124 * t2256) * t2023 + (t1961 * t1836 - t2091 * t1965 + t2142) * pkin(6);
t2201 = t1784 * t2311;
t2143 = t1861 * t2225;
t1785 = (t2143 * t2403 + t2124) * t2003 + ((t1848 * t1965 + t2095 * t1961) * t2353 - t2256 * t2395) * t2023 + (t1836 * t1965 + t1961 * t2091 - t2143) * pkin(6);
t2200 = t1785 * t2311;
t2174 = t2036 * t2285;
t1838 = t1942 * t2027 - t2026 * t2174;
t1935 = 0.2e1 * t1983 + pkin(1);
t1849 = t1935 * t2027 - t2174;
t2175 = t2027 * t2285;
t2090 = t1942 * t2036 + t2026 * t2175;
t2094 = t1935 * t2036 + t2175;
t2123 = t2146 * t1869;
t2138 = t1869 * t2224;
t1786 = (t2138 * t2404 + t2396) * t2006 + ((t1962 * t1849 - t2094 * t1966) * t2349 + t2123 * t2255) * t2035 + (t1962 * t1838 - t2090 * t1966 + t2138) * pkin(6);
t2199 = t1786 * t2309;
t2172 = t2039 * t2284;
t1839 = t1943 * t2030 - t2029 * t2172;
t1938 = 0.2e1 * t1984 + pkin(1);
t1850 = t1938 * t2030 - t2172;
t2173 = t2030 * t2284;
t2089 = t1943 * t2039 + t2029 * t2173;
t2093 = t1938 * t2039 + t2173;
t2122 = t2145 * t1870;
t2137 = t1870 * t2223;
t1787 = (t2137 * t2404 + t2397) * t2009 + ((t1963 * t1850 - t2093 * t1967) * t2348 + t2122 * t2254) * t2038 + (t1963 * t1839 - t2089 * t1967 + t2137) * pkin(6);
t2198 = t1787 * t2308;
t2170 = t2042 * t2283;
t1840 = t1944 * t2033 - t2032 * t2170;
t1941 = 0.2e1 * t1985 + pkin(1);
t1851 = t1941 * t2033 - t2170;
t2171 = t2033 * t2283;
t2088 = t1944 * t2042 + t2032 * t2171;
t2092 = t1941 * t2042 + t2171;
t2121 = t2144 * t1871;
t2136 = t1871 * t2222;
t1788 = (t2136 * t2404 + t2398) * t2012 + ((t1964 * t1851 - t2092 * t1968) * t2347 + t2121 * t2253) * t2041 + (t1964 * t1840 - t2088 * t1968 + t2136) * pkin(6);
t2197 = t1788 * t2307;
t2141 = t1866 * t2224;
t1789 = (t2141 * t2403 + t2123) * t2006 + ((t1849 * t1966 + t2094 * t1962) * t2349 - t2255 * t2396) * t2035 + (t1838 * t1966 + t1962 * t2090 - t2141) * pkin(6);
t2196 = t1789 * t2309;
t2140 = t1867 * t2223;
t1790 = (t2140 * t2403 + t2122) * t2009 + ((t1850 * t1967 + t2093 * t1963) * t2348 - t2254 * t2397) * t2038 + (t1839 * t1967 + t1963 * t2089 - t2140) * pkin(6);
t2195 = t1790 * t2308;
t2139 = t1868 * t2222;
t1791 = (t2139 * t2403 + t2121) * t2012 + ((t1851 * t1968 + t2092 * t1964) * t2347 - t2253 * t2398) * t2041 + (t1840 * t1968 + t1964 * t2088 - t2139) * pkin(6);
t2194 = t1791 * t2307;
t2193 = t1808 * t2310;
t1809 = (pkin(6) * t2168 + t1900 * t2000 + pkin(5) + (-t2003 + 0.1e1) * t2354) * t2022 - t2402 * t2252 + (t2156 * t2266 + t2252 * t2342) * t2020;
t2192 = t1809 * t2310;
t2191 = t1810 * t2306;
t2190 = t1811 * t2305;
t2189 = t1812 * t2304;
t1813 = (pkin(6) * t2167 + t1901 * t2000 + pkin(5) + (-t2006 + 0.1e1) * t2352) * t2034 - t2401 * t2249 + (t2155 * t2266 + t2249 * t2342) * t2026;
t2188 = t1813 * t2306;
t1814 = (pkin(6) * t2166 + t1902 * t2000 + pkin(5) + (-t2009 + 0.1e1) * t2351) * t2037 - t2400 * t2247 + (t2154 * t2266 + t2247 * t2342) * t2029;
t2187 = t1814 * t2305;
t1815 = (pkin(6) * t2165 + t1903 * t2000 + pkin(5) + (-t2012 + 0.1e1) * t2350) * t2040 - t2399 * t2245 + (t2153 * t2266 + t2245 * t2342) * t2032;
t2186 = t1815 * t2304;
t2134 = t2257 * t2311;
t2133 = t2257 * t2309;
t2132 = t2257 * t2308;
t2131 = t2257 * t2307;
t1804 = ((-t2278 - m(3) * (rSges(3,2) * t2346 + t2046)) * t2022 + (t2020 * t1952 - m(3) * (rSges(3,1) * t2346 - t2366)) * t2019 - t2023 * t1921) * t2013 - t2014 * (-Icges(3,5) * t2019 - t2022 * Icges(3,6) + (rSges(3,1) * t2019 + rSges(3,2) * t2022) * m(3) * (rSges(3,3) + t2362));
t2120 = t1804 * t2134;
t1805 = ((-t2277 - m(3) * (rSges(3,2) * t2345 + t2046)) * t2034 + (t2026 * t1952 - m(3) * (rSges(3,1) * t2345 - t2366)) * t2025 - t2035 * t1921) * t2013 - t2014 * (-Icges(3,5) * t2025 - t2034 * Icges(3,6) + (rSges(3,1) * t2025 + rSges(3,2) * t2034) * m(3) * (rSges(3,3) + t2361));
t2119 = t1805 * t2133;
t1806 = ((-t2276 - m(3) * (rSges(3,2) * t2344 + t2046)) * t2037 + (t2029 * t1952 - m(3) * (rSges(3,1) * t2344 - t2366)) * t2028 - t2038 * t1921) * t2013 - t2014 * (-Icges(3,5) * t2028 - t2037 * Icges(3,6) + (rSges(3,1) * t2028 + rSges(3,2) * t2037) * m(3) * (rSges(3,3) + t2360));
t2118 = t1806 * t2132;
t1807 = ((-t2275 - m(3) * (rSges(3,2) * t2343 + t2046)) * t2040 + (t2032 * t1952 - m(3) * (rSges(3,1) * t2343 - t2366)) * t2031 - t2041 * t1921) * t2013 - t2014 * (-Icges(3,5) * t2031 - t2040 * Icges(3,6) + (rSges(3,1) * t2031 + rSges(3,2) * t2040) * m(3) * (rSges(3,3) + t2359));
t2117 = t1807 * t2131;
t2116 = t1853 * t2134;
t2115 = t1921 * t2134;
t2114 = t1857 * t2133;
t2113 = t1921 * t2133;
t2112 = t1858 * t2132;
t2111 = t1921 * t2132;
t2110 = t1859 * t2131;
t2109 = t1921 * t2131;
t2100 = Icges(1,3) + (t2066 + (0.2e1 * pkin(5) + rSges(2,3)) * rSges(2,3) + t2236) * m(2) + (t2064 + t2067 + t2236) * m(3) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + t2334;
t2082 = -(rSges(2,1) + t2043) * (-rSges(2,1) + t2043) * m(2) + (-t2070 - t1975) * m(3) + t2333 - t2334;
t2063 = koppelP(1,1);
t2062 = koppelP(2,1);
t2061 = koppelP(3,1);
t2060 = koppelP(4,1);
t2059 = koppelP(1,2);
t2058 = koppelP(2,2);
t2057 = koppelP(3,2);
t2056 = koppelP(4,2);
t2050 = rSges(4,1);
t2049 = rSges(4,2);
t2048 = xP(4);
t1996 = cos(t2048);
t1995 = sin(t2048);
t1982 = 0.2e1 * m(3) * t2046;
t1981 = -0.2e1 * t2239;
t1980 = 0.2e1 * t2239;
t1896 = 0.2e1 * t1907;
t1881 = -t1995 * t2059 + t1996 * t2063;
t1880 = -t1995 * t2058 + t1996 * t2062;
t1879 = -t1995 * t2057 + t1996 * t2061;
t1878 = -t1995 * t2056 + t1996 * t2060;
t1877 = -t1995 * t2063 - t1996 * t2059;
t1876 = -t1995 * t2062 - t1996 * t2058;
t1875 = -t1995 * t2061 - t1996 * t2057;
t1874 = -t1995 * t2060 - t1996 * t2056;
t1873 = m(4) * (-t1995 * t2049 + t1996 * t2050);
t1872 = m(4) * (-t1995 * t2050 - t1996 * t2049);
t1803 = t1868 * t1985 + t1871 * t2288 + (t1868 * t2241 - t2178) * pkin(2);
t1802 = t1867 * t1984 + t1870 * t2289 + (t1867 * t2242 - t2179) * pkin(2);
t1801 = t1866 * t1983 + t1869 * t2290 + (t1866 * t2243 - t2180) * pkin(2);
t1800 = -t1871 * t1985 + t1868 * t2288 + (-t1871 * t2241 - t2181) * pkin(2);
t1799 = -t1870 * t1984 + t1867 * t2289 + (-t1870 * t2242 - t2182) * pkin(2);
t1798 = -t1869 * t1983 + t1866 * t2290 + (-t1869 * t2243 - t2183) * pkin(2);
t1797 = t1861 * t1972 + t1862 * t2291 + (t1861 * t2250 - t2184) * pkin(2);
t1796 = -t1862 * t1972 + t1861 * t2291 + (-t1862 * t2250 - t2185) * pkin(2);
t1783 = (t1800 * t1881 + t1803 * t1877) * t1832;
t1782 = (t1799 * t1880 + t1802 * t1876) * t1831;
t1781 = (t1798 * t1879 + t1801 * t1875) * t1830;
t1774 = (t1796 * t1878 + t1797 * t1874) * t1828;
t1771 = (0.2e1 * (-t1953 * t2040 + t2279 + t2332) * t2244 + t1896 * t2010 + (t1981 - 0.4e1 * t2271) * t2040 + t2082 + t2422) * t2000 - t2405 * t2230 + ((t1982 + 0.2e1 * t2275) * t2040 + (-t2238 - 0.2e1 * t2279) * t2032 + t1920 * t2373) * t2041 - t2414 + (t1980 + 0.2e1 * t2271) * t2040 + t2032 * t2235 + t2100 - t2422;
t1770 = (0.2e1 * (-t1953 * t2037 + t2280 + t2332) * t2246 + t1896 * t2007 + (t1981 - 0.4e1 * t2272) * t2037 + t2082 + t2423) * t2000 - t2406 * t2230 + ((t1982 + 0.2e1 * t2276) * t2037 + (-t2238 - 0.2e1 * t2280) * t2029 + t1919 * t2373) * t2038 - t2415 + (t1980 + 0.2e1 * t2272) * t2037 + t2029 * t2235 + t2100 - t2423;
t1769 = (0.2e1 * (-t1953 * t2034 + t2281 + t2332) * t2248 + t1896 * t2004 + (t1981 - 0.4e1 * t2273) * t2034 + t2082 + t2424) * t2000 - t2407 * t2230 + ((t1982 + 0.2e1 * t2277) * t2034 + (-t2238 - 0.2e1 * t2281) * t2026 + t1918 * t2373) * t2035 - t2416 + (t1980 + 0.2e1 * t2273) * t2034 + t2026 * t2235 + t2100 - t2424;
t1768 = (0.2e1 * (-t1953 * t2022 + t2282 + t2332) * t2251 + t1896 * t2001 + (t1981 - 0.4e1 * t2274) * t2022 + t2082 + t2425) * t2000 - t2408 * t2230 + ((t1982 + 0.2e1 * t2278) * t2022 + (-t2238 - 0.2e1 * t2282) * t2020 + t1917 * t2373) * t2023 - t2417 + (t1980 + 0.2e1 * t2274) * t2022 + t2020 * t2235 + t2100 - t2425;
t1767 = t1807 * t2300 + (t1812 * t2287 + t1815 * t1859) * t2304;
t1766 = t1806 * t2301 + (t1811 * t2287 + t1814 * t1858) * t2305;
t1765 = t1805 * t2302 + (t1810 * t2287 + t1813 * t1857) * t2306;
t1764 = t1804 * t2303 + (t1808 * t2287 + t1809 * t1853) * t2310;
t1763 = t1795 * t2300 + (t1815 * t1819 + t1859 * t2316) * t2304;
t1762 = t1794 * t2301 + (t1814 * t1818 + t1858 * t2317) * t2305;
t1761 = t1793 * t2302 + (t1813 * t1817 + t1857 * t2318) * t2306;
t1760 = t1792 * t2303 + (t1809 * t1816 + t1853 * t2319) * t2310;
t1759 = (t1788 * t1877 - t1791 * t1881) * t2131;
t1758 = (t1787 * t1876 - t1790 * t1880) * t2132;
t1757 = (t1786 * t1875 - t1789 * t1879) * t2133;
t1756 = (t1784 * t1874 - t1785 * t1878) * t2134;
t1755 = (t1777 * t1881 + t1780 * t1877) * t2011 * t1832;
t1754 = (t1776 * t1880 + t1779 * t1876) * t2008 * t1831;
t1753 = (t1775 * t1879 + t1778 * t1875) * t2005 * t1830;
t1752 = (t1772 * t1878 + t1773 * t1874) * t2002 * t1828;
t1751 = t1788 * t2109 - (t1780 * t2292 + t1803 * t1807) * t1832;
t1750 = t1787 * t2111 - (t1779 * t2293 + t1802 * t1806) * t1831;
t1749 = t1786 * t2113 - (t1778 * t2294 + t1801 * t1805) * t1830;
t1748 = -t1791 * t2109 - (t1777 * t2292 + t1800 * t1807) * t1832;
t1747 = -t1790 * t2111 - (t1776 * t2293 + t1799 * t1806) * t1831;
t1746 = -t1789 * t2113 - (t1775 * t2294 + t1798 * t1805) * t1830;
t1745 = t1784 * t2115 - (t1773 * t2295 + t1797 * t1804) * t1828;
t1744 = -t1785 * t2115 - (t1772 * t2295 + t1796 * t1804) * t1828;
t1743 = t1771 * t2300 + (t1795 * t1815 + t1807 * t2316) * t2304;
t1742 = t1770 * t2301 + (t1794 * t1814 + t1806 * t2317) * t2305;
t1741 = t1769 * t2302 + (t1793 * t1813 + t1805 * t2318) * t2306;
t1740 = t1768 * t2303 + (t1792 * t1809 + t1804 * t2319) * t2310;
t1739 = t1788 * t2110 - (t1780 * t2312 + t1795 * t1803) * t1832;
t1738 = t1787 * t2112 - (t1779 * t2313 + t1794 * t1802) * t1831;
t1737 = t1786 * t2114 - (t1778 * t2314 + t1793 * t1801) * t1830;
t1736 = -t1791 * t2110 - (t1777 * t2312 + t1795 * t1800) * t1832;
t1735 = -t1790 * t2112 - (t1776 * t2313 + t1794 * t1799) * t1831;
t1734 = -t1789 * t2114 - (t1775 * t2314 + t1793 * t1798) * t1830;
t1733 = t1784 * t2116 - (t1773 * t2315 + t1792 * t1797) * t1828;
t1732 = -t1785 * t2116 - (t1772 * t2315 + t1792 * t1796) * t1828;
t1731 = t1788 * t2117 - (t1771 * t1803 + t1780 * t2320) * t1832;
t1730 = t1787 * t2118 - (t1770 * t1802 + t1779 * t2321) * t1831;
t1729 = t1786 * t2119 - (t1769 * t1801 + t1778 * t2322) * t1830;
t1728 = -t1791 * t2117 - (t1771 * t1800 + t1777 * t2320) * t1832;
t1727 = -t1790 * t2118 - (t1770 * t1799 + t1776 * t2321) * t1831;
t1726 = -t1789 * t2119 - (t1769 * t1798 + t1775 * t2322) * t1830;
t1725 = t1784 * t2120 - (t1768 * t1797 + t1773 * t2323) * t1828;
t1724 = -t1785 * t2120 - (t1768 * t1796 + t1772 * t2323) * t1828;
t1723 = -t1755 * t1859 + t1759 * t1921 - t1783 * t1807;
t1722 = -t1754 * t1858 + t1758 * t1921 - t1782 * t1806;
t1721 = -t1753 * t1857 + t1757 * t1921 - t1781 * t1805;
t1720 = -t1752 * t1853 + t1756 * t1921 - t1774 * t1804;
t1719 = -t1755 * t1819 + t1759 * t1859 - t1783 * t1795;
t1718 = -t1754 * t1818 + t1758 * t1858 - t1782 * t1794;
t1717 = -t1753 * t1817 + t1757 * t1857 - t1781 * t1793;
t1716 = -t1752 * t1816 + t1756 * t1853 - t1774 * t1792;
t1715 = -t1755 * t1795 + t1759 * t1807 - t1771 * t1783;
t1714 = -t1754 * t1794 + t1758 * t1806 - t1770 * t1782;
t1713 = -t1753 * t1793 + t1757 * t1805 - t1769 * t1781;
t1712 = -t1752 * t1792 + t1756 * t1804 - t1768 * t1774;
t1 = [m(4) - (t1731 * t1803 + t1739 * t2324) * t1832 - (t1730 * t1802 + t1738 * t2325) * t1831 - (t1729 * t1801 + t1737 * t2326) * t1830 - (t1725 * t1797 + t1733 * t2330) * t1828 + (t1745 * t2201 + t1749 * t2199 + t1750 * t2198 + t1751 * t2197) * t2257, -(t1731 * t1800 + t1739 * t2327) * t1832 - (t1730 * t1799 + t1738 * t2328) * t1831 - (t1729 * t1798 + t1737 * t2329) * t1830 - (t1725 * t1796 + t1733 * t2331) * t1828 + (-t1745 * t2200 - t1749 * t2196 - t1750 * t2195 - t1751 * t2194) * t2257, t1733 * t2192 + t1737 * t2188 + t1738 * t2187 + t1739 * t2186 + t1725 * t2303 + t1729 * t2302 + t1730 * t2301 + t1731 * t2300 + (t1745 * t2193 + t1749 * t2191 + t1750 * t2190 + t1751 * t2189) * t2072, -t1725 * t1774 - t1729 * t1781 - t1730 * t1782 - t1731 * t1783 - t1733 * t1752 - t1737 * t1753 - t1738 * t1754 - t1739 * t1755 + t1745 * t1756 + t1749 * t1757 + t1750 * t1758 + t1751 * t1759 + t1872; -(t1728 * t1803 + t1736 * t2324) * t1832 - (t1727 * t1802 + t1735 * t2325) * t1831 - (t1726 * t1801 + t1734 * t2326) * t1830 - (t1724 * t1797 + t1732 * t2330) * t1828 + (t1744 * t2201 + t1746 * t2199 + t1747 * t2198 + t1748 * t2197) * t2257, m(4) - (t1728 * t1800 + t1736 * t2327) * t1832 - (t1727 * t1799 + t1735 * t2328) * t1831 - (t1726 * t1798 + t1734 * t2329) * t1830 - (t1724 * t1796 + t1732 * t2331) * t1828 + (-t1744 * t2200 - t1746 * t2196 - t1747 * t2195 - t1748 * t2194) * t2257, t1732 * t2192 + t1734 * t2188 + t1735 * t2187 + t1736 * t2186 + t1724 * t2303 + t1726 * t2302 + t1727 * t2301 + t1728 * t2300 + (t1744 * t2193 + t1746 * t2191 + t1747 * t2190 + t1748 * t2189) * t2072, -t1724 * t1774 - t1726 * t1781 - t1727 * t1782 - t1728 * t1783 - t1732 * t1752 - t1734 * t1753 - t1735 * t1754 - t1736 * t1755 + t1744 * t1756 + t1746 * t1757 + t1747 * t1758 + t1748 * t1759 + t1873; -(t1743 * t1803 + t1763 * t2324) * t1832 - (t1742 * t1802 + t1762 * t2325) * t1831 - (t1741 * t1801 + t1761 * t2326) * t1830 - (t1740 * t1797 + t1760 * t2330) * t1828 + (t1764 * t2201 + t1765 * t2199 + t1766 * t2198 + t1767 * t2197) * t2257, -(t1743 * t1800 + t1763 * t2327) * t1832 - (t1742 * t1799 + t1762 * t2328) * t1831 - (t1741 * t1798 + t1761 * t2329) * t1830 - (t1740 * t1796 + t1760 * t2331) * t1828 + (-t1764 * t2200 - t1765 * t2196 - t1766 * t2195 - t1767 * t2194) * t2257, t1760 * t2192 + t1761 * t2188 + t1762 * t2187 + t1763 * t2186 + t1740 * t2303 + t1741 * t2302 + t1742 * t2301 + t1743 * t2300 + m(4) + (t1764 * t2193 + t1765 * t2191 + t1766 * t2190 + t1767 * t2189) * t2072, -t1740 * t1774 - t1741 * t1781 - t1742 * t1782 - t1743 * t1783 - t1752 * t1760 - t1753 * t1761 - t1754 * t1762 - t1755 * t1763 + t1756 * t1764 + t1757 * t1765 + t1758 * t1766 + t1759 * t1767; t1872 - (t1715 * t1803 + t1719 * t2324) * t1832 - (t1714 * t1802 + t1718 * t2325) * t1831 - (t1713 * t1801 + t1717 * t2326) * t1830 - (t1712 * t1797 + t1716 * t2330) * t1828 + (t1720 * t2201 + t1721 * t2199 + t1722 * t2198 + t1723 * t2197) * t2257, t1873 - (t1715 * t1800 + t1719 * t2327) * t1832 - (t1714 * t1799 + t1718 * t2328) * t1831 - (t1713 * t1798 + t1717 * t2329) * t1830 - (t1712 * t1796 + t1716 * t2331) * t1828 + (-t1720 * t2200 - t1721 * t2196 - t1722 * t2195 - t1723 * t2194) * t2257, t1716 * t2192 + t1717 * t2188 + t1718 * t2187 + t1719 * t2186 + t1712 * t2303 + t1713 * t2302 + t1714 * t2301 + t1715 * t2300 + (t1720 * t2193 + t1721 * t2191 + t1722 * t2190 + t1723 * t2189) * t2072, -t1715 * t1783 - t1719 * t1755 + t1723 * t1759 - t1714 * t1782 - t1718 * t1754 + t1722 * t1758 - t1713 * t1781 - t1717 * t1753 + t1721 * t1757 - t1712 * t1774 - t1716 * t1752 + t1720 * t1756 + Icges(4,3) + m(4) * (t2049 ^ 2 + t2050 ^ 2);];
MX  = t1;