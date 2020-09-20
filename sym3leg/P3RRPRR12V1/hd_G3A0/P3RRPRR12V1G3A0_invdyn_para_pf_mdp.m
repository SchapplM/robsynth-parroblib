% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RRPRR12V1G3A0
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [15x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRPRR12V1G3A0_convert_par2_MPV_fixb.m

% Output:
% tauX [3x1]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:11
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRPRR12V1G3A0_invdyn_para_pf_mdp(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_mdp: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_mdp: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:10:50
% EndTime: 2020-08-06 19:11:04
% DurationCPUTime: 14.73s
% Computational Cost: add. (71853->585), mult. (123351->986), div. (7551->12), fcn. (78729->18), ass. (0->405)
t2029 = legFrame(3,2);
t2013 = cos(t2029);
t2274 = t2013 * g(1);
t2010 = sin(t2029);
t2277 = t2010 * g(2);
t1969 = t2274 - t2277;
t2036 = sin(qJ(1,3));
t2019 = g(3) * t2036;
t2042 = cos(qJ(1,3));
t2335 = t1969 * t2042 - t2019;
t2030 = legFrame(2,2);
t2014 = cos(t2030);
t2273 = t2014 * g(1);
t2011 = sin(t2030);
t2276 = t2011 * g(2);
t1970 = t2273 - t2276;
t2038 = sin(qJ(1,2));
t2020 = g(3) * t2038;
t2044 = cos(qJ(1,2));
t2334 = t1970 * t2044 - t2020;
t2031 = legFrame(1,2);
t2015 = cos(t2031);
t2272 = t2015 * g(1);
t2012 = sin(t2031);
t2275 = t2012 * g(2);
t1971 = t2272 - t2275;
t2040 = sin(qJ(1,1));
t2021 = g(3) * t2040;
t2046 = cos(qJ(1,1));
t2333 = t1971 * t2046 - t2021;
t2039 = sin(qJ(2,1));
t2045 = cos(qJ(2,1));
t2050 = pkin(1) + pkin(2);
t2174 = t2045 * t2050;
t1992 = qJ(3,1) * t2039 + t2174;
t1984 = 0.1e1 / t1992;
t2195 = t1984 * t2039;
t2047 = xDP(3);
t2048 = xDP(2);
t2049 = xDP(1);
t1938 = (t2012 * t2048 - t2015 * t2049) * t2040 - t2046 * t2047;
t2286 = pkin(4) * t1938;
t1923 = t2195 * t2286;
t2008 = t2048 * t2050;
t2009 = t2049 * t2050;
t2169 = t2047 * t2050;
t2266 = qJ(3,1) * t2049;
t2267 = qJ(3,1) * t2048;
t1908 = (t2046 * t2009 - t2267) * t2015 + (-t2046 * t2008 - t2266) * t2012 - t2040 * t2169;
t2282 = pkin(4) * t2048;
t2007 = t2040 * t2282;
t2027 = t2045 ^ 2;
t2281 = pkin(4) * t2049;
t2157 = t2040 * t2281;
t2018 = t2046 * pkin(4);
t2270 = qJ(3,1) * t2040;
t1976 = t2039 * t2270 + t2018;
t2170 = t2047 * t1976;
t2312 = -t2046 * t2267 + t2009;
t2315 = t2046 * t2266 + t2008;
t1890 = t1908 * t2027 + ((t2315 * t2039 - t2157) * t2015 + (t2312 * t2039 + t2007) * t2012 - t2170) * t2045 + qJ(3,1) * (t2012 * t2049 + t2015 * t2048);
t2058 = 0.1e1 / qJ(3,1);
t2234 = t1890 * t2058;
t2128 = t1984 * t2234;
t1878 = t2050 * t2128 + t1923;
t2268 = qJ(3,1) * t2046;
t2283 = pkin(4) * t2040;
t1977 = t2039 * t2268 - t2283;
t2173 = t2046 * t2050;
t2180 = t2039 * t2050;
t2248 = t2012 * qJ(3,1);
t1902 = (t2015 * t2173 - t2248) * t2027 + (t1977 * t2015 + t2012 * t2180) * t2045 + t2248;
t2245 = t2015 * qJ(3,1);
t1905 = (-t2012 * t2173 - t2245) * t2027 + (-t1977 * t2012 + t2015 * t2180) * t2045 + t2245;
t1953 = t2040 * t2174 + t1976;
t2032 = xDDP(3);
t2033 = xDDP(2);
t2034 = xDDP(1);
t1896 = -t2015 * t2157 + t1908 * t2045 - t2170 + t2007 * t2012 + (t2312 * t2012 + t2315 * t2015) * t2039;
t2059 = 0.1e1 / qJ(3,1) ^ 2;
t2230 = t1896 * t2059;
t2329 = ((t1953 * t2032 * t2045 - t1902 * t2034 - t1905 * t2033) * t2058 + (t1878 * t2045 + t1890 * t2195) * t2230) * t1984;
t2037 = sin(qJ(2,2));
t2043 = cos(qJ(2,2));
t2176 = t2043 * t2050;
t1991 = qJ(3,2) * t2037 + t2176;
t1981 = 0.1e1 / t1991;
t2197 = t1981 * t2037;
t1937 = (t2011 * t2048 - t2014 * t2049) * t2038 - t2044 * t2047;
t2287 = pkin(4) * t1937;
t1922 = t2197 * t2287;
t2260 = qJ(3,2) * t2049;
t2261 = qJ(3,2) * t2048;
t1907 = (t2044 * t2009 - t2261) * t2014 + (-t2044 * t2008 - t2260) * t2011 - t2038 * t2169;
t2006 = t2038 * t2282;
t2026 = t2043 ^ 2;
t2158 = t2038 * t2281;
t2017 = t2044 * pkin(4);
t2264 = qJ(3,2) * t2038;
t1974 = t2037 * t2264 + t2017;
t2171 = t2047 * t1974;
t2313 = -t2044 * t2261 + t2009;
t2316 = t2044 * t2260 + t2008;
t1889 = t1907 * t2026 + ((t2316 * t2037 - t2158) * t2014 + (t2313 * t2037 + t2006) * t2011 - t2171) * t2043 + qJ(3,2) * (t2011 * t2049 + t2014 * t2048);
t2055 = 0.1e1 / qJ(3,2);
t2236 = t1889 * t2055;
t2130 = t1981 * t2236;
t1877 = t2050 * t2130 + t1922;
t2262 = qJ(3,2) * t2044;
t2284 = pkin(4) * t2038;
t1975 = t2037 * t2262 - t2284;
t2175 = t2044 * t2050;
t2183 = t2037 * t2050;
t2249 = t2011 * qJ(3,2);
t1901 = (t2014 * t2175 - t2249) * t2026 + (t1975 * t2014 + t2011 * t2183) * t2043 + t2249;
t2246 = t2014 * qJ(3,2);
t1904 = (-t2011 * t2175 - t2246) * t2026 + (-t1975 * t2011 + t2014 * t2183) * t2043 + t2246;
t1952 = t2038 * t2176 + t1974;
t1895 = -t2014 * t2158 + t1907 * t2043 - t2171 + t2006 * t2011 + (t2313 * t2011 + t2316 * t2014) * t2037;
t2056 = 0.1e1 / qJ(3,2) ^ 2;
t2231 = t1895 * t2056;
t2328 = ((t1952 * t2032 * t2043 - t1901 * t2034 - t1904 * t2033) * t2055 + (t1877 * t2043 + t1889 * t2197) * t2231) * t1981;
t2035 = sin(qJ(2,3));
t2041 = cos(qJ(2,3));
t2178 = t2041 * t2050;
t1990 = qJ(3,3) * t2035 + t2178;
t1978 = 0.1e1 / t1990;
t2199 = t1978 * t2035;
t1936 = (t2010 * t2048 - t2013 * t2049) * t2036 - t2042 * t2047;
t2288 = pkin(4) * t1936;
t1921 = t2199 * t2288;
t2254 = qJ(3,3) * t2049;
t2255 = qJ(3,3) * t2048;
t1906 = (t2042 * t2009 - t2255) * t2013 + (-t2042 * t2008 - t2254) * t2010 - t2036 * t2169;
t2005 = t2036 * t2282;
t2025 = t2041 ^ 2;
t2159 = t2036 * t2281;
t2016 = t2042 * pkin(4);
t2258 = qJ(3,3) * t2036;
t1972 = t2035 * t2258 + t2016;
t2172 = t2047 * t1972;
t2314 = -t2042 * t2255 + t2009;
t2317 = t2042 * t2254 + t2008;
t1888 = t1906 * t2025 + ((t2317 * t2035 - t2159) * t2013 + (t2314 * t2035 + t2005) * t2010 - t2172) * t2041 + qJ(3,3) * (t2010 * t2049 + t2013 * t2048);
t2052 = 0.1e1 / qJ(3,3);
t2238 = t1888 * t2052;
t2132 = t1978 * t2238;
t1876 = t2050 * t2132 + t1921;
t2256 = qJ(3,3) * t2042;
t2285 = pkin(4) * t2036;
t1973 = t2035 * t2256 - t2285;
t2177 = t2042 * t2050;
t2186 = t2035 * t2050;
t2250 = t2010 * qJ(3,3);
t1900 = (t2013 * t2177 - t2250) * t2025 + (t1973 * t2013 + t2010 * t2186) * t2041 + t2250;
t2247 = t2013 * qJ(3,3);
t1903 = (-t2010 * t2177 - t2247) * t2025 + (-t1973 * t2010 + t2013 * t2186) * t2041 + t2247;
t1951 = t2036 * t2178 + t1972;
t1894 = -t2013 * t2159 + t1906 * t2041 - t2172 + t2005 * t2010 + (t2314 * t2010 + t2317 * t2013) * t2035;
t2053 = 0.1e1 / qJ(3,3) ^ 2;
t2232 = t1894 * t2053;
t2327 = ((t1951 * t2032 * t2041 - t1900 * t2034 - t1903 * t2033) * t2052 + (t1876 * t2041 + t1888 * t2199) * t2232) * t1978;
t2326 = t1990 * t2042 - t2285;
t2325 = t1991 * t2044 - t2284;
t2324 = t1992 * t2046 - t2283;
t1939 = t2334 * t2037;
t1967 = g(1) * t2011 + g(2) * t2014;
t2323 = -t1967 * t2043 + t1939;
t1940 = t2335 * t2035;
t1966 = g(1) * t2010 + g(2) * t2013;
t2322 = -t1966 * t2041 + t1940;
t1941 = t2333 * t2039;
t1968 = g(1) * t2012 + g(2) * t2015;
t2321 = -t1968 * t2045 + t1941;
t2022 = g(3) * t2042;
t2320 = t1969 * t2036 + t2022;
t2023 = g(3) * t2044;
t2319 = t1970 * t2038 + t2023;
t2024 = g(3) * t2046;
t2318 = t1971 * t2040 + t2024;
t1979 = 0.1e1 / t1990 ^ 2;
t2051 = qJ(3,3) ^ 2;
t2060 = pkin(4) ^ 2;
t2217 = t1936 * t1978;
t2122 = t2025 * t2217;
t2215 = t1936 * t2052;
t1880 = pkin(1) * t2132;
t1891 = t1894 * t2052;
t1862 = t1880 - t1891;
t1849 = pkin(2) * t2132 + t1862;
t2244 = t1849 * t2035;
t2257 = qJ(3,3) * t2041;
t2300 = -pkin(4) / 0.2e1;
t2144 = ((qJ(3,3) + t2050) * (-qJ(3,3) + t2050) * t2122 + 0.2e1 * (t1936 * t2186 + t2238 * t2300) * t1978 * t2257 + pkin(4) * t2244 + (t2051 + t2060) * t2217) * t2215;
t2237 = t1888 * t2053;
t2150 = (-(t1921 + t1849) * t2178 + (pkin(4) * t2122 - t2244) * qJ(3,3)) * t2237;
t1822 = (t2041 * t2144 - t2150) * t1979 - t2327;
t1980 = t1978 * t1979;
t2093 = -t2186 + t2257;
t2198 = t1979 * t2052;
t2120 = t1936 * t2198;
t2168 = t2050 * t2052;
t2190 = t2013 * t2036;
t2193 = t2010 * t2036;
t2216 = t1936 * t1979;
t1837 = -t2093 * t1980 * t1888 * t2215 + (-t1891 * t2216 - t1894 * t2120) * t2035 + (-t2034 * t2190 + t2033 * t2193 - t2042 * t2032 - (-t2288 + (-t2035 * t2168 + t2041) * t1888) * t2216) * t1978;
t2259 = qJ(3,3) * t1837;
t1834 = 0.2e1 * t2259;
t2131 = t1978 * t2237;
t1867 = 0.2e1 * t1894 * t2131;
t1885 = t1888 ^ 2;
t2062 = pkin(1) ^ 2;
t2111 = -t2062 + (-0.2e1 * pkin(1) - pkin(2)) * pkin(2);
t1960 = -t1990 * t2036 - t2016;
t2208 = t1960 * t2052;
t1927 = -t2010 * t2326 - t2093 * t2013;
t2226 = t1927 * t2052;
t1924 = -t2093 * t2010 + t2013 * t2326;
t2229 = t1924 * t2052;
t2089 = -t1978 * t2144 + (t1894 * t2168 + ((-t2051 + t2111) * t2238 + t2093 * t2288) * t1978) * t2131 + t1876 * t2232 - t2034 * t2229 - t2033 * t2226 - t2032 * t2208;
t2280 = t1822 * pkin(1);
t2086 = t2089 + t2280;
t2110 = t1888 * t2120;
t2121 = t1888 * t2216;
t2123 = t1862 * t2217;
t2241 = t1885 * t2053;
t2135 = t1979 * t2241;
t2138 = (0.4e1 * t1880 - 0.2e1 * t1891) * t2217;
t2162 = -t2051 + t2062;
t2187 = t2035 * t2041;
t2253 = t1822 * qJ(3,3);
t2294 = pkin(1) * t1837;
t2297 = 0.2e1 * t2025 - 0.1e1;
t2303 = 0.2e1 * t2041;
t2304 = -0.2e1 * qJ(3,3);
t2307 = -2 * MDP(5);
t2070 = (MDP(10) * t2035 - MDP(9) * t2041) * t2320 - t1837 * MDP(1) - ((0.4e1 * t2121 + 0.2e1 * t2294) * t2025 + ((t1834 - t2138) * t2035 + t2320) * t2041 - 0.2e1 * t2121) * MDP(11) - ((-pkin(1) * t2135 + t1867 + t2253) * t2041 + (-t1885 * t2198 - t2086) * t2035 - t2335) * MDP(12) - ((t2138 - 0.2e1 * t2259) * t2025 - 0.2e1 * t2123 + t1834 + ((0.2e1 * t2121 + t2294) * t2303 + t2320) * t2035) * MDP(13) - ((0.4e1 * (t1880 - t1891 / 0.2e1) * qJ(3,3) * t2217 + t2162 * t1837) * t2025 + (0.2e1 * (qJ(3,3) * t2121 + (-t2123 + t2259) * pkin(1)) * t2035 + t2320 * pkin(1)) * t2041 + (((-t2274 / 0.2e1 + t2277 / 0.2e1) * t2036 - t2022 / 0.2e1) * t2035 - t2259 / 0.2e1 + t2123) * t2304) * MDP(14) - t2320 * MDP(2) - t2335 * MDP(3) - (t1837 * t2035 + t2110 * t2303) * t2035 * MDP(4) - (t1822 * t2035 + t2041 * t2135) * MDP(6) - (t1822 * t2041 - t2035 * t2135) * MDP(7) + (t1837 * t2187 + t2297 * t2110) * t2307;
t1985 = 0.1e1 / t1992 ^ 2;
t2057 = qJ(3,1) ^ 2;
t2211 = t1938 * t1984;
t2114 = t2027 * t2211;
t2209 = t1938 * t2058;
t1884 = pkin(1) * t2128;
t1893 = t1896 * t2058;
t1866 = t1884 - t1893;
t1851 = pkin(2) * t2128 + t1866;
t2242 = t1851 * t2039;
t2269 = qJ(3,1) * t2045;
t2142 = ((qJ(3,1) + t2050) * (-qJ(3,1) + t2050) * t2114 + 0.2e1 * (t1938 * t2180 + t2234 * t2300) * t1984 * t2269 + pkin(4) * t2242 + (t2057 + t2060) * t2211) * t2209;
t2233 = t1890 * t2059;
t2148 = (-(t1923 + t1851) * t2174 + (pkin(4) * t2114 - t2242) * qJ(3,1)) * t2233;
t1824 = (t2045 * t2142 - t2148) * t1985 - t2329;
t1986 = t1984 * t1985;
t2095 = -t2180 + t2269;
t2194 = t1985 * t2058;
t2112 = t1938 * t2194;
t2166 = t2050 * t2058;
t2188 = t2015 * t2040;
t2191 = t2012 * t2040;
t2210 = t1938 * t1985;
t1839 = -t2095 * t1986 * t1890 * t2209 + (-t1893 * t2210 - t1896 * t2112) * t2039 + (-t2034 * t2188 + t2033 * t2191 - t2046 * t2032 - (-t2286 + (-t2039 * t2166 + t2045) * t1890) * t2210) * t1984;
t2271 = qJ(3,1) * t1839;
t1836 = 0.2e1 * t2271;
t2127 = t1984 * t2233;
t1869 = 0.2e1 * t1896 * t2127;
t1887 = t1890 ^ 2;
t1962 = -t1992 * t2040 - t2018;
t2206 = t1962 * t2058;
t1929 = -t2012 * t2324 - t2095 * t2015;
t2224 = t1929 * t2058;
t1926 = -t2095 * t2012 + t2015 * t2324;
t2227 = t1926 * t2058;
t2087 = -t1984 * t2142 + (t1896 * t2166 + ((-t2057 + t2111) * t2234 + t2095 * t2286) * t1984) * t2127 + t1878 * t2230 - t2034 * t2227 - t2033 * t2224 - t2032 * t2206;
t2278 = t1824 * pkin(1);
t2084 = t2087 + t2278;
t2108 = t1890 * t2112;
t2113 = t1890 * t2210;
t2115 = t1866 * t2211;
t2239 = t1887 * t2059;
t2133 = t1985 * t2239;
t2136 = (0.4e1 * t1884 - 0.2e1 * t1893) * t2211;
t2160 = -t2057 + t2062;
t2181 = t2039 * t2045;
t2251 = t1824 * qJ(3,1);
t2292 = pkin(1) * t1839;
t2295 = 0.2e1 * t2027 - 0.1e1;
t2301 = 0.2e1 * t2045;
t2306 = -0.2e1 * qJ(3,1);
t2069 = (MDP(10) * t2039 - MDP(9) * t2045) * t2318 - t1839 * MDP(1) - ((0.4e1 * t2113 + 0.2e1 * t2292) * t2027 + ((t1836 - t2136) * t2039 + t2318) * t2045 - 0.2e1 * t2113) * MDP(11) - ((-pkin(1) * t2133 + t1869 + t2251) * t2045 + (-t1887 * t2194 - t2084) * t2039 - t2333) * MDP(12) - ((t2136 - 0.2e1 * t2271) * t2027 - 0.2e1 * t2115 + t1836 + ((0.2e1 * t2113 + t2292) * t2301 + t2318) * t2039) * MDP(13) - ((0.4e1 * (t1884 - t1893 / 0.2e1) * qJ(3,1) * t2211 + t2160 * t1839) * t2027 + (0.2e1 * (qJ(3,1) * t2113 + (-t2115 + t2271) * pkin(1)) * t2039 + t2318 * pkin(1)) * t2045 + (((-t2272 / 0.2e1 + t2275 / 0.2e1) * t2040 - t2024 / 0.2e1) * t2039 - t2271 / 0.2e1 + t2115) * t2306) * MDP(14) - t2318 * MDP(2) - t2333 * MDP(3) - (t1839 * t2039 + t2108 * t2301) * t2039 * MDP(4) - (t1824 * t2039 + t2045 * t2133) * MDP(6) - (t1824 * t2045 - t2039 * t2133) * MDP(7) + (t1839 * t2181 + t2295 * t2108) * t2307;
t1982 = 0.1e1 / t1991 ^ 2;
t2054 = qJ(3,2) ^ 2;
t2214 = t1937 * t1981;
t2118 = t2026 * t2214;
t2212 = t1937 * t2055;
t1882 = pkin(1) * t2130;
t1892 = t1895 * t2055;
t1864 = t1882 - t1892;
t1850 = pkin(2) * t2130 + t1864;
t2243 = t1850 * t2037;
t2263 = qJ(3,2) * t2043;
t2143 = ((qJ(3,2) + t2050) * (-qJ(3,2) + t2050) * t2118 + 0.2e1 * (t1937 * t2183 + t2236 * t2300) * t1981 * t2263 + pkin(4) * t2243 + (t2054 + t2060) * t2214) * t2212;
t2235 = t1889 * t2056;
t2149 = (-(t1922 + t1850) * t2176 + (pkin(4) * t2118 - t2243) * qJ(3,2)) * t2235;
t1823 = (t2043 * t2143 - t2149) * t1982 - t2328;
t1983 = t1981 * t1982;
t2094 = -t2183 + t2263;
t2196 = t1982 * t2055;
t2116 = t1937 * t2196;
t2167 = t2050 * t2055;
t2189 = t2014 * t2038;
t2192 = t2011 * t2038;
t2213 = t1937 * t1982;
t1838 = -t2094 * t1983 * t1889 * t2212 + (-t1892 * t2213 - t1895 * t2116) * t2037 + (-t2034 * t2189 + t2033 * t2192 - t2044 * t2032 - (-t2287 + (-t2037 * t2167 + t2043) * t1889) * t2213) * t1981;
t2265 = qJ(3,2) * t1838;
t1835 = 0.2e1 * t2265;
t2129 = t1981 * t2235;
t1868 = 0.2e1 * t1895 * t2129;
t1886 = t1889 ^ 2;
t1961 = -t1991 * t2038 - t2017;
t2207 = t1961 * t2055;
t1928 = -t2011 * t2325 - t2094 * t2014;
t2225 = t1928 * t2055;
t1925 = -t2094 * t2011 + t2014 * t2325;
t2228 = t1925 * t2055;
t2088 = -t1981 * t2143 + (t1895 * t2167 + ((-t2054 + t2111) * t2236 + t2094 * t2287) * t1981) * t2129 + t1877 * t2231 - t2034 * t2228 - t2033 * t2225 - t2032 * t2207;
t2279 = t1823 * pkin(1);
t2085 = t2088 + t2279;
t2109 = t1889 * t2116;
t2117 = t1889 * t2213;
t2119 = t1864 * t2214;
t2240 = t1886 * t2056;
t2134 = t1982 * t2240;
t2137 = (0.4e1 * t1882 - 0.2e1 * t1892) * t2214;
t2161 = -t2054 + t2062;
t2184 = t2037 * t2043;
t2252 = t1823 * qJ(3,2);
t2293 = pkin(1) * t1838;
t2296 = 0.2e1 * t2026 - 0.1e1;
t2302 = 0.2e1 * t2043;
t2305 = -0.2e1 * qJ(3,2);
t2071 = (MDP(10) * t2037 - MDP(9) * t2043) * t2319 - t1838 * MDP(1) - ((0.4e1 * t2117 + 0.2e1 * t2293) * t2026 + ((t1835 - t2137) * t2037 + t2319) * t2043 - 0.2e1 * t2117) * MDP(11) - ((-pkin(1) * t2134 + t1868 + t2252) * t2043 + (-t1886 * t2196 - t2085) * t2037 - t2334) * MDP(12) - ((t2137 - 0.2e1 * t2265) * t2026 - 0.2e1 * t2119 + t1835 + ((0.2e1 * t2117 + t2293) * t2302 + t2319) * t2037) * MDP(13) - ((0.4e1 * (t1882 - t1892 / 0.2e1) * qJ(3,2) * t2214 + t2161 * t1838) * t2026 + (0.2e1 * (qJ(3,2) * t2117 + (-t2119 + t2265) * pkin(1)) * t2037 + t2319 * pkin(1)) * t2043 + (((-t2273 / 0.2e1 + t2276 / 0.2e1) * t2038 - t2023 / 0.2e1) * t2037 - t2265 / 0.2e1 + t2119) * t2305) * MDP(14) - t2319 * MDP(2) - t2334 * MDP(3) - (t1838 * t2037 + t2109 * t2302) * t2037 * MDP(4) - (t1823 * t2037 + t2043 * t2134) * MDP(6) - (t1823 * t2043 - t2037 * t2134) * MDP(7) + (t1838 * t2184 + t2296 * t2109) * t2307;
t2308 = 0.2e1 * pkin(1);
t2299 = pkin(1) * g(1);
t2298 = pkin(1) * g(2);
t2291 = pkin(1) * t2042;
t2290 = pkin(1) * t2044;
t2289 = pkin(1) * t2046;
t1933 = t1936 ^ 2;
t2223 = t1933 * t1979;
t2222 = t1933 * t2035;
t1934 = t1937 ^ 2;
t2221 = t1934 * t1982;
t2220 = t1934 * t2037;
t1935 = t1938 ^ 2;
t2219 = t1935 * t1985;
t2218 = t1935 * t2039;
t2205 = t1966 * t2035;
t2203 = t1967 * t2037;
t2201 = t1968 * t2039;
t2185 = t2035 * t2052;
t2182 = t2037 * t2055;
t2179 = t2039 * t2058;
t1930 = pkin(1) * t2223;
t2165 = t1930 + t1867;
t1931 = pkin(1) * t2221;
t2164 = t1931 + t1868;
t1932 = pkin(1) * t2219;
t2163 = t1932 + t1869;
t2147 = t1837 * t2185;
t2146 = t1838 * t2182;
t2145 = t1839 * t2179;
t2126 = t1979 * t2222;
t2125 = t1982 * t2220;
t2124 = t1985 * t2218;
t1897 = t2297 * t2223;
t1898 = t2296 * t2221;
t1899 = t2295 * t2219;
t2107 = t1933 * t1980 * t2185;
t2106 = t1934 * t1983 * t2182;
t2105 = t1935 * t1986 * t2179;
t2104 = -0.2e1 * t2025 * t1930;
t2103 = -0.2e1 * t2026 * t1931;
t2102 = -0.2e1 * t2027 * t1932;
t2098 = -MDP(12) * (-pkin(1) * t2035 + t2257) - MDP(6) * t2035;
t2097 = -MDP(12) * (-pkin(1) * t2037 + t2263) - MDP(6) * t2037;
t2096 = -MDP(12) * (-pkin(1) * t2039 + t2269) - MDP(6) * t2039;
t2092 = t2041 * t2107;
t2091 = t2043 * t2106;
t2090 = t2045 * t2105;
t2080 = -MDP(10) * (t2041 * t2335 + t2205) - MDP(11) * ((t2126 * t2308 - t1966) * t2041 + t1940 + 0.2e1 * t2280 - qJ(3,3) * t1897 + t2089) - MDP(13) * (t2104 + (t2126 * t2304 - t2335) * t2041 - t2205 + 0.2e1 * t2253 + t2165) - MDP(14) * (qJ(3,3) * t2104 + ((-g(1) * t2256 - t2298) * t2013 + (g(2) * t2256 - t2299) * t2010 + g(3) * t2258 + t2162 * t2126) * t2041 + ((g(1) * t2291 - g(2) * qJ(3,3)) * t2013 + (-g(1) * qJ(3,3) - g(2) * t2291) * t2010 - pkin(1) * t2019) * t2035 + t1822 * t2051 + t2165 * qJ(3,3) + t2086 * pkin(1)) + MDP(5) * t1897 - MDP(8) * t1822 - MDP(9) * t2322;
t2079 = -MDP(10) * (t2043 * t2334 + t2203) - MDP(11) * ((t2125 * t2308 - t1967) * t2043 + t1939 + 0.2e1 * t2279 - qJ(3,2) * t1898 + t2088) - MDP(13) * (t2103 + (t2125 * t2305 - t2334) * t2043 - t2203 + 0.2e1 * t2252 + t2164) - MDP(14) * (qJ(3,2) * t2103 + ((-g(1) * t2262 - t2298) * t2014 + (g(2) * t2262 - t2299) * t2011 + g(3) * t2264 + t2161 * t2125) * t2043 + ((g(1) * t2290 - g(2) * qJ(3,2)) * t2014 + (-g(1) * qJ(3,2) - g(2) * t2290) * t2011 - pkin(1) * t2020) * t2037 + t1823 * t2054 + t2164 * qJ(3,2) + t2085 * pkin(1)) + MDP(5) * t1898 - MDP(8) * t1823 - MDP(9) * t2323;
t2078 = -MDP(10) * (t2045 * t2333 + t2201) - MDP(11) * ((t2124 * t2308 - t1968) * t2045 + t1941 + 0.2e1 * t2278 - qJ(3,1) * t1899 + t2087) - MDP(13) * (t2102 + (t2124 * t2306 - t2333) * t2045 - t2201 + 0.2e1 * t2251 + t2163) - MDP(14) * (qJ(3,1) * t2102 + ((-g(1) * t2268 - t2298) * t2015 + (g(2) * t2268 - t2299) * t2012 + g(3) * t2270 + t2160 * t2124) * t2045 + ((g(1) * t2289 - g(2) * qJ(3,1)) * t2015 + (-g(1) * qJ(3,1) - g(2) * t2289) * t2012 - pkin(1) * t2021) * t2039 + t1824 * t2057 + t2163 * qJ(3,1) + t2084 * pkin(1)) + MDP(5) * t1899 - MDP(8) * t1824 - MDP(9) * t2321;
t2074 = t2052 * ((MDP(7) * t2041 - t2098) * t1837 - t2080);
t2073 = t2055 * ((MDP(7) * t2043 - t2097) * t1838 - t2079);
t2072 = t2058 * ((MDP(7) * t2045 - t2096) * t1839 - t2078);
t1872 = (t1935 * t2027 - t1935 - t2239) * t1985;
t1871 = (t1934 * t2026 - t1934 - t2240) * t1982;
t1870 = (t1933 * t2025 - t1933 - t2241) * t1979;
t1818 = (t2148 + (-t2142 - t2218) * t2045) * t1985 + t2329;
t1817 = (t2149 + (-t2143 - t2220) * t2043) * t1982 + t2328;
t1816 = (t2150 + (-t2144 - t2222) * t2041) * t1979 + t2327;
t1800 = (-t1887 * t2058 + (-pkin(1) * t2181 + qJ(3,1) * t2027 - qJ(3,1)) * t1935) * t1985 - t2084 - t2321;
t1799 = (-t1886 * t2055 + (-pkin(1) * t2184 + qJ(3,2) * t2026 - qJ(3,2)) * t1934) * t1982 - t2085 - t2323;
t1798 = (-t1885 * t2052 + (-pkin(1) * t2187 + qJ(3,3) * t2025 - qJ(3,3)) * t1933) * t1979 - t2086 - t2322;
t1 = [(-t1900 * t2092 - t1901 * t2091 - t1902 * t2090) * MDP(4) + (t1816 * t2229 + t1817 * t2228 + t1818 * t2227) * MDP(11) + (t1924 * t2147 + t1925 * t2146 + t1926 * t2145) * MDP(12) + (t1870 * t2229 + t1871 * t2228 + t1872 * t2227) * MDP(13) + (t1798 * t2229 + t1799 * t2228 + t1800 * t2227) * MDP(14) + (t2034 - g(1)) * MDP(15) + (t1902 * t2072 + t2069 * t2188) * t1984 + (t1901 * t2073 + t2071 * t2189) * t1981 + (t1900 * t2074 + t2070 * t2190) * t1978; (-t1903 * t2092 - t1904 * t2091 - t1905 * t2090) * MDP(4) + (t1816 * t2226 + t1817 * t2225 + t1818 * t2224) * MDP(11) + (t1927 * t2147 + t1928 * t2146 + t1929 * t2145) * MDP(12) + (t1870 * t2226 + t1871 * t2225 + t1872 * t2224) * MDP(13) + (t1798 * t2226 + t1799 * t2225 + t1800 * t2224) * MDP(14) + (t2033 - g(2)) * MDP(15) + (t1905 * t2072 - t2069 * t2191) * t1984 + (t1904 * t2073 - t2071 * t2192) * t1981 + (t1903 * t2074 - t2070 * t2193) * t1978; (t1951 * t2025 * t2107 + t1952 * t2026 * t2106 + t1953 * t2027 * t2105) * MDP(4) + (t1816 * t2208 + t1817 * t2207 + t1818 * t2206) * MDP(11) + (t1960 * t2147 + t1961 * t2146 + t1962 * t2145) * MDP(12) + (t1870 * t2208 + t1871 * t2207 + t1872 * t2206) * MDP(13) + (t1798 * t2208 + t1799 * t2207 + t1800 * t2206) * MDP(14) + (t2032 - g(3)) * MDP(15) + ((-MDP(7) * t1839 * t2027 + (t1839 * t2096 + t2078) * t2045) * t2058 * t1953 + t2069 * t2046) * t1984 + ((-MDP(7) * t1838 * t2026 + (t1838 * t2097 + t2079) * t2043) * t2055 * t1952 + t2071 * t2044) * t1981 + ((-MDP(7) * t1837 * t2025 + (t1837 * t2098 + t2080) * t2041) * t2052 * t1951 + t2070 * t2042) * t1978;];
tauX  = t1;
