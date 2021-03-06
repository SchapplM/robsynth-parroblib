% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RRPRR12V1G2A0
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
%   see P3RRPRR12V1G2A0_convert_par2_MPV_fixb.m

% Output:
% tauX [3x1]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:07
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRPRR12V1G2A0_invdyn_para_pf_mdp(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_mdp: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_mdp: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:06:26
% EndTime: 2020-08-06 19:06:39
% DurationCPUTime: 13.19s
% Computational Cost: add. (71853->585), mult. (123438->986), div. (7551->12), fcn. (78816->18), ass. (0->405)
t2050 = sin(qJ(1,3));
t2036 = g(3) * t2050;
t2056 = cos(qJ(1,3));
t2043 = legFrame(3,2);
t2030 = sin(t2043);
t2025 = t2030 * g(2);
t2033 = cos(t2043);
t2291 = t2033 * g(1);
t2335 = -t2291 + t2025;
t2182 = t2335 * t2056 + t2036;
t2052 = sin(qJ(1,2));
t2037 = g(3) * t2052;
t2058 = cos(qJ(1,2));
t2044 = legFrame(2,2);
t2031 = sin(t2044);
t2026 = t2031 * g(2);
t2034 = cos(t2044);
t2290 = t2034 * g(1);
t2336 = -t2290 + t2026;
t2181 = t2336 * t2058 + t2037;
t2054 = sin(qJ(1,1));
t2038 = g(3) * t2054;
t2060 = cos(qJ(1,1));
t2045 = legFrame(1,2);
t2032 = sin(t2045);
t2027 = t2032 * g(2);
t2035 = cos(t2045);
t2289 = t2035 * g(1);
t2337 = -t2289 + t2027;
t2180 = t2337 * t2060 + t2038;
t2297 = g(3) * t2056;
t2352 = -t2335 * t2050 + t2297;
t2296 = g(3) * t2058;
t2351 = -t2336 * t2052 + t2296;
t2295 = g(3) * t2060;
t2350 = -t2337 * t2054 + t2295;
t2053 = sin(qJ(2,1));
t2059 = cos(qJ(2,1));
t2064 = pkin(1) + pkin(2);
t2190 = t2059 * t2064;
t2012 = qJ(3,1) * t2053 + t2190;
t2004 = 0.1e1 / t2012;
t2212 = t2004 * t2053;
t2061 = xDP(3);
t2062 = xDP(2);
t2063 = xDP(1);
t1955 = (-t2032 * t2062 + t2035 * t2063) * t2060 - t2054 * t2061;
t2306 = pkin(4) * t1955;
t1946 = t2212 * t2306;
t2028 = t2062 * t2064;
t2029 = t2063 * t2064;
t2189 = t2061 * t2064;
t2283 = qJ(3,1) * t2063;
t2284 = qJ(3,1) * t2062;
t1931 = (t2054 * t2029 - t2284) * t2035 + (-t2054 * t2028 - t2283) * t2032 + t2060 * t2189;
t2285 = qJ(3,1) * t2060;
t2303 = pkin(4) * t2054;
t2125 = t2053 * t2285 - t2303;
t1988 = t2061 * t2125;
t2298 = pkin(4) * t2063;
t2024 = t2060 * t2298;
t2041 = t2059 ^ 2;
t2299 = pkin(4) * t2062;
t2174 = t2060 * t2299;
t2338 = -t2054 * t2284 + t2029;
t2341 = t2054 * t2283 + t2028;
t1907 = t1931 * t2041 + ((t2341 * t2053 + t2024) * t2035 + (t2338 * t2053 - t2174) * t2032 + t1988) * t2059 + qJ(3,1) * (t2032 * t2063 + t2035 * t2062);
t2072 = 0.1e1 / qJ(3,1);
t2251 = t1907 * t2072;
t2145 = t2004 * t2251;
t1895 = t2064 * t2145 + t1946;
t2287 = qJ(3,1) * t2054;
t2300 = pkin(4) * t2060;
t1997 = t2053 * t2287 + t2300;
t2193 = t2054 * t2064;
t2195 = t2053 * t2064;
t2265 = t2032 * qJ(3,1);
t1919 = (t2035 * t2193 - t2265) * t2041 + (t1997 * t2035 + t2032 * t2195) * t2059 + t2265;
t2262 = t2035 * qJ(3,1);
t1922 = (-t2032 * t2193 - t2262) * t2041 + (-t1997 * t2032 + t2035 * t2195) * t2059 + t2262;
t1973 = t2060 * t2190 + t2125;
t2046 = xDDP(3);
t2047 = xDDP(2);
t2048 = xDDP(1);
t1913 = -t2032 * t2174 + t1931 * t2059 + t2024 * t2035 + t1988 + (t2338 * t2032 + t2341 * t2035) * t2053;
t2073 = 0.1e1 / qJ(3,1) ^ 2;
t2247 = t1913 * t2073;
t2349 = ((t1973 * t2046 * t2059 + t1919 * t2048 + t1922 * t2047) * t2072 - (t1895 * t2059 + t1907 * t2212) * t2247) * t2004;
t2051 = sin(qJ(2,2));
t2057 = cos(qJ(2,2));
t2191 = t2057 * t2064;
t2011 = qJ(3,2) * t2051 + t2191;
t2001 = 0.1e1 / t2011;
t2214 = t2001 * t2051;
t1954 = (-t2031 * t2062 + t2034 * t2063) * t2058 - t2052 * t2061;
t2307 = pkin(4) * t1954;
t1945 = t2214 * t2307;
t2277 = qJ(3,2) * t2063;
t2278 = qJ(3,2) * t2062;
t1930 = (t2052 * t2029 - t2278) * t2034 + (-t2052 * t2028 - t2277) * t2031 + t2058 * t2189;
t2279 = qJ(3,2) * t2058;
t2304 = pkin(4) * t2052;
t2126 = t2051 * t2279 - t2304;
t1987 = t2061 * t2126;
t2023 = t2058 * t2298;
t2040 = t2057 ^ 2;
t2175 = t2058 * t2299;
t2339 = -t2052 * t2278 + t2029;
t2342 = t2052 * t2277 + t2028;
t1906 = t1930 * t2040 + ((t2342 * t2051 + t2023) * t2034 + (t2339 * t2051 - t2175) * t2031 + t1987) * t2057 + qJ(3,2) * (t2031 * t2063 + t2034 * t2062);
t2069 = 0.1e1 / qJ(3,2);
t2253 = t1906 * t2069;
t2147 = t2001 * t2253;
t1894 = t2064 * t2147 + t1945;
t2281 = qJ(3,2) * t2052;
t2301 = pkin(4) * t2058;
t1996 = t2051 * t2281 + t2301;
t2197 = t2052 * t2064;
t2199 = t2051 * t2064;
t2266 = t2031 * qJ(3,2);
t1918 = (t2034 * t2197 - t2266) * t2040 + (t1996 * t2034 + t2031 * t2199) * t2057 + t2266;
t2263 = t2034 * qJ(3,2);
t1921 = (-t2031 * t2197 - t2263) * t2040 + (-t1996 * t2031 + t2034 * t2199) * t2057 + t2263;
t1972 = t2058 * t2191 + t2126;
t1912 = -t2031 * t2175 + t1930 * t2057 + t2023 * t2034 + t1987 + (t2339 * t2031 + t2342 * t2034) * t2051;
t2070 = 0.1e1 / qJ(3,2) ^ 2;
t2248 = t1912 * t2070;
t2348 = ((t1972 * t2046 * t2057 + t1918 * t2048 + t1921 * t2047) * t2069 - (t1894 * t2057 + t1906 * t2214) * t2248) * t2001;
t2049 = sin(qJ(2,3));
t2055 = cos(qJ(2,3));
t2192 = t2055 * t2064;
t2010 = qJ(3,3) * t2049 + t2192;
t1998 = 0.1e1 / t2010;
t2216 = t1998 * t2049;
t1953 = (-t2030 * t2062 + t2033 * t2063) * t2056 - t2050 * t2061;
t2308 = pkin(4) * t1953;
t1944 = t2216 * t2308;
t2271 = qJ(3,3) * t2063;
t2272 = qJ(3,3) * t2062;
t1929 = (t2050 * t2029 - t2272) * t2033 + (-t2050 * t2028 - t2271) * t2030 + t2056 * t2189;
t2273 = qJ(3,3) * t2056;
t2305 = pkin(4) * t2050;
t2127 = t2049 * t2273 - t2305;
t1986 = t2061 * t2127;
t2022 = t2056 * t2298;
t2039 = t2055 ^ 2;
t2176 = t2056 * t2299;
t2340 = -t2050 * t2272 + t2029;
t2343 = t2050 * t2271 + t2028;
t1905 = t1929 * t2039 + ((t2343 * t2049 + t2022) * t2033 + (t2340 * t2049 - t2176) * t2030 + t1986) * t2055 + qJ(3,3) * (t2030 * t2063 + t2033 * t2062);
t2066 = 0.1e1 / qJ(3,3);
t2255 = t1905 * t2066;
t2149 = t1998 * t2255;
t1893 = t2064 * t2149 + t1944;
t2275 = qJ(3,3) * t2050;
t2302 = pkin(4) * t2056;
t1995 = t2049 * t2275 + t2302;
t2201 = t2050 * t2064;
t2203 = t2049 * t2064;
t2267 = t2030 * qJ(3,3);
t1917 = (t2033 * t2201 - t2267) * t2039 + (t1995 * t2033 + t2030 * t2203) * t2055 + t2267;
t2264 = t2033 * qJ(3,3);
t1920 = (-t2030 * t2201 - t2264) * t2039 + (-t1995 * t2030 + t2033 * t2203) * t2055 + t2264;
t1971 = t2056 * t2192 + t2127;
t1911 = -t2030 * t2176 + t1929 * t2055 + t2022 * t2033 + t1986 + (t2340 * t2030 + t2343 * t2033) * t2049;
t2067 = 0.1e1 / qJ(3,3) ^ 2;
t2249 = t1911 * t2067;
t2347 = (t2066 * (t1971 * t2046 * t2055 + t1917 * t2048 + t1920 * t2047) - (t1893 * t2055 + t1905 * t2216) * t2249) * t1998;
t1956 = t2352 * t2049;
t1989 = g(1) * t2030 + g(2) * t2033;
t2346 = -t1989 * t2055 + t1956;
t1957 = t2351 * t2051;
t1990 = g(1) * t2031 + g(2) * t2034;
t2345 = -t1990 * t2057 + t1957;
t1958 = t2350 * t2053;
t1991 = g(1) * t2032 + g(2) * t2035;
t2344 = -t1991 * t2059 + t1958;
t1999 = 0.1e1 / t2010 ^ 2;
t2065 = qJ(3,3) ^ 2;
t2074 = pkin(4) ^ 2;
t2234 = t1953 * t1998;
t2139 = t2039 * t2234;
t2232 = t1953 * t2066;
t1897 = pkin(1) * t2149;
t1908 = t1911 * t2066;
t1879 = t1897 - t1908;
t1866 = pkin(2) * t2149 + t1879;
t2261 = t1866 * t2049;
t2274 = qJ(3,3) * t2055;
t2320 = -pkin(4) / 0.2e1;
t2161 = ((qJ(3,3) + t2064) * (-qJ(3,3) + t2064) * t2139 + 0.2e1 * (t1953 * t2203 + t2255 * t2320) * t1998 * t2274 + pkin(4) * t2261 + (t2065 + t2074) * t2234) * t2232;
t2254 = t1905 * t2067;
t2167 = (-(t1944 + t1866) * t2192 + (pkin(4) * t2139 - t2261) * qJ(3,3)) * t2254;
t1839 = (t2055 * t2161 - t2167) * t1999 + t2347;
t2000 = t1998 * t1999;
t2107 = -t2203 + t2274;
t2215 = t1999 * t2066;
t2137 = t1953 * t2215;
t2188 = t2064 * t2066;
t2207 = t2033 * t2056;
t2210 = t2030 * t2056;
t2233 = t1953 * t1999;
t1854 = -t2000 * t2107 * t1905 * t2232 + (-t1908 * t2233 - t1911 * t2137) * t2049 + (t2048 * t2207 - t2047 * t2210 - t2050 * t2046 - (-t2308 + (-t2049 * t2188 + t2055) * t1905) * t2233) * t1998;
t2276 = qJ(3,3) * t1854;
t1851 = 0.2e1 * t2276;
t2148 = t1998 * t2254;
t1884 = 0.2e1 * t1911 * t2148;
t1902 = t1905 ^ 2;
t2076 = pkin(1) ^ 2;
t2128 = -t2076 + (-0.2e1 * pkin(1) - pkin(2)) * pkin(2);
t1980 = t2010 * t2056 - t2305;
t2225 = t1980 * t2066;
t1968 = t2010 * t2050 + t2302;
t1933 = -t1968 * t2030 - t2033 * t2107;
t2245 = t1933 * t2066;
t1932 = t1968 * t2033 - t2030 * t2107;
t2246 = t1932 * t2066;
t2103 = -t1998 * t2161 + (t1911 * t2188 + ((-t2065 + t2128) * t2255 + t2107 * t2308) * t1998) * t2148 + t1893 * t2249 - t2048 * t2246 - t2047 * t2245 - t2046 * t2225;
t2294 = t1839 * pkin(1);
t2100 = t2103 + t2294;
t2121 = t1905 * t2137;
t2138 = t1905 * t2233;
t2140 = t1879 * t2234;
t2258 = t1902 * t2067;
t2152 = t1999 * t2258;
t2155 = (0.4e1 * t1897 - 0.2e1 * t1908) * t2234;
t2179 = -t2065 + t2076;
t2204 = t2049 * t2055;
t2270 = t1839 * qJ(3,3);
t2314 = pkin(1) * t1854;
t2317 = 0.2e1 * t2039 - 0.1e1;
t2323 = 0.2e1 * t2055;
t2324 = -0.2e1 * qJ(3,3);
t2327 = -2 * MDP(5);
t2334 = -(MDP(10) * t2049 - MDP(2)) * t2182 + t2182 * t2055 * MDP(9) + t1854 * MDP(1) + ((0.4e1 * t2138 + 0.2e1 * t2314) * t2039 + ((t1851 - t2155) * t2049 + t2182) * t2055 - 0.2e1 * t2138) * MDP(11) + ((-pkin(1) * t2152 + t1884 + t2270) * t2055 + (-t1902 * t2215 - t2100) * t2049 - t2352) * MDP(12) + ((t2155 - 0.2e1 * t2276) * t2039 - 0.2e1 * t2140 + t1851 + ((0.2e1 * t2138 + t2314) * t2323 + t2182) * t2049) * MDP(13) + ((0.4e1 * (t1897 - t1908 / 0.2e1) * qJ(3,3) * t2234 + t2179 * t1854) * t2039 + (0.2e1 * (qJ(3,3) * t2138 + (-t2140 + t2276) * pkin(1)) * t2049 + t2182 * pkin(1)) * t2055 + (((t2291 / 0.2e1 - t2025 / 0.2e1) * t2056 - t2036 / 0.2e1) * t2049 - t2276 / 0.2e1 + t2140) * t2324) * MDP(14) + t2352 * MDP(3) + (t1854 * t2049 + t2121 * t2323) * t2049 * MDP(4) + (t1839 * t2049 + t2055 * t2152) * MDP(6) + (t1839 * t2055 - t2049 * t2152) * MDP(7) - (t1854 * t2204 + t2317 * t2121) * t2327;
t2002 = 0.1e1 / t2011 ^ 2;
t2068 = qJ(3,2) ^ 2;
t2231 = t1954 * t2001;
t2135 = t2040 * t2231;
t2229 = t1954 * t2069;
t1899 = pkin(1) * t2147;
t1909 = t1912 * t2069;
t1881 = t1899 - t1909;
t1867 = pkin(2) * t2147 + t1881;
t2260 = t1867 * t2051;
t2280 = qJ(3,2) * t2057;
t2160 = ((qJ(3,2) + t2064) * (-qJ(3,2) + t2064) * t2135 + 0.2e1 * (t1954 * t2199 + t2253 * t2320) * t2001 * t2280 + pkin(4) * t2260 + (t2068 + t2074) * t2231) * t2229;
t2252 = t1906 * t2070;
t2166 = (-(t1945 + t1867) * t2191 + (pkin(4) * t2135 - t2260) * qJ(3,2)) * t2252;
t1840 = (t2057 * t2160 - t2166) * t2002 + t2348;
t2003 = t2001 * t2002;
t2108 = -t2199 + t2280;
t2213 = t2002 * t2069;
t2133 = t1954 * t2213;
t2187 = t2064 * t2069;
t2206 = t2034 * t2058;
t2209 = t2031 * t2058;
t2230 = t1954 * t2002;
t1855 = -t2003 * t2108 * t1906 * t2229 + (-t1909 * t2230 - t1912 * t2133) * t2051 + (t2048 * t2206 - t2047 * t2209 - t2052 * t2046 - (-t2307 + (-t2051 * t2187 + t2057) * t1906) * t2230) * t2001;
t2282 = qJ(3,2) * t1855;
t1852 = 0.2e1 * t2282;
t2146 = t2001 * t2252;
t1885 = 0.2e1 * t1912 * t2146;
t1903 = t1906 ^ 2;
t1981 = t2011 * t2058 - t2304;
t2224 = t1981 * t2069;
t1969 = t2011 * t2052 + t2301;
t1935 = -t1969 * t2031 - t2034 * t2108;
t2243 = t1935 * t2069;
t1934 = t1969 * t2034 - t2031 * t2108;
t2244 = t1934 * t2069;
t2102 = -t2001 * t2160 + (t1912 * t2187 + ((-t2068 + t2128) * t2253 + t2108 * t2307) * t2001) * t2146 + t1894 * t2248 - t2048 * t2244 - t2047 * t2243 - t2046 * t2224;
t2293 = t1840 * pkin(1);
t2099 = t2102 + t2293;
t2120 = t1906 * t2133;
t2134 = t1906 * t2230;
t2136 = t1881 * t2231;
t2257 = t1903 * t2070;
t2151 = t2002 * t2257;
t2154 = (0.4e1 * t1899 - 0.2e1 * t1909) * t2231;
t2178 = -t2068 + t2076;
t2200 = t2051 * t2057;
t2269 = t1840 * qJ(3,2);
t2313 = pkin(1) * t1855;
t2316 = 0.2e1 * t2040 - 0.1e1;
t2322 = 0.2e1 * t2057;
t2325 = -0.2e1 * qJ(3,2);
t2333 = -(MDP(10) * t2051 - MDP(2)) * t2181 + t2181 * t2057 * MDP(9) + t1855 * MDP(1) + ((0.4e1 * t2134 + 0.2e1 * t2313) * t2040 + ((t1852 - t2154) * t2051 + t2181) * t2057 - 0.2e1 * t2134) * MDP(11) + ((-pkin(1) * t2151 + t1885 + t2269) * t2057 + (-t1903 * t2213 - t2099) * t2051 - t2351) * MDP(12) + ((t2154 - 0.2e1 * t2282) * t2040 - 0.2e1 * t2136 + t1852 + ((0.2e1 * t2134 + t2313) * t2322 + t2181) * t2051) * MDP(13) + ((0.4e1 * (t1899 - t1909 / 0.2e1) * qJ(3,2) * t2231 + t2178 * t1855) * t2040 + (0.2e1 * (qJ(3,2) * t2134 + (-t2136 + t2282) * pkin(1)) * t2051 + t2181 * pkin(1)) * t2057 + (((t2290 / 0.2e1 - t2026 / 0.2e1) * t2058 - t2037 / 0.2e1) * t2051 - t2282 / 0.2e1 + t2136) * t2325) * MDP(14) + t2351 * MDP(3) + (t1855 * t2051 + t2120 * t2322) * t2051 * MDP(4) + (t1840 * t2051 + t2057 * t2151) * MDP(6) + (t1840 * t2057 - t2051 * t2151) * MDP(7) - (t1855 * t2200 + t2316 * t2120) * t2327;
t2005 = 0.1e1 / t2012 ^ 2;
t2071 = qJ(3,1) ^ 2;
t2228 = t1955 * t2004;
t2131 = t2041 * t2228;
t2226 = t1955 * t2072;
t1901 = pkin(1) * t2145;
t1910 = t1913 * t2072;
t1883 = t1901 - t1910;
t1868 = pkin(2) * t2145 + t1883;
t2259 = t1868 * t2053;
t2286 = qJ(3,1) * t2059;
t2159 = ((qJ(3,1) + t2064) * (-qJ(3,1) + t2064) * t2131 + 0.2e1 * (t1955 * t2195 + t2251 * t2320) * t2004 * t2286 + pkin(4) * t2259 + (t2071 + t2074) * t2228) * t2226;
t2250 = t1907 * t2073;
t2165 = (-(t1946 + t1868) * t2190 + (pkin(4) * t2131 - t2259) * qJ(3,1)) * t2250;
t1841 = (t2059 * t2159 - t2165) * t2005 + t2349;
t2006 = t2004 * t2005;
t2109 = -t2195 + t2286;
t2211 = t2005 * t2072;
t2129 = t1955 * t2211;
t2186 = t2064 * t2072;
t2205 = t2035 * t2060;
t2208 = t2032 * t2060;
t2227 = t1955 * t2005;
t1856 = -t2006 * t2109 * t1907 * t2226 + (-t1910 * t2227 - t1913 * t2129) * t2053 + (t2048 * t2205 - t2047 * t2208 - t2054 * t2046 - (-t2306 + (-t2053 * t2186 + t2059) * t1907) * t2227) * t2004;
t2288 = qJ(3,1) * t1856;
t1853 = 0.2e1 * t2288;
t2144 = t2004 * t2250;
t1886 = 0.2e1 * t1913 * t2144;
t1904 = t1907 ^ 2;
t1982 = t2012 * t2060 - t2303;
t2223 = t1982 * t2072;
t1970 = t2012 * t2054 + t2300;
t1937 = -t1970 * t2032 - t2035 * t2109;
t2241 = t1937 * t2072;
t1936 = t1970 * t2035 - t2032 * t2109;
t2242 = t1936 * t2072;
t2101 = -t2004 * t2159 + (t1913 * t2186 + ((-t2071 + t2128) * t2251 + t2109 * t2306) * t2004) * t2144 + t1895 * t2247 - t2048 * t2242 - t2047 * t2241 - t2046 * t2223;
t2292 = t1841 * pkin(1);
t2098 = t2101 + t2292;
t2119 = t1907 * t2129;
t2130 = t1907 * t2227;
t2132 = t1883 * t2228;
t2256 = t1904 * t2073;
t2150 = t2005 * t2256;
t2153 = (0.4e1 * t1901 - 0.2e1 * t1910) * t2228;
t2177 = -t2071 + t2076;
t2196 = t2053 * t2059;
t2268 = t1841 * qJ(3,1);
t2312 = pkin(1) * t1856;
t2315 = 0.2e1 * t2041 - 0.1e1;
t2321 = 0.2e1 * t2059;
t2326 = -0.2e1 * qJ(3,1);
t2332 = -(MDP(10) * t2053 - MDP(2)) * t2180 + t2180 * t2059 * MDP(9) + t1856 * MDP(1) + ((0.4e1 * t2130 + 0.2e1 * t2312) * t2041 + ((t1853 - t2153) * t2053 + t2180) * t2059 - 0.2e1 * t2130) * MDP(11) + ((-pkin(1) * t2150 + t1886 + t2268) * t2059 + (-t1904 * t2211 - t2098) * t2053 - t2350) * MDP(12) + ((t2153 - 0.2e1 * t2288) * t2041 - 0.2e1 * t2132 + t1853 + ((0.2e1 * t2130 + t2312) * t2321 + t2180) * t2053) * MDP(13) + ((0.4e1 * (t1901 - t1910 / 0.2e1) * qJ(3,1) * t2228 + t2177 * t1856) * t2041 + (0.2e1 * (qJ(3,1) * t2130 + (-t2132 + t2288) * pkin(1)) * t2053 + t2180 * pkin(1)) * t2059 + (((t2289 / 0.2e1 - t2027 / 0.2e1) * t2060 - t2038 / 0.2e1) * t2053 - t2288 / 0.2e1 + t2132) * t2326) * MDP(14) + t2350 * MDP(3) + (t1856 * t2053 + t2119 * t2321) * t2053 * MDP(4) + (t1841 * t2053 + t2059 * t2150) * MDP(6) + (t1841 * t2059 - t2053 * t2150) * MDP(7) - (t1856 * t2196 + t2315 * t2119) * t2327;
t2328 = 0.2e1 * pkin(1);
t2319 = pkin(1) * g(1);
t2318 = pkin(1) * g(2);
t2311 = pkin(1) * t2050;
t2310 = pkin(1) * t2052;
t2309 = pkin(1) * t2054;
t1950 = t1953 ^ 2;
t2240 = t1950 * t1999;
t2239 = t1950 * t2049;
t1951 = t1954 ^ 2;
t2238 = t1951 * t2002;
t2237 = t1951 * t2051;
t1952 = t1955 ^ 2;
t2236 = t1952 * t2005;
t2235 = t1952 * t2053;
t2222 = t1989 * t2049;
t2220 = t1990 * t2051;
t2218 = t1991 * t2053;
t2202 = t2049 * t2066;
t2198 = t2051 * t2069;
t2194 = t2053 * t2072;
t1947 = pkin(1) * t2240;
t2185 = t1947 + t1884;
t1948 = pkin(1) * t2238;
t2184 = t1948 + t1885;
t1949 = pkin(1) * t2236;
t2183 = t1949 + t1886;
t2164 = t1854 * t2202;
t2163 = t1855 * t2198;
t2162 = t1856 * t2194;
t2143 = t1999 * t2239;
t2142 = t2002 * t2237;
t2141 = t2005 * t2235;
t1914 = t2317 * t2240;
t1915 = t2316 * t2238;
t1916 = t2315 * t2236;
t2118 = t1950 * t2000 * t2202;
t2117 = t1951 * t2003 * t2198;
t2116 = t1952 * t2006 * t2194;
t2115 = -0.2e1 * t2039 * t1947;
t2114 = -0.2e1 * t2040 * t1948;
t2113 = -0.2e1 * t2041 * t1949;
t2112 = MDP(12) * (-pkin(1) * t2049 + t2274) + MDP(6) * t2049;
t2111 = MDP(12) * (-pkin(1) * t2051 + t2280) + MDP(6) * t2051;
t2110 = MDP(12) * (-pkin(1) * t2053 + t2286) + MDP(6) * t2053;
t2106 = t2055 * t2118;
t2105 = t2057 * t2117;
t2104 = t2059 * t2116;
t2094 = MDP(10) * (t2055 * t2352 + t2222) + MDP(11) * ((t2143 * t2328 - t1989) * t2055 + t1956 + 0.2e1 * t2294 - qJ(3,3) * t1914 + t2103) + MDP(13) * (t2115 + (t2143 * t2324 - t2352) * t2055 - t2222 + 0.2e1 * t2270 + t2185) + MDP(14) * (qJ(3,3) * t2115 + ((-g(1) * t2275 - t2318) * t2033 + (g(2) * t2275 - t2319) * t2030 - g(3) * t2273 + t2179 * t2143) * t2055 + ((g(1) * t2311 - g(2) * qJ(3,3)) * t2033 + (-g(1) * qJ(3,3) - g(2) * t2311) * t2030 + pkin(1) * t2297) * t2049 + t1839 * t2065 + t2185 * qJ(3,3) + t2100 * pkin(1)) - MDP(5) * t1914 + MDP(8) * t1839 + MDP(9) * t2346;
t2093 = MDP(10) * (t2057 * t2351 + t2220) + MDP(11) * ((t2142 * t2328 - t1990) * t2057 + t1957 + 0.2e1 * t2293 - qJ(3,2) * t1915 + t2102) + MDP(13) * (t2114 + (t2142 * t2325 - t2351) * t2057 - t2220 + 0.2e1 * t2269 + t2184) + MDP(14) * (qJ(3,2) * t2114 + ((-g(1) * t2281 - t2318) * t2034 + (g(2) * t2281 - t2319) * t2031 - g(3) * t2279 + t2178 * t2142) * t2057 + ((g(1) * t2310 - g(2) * qJ(3,2)) * t2034 + (-g(1) * qJ(3,2) - g(2) * t2310) * t2031 + pkin(1) * t2296) * t2051 + t1840 * t2068 + t2184 * qJ(3,2) + t2099 * pkin(1)) - MDP(5) * t1915 + MDP(8) * t1840 + MDP(9) * t2345;
t2092 = MDP(10) * (t2059 * t2350 + t2218) + MDP(11) * ((t2141 * t2328 - t1991) * t2059 + t1958 + 0.2e1 * t2292 - qJ(3,1) * t1916 + t2101) + MDP(13) * (t2113 + (t2141 * t2326 - t2350) * t2059 - t2218 + 0.2e1 * t2268 + t2183) + MDP(14) * (qJ(3,1) * t2113 + ((-g(1) * t2287 - t2318) * t2035 + (g(2) * t2287 - t2319) * t2032 - g(3) * t2285 + t2177 * t2141) * t2059 + ((g(1) * t2309 - g(2) * qJ(3,1)) * t2035 + (-g(1) * qJ(3,1) - g(2) * t2309) * t2032 + pkin(1) * t2295) * t2053 + t1841 * t2071 + t2183 * qJ(3,1) + t2098 * pkin(1)) - MDP(5) * t1916 + MDP(8) * t1841 + MDP(9) * t2344;
t2088 = t2066 * ((MDP(7) * t2055 + t2112) * t1854 + t2094);
t2087 = t2069 * ((MDP(7) * t2057 + t2111) * t1855 + t2093);
t2086 = t2072 * ((MDP(7) * t2059 + t2110) * t1856 + t2092);
t1889 = (t1952 * t2041 - t1952 - t2256) * t2005;
t1888 = (t1951 * t2040 - t1951 - t2257) * t2002;
t1887 = (t1950 * t2039 - t1950 - t2258) * t1999;
t1835 = (t2165 + (-t2159 - t2235) * t2059) * t2005 - t2349;
t1834 = (t2166 + (-t2160 - t2237) * t2057) * t2002 - t2348;
t1833 = (t2167 + (-t2161 - t2239) * t2055) * t1999 - t2347;
t1817 = (-t1904 * t2072 + (-pkin(1) * t2196 + qJ(3,1) * t2041 - qJ(3,1)) * t1952) * t2005 - t2098 - t2344;
t1816 = (-t1903 * t2069 + (-pkin(1) * t2200 + qJ(3,2) * t2040 - qJ(3,2)) * t1951) * t2002 - t2099 - t2345;
t1815 = (-t1902 * t2066 + (-pkin(1) * t2204 + qJ(3,3) * t2039 - qJ(3,3)) * t1950) * t1999 - t2100 - t2346;
t1 = [(-t1917 * t2106 - t1918 * t2105 - t1919 * t2104) * MDP(4) + (t1833 * t2246 + t1834 * t2244 + t1835 * t2242) * MDP(11) + (t1932 * t2164 + t1934 * t2163 + t1936 * t2162) * MDP(12) + (t1887 * t2246 + t1888 * t2244 + t1889 * t2242) * MDP(13) + (t1815 * t2246 + t1816 * t2244 + t1817 * t2242) * MDP(14) + (t2048 - g(1)) * MDP(15) + (t1919 * t2086 + t2332 * t2205) * t2004 + (t1918 * t2087 + t2333 * t2206) * t2001 + (t1917 * t2088 + t2334 * t2207) * t1998; (-t1920 * t2106 - t1921 * t2105 - t1922 * t2104) * MDP(4) + (t1833 * t2245 + t1834 * t2243 + t1835 * t2241) * MDP(11) + (t1933 * t2164 + t1935 * t2163 + t1937 * t2162) * MDP(12) + (t1887 * t2245 + t1888 * t2243 + t1889 * t2241) * MDP(13) + (t1815 * t2245 + t1816 * t2243 + t1817 * t2241) * MDP(14) + (t2047 - g(2)) * MDP(15) + (t1922 * t2086 - t2208 * t2332) * t2004 + (t1921 * t2087 - t2209 * t2333) * t2001 + (t1920 * t2088 - t2210 * t2334) * t1998; (-t1971 * t2039 * t2118 - t1972 * t2040 * t2117 - t1973 * t2041 * t2116) * MDP(4) + (t1833 * t2225 + t1834 * t2224 + t1835 * t2223) * MDP(11) + (t1980 * t2164 + t1981 * t2163 + t1982 * t2162) * MDP(12) + (t1887 * t2225 + t1888 * t2224 + t1889 * t2223) * MDP(13) + (t1815 * t2225 + t1816 * t2224 + t1817 * t2223) * MDP(14) + (t2046 - g(3)) * MDP(15) + ((MDP(7) * t1856 * t2041 + (t1856 * t2110 + t2092) * t2059) * t2072 * t1973 - t2332 * t2054) * t2004 + ((MDP(7) * t1855 * t2040 + (t1855 * t2111 + t2093) * t2057) * t2069 * t1972 - t2333 * t2052) * t2001 + ((MDP(7) * t1854 * t2039 + (t1854 * t2112 + t2094) * t2055) * t2066 * t1971 - t2334 * t2050) * t1998;];
tauX  = t1;
