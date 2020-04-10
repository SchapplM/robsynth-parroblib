% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P4PRRR1G1P1A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [11x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P4PRRR1G1P1A0_convert_par2_MPV_fixb.m

% Output:
% taucX [4x1]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-02 20:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P4PRRR1G1P1A0_coriolisvec_para_pf_mdp(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(7,1),zeros(11,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRR1G1P1A0_coriolisvec_para_pf_mdp: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRR1G1P1A0_coriolisvec_para_pf_mdp: xDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRR1G1P1A0_coriolisvec_para_pf_mdp: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P4PRRR1G1P1A0_coriolisvec_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRR1G1P1A0_coriolisvec_para_pf_mdp: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRR1G1P1A0_coriolisvec_para_pf_mdp: Koppelpunkt has to be [4x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [11 1]), ...
  'P4PRRR1G1P1A0_coriolisvec_para_pf_mdp: MDP has to be [11x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-02 20:15:20
% EndTime: 2020-03-02 20:15:32
% DurationCPUTime: 13.26s
% Computational Cost: add. (46877->293), mult. (40660->539), div. (5824->11), fcn. (35532->34), ass. (0->264)
t2173 = xP(4);
t2151 = sin(t2173);
t2152 = cos(t2173);
t2174 = koppelP(4,2);
t2178 = koppelP(4,1);
t2103 = t2151 * t2178 + t2152 * t2174;
t2170 = xDP(4);
t2172 = xDP(1);
t2091 = t2103 * t2170 - t2172;
t2107 = -t2151 * t2174 + t2152 * t2178;
t2171 = xDP(2);
t2095 = t2107 * t2170 + t2171;
t2153 = pkin(7) + qJ(2,4);
t2133 = qJ(3,4) + t2153;
t2119 = sin(t2133);
t2120 = cos(t2133);
t2158 = legFrame(4,3);
t2143 = sin(t2158);
t2147 = cos(t2158);
t2028 = (-t2091 * t2147 + t2095 * t2143) * t2120 + t2119 * (t2091 * t2143 + t2095 * t2147);
t2131 = sin(t2153);
t2132 = cos(t2153);
t2067 = t2119 * t2132 - t2120 * t2131;
t2333 = 0.1e1 / t2067;
t2311 = t2028 * t2333;
t2175 = koppelP(3,2);
t2179 = koppelP(3,1);
t2104 = t2151 * t2179 + t2152 * t2175;
t2092 = t2104 * t2170 - t2172;
t2108 = -t2151 * t2175 + t2152 * t2179;
t2096 = t2108 * t2170 + t2171;
t2154 = pkin(7) + qJ(2,3);
t2140 = qJ(3,3) + t2154;
t2125 = sin(t2140);
t2128 = cos(t2140);
t2159 = legFrame(3,3);
t2144 = sin(t2159);
t2148 = cos(t2159);
t2032 = (-t2092 * t2148 + t2096 * t2144) * t2128 + t2125 * (t2092 * t2144 + t2096 * t2148);
t2134 = sin(t2154);
t2137 = cos(t2154);
t2080 = t2125 * t2137 - t2128 * t2134;
t2332 = 0.1e1 / t2080;
t2309 = t2032 * t2332;
t2176 = koppelP(2,2);
t2180 = koppelP(2,1);
t2105 = t2151 * t2180 + t2152 * t2176;
t2093 = t2105 * t2170 - t2172;
t2109 = -t2151 * t2176 + t2152 * t2180;
t2097 = t2109 * t2170 + t2171;
t2155 = pkin(7) + qJ(2,2);
t2141 = qJ(3,2) + t2155;
t2126 = sin(t2141);
t2129 = cos(t2141);
t2160 = legFrame(2,3);
t2145 = sin(t2160);
t2149 = cos(t2160);
t2033 = (-t2093 * t2149 + t2097 * t2145) * t2129 + t2126 * (t2093 * t2145 + t2097 * t2149);
t2135 = sin(t2155);
t2138 = cos(t2155);
t2081 = t2126 * t2138 - t2129 * t2135;
t2331 = 0.1e1 / t2081;
t2307 = t2033 * t2331;
t2177 = koppelP(1,2);
t2181 = koppelP(1,1);
t2106 = t2151 * t2181 + t2152 * t2177;
t2094 = t2106 * t2170 - t2172;
t2110 = -t2151 * t2177 + t2152 * t2181;
t2098 = t2110 * t2170 + t2171;
t2156 = pkin(7) + qJ(2,1);
t2142 = qJ(3,1) + t2156;
t2127 = sin(t2142);
t2130 = cos(t2142);
t2161 = legFrame(1,3);
t2146 = sin(t2161);
t2150 = cos(t2161);
t2034 = (-t2094 * t2150 + t2098 * t2146) * t2130 + t2127 * (t2094 * t2146 + t2098 * t2150);
t2136 = sin(t2156);
t2139 = cos(t2156);
t2082 = t2127 * t2139 - t2130 * t2136;
t2330 = 0.1e1 / t2082;
t2305 = t2034 * t2330;
t2341 = -t2127 * t2146 + t2150 * t2130;
t2340 = -t2126 * t2145 + t2149 * t2129;
t2339 = -t2125 * t2144 + t2148 * t2128;
t2338 = -t2119 * t2143 + t2147 * t2120;
t2044 = -pkin(2) * (t2131 * t2143 - t2132 * t2147) + t2338 * pkin(3);
t2183 = 0.1e1 / pkin(3);
t2184 = 0.1e1 / pkin(2);
t2278 = t2183 * t2184;
t2035 = t2044 * t2091 * t2333 * t2278;
t2084 = t2119 * t2147 + t2120 * t2143;
t2043 = pkin(2) * (t2131 * t2147 + t2132 * t2143) + t2084 * pkin(3);
t2249 = t2095 * t2043 * t2183;
t2299 = t2333 * t2184;
t2009 = t2035 + (t2028 - t2249) * t2299;
t2221 = t2119 * t2131 + t2120 * t2132;
t2185 = 0.1e1 / pkin(2) ^ 2;
t2256 = t2185 * t2311;
t2241 = (pkin(3) * t2009 + t2221 * t2311) * t2256;
t2049 = -pkin(2) * (t2134 * t2144 - t2137 * t2148) + t2339 * pkin(3);
t2036 = t2049 * t2092 * t2332 * t2278;
t2088 = t2125 * t2148 + t2128 * t2144;
t2046 = pkin(2) * (t2134 * t2148 + t2137 * t2144) + t2088 * pkin(3);
t2248 = t2096 * t2046 * t2183;
t2296 = t2332 * t2184;
t2014 = t2036 + (t2032 - t2248) * t2296;
t2220 = t2125 * t2134 + t2128 * t2137;
t2252 = t2185 * t2309;
t2240 = (pkin(3) * t2014 + t2220 * t2309) * t2252;
t2051 = -pkin(2) * (t2135 * t2145 - t2138 * t2149) + t2340 * pkin(3);
t2037 = t2051 * t2093 * t2331 * t2278;
t2089 = t2126 * t2149 + t2129 * t2145;
t2047 = pkin(2) * (t2135 * t2149 + t2138 * t2145) + t2089 * pkin(3);
t2247 = t2097 * t2047 * t2183;
t2295 = t2331 * t2184;
t2016 = t2037 + (t2033 - t2247) * t2295;
t2219 = t2126 * t2135 + t2129 * t2138;
t2251 = t2185 * t2307;
t2239 = (pkin(3) * t2016 + t2219 * t2307) * t2251;
t2053 = -pkin(2) * (t2136 * t2146 - t2139 * t2150) + t2341 * pkin(3);
t2038 = t2053 * t2094 * t2330 * t2278;
t2090 = t2127 * t2150 + t2130 * t2146;
t2048 = pkin(2) * (t2136 * t2150 + t2139 * t2146) + t2090 * pkin(3);
t2246 = t2098 * t2048 * t2183;
t2294 = t2330 * t2184;
t2018 = t2038 + (t2034 - t2246) * t2294;
t2218 = t2127 * t2136 + t2130 * t2139;
t2250 = t2185 * t2305;
t2238 = (pkin(3) * t2018 + t2218 * t2305) * t2250;
t2257 = t2028 ^ 2 / t2067 ^ 2 * t2333;
t2255 = t2032 ^ 2 / t2080 ^ 2 * t2332;
t2254 = t2033 ^ 2 / t2081 ^ 2 * t2331;
t2253 = t2034 ^ 2 / t2082 ^ 2 * t2330;
t2022 = -t2246 * t2294 + t2038;
t2262 = t2018 * t2022 * t2184;
t2234 = t2330 * t2262;
t2006 = (t2218 * pkin(2) + pkin(3)) * t2234;
t2182 = pkin(3) ^ 2;
t2157 = t2170 ^ 2;
t2279 = t2157 * t2184;
t2328 = 0.2e1 * pkin(3);
t2198 = (-(-t2018 * t2182 + (-t2305 - t2218 * (t2038 / 0.2e1 + (t2034 - t2246 / 0.2e1) * t2294) * t2328) * pkin(2)) * t2250 + (-t2048 * t2106 - t2053 * t2110) * t2279) * t2183;
t2206 = (-t2090 * t2106 - t2110 * t2341) * t2279;
t2202 = t2330 * t2206;
t2242 = pkin(3) * t2262;
t1987 = -t2006 + 0.2e1 * t2202 - (t2198 - 0.2e1 * t2238 - 0.2e1 * t2242) * t2330;
t2166 = sin(qJ(3,1));
t2169 = cos(qJ(3,1));
t2316 = (t2038 + (0.2e1 * t2034 - t2246) * t2294) * t2022;
t2337 = (-t1987 * t2169 + t2166 * t2316) * MDP(6) + (t1987 * t2166 + t2169 * t2316) * MDP(7);
t2021 = -t2247 * t2295 + t2037;
t2263 = t2016 * t2021 * t2184;
t2235 = t2331 * t2263;
t2005 = (t2219 * pkin(2) + pkin(3)) * t2235;
t2199 = (-(-t2016 * t2182 + (-t2307 - t2219 * (t2037 / 0.2e1 + (t2033 - t2247 / 0.2e1) * t2295) * t2328) * pkin(2)) * t2251 + (-t2047 * t2105 - t2051 * t2109) * t2279) * t2183;
t2207 = (-t2089 * t2105 - t2109 * t2340) * t2279;
t2203 = t2331 * t2207;
t2243 = pkin(3) * t2263;
t1986 = -t2005 + 0.2e1 * t2203 - (t2199 - 0.2e1 * t2239 - 0.2e1 * t2243) * t2331;
t2165 = sin(qJ(3,2));
t2168 = cos(qJ(3,2));
t2317 = (t2037 + (0.2e1 * t2033 - t2247) * t2295) * t2021;
t2336 = (-t1986 * t2168 + t2165 * t2317) * MDP(6) + (t1986 * t2165 + t2168 * t2317) * MDP(7);
t2020 = -t2248 * t2296 + t2036;
t2264 = t2014 * t2020 * t2184;
t2236 = t2332 * t2264;
t2004 = (t2220 * pkin(2) + pkin(3)) * t2236;
t2200 = (-(-t2014 * t2182 + (-t2309 - t2220 * (t2036 / 0.2e1 + (t2032 - t2248 / 0.2e1) * t2296) * t2328) * pkin(2)) * t2252 + (-t2046 * t2104 - t2049 * t2108) * t2279) * t2183;
t2208 = (-t2088 * t2104 - t2108 * t2339) * t2279;
t2204 = t2332 * t2208;
t2244 = pkin(3) * t2264;
t1985 = -t2004 + 0.2e1 * t2204 - (t2200 - 0.2e1 * t2240 - 0.2e1 * t2244) * t2332;
t2164 = sin(qJ(3,3));
t2167 = cos(qJ(3,3));
t2318 = (t2036 + (0.2e1 * t2032 - t2248) * t2296) * t2020;
t2335 = (-t1985 * t2167 + t2164 * t2318) * MDP(6) + (t1985 * t2164 + t2167 * t2318) * MDP(7);
t2019 = -t2249 * t2299 + t2035;
t2265 = t2009 * t2019 * t2184;
t2237 = t2333 * t2265;
t2003 = (t2221 * pkin(2) + pkin(3)) * t2237;
t2201 = (-(-t2009 * t2182 + (-t2311 - t2221 * (t2035 / 0.2e1 + (t2028 - t2249 / 0.2e1) * t2299) * t2328) * pkin(2)) * t2256 + (-t2043 * t2103 - t2044 * t2107) * t2279) * t2183;
t2209 = (-t2084 * t2103 - t2107 * t2338) * t2279;
t2205 = t2333 * t2209;
t2245 = pkin(3) * t2265;
t1983 = -t2003 + 0.2e1 * t2205 - (t2201 - 0.2e1 * t2241 - 0.2e1 * t2245) * t2333;
t2162 = sin(qJ(3,4));
t2163 = cos(qJ(3,4));
t2319 = (t2035 + (0.2e1 * t2028 - t2249) * t2299) * t2019;
t2334 = (-t1983 * t2163 + t2162 * t2319) * MDP(6) + (t1983 * t2162 + t2163 * t2319) * MDP(7);
t2329 = MDP(5) * t2184;
t1984 = -t2003 + t2205 - (t2201 - t2241 - t2245) * t2333;
t2327 = t1984 * t2333;
t1988 = -t2004 + t2204 - (t2200 - t2240 - t2244) * t2332;
t2326 = t1988 * t2332;
t1989 = -t2005 + t2203 - (t2199 - t2239 - t2243) * t2331;
t2325 = t1989 * t2331;
t1990 = -t2006 + t2202 - (t2198 - t2238 - t2242) * t2330;
t2324 = t1990 * t2330;
t1991 = pkin(3) * t2237 + (t2209 + t2241) * t2333;
t2323 = t1991 * t2333;
t1992 = pkin(3) * t2236 + (t2208 + t2240) * t2332;
t2322 = t1992 * t2332;
t1993 = pkin(3) * t2235 + (t2207 + t2239) * t2331;
t2321 = t1993 * t2331;
t1994 = pkin(3) * t2234 + (t2206 + t2238) * t2330;
t2320 = t1994 * t2330;
t2055 = t2103 * t2147 - t2107 * t2143;
t2059 = t2103 * t2143 + t2107 * t2147;
t2225 = t2055 * t2120 - t2059 * t2119;
t2023 = pkin(2) * (t2055 * t2132 - t2059 * t2131) + t2225 * pkin(3);
t2315 = t2023 * t2333;
t2056 = t2104 * t2148 - t2108 * t2144;
t2060 = t2104 * t2144 + t2108 * t2148;
t2224 = t2056 * t2128 - t2060 * t2125;
t2024 = pkin(2) * (t2056 * t2137 - t2060 * t2134) + t2224 * pkin(3);
t2314 = t2024 * t2332;
t2057 = t2105 * t2149 - t2109 * t2145;
t2061 = t2105 * t2145 + t2109 * t2149;
t2223 = t2057 * t2129 - t2061 * t2126;
t2025 = pkin(2) * (t2057 * t2138 - t2061 * t2135) + t2223 * pkin(3);
t2313 = t2025 * t2331;
t2058 = t2106 * t2150 - t2110 * t2146;
t2062 = t2106 * t2146 + t2110 * t2150;
t2222 = t2058 * t2130 - t2062 * t2127;
t2026 = pkin(2) * (t2058 * t2139 - t2062 * t2136) + t2222 * pkin(3);
t2312 = t2026 * t2330;
t2303 = t2225 * t2333;
t2302 = t2224 * t2332;
t2301 = t2223 * t2331;
t2300 = t2222 * t2330;
t2298 = t2333 * t2338;
t2297 = t2333 * t2084;
t2293 = t2332 * t2339;
t2292 = t2332 * t2088;
t2291 = t2331 * t2340;
t2290 = t2331 * t2089;
t2289 = t2330 * t2341;
t2288 = t2330 * t2090;
t2277 = t1991 * t2315;
t2276 = t2162 * t2323;
t2275 = t2163 * t2323;
t2274 = t1992 * t2314;
t2273 = t2164 * t2322;
t2272 = t2167 * t2322;
t2271 = t1993 * t2313;
t2270 = t2165 * t2321;
t2269 = t2168 * t2321;
t2268 = t1994 * t2312;
t2267 = t2166 * t2320;
t2266 = t2169 * t2320;
t2261 = t2023 * t2257;
t2260 = t2024 * t2255;
t2259 = t2025 * t2254;
t2258 = t2026 * t2253;
t2233 = t2162 * t2257;
t2232 = t2163 * t2257;
t2231 = t2164 * t2255;
t2230 = t2167 * t2255;
t2229 = t2165 * t2254;
t2228 = t2168 * t2254;
t2227 = t2166 * t2253;
t2226 = t2169 * t2253;
t1 = [(MDP(10) * t2151 - MDP(9) * t2152) * t2157 - t2337 * t2289 - t2336 * t2291 - t2335 * t2293 - t2334 * t2298 + ((t1991 * t2298 + t1992 * t2293 + t1993 * t2291 + t1994 * t2289) * MDP(2) + (t1984 * t2298 + t1988 * t2293 + t1989 * t2291 + t1990 * t2289) * MDP(5)) * t2184 + ((-t2044 * t2275 - t2049 * t2272 - t2051 * t2269 - t2053 * t2266) * MDP(6) + (t2044 * t2276 + t2049 * t2273 + t2051 * t2270 + t2053 * t2267) * MDP(7) + (-t2044 * t2327 - t2049 * t2326 - t2051 * t2325 - t2053 * t2324) * t2329 + ((-t2044 * t2233 - t2049 * t2231 - t2051 * t2229 - t2053 * t2227) * MDP(6) + (-t2044 * t2232 - t2049 * t2230 - t2051 * t2228 - t2053 * t2226) * MDP(7)) * t2185) * t2183; (-MDP(10) * t2152 - MDP(9) * t2151) * t2157 - t2337 * t2288 - t2336 * t2290 - t2335 * t2292 - t2334 * t2297 + ((t1991 * t2297 + t1992 * t2292 + t1993 * t2290 + t1994 * t2288) * MDP(2) + (t1984 * t2297 + t1988 * t2292 + t1989 * t2290 + t1990 * t2288) * MDP(5)) * t2184 + ((-t2043 * t2275 - t2046 * t2272 - t2047 * t2269 - t2048 * t2266) * MDP(6) + (t2043 * t2276 + t2046 * t2273 + t2047 * t2270 + t2048 * t2267) * MDP(7) + (-t2043 * t2327 - t2046 * t2326 - t2047 * t2325 - t2048 * t2324) * t2329 + ((-t2043 * t2233 - t2046 * t2231 - t2047 * t2229 - t2048 * t2227) * MDP(6) + (-t2043 * t2232 - t2046 * t2230 - t2047 * t2228 - t2048 * t2226) * MDP(7)) * t2185) * t2183; 0; t2337 * t2300 + t2336 * t2301 + t2335 * t2302 + t2334 * t2303 + ((-t1991 * t2303 - t1992 * t2302 - t1993 * t2301 - t1994 * t2300) * MDP(2) + (-t1984 * t2303 - t1988 * t2302 - t1989 * t2301 - t1990 * t2300) * MDP(5)) * t2184 + ((t2163 * t2277 + t2167 * t2274 + t2168 * t2271 + t2169 * t2268) * MDP(6) + (-t2162 * t2277 - t2164 * t2274 - t2165 * t2271 - t2166 * t2268) * MDP(7) + (t1984 * t2315 + t1988 * t2314 + t1989 * t2313 + t1990 * t2312) * t2329 + ((t2162 * t2261 + t2164 * t2260 + t2165 * t2259 + t2166 * t2258) * MDP(6) + (t2163 * t2261 + t2167 * t2260 + t2168 * t2259 + t2169 * t2258) * MDP(7)) * t2185) * t2183;];
taucX  = t1;
