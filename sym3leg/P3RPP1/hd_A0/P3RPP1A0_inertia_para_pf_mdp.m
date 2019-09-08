% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RPP1G1P1A0
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
%   pkin=[a2,a3,d1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [13x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPP1G1P1A0_convert_par2_MPV_fixb.m

% Output:
% MMX [3x3]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:53
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P3RPP1G1P1A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,1),zeros(13,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPP1G1P1A0_inertia_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPP1G1P1A0_inertia_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RPP1G1P1A0_inertia_para_pf_mdp: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPP1G1P1A0_inertia_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPP1G1P1A0_inertia_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'P3RPP1G1P1A0_inertia_para_pf_mdp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:52:54
% EndTime: 2019-05-03 14:52:58
% DurationCPUTime: 4.28s
% Computational Cost: add. (10503->423), mult. (14743->658), div. (642->6), fcn. (5450->14), ass. (0->306)
t2383 = MDP(4) - MDP(8);
t2254 = (qJ(2,3) ^ 2);
t2388 = 2 * pkin(1);
t2234 = qJ(3,3) * t2388;
t2265 = pkin(1) ^ 2;
t2319 = qJ(3,3) ^ 2 + t2265;
t2280 = t2234 + t2319;
t2211 = t2254 + t2280;
t2208 = 1 + t2211;
t2196 = 1 / t2208;
t2397 = MDP(6) * t2196;
t2256 = (qJ(2,2) ^ 2);
t2235 = qJ(3,2) * t2388;
t2318 = qJ(3,2) ^ 2 + t2265;
t2279 = t2235 + t2318;
t2212 = t2256 + t2279;
t2209 = 1 + t2212;
t2198 = 1 / t2209;
t2396 = MDP(6) * t2198;
t2258 = (qJ(2,1) ^ 2);
t2236 = qJ(3,1) * t2388;
t2317 = qJ(3,1) ^ 2 + t2265;
t2278 = t2236 + t2317;
t2213 = t2258 + t2278;
t2210 = 1 + t2213;
t2200 = 1 / t2210;
t2395 = MDP(6) * t2200;
t2394 = MDP(9) * t2196;
t2393 = MDP(9) * t2198;
t2392 = MDP(9) * t2200;
t2259 = koppelP(3,2);
t2316 = 1 + t2319;
t2262 = koppelP(3,1);
t2371 = (qJ(2,3) * t2262);
t2391 = t2259 * t2316 + pkin(1) * t2371 + (t2259 * t2388 + t2371) * qJ(3,3);
t2260 = koppelP(2,2);
t2315 = 1 + t2318;
t2263 = koppelP(2,1);
t2374 = (qJ(2,2) * t2263);
t2390 = t2260 * t2315 + pkin(1) * t2374 + (t2260 * t2388 + t2374) * qJ(3,2);
t2261 = koppelP(1,2);
t2314 = 1 + t2317;
t2264 = koppelP(1,1);
t2377 = (qJ(2,1) * t2264);
t2389 = t2261 * t2314 + pkin(1) * t2377 + (t2261 * t2388 + t2377) * qJ(3,1);
t2387 = 2 * MDP(5);
t2245 = pkin(1) + qJ(3,1);
t2244 = pkin(1) + qJ(3,2);
t2243 = pkin(1) + qJ(3,3);
t2249 = cos(qJ(1,3));
t2246 = sin(qJ(1,3));
t2334 = t2243 * t2246;
t2190 = qJ(2,3) * t2249 - t2334;
t2333 = t2243 * t2249;
t2193 = qJ(2,3) * t2246 + t2333;
t2240 = legFrame(3,3);
t2223 = sin(t2240);
t2226 = cos(t2240);
t2150 = t2190 * t2226 - t2193 * t2223;
t2153 = t2190 * t2223 + t2193 * t2226;
t2252 = xP(3);
t2232 = sin(t2252);
t2233 = cos(t2252);
t2177 = -t2232 * t2262 - t2233 * t2259;
t2180 = -t2232 * t2259 + t2233 * t2262;
t2051 = (t2150 * t2177 + t2153 * t2180) * t2196;
t2386 = t2051 * pkin(1);
t2250 = cos(qJ(1,2));
t2247 = sin(qJ(1,2));
t2330 = t2244 * t2247;
t2191 = qJ(2,2) * t2250 - t2330;
t2329 = t2244 * t2250;
t2194 = qJ(2,2) * t2247 + t2329;
t2241 = legFrame(2,3);
t2224 = sin(t2241);
t2227 = cos(t2241);
t2151 = t2191 * t2227 - t2194 * t2224;
t2154 = t2191 * t2224 + t2194 * t2227;
t2178 = -t2232 * t2263 - t2233 * t2260;
t2181 = -t2232 * t2260 + t2233 * t2263;
t2052 = (t2151 * t2178 + t2154 * t2181) * t2198;
t2385 = t2052 * pkin(1);
t2251 = cos(qJ(1,1));
t2248 = sin(qJ(1,1));
t2326 = t2245 * t2248;
t2192 = qJ(2,1) * t2251 - t2326;
t2325 = t2245 * t2251;
t2195 = qJ(2,1) * t2248 + t2325;
t2242 = legFrame(1,3);
t2225 = sin(t2242);
t2228 = cos(t2242);
t2152 = t2192 * t2228 - t2195 * t2225;
t2155 = t2192 * t2225 + t2195 * t2228;
t2179 = -t2232 * t2264 - t2233 * t2261;
t2182 = -t2232 * t2261 + t2233 * t2264;
t2053 = (t2152 * t2179 + t2155 * t2182) * t2200;
t2384 = t2053 * pkin(1);
t2237 = 1 + t2254;
t2299 = qJ(2,3) * t2333;
t2171 = t2237 * t2246 + t2299;
t2220 = qJ(2,3) * t2334;
t2174 = -t2237 * t2249 + t2220;
t2132 = t2171 * t2226 - t2174 * t2223;
t2120 = t2132 * t2196;
t2355 = t2150 * t2196;
t2313 = pkin(1) * t2355;
t2084 = t2120 - t2313;
t2382 = MDP(6) * t2084;
t2238 = 1 + t2256;
t2302 = qJ(2,2) * t2329;
t2172 = t2238 * t2247 + t2302;
t2221 = qJ(2,2) * t2330;
t2175 = -t2238 * t2250 + t2221;
t2133 = t2172 * t2227 - t2175 * t2224;
t2121 = t2133 * t2198;
t2353 = t2151 * t2198;
t2312 = pkin(1) * t2353;
t2086 = t2121 - t2312;
t2381 = MDP(6) * t2086;
t2239 = 1 + t2258;
t2305 = qJ(2,1) * t2325;
t2173 = t2239 * t2248 + t2305;
t2222 = qJ(2,1) * t2326;
t2176 = -t2239 * t2251 + t2222;
t2134 = t2173 * t2228 - t2176 * t2225;
t2122 = t2134 * t2200;
t2351 = t2152 * t2200;
t2311 = pkin(1) * t2351;
t2088 = t2122 - t2311;
t2380 = MDP(6) * t2088;
t2379 = qJ(2,1) * t2053;
t2378 = qJ(2,1) * t2261;
t2376 = qJ(2,2) * t2052;
t2375 = qJ(2,2) * t2260;
t2373 = qJ(2,3) * t2051;
t2372 = qJ(2,3) * t2259;
t2331 = t2243 * t2262;
t2183 = qJ(2,3) * t2331 - (t2254 * t2259) - t2259;
t2332 = t2243 * t2259;
t2184 = qJ(2,3) * t2332 + t2254 * t2262 + t2262;
t2144 = t2183 * t2249 + t2184 * t2246;
t2145 = t2183 * t2246 - t2184 * t2249;
t2042 = (-t2144 * t2232 + t2145 * t2233) * t2226 + (t2144 * t2233 + t2145 * t2232) * t2223;
t2370 = t2042 * t2051;
t2327 = t2244 * t2263;
t2185 = qJ(2,2) * t2327 - (t2256 * t2260) - t2260;
t2328 = t2244 * t2260;
t2186 = qJ(2,2) * t2328 + t2256 * t2263 + t2263;
t2146 = t2185 * t2250 + t2186 * t2247;
t2147 = t2185 * t2247 - t2186 * t2250;
t2043 = (-t2146 * t2232 + t2147 * t2233) * t2227 + (t2146 * t2233 + t2147 * t2232) * t2224;
t2369 = t2043 * t2052;
t2323 = t2245 * t2264;
t2187 = qJ(2,1) * t2323 - (t2258 * t2261) - t2261;
t2324 = t2245 * t2261;
t2188 = qJ(2,1) * t2324 + t2258 * t2264 + t2264;
t2148 = t2187 * t2251 + t2188 * t2248;
t2149 = t2187 * t2248 - t2188 * t2251;
t2044 = (-t2148 * t2232 + t2149 * t2233) * t2228 + (t2148 * t2233 + t2149 * t2232) * t2225;
t2368 = t2044 * t2053;
t2202 = t2331 - t2372;
t2205 = t2332 + t2371;
t2156 = t2202 * t2249 + t2205 * t2246;
t2159 = t2202 * t2246 - t2205 * t2249;
t2045 = (t2156 * t2233 + t2159 * t2232) * t2226 - (-t2156 * t2232 + t2159 * t2233) * t2223;
t2367 = t2045 * t2196;
t2203 = t2327 - t2375;
t2206 = t2328 + t2374;
t2157 = t2203 * t2250 + t2206 * t2247;
t2160 = t2203 * t2247 - t2206 * t2250;
t2046 = (t2157 * t2233 + t2160 * t2232) * t2227 - (-t2157 * t2232 + t2160 * t2233) * t2224;
t2366 = t2046 * t2198;
t2204 = t2323 - t2378;
t2207 = t2324 + t2377;
t2158 = t2204 * t2251 + t2207 * t2248;
t2161 = t2204 * t2248 - t2207 * t2251;
t2047 = (t2158 * t2233 + t2161 * t2232) * t2228 - (-t2158 * t2232 + t2161 * t2233) * t2225;
t2365 = t2047 * t2200;
t2364 = t2051 * t2243;
t2363 = t2052 * t2244;
t2362 = t2053 * t2245;
t2217 = 1 + t2280;
t2165 = t2217 * t2246 - t2299;
t2168 = t2217 * t2249 + t2220;
t2126 = t2165 * t2226 + t2168 * t2223;
t2108 = t2126 * t2196;
t2218 = 1 + t2279;
t2166 = t2218 * t2247 - t2302;
t2169 = t2218 * t2250 + t2221;
t2127 = t2166 * t2227 + t2169 * t2224;
t2109 = t2127 * t2198;
t2219 = 1 + t2278;
t2167 = t2219 * t2248 - t2305;
t2170 = t2219 * t2251 + t2222;
t2128 = t2167 * t2228 + t2170 * t2225;
t2110 = t2128 * t2200;
t2129 = -t2165 * t2223 + t2168 * t2226;
t2111 = t2129 * t2196;
t2130 = -t2166 * t2224 + t2169 * t2227;
t2112 = t2130 * t2198;
t2131 = -t2167 * t2225 + t2170 * t2228;
t2113 = t2131 * t2200;
t2135 = t2171 * t2223 + t2174 * t2226;
t2123 = t2135 * t2196;
t2136 = t2172 * t2224 + t2175 * t2227;
t2124 = t2136 * t2198;
t2137 = t2173 * t2225 + t2176 * t2228;
t2125 = t2137 * t2200;
t2197 = 1 / t2208 ^ 2;
t2361 = t2150 ^ 2 * t2197;
t2199 = 1 / t2209 ^ 2;
t2360 = t2151 ^ 2 * t2199;
t2201 = 1 / t2210 ^ 2;
t2359 = t2152 ^ 2 * t2201;
t2358 = t2153 ^ 2 * t2197;
t2357 = t2154 ^ 2 * t2199;
t2356 = t2155 ^ 2 * t2201;
t2354 = t2150 * t2197;
t2352 = t2151 * t2199;
t2350 = t2152 * t2201;
t2349 = t2153 * t2196;
t2348 = t2153 * t2197;
t2347 = t2154 * t2198;
t2346 = t2154 * t2199;
t2345 = t2155 * t2200;
t2344 = t2155 * t2201;
t2343 = t2177 * t2196;
t2342 = t2178 * t2198;
t2341 = t2179 * t2200;
t2340 = t2180 * t2196;
t2339 = t2181 * t2198;
t2338 = t2182 * t2200;
t2337 = t2196 * t2243;
t2336 = t2198 * t2244;
t2335 = t2200 * t2245;
t2322 = t2126 * t2340 + t2129 * t2343;
t2321 = t2127 * t2339 + t2130 * t2342;
t2320 = t2128 * t2338 + t2131 * t2341;
t2048 = t2132 * t2343 + t2135 * t2340;
t2049 = t2133 * t2342 + t2136 * t2339;
t2050 = t2134 * t2341 + t2137 * t2338;
t2310 = pkin(1) * t2349;
t2309 = pkin(1) * t2347;
t2308 = pkin(1) * t2345;
t2307 = qJ(2,1) * t2351;
t2306 = qJ(2,1) * t2345;
t2304 = qJ(2,2) * t2353;
t2303 = qJ(2,2) * t2347;
t2301 = qJ(2,3) * t2355;
t2300 = qJ(2,3) * t2349;
t2298 = t2051 * t2367;
t2297 = t2045 * t2354;
t2296 = t2045 * t2348;
t2295 = t2052 * t2366;
t2294 = t2046 * t2352;
t2293 = t2046 * t2346;
t2292 = t2053 * t2365;
t2291 = t2047 * t2350;
t2290 = t2047 * t2344;
t2289 = t2150 * t2348;
t2288 = t2150 * t2337;
t2287 = t2151 * t2346;
t2286 = t2151 * t2336;
t2285 = t2152 * t2344;
t2284 = t2152 * t2335;
t2283 = t2153 * t2337;
t2282 = t2154 * t2336;
t2281 = t2155 * t2335;
t2229 = t2254 + t2265;
t2271 = MDP(4) * (t2120 - 0.2e1 * t2313) + (-pkin(1) * t2132 + t2150 * t2229) * t2397 + MDP(7) * (t2111 + 0.2e1 * t2301) + MDP(8) * (-t2120 + 0.2e1 * t2288) + (qJ(2,3) * t2129 - t2132 * t2243 + t2150 * t2211) * t2394;
t2270 = MDP(4) * (t2123 - 0.2e1 * t2310) + (-pkin(1) * t2135 + t2153 * t2229) * t2397 + MDP(7) * (t2108 + 0.2e1 * t2300) + MDP(8) * (-t2123 + 0.2e1 * t2283) + (qJ(2,3) * t2126 - t2135 * t2243 + t2153 * t2211) * t2394;
t2230 = t2256 + t2265;
t2269 = MDP(4) * (t2121 - 0.2e1 * t2312) + (-pkin(1) * t2133 + t2151 * t2230) * t2396 + MDP(7) * (t2112 + 0.2e1 * t2304) + MDP(8) * (-t2121 + 0.2e1 * t2286) + (qJ(2,2) * t2130 - t2133 * t2244 + t2151 * t2212) * t2393;
t2268 = MDP(4) * (t2124 - 0.2e1 * t2309) + (-pkin(1) * t2136 + t2154 * t2230) * t2396 + MDP(7) * (t2109 + 0.2e1 * t2303) + MDP(8) * (-t2124 + 0.2e1 * t2282) + (qJ(2,2) * t2127 - t2136 * t2244 + t2154 * t2212) * t2393;
t2231 = t2258 + t2265;
t2267 = MDP(4) * (t2122 - 0.2e1 * t2311) + (-pkin(1) * t2134 + t2152 * t2231) * t2395 + MDP(7) * (t2113 + 0.2e1 * t2307) + MDP(8) * (-t2122 + 0.2e1 * t2284) + (qJ(2,1) * t2131 - t2134 * t2245 + t2152 * t2213) * t2392;
t2266 = MDP(4) * (t2125 - 0.2e1 * t2308) + (-pkin(1) * t2137 + t2155 * t2231) * t2395 + MDP(7) * (t2110 + 0.2e1 * t2306) + MDP(8) * (-t2125 + 0.2e1 * t2281) + (qJ(2,1) * t2128 - t2137 * t2245 + t2155 * t2213) * t2392;
t2189 = (t2232 ^ 2 + t2233 ^ 2) * MDP(13);
t2164 = -t2245 * t2378 + (t2314 + t2236) * t2264;
t2163 = -t2244 * t2375 + (t2315 + t2235) * t2263;
t2162 = -t2243 * t2372 + (t2316 + t2234) * t2262;
t2119 = t2164 * t2251 + t2248 * t2389;
t2118 = t2164 * t2248 - t2251 * t2389;
t2117 = t2163 * t2250 + t2247 * t2390;
t2116 = t2163 * t2247 - t2250 * t2390;
t2115 = t2162 * t2249 + t2246 * t2391;
t2114 = t2162 * t2246 - t2249 * t2391;
t2094 = t2125 - t2308;
t2092 = t2124 - t2309;
t2090 = t2123 - t2310;
t2083 = t2125 - t2281;
t2082 = t2122 - t2284;
t2081 = t2124 - t2282;
t2080 = t2121 - t2286;
t2079 = t2123 - t2283;
t2078 = t2120 - t2288;
t2065 = t2110 + t2306;
t2063 = t2113 + t2307;
t2061 = t2109 + t2303;
t2059 = t2112 + t2304;
t2057 = t2108 + t2300;
t2055 = t2111 + t2301;
t2035 = (t2118 * t2233 - t2119 * t2232) * t2228 + (t2118 * t2232 + t2119 * t2233) * t2225;
t2034 = (t2116 * t2233 - t2117 * t2232) * t2227 + (t2116 * t2232 + t2117 * t2233) * t2224;
t2033 = (t2114 * t2233 - t2115 * t2232) * t2226 + (t2114 * t2232 + t2115 * t2233) * t2223;
t1 = [(t2359 + t2360 + t2361) * MDP(1) + (qJ(2,1) * t2359 + qJ(2,2) * t2360 + qJ(2,3) * t2361) * t2387 + (t2084 * t2120 + t2086 * t2121 + t2088 * t2122) * MDP(6) + (t2055 * t2111 + t2059 * t2112 + t2063 * t2113 + t2078 * t2120 + t2080 * t2121 + t2082 * t2122) * MDP(9) + t2189 + ((t2131 * MDP(7) + t2134 * t2383) * t2201 + t2267 * t2200) * t2152 + ((t2130 * MDP(7) + t2133 * t2383) * t2199 + t2269 * t2198) * t2151 + ((t2129 * MDP(7) + t2132 * t2383) * t2197 + t2271 * t2196) * t2150; (t2285 + t2287 + t2289) * MDP(1) + (qJ(2,1) * t2285 + qJ(2,2) * t2287 + qJ(2,3) * t2289) * t2387 + (t2126 * t2354 + t2127 * t2352 + t2128 * t2350) * MDP(7) + (t2137 * t2380 + (t2063 * t2128 + t2082 * t2137) * MDP(9) + t2267 * t2155) * t2200 + (t2136 * t2381 + (t2059 * t2127 + t2080 * t2136) * MDP(9) + t2269 * t2154) * t2198 + (t2135 * t2382 + (t2055 * t2126 + t2078 * t2135) * MDP(9) + t2271 * t2153) * t2196 + t2383 * (t2135 * t2354 + t2136 * t2352 + t2137 * t2350); (t2356 + t2357 + t2358) * MDP(1) + (qJ(2,1) * t2356 + qJ(2,2) * t2357 + qJ(2,3) * t2358) * t2387 + (t2090 * t2123 + t2092 * t2124 + t2094 * t2125) * MDP(6) + (t2057 * t2108 + t2061 * t2109 + t2065 * t2110 + t2079 * t2123 + t2081 * t2124 + t2083 * t2125) * MDP(9) + t2189 + ((t2128 * MDP(7) + t2137 * t2383) * t2201 + t2266 * t2200) * t2155 + ((t2127 * MDP(7) + t2136 * t2383) * t2199 + t2268 * t2198) * t2154 + ((t2126 * MDP(7) + t2135 * t2383) * t2197 + t2270 * t2196) * t2153; (t2291 + t2294 + t2297) * MDP(1) + (qJ(2,1) * t2291 + qJ(2,2) * t2294 + qJ(2,3) * t2297) * t2387 + (t2033 * t2354 + t2034 * t2352 + t2035 * t2350) * MDP(7) - t2232 * MDP(11) - t2233 * MDP(12) + (t2044 * t2380 + (t2035 * t2063 + t2044 * t2082) * MDP(9) + t2267 * t2047) * t2200 + (t2043 * t2381 + (t2034 * t2059 + t2043 * t2080) * MDP(9) + t2269 * t2046) * t2198 + (t2042 * t2382 + (t2033 * t2055 + t2042 * t2078) * MDP(9) + t2271 * t2045) * t2196 + t2383 * (t2042 * t2354 + t2043 * t2352 + t2044 * t2350); (t2290 + t2293 + t2296) * MDP(1) + (qJ(2,1) * t2290 + qJ(2,2) * t2293 + qJ(2,3) * t2296) * t2387 + (t2033 * t2348 + t2034 * t2346 + t2035 * t2344) * MDP(7) + t2233 * MDP(11) - t2232 * MDP(12) + (t2044 * t2094 * MDP(6) + (t2035 * t2065 + t2044 * t2083) * MDP(9) + t2266 * t2047) * t2200 + (t2043 * t2092 * MDP(6) + (t2034 * t2061 + t2043 * t2081) * MDP(9) + t2268 * t2046) * t2198 + (t2042 * t2090 * MDP(6) + (t2033 * t2057 + t2042 * t2079) * MDP(9) + t2270 * t2045) * t2196 + t2383 * (t2042 * t2348 + t2043 * t2346 + t2044 * t2344); (t2292 + t2295 + t2298) * MDP(1) + ((t2047 * (t2050 - 0.2e1 * t2384) + t2368) * t2200 + (t2046 * (t2049 - 0.2e1 * t2385) + t2369) * t2198 + (t2045 * (t2048 - 0.2e1 * t2386) + t2370) * t2196) * MDP(4) + (qJ(2,1) * t2292 + qJ(2,2) * t2295 + qJ(2,3) * t2298) * t2387 + ((t2047 * (-pkin(1) * t2050 + t2053 * t2231) + t2044 * (t2050 - t2384)) * t2200 + (t2046 * (-pkin(1) * t2049 + t2052 * t2230) + t2043 * (t2049 - t2385)) * t2198 + (t2045 * (-pkin(1) * t2048 + t2051 * t2229) + t2042 * (t2048 - t2386)) * t2196) * MDP(6) + ((t2047 * (t2320 + 0.2e1 * t2379) + t2035 * t2053) * t2200 + (t2046 * (t2321 + 0.2e1 * t2376) + t2034 * t2052) * t2198 + (t2045 * (t2322 + 0.2e1 * t2373) + t2033 * t2051) * t2196) * MDP(7) + ((t2047 * (-t2050 + 0.2e1 * t2362) - t2368) * t2200 + (t2046 * (-t2049 + 0.2e1 * t2363) - t2369) * t2198 + (t2045 * (-t2048 + 0.2e1 * t2364) - t2370) * t2196) * MDP(8) + ((qJ(2,1) * t2320 - t2050 * t2245 + t2213 * t2053) * t2365 + (qJ(2,2) * t2321 - t2049 * t2244 + t2212 * t2052) * t2366 + (qJ(2,3) * t2322 - t2048 * t2243 + t2211 * t2051) * t2367 + (t2044 * (t2050 - t2362) + t2035 * (t2320 + t2379)) * t2200 + (t2043 * (t2049 - t2363) + t2034 * (t2321 + t2376)) * t2198 + (t2042 * (t2048 - t2364) + t2033 * (t2322 + t2373)) * t2196) * MDP(9) + MDP(10);];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t1(1), t1(2), t1(4); t1(2), t1(3), t1(5); t1(4), t1(5), t1(6);];
MMX  = res;
