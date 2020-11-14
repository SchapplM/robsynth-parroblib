% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RRPRR8V2G1A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [13x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRPRR8V2G1A0_convert_par2_MPV_fixb.m

% Output:
% MMX [3x3]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 21:04
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P3RRPRR8V2G1A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(13,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_mdp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:03:53
% EndTime: 2020-08-06 21:03:57
% DurationCPUTime: 3.81s
% Computational Cost: add. (4080->334), mult. (5294->611), div. (840->18), fcn. (4434->35), ass. (0->257)
t2404 = 2 * pkin(1);
t2403 = 2 * MDP(11);
t2244 = cos(qJ(2,3));
t2213 = t2244 * pkin(2);
t2195 = t2213 + pkin(1);
t2246 = cos(qJ(2,2));
t2214 = t2246 * pkin(2);
t2196 = t2214 + pkin(1);
t2248 = cos(qJ(2,1));
t2215 = t2248 * pkin(2);
t2197 = t2215 + pkin(1);
t2194 = cos(pkin(7)) * pkin(3) + pkin(2);
t2238 = sin(qJ(2,3));
t2381 = sin(pkin(7)) * pkin(3);
t2161 = t2194 * t2238 + t2244 * t2381;
t2164 = t2194 * t2244 - t2238 * t2381;
t2402 = t2161 ^ 2 / t2164 ^ 2;
t2240 = sin(qJ(2,2));
t2162 = t2194 * t2240 + t2246 * t2381;
t2165 = t2194 * t2246 - t2240 * t2381;
t2401 = t2162 ^ 2 / t2165 ^ 2;
t2242 = sin(qJ(2,1));
t2163 = t2194 * t2242 + t2248 * t2381;
t2166 = t2194 * t2248 - t2242 * t2381;
t2400 = t2163 ^ 2 / t2166 ^ 2;
t2237 = pkin(5) + qJ(3,1);
t2218 = -pkin(6) - t2237;
t2211 = 0.1e1 / t2218;
t2349 = t2211 * t2248;
t2391 = -0.2e1 * pkin(2);
t2299 = t2349 * t2391;
t2385 = pkin(1) * t2211;
t2399 = -t2299 / 0.2e1 + t2385;
t2236 = pkin(5) + qJ(3,2);
t2217 = -pkin(6) - t2236;
t2209 = 0.1e1 / t2217;
t2355 = t2209 * t2246;
t2300 = t2355 * t2391;
t2386 = pkin(1) * t2209;
t2398 = -t2300 / 0.2e1 + t2386;
t2235 = pkin(5) + qJ(3,3);
t2216 = -pkin(6) - t2235;
t2207 = 0.1e1 / t2216;
t2360 = t2207 * t2244;
t2301 = t2360 * t2391;
t2387 = pkin(1) * t2207;
t2397 = -t2301 / 0.2e1 + t2387;
t2344 = pkin(1) ^ 2 + pkin(5) ^ 2;
t2390 = 2 * pkin(5);
t2193 = (t2390 + qJ(3,1)) * qJ(3,1) + t2344;
t2255 = pkin(2) ^ 2;
t2306 = t2211 * t2248 ^ 2 * t2255;
t2396 = pkin(1) * t2299 - t2193 * t2211 - t2306;
t2192 = (t2390 + qJ(3,2)) * qJ(3,2) + t2344;
t2308 = t2209 * t2246 ^ 2 * t2255;
t2395 = pkin(1) * t2300 - t2192 * t2209 - t2308;
t2191 = (t2390 + qJ(3,3)) * qJ(3,3) + t2344;
t2310 = t2207 * t2244 ^ 2 * t2255;
t2394 = pkin(1) * t2301 - t2191 * t2207 - t2310;
t2393 = -2 * pkin(1);
t2389 = 2 * MDP(5);
t2222 = qJ(2,3) + pkin(7);
t2384 = pkin(3) * cos(t2222);
t2223 = qJ(2,2) + pkin(7);
t2383 = pkin(3) * cos(t2223);
t2224 = qJ(2,1) + pkin(7);
t2382 = pkin(3) * cos(t2224);
t2380 = 0.2e1 * pkin(2) * pkin(3);
t2188 = t2213 + t2384;
t2185 = 0.1e1 / t2188;
t2250 = 0.2e1 * qJ(2,3);
t2254 = pkin(3) ^ 2;
t2378 = (sin(t2250 + pkin(7)) * t2380 + t2255 * sin(t2250) + t2254 * sin(0.2e1 * t2222) + (sin(t2222) * pkin(3) + pkin(2) * t2238) * t2404) * t2185;
t2189 = t2214 + t2383;
t2186 = 0.1e1 / t2189;
t2251 = 0.2e1 * qJ(2,2);
t2377 = (sin(t2251 + pkin(7)) * t2380 + t2255 * sin(t2251) + t2254 * sin(0.2e1 * t2223) + (sin(t2223) * pkin(3) + pkin(2) * t2240) * t2404) * t2186;
t2252 = 0.2e1 * qJ(2,1);
t2151 = sin(t2252 + pkin(7)) * t2380 + t2255 * sin(t2252) + t2254 * sin(0.2e1 * t2224) + (sin(t2224) * pkin(3) + pkin(2) * t2242) * t2404;
t2190 = t2215 + t2382;
t2187 = 0.1e1 / t2190;
t2376 = t2151 * t2187;
t2375 = 0.1e1 / t2164 * t2161;
t2374 = t2162 / t2165;
t2373 = t2163 / t2166;
t2232 = legFrame(3,3);
t2201 = sin(t2232);
t2204 = cos(t2232);
t2239 = sin(qJ(1,3));
t2245 = cos(qJ(1,3));
t2179 = -t2201 * t2239 + t2204 * t2245;
t2372 = t2179 * t2207;
t2233 = legFrame(2,3);
t2202 = sin(t2233);
t2205 = cos(t2233);
t2241 = sin(qJ(1,2));
t2247 = cos(qJ(1,2));
t2180 = -t2202 * t2241 + t2205 * t2247;
t2371 = t2180 * t2209;
t2234 = legFrame(1,3);
t2203 = sin(t2234);
t2206 = cos(t2234);
t2243 = sin(qJ(1,1));
t2249 = cos(qJ(1,1));
t2181 = -t2203 * t2243 + t2206 * t2249;
t2370 = t2181 * t2211;
t2182 = t2201 * t2245 + t2204 * t2239;
t2369 = t2182 * t2207;
t2183 = t2202 * t2247 + t2205 * t2241;
t2368 = t2183 * t2209;
t2184 = t2203 * t2249 + t2206 * t2243;
t2367 = t2184 * t2211;
t2366 = t2185 * t2238;
t2365 = t2186 * t2240;
t2364 = t2187 * t2242;
t2208 = 0.1e1 / t2216 ^ 2;
t2225 = t2238 ^ 2;
t2359 = t2208 * t2225;
t2358 = t2208 * t2235;
t2357 = t2208 * t2238;
t2356 = t2208 * t2244;
t2210 = 0.1e1 / t2217 ^ 2;
t2226 = t2240 ^ 2;
t2354 = t2210 * t2226;
t2353 = t2210 * t2236;
t2352 = t2210 * t2240;
t2351 = t2210 * t2246;
t2350 = t2211 * t2187;
t2212 = 0.1e1 / t2218 ^ 2;
t2227 = t2242 ^ 2;
t2348 = t2212 * t2227;
t2347 = t2212 * t2237;
t2346 = t2212 * t2242;
t2345 = t2212 * t2248;
t2343 = 2 * MDP(9);
t2342 = 2 * MDP(10);
t2337 = t2208 * t2402;
t2336 = t2210 * t2401;
t2335 = t2212 * t2400;
t2334 = t2207 * t2375;
t2333 = t2208 * t2375;
t2332 = t2235 * t2375;
t2331 = t2238 * t2375;
t2330 = t2244 * t2375;
t2329 = t2209 * t2374;
t2328 = t2210 * t2374;
t2327 = t2236 * t2374;
t2326 = t2240 * t2374;
t2325 = t2246 * t2374;
t2324 = t2212 * t2373;
t2323 = t2211 * t2373;
t2322 = t2179 * t2182 * t2208;
t2321 = t2180 * t2183 * t2210;
t2320 = t2181 * t2184 * t2212;
t2319 = t2207 * t2366;
t2318 = t2185 * t2360;
t2317 = t2235 * t2366;
t2316 = t2209 * t2365;
t2315 = t2186 * t2355;
t2314 = t2236 * t2365;
t2313 = t2242 * t2350;
t2312 = t2187 * t2349;
t2311 = t2237 * t2364;
t2309 = t2238 * t2356;
t2307 = t2240 * t2351;
t2305 = t2242 * t2345;
t2304 = t2378 / 0.2e1;
t2303 = t2377 / 0.2e1;
t2302 = t2376 / 0.2e1;
t2298 = t2207 * t2331;
t2297 = t2207 * t2330;
t2296 = t2225 * t2333;
t2295 = t2208 * t2332;
t2294 = t2208 * t2331;
t2293 = t2208 * t2330;
t2292 = t2209 * t2326;
t2291 = t2209 * t2325;
t2290 = t2226 * t2328;
t2289 = t2210 * t2327;
t2288 = t2210 * t2326;
t2287 = t2210 * t2325;
t2286 = t2242 * t2323;
t2285 = t2248 * t2323;
t2284 = t2227 * t2324;
t2283 = t2237 * t2324;
t2282 = t2242 * t2324;
t2281 = t2248 * t2324;
t2280 = t2244 * t2322;
t2279 = t2246 * t2321;
t2278 = t2248 * t2320;
t2277 = t2207 * t2317;
t2276 = t2209 * t2314;
t2275 = t2211 * t2311;
t2274 = -t2207 * t2378 / 0.2e1;
t2273 = -t2209 * t2377 / 0.2e1;
t2272 = -t2151 * t2350 / 0.2e1;
t2271 = t2185 * t2298;
t2270 = t2238 * t2293;
t2269 = t2186 * t2292;
t2268 = t2240 * t2287;
t2267 = t2187 * t2286;
t2266 = t2242 * t2281;
t2261 = t2179 * t2318 + t2180 * t2315 + t2181 * t2312;
t2262 = t2179 * t2319 + t2180 * t2316 + t2181 * t2313;
t2265 = (t2179 * t2270 + t2180 * t2268 + t2181 * t2266) * t2389 + (t2179 * t2296 + t2180 * t2290 + t2181 * t2284) * MDP(4) + (t2179 * t2333 + t2180 * t2328 + t2181 * t2324) * MDP(1) - t2262 * MDP(6) - t2261 * MDP(7);
t2259 = t2182 * t2318 + t2183 * t2315 + t2184 * t2312;
t2260 = t2182 * t2319 + t2183 * t2316 + t2184 * t2313;
t2264 = (t2182 * t2270 + t2183 * t2268 + t2184 * t2266) * t2389 + (t2182 * t2296 + t2183 * t2290 + t2184 * t2284) * MDP(4) + (t2182 * t2333 + t2183 * t2328 + t2184 * t2324) * MDP(1) - t2260 * MDP(6) - t2259 * MDP(7);
t2263 = (t2238 * t2280 + t2240 * t2279 + t2242 * t2278) * t2389 + (t2235 * t2322 + t2236 * t2321 + t2237 * t2320) * t2403 + (t2225 * t2322 + t2226 * t2321 + t2227 * t2320) * MDP(4) + (t2320 + t2321 + t2322) * MDP(1) + ((-t2238 * t2322 - t2240 * t2321 - t2242 * t2320) * MDP(10) + (t2278 + t2279 + t2280) * MDP(9)) * t2404;
t2258 = t2267 + t2269 + t2271;
t2257 = t2185 * t2297 + t2186 * t2291 + t2187 * t2285;
t2178 = t2184 ^ 2;
t2177 = t2183 ^ 2;
t2176 = t2182 ^ 2;
t2175 = t2181 ^ 2;
t2174 = t2180 ^ 2;
t2173 = t2179 ^ 2;
t2172 = t2196 * t2247 - t2217 * t2241;
t2171 = t2195 * t2245 - t2216 * t2239;
t2170 = t2197 * t2249 - t2218 * t2243;
t2169 = t2197 * t2243 + t2218 * t2249;
t2168 = t2196 * t2241 + t2217 * t2247;
t2167 = t2195 * t2239 + t2216 * t2245;
t2148 = -pkin(2) * t2364 - 0.2e1 * t2237 * t2323;
t2147 = -pkin(2) * t2365 - 0.2e1 * t2209 * t2327;
t2146 = -pkin(2) * t2366 - 0.2e1 * t2207 * t2332;
t2145 = -pkin(5) * t2187 * t2248 + t2286 * t2404;
t2144 = -pkin(5) * t2186 * t2246 + t2292 * t2404;
t2143 = -pkin(5) * t2185 * t2244 + t2298 * t2404;
t2142 = -pkin(5) * t2364 + t2285 * t2393;
t2141 = -pkin(5) * t2365 + t2291 * t2393;
t2140 = -pkin(5) * t2366 + t2297 * t2393;
t2139 = t2169 * t2206 + t2170 * t2203 + t2184 * t2382;
t2138 = t2168 * t2205 + t2172 * t2202 + t2183 * t2383;
t2137 = t2167 * t2204 + t2171 * t2201 + t2182 * t2384;
t2136 = -t2169 * t2203 + t2170 * t2206 + t2181 * t2382;
t2135 = -t2168 * t2202 + t2172 * t2205 + t2180 * t2383;
t2134 = -t2167 * t2201 + t2171 * t2204 + t2179 * t2384;
t2128 = (-t2184 * t2197 + t2139) * t2211;
t2127 = (-t2183 * t2196 + t2138) * t2209;
t2126 = (-t2182 * t2195 + t2137) * t2207;
t2125 = (-t2181 * t2197 + t2136) * t2211;
t2124 = (-t2180 * t2196 + t2135) * t2209;
t2123 = (-t2179 * t2195 + t2134) * t2207;
t2117 = (-t2197 * t2373 + t2302) * t2211;
t2116 = (-t2196 * t2374 + t2303) * t2209;
t2115 = (-t2195 * t2375 + t2304) * t2207;
t2108 = t2399 * t2139 + t2396 * t2184;
t2107 = t2399 * t2136 + t2396 * t2181;
t2106 = t2398 * t2138 + t2395 * t2183;
t2105 = t2398 * t2135 + t2395 * t2180;
t2104 = t2397 * t2137 + t2394 * t2182;
t2103 = t2397 * t2134 + t2394 * t2179;
t2102 = -t2306 * t2373 + (pkin(1) * t2373 - t2376 / 0.4e1) * t2299 - pkin(2) * t2311 - t2193 * t2323 + t2302 * t2385;
t2101 = -t2308 * t2374 + (pkin(1) * t2374 - t2377 / 0.4e1) * t2300 - pkin(2) * t2314 - t2192 * t2329 + t2303 * t2386;
t2100 = -t2310 * t2375 + (pkin(1) * t2375 - t2378 / 0.4e1) * t2301 - pkin(2) * t2317 - t2191 * t2334 + t2304 * t2387;
t1 = [(t2173 * t2208 + t2174 * t2210 + t2175 * t2212) * MDP(1) + (t2173 * t2359 + t2174 * t2354 + t2175 * t2348) * MDP(4) + (-(t2107 * t2181 - t2125 * t2136) * t2211 - (t2105 * t2180 - t2124 * t2135) * t2209 - (t2103 * t2179 - t2123 * t2134) * t2207) * MDP(12) + MDP(13) + (t2173 * t2309 + t2174 * t2307 + t2175 * t2305) * t2389 + (t2173 * t2358 + t2174 * t2353 + t2175 * t2347) * t2403 + ((t2173 * t2356 + t2174 * t2351 + t2175 * t2345) * MDP(9) + (-t2173 * t2357 - t2174 * t2352 - t2175 * t2346) * MDP(10)) * t2404; (-(t2108 * t2181 - t2128 * t2136) * t2211 - (t2106 * t2180 - t2127 * t2135) * t2209 - (t2104 * t2179 - t2126 * t2134) * t2207) * MDP(12) + t2263; (-t2140 * t2372 - t2141 * t2371 - t2142 * t2370) * MDP(9) + (-t2143 * t2372 - t2144 * t2371 - t2145 * t2370) * MDP(10) + (-t2146 * t2372 - t2147 * t2371 - t2148 * t2370) * MDP(11) + (-(t2102 * t2181 - t2117 * t2136) * t2211 - (t2101 * t2180 - t2116 * t2135) * t2209 - (t2100 * t2179 - t2115 * t2134) * t2207) * MDP(12) + t2265; (-(t2107 * t2184 - t2125 * t2139) * t2211 - (t2105 * t2183 - t2124 * t2138) * t2209 - (t2103 * t2182 - t2123 * t2137) * t2207) * MDP(12) + t2263; (t2176 * t2208 + t2177 * t2210 + t2178 * t2212) * MDP(1) + (t2176 * t2359 + t2177 * t2354 + t2178 * t2348) * MDP(4) + (-(t2108 * t2184 - t2128 * t2139) * t2211 - (t2106 * t2183 - t2127 * t2138) * t2209 - (t2104 * t2182 - t2126 * t2137) * t2207) * MDP(12) + MDP(13) + (t2176 * t2309 + t2177 * t2307 + t2178 * t2305) * t2389 + (t2176 * t2358 + t2177 * t2353 + t2178 * t2347) * t2403 + ((t2176 * t2356 + t2177 * t2351 + t2178 * t2345) * MDP(9) + (-t2176 * t2357 - t2177 * t2352 - t2178 * t2346) * MDP(10)) * t2404; (-t2140 * t2369 - t2141 * t2368 - t2142 * t2367) * MDP(9) + (-t2143 * t2369 - t2144 * t2368 - t2145 * t2367) * MDP(10) + (-t2146 * t2369 - t2147 * t2368 - t2148 * t2367) * MDP(11) + (-(t2102 * t2184 - t2117 * t2139) * t2211 - (t2101 * t2183 - t2116 * t2138) * t2209 - (t2100 * t2182 - t2115 * t2137) * t2207) * MDP(12) + t2264; (t2179 * t2295 + t2180 * t2289 + t2181 * t2283) * t2403 + (-t2103 * t2334 - t2105 * t2329 - t2107 * t2323 - t2123 * t2274 - t2124 * t2273 - t2125 * t2272) * MDP(12) + (t2262 * MDP(11) + (t2179 * t2277 + t2180 * t2276 + t2181 * t2275) * MDP(12)) * pkin(2) + (MDP(10) * t2261 + MDP(9) * t2262) * pkin(5) + ((t2179 * t2293 + t2180 * t2287 + t2181 * t2281) * t2343 + (-t2179 * t2294 - t2180 * t2288 - t2181 * t2282) * t2342) * pkin(1) + t2265; (t2182 * t2295 + t2183 * t2289 + t2184 * t2283) * t2403 + (-t2104 * t2334 - t2106 * t2329 - t2108 * t2323 - t2126 * t2274 - t2127 * t2273 - t2128 * t2272) * MDP(12) + (t2260 * MDP(11) + (t2182 * t2277 + t2183 * t2276 + t2184 * t2275) * MDP(12)) * pkin(2) + (MDP(10) * t2259 + MDP(9) * t2260) * pkin(5) + ((t2182 * t2293 + t2183 * t2287 + t2184 * t2281) * t2343 + (-t2182 * t2294 - t2183 * t2288 - t2184 * t2282) * t2342) * pkin(1) + t2264; (t2335 + t2336 + t2337) * MDP(1) + (t2225 * t2337 + t2226 * t2336 + t2227 * t2335) * MDP(4) + (t2305 * t2400 + t2307 * t2401 + t2309 * t2402) * t2389 - 0.2e1 * t2258 * MDP(6) - 0.2e1 * t2257 * MDP(7) + (0.1e1 / t2190 ^ 2 + 0.1e1 / t2189 ^ 2 + 0.1e1 / t2188 ^ 2) * MDP(8) + (pkin(5) * t2258 - t2140 * t2334 - t2141 * t2329 - t2142 * t2323) * MDP(9) + (pkin(5) * t2257 - t2143 * t2334 - t2144 * t2329 - t2145 * t2323) * MDP(10) + (pkin(2) * t2258 - t2146 * t2334 - t2147 * t2329 - t2148 * t2323) * MDP(11) + (-(t2102 * t2373 - t2117 * t2302) * t2211 - (t2101 * t2374 - t2116 * t2303) * t2209 - (t2100 * t2375 - t2115 * t2304) * t2207 + (t2235 * t2271 + t2236 * t2269 + t2237 * t2267 + (t2185 ^ 2 + t2186 ^ 2 + t2187 ^ 2) * pkin(2)) * pkin(2)) * MDP(12) + MDP(13);];
%% Postprocessing: Reshape Output
% From vec2mat_3_matlab.m
res = [t1(1), t1(2), t1(3); t1(4), t1(5), t1(6); t1(7), t1(8), t1(9);];
MMX  = res;
