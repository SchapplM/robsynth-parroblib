% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3PRRRR8V2G3A0
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
%   see P3PRRRR8V2G3A0_convert_par2_MPV_fixb.m

% Output:
% MMX [3x3]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:05
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P3PRRRR8V2G3A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(8,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G3A0_inertia_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G3A0_inertia_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G3A0_inertia_para_pf_mdp: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G3A0_inertia_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G3A0_inertia_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3PRRRR8V2G3A0_inertia_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:05:35
% EndTime: 2020-08-06 18:05:43
% DurationCPUTime: 8.30s
% Computational Cost: add. (13434->490), mult. (31299->1095), div. (1272->8), fcn. (31404->22), ass. (0->395)
t2207 = sin(pkin(8));
t2209 = cos(pkin(8));
t2225 = cos(qJ(2,1));
t2210 = cos(pkin(4));
t2219 = sin(qJ(2,1));
t2371 = t2210 * t2219;
t2179 = t2207 * t2371 - t2209 * t2225;
t2218 = sin(qJ(3,1));
t2208 = sin(pkin(4));
t2224 = cos(qJ(3,1));
t2380 = t2208 * t2224;
t2168 = t2179 * t2218 + t2207 * t2380;
t2182 = t2207 * t2225 + t2209 * t2371;
t2173 = t2182 * t2218 + t2209 * t2380;
t2466 = t2168 * t2173;
t2221 = cos(qJ(2,3));
t2215 = sin(qJ(2,3));
t2375 = t2210 * t2215;
t2177 = t2207 * t2375 - t2209 * t2221;
t2214 = sin(qJ(3,3));
t2220 = cos(qJ(3,3));
t2384 = t2208 * t2220;
t2169 = t2177 * t2214 + t2207 * t2384;
t2180 = t2207 * t2221 + t2209 * t2375;
t2171 = t2180 * t2214 + t2209 * t2384;
t2465 = t2169 * t2171;
t2223 = cos(qJ(2,2));
t2217 = sin(qJ(2,2));
t2373 = t2210 * t2217;
t2178 = t2207 * t2373 - t2209 * t2223;
t2216 = sin(qJ(3,2));
t2222 = cos(qJ(3,2));
t2382 = t2208 * t2222;
t2170 = t2178 * t2216 + t2207 * t2382;
t2181 = t2207 * t2223 + t2209 * t2373;
t2172 = t2181 * t2216 + t2209 * t2382;
t2464 = t2170 * t2172;
t2226 = pkin(7) + pkin(6);
t2186 = pkin(2) * t2221 + t2215 * t2226;
t2183 = pkin(2) * t2215 - t2221 * t2226;
t2390 = t2208 * t2214;
t2267 = pkin(3) * t2390 - t2183 * t2210;
t2159 = t2186 * t2209 + t2267 * t2207;
t2376 = t2210 * t2214;
t2174 = pkin(3) * t2376 + t2183 * t2208;
t2211 = legFrame(3,2);
t2195 = sin(t2211);
t2198 = cos(t2211);
t2389 = t2208 * t2215;
t2392 = t2207 * t2208;
t2449 = pkin(3) * t2220 ^ 2;
t2452 = pkin(2) * t2214;
t2132 = -(t2177 * t2198 - t2195 * t2389) * t2449 + (t2159 * t2198 + t2174 * t2195) * t2220 + (t2195 * t2210 + t2198 * t2392) * t2452;
t2135 = (t2177 * t2195 + t2198 * t2389) * t2449 + (-t2159 * t2195 + t2174 * t2198) * t2220 + (-t2195 * t2392 + t2198 * t2210) * t2452;
t2463 = t2171 * (t2132 * t2195 - t2135 * t2198);
t2187 = pkin(2) * t2223 + t2217 * t2226;
t2184 = pkin(2) * t2217 - t2223 * t2226;
t2388 = t2208 * t2216;
t2266 = pkin(3) * t2388 - t2184 * t2210;
t2160 = t2187 * t2209 + t2266 * t2207;
t2374 = t2210 * t2216;
t2175 = pkin(3) * t2374 + t2184 * t2208;
t2212 = legFrame(2,2);
t2196 = sin(t2212);
t2199 = cos(t2212);
t2387 = t2208 * t2217;
t2448 = pkin(3) * t2222 ^ 2;
t2451 = pkin(2) * t2216;
t2133 = -(t2178 * t2199 - t2196 * t2387) * t2448 + (t2160 * t2199 + t2175 * t2196) * t2222 + (t2196 * t2210 + t2199 * t2392) * t2451;
t2136 = (t2178 * t2196 + t2199 * t2387) * t2448 + (-t2160 * t2196 + t2175 * t2199) * t2222 + (-t2196 * t2392 + t2199 * t2210) * t2451;
t2462 = t2172 * (t2133 * t2196 - t2136 * t2199);
t2188 = pkin(2) * t2225 + t2219 * t2226;
t2185 = pkin(2) * t2219 - t2225 * t2226;
t2386 = t2208 * t2218;
t2265 = pkin(3) * t2386 - t2185 * t2210;
t2161 = t2188 * t2209 + t2265 * t2207;
t2372 = t2210 * t2218;
t2176 = pkin(3) * t2372 + t2185 * t2208;
t2213 = legFrame(1,2);
t2197 = sin(t2213);
t2200 = cos(t2213);
t2385 = t2208 * t2219;
t2447 = pkin(3) * t2224 ^ 2;
t2450 = pkin(2) * t2218;
t2134 = -(t2179 * t2200 - t2197 * t2385) * t2447 + (t2161 * t2200 + t2176 * t2197) * t2224 + (t2197 * t2210 + t2200 * t2392) * t2450;
t2137 = (t2179 * t2197 + t2200 * t2385) * t2447 + (-t2161 * t2197 + t2176 * t2200) * t2224 + (-t2197 * t2392 + t2200 * t2210) * t2450;
t2461 = t2173 * (t2134 * t2197 - t2137 * t2200);
t2140 = -t2182 * t2447 - t2188 * t2207 * t2224 + (pkin(2) * t2386 + t2265 * t2224) * t2209;
t2158 = pkin(2) * t2372 + t2176 * t2224 + t2385 * t2447;
t2154 = 0.1e1 / t2158;
t2122 = (pkin(6) * t2168 + t2140 * t2385) * t2154;
t2327 = t2154 * t2210 * t2224;
t2078 = -t2122 * t2218 + t2140 * t2327;
t2328 = t2154 * t2372;
t2081 = -t2122 * t2224 - t2140 * t2328;
t2368 = t2210 * t2225;
t2391 = t2207 * t2210;
t2444 = pkin(3) * t2224;
t2149 = (t2207 * t2368 + t2209 * t2219) * t2444 + t2188 * t2391 + t2185 * t2209;
t2379 = t2208 * t2225;
t2457 = 0.2e1 * pkin(2);
t2268 = t2140 * t2379 + t2168 * t2457;
t2227 = 0.1e1 / pkin(3);
t2443 = pkin(6) * t2227;
t2361 = t2224 * t2443;
t2091 = (-t2149 * t2361 - t2268 * t2218) * t2154;
t2364 = t2218 * t2443;
t2092 = (-t2149 * t2364 + t2268 * t2224) * t2154;
t2377 = t2209 * t2210;
t2146 = (t2207 * t2219 - t2209 * t2368) * t2444 - t2188 * t2377 + t2185 * t2207;
t2426 = t2146 * t2227;
t2460 = MDP(10) * (t2078 * t2426 - t2092 * t2173) + MDP(11) * (t2081 * t2426 - t2091 * t2173);
t2139 = -t2181 * t2448 - t2187 * t2207 * t2222 + (pkin(2) * t2388 + t2266 * t2222) * t2209;
t2157 = pkin(2) * t2374 + t2175 * t2222 + t2387 * t2448;
t2152 = 0.1e1 / t2157;
t2121 = (pkin(6) * t2170 + t2139 * t2387) * t2152;
t2334 = t2152 * t2210 * t2222;
t2077 = -t2121 * t2216 + t2139 * t2334;
t2335 = t2152 * t2374;
t2080 = -t2121 * t2222 - t2139 * t2335;
t2369 = t2210 * t2223;
t2445 = pkin(3) * t2222;
t2148 = (t2207 * t2369 + t2209 * t2217) * t2445 + t2187 * t2391 + t2184 * t2209;
t2381 = t2208 * t2223;
t2269 = t2139 * t2381 + t2170 * t2457;
t2362 = t2222 * t2443;
t2089 = (-t2148 * t2362 - t2269 * t2216) * t2152;
t2365 = t2216 * t2443;
t2090 = (-t2148 * t2365 + t2269 * t2222) * t2152;
t2145 = (t2207 * t2217 - t2209 * t2369) * t2445 - t2187 * t2377 + t2184 * t2207;
t2429 = t2145 * t2227;
t2459 = MDP(10) * (t2077 * t2429 - t2090 * t2172) + MDP(11) * (t2080 * t2429 - t2089 * t2172);
t2138 = -t2180 * t2449 - t2186 * t2207 * t2220 + (pkin(2) * t2390 + t2267 * t2220) * t2209;
t2156 = pkin(2) * t2376 + t2174 * t2220 + t2389 * t2449;
t2150 = 0.1e1 / t2156;
t2120 = (pkin(6) * t2169 + t2138 * t2389) * t2150;
t2341 = t2150 * t2210 * t2220;
t2076 = -t2120 * t2214 + t2138 * t2341;
t2342 = t2150 * t2376;
t2079 = -t2120 * t2220 - t2138 * t2342;
t2370 = t2210 * t2221;
t2446 = pkin(3) * t2220;
t2147 = (t2207 * t2370 + t2209 * t2215) * t2446 + t2186 * t2391 + t2183 * t2209;
t2383 = t2208 * t2221;
t2270 = t2138 * t2383 + t2169 * t2457;
t2363 = t2220 * t2443;
t2087 = (-t2147 * t2363 - t2270 * t2214) * t2150;
t2366 = t2214 * t2443;
t2088 = (-t2147 * t2366 + t2270 * t2220) * t2150;
t2144 = (t2207 * t2215 - t2209 * t2370) * t2446 - t2186 * t2377 + t2183 * t2207;
t2432 = t2144 * t2227;
t2458 = MDP(10) * (t2076 * t2432 - t2088 * t2171) + MDP(11) * (t2079 * t2432 - t2087 * t2171);
t2456 = 2 * MDP(6);
t2455 = 2 * MDP(7);
t2454 = 2 * MDP(8);
t2453 = 0.2e1 * t2208;
t2442 = MDP(3) * t2208;
t2441 = MDP(4) * t2208;
t2440 = MDP(7) * t2227;
t2439 = MDP(8) * t2227;
t2438 = MDP(9) / pkin(3) ^ 2;
t2151 = 0.1e1 / t2156 ^ 2;
t2437 = t2138 * t2151;
t2153 = 0.1e1 / t2157 ^ 2;
t2436 = t2139 * t2153;
t2155 = 0.1e1 / t2158 ^ 2;
t2435 = t2140 * t2155;
t2434 = t2144 * t2150;
t2433 = t2144 * t2151;
t2431 = t2145 * t2152;
t2430 = t2145 * t2153;
t2428 = t2146 * t2154;
t2427 = t2146 * t2155;
t2425 = t2147 * t2150;
t2424 = t2147 * t2227;
t2423 = t2148 * t2152;
t2422 = t2148 * t2227;
t2421 = t2149 * t2154;
t2420 = t2149 * t2227;
t2419 = t2151 * t2169;
t2189 = t2195 ^ 2;
t2418 = t2151 * t2189;
t2192 = t2198 ^ 2;
t2417 = t2151 * t2192;
t2416 = t2151 * t2198;
t2201 = t2214 ^ 2;
t2415 = t2151 * t2201;
t2414 = t2151 * t2215;
t2413 = t2151 * t2221;
t2412 = t2153 * t2170;
t2190 = t2196 ^ 2;
t2411 = t2153 * t2190;
t2193 = t2199 ^ 2;
t2410 = t2153 * t2193;
t2409 = t2153 * t2199;
t2202 = t2216 ^ 2;
t2408 = t2153 * t2202;
t2407 = t2153 * t2217;
t2406 = t2153 * t2223;
t2405 = t2155 * t2168;
t2191 = t2197 ^ 2;
t2404 = t2155 * t2191;
t2194 = t2200 ^ 2;
t2403 = t2155 * t2194;
t2402 = t2155 * t2200;
t2203 = t2218 ^ 2;
t2401 = t2155 * t2203;
t2400 = t2155 * t2219;
t2399 = t2155 * t2225;
t2398 = t2171 * t2195;
t2397 = t2171 * t2198;
t2396 = t2172 * t2196;
t2395 = t2172 * t2199;
t2394 = t2173 * t2197;
t2393 = t2173 * t2200;
t2378 = t2208 * t2227;
t2367 = t2210 * t2227;
t2360 = t2138 * t2419;
t2359 = t2139 * t2412;
t2358 = t2140 * t2405;
t2357 = t2147 * t2433;
t2356 = t2195 * t2434;
t2355 = t2198 * t2434;
t2354 = t2171 * t2433;
t2353 = t2148 * t2430;
t2352 = t2196 * t2431;
t2351 = t2199 * t2431;
t2350 = t2172 * t2430;
t2349 = t2149 * t2427;
t2348 = t2197 * t2428;
t2347 = t2200 * t2428;
t2346 = t2173 * t2427;
t2345 = t2147 * t2419;
t2344 = t2148 * t2412;
t2343 = t2149 * t2405;
t2165 = t2171 ^ 2;
t2340 = t2165 * t2415;
t2339 = t2151 * t2398;
t2338 = t2151 * t2397;
t2337 = t2195 * t2416;
t2336 = t2151 * t2214 * t2220;
t2166 = t2172 ^ 2;
t2333 = t2166 * t2408;
t2332 = t2153 * t2396;
t2331 = t2153 * t2395;
t2330 = t2196 * t2409;
t2329 = t2153 * t2216 * t2222;
t2167 = t2173 ^ 2;
t2326 = t2167 * t2401;
t2325 = t2155 * t2394;
t2324 = t2155 * t2393;
t2323 = t2197 * t2402;
t2322 = t2155 * t2218 * t2224;
t2321 = t2214 * t2383;
t2320 = t2216 * t2381;
t2319 = t2218 * t2379;
t2318 = t2220 * t2383;
t2317 = t2222 * t2381;
t2316 = t2224 * t2379;
t2315 = t2132 * t2338;
t2314 = t2133 * t2331;
t2313 = t2134 * t2324;
t2312 = t2135 * t2339;
t2311 = t2136 * t2332;
t2310 = t2137 * t2325;
t2309 = t2214 * t2354;
t2308 = t2220 * t2354;
t2307 = t2216 * t2350;
t2306 = t2222 * t2350;
t2305 = t2218 * t2346;
t2304 = t2224 * t2346;
t2303 = t2150 * t2215 * t2378;
t2302 = t2165 * t2337;
t2301 = t2165 * t2336;
t2300 = t2415 * t2465;
t2299 = t2152 * t2217 * t2378;
t2298 = t2166 * t2330;
t2297 = t2166 * t2329;
t2296 = t2408 * t2464;
t2295 = t2154 * t2219 * t2378;
t2294 = t2167 * t2323;
t2293 = t2167 * t2322;
t2292 = t2401 * t2466;
t2291 = t2144 * t2171 * t2337;
t2290 = t2145 * t2172 * t2330;
t2289 = t2146 * t2173 * t2323;
t2288 = t2214 * t2303;
t2287 = t2220 * t2303;
t2286 = t2336 * t2465;
t2285 = t2216 * t2299;
t2284 = t2222 * t2299;
t2283 = t2329 * t2464;
t2282 = t2218 * t2295;
t2281 = t2224 * t2295;
t2280 = t2322 * t2466;
t2129 = (t2147 * t2367 + t2169 * t2383) * t2150;
t2105 = t2129 * t2220 - t2147 * t2288;
t2106 = -t2129 * t2214 - t2147 * t2287;
t2279 = t2105 * MDP(10) + t2106 * MDP(11);
t2130 = (t2148 * t2367 + t2170 * t2381) * t2152;
t2107 = t2130 * t2222 - t2148 * t2285;
t2108 = -t2130 * t2216 - t2148 * t2284;
t2278 = t2107 * MDP(10) + t2108 * MDP(11);
t2131 = (t2149 * t2367 + t2168 * t2379) * t2154;
t2109 = t2131 * t2224 - t2149 * t2282;
t2110 = -t2131 * t2218 - t2149 * t2281;
t2277 = t2109 * MDP(10) + t2110 * MDP(11);
t2273 = t2144 * t2169 - t2147 * t2171;
t2272 = t2145 * t2170 - t2148 * t2172;
t2271 = t2146 * t2168 - t2149 * t2173;
t2264 = t2144 * t2288;
t2263 = t2144 * t2287;
t2262 = t2145 * t2285;
t2261 = t2145 * t2284;
t2260 = t2146 * t2282;
t2259 = t2146 * t2281;
t2252 = -t2132 * t2169 + t2138 * t2397;
t2251 = -t2133 * t2170 + t2139 * t2395;
t2250 = -t2134 * t2168 + t2140 * t2393;
t2249 = t2135 * t2169 + t2138 * t2398;
t2248 = t2136 * t2170 + t2139 * t2396;
t2247 = t2137 * t2168 + t2140 * t2394;
t2246 = t2144 * t2363 - 0.2e1 * t2171 * t2452;
t2245 = t2171 * t2220 * t2457 + t2144 * t2366;
t2244 = t2145 * t2362 - 0.2e1 * t2172 * t2451;
t2243 = t2172 * t2222 * t2457 + t2145 * t2365;
t2242 = t2146 * t2361 - 0.2e1 * t2173 * t2450;
t2241 = t2173 * t2224 * t2457 + t2146 * t2364;
t2240 = t2151 * t2195 * t2273;
t2239 = t2273 * t2416;
t2238 = t2153 * t2196 * t2272;
t2237 = t2272 * t2409;
t2236 = t2155 * t2197 * t2271;
t2235 = t2271 * t2402;
t2234 = t2150 * (-t2144 * t2367 + t2171 * t2383);
t2233 = t2152 * (-t2145 * t2367 + t2172 * t2381);
t2232 = t2154 * (-t2146 * t2367 + t2173 * t2379);
t2141 = t2144 ^ 2;
t2142 = t2145 ^ 2;
t2143 = t2146 ^ 2;
t2231 = (-t2400 * t2461 - t2407 * t2462 - t2414 * t2463) * t2441 + (t2399 * t2461 + t2406 * t2462 + t2413 * t2463) * t2442 + (t2132 * t2135 * t2151 + t2133 * t2136 * t2153 + t2134 * t2137 * t2155) * MDP(1) + 0.2e1 * (t2214 * t2291 + t2216 * t2290 + t2218 * t2289) * t2440 + 0.2e1 * (t2220 * t2291 + t2222 * t2290 + t2224 * t2289) * t2439 + (-t2141 * t2337 - t2142 * t2330 - t2143 * t2323) * t2438 + (-t2195 * t2198 * t2301 - t2196 * t2199 * t2297 - t2197 * t2200 * t2293) * t2456 + (-t2201 * t2302 - t2202 * t2298 - t2203 * t2294) * MDP(5) + (-t2294 - t2298 - t2302) * MDP(2);
t2230 = (t2250 * t2400 + t2251 * t2407 + t2252 * t2414) * t2441 + (-t2250 * t2399 - t2251 * t2406 - t2252 * t2413) * t2442 + (t2214 * t2239 + t2216 * t2237 + t2218 * t2235) * t2440 + (t2220 * t2239 + t2222 * t2237 + t2224 * t2235) * t2439 + (t2132 * t2437 + t2133 * t2436 + t2134 * t2435) * MDP(1) + (t2198 * t2357 + t2199 * t2353 + t2200 * t2349) * t2438 + (-t2198 * t2286 - t2199 * t2283 - t2200 * t2280) * t2456 + (-t2198 * t2300 - t2199 * t2296 - t2200 * t2292) * MDP(5) + (-t2168 * t2324 - t2169 * t2338 - t2170 * t2331) * MDP(2);
t2229 = (-t2247 * t2400 - t2248 * t2407 - t2249 * t2414) * t2441 + (t2247 * t2399 + t2248 * t2406 + t2249 * t2413) * t2442 + (-t2214 * t2240 - t2216 * t2238 - t2218 * t2236) * t2440 + (-t2220 * t2240 - t2222 * t2238 - t2224 * t2236) * t2439 + (t2135 * t2437 + t2136 * t2436 + t2137 * t2435) * MDP(1) + (-t2195 * t2357 - t2196 * t2353 - t2197 * t2349) * t2438 + (t2195 * t2286 + t2196 * t2283 + t2197 * t2280) * t2456 + (t2195 * t2300 + t2196 * t2296 + t2197 * t2292) * MDP(5) + (t2168 * t2325 + t2169 * t2339 + t2170 * t2332) * MDP(2);
t2164 = t2170 ^ 2;
t2163 = t2169 ^ 2;
t2162 = t2168 ^ 2;
t2128 = t2200 * t2232;
t2127 = t2197 * t2232;
t2126 = t2199 * t2233;
t2125 = t2196 * t2233;
t2124 = t2198 * t2234;
t2123 = t2195 * t2234;
t2116 = (pkin(6) * t2394 + t2137 * t2385) * t2154;
t2115 = (-pkin(6) * t2393 + t2134 * t2385) * t2154;
t2114 = (pkin(6) * t2396 + t2136 * t2387) * t2152;
t2113 = (-pkin(6) * t2395 + t2133 * t2387) * t2152;
t2112 = (pkin(6) * t2398 + t2135 * t2389) * t2150;
t2111 = (-pkin(6) * t2397 + t2132 * t2389) * t2150;
t2104 = -t2128 * t2224 - t2200 * t2260;
t2103 = t2128 * t2218 - t2200 * t2259;
t2102 = -t2127 * t2218 + t2197 * t2259;
t2101 = t2127 * t2224 + t2197 * t2260;
t2100 = -t2126 * t2222 - t2199 * t2262;
t2099 = t2126 * t2216 - t2199 * t2261;
t2098 = -t2125 * t2216 + t2196 * t2261;
t2097 = t2125 * t2222 + t2196 * t2262;
t2096 = -t2124 * t2220 - t2198 * t2264;
t2095 = t2124 * t2214 - t2198 * t2263;
t2094 = -t2123 * t2214 + t2195 * t2263;
t2093 = t2123 * t2220 + t2195 * t2264;
t2073 = (t2137 * t2316 + t2241 * t2197) * t2154;
t2072 = (t2134 * t2316 - t2241 * t2200) * t2154;
t2071 = (-t2137 * t2319 + t2242 * t2197) * t2154;
t2070 = (-t2134 * t2319 - t2242 * t2200) * t2154;
t2069 = (t2136 * t2317 + t2243 * t2196) * t2152;
t2068 = (t2133 * t2317 - t2243 * t2199) * t2152;
t2067 = (-t2136 * t2320 + t2244 * t2196) * t2152;
t2066 = (-t2133 * t2320 - t2244 * t2199) * t2152;
t2065 = (t2135 * t2318 + t2245 * t2195) * t2150;
t2064 = (t2132 * t2318 - t2245 * t2198) * t2150;
t2063 = (-t2135 * t2321 + t2246 * t2195) * t2150;
t2062 = (-t2132 * t2321 - t2246 * t2198) * t2150;
t2059 = -t2116 * t2224 - t2137 * t2328;
t2058 = -t2115 * t2224 - t2134 * t2328;
t2057 = -t2114 * t2222 - t2136 * t2335;
t2056 = -t2113 * t2222 - t2133 * t2335;
t2055 = -t2112 * t2220 - t2135 * t2342;
t2054 = -t2111 * t2220 - t2132 * t2342;
t2053 = -t2116 * t2218 + t2137 * t2327;
t2052 = -t2115 * t2218 + t2134 * t2327;
t2051 = -t2114 * t2216 + t2136 * t2334;
t2050 = -t2113 * t2216 + t2133 * t2334;
t2049 = -t2112 * t2214 + t2135 * t2341;
t2048 = -t2111 * t2214 + t2132 * t2341;
t1 = [(t2132 ^ 2 * t2151 + t2133 ^ 2 * t2153 + t2134 ^ 2 * t2155) * MDP(1) + (t2165 * t2417 + t2166 * t2410 + t2167 * t2403) * MDP(2) + (t2192 * t2340 + t2193 * t2333 + t2194 * t2326) * MDP(5) + (t2192 * t2301 + t2193 * t2297 + t2194 * t2293) * t2456 + MDP(12) + ((-t2072 * t2393 + t2104 * t2134) * MDP(10) + (-t2070 * t2393 + t2103 * t2134) * MDP(11)) * t2154 + ((-t2068 * t2395 + t2100 * t2133) * MDP(10) + (-t2066 * t2395 + t2099 * t2133) * MDP(11)) * t2152 + ((-t2064 * t2397 + t2096 * t2132) * MDP(10) + (-t2062 * t2397 + t2095 * t2132) * MDP(11)) * t2150 + (t2141 * t2417 + t2142 * t2410 + t2143 * t2403) * t2438 + ((-t2221 * t2315 - t2223 * t2314 - t2225 * t2313) * MDP(3) + (t2215 * t2315 + t2217 * t2314 + t2219 * t2313) * MDP(4)) * t2453 + ((t2048 * t2355 + t2050 * t2351 + t2052 * t2347) * MDP(10) + (t2054 * t2355 + t2056 * t2351 + t2058 * t2347) * MDP(11) + (-t2192 * t2309 - t2193 * t2307 - t2194 * t2305) * t2455 + (-t2192 * t2308 - t2193 * t2306 - t2194 * t2304) * t2454) * t2227; ((t2101 * MDP(10) + t2102 * MDP(11)) * t2134 + ((t2053 * t2426 - t2073 * t2173) * MDP(10) + (t2059 * t2426 - t2071 * t2173) * MDP(11)) * t2200) * t2154 + ((t2097 * MDP(10) + t2098 * MDP(11)) * t2133 + ((t2051 * t2429 - t2069 * t2172) * MDP(10) + (t2057 * t2429 - t2067 * t2172) * MDP(11)) * t2199) * t2152 + ((t2093 * MDP(10) + t2094 * MDP(11)) * t2132 + ((t2049 * t2432 - t2065 * t2171) * MDP(10) + (t2055 * t2432 - t2063 * t2171) * MDP(11)) * t2198) * t2150 + t2231; (t2277 * t2134 + t2460 * t2200) * t2154 + (t2278 * t2133 + t2459 * t2199) * t2152 + (t2279 * t2132 + t2458 * t2198) * t2150 + t2230; ((MDP(10) * t2104 + MDP(11) * t2103) * t2137 + ((-t2052 * t2426 + t2072 * t2173) * MDP(10) + (-t2058 * t2426 + t2070 * t2173) * MDP(11)) * t2197) * t2154 + ((MDP(10) * t2100 + MDP(11) * t2099) * t2136 + ((-t2050 * t2429 + t2068 * t2172) * MDP(10) + (-t2056 * t2429 + t2066 * t2172) * MDP(11)) * t2196) * t2152 + ((MDP(10) * t2096 + MDP(11) * t2095) * t2135 + ((-t2048 * t2432 + t2064 * t2171) * MDP(10) + (-t2054 * t2432 + t2062 * t2171) * MDP(11)) * t2195) * t2150 + t2231; (t2135 ^ 2 * t2151 + t2136 ^ 2 * t2153 + t2137 ^ 2 * t2155) * MDP(1) + (t2165 * t2418 + t2166 * t2411 + t2167 * t2404) * MDP(2) + (t2189 * t2340 + t2190 * t2333 + t2191 * t2326) * MDP(5) + (t2189 * t2301 + t2190 * t2297 + t2191 * t2293) * t2456 + MDP(12) + ((t2073 * t2394 + t2101 * t2137) * MDP(10) + (t2071 * t2394 + t2102 * t2137) * MDP(11)) * t2154 + ((t2069 * t2396 + t2097 * t2136) * MDP(10) + (t2067 * t2396 + t2098 * t2136) * MDP(11)) * t2152 + ((t2065 * t2398 + t2093 * t2135) * MDP(10) + (t2063 * t2398 + t2094 * t2135) * MDP(11)) * t2150 + (t2141 * t2418 + t2142 * t2411 + t2143 * t2404) * t2438 + ((t2221 * t2312 + t2223 * t2311 + t2225 * t2310) * MDP(3) + (-t2215 * t2312 - t2217 * t2311 - t2219 * t2310) * MDP(4)) * t2453 + ((-t2049 * t2356 - t2051 * t2352 - t2053 * t2348) * MDP(10) + (-t2055 * t2356 - t2057 * t2352 - t2059 * t2348) * MDP(11) + (-t2189 * t2309 - t2190 * t2307 - t2191 * t2305) * t2455 + (-t2189 * t2308 - t2190 * t2306 - t2191 * t2304) * t2454) * t2227; (t2277 * t2137 - t2460 * t2197) * t2154 + (t2278 * t2136 - t2459 * t2196) * t2152 + (t2279 * t2135 - t2458 * t2195) * t2150 + t2229; ((t2052 * t2420 + t2072 * t2168 + t2104 * t2140) * MDP(10) + (t2058 * t2420 + t2070 * t2168 + t2103 * t2140) * MDP(11)) * t2154 + ((t2050 * t2422 + t2068 * t2170 + t2100 * t2139) * MDP(10) + (t2056 * t2422 + t2066 * t2170 + t2099 * t2139) * MDP(11)) * t2152 + ((t2048 * t2424 + t2064 * t2169 + t2096 * t2138) * MDP(10) + (t2054 * t2424 + t2062 * t2169 + t2095 * t2138) * MDP(11)) * t2150 + t2230; ((t2053 * t2420 + t2073 * t2168 + t2101 * t2140) * MDP(10) + (t2059 * t2420 + t2071 * t2168 + t2102 * t2140) * MDP(11)) * t2154 + ((t2051 * t2422 + t2069 * t2170 + t2097 * t2139) * MDP(10) + (t2057 * t2422 + t2067 * t2170 + t2098 * t2139) * MDP(11)) * t2152 + ((t2049 * t2424 + t2065 * t2169 + t2093 * t2138) * MDP(10) + (t2055 * t2424 + t2063 * t2169 + t2094 * t2138) * MDP(11)) * t2150 + t2229; (t2138 ^ 2 * t2151 + t2139 ^ 2 * t2153 + t2140 ^ 2 * t2155) * MDP(1) + (t2151 * t2163 + t2153 * t2164 + t2155 * t2162) * MDP(2) + (t2162 * t2401 + t2163 * t2415 + t2164 * t2408) * MDP(5) + (t2162 * t2322 + t2163 * t2336 + t2164 * t2329) * t2456 + MDP(12) + ((t2092 * t2168 + t2109 * t2140) * MDP(10) + (t2091 * t2168 + t2110 * t2140) * MDP(11)) * t2154 + ((t2090 * t2170 + t2107 * t2139) * MDP(10) + (t2089 * t2170 + t2108 * t2139) * MDP(11)) * t2152 + ((t2088 * t2169 + t2105 * t2138) * MDP(10) + (t2087 * t2169 + t2106 * t2138) * MDP(11)) * t2150 + (t2147 ^ 2 * t2151 + t2148 ^ 2 * t2153 + t2149 ^ 2 * t2155) * t2438 + ((t2221 * t2360 + t2223 * t2359 + t2225 * t2358) * MDP(3) + (-t2215 * t2360 - t2217 * t2359 - t2219 * t2358) * MDP(4)) * t2453 + ((t2076 * t2425 + t2077 * t2423 + t2078 * t2421) * MDP(10) + (t2079 * t2425 + t2080 * t2423 + t2081 * t2421) * MDP(11) + (t2214 * t2345 + t2216 * t2344 + t2218 * t2343) * t2455 + (t2220 * t2345 + t2222 * t2344 + t2224 * t2343) * t2454) * t2227;];
%% Postprocessing: Reshape Output
% From vec2mat_3_matlab.m
res = [t1(1), t1(2), t1(3); t1(4), t1(5), t1(6); t1(7), t1(8), t1(9);];
MMX  = res;