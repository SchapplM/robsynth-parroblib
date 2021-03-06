% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P4PRRRR1G3A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [15x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P4PRRRR1G3A0_convert_par2_MPV_fixb.m

% Output:
% taucX [4x1]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-02 19:06
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P4PRRRR1G3A0_coriolisvec_para_pf_mdp(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(2,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G3A0_coriolisvec_para_pf_mdp: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR1G3A0_coriolisvec_para_pf_mdp: xDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G3A0_coriolisvec_para_pf_mdp: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G3A0_coriolisvec_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G3A0_coriolisvec_para_pf_mdp: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G3A0_coriolisvec_para_pf_mdp: Koppelpunkt has to be [4x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P4PRRRR1G3A0_coriolisvec_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-02 19:06:10
% EndTime: 2020-03-02 19:06:23
% DurationCPUTime: 13.73s
% Computational Cost: add. (5832->409), mult. (15554->839), div. (3896->22), fcn. (11996->26), ass. (0->325)
t2266 = xP(4);
t2194 = sin(t2266);
t2195 = cos(t2266);
t2267 = koppelP(4,2);
t2271 = koppelP(4,1);
t2174 = t2194 * t2271 + t2195 * t2267;
t2178 = -t2194 * t2267 + t2195 * t2271;
t2246 = legFrame(4,2);
t2186 = sin(t2246);
t2190 = cos(t2246);
t2162 = t2174 * t2190 + t2186 * t2178;
t2242 = sin(qJ(3,4));
t2244 = cos(qJ(3,4));
t2262 = xDP(4);
t2264 = xDP(2);
t2265 = xDP(1);
t2245 = cos(qJ(2,4));
t2263 = xDP(3);
t2440 = t2245 * t2263;
t2155 = (t2162 * t2262 + t2186 * t2264 - t2265 * t2190) * t2244 - t2242 * t2440;
t2150 = t2155 ^ 2;
t2243 = sin(qJ(2,4));
t2199 = 0.1e1 / t2243 ^ 2;
t2499 = t2150 * t2199;
t2268 = koppelP(3,2);
t2272 = koppelP(3,1);
t2175 = t2194 * t2272 + t2195 * t2268;
t2179 = -t2194 * t2268 + t2195 * t2272;
t2247 = legFrame(3,2);
t2187 = sin(t2247);
t2191 = cos(t2247);
t2163 = t2175 * t2191 + t2187 * t2179;
t2250 = sin(qJ(3,3));
t2256 = cos(qJ(3,3));
t2257 = cos(qJ(2,3));
t2430 = t2257 * t2263;
t2159 = (t2163 * t2262 + t2187 * t2264 - t2265 * t2191) * t2256 - t2250 * t2430;
t2152 = t2159 ^ 2;
t2251 = sin(qJ(2,3));
t2210 = 0.1e1 / t2251 ^ 2;
t2498 = t2152 * t2210;
t2269 = koppelP(2,2);
t2273 = koppelP(2,1);
t2176 = t2194 * t2273 + t2195 * t2269;
t2180 = -t2194 * t2269 + t2195 * t2273;
t2248 = legFrame(2,2);
t2188 = sin(t2248);
t2192 = cos(t2248);
t2164 = t2176 * t2192 + t2188 * t2180;
t2252 = sin(qJ(3,2));
t2258 = cos(qJ(3,2));
t2259 = cos(qJ(2,2));
t2428 = t2259 * t2263;
t2160 = (t2164 * t2262 + t2188 * t2264 - t2265 * t2192) * t2258 - t2252 * t2428;
t2153 = t2160 ^ 2;
t2253 = sin(qJ(2,2));
t2215 = 0.1e1 / t2253 ^ 2;
t2497 = t2153 * t2215;
t2270 = koppelP(1,2);
t2274 = koppelP(1,1);
t2177 = t2194 * t2274 + t2195 * t2270;
t2181 = -t2194 * t2270 + t2195 * t2274;
t2249 = legFrame(1,2);
t2189 = sin(t2249);
t2193 = cos(t2249);
t2165 = t2177 * t2193 + t2189 * t2181;
t2254 = sin(qJ(3,1));
t2260 = cos(qJ(3,1));
t2261 = cos(qJ(2,1));
t2426 = t2261 * t2263;
t2161 = (t2165 * t2262 + t2189 * t2264 - t2265 * t2193) * t2260 - t2254 * t2426;
t2154 = t2161 ^ 2;
t2255 = sin(qJ(2,1));
t2220 = 0.1e1 / t2255 ^ 2;
t2496 = t2154 * t2220;
t2218 = t2254 ^ 2;
t2219 = 0.1e1 / t2255;
t2304 = t2260 ^ 2;
t2235 = 0.1e1 / t2304;
t2237 = 0.1e1 / t2304 ^ 2;
t2325 = (t2218 * t2237 + t2235) * t2219;
t2213 = t2252 ^ 2;
t2214 = 0.1e1 / t2253;
t2300 = t2258 ^ 2;
t2229 = 0.1e1 / t2300;
t2231 = 0.1e1 / t2300 ^ 2;
t2326 = (t2213 * t2231 + t2229) * t2214;
t2208 = t2250 ^ 2;
t2209 = 0.1e1 / t2251;
t2296 = t2256 ^ 2;
t2223 = 0.1e1 / t2296;
t2225 = 0.1e1 / t2296 ^ 2;
t2327 = (t2208 * t2225 + t2223) * t2209;
t2197 = t2242 ^ 2;
t2198 = 0.1e1 / t2243;
t2283 = t2244 ^ 2;
t2202 = 0.1e1 / t2283;
t2204 = 0.1e1 / t2283 ^ 2;
t2328 = (t2197 * t2204 + t2202) * t2198;
t2275 = 0.1e1 / pkin(2);
t2495 = 2 * MDP(6);
t2201 = 0.1e1 / t2244;
t2222 = 0.1e1 / t2256;
t2228 = 0.1e1 / t2258;
t2234 = 0.1e1 / t2260;
t2276 = 0.1e1 / pkin(2) ^ 2;
t2494 = MDP(2) * t2275;
t2166 = -t2186 * t2245 + t2190 * t2243;
t2167 = t2186 * t2243 + t2190 * t2245;
t2203 = t2201 * t2202;
t2240 = t2262 ^ 2;
t2241 = t2263 ^ 2;
t2445 = t2241 * t2275;
t2485 = t2198 * t2499;
t2138 = t2203 * t2275 * t2485 + (t2203 * t2445 + (-t2166 * t2174 - t2167 * t2178) * t2240) * t2198;
t2493 = t2138 * t2198;
t2492 = t2138 * t2201;
t2168 = -t2187 * t2257 + t2191 * t2251;
t2169 = t2187 * t2251 + t2191 * t2257;
t2224 = t2222 * t2223;
t2483 = t2209 * t2498;
t2139 = t2224 * t2275 * t2483 + (t2224 * t2445 + (-t2168 * t2175 - t2169 * t2179) * t2240) * t2209;
t2491 = t2139 * t2209;
t2490 = t2139 * t2222;
t2170 = -t2188 * t2259 + t2192 * t2253;
t2171 = t2188 * t2253 + t2192 * t2259;
t2230 = t2228 * t2229;
t2481 = t2214 * t2497;
t2140 = t2230 * t2275 * t2481 + (t2230 * t2445 + (-t2170 * t2176 - t2171 * t2180) * t2240) * t2214;
t2489 = t2140 * t2214;
t2488 = t2140 * t2228;
t2172 = -t2189 * t2261 + t2193 * t2255;
t2173 = t2189 * t2255 + t2193 * t2261;
t2236 = t2234 * t2235;
t2479 = t2219 * t2496;
t2141 = t2236 * t2275 * t2479 + (t2236 * t2445 + (-t2172 * t2177 - t2173 * t2181) * t2240) * t2219;
t2487 = t2141 * t2219;
t2486 = t2141 * t2234;
t2484 = t2150 * t2204;
t2482 = t2152 * t2225;
t2480 = t2153 * t2231;
t2478 = t2154 * t2237;
t2477 = t2155 * t2199;
t2476 = t2159 * t2210;
t2475 = t2160 * t2215;
t2474 = t2161 * t2220;
t2473 = t2194 * t2240;
t2472 = t2195 * t2240;
t2471 = t2198 * t2201;
t2470 = t2198 * t2245 ^ 2;
t2469 = t2198 * t2245;
t2468 = t2199 * t2276;
t2467 = t2201 * t2242;
t2466 = t2203 * t2242;
t2465 = t2204 * t2245;
t2464 = t2209 * t2222;
t2463 = t2209 * t2257 ^ 2;
t2462 = t2209 * t2257;
t2461 = t2210 * t2276;
t2460 = t2214 * t2228;
t2459 = t2214 * t2259 ^ 2;
t2458 = t2214 * t2259;
t2457 = t2215 * t2276;
t2456 = t2219 * t2234;
t2455 = t2219 * t2261 ^ 2;
t2454 = t2219 * t2261;
t2453 = t2220 * t2276;
t2452 = t2222 * t2250;
t2451 = t2224 * t2250;
t2450 = t2225 * t2257;
t2449 = t2230 * t2252;
t2448 = t2231 * t2259;
t2447 = t2236 * t2254;
t2446 = t2237 * t2261;
t2444 = t2241 * t2276;
t2443 = t2242 * t2243;
t2442 = t2242 * t2245;
t2383 = t2204 * t2468;
t2392 = t2155 * t2469;
t2425 = t2263 * t2275;
t2126 = -(-t2263 * t2443 + t2392) * t2155 * t2383 + (-(-t2155 * t2242 + t2440) * t2203 * t2425 + (-t2174 * t2186 + t2178 * t2190) * t2240) * t2275 * t2471;
t2441 = t2245 * t2126;
t2439 = t2250 * t2251;
t2438 = t2250 * t2257;
t2376 = t2231 * t2457;
t2388 = t2160 * t2458;
t2436 = t2252 * t2253;
t2128 = -(-t2263 * t2436 + t2388) * t2160 * t2376 + (-(-t2160 * t2252 + t2428) * t2230 * t2425 + (-t2176 * t2188 + t2180 * t2192) * t2240) * t2275 * t2460;
t2437 = t2252 * t2128;
t2435 = t2252 * t2259;
t2373 = t2237 * t2453;
t2386 = t2161 * t2454;
t2433 = t2254 * t2255;
t2129 = -(-t2263 * t2433 + t2386) * t2161 * t2373 + (-(-t2161 * t2254 + t2426) * t2236 * t2425 + (-t2177 * t2189 + t2181 * t2193) * t2240) * t2275 * t2456;
t2434 = t2254 * t2129;
t2432 = t2254 * t2261;
t2379 = t2225 * t2461;
t2390 = t2159 * t2462;
t2127 = -(-t2263 * t2439 + t2390) * t2159 * t2379 + (-(-t2159 * t2250 + t2430) * t2224 * t2425 + (-t2175 * t2187 + t2179 * t2191) * t2240) * t2275 * t2464;
t2431 = t2257 * t2127;
t2429 = t2259 * t2128;
t2427 = t2261 * t2129;
t2277 = t2275 * t2276;
t2424 = t2263 * t2277;
t2146 = t2150 * t2383;
t2142 = t2202 * t2444 + t2146;
t2415 = -0.2e1 * t2263 * t2276;
t2312 = t2392 * t2415;
t2369 = t2243 * t2444;
t2423 = -t2197 * t2203 * t2369 + t2312 * t2466 + (-t2142 * t2243 + t2441) * t2244;
t2147 = t2152 * t2379;
t2143 = t2223 * t2444 + t2147;
t2311 = t2390 * t2415;
t2368 = t2251 * t2444;
t2422 = -t2208 * t2224 * t2368 + t2311 * t2451 + (-t2143 * t2251 + t2431) * t2256;
t2148 = t2153 * t2376;
t2144 = t2229 * t2444 + t2148;
t2310 = t2388 * t2415;
t2367 = t2253 * t2444;
t2421 = -t2213 * t2230 * t2367 + t2310 * t2449 + (-t2144 * t2253 + t2429) * t2258;
t2149 = t2154 * t2373;
t2145 = t2235 * t2444 + t2149;
t2309 = t2386 * t2415;
t2366 = t2255 * t2444;
t2420 = -t2218 * t2236 * t2366 + t2309 * t2447 + (-t2145 * t2255 + t2427) * t2260;
t2419 = t2142 * t2443 + t2202 * t2312 + (-t2202 * t2369 - t2441) * t2242;
t2418 = t2143 * t2439 + t2223 * t2311 + (-t2223 * t2368 - t2431) * t2250;
t2417 = t2144 * t2436 + t2229 * t2310 + (-t2229 * t2367 - t2429) * t2252;
t2416 = t2145 * t2433 + t2235 * t2309 + (-t2235 * t2366 - t2427) * t2254;
t2414 = MDP(7) * t2241 * t2277;
t2413 = t2126 * t2471;
t2412 = t2126 * t2198 * t2242;
t2411 = t2127 * t2464;
t2410 = t2127 * t2209 * t2250;
t2409 = t2128 * t2460;
t2408 = t2214 * t2437;
t2407 = t2129 * t2456;
t2406 = t2219 * t2434;
t2405 = t2138 * t2469;
t2404 = t2139 * t2462;
t2403 = t2140 * t2458;
t2402 = t2141 * t2454;
t2401 = t2199 * t2484;
t2400 = t2242 * t2499;
t2399 = t2210 * t2482;
t2398 = t2250 * t2498;
t2397 = t2215 * t2480;
t2396 = t2252 * t2497;
t2395 = t2220 * t2478;
t2394 = t2254 * t2496;
t2393 = (-0.1e1 + 0.2e1 * t2283) * t2477;
t2391 = (-0.1e1 + 0.2e1 * t2296) * t2476;
t2389 = (-0.1e1 + 0.2e1 * t2300) * t2475;
t2387 = (-0.1e1 + 0.2e1 * t2304) * t2474;
t2385 = t2198 * t2441;
t2384 = t2199 * t2465;
t2205 = t2201 * t2204;
t2382 = t2205 * t2442;
t2381 = t2209 * t2431;
t2380 = t2210 * t2450;
t2378 = t2214 * t2429;
t2377 = t2215 * t2448;
t2375 = t2219 * t2427;
t2374 = t2220 * t2446;
t2226 = t2222 * t2225;
t2372 = t2226 * t2438;
t2232 = t2228 * t2231;
t2371 = t2232 * t2435;
t2238 = t2234 * t2237;
t2370 = t2238 * t2432;
t2365 = 0.2e1 * t2424;
t2364 = t2423 * t2198;
t2363 = t2422 * t2209;
t2362 = t2421 * t2214;
t2361 = t2420 * t2219;
t2360 = t2419 * t2198;
t2359 = t2418 * t2209;
t2358 = t2417 * t2214;
t2357 = t2416 * t2219;
t2352 = t2202 * t2385;
t2351 = t2223 * t2381;
t2350 = t2213 * t2409;
t2349 = t2229 * t2378;
t2348 = t2218 * t2407;
t2347 = t2235 * t2375;
t2346 = t2201 * t2405;
t2345 = t2222 * t2404;
t2344 = t2228 * t2403;
t2343 = t2234 * t2402;
t2342 = t2465 * t2485;
t2341 = t2450 * t2483;
t2340 = t2448 * t2481;
t2339 = t2446 * t2479;
t2338 = t2204 * t2393;
t2337 = t2466 * t2477;
t2336 = t2225 * t2391;
t2335 = t2451 * t2476;
t2334 = t2231 * t2389;
t2333 = t2449 * t2475;
t2332 = t2237 * t2387;
t2331 = t2447 * t2474;
t2330 = t2197 * t2413;
t2329 = t2208 * t2411;
t2324 = t2162 * t2346;
t2323 = t2163 * t2345;
t2322 = t2164 * t2344;
t2321 = t2165 * t2343;
t2320 = t2186 * t2346;
t2319 = t2187 * t2345;
t2318 = t2188 * t2344;
t2317 = t2189 * t2343;
t2316 = t2190 * t2346;
t2315 = t2191 * t2345;
t2314 = t2192 * t2344;
t2313 = t2193 * t2343;
t2308 = t2126 * t2467 + t2127 * t2452 + t2228 * t2437 + t2234 * t2434;
t2217 = t2254 * t2218;
t2212 = t2252 * t2213;
t2207 = t2250 * t2208;
t2196 = t2242 * t2197;
t2158 = (t2172 * t2181 - t2173 * t2177) * t2219;
t2157 = (t2170 * t2180 - t2171 * t2176) * t2214;
t2156 = (t2168 * t2179 - t2169 * t2175) * t2209;
t2151 = (t2166 * t2178 - t2167 * t2174) * t2198;
t1 = [(t2167 * t2493 + t2169 * t2491 + t2171 * t2489 + t2173 * t2487) * MDP(1) + (-t2190 * t2413 - t2191 * t2411 - t2192 * t2409 - t2193 * t2407) * t2494 + (t2167 * t2385 + t2169 * t2381 + t2171 * t2378 + t2173 * t2375 + (-t2167 * t2401 - t2169 * t2399 - t2171 * t2397 - t2173 * t2395) * t2276 + (-t2313 - t2314 - t2315 - t2316) * t2275) * MDP(3) + (-t2167 * t2126 - t2169 * t2127 - t2171 * t2128 - t2173 * t2129 + (-t2167 * t2342 - t2169 * t2341 - t2171 * t2340 - t2173 * t2339) * t2276 + (t2190 * t2492 + t2191 * t2490 + t2192 * t2488 + t2193 * t2486) * t2275) * MDP(4) + ((-t2190 * t2330 - t2191 * t2329 - t2192 * t2350 - t2193 * t2348) * t2275 + (-t2190 * t2337 - t2191 * t2335 - t2192 * t2333 - t2193 * t2331) * t2365) * MDP(5) + ((-t2190 * t2412 - t2191 * t2410 - t2192 * t2408 - t2193 * t2406) * t2275 + (-t2190 * t2338 - t2191 * t2336 - t2192 * t2334 - t2193 * t2332) * t2424) * t2495 + (-t2190 * t2328 - t2191 * t2327 - t2192 * t2326 - t2193 * t2325) * t2414 + (t2173 * t2361 + t2171 * t2362 + t2169 * t2363 + t2167 * t2364 + (-t2190 * t2405 - t2191 * t2404 - t2192 * t2403 - t2193 * t2402) * t2275) * MDP(10) + (t2173 * t2357 + t2171 * t2358 + t2169 * t2359 + t2167 * t2360 + (t2242 * t2316 + t2250 * t2315 + t2252 * t2314 + t2254 * t2313) * t2275) * MDP(11) - MDP(13) * t2472 + MDP(14) * t2473; (t2166 * t2493 + t2168 * t2491 + t2170 * t2489 + t2172 * t2487) * MDP(1) + (t2186 * t2413 + t2187 * t2411 + t2188 * t2409 + t2189 * t2407) * t2494 + (t2166 * t2385 + t2168 * t2381 + t2170 * t2378 + t2172 * t2375 + (-t2166 * t2401 - t2168 * t2399 - t2170 * t2397 - t2172 * t2395) * t2276 + (t2317 + t2318 + t2319 + t2320) * t2275) * MDP(3) + (-t2166 * t2126 - t2168 * t2127 - t2170 * t2128 - t2172 * t2129 + (-t2166 * t2342 - t2168 * t2341 - t2170 * t2340 - t2172 * t2339) * t2276 + (-t2186 * t2492 - t2187 * t2490 - t2188 * t2488 - t2189 * t2486) * t2275) * MDP(4) + ((t2186 * t2330 + t2187 * t2329 + t2188 * t2350 + t2189 * t2348) * t2275 + (t2186 * t2337 + t2187 * t2335 + t2188 * t2333 + t2189 * t2331) * t2365) * MDP(5) + ((t2186 * t2412 + t2187 * t2410 + t2188 * t2408 + t2189 * t2406) * t2275 + (t2186 * t2338 + t2187 * t2336 + t2188 * t2334 + t2189 * t2332) * t2424) * t2495 + (t2186 * t2328 + t2187 * t2327 + t2188 * t2326 + t2189 * t2325) * t2414 + (t2172 * t2361 + t2170 * t2362 + t2168 * t2363 + t2166 * t2364 + (t2186 * t2405 + t2187 * t2404 + t2188 * t2403 + t2189 * t2402) * t2275) * MDP(10) + (t2172 * t2357 + t2170 * t2358 + t2168 * t2359 + t2166 * t2360 + (-t2242 * t2320 - t2250 * t2319 - t2252 * t2318 - t2254 * t2317) * t2275) * MDP(11) - MDP(13) * t2473 - MDP(14) * t2472; -t2308 * MDP(4) + (t2141 * MDP(1) + t2420 * MDP(10) + t2416 * MDP(11) + MDP(3) * t2427) * t2254 * t2456 + (t2140 * MDP(1) + t2421 * MDP(10) + t2417 * MDP(11) + MDP(3) * t2429) * t2252 * t2460 + (t2139 * MDP(1) + t2422 * MDP(10) + t2418 * MDP(11) + MDP(3) * t2431) * t2209 * t2452 + (t2138 * MDP(1) + t2423 * MDP(10) + t2419 * MDP(11) + MDP(3) * t2441) * t2198 * t2467 + ((-t2205 * t2400 - t2226 * t2398 - t2232 * t2396 - t2238 * t2394) * MDP(3) + (-t2370 * t2479 - t2371 * t2481 - t2372 * t2483 - t2382 * t2485) * MDP(4)) * t2276 + ((-t2204 * t2400 - t2225 * t2398 - t2231 * t2396 - t2237 * t2394) * MDP(5) + 0.2e1 * ((-t2155 * t2197 * t2384 - t2159 * t2208 * t2380 - t2160 * t2213 * t2377 - t2161 * t2218 * t2374) * MDP(5) + (-t2370 * t2387 - t2371 * t2389 - t2372 * t2391 - t2382 * t2393) * MDP(6)) * t2263 + ((-t2196 * t2205 * t2469 - t2198 * t2203 * t2442 - t2207 * t2226 * t2462 - t2209 * t2224 * t2438 - t2212 * t2232 * t2458 - t2214 * t2230 * t2435 - t2217 * t2238 * t2454 - t2219 * t2236 * t2432) * MDP(7) + (t2204 * t2242 + t2225 * t2250 + t2231 * t2252 + t2237 * t2254) * MDP(9)) * t2241) * t2277 + ((-t2242 * t2352 - t2250 * t2351 - t2252 * t2349 - t2254 * t2347) * MDP(2) + (-t2196 * t2352 - t2207 * t2351 - t2212 * t2349 - t2217 * t2347) * MDP(5) + (t2234 * (-0.2e1 * t2154 * t2235 * t2453 + t2149) + t2228 * (-0.2e1 * t2153 * t2229 * t2457 + t2148) + t2222 * (-0.2e1 * t2152 * t2223 * t2461 + t2147) + t2201 * (-0.2e1 * t2150 * t2202 * t2468 + t2146) - 0.2e1 * t2261 * t2348 - 0.2e1 * t2259 * t2350 - 0.2e1 * t2257 * t2329 - 0.2e1 * t2245 * t2330) * MDP(6) + t2308 * MDP(7) + (t2126 + t2127 + t2128 + t2129) * MDP(8) + ((t2218 * t2235 * t2455 - t2255) * MDP(11) + ((-MDP(3) * t2455 + t2261 * MDP(4)) * t2235 + (-t2255 - t2455) * MDP(10) * t2234) * t2254) * t2141 + ((t2213 * t2229 * t2459 - t2253) * MDP(11) + ((-MDP(3) * t2459 + t2259 * MDP(4)) * t2229 + (-t2253 - t2459) * MDP(10) * t2228) * t2252) * t2140 + ((t2208 * t2223 * t2463 - t2251) * MDP(11) + ((-MDP(3) * t2463 + t2257 * MDP(4)) * t2223 + (-t2251 - t2463) * MDP(10) * t2222) * t2250) * t2139 + ((t2197 * t2202 * t2470 - t2243) * MDP(11) + ((-MDP(3) * t2470 + t2245 * MDP(4)) * t2202 + (-t2243 - t2470) * MDP(10) * t2201) * t2242) * t2138) * t2275; (t2138 * t2151 + t2139 * t2156 + t2140 * t2157 + t2141 * t2158) * MDP(1) + (t2162 * t2413 + t2163 * t2411 + t2164 * t2409 + t2165 * t2407) * t2494 + (t2151 * t2441 + t2156 * t2431 + t2157 * t2429 + t2158 * t2427 + (-t2151 * t2198 * t2484 - t2156 * t2209 * t2482 - t2157 * t2214 * t2480 - t2158 * t2219 * t2478) * t2276 + (t2321 + t2322 + t2323 + t2324) * t2275) * MDP(3) + (-t2151 * t2126 * t2243 - t2156 * t2127 * t2251 - t2157 * t2128 * t2253 - t2158 * t2129 * t2255 + (-t2151 * t2150 * t2384 - t2156 * t2152 * t2380 - t2157 * t2153 * t2377 - t2158 * t2154 * t2374) * t2276 + (-t2162 * t2492 - t2163 * t2490 - t2164 * t2488 - t2165 * t2486) * t2275) * MDP(4) + ((t2162 * t2330 + t2163 * t2329 + t2164 * t2350 + t2165 * t2348) * t2275 + (t2162 * t2337 + t2163 * t2335 + t2164 * t2333 + t2165 * t2331) * t2365) * MDP(5) + ((t2162 * t2412 + t2163 * t2410 + t2164 * t2408 + t2165 * t2406) * t2275 + (t2162 * t2338 + t2163 * t2336 + t2164 * t2334 + t2165 * t2332) * t2424) * t2495 + (t2162 * t2328 + t2163 * t2327 + t2164 * t2326 + t2165 * t2325) * t2414 + (t2420 * t2158 + t2421 * t2157 + t2422 * t2156 + t2423 * t2151 + (t2162 * t2405 + t2163 * t2404 + t2164 * t2403 + t2165 * t2402) * t2275) * MDP(10) + (t2416 * t2158 + t2417 * t2157 + t2418 * t2156 + t2419 * t2151 + (-t2242 * t2324 - t2250 * t2323 - t2252 * t2322 - t2254 * t2321) * t2275) * MDP(11);];
taucX  = t1;
