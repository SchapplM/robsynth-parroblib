% Calculate Gravitation load for parallel robot
% P3RRRRR7V2G2A0
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
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 10:08
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRRRR7V2G2A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 09:43:23
% EndTime: 2020-08-07 09:43:25
% DurationCPUTime: 2.89s
% Computational Cost: add. (1191->320), mult. (2187->448), div. (66->14), fcn. (1239->66), ass. (0->212)
t2463 = pkin(6) + pkin(5);
t2368 = (-rSges(2,3) - pkin(5)) * m(2) + m(1) * rSges(1,2) - (rSges(3,3) + t2463) * m(3);
t2474 = m(1) * rSges(1,1) + (m(2) + m(3)) * pkin(1);
t2326 = sin(qJ(3,1));
t2335 = cos(qJ(3,1));
t2393 = m(2) * rSges(2,1) + pkin(2) * m(3);
t2237 = (rSges(3,1) * t2335 - rSges(3,2) * t2326) * m(3) + t2393;
t2338 = m(2) * rSges(2,2);
t2243 = t2338 + (rSges(3,1) * t2326 + rSges(3,2) * t2335) * m(3);
t2327 = sin(qJ(2,1));
t2336 = cos(qJ(2,1));
t2471 = t2237 * t2336 - t2243 * t2327;
t2323 = sin(qJ(3,2));
t2332 = cos(qJ(3,2));
t2236 = (rSges(3,1) * t2332 - rSges(3,2) * t2323) * m(3) + t2393;
t2242 = t2338 + (rSges(3,1) * t2323 + rSges(3,2) * t2332) * m(3);
t2324 = sin(qJ(2,2));
t2333 = cos(qJ(2,2));
t2470 = t2236 * t2333 - t2242 * t2324;
t2320 = sin(qJ(3,3));
t2329 = cos(qJ(3,3));
t2235 = (rSges(3,1) * t2329 - rSges(3,2) * t2320) * m(3) + t2393;
t2241 = t2338 + (rSges(3,1) * t2320 + rSges(3,2) * t2329) * m(3);
t2321 = sin(qJ(2,3));
t2330 = cos(qJ(2,3));
t2469 = t2235 * t2330 - t2241 * t2321;
t2468 = 2 * pkin(3);
t2302 = (pkin(7) + t2463);
t2467 = 2 * t2302;
t2307 = t2330 ^ 2;
t2466 = 0.2e1 * t2307;
t2309 = t2333 ^ 2;
t2465 = 0.2e1 * t2309;
t2311 = t2336 ^ 2;
t2464 = 0.2e1 * t2311;
t2462 = m(3) / pkin(3);
t2306 = t2329 ^ 2;
t2294 = pkin(3) * t2306;
t2308 = t2332 ^ 2;
t2295 = pkin(3) * t2308;
t2310 = t2335 ^ 2;
t2296 = pkin(3) * t2310;
t2461 = t2320 * pkin(1);
t2460 = t2320 * pkin(3);
t2459 = t2323 * pkin(1);
t2458 = t2323 * pkin(3);
t2457 = t2326 * pkin(1);
t2456 = t2326 * pkin(3);
t2455 = t2329 * pkin(2);
t2291 = t2329 * pkin(3);
t2454 = t2332 * pkin(2);
t2292 = t2332 * pkin(3);
t2453 = t2335 * pkin(2);
t2293 = t2335 * pkin(3);
t2452 = -qJ(3,1) + qJ(1,1);
t2451 = qJ(3,1) + qJ(1,1);
t2450 = -qJ(3,2) + qJ(1,2);
t2449 = qJ(3,2) + qJ(1,2);
t2448 = -qJ(3,3) + qJ(1,3);
t2447 = qJ(3,3) + qJ(1,3);
t2446 = qJ(1,1) + 0.2e1 * qJ(2,1);
t2445 = qJ(1,1) - 0.2e1 * qJ(2,1);
t2444 = qJ(1,2) + 0.2e1 * qJ(2,2);
t2443 = qJ(1,2) - 0.2e1 * qJ(2,2);
t2442 = qJ(1,3) + 0.2e1 * qJ(2,3);
t2441 = qJ(1,3) - 0.2e1 * qJ(2,3);
t2440 = 0.2e1 * pkin(1);
t2322 = sin(qJ(1,3));
t2331 = cos(qJ(1,3));
t2366 = pkin(1) * t2322 - t2331 * t2302;
t2400 = t2320 * t2321;
t2220 = t2366 * t2400 + (t2306 - 0.1e1) * t2322 * pkin(3);
t2317 = legFrame(3,2);
t2285 = sin(t2317);
t2439 = t2220 * t2285;
t2288 = cos(t2317);
t2438 = t2220 * t2288;
t2325 = sin(qJ(1,2));
t2334 = cos(qJ(1,2));
t2365 = pkin(1) * t2325 - t2334 * t2302;
t2399 = t2323 * t2324;
t2221 = t2365 * t2399 + (t2308 - 0.1e1) * t2325 * pkin(3);
t2318 = legFrame(2,2);
t2286 = sin(t2318);
t2437 = t2221 * t2286;
t2289 = cos(t2318);
t2436 = t2221 * t2289;
t2328 = sin(qJ(1,1));
t2337 = cos(qJ(1,1));
t2364 = pkin(1) * t2328 - t2337 * t2302;
t2398 = t2326 * t2327;
t2222 = t2364 * t2398 + (t2310 - 0.1e1) * t2328 * pkin(3);
t2319 = legFrame(1,2);
t2287 = sin(t2319);
t2435 = t2222 * t2287;
t2290 = cos(t2319);
t2434 = t2222 * t2290;
t2385 = pkin(3) * t2400;
t2269 = t2291 + pkin(2);
t2417 = t2269 * t2330;
t2433 = 0.1e1 / (pkin(1) - t2385 + t2417) / t2320;
t2384 = pkin(3) * t2399;
t2270 = t2292 + pkin(2);
t2415 = t2270 * t2333;
t2432 = 0.1e1 / (pkin(1) - t2384 + t2415) / t2323;
t2383 = pkin(3) * t2398;
t2271 = t2293 + pkin(2);
t2413 = t2271 * t2336;
t2431 = 0.1e1 / (pkin(1) - t2383 + t2413) / t2326;
t2250 = t2288 * g(1) - t2285 * g(2);
t2312 = qJ(2,3) + qJ(3,3);
t2276 = cos(t2312);
t2367 = t2474 * g(3);
t2369 = t2368 * g(3);
t2424 = 0.1e1 / (t2330 * pkin(2) + pkin(3) * t2276 + pkin(1)) * (t2369 * t2331 + t2322 * (t2469 * g(3) + t2367) + ((-t2474 - t2469) * t2331 + t2322 * t2368) * t2250);
t2251 = t2289 * g(1) - t2286 * g(2);
t2313 = qJ(2,2) + qJ(3,2);
t2277 = cos(t2313);
t2423 = 0.1e1 / (t2333 * pkin(2) + pkin(3) * t2277 + pkin(1)) * (t2369 * t2334 + t2325 * (t2470 * g(3) + t2367) + ((-t2474 - t2470) * t2334 + t2325 * t2368) * t2251);
t2252 = t2290 * g(1) - t2287 * g(2);
t2314 = qJ(2,1) + qJ(3,1);
t2278 = cos(t2314);
t2422 = 0.1e1 / (t2336 * pkin(2) + pkin(3) * t2278 + pkin(1)) * (t2369 * t2337 + t2328 * (t2471 * g(3) + t2367) + ((-t2474 - t2471) * t2337 + t2328 * t2368) * t2252);
t2344 = pkin(2) / 0.2e1;
t2421 = (t2291 + t2344) * t2320;
t2420 = (t2292 + t2344) * t2323;
t2419 = (t2293 + t2344) * t2326;
t2418 = t2269 * t2285;
t2416 = t2270 * t2286;
t2414 = t2271 * t2287;
t2412 = t2285 * t2322;
t2411 = t2286 * t2325;
t2410 = t2287 * t2328;
t2409 = t2288 * t2269;
t2408 = t2288 * t2322;
t2407 = t2289 * t2270;
t2406 = t2289 * t2325;
t2405 = t2290 * t2271;
t2404 = t2290 * t2328;
t2354 = pkin(3) ^ 2;
t2403 = t2306 * t2354;
t2402 = t2308 * t2354;
t2401 = t2310 * t2354;
t2259 = pkin(1) * t2321 - t2460;
t2397 = t2329 * t2259;
t2260 = pkin(1) * t2324 - t2458;
t2396 = t2332 * t2260;
t2261 = pkin(1) * t2327 - t2456;
t2395 = t2335 * t2261;
t2356 = pkin(2) ^ 2;
t2394 = -t2354 / 0.2e1 + t2356 / 0.2e1;
t2392 = pkin(2) * t2291;
t2391 = pkin(2) * t2292;
t2390 = pkin(2) * t2293;
t2388 = t2269 * t2460;
t2387 = t2270 * t2458;
t2386 = t2271 * t2456;
t2247 = t2285 * g(1) + t2288 * g(2);
t2363 = g(3) * t2331 + t2250 * t2322;
t2214 = (t2363 * t2235 + t2247 * t2241) * t2321 + (-t2247 * t2235 + t2363 * t2241) * t2330;
t2381 = t2214 * t2433;
t2248 = t2286 * g(1) + t2289 * g(2);
t2362 = g(3) * t2334 + t2251 * t2325;
t2215 = (t2362 * t2236 + t2248 * t2242) * t2324 + (-t2248 * t2236 + t2362 * t2242) * t2333;
t2380 = t2215 * t2432;
t2249 = t2287 * g(1) + t2290 * g(2);
t2361 = g(3) * t2337 + t2252 * t2328;
t2216 = (t2361 * t2237 + t2249 * t2243) * t2327 + (-t2249 * t2237 + t2361 * t2243) * t2336;
t2379 = t2216 * t2431;
t2273 = sin(t2312);
t2378 = ((-rSges(3,1) * t2247 + t2363 * rSges(3,2)) * t2276 + t2273 * (t2363 * rSges(3,1) + rSges(3,2) * t2247)) * t2433;
t2274 = sin(t2313);
t2377 = ((-rSges(3,1) * t2248 + t2362 * rSges(3,2)) * t2277 + t2274 * (t2362 * rSges(3,1) + rSges(3,2) * t2248)) * t2432;
t2275 = sin(t2314);
t2376 = ((-rSges(3,1) * t2249 + t2361 * rSges(3,2)) * t2278 + t2275 * (t2361 * rSges(3,1) + rSges(3,2) * t2249)) * t2431;
t2375 = t2322 * t2400;
t2374 = t2325 * t2399;
t2373 = t2328 * t2398;
t2372 = t2331 * t2424;
t2371 = t2334 * t2423;
t2370 = t2337 * t2422;
t2229 = t2375 * t2468 - t2366;
t2360 = pkin(2) * t2375 + t2229 * t2329;
t2230 = t2374 * t2468 - t2365;
t2359 = pkin(2) * t2374 + t2230 * t2332;
t2231 = t2373 * t2468 - t2364;
t2358 = pkin(2) * t2373 + t2231 * t2335;
t2357 = 0.1e1 / pkin(2);
t2351 = 0.2e1 * qJ(3,1);
t2348 = 0.2e1 * qJ(3,2);
t2345 = 0.2e1 * qJ(3,3);
t2343 = -pkin(3) / 0.2e1;
t2300 = -t2354 + t2356;
t2284 = -qJ(2,1) + t2452;
t2283 = qJ(2,1) + t2451;
t2282 = -qJ(2,2) + t2450;
t2281 = qJ(2,2) + t2449;
t2280 = -qJ(2,3) + t2448;
t2279 = qJ(2,3) + t2447;
t2255 = t2296 + t2453 / 0.2e1 + t2343;
t2254 = t2295 + t2454 / 0.2e1 + t2343;
t2253 = t2294 + t2455 / 0.2e1 + t2343;
t2240 = t2390 + t2394 + t2401;
t2239 = t2391 + t2394 + t2402;
t2238 = t2392 + t2394 + t2403;
t2228 = t2457 + (-pkin(3) + t2453 + 0.2e1 * t2296) * t2327;
t2227 = t2459 + (-pkin(3) + t2454 + 0.2e1 * t2295) * t2324;
t2226 = t2461 + (-pkin(3) + t2455 + 0.2e1 * t2294) * t2321;
t2225 = pkin(1) * t2456 + (t2300 + 0.2e1 * t2390 + 0.2e1 * t2401) * t2327;
t2224 = pkin(1) * t2458 + (t2300 + 0.2e1 * t2391 + 0.2e1 * t2402) * t2324;
t2223 = pkin(1) * t2460 + (t2300 + 0.2e1 * t2392 + 0.2e1 * t2403) * t2321;
t1 = [t2288 * t2372 + t2289 * t2371 + t2290 * t2370 - m(4) * g(1) + (((t2255 * t2404 + t2287 * t2419) * t2464 + (t2287 * t2228 - t2358 * t2290) * t2336 - t2434 + t2287 * t2395) * t2379 + ((t2254 * t2406 + t2286 * t2420) * t2465 + (t2286 * t2227 - t2359 * t2289) * t2333 - t2436 + t2286 * t2396) * t2380 + ((t2253 * t2408 + t2285 * t2421) * t2466 + (t2285 * t2226 - t2360 * t2288) * t2330 - t2438 + t2285 * t2397) * t2381 + (((-t2240 * t2404 - t2287 * t2386) * t2464 + (-t2287 * t2225 + t2231 * t2405) * t2336 + pkin(3) * t2434 - t2261 * t2414) * t2376 + ((-t2239 * t2406 - t2286 * t2387) * t2465 + (-t2286 * t2224 + t2230 * t2407) * t2333 + pkin(3) * t2436 - t2260 * t2416) * t2377 + ((-t2238 * t2408 - t2285 * t2388) * t2466 + (-t2285 * t2223 + t2229 * t2409) * t2330 + pkin(3) * t2438 - t2259 * t2418) * t2378) * t2462) * t2357; -t2285 * t2372 - t2286 * t2371 - t2287 * t2370 - m(4) * g(2) + (((-t2255 * t2410 + t2290 * t2419) * t2464 + (t2290 * t2228 + t2358 * t2287) * t2336 + t2435 + t2290 * t2395) * t2379 + ((-t2254 * t2411 + t2289 * t2420) * t2465 + (t2289 * t2227 + t2359 * t2286) * t2333 + t2437 + t2289 * t2396) * t2380 + ((-t2253 * t2412 + t2288 * t2421) * t2466 + (t2288 * t2226 + t2360 * t2285) * t2330 + t2439 + t2288 * t2397) * t2381 + (((t2240 * t2410 - t2290 * t2386) * t2464 + (-t2225 * t2290 - t2231 * t2414) * t2336 - pkin(3) * t2435 - t2261 * t2405) * t2376 + ((t2239 * t2411 - t2289 * t2387) * t2465 + (-t2224 * t2289 - t2230 * t2416) * t2333 - pkin(3) * t2437 - t2260 * t2407) * t2377 + ((t2238 * t2412 - t2288 * t2388) * t2466 + (-t2223 * t2288 - t2229 * t2418) * t2330 - pkin(3) * t2439 - t2259 * t2409) * t2378) * t2462) * t2357; -m(4) * g(3) - t2322 * t2424 - t2325 * t2423 - t2328 * t2422 + ((-0.2e1 * t2240 * t2337 * t2311 - ((pkin(1) - 0.2e1 * t2383) * t2337 + t2328 * t2302) * t2413 + pkin(3) * ((pkin(1) * t2398 - pkin(3) + t2296) * t2337 + t2302 * t2373)) * t2376 + (-0.2e1 * t2239 * t2334 * t2309 - ((pkin(1) - 0.2e1 * t2384) * t2334 + t2325 * t2302) * t2415 + pkin(3) * ((pkin(1) * t2399 - pkin(3) + t2295) * t2334 + t2302 * t2374)) * t2377 + (-0.2e1 * t2238 * t2331 * t2307 - ((pkin(1) - 0.2e1 * t2385) * t2331 + t2322 * t2302) * t2417 + pkin(3) * ((pkin(1) * t2400 - pkin(3) + t2294) * t2331 + t2302 * t2375)) * t2378) * t2357 * t2462 + ((cos(t2284) + cos(t2283)) * t2440 + (sin(t2284) + sin(t2283)) * t2467 + (cos(-0.2e1 * qJ(3,1) + t2445) + cos(t2351 + t2446) + 0.2e1 * t2337) * pkin(3) + (cos(-qJ(3,1) + t2445) + cos(qJ(3,1) + t2446) + cos(t2452) + cos(t2451)) * pkin(2)) / (-t2356 * sin(qJ(2,1) - qJ(3,1)) + pkin(2) * (0.2e1 * t2457 + pkin(2) * t2275 + (sin(t2351 + qJ(2,1)) - t2327) * pkin(3))) * t2216 / 0.2e1 + ((cos(t2282) + cos(t2281)) * t2440 + (sin(t2282) + sin(t2281)) * t2467 + (cos(-0.2e1 * qJ(3,2) + t2443) + cos(t2348 + t2444) + 0.2e1 * t2334) * pkin(3) + (cos(-qJ(3,2) + t2443) + cos(qJ(3,2) + t2444) + cos(t2450) + cos(t2449)) * pkin(2)) / (-t2356 * sin(qJ(2,2) - qJ(3,2)) + pkin(2) * (0.2e1 * t2459 + pkin(2) * t2274 + (sin(t2348 + qJ(2,2)) - t2324) * pkin(3))) * t2215 / 0.2e1 + ((cos(t2280) + cos(t2279)) * t2440 + (sin(t2280) + sin(t2279)) * t2467 + (cos(-0.2e1 * qJ(3,3) + t2441) + cos(t2345 + t2442) + 0.2e1 * t2331) * pkin(3) + (cos(-qJ(3,3) + t2441) + cos(qJ(3,3) + t2442) + cos(t2448) + cos(t2447)) * pkin(2)) / (-t2356 * sin(qJ(2,3) - qJ(3,3)) + pkin(2) * (0.2e1 * t2461 + pkin(2) * t2273 + (sin(t2345 + qJ(2,3)) - t2321) * pkin(3))) * t2214 / 0.2e1;];
taugX  = t1;
