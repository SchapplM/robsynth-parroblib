% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RPRRR6V1G3A0
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [12x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPRRR6V1G3A0_convert_par2_MPV_fixb.m

% Output:
% MMX [3x3]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:43
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P3RPRRR6V1G3A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G3A0_inertia_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G3A0_inertia_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G3A0_inertia_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G3A0_inertia_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G3A0_inertia_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3RPRRR6V1G3A0_inertia_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:42:47
% EndTime: 2020-08-06 18:42:55
% DurationCPUTime: 7.90s
% Computational Cost: add. (8262->448), mult. (9696->903), div. (1260->20), fcn. (8127->53), ass. (0->377)
t2338 = sin(pkin(7));
t2357 = -pkin(6) - pkin(5);
t2294 = t2338 * t2357 - pkin(1);
t2350 = cos(qJ(1,3));
t2279 = t2294 * t2350;
t2339 = cos(pkin(7));
t2344 = sin(qJ(1,3));
t2496 = t2344 * t2357;
t2499 = t2338 * t2344;
t2605 = t2279 + pkin(2) * t2499 - (pkin(2) * t2350 - t2496) * t2339;
t2352 = cos(qJ(1,2));
t2280 = t2294 * t2352;
t2346 = sin(qJ(1,2));
t2495 = t2346 * t2357;
t2498 = t2338 * t2346;
t2604 = t2280 + pkin(2) * t2498 - (pkin(2) * t2352 - t2495) * t2339;
t2354 = cos(qJ(1,1));
t2281 = t2294 * t2354;
t2348 = sin(qJ(1,1));
t2494 = t2348 * t2357;
t2497 = t2338 * t2348;
t2603 = t2281 + pkin(2) * t2497 - (pkin(2) * t2354 - t2494) * t2339;
t2486 = 2 * MDP(6);
t2602 = MDP(10) / 0.2e1;
t2601 = MDP(11) / 0.2e1;
t2349 = cos(qJ(3,3));
t2304 = t2349 * pkin(3) + pkin(2);
t2313 = t2339 * pkin(1);
t2291 = t2313 + t2304;
t2283 = 0.1e1 / t2291 ^ 2;
t2600 = t2283 / 0.4e1;
t2351 = cos(qJ(3,2));
t2305 = t2351 * pkin(3) + pkin(2);
t2292 = t2313 + t2305;
t2285 = 0.1e1 / t2292 ^ 2;
t2599 = t2285 / 0.4e1;
t2353 = cos(qJ(3,1));
t2306 = t2353 * pkin(3) + pkin(2);
t2293 = t2313 + t2306;
t2287 = 0.1e1 / t2293 ^ 2;
t2598 = t2287 / 0.4e1;
t2286 = 0.1e1 / t2293;
t2325 = qJ(1,1) + pkin(7);
t2312 = cos(t2325);
t2489 = pkin(7) + qJ(3,1);
t2492 = -pkin(7) + qJ(3,1);
t2574 = 0.2e1 * t2357;
t2578 = 0.2e1 * pkin(2);
t2579 = 0.2e1 * pkin(1);
t2257 = t2312 * t2574 + sin(t2325) * t2578 + t2348 * t2579 + (sin(qJ(1,1) - t2492) + sin(qJ(1,1) + t2489)) * pkin(3);
t2347 = sin(qJ(3,1));
t2575 = 0.2e1 * t2347;
t2263 = pkin(3) * sin(0.2e1 * qJ(3,1)) + pkin(2) * t2575 + (sin(t2489) + sin(t2492)) * pkin(1);
t2557 = t2257 / t2263;
t2450 = t2353 * t2557;
t2397 = t2286 * t2450;
t2597 = t2397 / 0.2e1;
t2451 = t2286 * t2557;
t2398 = t2347 * t2451;
t2596 = t2398 / 0.2e1;
t2284 = 0.1e1 / t2292;
t2324 = qJ(1,2) + pkin(7);
t2311 = cos(t2324);
t2488 = pkin(7) + qJ(3,2);
t2491 = -pkin(7) + qJ(3,2);
t2256 = t2311 * t2574 + sin(t2324) * t2578 + t2346 * t2579 + (sin(qJ(1,2) - t2491) + sin(qJ(1,2) + t2488)) * pkin(3);
t2345 = sin(qJ(3,2));
t2576 = 0.2e1 * t2345;
t2262 = pkin(3) * sin(0.2e1 * qJ(3,2)) + pkin(2) * t2576 + (sin(t2488) + sin(t2491)) * pkin(1);
t2558 = t2256 / t2262;
t2452 = t2351 * t2558;
t2400 = t2284 * t2452;
t2595 = t2400 / 0.2e1;
t2453 = t2284 * t2558;
t2401 = t2345 * t2453;
t2594 = t2401 / 0.2e1;
t2282 = 0.1e1 / t2291;
t2323 = qJ(1,3) + pkin(7);
t2310 = cos(t2323);
t2487 = pkin(7) + qJ(3,3);
t2490 = -pkin(7) + qJ(3,3);
t2255 = t2310 * t2574 + sin(t2323) * t2578 + t2344 * t2579 + (sin(qJ(1,3) - t2490) + sin(qJ(1,3) + t2487)) * pkin(3);
t2343 = sin(qJ(3,3));
t2577 = 0.2e1 * t2343;
t2261 = pkin(3) * sin(0.2e1 * qJ(3,3)) + pkin(2) * t2577 + (sin(t2487) + sin(t2490)) * pkin(1);
t2559 = t2255 / t2261;
t2454 = t2349 * t2559;
t2403 = t2282 * t2454;
t2593 = t2403 / 0.2e1;
t2455 = t2282 * t2559;
t2404 = t2343 * t2455;
t2592 = t2404 / 0.2e1;
t2591 = -t2287 * t2312 / 0.2e1;
t2590 = -t2285 * t2311 / 0.2e1;
t2589 = -t2283 * t2310 / 0.2e1;
t2342 = legFrame(1,2);
t2300 = t2342 + t2325;
t2301 = -t2342 + t2325;
t2275 = cos(t2301) - cos(t2300);
t2551 = t2275 * t2286;
t2588 = t2551 / 0.2e1;
t2341 = legFrame(2,2);
t2298 = t2341 + t2324;
t2299 = -t2341 + t2324;
t2274 = cos(t2299) - cos(t2298);
t2552 = t2274 * t2284;
t2587 = t2552 / 0.2e1;
t2340 = legFrame(3,2);
t2296 = t2340 + t2323;
t2297 = -t2340 + t2323;
t2273 = cos(t2297) - cos(t2296);
t2553 = t2273 * t2282;
t2586 = t2553 / 0.2e1;
t2272 = -sin(t2300) - sin(t2301);
t2554 = t2272 * t2286;
t2585 = t2554 / 0.2e1;
t2271 = -sin(t2298) - sin(t2299);
t2555 = t2271 * t2284;
t2584 = t2555 / 0.2e1;
t2270 = -sin(t2296) - sin(t2297);
t2556 = t2270 * t2282;
t2583 = t2556 / 0.2e1;
t2582 = t2270 * t2273;
t2581 = t2271 * t2274;
t2580 = t2272 * t2275;
t2573 = MDP(5) / 0.2e1;
t2572 = MDP(5) / 0.4e1;
t2358 = 0.1e1 / pkin(3);
t2571 = MDP(7) * t2358;
t2570 = MDP(8) * t2358;
t2569 = MDP(9) / pkin(3) ^ 2;
t2511 = t2304 * t2350;
t2250 = (-t2496 + t2511) * t2339 - t2279 - t2304 * t2499;
t2314 = sin(t2340);
t2568 = t2250 * t2314;
t2317 = cos(t2340);
t2567 = t2250 * t2317;
t2327 = 0.1e1 / t2343;
t2566 = t2250 * t2327;
t2510 = t2305 * t2352;
t2252 = (-t2495 + t2510) * t2339 - t2280 - t2305 * t2498;
t2315 = sin(t2341);
t2565 = t2252 * t2315;
t2318 = cos(t2341);
t2564 = t2252 * t2318;
t2330 = 0.1e1 / t2345;
t2563 = t2252 * t2330;
t2509 = t2306 * t2354;
t2254 = (-t2494 + t2509) * t2339 - t2281 - t2306 * t2497;
t2316 = sin(t2342);
t2562 = t2254 * t2316;
t2319 = cos(t2342);
t2561 = t2254 * t2319;
t2333 = 0.1e1 / t2347;
t2560 = t2254 * t2333;
t2550 = t2282 * t2310;
t2549 = t2282 * t2314;
t2548 = t2282 * t2317;
t2547 = t2282 * t2327;
t2546 = t2282 * t2349;
t2544 = t2283 * t2314;
t2543 = t2283 * t2317;
t2326 = t2343 ^ 2;
t2542 = t2283 * t2326;
t2328 = 0.1e1 / t2343 ^ 2;
t2541 = t2283 * t2328;
t2540 = t2283 * t2349;
t2539 = t2284 * t2311;
t2538 = t2284 * t2315;
t2537 = t2284 * t2318;
t2536 = t2284 * t2330;
t2535 = t2284 * t2351;
t2533 = t2285 * t2315;
t2532 = t2285 * t2318;
t2329 = t2345 ^ 2;
t2531 = t2285 * t2329;
t2331 = 0.1e1 / t2345 ^ 2;
t2530 = t2285 * t2331;
t2529 = t2285 * t2351;
t2528 = t2286 * t2312;
t2527 = t2286 * t2316;
t2526 = t2286 * t2319;
t2525 = t2286 * t2333;
t2524 = t2286 * t2353;
t2522 = t2287 * t2316;
t2521 = t2287 * t2319;
t2332 = t2347 ^ 2;
t2520 = t2287 * t2332;
t2334 = 0.1e1 / t2347 ^ 2;
t2519 = t2287 * t2334;
t2518 = t2287 * t2353;
t2302 = pkin(1) * t2338 + pkin(5);
t2517 = t2302 * t2310;
t2516 = t2302 * t2311;
t2515 = t2302 * t2312;
t2514 = t2302 / 0.2e1;
t2513 = t2302 * t2358;
t2303 = t2313 + pkin(2);
t2512 = t2303 / 0.2e1;
t2508 = t2314 * t2343;
t2507 = t2315 * t2345;
t2506 = t2316 * t2347;
t2505 = t2317 * t2343;
t2504 = t2318 * t2345;
t2503 = t2319 * t2347;
t2502 = t2327 * t2349;
t2501 = t2330 * t2351;
t2500 = t2333 * t2353;
t2485 = 0.2e1 * t2358;
t2484 = -0.2e1 * t2512;
t2483 = 0.2e1 * t2512;
t2335 = t2349 ^ 2;
t2482 = pkin(3) * (-t2339 * t2350 + t2499) * t2335;
t2336 = t2351 ^ 2;
t2481 = pkin(3) * (-t2339 * t2352 + t2498) * t2336;
t2337 = t2353 ^ 2;
t2480 = pkin(3) * (-t2339 * t2354 + t2497) * t2337;
t2249 = (t2304 * t2344 + t2350 * t2357) * t2339 - t2294 * t2344 + t2338 * t2511;
t2243 = (t2249 + t2517) * t2546;
t2479 = MDP(11) * t2243 * t2250;
t2251 = (t2305 * t2346 + t2352 * t2357) * t2339 - t2294 * t2346 + t2338 * t2510;
t2244 = (t2251 + t2516) * t2535;
t2478 = MDP(11) * t2244 * t2252;
t2253 = (t2306 * t2348 + t2354 * t2357) * t2339 - t2294 * t2348 + t2338 * t2509;
t2245 = (t2253 + t2515) * t2524;
t2477 = MDP(11) * t2245 * t2254;
t2476 = MDP(11) * t2559;
t2475 = MDP(11) * t2558;
t2474 = MDP(11) * t2557;
t2473 = t2250 ^ 2 * t2541;
t2472 = t2252 ^ 2 * t2530;
t2471 = t2254 ^ 2 * t2519;
t2470 = t2249 * t2327 * t2335;
t2469 = t2250 * t2544;
t2468 = t2250 * t2543;
t2467 = t2250 * t2513;
t2466 = t2250 * t2502;
t2465 = t2251 * t2330 * t2336;
t2464 = t2252 * t2533;
t2463 = t2252 * t2532;
t2462 = t2252 * t2513;
t2461 = t2252 * t2501;
t2460 = t2253 * t2333 * t2337;
t2459 = t2254 * t2522;
t2458 = t2254 * t2521;
t2457 = t2254 * t2513;
t2456 = t2254 * t2500;
t2449 = t2283 * t2582;
t2448 = t2285 * t2581;
t2447 = t2287 * t2580;
t2446 = t2303 * t2550;
t2445 = t2310 * t2542;
t2444 = t2335 * t2541;
t2443 = t2328 * t2540;
t2442 = t2343 * t2540;
t2441 = t2303 * t2539;
t2440 = t2311 * t2531;
t2439 = t2336 * t2530;
t2438 = t2331 * t2529;
t2437 = t2345 * t2529;
t2436 = t2303 * t2528;
t2435 = t2312 * t2520;
t2434 = t2337 * t2519;
t2433 = t2334 * t2518;
t2432 = t2347 * t2518;
t2431 = t2343 * t2514;
t2430 = t2345 * t2514;
t2429 = t2347 * t2514;
t2428 = t2349 * t2514;
t2427 = t2351 * t2514;
t2426 = t2353 * t2514;
t2425 = t2486 / 0.2e1;
t2424 = t2486 / 0.4e1;
t2423 = t2485 / 0.2e1;
t2422 = t2343 * t2484;
t2421 = t2345 * t2484;
t2420 = t2347 * t2484;
t2419 = t2349 * t2483;
t2418 = t2351 * t2483;
t2417 = t2353 * t2483;
t2234 = t2314 * t2482 + (pkin(3) * t2505 + t2605 * t2314) * t2349 + t2303 * t2505;
t2416 = t2234 * t2443;
t2235 = t2315 * t2481 + (pkin(3) * t2504 + t2604 * t2315) * t2351 + t2303 * t2504;
t2415 = t2235 * t2438;
t2236 = t2316 * t2480 + (pkin(3) * t2503 + t2603 * t2316) * t2353 + t2303 * t2503;
t2414 = t2236 * t2433;
t2237 = -t2317 * t2482 + (pkin(3) * t2508 - t2605 * t2317) * t2349 + t2303 * t2508;
t2413 = t2237 * t2443;
t2238 = -t2318 * t2481 + (pkin(3) * t2507 - t2604 * t2318) * t2351 + t2303 * t2507;
t2412 = t2238 * t2438;
t2239 = -t2319 * t2480 + (pkin(3) * t2506 - t2603 * t2319) * t2353 + t2303 * t2506;
t2411 = t2239 * t2433;
t2410 = t2249 * t2443;
t2409 = t2283 * t2466;
t2408 = t2251 * t2438;
t2407 = t2285 * t2461;
t2406 = t2253 * t2433;
t2405 = t2287 * t2456;
t2402 = t2513 * t2559;
t2399 = t2513 * t2558;
t2396 = t2513 * t2557;
t2395 = t2310 * t2442;
t2394 = t2311 * t2437;
t2393 = t2312 * t2432;
t2392 = MDP(7) * t2423;
t2391 = MDP(8) * t2423;
t2390 = t2249 * t2250 * t2444;
t2389 = t2455 * t2566;
t2388 = t2314 * t2409;
t2387 = t2317 * t2409;
t2386 = t2466 * t2513;
t2385 = t2251 * t2252 * t2439;
t2384 = t2453 * t2563;
t2383 = t2315 * t2407;
t2382 = t2318 * t2407;
t2381 = t2461 * t2513;
t2380 = t2253 * t2254 * t2434;
t2379 = t2451 * t2560;
t2378 = t2316 * t2405;
t2377 = t2319 * t2405;
t2376 = t2456 * t2513;
t2360 = pkin(1) ^ 2;
t2363 = t2447 / 0.4e1 + t2448 / 0.4e1 + t2449 / 0.4e1;
t2367 = (t2272 * t2316 - t2275 * t2319) * t2287 * t2254;
t2368 = (t2271 * t2315 - t2274 * t2318) * t2285 * t2252;
t2369 = (t2270 * t2314 - t2273 * t2317) * t2283 * t2250;
t2372 = (t2234 * t2237 * t2541 + t2235 * t2238 * t2530 + t2236 * t2239 * t2519 + t2360 * t2363) * MDP(4) + (-t2314 * t2317 * t2473 - t2315 * t2318 * t2472 - t2316 * t2319 * t2471) * t2569 + (t2432 * t2580 + t2437 * t2581 + t2442 * t2582) * t2424 + (t2326 * t2449 + t2329 * t2448 + t2332 * t2447) * t2572 + MDP(1) * t2363 + ((t2367 * t2500 + t2368 * t2501 + t2369 * t2502) * MDP(8) + (t2367 + t2368 + t2369) * MDP(7)) * t2358 / 0.2e1;
t2362 = t2270 * t2589 + t2271 * t2590 + t2272 * t2591;
t2371 = (-t2237 * t2410 - t2238 * t2408 - t2239 * t2406 + t2360 * t2362) * MDP(4) + (t2270 * t2593 + t2271 * t2595 + t2272 * t2597 + t2310 * t2387 + t2311 * t2382 + t2312 * t2377) * t2570 + (t2270 * t2592 + t2271 * t2594 + t2272 * t2596 + t2310 * t2468 + t2311 * t2463 + t2312 * t2458) * t2571 + (-t2317 * t2389 - t2318 * t2384 - t2319 * t2379) * t2569 + (-t2270 * t2395 - t2271 * t2394 - t2272 * t2393) * t2425 + (-t2270 * t2445 - t2271 * t2440 - t2272 * t2435) * t2573 + MDP(1) * t2362;
t2361 = t2273 * t2589 + t2274 * t2590 + t2275 * t2591;
t2370 = (-t2234 * t2410 - t2235 * t2408 - t2236 * t2406 + t2360 * t2361) * MDP(4) + (t2273 * t2593 + t2274 * t2595 + t2275 * t2597 - t2310 * t2388 - t2311 * t2383 - t2312 * t2378) * t2570 + (t2273 * t2592 + t2274 * t2594 + t2275 * t2596 - t2310 * t2469 - t2311 * t2464 - t2312 * t2459) * t2571 + (t2314 * t2389 + t2315 * t2384 + t2316 * t2379) * t2569 + (-t2273 * t2395 - t2274 * t2394 - t2275 * t2393) * t2425 + (-t2273 * t2445 - t2274 * t2440 - t2275 * t2435) * t2573 + MDP(1) * t2361;
t2307 = t2310 ^ 2;
t2308 = t2311 ^ 2;
t2309 = t2312 ^ 2;
t2366 = t2283 * t2307 + t2285 * t2308 + t2287 * t2309;
t2264 = t2270 ^ 2;
t2265 = t2271 ^ 2;
t2266 = t2272 ^ 2;
t2365 = t2264 * t2600 + t2265 * t2599 + t2266 * t2598;
t2267 = t2273 ^ 2;
t2268 = t2274 ^ 2;
t2269 = t2275 ^ 2;
t2364 = t2267 * t2600 + t2268 * t2599 + t2269 * t2598;
t2242 = (t2347 * t2515 - t2460) * t2286;
t2241 = (t2345 * t2516 - t2465) * t2284;
t2240 = (t2343 * t2517 - t2470) * t2282;
t2233 = -t2343 * t2402 - 0.2e1 * t2349 * t2446;
t2232 = -t2347 * t2396 - 0.2e1 * t2353 * t2436;
t2231 = -t2345 * t2399 - 0.2e1 * t2351 * t2441;
t2230 = -t2353 * t2396 + t2436 * t2575;
t2229 = -t2351 * t2399 + t2441 * t2576;
t2228 = -t2349 * t2402 + t2446 * t2577;
t2227 = (t2273 * t2419 - t2314 * t2467) * t2282;
t2226 = (t2270 * t2419 + t2317 * t2467) * t2282;
t2225 = (t2275 * t2417 - t2316 * t2457) * t2286;
t2224 = (t2274 * t2418 - t2315 * t2462) * t2284;
t2223 = (t2272 * t2417 + t2319 * t2457) * t2286;
t2222 = (t2271 * t2418 + t2318 * t2462) * t2284;
t2219 = (t2272 * t2420 + t2319 * t2376) * t2286;
t2218 = (t2275 * t2420 - t2316 * t2376) * t2286;
t2217 = (t2271 * t2421 + t2318 * t2381) * t2284;
t2216 = (t2274 * t2421 - t2315 * t2381) * t2284;
t2215 = (t2270 * t2422 + t2317 * t2386) * t2282;
t2214 = (t2273 * t2422 - t2314 * t2386) * t2282;
t2209 = (-t2275 * t2426 - t2236) * t2286;
t2208 = (-t2274 * t2427 - t2235) * t2284;
t2207 = (-t2273 * t2428 - t2234) * t2282;
t2206 = (-t2272 * t2426 - t2239) * t2286;
t2205 = (-t2271 * t2427 - t2238) * t2284;
t2204 = (-t2270 * t2428 - t2237) * t2282;
t2202 = (t2236 * t2500 - t2275 * t2429) * t2286;
t2201 = (t2235 * t2501 - t2274 * t2430) * t2284;
t2200 = (t2234 * t2502 - t2273 * t2431) * t2282;
t2199 = (t2239 * t2500 - t2272 * t2429) * t2286;
t2198 = (t2238 * t2501 - t2271 * t2430) * t2284;
t2197 = (t2237 * t2502 - t2270 * t2431) * t2282;
t1 = [MDP(1) * t2365 + (t2237 ^ 2 * t2541 + t2238 ^ 2 * t2530 + t2239 ^ 2 * t2519 + t2360 * t2365) * MDP(4) + (t2264 * t2542 + t2265 * t2531 + t2266 * t2520) * t2572 + (t2264 * t2442 + t2265 * t2437 + t2266 * t2432) * t2424 + (-t2270 * t2468 - t2271 * t2463 - t2272 * t2458) * t2392 + (-t2270 * t2387 - t2271 * t2382 - t2272 * t2377) * t2391 + (t2317 ^ 2 * t2473 + t2318 ^ 2 * t2472 + t2319 ^ 2 * t2471) * t2569 + (t2222 * t2584 + t2223 * t2585 + t2226 * t2583 + ((-t2199 * t2525 - t2411) * t2561 + (-t2198 * t2536 - t2412) * t2564 + (-t2197 * t2547 - t2413) * t2567) * t2358) * MDP(10) + (t2215 * t2583 + t2217 * t2584 + t2219 * t2585 + ((-t2206 * t2286 + t2239 * t2287) * t2319 * t2560 + (-t2205 * t2284 + t2238 * t2285) * t2318 * t2563 + (-t2204 * t2282 + t2237 * t2283) * t2317 * t2566) * t2358) * MDP(11) + MDP(12); (t2224 * t2555 + t2225 * t2554 + t2227 * t2556) * t2602 + (t2214 * t2556 + t2216 * t2555 + t2218 * t2554) * t2601 + ((t2316 * MDP(10) * t2411 + (-t2202 * MDP(10) * t2526 + (-t2209 * t2526 - t2239 * t2522) * MDP(11)) * t2333) * t2254 + (t2315 * MDP(10) * t2412 + (-t2201 * MDP(10) * t2537 + (-t2208 * t2537 - t2238 * t2533) * MDP(11)) * t2330) * t2252 + (t2314 * MDP(10) * t2413 + (-t2200 * MDP(10) * t2548 + (-t2207 * t2548 - t2237 * t2544) * MDP(11)) * t2327) * t2250) * t2358 + t2372; (t2231 * t2555 + t2232 * t2554 + t2233 * t2556) * t2602 + (t2228 * t2556 + t2229 * t2555 + t2230 * t2554) * t2601 + ((-t2239 * t2474 + ((t2239 * t2450 - t2242 * t2561) * MDP(10) - t2319 * t2477) * t2333) * t2286 + (-t2238 * t2475 + ((t2238 * t2452 - t2241 * t2564) * MDP(10) - t2318 * t2478) * t2330) * t2284 + (-t2237 * t2476 + ((t2237 * t2454 - t2240 * t2567) * MDP(10) - t2317 * t2479) * t2327) * t2282) * t2358 + t2371; (t2222 * t2552 + t2223 * t2551 + t2226 * t2553) * t2602 + (t2215 * t2553 + t2217 * t2552 + t2219 * t2551) * t2601 + ((-t2319 * MDP(10) * t2414 + (t2199 * MDP(10) * t2527 + (t2206 * t2527 + t2236 * t2521) * MDP(11)) * t2333) * t2254 + (-t2318 * MDP(10) * t2415 + (t2198 * MDP(10) * t2538 + (t2205 * t2538 + t2235 * t2532) * MDP(11)) * t2330) * t2252 + (-t2317 * MDP(10) * t2416 + (t2197 * MDP(10) * t2549 + (t2204 * t2549 + t2234 * t2543) * MDP(11)) * t2327) * t2250) * t2358 + t2372; MDP(1) * t2364 + (t2234 ^ 2 * t2541 + t2235 ^ 2 * t2530 + t2236 ^ 2 * t2519 + t2360 * t2364) * MDP(4) + (t2267 * t2542 + t2268 * t2531 + t2269 * t2520) * t2572 + (t2267 * t2442 + t2268 * t2437 + t2269 * t2432) * t2424 + (t2273 * t2469 + t2274 * t2464 + t2275 * t2459) * t2392 + (t2273 * t2388 + t2274 * t2383 + t2275 * t2378) * t2391 + (t2314 ^ 2 * t2473 + t2315 ^ 2 * t2472 + t2316 ^ 2 * t2471) * t2569 + (t2224 * t2587 + t2225 * t2588 + t2227 * t2586 + ((t2202 * t2525 + t2414) * t2562 + (t2201 * t2536 + t2415) * t2565 + (t2200 * t2547 + t2416) * t2568) * t2358) * MDP(10) + (t2214 * t2586 + t2216 * t2587 + t2218 * t2588 + ((t2209 * t2286 - t2236 * t2287) * t2316 * t2560 + (t2208 * t2284 - t2235 * t2285) * t2315 * t2563 + (t2207 * t2282 - t2234 * t2283) * t2314 * t2566) * t2358) * MDP(11) + MDP(12); (t2231 * t2552 + t2232 * t2551 + t2233 * t2553) * t2602 + (t2228 * t2553 + t2229 * t2552 + t2230 * t2551) * t2601 + ((-t2236 * t2474 + ((t2236 * t2450 + t2242 * t2562) * MDP(10) + t2316 * t2477) * t2333) * t2286 + (-t2235 * t2475 + ((t2235 * t2452 + t2241 * t2565) * MDP(10) + t2315 * t2478) * t2330) * t2284 + (-t2234 * t2476 + ((t2234 * t2454 + t2240 * t2568) * MDP(10) + t2314 * t2479) * t2327) * t2282) * t2358 + t2370; (-t2222 * t2539 - t2223 * t2528 - t2226 * t2550) * MDP(10) + (-t2215 * t2550 - t2217 * t2539 - t2219 * t2528) * MDP(11) + ((t2197 * t2559 + t2198 * t2558 + t2199 * t2557 + t2317 * t2390 + t2318 * t2385 + t2319 * t2380) * MDP(10) + (t2204 * t2559 + t2205 * t2558 + t2206 * t2557 - t2249 * t2387 - t2251 * t2382 - t2253 * t2377) * MDP(11)) * t2358 + t2371; (-t2224 * t2539 - t2225 * t2528 - t2227 * t2550) * MDP(10) + (-t2214 * t2550 - t2216 * t2539 - t2218 * t2528) * MDP(11) + ((t2200 * t2559 + t2201 * t2558 + t2202 * t2557 - t2314 * t2390 - t2315 * t2385 - t2316 * t2380) * MDP(10) + (t2207 * t2559 + t2208 * t2558 + t2209 * t2557 + t2249 * t2388 + t2251 * t2383 + t2253 * t2378) * MDP(11)) * t2358 + t2370; t2366 * MDP(1) + (t2249 ^ 2 * t2444 + t2251 ^ 2 * t2439 + t2253 ^ 2 * t2434 + t2360 * t2366) * MDP(4) + (t2307 * t2542 + t2308 * t2531 + t2309 * t2520) * MDP(5) + (t2307 * t2442 + t2308 * t2437 + t2309 * t2432) * t2486 + (t2257 ^ 2 / t2263 ^ 2 + t2256 ^ 2 / t2262 ^ 2 + t2255 ^ 2 / t2261 ^ 2) * t2569 + (-t2231 * t2539 - t2232 * t2528 - t2233 * t2550 + ((-t2286 * t2460 + t2242) * t2557 + (-t2284 * t2465 + t2241) * t2558 + (-t2282 * t2470 + t2240) * t2559) * t2358) * MDP(10) + (-t2228 * t2550 - t2229 * t2539 - t2230 * t2528 + ((t2253 * t2524 + t2245) * t2557 + (t2251 * t2535 + t2244) * t2558 + (t2249 * t2546 + t2243) * t2559) * t2358) * MDP(11) + MDP(12) + ((-t2310 * t2404 - t2311 * t2401 - t2312 * t2398) * MDP(7) + (-t2310 * t2403 - t2311 * t2400 - t2312 * t2397) * MDP(8)) * t2485;];
%% Postprocessing: Reshape Output
% From vec2mat_3_matlab.m
res = [t1(1), t1(2), t1(3); t1(4), t1(5), t1(6); t1(7), t1(8), t1(9);];
MMX  = res;
