% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P4PRRRR8V2G1A0
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
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [4x15]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:13
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P4PRRRR8V2G1A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_regmin: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_regmin: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_regmin: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_regmin: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_regmin: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:11:43
% EndTime: 2020-08-07 11:11:48
% DurationCPUTime: 4.64s
% Computational Cost: add. (2194->409), mult. (4645->792), div. (112->17), fcn. (4976->30), ass. (0->300)
t2420 = cos(qJ(2,1));
t2414 = sin(qJ(2,1));
t2421 = pkin(7) + pkin(6);
t2467 = t2414 * t2421;
t2368 = pkin(2) * t2420 + t2467;
t2397 = sin(pkin(8));
t2399 = cos(pkin(8));
t2378 = t2421 * t2420;
t2365 = pkin(2) * t2414 - t2378;
t2400 = cos(pkin(4));
t2398 = sin(pkin(4));
t2413 = sin(qJ(3,1));
t2501 = t2398 * t2413;
t2432 = pkin(3) * t2501 - t2365 * t2400;
t2553 = t2368 * t2399 + t2432 * t2397;
t2418 = cos(qJ(2,2));
t2412 = sin(qJ(2,2));
t2471 = t2412 * t2421;
t2367 = pkin(2) * t2418 + t2471;
t2377 = t2421 * t2418;
t2364 = pkin(2) * t2412 - t2377;
t2411 = sin(qJ(3,2));
t2503 = t2398 * t2411;
t2433 = pkin(3) * t2503 - t2364 * t2400;
t2552 = t2367 * t2399 + t2433 * t2397;
t2416 = cos(qJ(2,3));
t2410 = sin(qJ(2,3));
t2475 = t2410 * t2421;
t2366 = pkin(2) * t2416 + t2475;
t2376 = t2421 * t2416;
t2363 = pkin(2) * t2410 - t2376;
t2409 = sin(qJ(3,3));
t2505 = t2398 * t2409;
t2434 = pkin(3) * t2505 - t2363 * t2400;
t2551 = t2366 * t2399 + t2434 * t2397;
t2408 = cos(qJ(2,4));
t2406 = sin(qJ(2,4));
t2480 = t2406 * t2421;
t2362 = pkin(2) * t2408 + t2480;
t2372 = t2421 * t2408;
t2361 = pkin(2) * t2406 - t2372;
t2405 = sin(qJ(3,4));
t2508 = t2398 * t2405;
t2435 = pkin(3) * t2508 - t2361 * t2400;
t2550 = t2362 * t2399 + t2435 * t2397;
t2407 = cos(qJ(3,4));
t2549 = pkin(3) * t2407 ^ 2;
t2415 = cos(qJ(3,3));
t2548 = pkin(3) * t2415 ^ 2;
t2417 = cos(qJ(3,2));
t2547 = pkin(3) * t2417 ^ 2;
t2419 = cos(qJ(3,1));
t2546 = pkin(3) * t2419 ^ 2;
t2545 = g(3) * t2398;
t2423 = koppelP(4,2);
t2427 = koppelP(4,1);
t2351 = t2397 * t2423 + t2399 * t2427;
t2352 = -t2397 * t2427 + t2399 * t2423;
t2422 = xP(4);
t2391 = sin(t2422);
t2392 = cos(t2422);
t2291 = t2351 * t2391 + t2352 * t2392;
t2401 = legFrame(4,3);
t2379 = sin(t2401);
t2383 = cos(t2401);
t2439 = t2351 * t2392 - t2352 * t2391;
t2235 = t2291 * t2383 - t2439 * t2379;
t2495 = t2400 * t2406;
t2331 = t2408 * t2423 - t2427 * t2495;
t2332 = t2408 * t2427 + t2423 * t2495;
t2275 = t2331 * t2399 - t2332 * t2397;
t2276 = t2331 * t2397 + t2332 * t2399;
t2371 = pkin(3) * t2407 + pkin(2);
t2496 = t2400 * t2405;
t2347 = t2371 * t2496;
t2481 = t2406 * t2407;
t2506 = t2398 * t2407;
t2544 = (((-t2275 * t2391 + t2276 * t2392) * t2383 + (t2275 * t2392 + t2276 * t2391) * t2379) * t2405 + t2235 * t2506) / (t2347 + (pkin(3) * t2481 + t2361) * t2506);
t2424 = koppelP(3,2);
t2428 = koppelP(3,1);
t2353 = t2397 * t2424 + t2399 * t2428;
t2354 = -t2397 * t2428 + t2399 * t2424;
t2292 = t2353 * t2391 + t2354 * t2392;
t2402 = legFrame(3,3);
t2380 = sin(t2402);
t2384 = cos(t2402);
t2438 = t2353 * t2392 - t2354 * t2391;
t2236 = t2292 * t2384 - t2438 * t2380;
t2492 = t2400 * t2410;
t2333 = t2416 * t2424 - t2428 * t2492;
t2336 = t2416 * t2428 + t2424 * t2492;
t2277 = t2333 * t2399 - t2336 * t2397;
t2280 = t2333 * t2397 + t2336 * t2399;
t2373 = pkin(3) * t2415 + pkin(2);
t2493 = t2400 * t2409;
t2348 = t2373 * t2493;
t2476 = t2410 * t2415;
t2499 = t2398 * t2415;
t2543 = (((-t2277 * t2391 + t2280 * t2392) * t2384 + (t2277 * t2392 + t2280 * t2391) * t2380) * t2409 + t2236 * t2499) / (t2348 + (pkin(3) * t2476 + t2363) * t2499);
t2425 = koppelP(2,2);
t2429 = koppelP(2,1);
t2355 = t2397 * t2425 + t2399 * t2429;
t2356 = -t2397 * t2429 + t2399 * t2425;
t2293 = t2355 * t2391 + t2356 * t2392;
t2403 = legFrame(2,3);
t2381 = sin(t2403);
t2385 = cos(t2403);
t2437 = t2355 * t2392 - t2356 * t2391;
t2237 = t2293 * t2385 - t2437 * t2381;
t2490 = t2400 * t2412;
t2334 = t2418 * t2425 - t2429 * t2490;
t2337 = t2418 * t2429 + t2425 * t2490;
t2278 = t2334 * t2399 - t2337 * t2397;
t2281 = t2334 * t2397 + t2337 * t2399;
t2374 = pkin(3) * t2417 + pkin(2);
t2491 = t2400 * t2411;
t2349 = t2374 * t2491;
t2472 = t2412 * t2417;
t2498 = t2398 * t2417;
t2542 = (((-t2278 * t2391 + t2281 * t2392) * t2385 + t2381 * (t2278 * t2392 + t2281 * t2391)) * t2411 + t2237 * t2498) / (t2349 + (pkin(3) * t2472 + t2364) * t2498);
t2426 = koppelP(1,2);
t2430 = koppelP(1,1);
t2357 = t2397 * t2426 + t2399 * t2430;
t2358 = -t2397 * t2430 + t2399 * t2426;
t2294 = t2357 * t2391 + t2358 * t2392;
t2404 = legFrame(1,3);
t2382 = sin(t2404);
t2386 = cos(t2404);
t2436 = t2357 * t2392 - t2358 * t2391;
t2238 = t2294 * t2386 - t2382 * t2436;
t2488 = t2400 * t2414;
t2335 = t2420 * t2426 - t2430 * t2488;
t2338 = t2420 * t2430 + t2426 * t2488;
t2279 = t2335 * t2399 - t2338 * t2397;
t2282 = t2335 * t2397 + t2338 * t2399;
t2375 = pkin(3) * t2419 + pkin(2);
t2489 = t2400 * t2413;
t2350 = t2375 * t2489;
t2468 = t2414 * t2419;
t2497 = t2398 * t2419;
t2541 = (((-t2279 * t2391 + t2282 * t2392) * t2386 + (t2279 * t2392 + t2282 * t2391) * t2382) * t2413 + t2238 * t2497) / (t2350 + (pkin(3) * t2468 + t2365) * t2497);
t2315 = t2371 * t2406 - t2372;
t2283 = 0.1e1 / (t2315 * t2506 + t2347);
t2516 = (t2371 * t2408 + t2480) * t2400;
t2540 = (t2235 * t2516 - t2315 * (t2291 * t2379 + t2439 * t2383)) * t2283;
t2328 = t2373 * t2410 - t2376;
t2285 = 0.1e1 / (t2328 * t2499 + t2348);
t2515 = (t2373 * t2416 + t2475) * t2400;
t2539 = (t2236 * t2515 - t2328 * (t2292 * t2380 + t2438 * t2384)) * t2285;
t2329 = t2374 * t2412 - t2377;
t2286 = 0.1e1 / (t2329 * t2498 + t2349);
t2514 = (t2374 * t2418 + t2471) * t2400;
t2538 = (t2237 * t2514 - (t2293 * t2381 + t2437 * t2385) * t2329) * t2286;
t2330 = t2375 * t2414 - t2378;
t2287 = 0.1e1 / (t2330 * t2497 + t2350);
t2513 = (t2375 * t2420 + t2467) * t2400;
t2537 = (t2238 * t2513 - (t2382 * t2294 + t2436 * t2386) * t2330) * t2287;
t2369 = g(1) * t2397 - g(2) * t2399;
t2370 = g(1) * t2399 + g(2) * t2397;
t2231 = (-t2545 + (t2369 * t2383 + t2370 * t2379) * t2400) * t2408 + (-t2369 * t2379 + t2370 * t2383) * t2406;
t2507 = t2398 * t2406;
t2263 = 0.1e1 / (t2507 * t2549 + (pkin(3) * t2496 + t2361 * t2398) * t2407 + pkin(2) * t2496);
t2536 = t2231 * t2263;
t2232 = (-t2545 + (t2369 * t2384 + t2370 * t2380) * t2400) * t2416 + (-t2369 * t2380 + t2370 * t2384) * t2410;
t2504 = t2398 * t2410;
t2264 = 0.1e1 / (t2504 * t2548 + (pkin(3) * t2493 + t2363 * t2398) * t2415 + pkin(2) * t2493);
t2535 = t2232 * t2264;
t2233 = (-t2545 + (t2369 * t2385 + t2370 * t2381) * t2400) * t2418 + (-t2369 * t2381 + t2370 * t2385) * t2412;
t2502 = t2398 * t2412;
t2265 = 0.1e1 / (t2502 * t2547 + (pkin(3) * t2491 + t2364 * t2398) * t2417 + pkin(2) * t2491);
t2534 = t2233 * t2265;
t2234 = (-t2545 + (t2369 * t2386 + t2370 * t2382) * t2400) * t2420 + (-t2369 * t2382 + t2370 * t2386) * t2414;
t2500 = t2398 * t2414;
t2266 = 0.1e1 / (t2500 * t2546 + (pkin(3) * t2489 + t2365 * t2398) * t2419 + pkin(2) * t2489);
t2533 = t2234 * t2266;
t2295 = -t2379 * t2397 + t2383 * t2399;
t2299 = t2379 * t2399 + t2383 * t2397;
t2532 = (-t2295 * t2516 + t2299 * t2315) * t2283;
t2531 = (-t2295 * t2315 - t2299 * t2516) * t2283;
t2296 = -t2380 * t2397 + t2384 * t2399;
t2300 = t2380 * t2399 + t2384 * t2397;
t2530 = (-t2296 * t2515 + t2300 * t2328) * t2285;
t2297 = -t2381 * t2397 + t2385 * t2399;
t2301 = t2381 * t2399 + t2385 * t2397;
t2529 = (-t2297 * t2514 + t2301 * t2329) * t2286;
t2298 = -t2382 * t2397 + t2386 * t2399;
t2302 = t2382 * t2399 + t2386 * t2397;
t2528 = (-t2298 * t2513 + t2302 * t2330) * t2287;
t2527 = (-t2296 * t2328 - t2300 * t2515) * t2285;
t2526 = (-t2297 * t2329 - t2301 * t2514) * t2286;
t2525 = (-t2298 * t2330 - t2302 * t2513) * t2287;
t2339 = g(1) * t2379 - g(2) * t2383;
t2343 = g(1) * t2383 + g(2) * t2379;
t2494 = t2400 * t2408;
t2255 = t2343 * (t2397 * t2494 + t2399 * t2406) + t2339 * (-t2397 * t2406 + t2399 * t2494) - t2408 * t2545;
t2524 = t2255 * t2263;
t2303 = t2397 * t2495 - t2399 * t2408;
t2304 = t2397 * t2408 + t2399 * t2495;
t2256 = g(3) * t2507 - t2303 * t2343 - t2304 * t2339;
t2523 = t2256 * t2263;
t2340 = g(1) * t2380 - g(2) * t2384;
t2344 = g(1) * t2384 + g(2) * t2380;
t2487 = t2400 * t2416;
t2257 = t2344 * (t2397 * t2487 + t2399 * t2410) + t2340 * (-t2397 * t2410 + t2399 * t2487) - t2416 * t2545;
t2522 = t2257 * t2264;
t2308 = t2397 * t2492 - t2399 * t2416;
t2311 = t2397 * t2416 + t2399 * t2492;
t2258 = g(3) * t2504 - t2308 * t2344 - t2311 * t2340;
t2521 = t2258 * t2264;
t2341 = g(1) * t2381 - g(2) * t2385;
t2345 = g(1) * t2385 + g(2) * t2381;
t2486 = t2400 * t2418;
t2259 = t2345 * (t2397 * t2486 + t2399 * t2412) + t2341 * (-t2397 * t2412 + t2399 * t2486) - t2418 * t2545;
t2520 = t2259 * t2265;
t2309 = t2397 * t2490 - t2399 * t2418;
t2312 = t2397 * t2418 + t2399 * t2490;
t2260 = g(3) * t2502 - t2309 * t2345 - t2312 * t2341;
t2519 = t2260 * t2265;
t2342 = g(1) * t2382 - g(2) * t2386;
t2346 = g(1) * t2386 + g(2) * t2382;
t2485 = t2400 * t2420;
t2261 = t2346 * (t2397 * t2485 + t2399 * t2414) + t2342 * (-t2397 * t2414 + t2399 * t2485) - t2420 * t2545;
t2518 = t2261 * t2266;
t2310 = t2397 * t2488 - t2399 * t2420;
t2313 = t2397 * t2420 + t2399 * t2488;
t2262 = g(3) * t2500 - t2310 * t2346 - t2313 * t2342;
t2517 = t2262 * t2266;
t2484 = t2400 * t2421;
t2483 = t2405 * t2406;
t2482 = t2405 * t2408;
t2479 = t2407 * t2408;
t2478 = t2409 * t2410;
t2477 = t2409 * t2416;
t2474 = t2411 * t2412;
t2473 = t2411 * t2418;
t2470 = t2413 * t2414;
t2469 = t2413 * t2420;
t2466 = t2415 * t2416;
t2465 = t2417 * t2418;
t2464 = t2419 * t2420;
t2463 = pkin(2) * t2508;
t2462 = pkin(2) * t2505;
t2461 = pkin(2) * t2503;
t2460 = pkin(2) * t2501;
t2459 = t2231 * t2544;
t2458 = t2232 * t2543;
t2457 = t2233 * t2542;
t2456 = t2234 * t2541;
t2455 = t2405 * t2536;
t2454 = t2407 * t2536;
t2453 = t2409 * t2535;
t2452 = t2415 * t2535;
t2451 = t2411 * t2534;
t2450 = t2417 * t2534;
t2449 = t2413 * t2533;
t2448 = t2419 * t2533;
t2447 = t2371 * t2508;
t2446 = t2373 * t2505;
t2445 = t2374 * t2503;
t2444 = t2375 * t2501;
t2443 = t2371 * t2400;
t2442 = t2373 * t2400;
t2441 = t2374 * t2400;
t2440 = t2375 * t2400;
t2431 = 0.1e1 / pkin(3);
t2360 = g(1) * t2392 + g(2) * t2391;
t2359 = g(1) * t2391 - g(2) * t2392;
t2324 = t2400 * t2468 - t2501;
t2323 = t2400 * t2472 - t2503;
t2322 = t2400 * t2476 - t2505;
t2321 = t2398 * t2468 + t2489;
t2320 = t2400 * t2470 + t2497;
t2319 = t2398 * t2472 + t2491;
t2318 = t2400 * t2474 + t2498;
t2317 = t2398 * t2476 + t2493;
t2316 = t2400 * t2478 + t2499;
t2307 = t2400 * t2481 - t2508;
t2306 = t2398 * t2481 + t2496;
t2305 = t2400 * t2483 + t2506;
t2270 = t2368 * t2397 - t2432 * t2399;
t2269 = t2367 * t2397 - t2433 * t2399;
t2268 = t2366 * t2397 - t2434 * t2399;
t2267 = t2362 * t2397 - t2435 * t2399;
t2246 = -t2302 * t2497 - (-t2298 * t2420 + t2302 * t2488) * t2413;
t2245 = -t2301 * t2498 - (-t2297 * t2418 + t2301 * t2490) * t2411;
t2244 = -t2300 * t2499 - (-t2296 * t2416 + t2300 * t2492) * t2409;
t2243 = -t2298 * t2497 - (t2298 * t2488 + t2302 * t2420) * t2413;
t2242 = -t2297 * t2498 - (t2297 * t2490 + t2301 * t2418) * t2411;
t2241 = -t2296 * t2499 - (t2296 * t2492 + t2300 * t2416) * t2409;
t2240 = -t2299 * t2506 - (-t2295 * t2408 + t2299 * t2495) * t2405;
t2239 = -t2295 * t2506 - (t2295 * t2495 + t2299 * t2408) * t2405;
t2230 = (-t2324 * t2397 + t2399 * t2464) * t2346 - t2342 * (t2324 * t2399 + t2397 * t2464) + g(3) * t2321;
t2229 = (-t2323 * t2397 + t2399 * t2465) * t2345 - (t2323 * t2399 + t2397 * t2465) * t2341 + g(3) * t2319;
t2228 = t2346 * (-t2320 * t2397 + t2399 * t2469) - t2342 * (t2320 * t2399 + t2397 * t2469) - g(3) * (-t2398 * t2470 + t2400 * t2419);
t2227 = t2345 * (-t2318 * t2397 + t2399 * t2473) - (t2318 * t2399 + t2397 * t2473) * t2341 - g(3) * (-t2398 * t2474 + t2400 * t2417);
t2226 = t2344 * (-t2322 * t2397 + t2399 * t2466) - (t2322 * t2399 + t2397 * t2466) * t2340 + g(3) * t2317;
t2225 = t2344 * (-t2316 * t2397 + t2399 * t2477) - (t2316 * t2399 + t2397 * t2477) * t2340 - g(3) * (-t2398 * t2478 + t2400 * t2415);
t2224 = t2343 * (-t2307 * t2397 + t2399 * t2479) - (t2307 * t2399 + t2397 * t2479) * t2339 + g(3) * t2306;
t2223 = t2343 * (-t2305 * t2397 + t2399 * t2482) - (t2305 * t2399 + t2397 * t2482) * t2339 - g(3) * (-t2398 * t2483 + t2400 * t2407);
t1 = [(-(-(t2310 * t2386 + t2313 * t2382) * t2546 + (-t2382 * t2270 + t2553 * t2386) * t2419 + t2302 * t2460) * t2266 - (-(t2309 * t2385 + t2312 * t2381) * t2547 + (-t2381 * t2269 + t2552 * t2385) * t2417 + t2301 * t2461) * t2265 - (-(t2308 * t2384 + t2311 * t2380) * t2548 + (-t2380 * t2268 + t2551 * t2384) * t2415 + t2300 * t2462) * t2264 - (-(t2303 * t2383 + t2304 * t2379) * t2549 + (-t2379 * t2267 + t2550 * t2383) * t2407 + t2299 * t2463) * t2263) * g(3), 0, t2239 * t2524 + t2241 * t2522 + t2242 * t2520 + t2243 * t2518, t2239 * t2523 + t2241 * t2521 + t2242 * t2519 + t2243 * t2517, 0, 0, 0, 0, 0, t2239 * t2454 + t2241 * t2452 + t2242 * t2450 + t2243 * t2448 + (t2223 * t2532 + t2225 * t2530 + t2227 * t2529 + t2228 * t2528) * t2431, -t2239 * t2455 - t2241 * t2453 - t2242 * t2451 - t2243 * t2449 + (t2224 * t2532 + t2226 * t2530 + t2229 * t2529 + t2230 * t2528) * t2431, 0, 0, 0, -t2359 * t2391 - t2360 * t2392; (-((-t2310 * t2382 + t2313 * t2386) * t2546 + (t2270 * t2386 + t2553 * t2382) * t2419 - t2298 * t2460) * t2266 - ((-t2309 * t2381 + t2312 * t2385) * t2547 + (t2269 * t2385 + t2552 * t2381) * t2417 - t2297 * t2461) * t2265 - ((-t2308 * t2380 + t2311 * t2384) * t2548 + (t2268 * t2384 + t2551 * t2380) * t2415 - t2296 * t2462) * t2264 - ((-t2303 * t2379 + t2304 * t2383) * t2549 + (t2267 * t2383 + t2550 * t2379) * t2407 - t2295 * t2463) * t2263) * g(3), 0, t2240 * t2524 + t2244 * t2522 + t2245 * t2520 + t2246 * t2518, t2240 * t2523 + t2244 * t2521 + t2245 * t2519 + t2246 * t2517, 0, 0, 0, 0, 0, t2240 * t2454 + t2244 * t2452 + t2245 * t2450 + t2246 * t2448 + (t2223 * t2531 + t2225 * t2527 + t2227 * t2526 + t2228 * t2525) * t2431, -t2240 * t2455 - t2244 * t2453 - t2245 * t2451 - t2246 * t2449 + (t2224 * t2531 + t2226 * t2527 + t2229 * t2526 + t2230 * t2525) * t2431, 0, 0, 0, t2359 * t2392 - t2360 * t2391; -0.4e1 * g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); (-(-(t2391 * t2430 + t2392 * t2426) * ((t2298 * t2375 + t2302 * t2484) * t2464 - (-t2298 * t2421 + t2302 * t2440) * t2468 + t2302 * t2444) + (-t2391 * t2426 + t2392 * t2430) * ((-t2298 * t2484 + t2302 * t2375) * t2464 + (t2298 * t2440 + t2302 * t2421) * t2468 - t2298 * t2444)) / (t2321 * t2375 - t2378 * t2497) - (-(t2391 * t2429 + t2392 * t2425) * ((t2297 * t2374 + t2301 * t2484) * t2465 - (-t2297 * t2421 + t2301 * t2441) * t2472 + t2301 * t2445) + (-t2391 * t2425 + t2392 * t2429) * ((-t2297 * t2484 + t2301 * t2374) * t2465 + (t2297 * t2441 + t2301 * t2421) * t2472 - t2297 * t2445)) / (t2319 * t2374 - t2377 * t2498) - (-(t2391 * t2428 + t2392 * t2424) * ((t2296 * t2373 + t2300 * t2484) * t2466 - (-t2296 * t2421 + t2300 * t2442) * t2476 + t2300 * t2446) + (-t2391 * t2424 + t2392 * t2428) * ((-t2296 * t2484 + t2300 * t2373) * t2466 + (t2296 * t2442 + t2300 * t2421) * t2476 - t2296 * t2446)) / (t2317 * t2373 - t2376 * t2499) - (-(t2391 * t2427 + t2392 * t2423) * ((t2295 * t2371 + t2299 * t2484) * t2479 - (-t2295 * t2421 + t2299 * t2443) * t2481 + t2299 * t2447) + (-t2391 * t2423 + t2392 * t2427) * ((-t2295 * t2484 + t2299 * t2371) * t2479 + (t2295 * t2443 + t2299 * t2421) * t2481 - t2295 * t2447)) / (t2306 * t2371 - t2372 * t2506)) * g(3), 0, t2255 * t2544 + t2257 * t2543 + t2259 * t2542 + t2261 * t2541, t2256 * t2544 + t2258 * t2543 + t2260 * t2542 + t2262 * t2541, 0, 0, 0, 0, 0, t2407 * t2459 + t2415 * t2458 + t2417 * t2457 + t2419 * t2456 + (t2223 * t2540 + t2225 * t2539 + t2227 * t2538 + t2228 * t2537) * t2431, -t2405 * t2459 - t2409 * t2458 - t2411 * t2457 - t2413 * t2456 + (t2224 * t2540 + t2226 * t2539 + t2229 * t2538 + t2230 * t2537) * t2431, 0, t2359, t2360, 0;];
tau_reg  = t1;