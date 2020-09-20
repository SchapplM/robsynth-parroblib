% Calculate Gravitation load for parallel robot
% P3RRRRR7V2G3A0
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
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 10:47
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRRRR7V2G3A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 10:23:37
% EndTime: 2020-08-07 10:23:39
% DurationCPUTime: 2.17s
% Computational Cost: add. (1173->320), mult. (1818->432), div. (66->14), fcn. (1203->66), ass. (0->208)
t2447 = sin(qJ(3,1));
t2456 = cos(qJ(3,1));
t2568 = -m(3) * pkin(2) - mrSges(2,1);
t2357 = -mrSges(3,1) * t2456 + mrSges(3,2) * t2447 + t2568;
t2373 = t2447 * mrSges(3,1) + t2456 * mrSges(3,2) + mrSges(2,2);
t2448 = sin(qJ(2,1));
t2457 = cos(qJ(2,1));
t2585 = -t2457 * t2357 - t2373 * t2448;
t2444 = sin(qJ(3,2));
t2453 = cos(qJ(3,2));
t2356 = -mrSges(3,1) * t2453 + mrSges(3,2) * t2444 + t2568;
t2372 = t2444 * mrSges(3,1) + t2453 * mrSges(3,2) + mrSges(2,2);
t2445 = sin(qJ(2,2));
t2454 = cos(qJ(2,2));
t2584 = -t2454 * t2356 - t2372 * t2445;
t2441 = sin(qJ(3,3));
t2450 = cos(qJ(3,3));
t2355 = -mrSges(3,1) * t2450 + mrSges(3,2) * t2441 + t2568;
t2371 = t2441 * mrSges(3,1) + t2450 * mrSges(3,2) + mrSges(2,2);
t2442 = sin(qJ(2,3));
t2451 = cos(qJ(2,3));
t2583 = -t2451 * t2355 - t2371 * t2442;
t2582 = 2 * pkin(3);
t2460 = pkin(6) + pkin(5);
t2423 = (pkin(7) + t2460);
t2581 = 2 * t2423;
t2580 = 0.2e1 * t2451 ^ 2;
t2579 = 0.2e1 * t2454 ^ 2;
t2578 = 0.2e1 * t2457 ^ 2;
t2427 = t2450 ^ 2;
t2416 = pkin(3) * t2427;
t2429 = t2453 ^ 2;
t2417 = pkin(3) * t2429;
t2431 = t2456 ^ 2;
t2418 = pkin(3) * t2431;
t2577 = t2441 * pkin(1);
t2576 = t2441 * pkin(3);
t2575 = t2444 * pkin(1);
t2574 = t2444 * pkin(3);
t2573 = t2447 * pkin(1);
t2572 = t2447 * pkin(3);
t2571 = t2450 * pkin(2);
t2413 = t2450 * pkin(3);
t2570 = t2453 * pkin(2);
t2414 = t2453 * pkin(3);
t2569 = t2456 * pkin(2);
t2415 = t2456 * pkin(3);
t2567 = -qJ(3,1) + qJ(1,1);
t2566 = qJ(3,1) + qJ(1,1);
t2565 = -qJ(3,2) + qJ(1,2);
t2564 = qJ(3,2) + qJ(1,2);
t2563 = -qJ(3,3) + qJ(1,3);
t2562 = qJ(3,3) + qJ(1,3);
t2561 = qJ(1,1) + 0.2e1 * qJ(2,1);
t2560 = qJ(1,1) - 0.2e1 * qJ(2,1);
t2559 = qJ(1,2) + 0.2e1 * qJ(2,2);
t2558 = qJ(1,2) - 0.2e1 * qJ(2,2);
t2557 = qJ(1,3) + 0.2e1 * qJ(2,3);
t2556 = -0.2e1 * qJ(2,3) + qJ(1,3);
t2555 = 0.2e1 * pkin(1);
t2452 = cos(qJ(1,3));
t2443 = sin(qJ(1,3));
t2508 = pkin(1) * t2452 + t2443 * t2423;
t2518 = t2441 * t2442;
t2337 = t2508 * t2518 + (t2427 - 0.1e1) * t2452 * pkin(3);
t2438 = legFrame(3,2);
t2407 = sin(t2438);
t2554 = t2337 * t2407;
t2410 = cos(t2438);
t2553 = t2337 * t2410;
t2455 = cos(qJ(1,2));
t2446 = sin(qJ(1,2));
t2507 = pkin(1) * t2455 + t2446 * t2423;
t2517 = t2444 * t2445;
t2338 = t2507 * t2517 + (t2429 - 0.1e1) * t2455 * pkin(3);
t2439 = legFrame(2,2);
t2408 = sin(t2439);
t2552 = t2338 * t2408;
t2411 = cos(t2439);
t2551 = t2338 * t2411;
t2458 = cos(qJ(1,1));
t2449 = sin(qJ(1,1));
t2506 = pkin(1) * t2458 + t2449 * t2423;
t2516 = t2447 * t2448;
t2339 = t2506 * t2516 + (t2431 - 0.1e1) * t2458 * pkin(3);
t2440 = legFrame(1,2);
t2409 = sin(t2440);
t2550 = t2339 * t2409;
t2412 = cos(t2440);
t2549 = t2339 * t2412;
t2498 = pkin(3) * t2518;
t2389 = t2413 + pkin(2);
t2535 = t2389 * t2451;
t2548 = 0.1e1 / (pkin(1) - t2498 + t2535) / t2441;
t2497 = pkin(3) * t2517;
t2390 = t2414 + pkin(2);
t2533 = t2390 * t2454;
t2547 = 0.1e1 / (pkin(1) - t2497 + t2533) / t2444;
t2496 = pkin(3) * t2516;
t2391 = t2415 + pkin(2);
t2531 = t2391 * t2457;
t2546 = 0.1e1 / (pkin(1) - t2496 + t2531) / t2447;
t2364 = t2410 * g(1) - t2407 * g(2);
t2480 = m(2) * pkin(5) + t2460 * m(3) - mrSges(1,2) + mrSges(2,3) + mrSges(3,3);
t2367 = t2480 * g(3);
t2382 = mrSges(1,1) + (m(2) + m(3)) * pkin(1);
t2381 = g(3) * t2382;
t2433 = qJ(2,3) + qJ(3,3);
t2395 = cos(t2433);
t2545 = 0.1e1 / (t2451 * pkin(2) + pkin(3) * t2395 + pkin(1)) * ((t2583 * g(3) + t2381) * t2452 + t2367 * t2443 + (-t2480 * t2452 + (t2382 + t2583) * t2443) * t2364);
t2365 = t2411 * g(1) - t2408 * g(2);
t2434 = qJ(2,2) + qJ(3,2);
t2396 = cos(t2434);
t2544 = 0.1e1 / (t2454 * pkin(2) + pkin(3) * t2396 + pkin(1)) * ((t2584 * g(3) + t2381) * t2455 + t2367 * t2446 + (-t2480 * t2455 + (t2382 + t2584) * t2446) * t2365);
t2366 = t2412 * g(1) - t2409 * g(2);
t2435 = qJ(2,1) + qJ(3,1);
t2397 = cos(t2435);
t2543 = 0.1e1 / (t2457 * pkin(2) + pkin(3) * t2397 + pkin(1)) * ((t2585 * g(3) + t2381) * t2458 + t2367 * t2449 + (-t2480 * t2458 + (t2382 + t2585) * t2449) * t2366);
t2463 = pkin(2) / 0.2e1;
t2539 = (t2413 + t2463) * t2441;
t2538 = (t2414 + t2463) * t2444;
t2537 = (t2415 + t2463) * t2447;
t2536 = t2389 * t2407;
t2534 = t2390 * t2408;
t2532 = t2391 * t2409;
t2530 = t2407 * t2452;
t2529 = t2408 * t2455;
t2528 = t2409 * t2458;
t2527 = t2410 * t2389;
t2526 = t2410 * t2452;
t2525 = t2411 * t2390;
t2524 = t2411 * t2455;
t2523 = t2412 * t2391;
t2522 = t2412 * t2458;
t2473 = pkin(3) ^ 2;
t2521 = t2427 * t2473;
t2520 = t2429 * t2473;
t2519 = t2431 * t2473;
t2378 = pkin(1) * t2442 - t2576;
t2515 = t2450 * t2378;
t2379 = pkin(1) * t2445 - t2574;
t2513 = t2453 * t2379;
t2380 = pkin(1) * t2448 - t2572;
t2511 = t2456 * t2380;
t2475 = pkin(2) ^ 2;
t2505 = -t2473 / 0.2e1 + t2475 / 0.2e1;
t2504 = pkin(2) * t2413;
t2503 = pkin(2) * t2414;
t2502 = pkin(2) * t2415;
t2501 = t2389 * t2576;
t2500 = t2390 * t2574;
t2499 = t2391 * t2572;
t2361 = t2407 * g(1) + t2410 * g(2);
t2483 = g(3) * t2443 - t2364 * t2452;
t2331 = (t2355 * t2483 + t2361 * t2371) * t2442 - t2451 * (-t2361 * t2355 + t2371 * t2483);
t2495 = t2331 * t2548;
t2362 = t2408 * g(1) + t2411 * g(2);
t2482 = g(3) * t2446 - t2365 * t2455;
t2332 = (t2356 * t2482 + t2362 * t2372) * t2445 - t2454 * (-t2362 * t2356 + t2372 * t2482);
t2494 = t2332 * t2547;
t2363 = t2409 * g(1) + t2412 * g(2);
t2481 = g(3) * t2449 - t2366 * t2458;
t2333 = (t2357 * t2481 + t2363 * t2373) * t2448 - t2457 * (-t2363 * t2357 + t2373 * t2481);
t2493 = t2333 * t2546;
t2392 = sin(t2433);
t2492 = ((-mrSges(3,1) * t2361 - mrSges(3,2) * t2483) * t2395 + t2392 * (-mrSges(3,1) * t2483 + mrSges(3,2) * t2361)) * t2548;
t2393 = sin(t2434);
t2491 = ((-mrSges(3,1) * t2362 - mrSges(3,2) * t2482) * t2396 + t2393 * (-mrSges(3,1) * t2482 + mrSges(3,2) * t2362)) * t2547;
t2394 = sin(t2435);
t2490 = ((-mrSges(3,1) * t2363 - mrSges(3,2) * t2481) * t2397 + t2394 * (-mrSges(3,1) * t2481 + mrSges(3,2) * t2363)) * t2546;
t2489 = t2443 * t2545;
t2488 = t2446 * t2544;
t2487 = t2449 * t2543;
t2486 = t2452 * t2518;
t2485 = t2455 * t2517;
t2484 = t2458 * t2516;
t2346 = t2486 * t2582 - t2508;
t2479 = pkin(2) * t2486 + t2346 * t2450;
t2347 = t2485 * t2582 - t2507;
t2478 = pkin(2) * t2485 + t2347 * t2453;
t2348 = t2484 * t2582 - t2506;
t2477 = pkin(2) * t2484 + t2348 * t2456;
t2476 = 0.1e1 / pkin(2);
t2474 = 1 / pkin(3);
t2470 = 0.2e1 * qJ(3,1);
t2467 = 0.2e1 * qJ(3,2);
t2464 = 0.2e1 * qJ(3,3);
t2462 = -pkin(3) / 0.2e1;
t2422 = -t2473 + t2475;
t2403 = -qJ(2,1) + t2567;
t2402 = qJ(2,1) + t2566;
t2401 = -qJ(2,2) + t2565;
t2400 = qJ(2,2) + t2564;
t2399 = -qJ(2,3) + t2563;
t2398 = qJ(2,3) + t2562;
t2370 = t2418 + t2569 / 0.2e1 + t2462;
t2369 = t2417 + t2570 / 0.2e1 + t2462;
t2368 = t2416 + t2571 / 0.2e1 + t2462;
t2354 = t2502 + t2505 + t2519;
t2353 = t2503 + t2505 + t2520;
t2352 = t2504 + t2505 + t2521;
t2345 = t2573 + (-pkin(3) + t2569 + 0.2e1 * t2418) * t2448;
t2344 = t2575 + (-pkin(3) + t2570 + 0.2e1 * t2417) * t2445;
t2343 = t2577 + (-pkin(3) + t2571 + 0.2e1 * t2416) * t2442;
t2342 = pkin(1) * t2572 + (t2422 + 0.2e1 * t2502 + 0.2e1 * t2519) * t2448;
t2341 = pkin(1) * t2574 + (t2422 + 0.2e1 * t2503 + 0.2e1 * t2520) * t2445;
t2340 = pkin(1) * t2576 + (t2422 + 0.2e1 * t2504 + 0.2e1 * t2521) * t2442;
t1 = [-t2410 * t2489 - t2411 * t2488 - t2412 * t2487 - g(1) * m(4) + (((t2370 * t2522 + t2409 * t2537) * t2578 + (t2409 * t2345 - t2412 * t2477) * t2457 - t2549 + t2409 * t2511) * t2493 + ((t2369 * t2524 + t2408 * t2538) * t2579 + (t2408 * t2344 - t2411 * t2478) * t2454 - t2551 + t2408 * t2513) * t2494 + ((t2368 * t2526 + t2407 * t2539) * t2580 + (t2407 * t2343 - t2410 * t2479) * t2451 - t2553 + t2407 * t2515) * t2495 + (((-t2354 * t2522 - t2409 * t2499) * t2578 + (-t2409 * t2342 + t2348 * t2523) * t2457 + pkin(3) * t2549 - t2380 * t2532) * t2490 + ((-t2353 * t2524 - t2408 * t2500) * t2579 + (-t2408 * t2341 + t2347 * t2525) * t2454 + pkin(3) * t2551 - t2379 * t2534) * t2491 + ((-t2352 * t2526 - t2407 * t2501) * t2580 + (-t2407 * t2340 + t2346 * t2527) * t2451 + pkin(3) * t2553 - t2378 * t2536) * t2492) * t2474) * t2476; t2407 * t2489 + t2408 * t2488 + t2409 * t2487 - g(2) * m(4) + (((-t2370 * t2528 + t2412 * t2537) * t2578 + (t2412 * t2345 + t2409 * t2477) * t2457 + t2550 + t2412 * t2511) * t2493 + ((-t2369 * t2529 + t2411 * t2538) * t2579 + (t2411 * t2344 + t2408 * t2478) * t2454 + t2552 + t2411 * t2513) * t2494 + ((-t2368 * t2530 + t2410 * t2539) * t2580 + (t2410 * t2343 + t2407 * t2479) * t2451 + t2554 + t2410 * t2515) * t2495 + (((t2354 * t2528 - t2412 * t2499) * t2578 + (-t2412 * t2342 - t2348 * t2532) * t2457 - pkin(3) * t2550 - t2380 * t2523) * t2490 + ((t2353 * t2529 - t2411 * t2500) * t2579 + (-t2411 * t2341 - t2347 * t2534) * t2454 - pkin(3) * t2552 - t2379 * t2525) * t2491 + ((t2352 * t2530 - t2410 * t2501) * t2580 + (-t2410 * t2340 - t2346 * t2536) * t2451 - pkin(3) * t2554 - t2378 * t2527) * t2492) * t2474) * t2476; -g(3) * m(4) - t2452 * t2545 - t2455 * t2544 - t2458 * t2543 + ((t2354 * t2449 * t2578 + ((pkin(1) - 0.2e1 * t2496) * t2449 - t2458 * t2423) * t2531 - pkin(3) * ((pkin(1) * t2516 - pkin(3) + t2418) * t2449 - t2423 * t2484)) * t2490 + (t2353 * t2446 * t2579 + ((pkin(1) - 0.2e1 * t2497) * t2446 - t2455 * t2423) * t2533 - pkin(3) * ((pkin(1) * t2517 - pkin(3) + t2417) * t2446 - t2423 * t2485)) * t2491 + (t2352 * t2443 * t2580 + ((pkin(1) - 0.2e1 * t2498) * t2443 - t2452 * t2423) * t2535 - pkin(3) * ((pkin(1) * t2518 - pkin(3) + t2416) * t2443 - t2423 * t2486)) * t2492) * t2474 * t2476 + ((-sin(t2403) - sin(t2402)) * t2555 + (cos(t2403) + cos(t2402)) * t2581 + (-sin(-0.2e1 * qJ(3,1) + t2560) - sin(t2470 + t2561) - 0.2e1 * t2449) * pkin(3) + (-sin(-qJ(3,1) + t2560) - sin(qJ(3,1) + t2561) - sin(t2567) - sin(t2566)) * pkin(2)) / (-t2475 * sin(qJ(2,1) - qJ(3,1)) + pkin(2) * (0.2e1 * t2573 + pkin(2) * t2394 + (sin(t2470 + qJ(2,1)) - t2448) * pkin(3))) * t2333 / 0.2e1 + ((-sin(t2401) - sin(t2400)) * t2555 + (cos(t2401) + cos(t2400)) * t2581 + (-sin(-0.2e1 * qJ(3,2) + t2558) - sin(t2467 + t2559) - 0.2e1 * t2446) * pkin(3) + (-sin(-qJ(3,2) + t2558) - sin(qJ(3,2) + t2559) - sin(t2565) - sin(t2564)) * pkin(2)) / (-t2475 * sin(qJ(2,2) - qJ(3,2)) + pkin(2) * (0.2e1 * t2575 + pkin(2) * t2393 + (sin(t2467 + qJ(2,2)) - t2445) * pkin(3))) * t2332 / 0.2e1 + ((-sin(t2399) - sin(t2398)) * t2555 + (cos(t2399) + cos(t2398)) * t2581 + (-sin(-0.2e1 * qJ(3,3) + t2556) - sin(t2464 + t2557) - 0.2e1 * t2443) * pkin(3) + (-sin(-qJ(3,3) + t2556) - sin(qJ(3,3) + t2557) - sin(t2563) - sin(t2562)) * pkin(2)) / (-t2475 * sin(qJ(2,3) - qJ(3,3)) + pkin(2) * (0.2e1 * t2577 + pkin(2) * t2392 + (sin(t2464 + qJ(2,3)) - t2442) * pkin(3))) * t2331 / 0.2e1;];
taugX  = t1;
