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
% Datum: 2020-08-07 10:08
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRRRR7V2G2A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 09:45:43
% EndTime: 2020-08-07 09:45:45
% DurationCPUTime: 2.18s
% Computational Cost: add. (1173->320), mult. (1818->438), div. (66->14), fcn. (1203->66), ass. (0->211)
t2438 = sin(qJ(3,1));
t2447 = cos(qJ(3,1));
t2558 = -m(3) * pkin(2) - mrSges(2,1);
t2351 = -mrSges(3,1) * t2447 + mrSges(3,2) * t2438 + t2558;
t2367 = t2438 * mrSges(3,1) + t2447 * mrSges(3,2) + mrSges(2,2);
t2439 = sin(qJ(2,1));
t2448 = cos(qJ(2,1));
t2575 = -t2351 * t2448 - t2367 * t2439;
t2435 = sin(qJ(3,2));
t2444 = cos(qJ(3,2));
t2350 = -mrSges(3,1) * t2444 + mrSges(3,2) * t2435 + t2558;
t2366 = t2435 * mrSges(3,1) + t2444 * mrSges(3,2) + mrSges(2,2);
t2436 = sin(qJ(2,2));
t2445 = cos(qJ(2,2));
t2574 = -t2350 * t2445 - t2366 * t2436;
t2432 = sin(qJ(3,3));
t2441 = cos(qJ(3,3));
t2349 = -mrSges(3,1) * t2441 + mrSges(3,2) * t2432 + t2558;
t2365 = t2432 * mrSges(3,1) + t2441 * mrSges(3,2) + mrSges(2,2);
t2433 = sin(qJ(2,3));
t2442 = cos(qJ(2,3));
t2573 = -t2349 * t2442 - t2365 * t2433;
t2572 = 2 * pkin(3);
t2451 = pkin(6) + pkin(5);
t2414 = (pkin(7) + t2451);
t2571 = 2 * t2414;
t2419 = t2442 ^ 2;
t2570 = 0.2e1 * t2419;
t2421 = t2445 ^ 2;
t2569 = 0.2e1 * t2421;
t2423 = t2448 ^ 2;
t2568 = 0.2e1 * t2423;
t2418 = t2441 ^ 2;
t2407 = pkin(3) * t2418;
t2420 = t2444 ^ 2;
t2408 = pkin(3) * t2420;
t2422 = t2447 ^ 2;
t2409 = pkin(3) * t2422;
t2567 = t2432 * pkin(1);
t2566 = t2432 * pkin(3);
t2565 = t2435 * pkin(1);
t2564 = t2435 * pkin(3);
t2563 = t2438 * pkin(1);
t2562 = t2438 * pkin(3);
t2561 = t2441 * pkin(2);
t2404 = t2441 * pkin(3);
t2560 = t2444 * pkin(2);
t2405 = t2444 * pkin(3);
t2559 = t2447 * pkin(2);
t2406 = t2447 * pkin(3);
t2557 = -qJ(3,1) + qJ(1,1);
t2556 = qJ(3,1) + qJ(1,1);
t2555 = -qJ(3,2) + qJ(1,2);
t2554 = qJ(3,2) + qJ(1,2);
t2553 = -qJ(3,3) + qJ(1,3);
t2552 = qJ(3,3) + qJ(1,3);
t2551 = qJ(1,1) + 0.2e1 * qJ(2,1);
t2550 = qJ(1,1) - 0.2e1 * qJ(2,1);
t2549 = qJ(1,2) + 0.2e1 * qJ(2,2);
t2548 = qJ(1,2) - 0.2e1 * qJ(2,2);
t2547 = 0.2e1 * qJ(2,3) + qJ(1,3);
t2546 = qJ(1,3) - 0.2e1 * qJ(2,3);
t2545 = 0.2e1 * pkin(1);
t2434 = sin(qJ(1,3));
t2443 = cos(qJ(1,3));
t2476 = pkin(1) * t2434 - t2443 * t2414;
t2503 = t2432 * t2433;
t2331 = t2476 * t2503 + (t2418 - 0.1e1) * t2434 * pkin(3);
t2429 = legFrame(3,2);
t2398 = sin(t2429);
t2544 = t2331 * t2398;
t2401 = cos(t2429);
t2543 = t2331 * t2401;
t2437 = sin(qJ(1,2));
t2446 = cos(qJ(1,2));
t2475 = pkin(1) * t2437 - t2446 * t2414;
t2502 = t2435 * t2436;
t2332 = t2475 * t2502 + (t2420 - 0.1e1) * t2437 * pkin(3);
t2430 = legFrame(2,2);
t2399 = sin(t2430);
t2542 = t2332 * t2399;
t2402 = cos(t2430);
t2541 = t2332 * t2402;
t2440 = sin(qJ(1,1));
t2449 = cos(qJ(1,1));
t2474 = pkin(1) * t2440 - t2449 * t2414;
t2501 = t2438 * t2439;
t2333 = t2474 * t2501 + (t2422 - 0.1e1) * t2440 * pkin(3);
t2431 = legFrame(1,2);
t2400 = sin(t2431);
t2540 = t2333 * t2400;
t2403 = cos(t2431);
t2539 = t2333 * t2403;
t2491 = pkin(3) * t2503;
t2380 = t2404 + pkin(2);
t2520 = t2380 * t2442;
t2538 = 0.1e1 / (pkin(1) - t2491 + t2520) / t2432;
t2490 = pkin(3) * t2502;
t2381 = t2405 + pkin(2);
t2518 = t2381 * t2445;
t2537 = 0.1e1 / (pkin(1) - t2490 + t2518) / t2435;
t2489 = pkin(3) * t2501;
t2382 = t2406 + pkin(2);
t2516 = t2382 * t2448;
t2536 = 0.1e1 / (pkin(1) - t2489 + t2516) / t2438;
t2358 = t2401 * g(1) - t2398 * g(2);
t2368 = m(2) * pkin(5) + t2451 * m(3) - mrSges(1,2) + mrSges(2,3) + mrSges(3,3);
t2361 = t2368 * g(3);
t2376 = mrSges(1,1) + (m(2) + m(3)) * pkin(1);
t2375 = t2376 * g(3);
t2424 = qJ(2,3) + qJ(3,3);
t2386 = cos(t2424);
t2532 = 0.1e1 / (t2442 * pkin(2) + pkin(3) * t2386 + pkin(1)) * (-t2361 * t2443 + t2434 * (t2573 * g(3) + t2375) + ((-t2376 - t2573) * t2443 - t2434 * t2368) * t2358);
t2359 = t2402 * g(1) - t2399 * g(2);
t2425 = qJ(2,2) + qJ(3,2);
t2387 = cos(t2425);
t2531 = 0.1e1 / (t2445 * pkin(2) + pkin(3) * t2387 + pkin(1)) * (-t2361 * t2446 + t2437 * (t2574 * g(3) + t2375) + ((-t2376 - t2574) * t2446 - t2437 * t2368) * t2359);
t2360 = t2403 * g(1) - t2400 * g(2);
t2426 = qJ(2,1) + qJ(3,1);
t2388 = cos(t2426);
t2530 = 0.1e1 / (t2448 * pkin(2) + pkin(3) * t2388 + pkin(1)) * (-t2361 * t2449 + t2440 * (t2575 * g(3) + t2375) + ((-t2376 - t2575) * t2449 - t2440 * t2368) * t2360);
t2373 = pkin(1) * t2436 - t2564;
t2526 = t2373 * t2444;
t2374 = pkin(1) * t2439 - t2562;
t2525 = t2374 * t2447;
t2454 = pkin(2) / 0.2e1;
t2524 = (t2404 + t2454) * t2432;
t2523 = (t2405 + t2454) * t2435;
t2522 = (t2406 + t2454) * t2438;
t2521 = t2380 * t2398;
t2519 = t2381 * t2399;
t2517 = t2382 * t2400;
t2515 = t2398 * t2434;
t2514 = t2399 * t2437;
t2513 = t2400 * t2440;
t2512 = t2401 * t2380;
t2511 = t2401 * t2434;
t2510 = t2402 * t2381;
t2509 = t2402 * t2437;
t2508 = t2403 * t2382;
t2507 = t2403 * t2440;
t2464 = pkin(3) ^ 2;
t2506 = t2418 * t2464;
t2505 = t2420 * t2464;
t2504 = t2422 * t2464;
t2372 = pkin(1) * t2433 - t2566;
t2500 = t2441 * t2372;
t2466 = pkin(2) ^ 2;
t2498 = -t2464 / 0.2e1 + t2466 / 0.2e1;
t2497 = pkin(2) * t2404;
t2496 = pkin(2) * t2405;
t2495 = pkin(2) * t2406;
t2494 = t2380 * t2566;
t2493 = t2381 * t2564;
t2492 = t2382 * t2562;
t2355 = t2398 * g(1) + t2401 * g(2);
t2473 = g(3) * t2443 + t2358 * t2434;
t2325 = (-t2473 * t2349 + t2355 * t2365) * t2433 - (-t2355 * t2349 - t2473 * t2365) * t2442;
t2488 = t2325 * t2538;
t2356 = t2399 * g(1) + t2402 * g(2);
t2472 = g(3) * t2446 + t2359 * t2437;
t2326 = (-t2472 * t2350 + t2356 * t2366) * t2436 - (-t2356 * t2350 - t2472 * t2366) * t2445;
t2487 = t2326 * t2537;
t2357 = t2400 * g(1) + t2403 * g(2);
t2471 = g(3) * t2449 + t2360 * t2440;
t2327 = (-t2471 * t2351 + t2357 * t2367) * t2439 - (-t2357 * t2351 - t2471 * t2367) * t2448;
t2486 = t2327 * t2536;
t2383 = sin(t2424);
t2485 = ((-mrSges(3,1) * t2355 + t2473 * mrSges(3,2)) * t2386 + t2383 * (t2473 * mrSges(3,1) + mrSges(3,2) * t2355)) * t2538;
t2384 = sin(t2425);
t2484 = ((-mrSges(3,1) * t2356 + t2472 * mrSges(3,2)) * t2387 + t2384 * (t2472 * mrSges(3,1) + mrSges(3,2) * t2356)) * t2537;
t2385 = sin(t2426);
t2483 = ((-mrSges(3,1) * t2357 + t2471 * mrSges(3,2)) * t2388 + t2385 * (t2471 * mrSges(3,1) + mrSges(3,2) * t2357)) * t2536;
t2482 = t2434 * t2503;
t2481 = t2437 * t2502;
t2480 = t2440 * t2501;
t2479 = t2443 * t2532;
t2478 = t2446 * t2531;
t2477 = t2449 * t2530;
t2340 = t2482 * t2572 - t2476;
t2470 = pkin(2) * t2482 + t2340 * t2441;
t2341 = t2481 * t2572 - t2475;
t2469 = pkin(2) * t2481 + t2341 * t2444;
t2342 = t2480 * t2572 - t2474;
t2468 = pkin(2) * t2480 + t2342 * t2447;
t2467 = 0.1e1 / pkin(2);
t2465 = 1 / pkin(3);
t2461 = 0.2e1 * qJ(3,1);
t2458 = 0.2e1 * qJ(3,2);
t2455 = 0.2e1 * qJ(3,3);
t2453 = -pkin(3) / 0.2e1;
t2413 = -t2464 + t2466;
t2394 = -qJ(2,1) + t2557;
t2393 = qJ(2,1) + t2556;
t2392 = -qJ(2,2) + t2555;
t2391 = qJ(2,2) + t2554;
t2390 = -qJ(2,3) + t2553;
t2389 = qJ(2,3) + t2552;
t2364 = t2409 + t2559 / 0.2e1 + t2453;
t2363 = t2408 + t2560 / 0.2e1 + t2453;
t2362 = t2407 + t2561 / 0.2e1 + t2453;
t2348 = t2495 + t2498 + t2504;
t2347 = t2496 + t2498 + t2505;
t2346 = t2497 + t2498 + t2506;
t2339 = t2563 + (-pkin(3) + t2559 + 0.2e1 * t2409) * t2439;
t2338 = t2565 + (-pkin(3) + t2560 + 0.2e1 * t2408) * t2436;
t2337 = t2567 + (-pkin(3) + t2561 + 0.2e1 * t2407) * t2433;
t2336 = pkin(1) * t2562 + (t2413 + 0.2e1 * t2495 + 0.2e1 * t2504) * t2439;
t2335 = pkin(1) * t2564 + (t2413 + 0.2e1 * t2496 + 0.2e1 * t2505) * t2436;
t2334 = pkin(1) * t2566 + (t2413 + 0.2e1 * t2497 + 0.2e1 * t2506) * t2433;
t1 = [t2401 * t2479 + t2402 * t2478 + t2403 * t2477 - g(1) * m(4) + (((t2364 * t2507 + t2400 * t2522) * t2568 + (t2400 * t2339 - t2468 * t2403) * t2448 - t2539 + t2400 * t2525) * t2486 + ((t2363 * t2509 + t2399 * t2523) * t2569 + (t2399 * t2338 - t2469 * t2402) * t2445 - t2541 + t2399 * t2526) * t2487 + ((t2362 * t2511 + t2398 * t2524) * t2570 + (t2398 * t2337 - t2470 * t2401) * t2442 - t2543 + t2398 * t2500) * t2488 + (((-t2348 * t2507 - t2400 * t2492) * t2568 + (-t2400 * t2336 + t2342 * t2508) * t2448 + pkin(3) * t2539 - t2374 * t2517) * t2483 + ((-t2347 * t2509 - t2399 * t2493) * t2569 + (-t2399 * t2335 + t2341 * t2510) * t2445 + pkin(3) * t2541 - t2373 * t2519) * t2484 + ((-t2346 * t2511 - t2398 * t2494) * t2570 + (-t2398 * t2334 + t2340 * t2512) * t2442 + pkin(3) * t2543 - t2372 * t2521) * t2485) * t2465) * t2467; -t2398 * t2479 - t2399 * t2478 - t2400 * t2477 - g(2) * m(4) + (((-t2364 * t2513 + t2403 * t2522) * t2568 + (t2403 * t2339 + t2468 * t2400) * t2448 + t2540 + t2403 * t2525) * t2486 + ((-t2363 * t2514 + t2402 * t2523) * t2569 + (t2402 * t2338 + t2469 * t2399) * t2445 + t2542 + t2402 * t2526) * t2487 + ((-t2362 * t2515 + t2401 * t2524) * t2570 + (t2401 * t2337 + t2470 * t2398) * t2442 + t2544 + t2401 * t2500) * t2488 + (((t2348 * t2513 - t2403 * t2492) * t2568 + (-t2336 * t2403 - t2342 * t2517) * t2448 - pkin(3) * t2540 - t2374 * t2508) * t2483 + ((t2347 * t2514 - t2402 * t2493) * t2569 + (-t2335 * t2402 - t2341 * t2519) * t2445 - pkin(3) * t2542 - t2373 * t2510) * t2484 + ((t2346 * t2515 - t2401 * t2494) * t2570 + (-t2334 * t2401 - t2340 * t2521) * t2442 - pkin(3) * t2544 - t2372 * t2512) * t2485) * t2465) * t2467; -g(3) * m(4) - t2434 * t2532 - t2437 * t2531 - t2440 * t2530 + ((-0.2e1 * t2348 * t2449 * t2423 - ((pkin(1) - 0.2e1 * t2489) * t2449 + t2440 * t2414) * t2516 + pkin(3) * ((pkin(1) * t2501 - pkin(3) + t2409) * t2449 + t2414 * t2480)) * t2483 + (-0.2e1 * t2347 * t2446 * t2421 - ((pkin(1) - 0.2e1 * t2490) * t2446 + t2437 * t2414) * t2518 + pkin(3) * ((pkin(1) * t2502 - pkin(3) + t2408) * t2446 + t2414 * t2481)) * t2484 + (-0.2e1 * t2346 * t2443 * t2419 - ((pkin(1) - 0.2e1 * t2491) * t2443 + t2434 * t2414) * t2520 + pkin(3) * ((pkin(1) * t2503 - pkin(3) + t2407) * t2443 + t2414 * t2482)) * t2485) * t2465 * t2467 + ((cos(t2394) + cos(t2393)) * t2545 + (sin(t2394) + sin(t2393)) * t2571 + (cos(-0.2e1 * qJ(3,1) + t2550) + cos(t2461 + t2551) + 0.2e1 * t2449) * pkin(3) + (cos(-qJ(3,1) + t2550) + cos(qJ(3,1) + t2551) + cos(t2557) + cos(t2556)) * pkin(2)) / (-t2466 * sin(qJ(2,1) - qJ(3,1)) + pkin(2) * (0.2e1 * t2563 + pkin(2) * t2385 + (sin(t2461 + qJ(2,1)) - t2439) * pkin(3))) * t2327 / 0.2e1 + ((cos(t2392) + cos(t2391)) * t2545 + (sin(t2392) + sin(t2391)) * t2571 + (cos(-0.2e1 * qJ(3,2) + t2548) + cos(t2458 + t2549) + 0.2e1 * t2446) * pkin(3) + (cos(-qJ(3,2) + t2548) + cos(qJ(3,2) + t2549) + cos(t2555) + cos(t2554)) * pkin(2)) / (-t2466 * sin(qJ(2,2) - qJ(3,2)) + pkin(2) * (0.2e1 * t2565 + pkin(2) * t2384 + (sin(t2458 + qJ(2,2)) - t2436) * pkin(3))) * t2326 / 0.2e1 + ((cos(t2390) + cos(t2389)) * t2545 + (sin(t2390) + sin(t2389)) * t2571 + (cos(-0.2e1 * qJ(3,3) + t2546) + cos(t2455 + t2547) + 0.2e1 * t2443) * pkin(3) + (cos(-qJ(3,3) + t2546) + cos(qJ(3,3) + t2547) + cos(t2553) + cos(t2552)) * pkin(2)) / (-t2466 * sin(qJ(2,3) - qJ(3,3)) + pkin(2) * (0.2e1 * t2567 + pkin(2) * t2383 + (sin(t2455 + qJ(2,3)) - t2433) * pkin(3))) * t2325 / 0.2e1;];
taugX  = t1;
