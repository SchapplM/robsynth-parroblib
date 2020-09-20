% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P4PRRRR1G2A0
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
%   see P4PRRRR1G2A0_convert_par2_MPV_fixb.m

% Output:
% MMX [4x4]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:00
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P4PRRRR1G2A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(2,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G2A0_inertia_para_pf_mdp: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G2A0_inertia_para_pf_mdp: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G2A0_inertia_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G2A0_inertia_para_pf_mdp: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G2A0_inertia_para_pf_mdp: Koppelpunkt has to be [4x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P4PRRRR1G2A0_inertia_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 10:59:31
% EndTime: 2020-08-07 10:59:50
% DurationCPUTime: 19.95s
% Computational Cost: add. (3594->761), mult. (11709->1411), div. (4620->22), fcn. (13296->26), ass. (0->550)
t2250 = xP(4);
t2184 = sin(t2250);
t2185 = cos(t2250);
t2254 = koppelP(1,2);
t2258 = koppelP(1,1);
t2161 = t2184 * t2258 + t2185 * t2254;
t2165 = -t2184 * t2254 + t2185 * t2258;
t2237 = legFrame(1,2);
t2179 = sin(t2237);
t2183 = cos(t2237);
t2242 = sin(qJ(3,1));
t2243 = sin(qJ(2,1));
t2248 = cos(qJ(3,1));
t2563 = t2243 * t2248;
t2157 = t2179 * t2242 + t2183 * t2563;
t2210 = 0.1e1 / t2243;
t2224 = 0.1e1 / t2248;
t2595 = t2210 * t2224;
t2535 = t2157 * t2595;
t2154 = t2179 * t2563 - t2183 * t2242;
t2538 = t2154 * t2595;
t2355 = -t2161 * t2538 + t2165 * t2535;
t2253 = koppelP(2,2);
t2257 = koppelP(2,1);
t2160 = t2184 * t2257 + t2185 * t2253;
t2164 = -t2184 * t2253 + t2185 * t2257;
t2236 = legFrame(2,2);
t2178 = sin(t2236);
t2182 = cos(t2236);
t2240 = sin(qJ(3,2));
t2241 = sin(qJ(2,2));
t2246 = cos(qJ(3,2));
t2565 = t2241 * t2246;
t2156 = t2178 * t2240 + t2182 * t2565;
t2205 = 0.1e1 / t2241;
t2218 = 0.1e1 / t2246;
t2602 = t2205 * t2218;
t2536 = t2156 * t2602;
t2153 = t2178 * t2565 - t2182 * t2240;
t2539 = t2153 * t2602;
t2356 = -t2160 * t2539 + t2164 * t2536;
t2252 = koppelP(3,2);
t2256 = koppelP(3,1);
t2159 = t2184 * t2256 + t2185 * t2252;
t2163 = -t2184 * t2252 + t2185 * t2256;
t2235 = legFrame(3,2);
t2177 = sin(t2235);
t2181 = cos(t2235);
t2238 = sin(qJ(3,3));
t2239 = sin(qJ(2,3));
t2244 = cos(qJ(3,3));
t2567 = t2239 * t2244;
t2155 = t2177 * t2238 + t2181 * t2567;
t2200 = 0.1e1 / t2239;
t2212 = 0.1e1 / t2244;
t2609 = t2200 * t2212;
t2537 = t2155 * t2609;
t2152 = t2177 * t2567 - t2181 * t2238;
t2540 = t2152 * t2609;
t2357 = -t2159 * t2540 + t2163 * t2537;
t2251 = koppelP(4,2);
t2255 = koppelP(4,1);
t2158 = t2184 * t2255 + t2185 * t2251;
t2162 = -t2184 * t2251 + t2185 * t2255;
t2234 = legFrame(4,2);
t2176 = sin(t2234);
t2180 = cos(t2234);
t2230 = sin(qJ(3,4));
t2231 = sin(qJ(2,4));
t2232 = cos(qJ(3,4));
t2569 = t2231 * t2232;
t2151 = t2176 * t2230 + t2180 * t2569;
t2189 = 0.1e1 / t2231;
t2191 = 0.1e1 / t2232;
t2622 = t2189 * t2191;
t2541 = t2151 * t2622;
t2150 = t2176 * t2569 - t2180 * t2230;
t2542 = t2150 * t2622;
t2358 = -t2158 * t2542 + t2162 * t2541;
t2249 = cos(qJ(2,1));
t2592 = t2210 * t2249;
t2247 = cos(qJ(2,2));
t2599 = t2205 * t2247;
t2245 = cos(qJ(2,3));
t2606 = t2200 * t2245;
t2233 = cos(qJ(2,4));
t2619 = t2189 * t2233;
t2691 = t2355 * t2592 + t2356 * t2599 + t2357 * t2606 + t2358 * t2619;
t2719 = t2691 * MDP(1);
t2188 = t2230 ^ 2;
t2263 = t2232 ^ 2;
t2192 = 0.1e1 / t2263;
t2193 = t2191 * t2192;
t2190 = 0.1e1 / t2231 ^ 2;
t2196 = t2233 ^ 2;
t2611 = t2196 * t2190;
t2528 = t2193 * t2611;
t2353 = t2188 * t2528 - t2191;
t2614 = t2192 * t2230;
t2362 = (0.1e1 + t2611) * t2614;
t2718 = -MDP(10) * t2362 + t2353 * MDP(11);
t2199 = t2238 ^ 2;
t2273 = t2244 ^ 2;
t2213 = 0.1e1 / t2273;
t2214 = t2212 * t2213;
t2201 = 0.1e1 / t2239 ^ 2;
t2217 = t2245 ^ 2;
t2603 = t2201 * t2217;
t2517 = t2214 * t2603;
t2350 = t2199 * t2517 - t2212;
t2586 = t2213 * t2238;
t2361 = (0.1e1 + t2603) * t2586;
t2717 = -MDP(10) * t2361 + t2350 * MDP(11);
t2204 = t2240 ^ 2;
t2277 = t2246 ^ 2;
t2219 = 0.1e1 / t2277;
t2220 = t2218 * t2219;
t2206 = 0.1e1 / t2241 ^ 2;
t2223 = t2247 ^ 2;
t2596 = t2206 * t2223;
t2508 = t2220 * t2596;
t2347 = t2204 * t2508 - t2218;
t2580 = t2219 * t2240;
t2360 = (0.1e1 + t2596) * t2580;
t2716 = -MDP(10) * t2360 + t2347 * MDP(11);
t2209 = t2242 ^ 2;
t2281 = t2248 ^ 2;
t2225 = 0.1e1 / t2281;
t2226 = t2224 * t2225;
t2211 = 0.1e1 / t2243 ^ 2;
t2229 = t2249 ^ 2;
t2589 = t2211 * t2229;
t2499 = t2226 * t2589;
t2344 = t2209 * t2499 - t2224;
t2574 = t2225 * t2242;
t2359 = (0.1e1 + t2589) * t2574;
t2715 = -MDP(10) * t2359 + t2344 * MDP(11);
t2616 = t2191 * t2230;
t2195 = t2233 * t2196;
t2617 = t2190 * t2195;
t2341 = (t2233 + t2617) * t2616;
t2529 = t2192 * t2617;
t2352 = t2188 * t2529 - t2233;
t2714 = -MDP(10) * t2341 + MDP(11) * t2352;
t2588 = t2212 * t2238;
t2216 = t2245 * t2217;
t2604 = t2201 * t2216;
t2339 = (t2245 + t2604) * t2588;
t2518 = t2213 * t2604;
t2349 = t2199 * t2518 - t2245;
t2713 = -MDP(10) * t2339 + MDP(11) * t2349;
t2582 = t2218 * t2240;
t2222 = t2247 * t2223;
t2597 = t2206 * t2222;
t2337 = (t2247 + t2597) * t2582;
t2509 = t2219 * t2597;
t2346 = t2204 * t2509 - t2247;
t2712 = -MDP(10) * t2337 + MDP(11) * t2346;
t2576 = t2224 * t2242;
t2228 = t2249 * t2229;
t2590 = t2211 * t2228;
t2335 = (t2249 + t2590) * t2576;
t2500 = t2225 * t2590;
t2343 = t2209 * t2500 - t2249;
t2711 = -MDP(10) * t2335 + MDP(11) * t2343;
t2593 = t2210 * t2229;
t2336 = (t2243 + t2593) * t2576;
t2600 = t2205 * t2223;
t2338 = (t2241 + t2600) * t2582;
t2607 = t2200 * t2217;
t2340 = (t2239 + t2607) * t2588;
t2620 = t2189 * t2196;
t2342 = (t2231 + t2620) * t2616;
t2694 = t2191 * t2358 + t2212 * t2357 + t2218 * t2356 + t2224 * t2355;
t2575 = t2224 * t2249;
t2504 = t2210 * t2575;
t2427 = t2242 * t2504;
t2581 = t2218 * t2247;
t2513 = t2205 * t2581;
t2439 = t2240 * t2513;
t2587 = t2212 * t2245;
t2522 = t2200 * t2587;
t2451 = t2238 * t2522;
t2615 = t2191 * t2233;
t2533 = t2189 * t2615;
t2463 = t2230 * t2533;
t2693 = t2355 * t2427 + t2356 * t2439 + t2357 * t2451 + t2358 * t2463;
t2692 = -t2355 * t2504 - t2356 * t2513 - t2357 * t2522 - t2358 * t2533;
t2690 = -t2184 * MDP(13) - t2185 * MDP(14) + (t2355 * t2538 + t2356 * t2539 + t2357 * t2540 + t2358 * t2542) * MDP(1);
t2689 = t2185 * MDP(13) - t2184 * MDP(14) + (t2355 * t2535 + t2356 * t2536 + t2357 * t2537 + t2358 * t2541) * MDP(1);
t2688 = 2 * MDP(3);
t2687 = 2 * MDP(4);
t2686 = 2 * MDP(6);
t2685 = 2 * MDP(7);
t2684 = 2 * MDP(8);
t2259 = 1 / pkin(2);
t2683 = 2 * t2259;
t2682 = MDP(2) * t2259;
t2260 = 1 / (pkin(2) ^ 2);
t2681 = MDP(2) * t2260;
t2680 = MDP(3) * t2259;
t2679 = MDP(4) * t2259;
t2678 = MDP(5) * t2259;
t2677 = MDP(5) * t2260;
t2676 = MDP(7) * t2259;
t2675 = MDP(7) * t2260;
t2674 = MDP(8) * t2259;
t2673 = MDP(8) * t2260;
t2672 = MDP(9) * t2259;
t2671 = MDP(9) * t2260;
t2119 = t2158 * t2180 + t2162 * t2176;
t2670 = t2358 * t2119;
t2120 = t2159 * t2181 + t2163 * t2177;
t2669 = t2357 * t2120;
t2121 = t2160 * t2182 + t2164 * t2178;
t2668 = t2356 * t2121;
t2118 = t2161 * t2183 + t2165 * t2179;
t2667 = t2355 * t2118;
t2570 = t2230 * t2233;
t2525 = t2192 * t2570;
t2461 = t2189 * t2525;
t2402 = t2119 * t2461;
t2105 = t2259 * t2402;
t2666 = t2105 * t2191;
t2665 = t2105 * t2233;
t2568 = t2238 * t2245;
t2496 = t2213 * t2568;
t2449 = t2200 * t2496;
t2399 = t2120 * t2449;
t2106 = t2259 * t2399;
t2664 = t2106 * t2212;
t2663 = t2106 * t2245;
t2566 = t2240 * t2247;
t2494 = t2219 * t2566;
t2437 = t2205 * t2494;
t2396 = t2121 * t2437;
t2107 = t2259 * t2396;
t2662 = t2107 * t2218;
t2661 = t2107 * t2247;
t2564 = t2242 * t2249;
t2492 = t2225 * t2564;
t2425 = t2210 * t2492;
t2405 = t2118 * t2425;
t2108 = t2259 * t2405;
t2660 = t2108 * t2224;
t2659 = t2108 * t2249;
t2646 = t2119 * t2191;
t2114 = t2259 * t2646;
t2658 = t2114 * t2189;
t2657 = t2114 * t2191;
t2643 = t2120 * t2212;
t2115 = t2259 * t2643;
t2656 = t2115 * t2200;
t2655 = t2115 * t2212;
t2640 = t2121 * t2218;
t2116 = t2259 * t2640;
t2654 = t2116 * t2205;
t2653 = t2116 * t2218;
t2649 = t2118 * t2224;
t2117 = t2259 * t2649;
t2652 = t2117 * t2210;
t2651 = t2117 * t2224;
t2650 = t2118 * t2179;
t2648 = t2118 * t2226;
t2647 = t2119 * t2176;
t2645 = t2119 * t2193;
t2644 = t2120 * t2177;
t2642 = t2120 * t2214;
t2641 = t2121 * t2178;
t2639 = t2121 * t2220;
t2638 = t2150 * t2180;
t2637 = t2151 * t2176;
t2636 = t2152 * t2181;
t2635 = t2153 * t2182;
t2634 = t2154 * t2183;
t2633 = t2155 * t2177;
t2632 = t2156 * t2178;
t2631 = t2157 * t2179;
t2630 = t2176 * t2192;
t2629 = t2177 * t2213;
t2628 = t2178 * t2219;
t2627 = t2179 * t2225;
t2626 = t2180 * t2192;
t2625 = t2181 * t2213;
t2624 = t2182 * t2219;
t2623 = t2183 * t2225;
t2621 = t2189 * t2192;
t2618 = t2190 * t2192;
t2613 = t2193 * t2233;
t2612 = t2195 * t2230;
t2610 = t2196 * t2230;
t2608 = t2200 * t2213;
t2605 = t2201 * t2213;
t2601 = t2205 * t2219;
t2598 = t2206 * t2219;
t2594 = t2210 * t2225;
t2591 = t2211 * t2225;
t2585 = t2214 * t2245;
t2584 = t2216 * t2238;
t2583 = t2217 * t2238;
t2579 = t2220 * t2247;
t2578 = t2222 * t2240;
t2577 = t2223 * t2240;
t2573 = t2226 * t2249;
t2572 = t2228 * t2242;
t2571 = t2229 * t2242;
t2562 = MDP(6) * t2683;
t2561 = t2260 * t2686;
t2560 = 2 * t2675;
t2559 = 2 * t2673;
t2558 = t2358 * t2621;
t2557 = t2357 * t2608;
t2556 = t2356 * t2601;
t2555 = t2355 * t2594;
t2554 = t2105 * t2622;
t2553 = t2106 * t2609;
t2552 = t2107 * t2602;
t2551 = t2108 * t2595;
t2550 = t2210 * t2648;
t2549 = t2211 * t2648;
t2548 = t2189 * t2645;
t2547 = t2190 * t2645;
t2546 = t2200 * t2642;
t2545 = t2201 * t2642;
t2544 = t2205 * t2639;
t2543 = t2206 * t2639;
t2534 = t2188 * t2618;
t2532 = t2189 * t2614;
t2531 = t2192 * t2619;
t2530 = t2190 * t2615;
t2527 = 0.1e1 / t2263 ^ 2 * t2611;
t2526 = t2190 * t2610;
t2524 = t2193 * t2570;
t2523 = t2199 * t2605;
t2521 = t2200 * t2586;
t2520 = t2213 * t2606;
t2519 = t2201 * t2587;
t2516 = 0.1e1 / t2273 ^ 2 * t2603;
t2515 = t2201 * t2583;
t2514 = t2204 * t2598;
t2512 = t2205 * t2580;
t2511 = t2219 * t2599;
t2510 = t2206 * t2581;
t2507 = 0.1e1 / t2277 ^ 2 * t2596;
t2506 = t2206 * t2577;
t2505 = t2209 * t2591;
t2503 = t2210 * t2574;
t2502 = t2225 * t2592;
t2501 = t2211 * t2575;
t2498 = 0.1e1 / t2281 ^ 2 * t2589;
t2497 = t2211 * t2571;
t2495 = t2214 * t2568;
t2493 = t2220 * t2566;
t2491 = t2226 * t2564;
t2477 = t2358 * t2525;
t2476 = t2357 * t2496;
t2475 = t2356 * t2494;
t2474 = t2355 * t2492;
t2473 = t2105 * t2533;
t2472 = t2106 * t2522;
t2471 = t2107 * t2513;
t2470 = t2108 * t2504;
t2186 = t2188 ^ 2;
t2469 = t2186 * t2527;
t2187 = t2230 * t2188;
t2468 = t2187 * t2528;
t2467 = t2187 * t2190 * t2613;
t2466 = t2188 * t2189 * t2613;
t2465 = t2233 * t2534;
t2464 = t2188 * t2527;
t2462 = t2196 * t2532;
t2460 = t2189 * t2524;
t2459 = t2193 * t2526;
t2458 = t2190 * t2524;
t2197 = t2199 ^ 2;
t2457 = t2197 * t2516;
t2198 = t2238 * t2199;
t2456 = t2198 * t2517;
t2455 = t2198 * t2201 * t2585;
t2454 = t2199 * t2200 * t2585;
t2453 = t2245 * t2523;
t2452 = t2199 * t2516;
t2450 = t2217 * t2521;
t2448 = t2200 * t2495;
t2447 = t2214 * t2515;
t2446 = t2201 * t2495;
t2202 = t2204 ^ 2;
t2445 = t2202 * t2507;
t2203 = t2240 * t2204;
t2444 = t2203 * t2508;
t2443 = t2203 * t2206 * t2579;
t2442 = t2204 * t2205 * t2579;
t2441 = t2247 * t2514;
t2440 = t2204 * t2507;
t2438 = t2223 * t2512;
t2436 = t2205 * t2493;
t2435 = t2220 * t2506;
t2434 = t2206 * t2493;
t2207 = t2209 ^ 2;
t2433 = t2207 * t2498;
t2208 = t2242 * t2209;
t2432 = t2208 * t2499;
t2431 = t2208 * t2211 * t2573;
t2430 = t2209 * t2210 * t2573;
t2429 = t2249 * t2505;
t2428 = t2209 * t2498;
t2426 = t2229 * t2503;
t2424 = t2210 * t2491;
t2423 = t2226 * t2497;
t2422 = t2211 * t2491;
t2368 = t2211 * t2259 * t2492;
t2369 = t2206 * t2259 * t2494;
t2370 = t2201 * t2259 * t2496;
t2371 = t2190 * t2259 * t2525;
t2421 = t2150 * t2371 + t2152 * t2370 + t2153 * t2369 + t2154 * t2368;
t2420 = t2151 * t2371 + t2155 * t2370 + t2156 * t2369 + t2157 * t2368;
t2419 = t2358 * t2462;
t2418 = t2357 * t2450;
t2417 = t2356 * t2438;
t2416 = t2355 * t2426;
t2415 = t2105 * t2187 * t2531;
t2414 = t2188 * t2473;
t2413 = t2106 * t2198 * t2520;
t2412 = t2199 * t2472;
t2411 = t2107 * t2203 * t2511;
t2410 = t2204 * t2471;
t2409 = t2108 * t2208 * t2502;
t2408 = t2209 * t2470;
t2407 = t2498 * t2650;
t2406 = t2118 * t2426;
t2404 = t2527 * t2647;
t2403 = t2119 * t2462;
t2401 = t2516 * t2644;
t2400 = t2120 * t2450;
t2398 = t2507 * t2641;
t2397 = t2121 * t2438;
t2395 = t2176 * t2461;
t2394 = t2177 * t2449;
t2393 = t2178 * t2437;
t2392 = t2179 * t2425;
t2391 = t2180 * t2469;
t2390 = t2180 * t2468;
t2389 = t2180 * t2466;
t2388 = t2180 * t2464;
t2387 = t2180 * t2461;
t2386 = t2181 * t2457;
t2385 = t2181 * t2456;
t2384 = t2181 * t2454;
t2383 = t2181 * t2452;
t2382 = t2181 * t2449;
t2381 = t2182 * t2445;
t2380 = t2182 * t2444;
t2379 = t2182 * t2442;
t2378 = t2182 * t2440;
t2377 = t2182 * t2437;
t2376 = t2183 * t2433;
t2375 = t2183 * t2432;
t2374 = t2183 * t2430;
t2373 = t2183 * t2428;
t2372 = t2183 * t2425;
t2366 = t2150 * t2176 - t2151 * t2180;
t2365 = t2152 * t2177 - t2155 * t2181;
t2364 = t2153 * t2178 - t2156 * t2182;
t2363 = t2154 * t2179 - t2157 * t2183;
t2354 = t2188 * t2192 * t2620 - t2231;
t2351 = t2199 * t2213 * t2607 - t2239;
t2348 = t2204 * t2219 * t2600 - t2241;
t2345 = t2209 * t2225 * t2593 - t2243;
t2334 = t2118 * t2359;
t2333 = t2119 * t2362;
t2332 = t2120 * t2361;
t2331 = t2121 * t2360;
t2330 = t2114 * t2463 + t2105;
t2329 = t2115 * t2451 + t2106;
t2328 = t2116 * t2439 + t2107;
t2327 = t2117 * t2427 + t2108;
t2326 = t2118 * t2344;
t2325 = t2119 * t2353;
t2324 = t2120 * t2350;
t2323 = t2121 * t2347;
t2322 = t2176 * t2354;
t2321 = t2177 * t2351;
t2320 = t2178 * t2348;
t2319 = t2179 * t2345;
t2318 = t2180 * t2354;
t2317 = t2181 * t2351;
t2316 = t2182 * t2348;
t2315 = t2183 * t2345;
t2314 = t2176 * t2342;
t2313 = t2177 * t2340;
t2312 = t2178 * t2338;
t2311 = t2179 * t2336;
t2310 = t2180 * t2342;
t2309 = t2181 * t2340;
t2308 = t2182 * t2338;
t2307 = t2183 * t2336;
t2288 = t2150 * t2530 + t2152 * t2519 + t2153 * t2510 + t2154 * t2501;
t2306 = ((-t2154 * t2249 + t2183 * t2572) * t2591 + (-t2153 * t2247 + t2182 * t2578) * t2598 + (-t2152 * t2245 + t2181 * t2584) * t2605 + (-t2150 * t2233 + t2180 * t2612) * t2618) * t2680 + ((-t2183 * t2571 + t2154) * t2594 + (-t2182 * t2577 + t2153) * t2601 + (-t2181 * t2583 + t2152) * t2608 + (-t2180 * t2610 + t2150) * t2621) * t2679 + t2288 * MDP(1) + (-t2180 * t2465 - t2181 * t2453 - t2182 * t2441 - t2183 * t2429) * t2561 + (-t2180 * t2467 - t2181 * t2455 - t2182 * t2443 - t2183 * t2431) * t2677 + (-t2180 * t2458 - t2181 * t2446 - t2182 * t2434 - t2183 * t2422) * t2681 + (t2180 * t2532 + t2181 * t2521 + t2182 * t2512 + t2183 * t2503) * t2675 + (t2180 * t2622 + t2181 * t2609 + t2182 * t2602 + t2183 * t2595) * t2673;
t2287 = t2151 * t2530 + t2155 * t2519 + t2156 * t2510 + t2157 * t2501;
t2305 = ((-t2157 * t2249 - t2179 * t2572) * t2591 + (-t2156 * t2247 - t2178 * t2578) * t2598 + (-t2155 * t2245 - t2177 * t2584) * t2605 + (-t2151 * t2233 - t2176 * t2612) * t2618) * t2680 + ((t2179 * t2571 + t2157) * t2594 + (t2178 * t2577 + t2156) * t2601 + (t2177 * t2583 + t2155) * t2608 + (t2176 * t2610 + t2151) * t2621) * t2679 + t2287 * MDP(1) + (t2176 * t2465 + t2177 * t2453 + t2178 * t2441 + t2179 * t2429) * t2561 + (t2176 * t2467 + t2177 * t2455 + t2178 * t2443 + t2179 * t2431) * t2677 + (t2176 * t2458 + t2177 * t2446 + t2178 * t2434 + t2179 * t2422) * t2681 + (-t2176 * t2532 - t2177 * t2521 - t2178 * t2512 - t2179 * t2503) * t2675 + (-t2176 * t2622 - t2177 * t2609 - t2178 * t2602 - t2179 * t2595) * t2673;
t2304 = (-t2363 * t2423 - t2364 * t2435 - t2365 * t2447 - t2366 * t2459) * t2680 + (t2363 * t2424 + t2364 * t2436 + t2365 * t2448 + t2366 * t2460) * t2679 + (t2150 * t2151 * t2618 + t2152 * t2155 * t2605 + t2153 * t2156 * t2598 + t2154 * t2157 * t2591) * MDP(1) + (-t2176 * t2390 - t2177 * t2385 - t2178 * t2380 - t2179 * t2375) * t2561 + (t2176 * t2389 + t2177 * t2384 + t2178 * t2379 + t2179 * t2374) * t2560 + (t2176 * t2387 + t2177 * t2382 + t2178 * t2377 + t2179 * t2372) * t2559 + (-t2176 * t2391 - t2177 * t2386 - t2178 * t2381 - t2179 * t2376) * t2677 + (-t2176 * t2388 - t2177 * t2383 - t2178 * t2378 - t2179 * t2373) * t2681 + (-t2176 * t2626 - t2177 * t2625 - t2178 * t2624 - t2179 * t2623) * t2671;
t2303 = t2114 * t2188 * t2531 + t2105 * t2616;
t2302 = t2199 * t2115 * t2520 + t2106 * t2588;
t2301 = t2204 * t2116 * t2511 + t2107 * t2582;
t2300 = t2209 * t2117 * t2502 + t2108 * t2576;
t2299 = t2589 + t2596 + t2603 + t2611;
t2296 = t2176 * t2718;
t2295 = t2177 * t2717;
t2294 = t2178 * t2716;
t2293 = t2179 * t2715;
t2292 = t2180 * t2718;
t2291 = t2181 * t2717;
t2290 = t2182 * t2716;
t2289 = t2183 * t2715;
t2286 = t2288 * MDP(10);
t2285 = t2287 * MDP(10);
t2175 = t2183 ^ 2;
t2174 = t2182 ^ 2;
t2173 = t2181 ^ 2;
t2172 = t2180 ^ 2;
t2171 = t2179 ^ 2;
t2170 = t2178 ^ 2;
t2169 = t2177 ^ 2;
t2168 = t2176 ^ 2;
t2166 = (t2184 ^ 2 + t2185 ^ 2) * MDP(15);
t2149 = t2259 * t2315;
t2148 = t2259 * t2319;
t2147 = t2259 * t2316;
t2146 = t2259 * t2320;
t2145 = t2259 * t2317;
t2144 = t2259 * t2321;
t2143 = t2259 * t2318;
t2142 = t2259 * t2322;
t2135 = t2259 * t2307;
t2134 = t2259 * t2311;
t2133 = t2259 * t2308;
t2132 = t2259 * t2312;
t2131 = t2259 * t2309;
t2130 = t2259 * t2313;
t2127 = t2259 * t2310;
t2126 = t2259 * t2314;
t2083 = t2108 * t2564 - t2117 * t2563;
t2082 = t2107 * t2566 - t2116 * t2565;
t2081 = t2106 * t2568 - t2115 * t2567;
t2080 = -t2117 * t2242 * t2243 - t2248 * t2659;
t2079 = -t2116 * t2240 * t2241 - t2246 * t2661;
t2078 = -t2115 * t2238 * t2239 - t2244 * t2663;
t2077 = t2105 * t2570 - t2114 * t2569;
t2076 = -t2114 * t2230 * t2231 - t2232 * t2665;
t1 = [(t2150 ^ 2 * t2618 + t2152 ^ 2 * t2605 + t2153 ^ 2 * t2598 + t2154 ^ 2 * t2591) * MDP(1) + (t2127 * t2542 + t2131 * t2540 + t2133 * t2539 + t2135 * t2538) * MDP(10) + (-t2143 * t2542 - t2145 * t2540 - t2147 * t2539 - t2149 * t2538) * MDP(11) + t2166 + (-t2154 * t2289 - t2153 * t2290 - t2152 * t2291 - t2150 * t2292 + (t2423 * t2634 + t2435 * t2635 + t2447 * t2636 + t2459 * t2638) * t2688 + (-t2424 * t2634 - t2436 * t2635 - t2448 * t2636 - t2460 * t2638) * t2687) * t2259 + ((t2172 * t2464 + t2173 * t2452 + t2174 * t2440 + t2175 * t2428) * MDP(2) + (t2172 * t2469 + t2173 * t2457 + t2174 * t2445 + t2175 * t2433) * MDP(5) + (t2172 * t2192 + t2173 * t2213 + t2174 * t2219 + t2175 * t2225) * MDP(9) + (t2172 * t2468 + t2173 * t2456 + t2174 * t2444 + t2175 * t2432) * t2686 + (-t2172 * t2466 - t2173 * t2454 - t2174 * t2442 - t2175 * t2430) * t2685 + (-t2172 * t2461 - t2173 * t2449 - t2174 * t2437 - t2175 * t2425) * t2684) * t2260; (-t2126 * t2542 - t2130 * t2540 - t2132 * t2539 - t2134 * t2538) * MDP(10) + (t2142 * t2542 + t2144 * t2540 + t2146 * t2539 + t2148 * t2538) * MDP(11) + (-t2151 * t2292 - t2155 * t2291 - t2156 * t2290 - t2157 * t2289) * t2259 + t2304; t2421 * MDP(11) + (-t2714 * t2180 - t2713 * t2181 - t2712 * t2182 - t2711 * t2183 - t2286) * t2259 + t2306; (-t2105 * t2387 - t2106 * t2382 - t2107 * t2377 - t2108 * t2372) * t2682 + (-t2150 * t2473 - t2152 * t2472 - t2153 * t2471 - t2154 * t2470 + (t2180 * t2419 + t2181 * t2418 + t2182 * t2417 + t2183 * t2416) * t2259) * MDP(3) + (t2150 * t2666 + t2152 * t2664 + t2153 * t2662 + t2154 * t2660 + (-t2180 * t2477 - t2181 * t2476 - t2182 * t2475 - t2183 * t2474) * t2259) * MDP(4) + (-t2180 * t2415 - t2181 * t2413 - t2182 * t2411 - t2183 * t2409) * t2678 + (-t2180 * t2414 - t2181 * t2412 - t2182 * t2410 - t2183 * t2408) * t2562 + (t2180 * t2303 + t2181 * t2302 + t2182 * t2301 + t2183 * t2300) * t2676 + (t2180 * t2330 + t2181 * t2329 + t2182 * t2328 + t2183 * t2327) * t2674 + (-t2180 * t2657 - t2181 * t2655 - t2182 * t2653 - t2183 * t2651) * t2672 + (t2076 * t2542 + t2078 * t2540 + t2079 * t2539 + t2080 * t2538 + (t2307 * t2355 + t2308 * t2356 + t2309 * t2357 + t2310 * t2358) * t2259) * MDP(10) + (t2077 * t2542 + t2081 * t2540 + t2082 * t2539 + t2083 * t2538 + (-t2315 * t2355 - t2316 * t2356 - t2317 * t2357 - t2318 * t2358) * t2259) * MDP(11) + t2690; (t2127 * t2541 + t2131 * t2537 + t2133 * t2536 + t2135 * t2535) * MDP(10) + (-t2143 * t2541 - t2145 * t2537 - t2147 * t2536 - t2149 * t2535) * MDP(11) + (t2150 * t2296 + t2152 * t2295 + t2153 * t2294 + t2154 * t2293) * t2259 + t2304; (t2151 ^ 2 * t2618 + t2155 ^ 2 * t2605 + t2156 ^ 2 * t2598 + t2157 ^ 2 * t2591) * MDP(1) + (-t2126 * t2541 - t2130 * t2537 - t2132 * t2536 - t2134 * t2535) * MDP(10) + (t2142 * t2541 + t2144 * t2537 + t2146 * t2536 + t2148 * t2535) * MDP(11) + t2166 + (t2157 * t2293 + t2156 * t2294 + t2155 * t2295 + t2151 * t2296 + (-t2423 * t2631 - t2435 * t2632 - t2447 * t2633 - t2459 * t2637) * t2688 + (t2424 * t2631 + t2436 * t2632 + t2448 * t2633 + t2460 * t2637) * t2687) * t2259 + ((t2168 * t2464 + t2169 * t2452 + t2170 * t2440 + t2171 * t2428) * MDP(2) + (t2168 * t2469 + t2169 * t2457 + t2170 * t2445 + t2171 * t2433) * MDP(5) + (t2168 * t2192 + t2169 * t2213 + t2170 * t2219 + t2171 * t2225) * MDP(9) + (t2168 * t2468 + t2169 * t2456 + t2170 * t2444 + t2171 * t2432) * t2686 + (-t2168 * t2466 - t2169 * t2454 - t2170 * t2442 - t2171 * t2430) * t2685 + (-t2168 * t2461 - t2169 * t2449 - t2170 * t2437 - t2171 * t2425) * t2684) * t2260; t2420 * MDP(11) + (t2714 * t2176 + t2713 * t2177 + t2712 * t2178 + t2711 * t2179 - t2285) * t2259 + t2305; (t2105 * t2395 + t2106 * t2394 + t2107 * t2393 + t2108 * t2392) * t2682 + (-t2151 * t2473 - t2155 * t2472 - t2156 * t2471 - t2157 * t2470 + (-t2176 * t2419 - t2177 * t2418 - t2178 * t2417 - t2179 * t2416) * t2259) * MDP(3) + (t2151 * t2666 + t2155 * t2664 + t2156 * t2662 + t2157 * t2660 + (t2176 * t2477 + t2177 * t2476 + t2178 * t2475 + t2179 * t2474) * t2259) * MDP(4) + (t2176 * t2415 + t2177 * t2413 + t2178 * t2411 + t2179 * t2409) * t2678 + (t2176 * t2414 + t2177 * t2412 + t2178 * t2410 + t2179 * t2408) * t2562 + (-t2176 * t2303 - t2177 * t2302 - t2178 * t2301 - t2179 * t2300) * t2676 + (-t2176 * t2330 - t2177 * t2329 - t2178 * t2328 - t2179 * t2327) * t2674 + (t2176 * t2657 + t2177 * t2655 + t2178 * t2653 + t2179 * t2651) * t2672 + (t2076 * t2541 + t2078 * t2537 + t2079 * t2536 + t2080 * t2535 + (-t2311 * t2355 - t2312 * t2356 - t2313 * t2357 - t2314 * t2358) * t2259) * MDP(10) + (t2077 * t2541 + t2081 * t2537 + t2082 * t2536 + t2083 * t2535 + (t2319 * t2355 + t2320 * t2356 + t2321 * t2357 + t2322 * t2358) * t2259) * MDP(11) + t2689; (t2127 * t2619 + t2131 * t2606 + t2133 * t2599 + t2135 * t2592) * MDP(10) + (-t2143 * t2619 - t2145 * t2606 - t2147 * t2599 - t2149 * t2592 + t2421) * MDP(11) - t2259 * t2286 + t2306; (-t2126 * t2619 - t2130 * t2606 - t2132 * t2599 - t2134 * t2592) * MDP(10) + (t2142 * t2619 + t2144 * t2606 + t2146 * t2599 + t2148 * t2592 + t2420) * MDP(11) - t2259 * t2285 + t2305; t2299 * MDP(1) + MDP(15) + ((t2591 + t2598 + t2605 + t2618) * MDP(2) + (t2505 + t2514 + t2523 + t2534) * MDP(5) + (t2190 * t2616 + t2201 * t2588 + t2206 * t2582 + t2211 * t2576) * t2686) * t2260 + ((-t2191 * t2611 - t2212 * t2603 - t2218 * t2596 - t2224 * t2589) * MDP(3) + (t2504 + t2513 + t2522 + t2533) * MDP(4) - t2299 * MDP(10) + (t2191 * t2526 + t2212 * t2515 + t2218 * t2506 + t2224 * t2497) * MDP(11)) * t2683; t2719 + (-t2105 * t2620 - t2106 * t2607 - t2107 * t2600 - t2108 * t2593) * MDP(3) + (t2659 + t2661 + t2663 + t2665) * MDP(4) + (t2076 * t2619 + t2078 * t2606 + t2079 * t2599 + t2080 * t2592) * MDP(10) + (t2077 * t2619 + t2081 * t2606 + t2082 * t2599 + t2083 * t2592) * MDP(11) + ((t2551 + t2552 + t2553 + t2554) * MDP(2) + t2692 * MDP(3) + t2694 * MDP(4) + (t2188 * t2554 + t2199 * t2553 + t2204 * t2552 + t2209 * t2551) * MDP(5) + (t2105 * t2189 * t2230 + t2106 * t2200 * t2238 + t2107 * t2205 * t2240 + t2108 * t2210 * t2242) * t2686 + (-t2576 * t2652 - t2582 * t2654 - t2588 * t2656 - t2616 * t2658) * MDP(7) + (-t2652 - t2654 - t2656 - t2658) * MDP(8) - t2691 * MDP(10) + t2693 * MDP(11)) * t2259; (-t2118 * t2373 - t2119 * t2388 - t2120 * t2383 - t2121 * t2378) * t2681 + ((-t2154 * t2549 + t2183 * t2555) * t2571 + (-t2153 * t2543 + t2182 * t2556) * t2577 + (-t2152 * t2545 + t2181 * t2557) * t2583 + (-t2150 * t2547 + t2180 * t2558) * t2610) * t2680 + ((t2154 * t2550 - t2355 * t2623) * t2564 + (t2153 * t2544 - t2356 * t2624) * t2566 + (t2152 * t2546 - t2357 * t2625) * t2568 + (t2150 * t2548 - t2358 * t2626) * t2570) * t2679 + (-t2118 * t2376 - t2119 * t2391 - t2120 * t2386 - t2121 * t2381) * t2677 + (-t2118 * t2375 - t2119 * t2390 - t2120 * t2385 - t2121 * t2380) * t2561 + (t2118 * t2374 + t2119 * t2389 + t2120 * t2384 + t2121 * t2379) * t2560 + (t2118 * t2372 + t2119 * t2387 + t2120 * t2382 + t2121 * t2377) * t2559 + (-t2118 * t2623 - t2119 * t2626 - t2120 * t2625 - t2121 * t2624) * t2671 + (t2358 * t2127 + t2357 * t2131 + t2356 * t2133 + t2355 * t2135 + (-t2150 * t2333 - t2152 * t2332 - t2153 * t2331 - t2154 * t2334) * t2259) * MDP(10) + (-t2358 * t2143 - t2357 * t2145 - t2356 * t2147 - t2355 * t2149 + (t2150 * t2325 + t2152 * t2324 + t2153 * t2323 + t2154 * t2326) * t2259) * MDP(11) + t2690; (t2188 * t2404 + t2199 * t2401 + t2204 * t2398 + t2209 * t2407) * t2681 + ((-t2157 * t2549 - t2179 * t2555) * t2571 + (-t2156 * t2543 - t2178 * t2556) * t2577 + (-t2155 * t2545 - t2177 * t2557) * t2583 + (-t2151 * t2547 - t2176 * t2558) * t2610) * t2680 + ((t2157 * t2550 + t2355 * t2627) * t2564 + (t2156 * t2544 + t2356 * t2628) * t2566 + (t2155 * t2546 + t2357 * t2629) * t2568 + (t2151 * t2548 + t2358 * t2630) * t2570) * t2679 + (t2186 * t2404 + t2197 * t2401 + t2202 * t2398 + t2207 * t2407) * t2677 + (t2432 * t2650 + t2444 * t2641 + t2456 * t2644 + t2468 * t2647) * t2561 + (-t2430 * t2650 - t2442 * t2641 - t2454 * t2644 - t2466 * t2647) * t2560 + (-t2118 * t2392 - t2119 * t2395 - t2120 * t2394 - t2121 * t2393) * t2559 + (t2118 * t2627 + t2119 * t2630 + t2120 * t2629 + t2121 * t2628) * t2671 + (-t2358 * t2126 - t2357 * t2130 - t2356 * t2132 - t2355 * t2134 + (-t2151 * t2333 - t2155 * t2332 - t2156 * t2331 - t2157 * t2334) * t2259) * MDP(10) + (t2358 * t2142 + t2357 * t2144 + t2356 * t2146 + t2355 * t2148 + (t2151 * t2325 + t2155 * t2324 + t2156 * t2323 + t2157 * t2326) * t2259) * MDP(11) + t2689; t2719 + (t2118 * t2422 + t2119 * t2458 + t2120 * t2446 + t2121 * t2434) * t2681 + (-t2118 * t2242 * t2500 - t2119 * t2230 * t2529 - t2120 * t2238 * t2518 - t2121 * t2240 * t2509 + t2692) * t2680 + (t2397 + t2400 + t2403 + t2406 + t2694) * t2679 + (t2118 * t2431 + t2119 * t2467 + t2120 * t2455 + t2121 * t2443) * t2677 + (t2118 * t2429 + t2119 * t2465 + t2120 * t2453 + t2121 * t2441) * t2561 + (-t2118 * t2503 - t2119 * t2532 - t2120 * t2521 - t2121 * t2512) * t2675 + (-t2118 * t2595 - t2119 * t2622 - t2120 * t2609 - t2121 * t2602) * t2673 + ((-t2118 * t2335 - t2119 * t2341 - t2120 * t2339 - t2121 * t2337 - t2691) * MDP(10) + (t2118 * t2343 + t2119 * t2352 + t2120 * t2349 + t2121 * t2346 + t2693) * MDP(11)) * t2259; (t2355 ^ 2 + t2356 ^ 2 + t2357 ^ 2 + t2358 ^ 2) * MDP(1) + (t2105 * t2402 + t2106 * t2399 + t2107 * t2396 + t2108 * t2405) * t2682 + (-t2358 * t2665 - t2357 * t2663 - t2356 * t2661 - t2355 * t2659 + (-t2355 * t2406 - t2356 * t2397 - t2357 * t2400 - t2358 * t2403) * t2259) * MDP(3) + (t2358 * t2105 * t2231 + t2357 * t2106 * t2239 + t2356 * t2107 * t2241 + t2355 * t2108 * t2243 + (t2118 * t2474 + t2119 * t2477 + t2120 * t2476 + t2121 * t2475) * t2259) * MDP(4) + (t2118 * t2409 + t2119 * t2415 + t2120 * t2413 + t2121 * t2411) * t2678 + (t2118 * t2408 + t2119 * t2414 + t2120 * t2412 + t2121 * t2410) * t2562 + (-t2118 * t2300 - t2119 * t2303 - t2120 * t2302 - t2121 * t2301) * t2676 + (-t2118 * t2327 - t2119 * t2330 - t2120 * t2329 - t2121 * t2328) * t2674 + (t2114 * t2646 + t2115 * t2643 + t2116 * t2640 + t2117 * t2649) * t2672 + (t2358 * t2076 + t2357 * t2078 + t2356 * t2079 + t2355 * t2080 + (-t2336 * t2667 - t2338 * t2668 - t2340 * t2669 - t2342 * t2670) * t2259) * MDP(10) + (t2358 * t2077 + t2357 * t2081 + t2356 * t2082 + t2355 * t2083 + (t2345 * t2667 + t2348 * t2668 + t2351 * t2669 + t2354 * t2670) * t2259) * MDP(11) + MDP(12);];
%% Postprocessing: Reshape Output
% From vec2mat_4_matlab.m
res = [t1(1), t1(2), t1(3), t1(4); t1(5), t1(6), t1(7), t1(8); t1(9), t1(10), t1(11), t1(12); t1(13), t1(14), t1(15), t1(16);];
MMX  = res;
