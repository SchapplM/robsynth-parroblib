% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P4PRRR1G1P1A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [11x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P4PRRR1G1P1A0_convert_par2_MPV_fixb.m

% Output:
% taugX [4x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-02 20:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P4PRRR1G1P1A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(7,1),zeros(11,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRR1G1P1A0_gravload_para_pf_mdp: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRR1G1P1A0_gravload_para_pf_mdp: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRR1G1P1A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P4PRRR1G1P1A0_gravload_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRR1G1P1A0_gravload_para_pf_mdp: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRR1G1P1A0_gravload_para_pf_mdp: Koppelpunkt has to be [4x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [11 1]), ...
  'P4PRRR1G1P1A0_gravload_para_pf_mdp: MDP has to be [11x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-02 20:15:32
% EndTime: 2020-03-02 20:15:35
% DurationCPUTime: 2.34s
% Computational Cost: add. (1690->166), mult. (1331->290), div. (168->6), fcn. (1440->26), ass. (0->140)
t2462 = legFrame(4,3);
t2448 = sin(t2462);
t2452 = cos(t2462);
t2410 = t2448 * g(1) - t2452 * g(2);
t2414 = t2452 * g(1) + t2448 * g(2);
t2458 = pkin(7) + qJ(2,4);
t2438 = qJ(3,4) + t2458;
t2420 = sin(t2438);
t2421 = cos(t2438);
t2358 = t2410 * t2421 + t2414 * t2420;
t2436 = sin(t2458);
t2437 = cos(t2458);
t2512 = 0.1e1 / (t2420 * t2437 - t2436 * t2421);
t2524 = t2358 * t2512;
t2359 = -t2410 * t2420 + t2414 * t2421;
t2523 = t2359 * t2512;
t2463 = legFrame(3,3);
t2449 = sin(t2463);
t2453 = cos(t2463);
t2411 = t2449 * g(1) - t2453 * g(2);
t2415 = t2453 * g(1) + t2449 * g(2);
t2459 = pkin(7) + qJ(2,3);
t2445 = qJ(3,3) + t2459;
t2426 = sin(t2445);
t2429 = cos(t2445);
t2360 = t2411 * t2429 + t2415 * t2426;
t2439 = sin(t2459);
t2442 = cos(t2459);
t2511 = 0.1e1 / (t2426 * t2442 - t2439 * t2429);
t2522 = t2360 * t2511;
t2464 = legFrame(2,3);
t2450 = sin(t2464);
t2454 = cos(t2464);
t2412 = t2450 * g(1) - t2454 * g(2);
t2416 = t2454 * g(1) + t2450 * g(2);
t2460 = pkin(7) + qJ(2,2);
t2446 = qJ(3,2) + t2460;
t2427 = sin(t2446);
t2430 = cos(t2446);
t2361 = t2412 * t2430 + t2416 * t2427;
t2440 = sin(t2460);
t2443 = cos(t2460);
t2510 = 0.1e1 / (t2427 * t2443 - t2440 * t2430);
t2521 = t2361 * t2510;
t2465 = legFrame(1,3);
t2451 = sin(t2465);
t2455 = cos(t2465);
t2413 = t2451 * g(1) - t2455 * g(2);
t2417 = t2455 * g(1) + t2451 * g(2);
t2461 = pkin(7) + qJ(2,1);
t2447 = qJ(3,1) + t2461;
t2428 = sin(t2447);
t2431 = cos(t2447);
t2362 = t2413 * t2431 + t2417 * t2428;
t2441 = sin(t2461);
t2444 = cos(t2461);
t2509 = 0.1e1 / (t2428 * t2444 - t2441 * t2431);
t2520 = t2362 * t2509;
t2363 = -t2411 * t2426 + t2415 * t2429;
t2519 = t2363 * t2511;
t2364 = -t2412 * t2427 + t2416 * t2430;
t2518 = t2364 * t2510;
t2365 = -t2413 * t2428 + t2417 * t2431;
t2517 = t2365 * t2509;
t2516 = koppelP(1,2);
t2515 = koppelP(2,2);
t2514 = koppelP(3,2);
t2513 = koppelP(4,2);
t2466 = xP(4);
t2456 = sin(t2466);
t2457 = cos(t2466);
t2467 = koppelP(4,1);
t2402 = t2456 * t2467 + t2457 * t2513;
t2406 = -t2456 * t2513 + t2457 * t2467;
t2366 = t2402 * t2452 - t2448 * t2406;
t2370 = t2448 * t2402 + t2406 * t2452;
t2480 = t2366 * t2421 - t2370 * t2420;
t2508 = t2480 * t2512;
t2468 = koppelP(3,1);
t2403 = t2456 * t2468 + t2457 * t2514;
t2407 = -t2456 * t2514 + t2457 * t2468;
t2367 = t2403 * t2453 - t2449 * t2407;
t2371 = t2449 * t2403 + t2407 * t2453;
t2479 = t2367 * t2429 - t2371 * t2426;
t2507 = t2479 * t2511;
t2469 = koppelP(2,1);
t2404 = t2456 * t2469 + t2457 * t2515;
t2408 = -t2456 * t2515 + t2457 * t2469;
t2368 = t2404 * t2454 - t2450 * t2408;
t2372 = t2450 * t2404 + t2408 * t2454;
t2478 = t2368 * t2430 - t2372 * t2427;
t2506 = t2478 * t2510;
t2470 = koppelP(1,1);
t2405 = t2456 * t2470 + t2457 * t2516;
t2409 = -t2456 * t2516 + t2457 * t2470;
t2369 = t2405 * t2455 - t2451 * t2409;
t2373 = t2451 * t2405 + t2409 * t2455;
t2477 = t2369 * t2431 - t2373 * t2428;
t2505 = t2477 * t2509;
t2394 = t2420 * t2448 - t2452 * t2421;
t2488 = t2512 * t2394;
t2476 = t2452 * t2420 + t2448 * t2421;
t2487 = t2512 * t2476;
t2396 = t2426 * t2449 - t2453 * t2429;
t2486 = t2511 * t2396;
t2475 = t2453 * t2426 + t2449 * t2429;
t2485 = t2511 * t2475;
t2397 = t2427 * t2450 - t2454 * t2430;
t2484 = t2510 * t2397;
t2474 = t2454 * t2427 + t2450 * t2430;
t2483 = t2510 * t2474;
t2398 = t2428 * t2451 - t2455 * t2431;
t2482 = t2509 * t2398;
t2473 = t2455 * t2428 + t2451 * t2431;
t2481 = t2509 * t2473;
t2472 = 0.1e1 / pkin(2);
t2471 = 0.1e1 / pkin(3);
t2419 = t2457 * g(1) + t2456 * g(2);
t2418 = t2456 * g(1) - t2457 * g(2);
t2381 = -t2413 * t2441 + t2417 * t2444;
t2380 = -t2412 * t2440 + t2416 * t2443;
t2379 = -t2411 * t2439 + t2415 * t2442;
t2378 = t2413 * t2444 + t2417 * t2441;
t2377 = t2412 * t2443 + t2416 * t2440;
t2376 = t2411 * t2442 + t2415 * t2439;
t2375 = -t2410 * t2436 + t2414 * t2437;
t2374 = t2410 * t2437 + t2414 * t2436;
t2357 = -pkin(2) * (t2441 * t2451 - t2455 * t2444) - t2398 * pkin(3);
t2356 = -pkin(2) * (t2440 * t2450 - t2454 * t2443) - t2397 * pkin(3);
t2355 = -pkin(2) * (t2439 * t2449 - t2453 * t2442) - t2396 * pkin(3);
t2354 = pkin(2) * (t2455 * t2441 + t2451 * t2444) + t2473 * pkin(3);
t2353 = pkin(2) * (t2454 * t2440 + t2450 * t2443) + t2474 * pkin(3);
t2352 = pkin(2) * (t2453 * t2439 + t2449 * t2442) + t2475 * pkin(3);
t2351 = -pkin(2) * (t2436 * t2448 - t2452 * t2437) - t2394 * pkin(3);
t2350 = pkin(2) * (t2452 * t2436 + t2448 * t2437) + t2476 * pkin(3);
t2345 = pkin(2) * (t2369 * t2444 - t2373 * t2441) + t2477 * pkin(3);
t2344 = pkin(2) * (t2368 * t2443 - t2372 * t2440) + t2478 * pkin(3);
t2343 = pkin(2) * (t2367 * t2442 - t2371 * t2439) + t2479 * pkin(3);
t2342 = pkin(2) * (t2366 * t2437 - t2370 * t2436) + t2480 * pkin(3);
t1 = [(-t2456 * t2418 - t2457 * t2419) * MDP(11) + ((-t2374 * t2488 - t2376 * t2486 - t2377 * t2484 - t2378 * t2482) * MDP(3) + (-t2375 * t2488 - t2379 * t2486 - t2380 * t2484 - t2381 * t2482) * MDP(4) + (-t2358 * t2488 - t2360 * t2486 - t2361 * t2484 - t2362 * t2482) * MDP(6) + (-t2359 * t2488 - t2363 * t2486 - t2364 * t2484 - t2365 * t2482) * MDP(7) + ((-t2351 * t2524 - t2355 * t2522 - t2356 * t2521 - t2357 * t2520) * MDP(6) + (-t2351 * t2523 - t2355 * t2519 - t2356 * t2518 - t2357 * t2517) * MDP(7)) * t2471) * t2472; (t2457 * t2418 - t2456 * t2419) * MDP(11) + ((t2374 * t2487 + t2376 * t2485 + t2377 * t2483 + t2378 * t2481) * MDP(3) + (t2375 * t2487 + t2379 * t2485 + t2380 * t2483 + t2381 * t2481) * MDP(4) + (t2358 * t2487 + t2360 * t2485 + t2361 * t2483 + t2362 * t2481) * MDP(6) + (t2359 * t2487 + t2363 * t2485 + t2364 * t2483 + t2365 * t2481) * MDP(7) + ((-t2350 * t2524 - t2352 * t2522 - t2353 * t2521 - t2354 * t2520) * MDP(6) + (-t2350 * t2523 - t2352 * t2519 - t2353 * t2518 - t2354 * t2517) * MDP(7)) * t2471) * t2472; (-(4 * MDP(1)) - MDP(11)) * g(3); t2419 * MDP(10) + t2418 * MDP(9) + ((-t2374 * t2508 - t2376 * t2507 - t2377 * t2506 - t2378 * t2505) * MDP(3) + (-t2375 * t2508 - t2379 * t2507 - t2380 * t2506 - t2381 * t2505) * MDP(4) + (-t2477 * t2520 - t2478 * t2521 - t2479 * t2522 - t2480 * t2524) * MDP(6) + (-t2477 * t2517 - t2478 * t2518 - t2479 * t2519 - t2480 * t2523) * MDP(7) + ((t2342 * t2524 + t2343 * t2522 + t2344 * t2521 + t2345 * t2520) * MDP(6) + (t2342 * t2523 + t2343 * t2519 + t2344 * t2518 + t2345 * t2517) * MDP(7)) * t2471) * t2472;];
taugX  = t1;
