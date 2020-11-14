% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRPRR8V1G2A0
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4,theta3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [13x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRPRR8V1G2A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:59
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR8V1G2A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1),zeros(13,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_mdp: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_mdp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:59:55
% EndTime: 2020-08-06 19:59:56
% DurationCPUTime: 0.93s
% Computational Cost: add. (702->145), mult. (1308->257), div. (117->9), fcn. (1218->23), ass. (0->133)
t2474 = MDP(3) - MDP(11);
t2390 = sin(qJ(1,3));
t2396 = cos(qJ(1,3));
t2386 = legFrame(3,2);
t2373 = sin(t2386);
t2376 = cos(t2386);
t2414 = t2376 * g(1) - t2373 * g(2);
t2345 = g(3) * t2390 - t2396 * t2414;
t2383 = pkin(4) + qJ(3,3);
t2379 = 0.1e1 / t2383;
t2462 = t2345 * t2379;
t2392 = sin(qJ(1,2));
t2398 = cos(qJ(1,2));
t2387 = legFrame(2,2);
t2374 = sin(t2387);
t2377 = cos(t2387);
t2413 = t2377 * g(1) - t2374 * g(2);
t2346 = g(3) * t2392 - t2398 * t2413;
t2384 = pkin(4) + qJ(3,2);
t2380 = 0.1e1 / t2384;
t2461 = t2346 * t2380;
t2394 = sin(qJ(1,1));
t2400 = cos(qJ(1,1));
t2388 = legFrame(1,2);
t2375 = sin(t2388);
t2378 = cos(t2388);
t2412 = t2378 * g(1) - t2375 * g(2);
t2347 = g(3) * t2394 - t2400 * t2412;
t2385 = pkin(4) + qJ(3,1);
t2381 = 0.1e1 / t2385;
t2460 = t2347 * t2381;
t2473 = MDP(12) * pkin(1);
t2472 = pkin(2) * cos(qJ(2,3) + pkin(5));
t2471 = pkin(2) * cos(qJ(2,2) + pkin(5));
t2470 = pkin(2) * cos(qJ(2,1) + pkin(5));
t2469 = pkin(2) * sin(pkin(5));
t2395 = cos(qJ(2,3));
t2465 = t2395 * pkin(1);
t2397 = cos(qJ(2,2));
t2464 = t2397 * pkin(1);
t2399 = cos(qJ(2,1));
t2463 = t2399 * pkin(1);
t2369 = cos(pkin(5)) * pkin(2) + pkin(1);
t2389 = sin(qJ(2,3));
t2411 = t2369 * t2395 - t2389 * t2469;
t2351 = 0.1e1 / t2411;
t2459 = t2351 * t2379;
t2391 = sin(qJ(2,2));
t2410 = t2369 * t2397 - t2391 * t2469;
t2352 = 0.1e1 / t2410;
t2458 = t2352 * t2380;
t2393 = sin(qJ(2,1));
t2409 = t2369 * t2399 - t2393 * t2469;
t2353 = 0.1e1 / t2409;
t2457 = t2353 * t2381;
t2363 = 0.1e1 / (t2465 + t2472);
t2456 = t2363 * t2373;
t2455 = t2363 * t2376;
t2364 = 0.1e1 / (t2464 + t2471);
t2454 = t2364 * t2374;
t2453 = t2364 * t2377;
t2365 = 0.1e1 / (t2463 + t2470);
t2452 = t2365 * t2375;
t2451 = t2365 * t2378;
t2450 = t2369 * t2376;
t2449 = t2369 * t2377;
t2448 = t2369 * t2378;
t2447 = t2373 * t2369;
t2446 = t2374 * t2369;
t2445 = t2375 * t2369;
t2444 = t2379 * t2396;
t2443 = t2380 * t2398;
t2442 = t2381 * t2400;
t2441 = t2376 * t2469;
t2440 = t2377 * t2469;
t2439 = t2378 * t2469;
t2438 = t2373 * t2469;
t2437 = t2374 * t2469;
t2436 = t2375 * t2469;
t2435 = t2395 * t2462;
t2434 = t2397 * t2461;
t2433 = t2399 * t2460;
t2432 = t2345 * t2459;
t2431 = t2346 * t2458;
t2430 = t2347 * t2457;
t2408 = g(3) * t2396 + t2390 * t2414;
t2429 = t2408 * t2459;
t2407 = g(3) * t2398 + t2392 * t2413;
t2428 = t2407 * t2458;
t2406 = g(3) * t2400 + t2394 * t2412;
t2427 = t2406 * t2457;
t2366 = t2396 * t2465;
t2333 = g(3) * (-t2396 * qJ(3,3) + t2390 * t2465) - t2414 * (t2390 * qJ(3,3) + t2366);
t2426 = t2333 * t2459;
t2367 = t2398 * t2464;
t2334 = g(3) * (-t2398 * qJ(3,2) + t2392 * t2464) - t2413 * (t2392 * qJ(3,2) + t2367);
t2425 = t2334 * t2458;
t2368 = t2400 * t2463;
t2335 = g(3) * (-t2400 * qJ(3,1) + t2394 * t2463) - t2412 * (t2394 * qJ(3,1) + t2368);
t2424 = t2335 * t2457;
t2423 = t2345 * t2444;
t2422 = t2346 * t2443;
t2421 = t2347 * t2442;
t2420 = t2351 * t2435;
t2419 = t2352 * t2434;
t2418 = t2353 * t2433;
t2417 = t2389 * t2432;
t2416 = t2391 * t2431;
t2415 = t2393 * t2430;
t2360 = t2373 * g(1) + t2376 * g(2);
t2327 = -t2360 * t2395 + t2389 * t2408;
t2361 = t2374 * g(1) + t2377 * g(2);
t2329 = -t2361 * t2397 + t2391 * t2407;
t2362 = t2375 * g(1) + t2378 * g(2);
t2331 = -t2362 * t2399 + t2393 * t2406;
t2405 = t2327 * t2456 + t2329 * t2454 + t2331 * t2452;
t2404 = t2327 * t2455 + t2329 * t2453 + t2331 * t2451;
t2356 = t2393 * t2369 + t2399 * t2469;
t2355 = t2391 * t2369 + t2397 * t2469;
t2354 = t2389 * t2369 + t2395 * t2469;
t2338 = -t2400 * t2385 + t2394 * t2409;
t2337 = -t2398 * t2384 + t2392 * t2410;
t2336 = -t2396 * t2383 + t2390 * t2411;
t2332 = t2362 * t2393 + t2399 * t2406;
t2330 = t2361 * t2391 + t2397 * t2407;
t2328 = t2360 * t2389 + t2395 * t2408;
t2326 = (t2394 * t2448 + t2436) * t2399 + (-t2394 * t2439 + t2445) * t2393;
t2325 = (t2392 * t2449 + t2437) * t2397 + (-t2392 * t2440 + t2446) * t2391;
t2324 = (t2390 * t2450 + t2438) * t2395 + (-t2390 * t2441 + t2447) * t2389;
t2323 = (-t2394 * t2445 + t2439) * t2399 + t2393 * (t2394 * t2436 + t2448);
t2322 = (-t2392 * t2446 + t2440) * t2397 + t2391 * (t2392 * t2437 + t2449);
t2321 = (-t2390 * t2447 + t2441) * t2395 + t2389 * (t2390 * t2438 + t2450);
t1 = [(t2324 * t2432 + t2325 * t2431 + t2326 * t2430) * MDP(2) + (t2324 * t2420 + t2325 * t2419 + t2326 * t2418 + t2405) * MDP(9) + (-t2324 * t2417 - t2325 * t2416 - t2326 * t2415 + t2328 * t2456 + t2330 * t2454 + t2332 * t2452) * MDP(10) + (t2326 * t2424 - (t2338 * t2378 + t2356 * t2375) * t2460 + t2325 * t2425 - (t2337 * t2377 + t2355 * t2374) * t2461 + t2324 * t2426 - (t2336 * t2376 + t2354 * t2373) * t2462) * MDP(12) - g(1) * MDP(13) + t2405 * t2473 + t2474 * (t2324 * t2429 + t2325 * t2428 + t2326 * t2427); (t2321 * t2432 + t2322 * t2431 + t2323 * t2430) * MDP(2) + (t2321 * t2420 + t2322 * t2419 + t2323 * t2418 + t2404) * MDP(9) + (-t2321 * t2417 - t2322 * t2416 - t2323 * t2415 + t2328 * t2455 + t2330 * t2453 + t2332 * t2451) * MDP(10) + (t2323 * t2424 - (-t2338 * t2375 + t2356 * t2378) * t2460 + t2322 * t2425 - (-t2337 * t2374 + t2355 * t2377) * t2461 + t2321 * t2426 - (-t2336 * t2373 + t2354 * t2376) * t2462) * MDP(12) - g(2) * MDP(13) + t2404 * t2473 + t2474 * (t2321 * t2429 + t2322 * t2428 + t2323 * t2427); (t2421 + t2422 + t2423) * MDP(2) + (t2396 * t2435 + t2398 * t2434 + t2400 * t2433) * MDP(9) + (-t2389 * t2423 - t2391 * t2422 - t2393 * t2421) * MDP(10) + ((t2400 * t2335 - (t2394 * t2385 + t2400 * t2470 + t2368) * t2347) * t2381 + (t2398 * t2334 - (t2392 * t2384 + t2398 * t2471 + t2367) * t2346) * t2380 + (t2396 * t2333 - (t2390 * t2383 + t2396 * t2472 + t2366) * t2345) * t2379) * MDP(12) - g(3) * MDP(13) + t2474 * (t2406 * t2442 + t2407 * t2443 + t2408 * t2444);];
taugX  = t1;
