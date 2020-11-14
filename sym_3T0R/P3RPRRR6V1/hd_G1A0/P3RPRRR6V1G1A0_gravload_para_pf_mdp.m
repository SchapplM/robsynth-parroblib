% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RPRRR6V1G1A0
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
%   see P3RPRRR6V1G1A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:32
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRRR6V1G1A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:32:39
% EndTime: 2020-08-06 18:32:40
% DurationCPUTime: 1.24s
% Computational Cost: add. (841->173), mult. (783->213), div. (60->7), fcn. (576->82), ass. (0->132)
t2450 = MDP(4) * pkin(1) + MDP(2);
t2446 = 2 * pkin(2);
t2385 = sin(qJ(3,3));
t2399 = 0.2e1 * qJ(3,3);
t2312 = 0.1e1 / (pkin(3) * sin(t2399) + t2385 * t2446 + (sin(pkin(7) + qJ(3,3)) + sin(-pkin(7) + qJ(3,3))) * pkin(1));
t2387 = sin(qJ(3,2));
t2400 = 0.2e1 * qJ(3,2);
t2313 = 0.1e1 / (pkin(3) * sin(t2400) + t2387 * t2446 + (sin(pkin(7) + qJ(3,2)) + sin(-pkin(7) + qJ(3,2))) * pkin(1));
t2389 = sin(qJ(3,1));
t2401 = 0.2e1 * qJ(3,1);
t2314 = 0.1e1 / (pkin(3) * sin(t2401) + t2389 * t2446 + (sin(pkin(7) + qJ(3,1)) + sin(-pkin(7) + qJ(3,1))) * pkin(1));
t2381 = qJ(1,1) + pkin(7);
t2384 = legFrame(1,3);
t2356 = t2384 + t2381;
t2351 = qJ(3,1) + t2356;
t2352 = -qJ(3,1) + t2356;
t2424 = sin(t2351) + sin(t2352);
t2380 = qJ(1,2) + pkin(7);
t2383 = legFrame(2,3);
t2355 = t2383 + t2380;
t2347 = qJ(3,2) + t2355;
t2348 = -qJ(3,2) + t2355;
t2425 = sin(t2347) + sin(t2348);
t2379 = qJ(1,3) + pkin(7);
t2382 = legFrame(3,3);
t2354 = t2382 + t2379;
t2343 = qJ(3,3) + t2354;
t2344 = -qJ(3,3) + t2354;
t2426 = sin(t2343) + sin(t2344);
t2449 = -t2426 * t2312 - t2425 * t2313 - t2424 * t2314;
t2448 = -0.2e1 * pkin(1);
t2447 = -2 * pkin(2);
t2398 = (pkin(6) + pkin(5));
t2445 = -2 * t2398;
t2444 = 2 * t2398;
t2442 = cos(pkin(7)) * pkin(1) + pkin(2);
t2371 = qJ(1,1) + t2384;
t2370 = qJ(1,2) + t2383;
t2369 = qJ(1,3) + t2382;
t2391 = cos(qJ(3,3));
t2373 = sin(t2382);
t2376 = cos(t2382);
t2315 = t2373 * g(1) - t2376 * g(2);
t2318 = t2376 * g(1) + t2373 * g(2);
t2363 = sin(t2379);
t2366 = cos(t2379);
t2407 = -t2315 * t2363 + t2318 * t2366;
t2441 = (-g(3) * t2391 + t2407 * t2385) * t2312;
t2440 = (g(3) * t2385 + t2407 * t2391) * t2312;
t2393 = cos(qJ(3,2));
t2374 = sin(t2383);
t2377 = cos(t2383);
t2316 = t2374 * g(1) - t2377 * g(2);
t2319 = t2377 * g(1) + t2374 * g(2);
t2364 = sin(t2380);
t2367 = cos(t2380);
t2406 = -t2316 * t2364 + t2319 * t2367;
t2439 = (-g(3) * t2393 + t2406 * t2387) * t2313;
t2438 = (g(3) * t2387 + t2406 * t2393) * t2313;
t2395 = cos(qJ(3,1));
t2375 = sin(t2384);
t2378 = cos(t2384);
t2317 = t2375 * g(1) - t2378 * g(2);
t2320 = t2378 * g(1) + t2375 * g(2);
t2365 = sin(t2381);
t2368 = cos(t2381);
t2405 = -t2317 * t2365 + t2320 * t2368;
t2437 = (-g(3) * t2395 + t2405 * t2389) * t2314;
t2436 = (g(3) * t2389 + t2405 * t2395) * t2314;
t2321 = 0.1e1 / (t2391 * pkin(3) + t2442);
t2435 = (t2315 * t2366 + t2363 * t2318) * t2321;
t2322 = 0.1e1 / (t2393 * pkin(3) + t2442);
t2434 = (t2316 * t2367 + t2364 * t2319) * t2322;
t2323 = 0.1e1 / (t2395 * pkin(3) + t2442);
t2433 = (t2317 * t2368 + t2365 * t2320) * t2323;
t2336 = sin(t2354);
t2432 = t2321 * t2336;
t2339 = cos(t2354);
t2431 = t2321 * t2339;
t2337 = sin(t2355);
t2430 = t2322 * t2337;
t2340 = cos(t2355);
t2429 = t2322 * t2340;
t2338 = sin(t2356);
t2428 = t2323 * t2338;
t2341 = cos(t2356);
t2427 = t2323 * t2341;
t2423 = cos(t2343) + cos(t2344);
t2422 = cos(t2347) + cos(t2348);
t2421 = cos(t2351) + cos(t2352);
t2420 = 0.2e1 * t2312;
t2419 = 0.2e1 * t2313;
t2418 = 0.2e1 * t2314;
t2417 = MDP(4) * g(3) / 0.2e1;
t2416 = t2385 * t2435;
t2415 = t2391 * t2435;
t2414 = t2387 * t2434;
t2413 = t2393 * t2434;
t2412 = t2389 * t2433;
t2411 = t2395 * t2433;
t2402 = 0.1e1 / pkin(3);
t2396 = cos(qJ(1,1));
t2394 = cos(qJ(1,2));
t2392 = cos(qJ(1,3));
t2390 = sin(qJ(1,1));
t2388 = sin(qJ(1,2));
t2386 = sin(qJ(1,3));
t2362 = -qJ(3,1) + t2371;
t2361 = qJ(3,1) + t2371;
t2360 = -qJ(3,2) + t2370;
t2359 = qJ(3,2) + t2370;
t2358 = -qJ(3,3) + t2369;
t2357 = qJ(3,3) + t2369;
t2353 = -0.2e1 * qJ(3,1) + t2356;
t2350 = t2401 + t2356;
t2349 = -0.2e1 * qJ(3,2) + t2355;
t2346 = t2400 + t2355;
t2345 = -0.2e1 * qJ(3,3) + t2354;
t2342 = t2399 + t2354;
t2311 = -t2317 * t2390 + t2320 * t2396;
t2310 = -t2316 * t2388 + t2319 * t2394;
t2309 = -t2315 * t2386 + t2318 * t2392;
t2308 = t2317 * t2396 + t2320 * t2390;
t2307 = t2316 * t2394 + t2319 * t2388;
t2306 = t2315 * t2392 + t2318 * t2386;
t2296 = t2338 * t2445 + cos(t2371) * t2448 + t2341 * t2447 - t2421 * pkin(3);
t2295 = t2337 * t2445 + cos(t2370) * t2448 + t2340 * t2447 - t2422 * pkin(3);
t2294 = t2336 * t2445 + cos(t2369) * t2448 + t2339 * t2447 - t2423 * pkin(3);
t2293 = t2341 * t2444 + sin(t2371) * t2448 + t2338 * t2447 - t2424 * pkin(3);
t2292 = t2340 * t2444 + sin(t2370) * t2448 + t2337 * t2447 - t2425 * pkin(3);
t2291 = t2339 * t2444 + sin(t2369) * t2448 + t2336 * t2447 - t2426 * pkin(3);
t1 = [(-t2309 * t2432 - t2310 * t2430 - t2311 * t2428) * MDP(3) + (-t2336 * t2415 - t2337 * t2413 - t2338 * t2411) * MDP(10) + (t2336 * t2416 + t2337 * t2414 + t2338 * t2412) * MDP(11) - g(1) * MDP(12) + ((-t2421 * t2418 - t2422 * t2419 - t2423 * t2420) * pkin(2) + (-(cos(t2362) + cos(t2361)) * t2418 - (cos(t2360) + cos(t2359)) * t2419 - (cos(t2358) + cos(t2357)) * t2420) * pkin(1) + t2449 * t2444 + (-(cos(t2353) + cos(t2350) + 0.2e1 * t2341) * t2314 - (cos(t2349) + cos(t2346) + 0.2e1 * t2340) * t2313 - (cos(t2345) + cos(t2342) + 0.2e1 * t2339) * t2312) * pkin(3)) * t2417 + ((t2294 * t2441 + t2295 * t2439 + t2296 * t2437) * MDP(10) + (t2294 * t2440 + t2295 * t2438 + t2296 * t2436) * MDP(11)) * t2402 + t2450 * (-t2306 * t2432 - t2307 * t2430 - t2308 * t2428); (t2309 * t2431 + t2310 * t2429 + t2311 * t2427) * MDP(3) + (t2339 * t2415 + t2340 * t2413 + t2341 * t2411) * MDP(10) + (-t2339 * t2416 - t2340 * t2414 - t2341 * t2412) * MDP(11) - g(2) * MDP(12) + ((-(sin(t2362) + sin(t2361)) * t2418 - (sin(t2360) + sin(t2359)) * t2419 - (sin(t2358) + sin(t2357)) * t2420) * pkin(1) + (t2423 * t2312 + t2422 * t2313 + t2421 * t2314) * t2444 + (-(sin(t2353) + sin(t2350) + 0.2e1 * t2338) * t2314 - (sin(t2349) + sin(t2346) + 0.2e1 * t2337) * t2313 - (sin(t2345) + sin(t2342) + 0.2e1 * t2336) * t2312) * pkin(3) + t2449 * t2446) * t2417 + ((t2291 * t2441 + t2292 * t2439 + t2293 * t2437) * MDP(10) + (t2291 * t2440 + t2292 * t2438 + t2293 * t2436) * MDP(11)) * t2402 + t2450 * (t2306 * t2431 + t2307 * t2429 + t2308 * t2427); (-MDP(12) - 0.3e1 * MDP(4)) * g(3);];
taugX  = t1;
