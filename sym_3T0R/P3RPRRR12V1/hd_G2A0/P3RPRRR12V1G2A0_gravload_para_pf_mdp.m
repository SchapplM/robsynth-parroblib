% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RPRRR12V1G2A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [14x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPRRR12V1G2A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:25
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRRR12V1G2A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(14,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_mdp: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_mdp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:25:10
% EndTime: 2020-08-06 18:25:11
% DurationCPUTime: 1.04s
% Computational Cost: add. (423->132), mult. (771->247), div. (102->7), fcn. (735->18), ass. (0->107)
t2414 = MDP(2) - MDP(4);
t2413 = MDP(3) - MDP(5);
t2339 = sin(qJ(1,1));
t2345 = cos(qJ(1,1));
t2333 = legFrame(1,2);
t2323 = sin(t2333);
t2326 = cos(t2333);
t2356 = t2326 * g(1) - t2323 * g(2);
t2353 = g(3) * t2345 + t2356 * t2339;
t2338 = sin(qJ(3,1));
t2320 = t2338 * pkin(3) + qJ(2,1);
t2317 = 0.1e1 / t2320;
t2377 = t2345 * t2317;
t2371 = t2353 * t2377;
t2337 = sin(qJ(1,2));
t2343 = cos(qJ(1,2));
t2332 = legFrame(2,2);
t2322 = sin(t2332);
t2325 = cos(t2332);
t2357 = t2325 * g(1) - t2322 * g(2);
t2354 = g(3) * t2343 + t2357 * t2337;
t2336 = sin(qJ(3,2));
t2319 = t2336 * pkin(3) + qJ(2,2);
t2316 = 0.1e1 / t2319;
t2378 = t2343 * t2316;
t2373 = t2354 * t2378;
t2335 = sin(qJ(1,3));
t2341 = cos(qJ(1,3));
t2331 = legFrame(3,2);
t2321 = sin(t2331);
t2324 = cos(t2331);
t2358 = t2324 * g(1) - t2321 * g(2);
t2355 = g(3) * t2341 + t2358 * t2335;
t2334 = sin(qJ(3,3));
t2318 = t2334 * pkin(3) + qJ(2,3);
t2315 = 0.1e1 / t2318;
t2379 = t2341 * t2315;
t2375 = t2355 * t2379;
t2412 = g(3) * t2335 - t2358 * t2341;
t2411 = g(3) * t2337 - t2357 * t2343;
t2410 = g(3) * t2339 - t2356 * t2345;
t2409 = pkin(3) * t2341;
t2408 = pkin(3) * t2343;
t2407 = pkin(3) * t2345;
t2403 = qJ(2,2) * t2343;
t2402 = qJ(2,3) * t2341;
t2401 = t2345 * qJ(2,1);
t2328 = 0.1e1 / t2334;
t2400 = t2412 * t2328;
t2329 = 0.1e1 / t2336;
t2399 = t2411 * t2329;
t2330 = 0.1e1 / t2338;
t2398 = t2410 * t2330;
t2397 = t2315 * t2335;
t2396 = t2316 * t2337;
t2395 = t2317 * t2339;
t2394 = t2321 * t2328;
t2340 = cos(qJ(3,3));
t2393 = t2321 * t2340;
t2392 = t2322 * t2329;
t2342 = cos(qJ(3,2));
t2391 = t2322 * t2342;
t2390 = t2323 * t2330;
t2344 = cos(qJ(3,1));
t2389 = t2323 * t2344;
t2388 = t2324 * t2328;
t2387 = t2324 * t2340;
t2386 = t2325 * t2329;
t2385 = t2325 * t2342;
t2384 = t2326 * t2330;
t2383 = t2326 * t2344;
t2291 = -g(3) * (-t2335 * pkin(1) + t2402) - t2358 * (t2341 * pkin(1) + t2335 * qJ(2,3));
t2382 = t2335 * t2291;
t2292 = -g(3) * (-t2337 * pkin(1) + t2403) - t2357 * (t2343 * pkin(1) + t2337 * qJ(2,2));
t2381 = t2337 * t2292;
t2293 = -g(3) * (-t2339 * pkin(1) + t2401) - t2356 * (t2345 * pkin(1) + t2339 * qJ(2,1));
t2380 = t2339 * t2293;
t2376 = t2355 * t2397;
t2374 = t2354 * t2396;
t2372 = t2353 * t2395;
t2370 = t2321 * t2397;
t2369 = t2324 * t2397;
t2368 = t2322 * t2396;
t2367 = t2325 * t2396;
t2366 = t2323 * t2395;
t2365 = t2326 * t2395;
t2364 = t2334 * t2376;
t2363 = t2340 * t2376;
t2362 = t2336 * t2374;
t2361 = t2342 * t2374;
t2360 = t2338 * t2372;
t2359 = t2344 * t2372;
t2346 = 0.1e1 / pkin(3);
t2327 = pkin(1) + pkin(5) + pkin(6);
t2314 = -t2327 * t2337 + t2403;
t2313 = -t2327 * t2335 + t2402;
t2312 = -t2327 * t2339 + t2401;
t2311 = t2323 * g(1) + t2326 * g(2);
t2310 = t2322 * g(1) + t2325 * g(2);
t2309 = t2321 * g(1) + t2324 * g(2);
t2290 = t2311 * t2338 - t2344 * t2410;
t2289 = t2311 * t2344 + t2338 * t2410;
t2288 = t2310 * t2336 - t2342 * t2411;
t2287 = t2310 * t2342 + t2336 * t2411;
t2286 = t2309 * t2334 - t2340 * t2412;
t2285 = t2309 * t2340 + t2334 * t2412;
t1 = [((t2326 * t2380 - ((pkin(3) * t2389 - t2312 * t2326) * t2338 + (t2344 - 0.1e1) * (t2344 + 0.1e1) * t2326 * t2407 + qJ(2,1) * t2389) * t2398) * t2317 + (t2325 * t2381 - ((pkin(3) * t2391 - t2314 * t2325) * t2336 + (t2342 - 0.1e1) * (t2342 + 0.1e1) * t2325 * t2408 + qJ(2,2) * t2391) * t2399) * t2316 + (t2324 * t2382 - ((pkin(3) * t2393 - t2313 * t2324) * t2334 + (t2340 - 0.1e1) * (t2340 + 0.1e1) * t2324 * t2409 + qJ(2,3) * t2393) * t2400) * t2315) * MDP(6) + (-t2324 * t2364 - t2325 * t2362 - t2326 * t2360 + (-t2286 * t2394 - t2288 * t2392 - t2290 * t2390) * t2346) * MDP(12) + (-t2324 * t2363 - t2325 * t2361 - t2326 * t2359 + (-t2285 * t2394 - t2287 * t2392 - t2289 * t2390) * t2346) * MDP(13) - g(1) * MDP(14) + t2414 * (t2365 * t2410 + t2367 * t2411 + t2369 * t2412) + t2413 * (t2353 * t2365 + t2354 * t2367 + t2355 * t2369); ((-t2323 * t2380 - ((pkin(3) * t2383 + t2312 * t2323) * t2338 + (-t2344 ^ 2 + 0.1e1) * t2323 * t2407 + qJ(2,1) * t2383) * t2398) * t2317 + (-t2322 * t2381 - ((pkin(3) * t2385 + t2314 * t2322) * t2336 + (-t2342 ^ 2 + 0.1e1) * t2322 * t2408 + qJ(2,2) * t2385) * t2399) * t2316 + (-t2321 * t2382 - ((pkin(3) * t2387 + t2313 * t2321) * t2334 + (-t2340 ^ 2 + 0.1e1) * t2321 * t2409 + qJ(2,3) * t2387) * t2400) * t2315) * MDP(6) + (t2321 * t2364 + t2322 * t2362 + t2323 * t2360 + (-t2286 * t2388 - t2288 * t2386 - t2290 * t2384) * t2346) * MDP(12) + (t2321 * t2363 + t2322 * t2361 + t2323 * t2359 + (-t2285 * t2388 - t2287 * t2386 - t2289 * t2384) * t2346) * MDP(13) - g(2) * MDP(14) - t2414 * (t2366 * t2410 + t2368 * t2411 + t2370 * t2412) - t2413 * (t2353 * t2366 + t2354 * t2368 + t2355 * t2370); ((t2345 * t2293 - (t2339 * t2320 + t2327 * t2345) * t2410) * t2317 + (t2343 * t2292 - (t2337 * t2319 + t2327 * t2343) * t2411) * t2316 + (t2341 * t2291 - (t2335 * t2318 + t2327 * t2341) * t2412) * t2315) * MDP(6) + (-t2334 * t2375 - t2336 * t2373 - t2338 * t2371) * MDP(12) + (-t2340 * t2375 - t2342 * t2373 - t2344 * t2371) * MDP(13) - g(3) * MDP(14) + t2414 * (t2377 * t2410 + t2378 * t2411 + t2379 * t2412) + t2413 * (t2373 + t2375 + t2371);];
taugX  = t1;
