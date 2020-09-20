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
% Datum: 2020-08-07 10:47
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRRRR7V2G3A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 10:21:12
% EndTime: 2020-08-07 10:21:15
% DurationCPUTime: 2.94s
% Computational Cost: add. (1191->320), mult. (2187->442), div. (66->14), fcn. (1239->66), ass. (0->209)
t2471 = pkin(6) + pkin(5);
t2374 = (-rSges(2,3) - pkin(5)) * m(2) + m(1) * rSges(1,2) - (rSges(3,3) + t2471) * m(3);
t2387 = m(1) * rSges(1,1) + (m(2) + m(3)) * pkin(1);
t2334 = sin(qJ(3,1));
t2343 = cos(qJ(3,1));
t2398 = m(2) * rSges(2,1) + pkin(2) * m(3);
t2243 = (rSges(3,1) * t2343 - rSges(3,2) * t2334) * m(3) + t2398;
t2346 = m(2) * rSges(2,2);
t2249 = t2346 + (rSges(3,1) * t2334 + rSges(3,2) * t2343) * m(3);
t2335 = sin(qJ(2,1));
t2344 = cos(qJ(2,1));
t2479 = t2243 * t2344 - t2249 * t2335;
t2331 = sin(qJ(3,2));
t2340 = cos(qJ(3,2));
t2242 = (rSges(3,1) * t2340 - rSges(3,2) * t2331) * m(3) + t2398;
t2248 = t2346 + (rSges(3,1) * t2331 + rSges(3,2) * t2340) * m(3);
t2332 = sin(qJ(2,2));
t2341 = cos(qJ(2,2));
t2478 = t2242 * t2341 - t2248 * t2332;
t2328 = sin(qJ(3,3));
t2337 = cos(qJ(3,3));
t2241 = (rSges(3,1) * t2337 - rSges(3,2) * t2328) * m(3) + t2398;
t2247 = t2346 + (rSges(3,1) * t2328 + rSges(3,2) * t2337) * m(3);
t2329 = sin(qJ(2,3));
t2338 = cos(qJ(2,3));
t2477 = t2241 * t2338 - t2247 * t2329;
t2476 = 2 * pkin(3);
t2310 = (pkin(7) + t2471);
t2475 = 2 * t2310;
t2474 = 0.2e1 * t2338 ^ 2;
t2473 = 0.2e1 * t2341 ^ 2;
t2472 = 0.2e1 * t2344 ^ 2;
t2470 = m(3) / pkin(3);
t2314 = t2337 ^ 2;
t2302 = pkin(3) * t2314;
t2316 = t2340 ^ 2;
t2303 = pkin(3) * t2316;
t2318 = t2343 ^ 2;
t2304 = pkin(3) * t2318;
t2469 = t2328 * pkin(1);
t2468 = t2328 * pkin(3);
t2467 = t2331 * pkin(1);
t2466 = t2331 * pkin(3);
t2465 = t2334 * pkin(1);
t2464 = t2334 * pkin(3);
t2463 = t2337 * pkin(2);
t2299 = t2337 * pkin(3);
t2462 = t2340 * pkin(2);
t2300 = t2340 * pkin(3);
t2461 = t2343 * pkin(2);
t2301 = t2343 * pkin(3);
t2460 = -qJ(3,1) + qJ(1,1);
t2459 = qJ(3,1) + qJ(1,1);
t2458 = -qJ(3,2) + qJ(1,2);
t2457 = qJ(3,2) + qJ(1,2);
t2456 = -qJ(3,3) + qJ(1,3);
t2455 = qJ(3,3) + qJ(1,3);
t2454 = qJ(1,1) + 0.2e1 * qJ(2,1);
t2453 = qJ(1,1) - 0.2e1 * qJ(2,1);
t2452 = qJ(1,2) + 0.2e1 * qJ(2,2);
t2451 = qJ(1,2) - 0.2e1 * qJ(2,2);
t2450 = 0.2e1 * qJ(2,3) + qJ(1,3);
t2449 = qJ(1,3) - 0.2e1 * qJ(2,3);
t2448 = 0.2e1 * pkin(1);
t2339 = cos(qJ(1,3));
t2330 = sin(qJ(1,3));
t2402 = pkin(1) * t2339 + t2330 * t2310;
t2408 = t2328 * t2329;
t2226 = t2402 * t2408 + (t2314 - 0.1e1) * t2339 * pkin(3);
t2325 = legFrame(3,2);
t2293 = sin(t2325);
t2447 = t2226 * t2293;
t2296 = cos(t2325);
t2446 = t2226 * t2296;
t2342 = cos(qJ(1,2));
t2333 = sin(qJ(1,2));
t2401 = pkin(1) * t2342 + t2333 * t2310;
t2407 = t2331 * t2332;
t2227 = t2401 * t2407 + (t2316 - 0.1e1) * t2342 * pkin(3);
t2326 = legFrame(2,2);
t2294 = sin(t2326);
t2445 = t2227 * t2294;
t2297 = cos(t2326);
t2444 = t2227 * t2297;
t2345 = cos(qJ(1,1));
t2336 = sin(qJ(1,1));
t2400 = pkin(1) * t2345 + t2336 * t2310;
t2406 = t2334 * t2335;
t2228 = t2400 * t2406 + (t2318 - 0.1e1) * t2345 * pkin(3);
t2327 = legFrame(1,2);
t2295 = sin(t2327);
t2443 = t2228 * t2295;
t2298 = cos(t2327);
t2442 = t2228 * t2298;
t2390 = pkin(3) * t2408;
t2277 = t2299 + pkin(2);
t2425 = t2277 * t2338;
t2441 = 0.1e1 / (pkin(1) - t2390 + t2425) / t2328;
t2389 = pkin(3) * t2407;
t2278 = t2300 + pkin(2);
t2423 = t2278 * t2341;
t2440 = 0.1e1 / (pkin(1) - t2389 + t2423) / t2331;
t2388 = pkin(3) * t2406;
t2279 = t2301 + pkin(2);
t2421 = t2279 * t2344;
t2439 = 0.1e1 / (pkin(1) - t2388 + t2421) / t2334;
t2256 = t2296 * g(1) - t2293 * g(2);
t2320 = qJ(2,3) + qJ(3,3);
t2284 = cos(t2320);
t2372 = t2374 * g(3);
t2373 = t2387 * g(3);
t2432 = 0.1e1 / (t2338 * pkin(2) + pkin(3) * t2284 + pkin(1)) * ((t2477 * g(3) + t2373) * t2339 - t2330 * t2372 + (t2374 * t2339 + t2330 * (t2387 + t2477)) * t2256);
t2257 = t2297 * g(1) - t2294 * g(2);
t2321 = qJ(2,2) + qJ(3,2);
t2285 = cos(t2321);
t2431 = 0.1e1 / (t2341 * pkin(2) + pkin(3) * t2285 + pkin(1)) * ((t2478 * g(3) + t2373) * t2342 - t2333 * t2372 + (t2374 * t2342 + t2333 * (t2387 + t2478)) * t2257);
t2258 = t2298 * g(1) - t2295 * g(2);
t2322 = qJ(2,1) + qJ(3,1);
t2286 = cos(t2322);
t2430 = 0.1e1 / (t2344 * pkin(2) + pkin(3) * t2286 + pkin(1)) * ((t2479 * g(3) + t2373) * t2345 - t2336 * t2372 + (t2374 * t2345 + t2336 * (t2387 + t2479)) * t2258);
t2352 = pkin(2) / 0.2e1;
t2429 = (t2299 + t2352) * t2328;
t2428 = (t2300 + t2352) * t2331;
t2427 = (t2301 + t2352) * t2334;
t2426 = t2277 * t2293;
t2424 = t2278 * t2294;
t2422 = t2279 * t2295;
t2420 = t2293 * t2339;
t2419 = t2294 * t2342;
t2418 = t2295 * t2345;
t2417 = t2296 * t2277;
t2416 = t2296 * t2339;
t2415 = t2297 * t2278;
t2414 = t2297 * t2342;
t2413 = t2298 * t2279;
t2412 = t2298 * t2345;
t2362 = pkin(3) ^ 2;
t2411 = t2314 * t2362;
t2410 = t2316 * t2362;
t2409 = t2318 * t2362;
t2265 = pkin(1) * t2329 - t2468;
t2405 = t2337 * t2265;
t2266 = pkin(1) * t2332 - t2466;
t2404 = t2340 * t2266;
t2267 = pkin(1) * t2335 - t2464;
t2403 = t2343 * t2267;
t2364 = pkin(2) ^ 2;
t2399 = -t2362 / 0.2e1 + t2364 / 0.2e1;
t2397 = pkin(2) * t2299;
t2396 = pkin(2) * t2300;
t2395 = pkin(2) * t2301;
t2393 = t2277 * t2468;
t2392 = t2278 * t2466;
t2391 = t2279 * t2464;
t2253 = t2293 * g(1) + t2296 * g(2);
t2371 = -t2330 * g(3) + t2256 * t2339;
t2220 = (t2371 * t2241 + t2253 * t2247) * t2329 + (-t2253 * t2241 + t2371 * t2247) * t2338;
t2386 = t2220 * t2441;
t2254 = t2294 * g(1) + t2297 * g(2);
t2370 = -t2333 * g(3) + t2257 * t2342;
t2221 = (t2370 * t2242 + t2254 * t2248) * t2332 + (-t2254 * t2242 + t2370 * t2248) * t2341;
t2385 = t2221 * t2440;
t2255 = t2295 * g(1) + t2298 * g(2);
t2369 = -t2336 * g(3) + t2258 * t2345;
t2222 = (t2369 * t2243 + t2255 * t2249) * t2335 + (-t2255 * t2243 + t2369 * t2249) * t2344;
t2384 = t2222 * t2439;
t2281 = sin(t2320);
t2383 = ((-rSges(3,1) * t2253 + t2371 * rSges(3,2)) * t2284 + t2281 * (t2371 * rSges(3,1) + rSges(3,2) * t2253)) * t2441;
t2282 = sin(t2321);
t2382 = ((-rSges(3,1) * t2254 + t2370 * rSges(3,2)) * t2285 + t2282 * (t2370 * rSges(3,1) + rSges(3,2) * t2254)) * t2440;
t2283 = sin(t2322);
t2381 = ((-rSges(3,1) * t2255 + t2369 * rSges(3,2)) * t2286 + t2283 * (t2369 * rSges(3,1) + rSges(3,2) * t2255)) * t2439;
t2380 = t2330 * t2432;
t2379 = t2333 * t2431;
t2378 = t2336 * t2430;
t2377 = t2339 * t2408;
t2376 = t2342 * t2407;
t2375 = t2345 * t2406;
t2235 = t2377 * t2476 - t2402;
t2368 = pkin(2) * t2377 + t2235 * t2337;
t2236 = t2376 * t2476 - t2401;
t2367 = pkin(2) * t2376 + t2236 * t2340;
t2237 = t2375 * t2476 - t2400;
t2366 = pkin(2) * t2375 + t2237 * t2343;
t2365 = 0.1e1 / pkin(2);
t2359 = 0.2e1 * qJ(3,1);
t2356 = 0.2e1 * qJ(3,2);
t2353 = 0.2e1 * qJ(3,3);
t2351 = -pkin(3) / 0.2e1;
t2308 = -t2362 + t2364;
t2292 = -qJ(2,1) + t2460;
t2291 = qJ(2,1) + t2459;
t2290 = -qJ(2,2) + t2458;
t2289 = qJ(2,2) + t2457;
t2288 = -qJ(2,3) + t2456;
t2287 = qJ(2,3) + t2455;
t2261 = t2304 + t2461 / 0.2e1 + t2351;
t2260 = t2303 + t2462 / 0.2e1 + t2351;
t2259 = t2302 + t2463 / 0.2e1 + t2351;
t2246 = t2395 + t2399 + t2409;
t2245 = t2396 + t2399 + t2410;
t2244 = t2397 + t2399 + t2411;
t2234 = t2465 + (-pkin(3) + t2461 + 0.2e1 * t2304) * t2335;
t2233 = t2467 + (-pkin(3) + t2462 + 0.2e1 * t2303) * t2332;
t2232 = t2469 + (-pkin(3) + t2463 + 0.2e1 * t2302) * t2329;
t2231 = pkin(1) * t2464 + (t2308 + 0.2e1 * t2395 + 0.2e1 * t2409) * t2335;
t2230 = pkin(1) * t2466 + (t2308 + 0.2e1 * t2396 + 0.2e1 * t2410) * t2332;
t2229 = pkin(1) * t2468 + (t2308 + 0.2e1 * t2397 + 0.2e1 * t2411) * t2329;
t1 = [-t2296 * t2380 - t2297 * t2379 - t2298 * t2378 - m(4) * g(1) + (((t2261 * t2412 + t2295 * t2427) * t2472 + (t2295 * t2234 - t2366 * t2298) * t2344 - t2442 + t2295 * t2403) * t2384 + ((t2260 * t2414 + t2294 * t2428) * t2473 + (t2294 * t2233 - t2367 * t2297) * t2341 - t2444 + t2294 * t2404) * t2385 + ((t2259 * t2416 + t2293 * t2429) * t2474 + (t2293 * t2232 - t2368 * t2296) * t2338 - t2446 + t2293 * t2405) * t2386 + (((-t2246 * t2412 - t2295 * t2391) * t2472 + (-t2295 * t2231 + t2237 * t2413) * t2344 + pkin(3) * t2442 - t2267 * t2422) * t2381 + ((-t2245 * t2414 - t2294 * t2392) * t2473 + (-t2294 * t2230 + t2236 * t2415) * t2341 + pkin(3) * t2444 - t2266 * t2424) * t2382 + ((-t2244 * t2416 - t2293 * t2393) * t2474 + (-t2293 * t2229 + t2235 * t2417) * t2338 + pkin(3) * t2446 - t2265 * t2426) * t2383) * t2470) * t2365; t2293 * t2380 + t2294 * t2379 + t2295 * t2378 - m(4) * g(2) + (((-t2261 * t2418 + t2298 * t2427) * t2472 + (t2298 * t2234 + t2366 * t2295) * t2344 + t2443 + t2298 * t2403) * t2384 + ((-t2260 * t2419 + t2297 * t2428) * t2473 + (t2297 * t2233 + t2367 * t2294) * t2341 + t2445 + t2297 * t2404) * t2385 + ((-t2259 * t2420 + t2296 * t2429) * t2474 + (t2296 * t2232 + t2368 * t2293) * t2338 + t2447 + t2296 * t2405) * t2386 + (((t2246 * t2418 - t2298 * t2391) * t2472 + (-t2298 * t2231 - t2237 * t2422) * t2344 - pkin(3) * t2443 - t2267 * t2413) * t2381 + ((t2245 * t2419 - t2297 * t2392) * t2473 + (-t2297 * t2230 - t2236 * t2424) * t2341 - pkin(3) * t2445 - t2266 * t2415) * t2382 + ((t2244 * t2420 - t2296 * t2393) * t2474 + (-t2296 * t2229 - t2235 * t2426) * t2338 - pkin(3) * t2447 - t2265 * t2417) * t2383) * t2470) * t2365; -m(4) * g(3) - t2339 * t2432 - t2342 * t2431 - t2345 * t2430 + ((t2246 * t2336 * t2472 + ((pkin(1) - 0.2e1 * t2388) * t2336 - t2345 * t2310) * t2421 - pkin(3) * ((pkin(1) * t2406 - pkin(3) + t2304) * t2336 - t2310 * t2375)) * t2381 + (t2245 * t2333 * t2473 + ((pkin(1) - 0.2e1 * t2389) * t2333 - t2342 * t2310) * t2423 - pkin(3) * ((pkin(1) * t2407 - pkin(3) + t2303) * t2333 - t2310 * t2376)) * t2382 + (t2244 * t2330 * t2474 + ((pkin(1) - 0.2e1 * t2390) * t2330 - t2339 * t2310) * t2425 - pkin(3) * ((pkin(1) * t2408 - pkin(3) + t2302) * t2330 - t2310 * t2377)) * t2383) * t2365 * t2470 + ((-sin(t2292) - sin(t2291)) * t2448 + (cos(t2292) + cos(t2291)) * t2475 + (-sin(-0.2e1 * qJ(3,1) + t2453) - sin(t2359 + t2454) - 0.2e1 * t2336) * pkin(3) + (-sin(-qJ(3,1) + t2453) - sin(qJ(3,1) + t2454) - sin(t2460) - sin(t2459)) * pkin(2)) / (-t2364 * sin(qJ(2,1) - qJ(3,1)) + pkin(2) * (0.2e1 * t2465 + pkin(2) * t2283 + (sin(t2359 + qJ(2,1)) - t2335) * pkin(3))) * t2222 / 0.2e1 + ((-sin(t2290) - sin(t2289)) * t2448 + (cos(t2290) + cos(t2289)) * t2475 + (-sin(-0.2e1 * qJ(3,2) + t2451) - sin(t2356 + t2452) - 0.2e1 * t2333) * pkin(3) + (-sin(-qJ(3,2) + t2451) - sin(qJ(3,2) + t2452) - sin(t2458) - sin(t2457)) * pkin(2)) / (-t2364 * sin(qJ(2,2) - qJ(3,2)) + pkin(2) * (0.2e1 * t2467 + pkin(2) * t2282 + (sin(t2356 + qJ(2,2)) - t2332) * pkin(3))) * t2221 / 0.2e1 + ((-sin(t2288) - sin(t2287)) * t2448 + (cos(t2288) + cos(t2287)) * t2475 + (-sin(-0.2e1 * qJ(3,3) + t2449) - sin(t2353 + t2450) - 0.2e1 * t2330) * pkin(3) + (-sin(-qJ(3,3) + t2449) - sin(qJ(3,3) + t2450) - sin(t2456) - sin(t2455)) * pkin(2)) / (-t2364 * sin(qJ(2,3) - qJ(3,3)) + pkin(2) * (0.2e1 * t2469 + pkin(2) * t2281 + (sin(t2353 + qJ(2,3)) - t2329) * pkin(3))) * t2220 / 0.2e1;];
taugX  = t1;
