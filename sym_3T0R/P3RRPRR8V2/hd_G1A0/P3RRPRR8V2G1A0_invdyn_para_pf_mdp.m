% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RRPRR8V2G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [15x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRPRR8V2G1A0_convert_par2_MPV_fixb.m

% Output:
% tauX [3x1]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2022-11-07 13:12
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRPRR8V2G1A0_invdyn_para_pf_mdp(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_mdp: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-07 13:12:14
% EndTime: 2022-11-07 13:12:29
% DurationCPUTime: 16.58s
% Computational Cost: add. (214803->811), mult. (252940->1144), div. (21354->15), fcn. (112311->140), ass. (0->544)
t2318 = (pkin(3) ^ 2);
t2320 = (pkin(2) ^ 2);
t2630 = 3 * t2320 + 6 * t2318;
t2629 = 6 * t2320 + 3 * t2318;
t2293 = cos(qJ(2,3));
t2241 = t2293 * pkin(2);
t2256 = qJ(2,3) + pkin(7);
t2217 = cos(t2256);
t2551 = pkin(3) * t2217;
t2615 = t2241 + t2551;
t2295 = cos(qJ(2,2));
t2242 = t2295 * pkin(2);
t2259 = qJ(2,2) + pkin(7);
t2220 = cos(t2259);
t2550 = pkin(3) * t2220;
t2614 = t2242 + t2550;
t2297 = cos(qJ(2,1));
t2243 = t2297 * pkin(2);
t2262 = qJ(2,1) + pkin(7);
t2223 = cos(t2262);
t2549 = pkin(3) * t2223;
t2613 = t2243 + t2549;
t2277 = legFrame(1,3);
t2229 = sin(t2277);
t2232 = cos(t2277);
t2292 = sin(qJ(1,1));
t2298 = cos(qJ(1,1));
t2080 = t2298 * t2229 + t2292 * t2232;
t2081 = -t2292 * t2229 + t2298 * t2232;
t2274 = cos(pkin(7));
t2547 = pkin(3) * t2274;
t2195 = pkin(2) + t2547;
t2280 = (qJ(3,1) + pkin(5));
t2250 = -pkin(6) - t2280;
t2237 = 0.1e1 / t2250;
t2291 = sin(qJ(2,1));
t2299 = xDP(3);
t2300 = xDP(2);
t2301 = xDP(1);
t2273 = sin(pkin(7));
t2548 = pkin(3) * t2273;
t2039 = (t2081 * t2301 + t2080 * t2300 + t2299 * (t2195 * t2291 + t2297 * t2548) / (t2195 * t2297 - t2291 * t2548)) * t2237;
t2028 = t2039 * pkin(1);
t2119 = 0.1e1 / t2613;
t2120 = 0.1e1 / t2613 ^ 2;
t2121 = t2119 * t2120;
t2208 = sin(t2262);
t2553 = pkin(2) * t2291;
t2133 = pkin(3) * t2208 + t2553;
t2285 = xDDP(2);
t2286 = xDDP(1);
t2247 = t2318 + t2320;
t2304 = 0.2e1 * pkin(7);
t2313 = 0.2e1 * qJ(2,1);
t2260 = t2313 + t2304;
t2611 = cos(t2260) * t2318 + cos(t2313) * t2320;
t2357 = -t2247 - t2611;
t2226 = t2277 + qJ(1,1);
t2179 = pkin(7) + t2226;
t2173 = qJ(2,1) + t2179;
t2139 = sin(t2173);
t2180 = -pkin(7) + t2226;
t2174 = -qJ(2,1) + t2180;
t2140 = sin(t2174);
t2145 = cos(t2173);
t2146 = cos(t2174);
t2161 = t2313 + t2179;
t2314 = -0.2e1 * qJ(2,1);
t2162 = t2314 + t2180;
t2193 = t2313 + t2226;
t2175 = t2304 + t2193;
t2194 = t2314 + t2226;
t2305 = -0.2e1 * pkin(7);
t2176 = t2305 + t2194;
t2191 = qJ(2,1) + t2226;
t2192 = -qJ(2,1) + t2226;
t2261 = t2313 + pkin(7);
t2207 = sin(t2261);
t2265 = sin(t2313);
t2472 = t2318 * sin(t2260);
t2434 = 0.2e1 * t2472;
t2437 = 0.2e1 * t2301;
t2438 = 0.2e1 * t2300;
t2246 = pkin(1) * t2301;
t2439 = 0.2e1 * t2250 * t2300 + 0.2e1 * t2246;
t2245 = pkin(1) * t2300;
t2442 = -0.2e1 * t2250 * t2301 + 0.2e1 * t2245;
t2583 = 0.4e1 * t2299;
t2453 = pkin(1) * t2583;
t2586 = 2 * t2247;
t2596 = 0.2e1 * pkin(2);
t2598 = 0.4e1 * pkin(1);
t2011 = t2299 * t2434 + ((t2265 * t2596 + t2291 * t2598) * t2299 + (cos(t2191) + cos(t2192)) * t2439 + (sin(t2191) + sin(t2192)) * t2442) * pkin(2) + (cos(t2226) * t2586 + (cos(t2193) + cos(t2194)) * t2320 + (cos(t2175) + cos(t2176)) * t2318) * t2301 + (sin(t2226) * t2586 + (sin(t2193) + sin(t2194)) * t2320 + (sin(t2176) + sin(t2175)) * t2318) * t2300 + (t2208 * t2453 + (t2145 + t2146) * t2439 + (t2139 + t2140) * t2442 + (t2207 * t2583 + (cos(t2162) + cos(t2179) + cos(t2180) + cos(t2161)) * t2437 + (sin(t2161) + sin(t2162) + sin(t2179) + sin(t2180)) * t2438) * pkin(2)) * pkin(3);
t2491 = t2119 * t2237;
t2427 = t2011 * t2491;
t2388 = t2427 / 0.4e1;
t2332 = t2250 ^ 2;
t2510 = t2011 / t2332;
t2412 = -t2510 / 0.4e1;
t2527 = 0.2e1 * pkin(3);
t2433 = pkin(2) * t2527;
t2222 = cos(t2261);
t2454 = t2222 + t2274;
t2272 = t2299 ^ 2;
t2445 = pkin(2) * t2547;
t2484 = (0.2e1 * t2445 + t2247) * t2272;
t2284 = xDDP(3);
t2490 = t2119 * t2284;
t1994 = -t2613 * t2039 * t2120 * t2412 - (t2133 * t2490 + t2121 * t2484 + t2080 * t2285 + t2081 * t2286 - ((-t2433 * t2454 + t2357) * t2039 / 0.2e1 + t2613 * (-t2028 + t2388)) * t2119 * t2039) * t2237;
t1991 = t1994 * pkin(1);
t2244 = g(2) * t2298;
t2542 = g(1) * t2292;
t2151 = -t2244 + t2542;
t2540 = g(2) * t2292;
t2541 = g(1) * t2298;
t2152 = t2540 + t2541;
t2463 = -t2151 * t2232 - t2152 * t2229;
t2488 = t2120 * t2272;
t2628 = pkin(5) * t2488 - 0.2e1 * t1991 + t2463;
t2276 = legFrame(2,3);
t2228 = sin(t2276);
t2231 = cos(t2276);
t2290 = sin(qJ(1,2));
t2296 = cos(qJ(1,2));
t2078 = t2296 * t2228 + t2290 * t2231;
t2079 = -t2290 * t2228 + t2296 * t2231;
t2279 = qJ(3,2) + pkin(5);
t2249 = -pkin(6) - t2279;
t2235 = 0.1e1 / t2249;
t2289 = sin(qJ(2,2));
t2038 = (t2079 * t2301 + t2078 * t2300 + t2299 * (t2195 * t2289 + t2295 * t2548) / (t2195 * t2295 - t2289 * t2548)) * t2235;
t2032 = pkin(1) * t2038;
t2116 = 0.1e1 / t2614;
t2117 = 0.1e1 / t2614 ^ 2;
t2118 = t2116 * t2117;
t2205 = sin(t2259);
t2554 = pkin(2) * t2289;
t2132 = pkin(3) * t2205 + t2554;
t2310 = 0.2e1 * qJ(2,2);
t2257 = t2310 + t2304;
t2610 = cos(t2257) * t2318 + cos(t2310) * t2320;
t2358 = -t2247 - t2610;
t2225 = t2276 + qJ(1,2);
t2181 = pkin(7) + t2225;
t2170 = qJ(2,2) + t2181;
t2137 = sin(t2170);
t2182 = -pkin(7) + t2225;
t2171 = -qJ(2,2) + t2182;
t2138 = sin(t2171);
t2143 = cos(t2170);
t2144 = cos(t2171);
t2163 = t2310 + t2181;
t2311 = -0.2e1 * qJ(2,2);
t2164 = t2311 + t2182;
t2189 = t2310 + t2225;
t2169 = t2304 + t2189;
t2190 = t2311 + t2225;
t2172 = t2305 + t2190;
t2187 = qJ(2,2) + t2225;
t2188 = -qJ(2,2) + t2225;
t2258 = t2310 + pkin(7);
t2204 = sin(t2258);
t2264 = sin(t2310);
t2473 = t2318 * sin(t2257);
t2435 = 0.2e1 * t2473;
t2440 = 0.2e1 * t2249 * t2300 + 0.2e1 * t2246;
t2443 = -0.2e1 * t2249 * t2301 + 0.2e1 * t2245;
t2012 = t2299 * t2435 + ((t2264 * t2596 + t2289 * t2598) * t2299 + (cos(t2187) + cos(t2188)) * t2440 + (sin(t2187) + sin(t2188)) * t2443) * pkin(2) + (cos(t2225) * t2586 + (cos(t2189) + cos(t2190)) * t2320 + (cos(t2169) + cos(t2172)) * t2318) * t2301 + (sin(t2225) * t2586 + (sin(t2189) + sin(t2190)) * t2320 + (sin(t2169) + sin(t2172)) * t2318) * t2300 + (t2205 * t2453 + (t2143 + t2144) * t2440 + (t2137 + t2138) * t2443 + (t2204 * t2583 + (cos(t2163) + cos(t2164) + cos(t2181) + cos(t2182)) * t2437 + (sin(t2163) + sin(t2164) + sin(t2181) + sin(t2182)) * t2438) * pkin(2)) * pkin(3);
t2495 = t2116 * t2235;
t2426 = t2012 * t2495;
t2387 = t2426 / 0.4e1;
t2331 = t2249 ^ 2;
t2509 = t2012 / t2331;
t2411 = -t2509 / 0.4e1;
t2219 = cos(t2258);
t2455 = t2219 + t2274;
t2494 = t2116 * t2284;
t1993 = -t2614 * t2038 * t2117 * t2411 - (t2132 * t2494 + t2118 * t2484 + t2078 * t2285 + t2079 * t2286 - ((-t2433 * t2455 + t2358) * t2038 / 0.2e1 + t2614 * (-t2032 + t2387)) * t2116 * t2038) * t2235;
t1990 = t1993 * pkin(1);
t2240 = t2290 * g(1);
t2531 = t2296 * g(2);
t2149 = -t2240 + t2531;
t2532 = t2296 * g(1);
t2535 = t2290 * g(2);
t2150 = t2532 + t2535;
t2464 = t2149 * t2231 - t2150 * t2228;
t2492 = t2117 * t2272;
t2627 = pkin(5) * t2492 - 0.2e1 * t1990 + t2464;
t2275 = legFrame(3,3);
t2227 = sin(t2275);
t2230 = cos(t2275);
t2288 = sin(qJ(1,3));
t2294 = cos(qJ(1,3));
t2076 = t2294 * t2227 + t2288 * t2230;
t2077 = -t2288 * t2227 + t2294 * t2230;
t2278 = qJ(3,3) + pkin(5);
t2248 = -pkin(6) - t2278;
t2233 = 0.1e1 / t2248;
t2287 = sin(qJ(2,3));
t2037 = (t2077 * t2301 + t2076 * t2300 + t2299 * (t2195 * t2287 + t2293 * t2548) / (t2195 * t2293 - t2287 * t2548)) * t2233;
t2030 = pkin(1) * t2037;
t2113 = 0.1e1 / t2615;
t2114 = 0.1e1 / t2615 ^ 2;
t2115 = t2113 * t2114;
t2202 = sin(t2256);
t2555 = pkin(2) * t2287;
t2131 = pkin(3) * t2202 + t2555;
t2307 = 0.2e1 * qJ(2,3);
t2254 = t2307 + t2304;
t2609 = cos(t2254) * t2318 + cos(t2307) * t2320;
t2359 = -t2247 - t2609;
t2224 = t2275 + qJ(1,3);
t2177 = pkin(7) + t2224;
t2165 = qJ(2,3) + t2177;
t2135 = sin(t2165);
t2178 = -pkin(7) + t2224;
t2166 = -qJ(2,3) + t2178;
t2136 = sin(t2166);
t2141 = cos(t2165);
t2142 = cos(t2166);
t2159 = t2307 + t2177;
t2308 = -0.2e1 * qJ(2,3);
t2160 = t2308 + t2178;
t2185 = t2307 + t2224;
t2167 = t2304 + t2185;
t2186 = t2308 + t2224;
t2168 = t2305 + t2186;
t2183 = qJ(2,3) + t2224;
t2184 = -qJ(2,3) + t2224;
t2255 = t2307 + pkin(7);
t2201 = sin(t2255);
t2263 = sin(t2307);
t2474 = t2318 * sin(t2254);
t2436 = 0.2e1 * t2474;
t2441 = 0.2e1 * t2248 * t2300 + 0.2e1 * t2246;
t2444 = -0.2e1 * t2248 * t2301 + 0.2e1 * t2245;
t2010 = t2299 * t2436 + ((t2263 * t2596 + t2287 * t2598) * t2299 + (cos(t2183) + cos(t2184)) * t2441 + (sin(t2183) + sin(t2184)) * t2444) * pkin(2) + (cos(t2224) * t2586 + (cos(t2185) + cos(t2186)) * t2320 + (cos(t2167) + cos(t2168)) * t2318) * t2301 + (sin(t2224) * t2586 + (sin(t2185) + sin(t2186)) * t2320 + (sin(t2167) + sin(t2168)) * t2318) * t2300 + (t2202 * t2453 + (t2141 + t2142) * t2441 + (t2135 + t2136) * t2444 + (t2201 * t2583 + (cos(t2177) + cos(t2178) + cos(t2159) + cos(t2160)) * t2437 + (sin(t2177) + sin(t2178) + sin(t2159) + sin(t2160)) * t2438) * pkin(2)) * pkin(3);
t2499 = t2113 * t2233;
t2428 = t2010 * t2499;
t2389 = t2428 / 0.4e1;
t2330 = t2248 ^ 2;
t2511 = t2010 / t2330;
t2413 = -t2511 / 0.4e1;
t2215 = cos(t2255);
t2456 = t2215 + t2274;
t2498 = t2113 * t2284;
t1992 = -t2615 * t2037 * t2114 * t2413 - (t2131 * t2498 + t2115 * t2484 + t2076 * t2285 + t2077 * t2286 - ((-t2433 * t2456 + t2359) * t2037 / 0.2e1 + t2615 * (-t2030 + t2389)) * t2113 * t2037) * t2233;
t1989 = t1992 * pkin(1);
t2239 = t2288 * g(1);
t2533 = t2294 * g(2);
t2147 = -t2239 + t2533;
t2534 = t2294 * g(1);
t2536 = t2288 * g(2);
t2148 = t2534 + t2536;
t2465 = t2147 * t2230 - t2148 * t2227;
t2496 = t2114 * t2272;
t2626 = pkin(5) * t2496 - 0.2e1 * t1989 + t2465;
t2211 = cos(-pkin(7) + qJ(2,3));
t2216 = cos(qJ(2,3) + t2304);
t2558 = pkin(2) * t2201;
t2448 = pkin(3) * t2558;
t2471 = t2320 * t2263;
t2347 = 0.2e1 * t2448 + t2471 + t2474;
t2317 = pkin(3) * t2318;
t2546 = pkin(3) * t2320;
t2399 = -0.2e1 * t2317 - 0.4e1 * t2546;
t2449 = -0.2e1 * t2546;
t2552 = pkin(2) * t2318;
t2451 = -0.2e1 * t2552;
t2497 = t2113 * t2299;
t2319 = pkin(2) * t2320;
t2588 = -0.2e1 * t2319 - 0.4e1 * t2552;
t2589 = -0.4e1 * pkin(1) * (t2445 + t2320 / 0.2e1 + t2318 / 0.2e1);
t2425 = (t2347 * t2248 * t2037 + (t2211 * t2449 + t2216 * t2451 + t2217 * t2399 + t2293 * t2588 + t2589) * t2497) * t2114 * t2299;
t2619 = 0.2e1 * pkin(1);
t2530 = pkin(3) * t2619;
t1995 = (t2217 * t2530 + (t2293 * t2619 + t2456 * t2527) * pkin(2) - t2359) * t2037;
t2431 = t1995 * t2511;
t2067 = t2131 * t2619 + t2347;
t2502 = t2067 * t2284;
t2196 = t2241 + pkin(1);
t2070 = t2196 * t2288 + t2248 * t2294;
t2074 = t2196 * t2294 - t2248 * t2288;
t2041 = -t2070 * t2227 + t2074 * t2230 + t2077 * t2551;
t2507 = t2041 * t2286;
t2040 = t2070 * t2230 + t2074 * t2227 + t2076 * t2551;
t2508 = t2040 * t2285;
t2303 = 0.3e1 * pkin(7);
t2306 = 0.3e1 * qJ(2,3);
t2321 = pkin(1) ^ 2;
t2450 = -0.3e1 * t2546;
t2452 = -0.3e1 * t2552;
t2529 = -0.8e1 * pkin(2) * pkin(3);
t2612 = -0.8e1 * t2445 - (4 * t2247);
t2526 = (-(t2436 + 0.4e1 * t2448 + 0.2e1 * t2471) * t2248 * t2497 - 0.4e1 * t2615 * (pkin(1) * t2389 - (t2321 + t2330) * t2037) + (t2317 * cos(t2306 + t2303) + t2551 * t2629 + t2319 * cos(t2306) + t2241 * t2630 - (cos(t2304 + t2306) + t2216) * t2452 - (cos(t2306 + pkin(7)) + t2211) * t2450) * t2037 + (t2215 * t2529 - 0.4e1 * t2609 + t2612) * (-t2030 + t2428 / 0.8e1)) * t2037;
t2625 = (t2113 * (t2502 / 0.2e1 + t2526 / 0.4e1) + t2508 + t2507 - t2425 / 0.2e1) * t2233 + t1989 - t2114 * t2431 / 0.8e1;
t2281 = 0.2e1 * t2278;
t2624 = (t2113 * (t2502 + t2526 / 0.2e1) - t2425 + 0.2e1 * t2507 + 0.2e1 * t2508) * t2233 - (-t1995 * t2413 + t2281 * t2272) * t2114 + 0.4e1 * t1989;
t2209 = cos(t2304 + qJ(2,2));
t2212 = cos(-pkin(7) + qJ(2,2));
t2557 = pkin(2) * t2204;
t2447 = pkin(3) * t2557;
t2470 = t2320 * t2264;
t2346 = 0.2e1 * t2447 + t2470 + t2473;
t2493 = t2116 * t2299;
t2424 = (t2346 * t2249 * t2038 + (t2209 * t2451 + t2212 * t2449 + t2220 * t2399 + t2295 * t2588 + t2589) * t2493) * t2117 * t2299;
t1996 = (t2220 * t2530 + (t2295 * t2619 + t2455 * t2527) * pkin(2) - t2358) * t2038;
t2430 = t1996 * t2509;
t2068 = t2132 * t2619 + t2346;
t2501 = t2068 * t2284;
t2197 = t2242 + pkin(1);
t2071 = t2197 * t2290 + t2249 * t2296;
t2075 = t2197 * t2296 - t2249 * t2290;
t2043 = -t2071 * t2228 + t2075 * t2231 + t2079 * t2550;
t2505 = t2043 * t2286;
t2042 = t2071 * t2231 + t2075 * t2228 + t2078 * t2550;
t2506 = t2042 * t2285;
t2309 = 0.3e1 * qJ(2,2);
t2525 = (-(t2435 + 0.4e1 * t2447 + 0.2e1 * t2470) * t2249 * t2493 - 0.4e1 * t2614 * (pkin(1) * t2387 - (t2321 + t2331) * t2038) + (t2317 * cos(t2309 + t2303) + t2550 * t2629 + t2319 * cos(t2309) + t2242 * t2630 - (cos(t2304 + t2309) + t2209) * t2452 - (cos(t2309 + pkin(7)) + t2212) * t2450) * t2038 + (t2219 * t2529 - 0.4e1 * t2610 + t2612) * (-t2032 + t2426 / 0.8e1)) * t2038;
t2623 = (t2116 * (t2501 / 0.2e1 + t2525 / 0.4e1) + t2506 + t2505 - t2424 / 0.2e1) * t2235 + t1990 - t2117 * t2430 / 0.8e1;
t2282 = 0.2e1 * t2279;
t2622 = (t2116 * (t2501 + t2525 / 0.2e1) - t2424 + 0.2e1 * t2505 + 0.2e1 * t2506) * t2235 - (-t1996 * t2411 + t2282 * t2272) * t2117 + 0.4e1 * t1990;
t2210 = cos(t2304 + qJ(2,1));
t2213 = cos(-pkin(7) + qJ(2,1));
t2556 = pkin(2) * t2207;
t2446 = pkin(3) * t2556;
t2469 = t2320 * t2265;
t2345 = 0.2e1 * t2446 + t2469 + t2472;
t2489 = t2119 * t2299;
t2423 = (t2345 * t2250 * t2039 + (t2210 * t2451 + t2213 * t2449 + t2223 * t2399 + t2297 * t2588 + t2589) * t2489) * t2120 * t2299;
t1997 = (t2223 * t2530 + (t2297 * t2619 + t2454 * t2527) * pkin(2) - t2357) * t2039;
t2429 = t1997 * t2510;
t2069 = t2133 * t2619 + t2345;
t2500 = t2069 * t2284;
t2198 = t2243 + pkin(1);
t2072 = t2198 * t2292 + t2250 * t2298;
t2073 = t2198 * t2298 - t2250 * t2292;
t2045 = -t2072 * t2229 + t2073 * t2232 + t2081 * t2549;
t2503 = t2045 * t2286;
t2044 = t2072 * t2232 + t2073 * t2229 + t2080 * t2549;
t2504 = t2044 * t2285;
t2312 = 0.3e1 * qJ(2,1);
t2524 = (-(t2434 + 0.4e1 * t2446 + 0.2e1 * t2469) * t2250 * t2489 - 0.4e1 * t2613 * (pkin(1) * t2388 - (t2321 + t2332) * t2039) + (t2317 * cos(t2303 + t2312) + t2549 * t2629 + t2319 * cos(t2312) + t2243 * t2630 - (cos(t2312 + t2304) + t2210) * t2452 - (t2213 + cos(pkin(7) + t2312)) * t2450) * t2039 + (t2222 * t2529 - 0.4e1 * t2611 + t2612) * (-t2028 + t2427 / 0.8e1)) * t2039;
t2621 = (t2119 * (t2500 / 0.2e1 + t2524 / 0.4e1) + t2504 + t2503 - t2423 / 0.2e1) * t2237 + t1991 - t2120 * t2429 / 0.8e1;
t2283 = 2 * t2280;
t2620 = (t2119 * (t2500 + t2524 / 0.2e1) - t2423 + 0.2e1 * t2503 + 0.2e1 * t2504) * t2237 - (-t1997 * t2412 + t2283 * t2272) * t2120 + 0.4e1 * t1991;
t2618 = t2147 * t2227 + t2148 * t2230;
t2617 = t2149 * t2228 + t2150 * t2231;
t2616 = -t2151 * t2229 + t2152 * t2232;
t2597 = -0.2e1 * pkin(2);
t2419 = t2037 * t2497;
t2022 = pkin(1) * t2419;
t2058 = t2115 * t2131 * t2272 + t2498;
t2576 = t2058 / 0.2e1;
t2019 = pkin(5) * t2576 - t2022;
t2592 = -0.2e1 * t2019;
t2418 = t2038 * t2493;
t2023 = pkin(1) * t2418;
t2059 = t2118 * t2132 * t2272 + t2494;
t2575 = t2059 / 0.2e1;
t2020 = pkin(5) * t2575 - t2023;
t2591 = -0.2e1 * t2020;
t2417 = t2039 * t2489;
t2024 = pkin(1) * t2417;
t2060 = t2121 * t2133 * t2272 + t2490;
t2574 = t2060 / 0.2e1;
t2021 = pkin(5) * t2574 - t2024;
t2590 = -0.2e1 * t2021;
t2585 = -0.2e1 * t2273;
t2584 = 0.2e1 * t2274;
t2582 = -pkin(5) / 0.2e1;
t2581 = -g(1) / 0.2e1;
t2580 = g(1) / 0.2e1;
t2579 = -g(2) / 0.2e1;
t2578 = g(2) / 0.2e1;
t2577 = pkin(1) * g(2);
t2573 = -t2067 / 0.2e1;
t2572 = -t2068 / 0.2e1;
t2571 = -t2069 / 0.2e1;
t2570 = t2135 / 0.2e1;
t2569 = t2137 / 0.2e1;
t2568 = t2139 / 0.2e1;
t2567 = t2142 / 0.2e1;
t2566 = t2144 / 0.2e1;
t2565 = t2146 / 0.2e1;
t2564 = t2202 / 0.2e1;
t2563 = t2205 / 0.2e1;
t2562 = t2208 / 0.2e1;
t2561 = t2217 / 0.2e1;
t2560 = t2220 / 0.2e1;
t2559 = t2223 / 0.2e1;
t2545 = pkin(5) * t1992;
t2544 = pkin(5) * t1993;
t2543 = pkin(5) * t1994;
t2539 = g(3) * t2293;
t2538 = g(3) * t2295;
t2537 = g(3) * t2297;
t2523 = t1992 * t2273;
t2522 = t1992 * t2274;
t2521 = t1992 * t2278;
t2520 = t1992 * t2287;
t2519 = t1993 * t2273;
t2518 = t1993 * t2274;
t2517 = t1993 * t2279;
t2516 = t1993 * t2289;
t2515 = t1994 * t2273;
t2514 = t1994 * t2274;
t2513 = t1994 * t2280;
t2512 = t1994 * t2291;
t2487 = t2131 * t2233;
t2486 = t2132 * t2235;
t2485 = t2133 * t2237;
t2477 = t2293 * t1992;
t2476 = t2295 * t1993;
t2475 = t2297 * t1994;
t2462 = g(2) * t2570 + t2141 * t2580;
t2461 = g(2) * t2569 + t2143 * t2580;
t2460 = g(2) * t2568 + t2145 * t2580;
t2459 = g(1) * t2570 + t2141 * t2579;
t2458 = g(1) * t2569 + t2143 * t2579;
t2457 = g(1) * t2568 + t2145 * t2579;
t2034 = t2037 ^ 2;
t2422 = t2034 * t2293 * t2287;
t2035 = t2038 ^ 2;
t2421 = t2035 * t2295 * t2289;
t2036 = t2039 ^ 2;
t2420 = t2036 * t2297 * t2291;
t2416 = t2113 * t2487;
t2415 = t2116 * t2486;
t2414 = t2119 * t2485;
t2029 = pkin(1) * t2034;
t2408 = t2029 + t2618;
t2031 = pkin(1) * t2035;
t2407 = t2031 + t2617;
t2033 = pkin(1) * t2036;
t2406 = t2033 + t2616;
t2398 = -0.2e1 * t2419;
t2397 = -0.2e1 * t2418;
t2396 = -0.2e1 * t2417;
t2395 = t2037 * t2428;
t2394 = t2038 * t2426;
t2393 = t2039 * t2427;
t2392 = t2287 * t2419;
t2391 = t2289 * t2418;
t2390 = t2291 * t2417;
t2386 = t2499 * t2573;
t2385 = t2495 * t2572;
t2384 = t2491 * t2571;
t2383 = t2419 * t2597;
t2382 = t2418 * t2597;
t2381 = t2417 * t2597;
t2004 = t2395 / 0.2e1;
t2005 = t2394 / 0.2e1;
t2006 = t2393 / 0.2e1;
t2374 = g(2) * t2567 + t2136 * t2581;
t2373 = g(2) * t2566 + t2138 * t2581;
t2372 = g(2) * t2565 + t2140 * t2581;
t2371 = g(1) * t2567 + t2136 * t2578;
t2370 = g(1) * t2566 + t2138 * t2578;
t2369 = g(1) * t2565 + t2140 * t2578;
t2046 = t2058 * t2287 + t2293 * t2496;
t2048 = t2059 * t2289 + t2295 * t2492;
t2050 = t2060 * t2291 + t2297 * t2488;
t2353 = t2274 * t2398 + t2523;
t2352 = t2273 * t2398 - t2522;
t2351 = t2274 * t2397 + t2519;
t2350 = t2273 * t2397 - t2518;
t2349 = t2274 * t2396 + t2515;
t2348 = t2273 * t2396 - t2514;
t1956 = -t2278 * t2034 + (-0.2e1 * t2392 - t2477) * pkin(2) + t2465 - t2625;
t1959 = t2287 * t2353 + t2293 * t2352;
t1964 = -t2287 * t2352 + t2293 * t2353;
t2338 = t1959 * MDP(11) + t1964 * MDP(12) - t2034 * MDP(13) + t1956 * MDP(14);
t1957 = -t2279 * t2035 + (-0.2e1 * t2391 - t2476) * pkin(2) + t2464 - t2623;
t1960 = t2289 * t2351 + t2295 * t2350;
t1962 = -t2289 * t2350 + t2295 * t2351;
t2337 = t1960 * MDP(11) + t1962 * MDP(12) - t2035 * MDP(13) + t1957 * MDP(14);
t1958 = -t2280 * t2036 + (-0.2e1 * t2390 - t2475) * pkin(2) + t2463 - t2621;
t1961 = t2291 * t2349 + t2297 * t2348;
t1963 = -t2291 * t2348 + t2297 * t2349;
t2336 = t1961 * MDP(11) + t1963 * MDP(12) - t2036 * MDP(13) + t1958 * MDP(14);
t2302 = pkin(1) * g(1);
t2153 = -t2278 * g(2) + t2302;
t2154 = t2278 * g(1) + t2577;
t2269 = t2293 ^ 2;
t2316 = pkin(5) ^ 2;
t1947 = 0.2e1 * ((-t2533 / 0.2e1 + t2239 / 0.2e1) * t2230 + (t2534 / 0.2e1 + t2536 / 0.2e1) * t2227 + pkin(2) * t2392 + t1989 + (-t2431 / 0.16e2 + (t2582 - qJ(3,3) / 0.2e1) * t2272) * t2114 - (-t2507 / 0.2e1 - t2508 / 0.2e1 + t2425 / 0.4e1 + (-t2502 / 0.4e1 - t2526 / 0.8e1) * t2113) * t2233) * t2241 + (t2153 * t2288 - t2154 * t2294) * t2230 + (t2153 * t2294 + t2154 * t2288) * t2227 - 0.2e1 * (qJ(3,3) * t2576 + t2019) * t2555 + (t2004 + 0.2e1 * t2545) * qJ(3,3) + pkin(5) * t2004 + pkin(1) * t2625 + (qJ(3,3) ^ 2 + t2269 * t2320 + t2316) * t1992;
t2016 = -0.2e1 * t2058 * t2278 + 0.4e1 * t2022;
t1950 = -t2201 * t2383 + t2624 * t2561 + t2016 * t2564 - t2374 + t2459 + (t1992 * t2215 + t2522) * pkin(2);
t1953 = -pkin(2) * t2523 - t1992 * t2558 + t2016 * t2561 - t2215 * t2383 - t2624 * t2564 - t2371 + t2462;
t1968 = -pkin(2) * t2046 + t1992 * t2281 + t2004 - t2618;
t1974 = 0.2e1 * t2287 * t2477 - (0.4e1 * t2269 - 0.2e1) * t2419;
t1977 = t2626 * t2287 + t2293 * t2592;
t1978 = t2287 * t2592 - t2626 * t2293;
t1986 = (t2293 * t2398 + t2520) * t2287;
t2047 = t2058 * t2293 - t2287 * t2496;
t2094 = g(1) * t2227 - g(2) * t2230;
t2097 = g(1) * t2230 + g(2) * t2227;
t2061 = t2094 * t2294 + t2097 * t2288;
t2064 = -t2094 * t2288 + t2097 * t2294;
t2335 = t1992 * MDP(1) + t1977 * MDP(10) + t1950 * MDP(11) + t1953 * MDP(12) + t1968 * MDP(13) + t1947 * MDP(14) + t2061 * MDP(2) + t2064 * MDP(3) + t1986 * MDP(4) + t1974 * MDP(5) + t2046 * MDP(6) + t2047 * MDP(7) + t1978 * MDP(9);
t2155 = -t2279 * g(2) + t2302;
t2156 = t2279 * g(1) + t2577;
t2270 = t2295 ^ 2;
t1948 = 0.2e1 * ((-t2531 / 0.2e1 + t2240 / 0.2e1) * t2231 + (t2532 / 0.2e1 + t2535 / 0.2e1) * t2228 + pkin(2) * t2391 + t1990 + (-t2430 / 0.16e2 + (t2582 - qJ(3,2) / 0.2e1) * t2272) * t2117 - (-t2505 / 0.2e1 - t2506 / 0.2e1 + t2424 / 0.4e1 + (-t2501 / 0.4e1 - t2525 / 0.8e1) * t2116) * t2235) * t2242 + (t2155 * t2290 - t2156 * t2296) * t2231 + (t2155 * t2296 + t2156 * t2290) * t2228 - 0.2e1 * (qJ(3,2) * t2575 + t2020) * t2554 + (t2005 + 0.2e1 * t2544) * qJ(3,2) + pkin(5) * t2005 + pkin(1) * t2623 + (qJ(3,2) ^ 2 + t2270 * t2320 + t2316) * t1993;
t2017 = -0.2e1 * t2059 * t2279 + 0.4e1 * t2023;
t1951 = -t2204 * t2382 + t2622 * t2560 + t2017 * t2563 - t2373 + t2458 + (t1993 * t2219 + t2518) * pkin(2);
t1954 = -pkin(2) * t2519 - t1993 * t2557 + t2017 * t2560 - t2219 * t2382 - t2622 * t2563 - t2370 + t2461;
t1969 = -pkin(2) * t2048 + t1993 * t2282 + t2005 - t2617;
t1975 = 0.2e1 * t2289 * t2476 - (0.4e1 * t2270 - 0.2e1) * t2418;
t1979 = t2627 * t2289 + t2295 * t2591;
t1980 = t2289 * t2591 - t2627 * t2295;
t1987 = (t2295 * t2397 + t2516) * t2289;
t2049 = t2059 * t2295 - t2289 * t2492;
t2095 = g(1) * t2228 - g(2) * t2231;
t2098 = g(1) * t2231 + g(2) * t2228;
t2062 = t2095 * t2296 + t2098 * t2290;
t2065 = -t2095 * t2290 + t2098 * t2296;
t2334 = t1993 * MDP(1) + t1979 * MDP(10) + t1951 * MDP(11) + t1954 * MDP(12) + t1969 * MDP(13) + t1948 * MDP(14) + t2062 * MDP(2) + t2065 * MDP(3) + t1987 * MDP(4) + t1975 * MDP(5) + t2048 * MDP(6) + t2049 * MDP(7) + t1980 * MDP(9);
t2157 = -t2280 * g(2) + t2302;
t2158 = t2280 * g(1) + t2577;
t2271 = t2297 ^ 2;
t1949 = 0.2e1 * ((-t2244 / 0.2e1 + t2542 / 0.2e1) * t2232 + (t2541 / 0.2e1 + t2540 / 0.2e1) * t2229 + pkin(2) * t2390 + t1991 + (-t2429 / 0.16e2 + (t2582 - qJ(3,1) / 0.2e1) * t2272) * t2120 - (-t2503 / 0.2e1 - t2504 / 0.2e1 + t2423 / 0.4e1 + (-t2500 / 0.4e1 - t2524 / 0.8e1) * t2119) * t2237) * t2243 + (t2157 * t2292 - t2158 * t2298) * t2232 + (t2157 * t2298 + t2158 * t2292) * t2229 - 0.2e1 * (qJ(3,1) * t2574 + t2021) * t2553 + (t2006 + 0.2e1 * t2543) * qJ(3,1) + pkin(5) * t2006 + pkin(1) * t2621 + (qJ(3,1) ^ 2 + t2271 * t2320 + t2316) * t1994;
t2018 = -0.2e1 * t2060 * t2280 + 0.4e1 * t2024;
t1952 = -t2207 * t2381 + t2620 * t2559 + t2018 * t2562 - t2372 + t2457 + (t1994 * t2222 + t2514) * pkin(2);
t1955 = -pkin(2) * t2515 - t1994 * t2556 + t2018 * t2559 - t2222 * t2381 - t2620 * t2562 - t2369 + t2460;
t1970 = -pkin(2) * t2050 + t1994 * t2283 + t2006 - t2616;
t1976 = 0.2e1 * t2291 * t2475 - (0.4e1 * t2271 - 0.2e1) * t2417;
t1981 = t2291 * t2590 - t2628 * t2297;
t1982 = t2628 * t2291 + t2297 * t2590;
t1988 = (t2297 * t2396 + t2512) * t2291;
t2051 = t2060 * t2297 - t2291 * t2488;
t2096 = g(1) * t2229 - g(2) * t2232;
t2099 = g(1) * t2232 + g(2) * t2229;
t2063 = t2096 * t2298 + t2099 * t2292;
t2066 = -t2096 * t2292 + t2099 * t2298;
t2333 = t1994 * MDP(1) + t1982 * MDP(10) + t1952 * MDP(11) + t1955 * MDP(12) + t1970 * MDP(13) + t1949 * MDP(14) + t2063 * MDP(2) + t2066 * MDP(3) + t1988 * MDP(4) + t1976 * MDP(5) + t2050 * MDP(6) + t2051 * MDP(7) + t1981 * MDP(9);
t1985 = t2406 - t2543;
t1984 = t2407 - t2544;
t1983 = t2408 - t2545;
t1973 = 0.2e1 * t2033 - t2393 - 0.2e1 * t2513;
t1972 = 0.2e1 * t2031 - t2394 - 0.2e1 * t2517;
t1971 = 0.2e1 * t2029 - t2395 - 0.2e1 * t2521;
t1 = [(t2286 - g(1)) * MDP(15) - (t2045 * t2336 + t2081 * t2333) * t2237 - (t2043 * t2337 + t2079 * t2334) * t2235 - (t2041 * t2338 + t2077 * t2335) * t2233; (t2285 - g(2)) * MDP(15) - (t2044 * t2336 + t2080 * t2333) * t2237 - (t2042 * t2337 + t2078 * t2334) * t2235 - (t2040 * t2338 + t2076 * t2335) * t2233; (-t1992 * t2416 - t1993 * t2415 - t1994 * t2414) * MDP(1) + (-t2061 * t2416 - t2062 * t2415 - t2063 * t2414) * MDP(2) + (-t2064 * t2416 - t2065 * t2415 - t2066 * t2414) * MDP(3) + ((-t1988 * t2485 - t2420) * t2119 + (-t1987 * t2486 - t2421) * t2116 + (-t1986 * t2487 - t2422) * t2113) * MDP(4) + ((-t1976 * t2485 - 0.2e1 * t2036 * t2271 + t2036) * t2119 + (-t1975 * t2486 - 0.2e1 * t2035 * t2270 + t2035) * t2116 + (-t1974 * t2487 - 0.2e1 * t2034 * t2269 + t2034) * t2113) * MDP(5) + ((-t2050 * t2485 + t2512) * t2119 + (-t2048 * t2486 + t2516) * t2116 + (-t2046 * t2487 + t2520) * t2113) * MDP(6) + ((-t2051 * t2485 + t2475) * t2119 + (-t2049 * t2486 + t2476) * t2116 + (-t2047 * t2487 + t2477) * t2113) * MDP(7) + (t2113 * t2058 + t2116 * t2059 + t2119 * t2060) * MDP(8) + ((-t1981 * t2485 + t1985 * t2291 - t2537) * t2119 + (-t1980 * t2486 + t1984 * t2289 - t2538) * t2116 + (-t1978 * t2487 + t1983 * t2287 - t2539) * t2113) * MDP(9) + ((g(3) * t2291 - t1982 * t2485 + t1985 * t2297) * t2119 + (g(3) * t2289 - t1979 * t2486 + t1984 * t2295) * t2116 + (g(3) * t2287 - t1977 * t2487 + t1983 * t2293) * t2113) * MDP(10) + (-t1952 * t2414 + t2119 * (t1973 * t2562 - g(3) * t2223 + (t2036 * t2207 + t2060 * t2584) * pkin(2) + t2372 + t2457) + t1961 * t2384 - t1951 * t2415 + t2116 * (t1972 * t2563 - g(3) * t2220 + (t2035 * t2204 + t2059 * t2584) * pkin(2) + t2373 + t2458) + t1960 * t2385 - t1950 * t2416 + t2113 * (t1971 * t2564 - g(3) * t2217 + (t2034 * t2201 + t2058 * t2584) * pkin(2) + t2374 + t2459) + t1959 * t2386) * MDP(11) + (-t1955 * t2414 + t2119 * (t1973 * t2559 + g(3) * t2208 + (t2036 * t2222 + t2060 * t2585) * pkin(2) + t2369 + t2460) + t1963 * t2384 - t1954 * t2415 + t2116 * (t1972 * t2560 + g(3) * t2205 + (t2035 * t2219 + t2059 * t2585) * pkin(2) + t2370 + t2461) + t1962 * t2385 - t1953 * t2416 + t2113 * (t1971 * t2561 + g(3) * t2202 + (t2034 * t2215 + t2058 * t2585) * pkin(2) + t2371 + t2462) + t1964 * t2386) * MDP(12) + (-(t2133 * t1970 + t2036 * t2571) * t2491 - (t2132 * t1969 + t2035 * t2572) * t2495 - (t2131 * t1968 + t2034 * t2573) * t2499 + (-t2113 * t2520 - t2116 * t2516 - t2119 * t2512) * pkin(2)) * MDP(13) + (-(t2133 * t1949 + t2069 * t1958 / 0.2e1) * t2491 - (t2132 * t1948 + t2068 * t1957 / 0.2e1) * t2495 - (t2131 * t1947 + t2067 * t1956 / 0.2e1) * t2499 + (t2119 * ((-t2393 / 0.2e1 - t2513 + t2406) * t2291 - t2537) + t2116 * ((-t2394 / 0.2e1 - t2517 + t2407) * t2289 - t2538) + t2113 * ((-t2395 / 0.2e1 - t2521 + t2408) * t2287 - t2539) + (t2119 * (t2060 + t2420) + t2116 * (t2059 + t2421) + t2113 * (t2058 + t2422)) * pkin(2)) * pkin(2)) * MDP(14) + (t2284 - g(3)) * MDP(15);];
tauX  = t1;
