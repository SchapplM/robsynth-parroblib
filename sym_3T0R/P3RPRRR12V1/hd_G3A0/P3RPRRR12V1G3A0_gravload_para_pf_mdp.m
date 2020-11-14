% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RPRRR12V1G3A0
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
%   see P3RPRRR12V1G3A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:28
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRRR12V1G3A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(14,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_mdp: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'P3RPRRR12V1G3A0_gravload_para_pf_mdp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:28:44
% EndTime: 2020-08-06 18:28:45
% DurationCPUTime: 1.10s
% Computational Cost: add. (423->132), mult. (774->250), div. (102->7), fcn. (738->18), ass. (0->104)
t2291 = legFrame(3,2);
t2278 = sin(t2291);
t2281 = cos(t2291);
t2266 = t2281 * g(1) - t2278 * g(2);
t2295 = sin(qJ(1,3));
t2301 = cos(qJ(1,3));
t2258 = -g(3) * t2295 + t2266 * t2301;
t2294 = sin(qJ(3,3));
t2272 = t2294 * pkin(3) + qJ(2,3);
t2269 = 0.1e1 / t2272;
t2339 = t2295 * t2269;
t2333 = t2258 * t2339;
t2292 = legFrame(2,2);
t2279 = sin(t2292);
t2282 = cos(t2292);
t2267 = t2282 * g(1) - t2279 * g(2);
t2297 = sin(qJ(1,2));
t2303 = cos(qJ(1,2));
t2260 = -g(3) * t2297 + t2267 * t2303;
t2296 = sin(qJ(3,2));
t2273 = t2296 * pkin(3) + qJ(2,2);
t2270 = 0.1e1 / t2273;
t2338 = t2297 * t2270;
t2331 = t2260 * t2338;
t2293 = legFrame(1,2);
t2280 = sin(t2293);
t2283 = cos(t2293);
t2268 = t2283 * g(1) - t2280 * g(2);
t2299 = sin(qJ(1,1));
t2305 = cos(qJ(1,1));
t2262 = -g(3) * t2299 + t2268 * t2305;
t2298 = sin(qJ(3,1));
t2274 = t2298 * pkin(3) + qJ(2,1);
t2271 = 0.1e1 / t2274;
t2337 = t2299 * t2271;
t2329 = t2262 * t2337;
t2365 = MDP(2) - MDP(4);
t2364 = MDP(3) - MDP(5);
t2261 = g(3) * t2305 + t2268 * t2299;
t2259 = g(3) * t2303 + t2267 * t2297;
t2257 = g(3) * t2301 + t2266 * t2295;
t2288 = 0.1e1 / t2294;
t2360 = t2257 * t2288;
t2289 = 0.1e1 / t2296;
t2359 = t2259 * t2289;
t2290 = 0.1e1 / t2298;
t2358 = t2261 * t2290;
t2354 = t2269 * t2301;
t2353 = t2270 * t2303;
t2352 = t2271 * t2305;
t2351 = t2278 * t2288;
t2300 = cos(qJ(3,3));
t2350 = t2278 * t2300;
t2349 = t2279 * t2289;
t2302 = cos(qJ(3,2));
t2348 = t2279 * t2302;
t2347 = t2280 * t2290;
t2304 = cos(qJ(3,1));
t2346 = t2280 * t2304;
t2345 = t2281 * t2288;
t2344 = t2281 * t2300;
t2343 = t2282 * t2289;
t2342 = t2282 * t2302;
t2341 = t2283 * t2290;
t2340 = t2283 * t2304;
t2336 = t2301 * t2278;
t2335 = t2303 * t2279;
t2334 = t2305 * t2280;
t2332 = t2258 * t2354;
t2330 = t2260 * t2353;
t2328 = t2262 * t2352;
t2327 = t2269 * t2336;
t2326 = t2281 * t2354;
t2325 = t2270 * t2335;
t2324 = t2282 * t2353;
t2323 = t2271 * t2334;
t2322 = t2283 * t2352;
t2275 = t2299 * qJ(2,1);
t2287 = pkin(1) + pkin(5) + pkin(6);
t2321 = t2287 * t2305 + t2275;
t2276 = qJ(2,3) * t2295;
t2320 = t2287 * t2301 + t2276;
t2277 = qJ(2,2) * t2297;
t2319 = t2287 * t2303 + t2277;
t2318 = t2294 * t2332;
t2317 = t2300 * t2332;
t2316 = t2296 * t2330;
t2315 = t2302 * t2330;
t2314 = t2298 * t2328;
t2313 = t2304 * t2328;
t2306 = 0.1e1 / pkin(3);
t2265 = t2280 * g(1) + t2283 * g(2);
t2264 = t2279 * g(1) + t2282 * g(2);
t2263 = t2278 * g(1) + t2281 * g(2);
t2250 = t2268 * (t2299 * pkin(1) - t2305 * qJ(2,1)) + g(3) * (t2305 * pkin(1) + t2275);
t2249 = t2267 * (t2297 * pkin(1) - t2303 * qJ(2,2)) + g(3) * (t2303 * pkin(1) + t2277);
t2248 = t2266 * (t2295 * pkin(1) - t2301 * qJ(2,3)) + g(3) * (t2301 * pkin(1) + t2276);
t2247 = t2261 * t2298 + t2265 * t2304;
t2246 = -t2261 * t2304 + t2265 * t2298;
t2245 = t2259 * t2296 + t2264 * t2302;
t2244 = -t2259 * t2302 + t2264 * t2296;
t2243 = t2257 * t2294 + t2263 * t2300;
t2242 = -t2257 * t2300 + t2263 * t2294;
t1 = [(t2250 * t2322 - (t2321 * t2283 * t2298 + qJ(2,1) * t2346 + (t2298 * t2346 + (-t2304 ^ 2 + 0.1e1) * t2283 * t2299) * pkin(3)) * t2271 * t2358 + t2249 * t2324 - (t2319 * t2282 * t2296 + qJ(2,2) * t2348 + (t2296 * t2348 + (-t2302 ^ 2 + 0.1e1) * t2282 * t2297) * pkin(3)) * t2270 * t2359 + t2248 * t2326 - (t2320 * t2281 * t2294 + qJ(2,3) * t2350 + (t2294 * t2350 + (-t2300 ^ 2 + 0.1e1) * t2281 * t2295) * pkin(3)) * t2269 * t2360) * MDP(6) + (-t2281 * t2318 - t2282 * t2316 - t2283 * t2314 + (-t2242 * t2351 - t2244 * t2349 - t2246 * t2347) * t2306) * MDP(12) + (-t2281 * t2317 - t2282 * t2315 - t2283 * t2313 + (-t2243 * t2351 - t2245 * t2349 - t2247 * t2347) * t2306) * MDP(13) - g(1) * MDP(14) + t2365 * (t2257 * t2326 + t2259 * t2324 + t2261 * t2322) + t2364 * (t2258 * t2326 + t2260 * t2324 + t2262 * t2322); ((-t2250 * t2334 - ((pkin(3) * t2340 - t2321 * t2280) * t2298 + t2299 * pkin(3) * (t2304 - 0.1e1) * (t2304 + 0.1e1) * t2280 + qJ(2,1) * t2340) * t2358) * t2271 + (-t2249 * t2335 - ((pkin(3) * t2342 - t2319 * t2279) * t2296 + t2297 * pkin(3) * (t2302 - 0.1e1) * (t2302 + 0.1e1) * t2279 + qJ(2,2) * t2342) * t2359) * t2270 + (-t2248 * t2336 - ((pkin(3) * t2344 - t2320 * t2278) * t2294 + t2295 * pkin(3) * (t2300 - 0.1e1) * (t2300 + 0.1e1) * t2278 + qJ(2,3) * t2344) * t2360) * t2269) * MDP(6) + (t2278 * t2318 + t2279 * t2316 + t2280 * t2314 + (-t2242 * t2345 - t2244 * t2343 - t2246 * t2341) * t2306) * MDP(12) + (t2278 * t2317 + t2279 * t2315 + t2280 * t2313 + (-t2243 * t2345 - t2245 * t2343 - t2247 * t2341) * t2306) * MDP(13) - g(2) * MDP(14) - t2365 * (t2257 * t2327 + t2259 * t2325 + t2261 * t2323) - t2364 * (t2258 * t2327 + t2260 * t2325 + t2262 * t2323); ((-t2299 * t2250 - (t2274 * t2305 - t2287 * t2299) * t2261) * t2271 + (-t2297 * t2249 - (t2273 * t2303 - t2287 * t2297) * t2259) * t2270 + (-t2295 * t2248 - (t2272 * t2301 - t2287 * t2295) * t2257) * t2269) * MDP(6) + (t2294 * t2333 + t2296 * t2331 + t2298 * t2329) * MDP(12) + (t2300 * t2333 + t2302 * t2331 + t2304 * t2329) * MDP(13) - g(3) * MDP(14) - t2365 * (t2257 * t2339 + t2259 * t2338 + t2261 * t2337) - t2364 * (t2329 + t2331 + t2333);];
taugX  = t1;
