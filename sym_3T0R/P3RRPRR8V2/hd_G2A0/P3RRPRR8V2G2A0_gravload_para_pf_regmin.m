% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRPRR8V2G2A0
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x15]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2022-11-07 13:08
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRPRR8V2G2A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-07 13:07:13
% EndTime: 2022-11-07 13:07:14
% DurationCPUTime: 1.55s
% Computational Cost: add. (1206->211), mult. (2040->354), div. (165->9), fcn. (1815->26), ass. (0->181)
t2280 = qJ(3,1) + pkin(5);
t2268 = pkin(6) + t2280;
t2252 = 0.1e1 / t2268;
t2295 = cos(qJ(1,1));
t2363 = t2252 * t2295;
t2289 = sin(qJ(1,1));
t2283 = legFrame(1,2);
t2255 = sin(t2283);
t2258 = cos(t2283);
t2307 = t2258 * g(1) - t2255 * g(2);
t2392 = g(3) * t2289 - t2307 * t2295;
t2400 = t2392 * t2363;
t2277 = cos(pkin(7));
t2378 = t2277 * pkin(3);
t2239 = pkin(2) + t2378;
t2294 = cos(qJ(2,1));
t2276 = sin(pkin(7));
t2288 = sin(qJ(2,1));
t2351 = t2276 * t2288;
t2340 = pkin(3) * t2351;
t2375 = 0.1e1 / (t2239 * t2294 - t2340) * t2252;
t2399 = t2392 * t2375;
t2279 = qJ(3,2) + pkin(5);
t2267 = pkin(6) + t2279;
t2251 = 0.1e1 / t2267;
t2293 = cos(qJ(1,2));
t2364 = t2251 * t2293;
t2287 = sin(qJ(1,2));
t2282 = legFrame(2,2);
t2254 = sin(t2282);
t2257 = cos(t2282);
t2308 = t2257 * g(1) - t2254 * g(2);
t2393 = g(3) * t2287 - t2308 * t2293;
t2398 = t2393 * t2364;
t2292 = cos(qJ(2,2));
t2286 = sin(qJ(2,2));
t2352 = t2276 * t2286;
t2341 = pkin(3) * t2352;
t2376 = 0.1e1 / (t2239 * t2292 - t2341) * t2251;
t2397 = t2393 * t2376;
t2278 = qJ(3,3) + pkin(5);
t2266 = pkin(6) + t2278;
t2250 = 0.1e1 / t2266;
t2291 = cos(qJ(1,3));
t2365 = t2250 * t2291;
t2285 = sin(qJ(1,3));
t2281 = legFrame(3,2);
t2253 = sin(t2281);
t2256 = cos(t2281);
t2309 = t2256 * g(1) - t2253 * g(2);
t2394 = g(3) * t2285 - t2309 * t2291;
t2396 = t2394 * t2365;
t2290 = cos(qJ(2,3));
t2284 = sin(qJ(2,3));
t2353 = t2276 * t2284;
t2342 = pkin(3) * t2353;
t2377 = 0.1e1 / (t2239 * t2290 - t2342) * t2250;
t2395 = t2394 * t2377;
t2391 = 0.2e1 * t2290 ^ 2;
t2390 = 0.2e1 * t2292 ^ 2;
t2389 = 0.2e1 * t2294 ^ 2;
t2270 = qJ(2,3) + pkin(7);
t2247 = cos(t2270);
t2388 = pkin(3) * t2247;
t2271 = qJ(2,2) + pkin(7);
t2248 = cos(t2271);
t2387 = pkin(3) * t2248;
t2272 = qJ(2,1) + pkin(7);
t2249 = cos(t2272);
t2386 = pkin(3) * t2249;
t2324 = pkin(1) * t2285 - t2291 * t2266;
t2269 = t2277 ^ 2;
t2349 = pkin(3) * (t2269 - 0.1e1);
t2382 = (t2285 * t2349 + t2324 * t2353) * pkin(3);
t2323 = pkin(1) * t2287 - t2293 * t2267;
t2381 = (t2287 * t2349 + t2323 * t2352) * pkin(3);
t2322 = pkin(1) * t2289 - t2295 * t2268;
t2380 = (t2289 * t2349 + t2322 * t2351) * pkin(3);
t2379 = t2276 * pkin(3);
t2259 = t2290 * pkin(2);
t2230 = 0.1e1 / (t2259 + t2388);
t2374 = t2230 * t2253;
t2373 = t2230 * t2256;
t2260 = t2292 * pkin(2);
t2231 = 0.1e1 / (t2260 + t2387);
t2372 = t2231 * t2254;
t2371 = t2231 * t2257;
t2261 = t2294 * pkin(2);
t2232 = 0.1e1 / (t2261 + t2386);
t2370 = t2232 * t2255;
t2369 = t2232 * t2258;
t2368 = t2239 * t2256;
t2367 = t2239 * t2257;
t2366 = t2239 * t2258;
t2362 = t2253 * t2239;
t2361 = t2253 * t2285;
t2360 = t2254 * t2239;
t2359 = t2254 * t2287;
t2358 = t2255 * t2239;
t2357 = t2255 * t2289;
t2356 = t2256 * t2285;
t2355 = t2257 * t2287;
t2354 = t2258 * t2289;
t2350 = pkin(2) * t2378;
t2348 = t2256 * t2379;
t2347 = t2257 * t2379;
t2346 = t2258 * t2379;
t2345 = t2253 * t2379;
t2344 = t2254 * t2379;
t2343 = t2255 * t2379;
t2305 = g(3) * t2291 + t2309 * t2285;
t2330 = t2305 * t2377;
t2304 = g(3) * t2293 + t2308 * t2287;
t2329 = t2304 * t2376;
t2303 = g(3) * t2295 + t2307 * t2289;
t2328 = t2303 * t2375;
t2321 = t2247 * t2395;
t2320 = t2290 * t2395;
t2319 = t2248 * t2397;
t2318 = t2292 * t2397;
t2317 = t2249 * t2399;
t2316 = t2294 * t2399;
t2244 = sin(t2270);
t2315 = t2244 * t2395;
t2314 = t2284 * t2395;
t2245 = sin(t2271);
t2313 = t2245 * t2397;
t2312 = t2286 * t2397;
t2246 = sin(t2272);
t2311 = t2246 * t2399;
t2310 = t2288 * t2399;
t2296 = pkin(3) ^ 2;
t2297 = pkin(2) ^ 2;
t2306 = 0.2e1 * t2269 * t2296 - t2296 + t2297 + 0.2e1 * t2350;
t2227 = t2253 * g(1) + t2256 * g(2);
t2193 = -t2227 * t2290 + t2284 * t2305;
t2228 = t2254 * g(1) + t2257 * g(2);
t2195 = -t2228 * t2292 + t2286 * t2304;
t2229 = t2255 * g(1) + t2258 * g(2);
t2197 = -t2229 * t2294 + t2288 * t2303;
t2302 = t2193 * t2374 + t2195 * t2372 + t2197 * t2370;
t2301 = t2193 * t2373 + t2195 * t2371 + t2197 * t2369;
t2300 = t2303 * t2363 + t2304 * t2364 + t2305 * t2365;
t2178 = (-t2239 * t2361 + t2348) * t2290 + (t2285 * t2345 + t2368) * t2284;
t2179 = (-t2239 * t2359 + t2347) * t2292 + (t2287 * t2344 + t2367) * t2286;
t2180 = (-t2239 * t2357 + t2346) * t2294 + (t2289 * t2343 + t2366) * t2288;
t2299 = t2178 * t2330 + t2179 * t2329 + t2180 * t2328;
t2181 = (t2239 * t2356 + t2345) * t2290 + t2284 * (-t2285 * t2348 + t2362);
t2182 = (t2239 * t2355 + t2344) * t2292 + t2286 * (-t2287 * t2347 + t2360);
t2183 = (t2239 * t2354 + t2343) * t2294 + t2288 * (-t2289 * t2346 + t2358);
t2298 = t2181 * t2330 + t2182 * t2329 + t2183 * t2328;
t2243 = pkin(1) * t2379;
t2242 = t2261 + pkin(1);
t2241 = t2260 + pkin(1);
t2240 = t2259 + pkin(1);
t2238 = t2241 * t2293;
t2237 = t2240 * t2291;
t2236 = t2295 * t2242;
t2235 = pkin(1) * t2288 - t2379;
t2234 = pkin(1) * t2286 - t2379;
t2233 = pkin(1) * t2284 - t2379;
t2226 = t2350 + t2297 / 0.2e1 + (t2269 - 0.1e1 / 0.2e1) * t2296;
t2207 = -0.2e1 * t2289 * t2340 + t2322;
t2206 = -0.2e1 * t2287 * t2341 + t2323;
t2205 = -0.2e1 * t2285 * t2342 + t2324;
t2204 = t2306 * t2288 + t2243;
t2203 = t2306 * t2286 + t2243;
t2202 = t2306 * t2284 + t2243;
t2198 = t2229 * t2288 + t2294 * t2303;
t2196 = t2228 * t2286 + t2292 * t2304;
t2194 = t2227 * t2284 + t2290 * t2305;
t2192 = t2229 * t2246 + t2249 * t2303;
t2191 = -t2229 * t2249 + t2246 * t2303;
t2190 = t2228 * t2245 + t2248 * t2304;
t2189 = -t2228 * t2248 + t2245 * t2304;
t2188 = t2227 * t2244 + t2247 * t2305;
t2187 = -t2227 * t2247 + t2244 * t2305;
t2186 = (t2289 * t2242 - t2280 * t2295) * g(3) - (t2289 * t2280 + t2236) * t2307;
t2185 = (t2287 * t2241 - t2279 * t2293) * g(3) - (t2287 * t2279 + t2238) * t2308;
t2184 = (t2285 * t2240 - t2278 * t2291) * g(3) - t2309 * (t2285 * t2278 + t2237);
t1 = [0, t2181 * t2395 + t2182 * t2397 + t2183 * t2399, t2298, 0, 0, 0, 0, 0, t2181 * t2320 + t2182 * t2318 + t2183 * t2316 + t2302, -t2181 * t2314 - t2182 * t2312 - t2183 * t2310 + t2194 * t2374 + t2196 * t2372 + t2198 * t2370, t2181 * t2321 + t2182 * t2319 + t2183 * t2317 + t2187 * t2374 + t2189 * t2372 + t2191 * t2370, -t2181 * t2315 - t2182 * t2313 - t2183 * t2311 + t2188 * t2374 + t2190 * t2372 + t2192 * t2370, -t2298, (t2183 * t2186 - ((t2226 * t2354 + t2239 * t2343) * t2389 + (t2204 * t2255 + t2207 * t2366) * t2294 - t2258 * t2380 + t2235 * t2358) * t2392) * t2375 + (t2182 * t2185 - ((t2226 * t2355 + t2239 * t2344) * t2390 + (t2203 * t2254 + t2206 * t2367) * t2292 - t2257 * t2381 + t2234 * t2360) * t2393) * t2376 + (t2181 * t2184 - ((t2226 * t2356 + t2239 * t2345) * t2391 + (t2202 * t2253 + t2205 * t2368) * t2290 - t2256 * t2382 + t2233 * t2362) * t2394) * t2377 + t2302 * pkin(2), -g(1); 0, t2178 * t2395 + t2179 * t2397 + t2180 * t2399, t2299, 0, 0, 0, 0, 0, t2178 * t2320 + t2179 * t2318 + t2180 * t2316 + t2301, -t2178 * t2314 - t2179 * t2312 - t2180 * t2310 + t2194 * t2373 + t2196 * t2371 + t2198 * t2369, t2178 * t2321 + t2179 * t2319 + t2180 * t2317 + t2187 * t2373 + t2189 * t2371 + t2191 * t2369, -t2178 * t2315 - t2179 * t2313 - t2180 * t2311 + t2188 * t2373 + t2190 * t2371 + t2192 * t2369, -t2299, (t2180 * t2186 - ((-t2226 * t2357 + t2239 * t2346) * t2389 + (t2258 * t2204 - t2207 * t2358) * t2294 + t2255 * t2380 + t2235 * t2366) * t2392) * t2375 + (t2179 * t2185 - ((-t2226 * t2359 + t2239 * t2347) * t2390 + (t2257 * t2203 - t2206 * t2360) * t2292 + t2254 * t2381 + t2234 * t2367) * t2393) * t2376 + (t2178 * t2184 - ((-t2226 * t2361 + t2239 * t2348) * t2391 + (t2256 * t2202 - t2205 * t2362) * t2290 + t2253 * t2382 + t2233 * t2368) * t2394) * t2377 + t2301 * pkin(2), -g(2); 0, t2400 + t2398 + t2396, t2300, 0, 0, 0, 0, 0, t2290 * t2396 + t2292 * t2398 + t2294 * t2400, -t2284 * t2396 - t2286 * t2398 - t2288 * t2400, t2247 * t2396 + t2248 * t2398 + t2249 * t2400, -t2244 * t2396 - t2245 * t2398 - t2246 * t2400, -t2300, (t2295 * t2186 - (t2289 * t2268 + t2295 * t2386 + t2236) * t2392) * t2252 + (t2293 * t2185 - (t2287 * t2267 + t2293 * t2387 + t2238) * t2393) * t2251 + (t2291 * t2184 - (t2285 * t2266 + t2291 * t2388 + t2237) * t2394) * t2250, -g(3);];
tau_reg  = t1;
