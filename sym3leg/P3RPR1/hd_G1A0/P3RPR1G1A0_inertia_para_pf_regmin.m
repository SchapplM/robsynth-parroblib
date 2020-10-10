% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RPR1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [2x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3*(3+1)/2x10]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:58
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RPR1G1A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(4,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPR1G1A0_inertia_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RPR1G1A0_inertia_para_pf_regmin: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RPR1G1A0_inertia_para_pf_regmin: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPR1G1A0_inertia_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPR1G1A0_inertia_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:58:11
% EndTime: 2019-05-03 14:58:12
% DurationCPUTime: 0.79s
% Computational Cost: add. (1463->185), mult. (2278->369), div. (249->6), fcn. (1958->14), ass. (0->145)
t320 = legFrame(3,3);
t309 = sin(t320);
t312 = cos(t320);
t323 = sin(qJ(1,3));
t326 = cos(qJ(1,3));
t278 = t309 * t326 + t312 * t323;
t279 = -t309 * t323 + t312 * t326;
t330 = xP(3);
t318 = sin(t330);
t319 = cos(t330);
t337 = koppelP(3,2);
t340 = koppelP(3,1);
t284 = -t318 * t340 - t319 * t337;
t287 = -t318 * t337 + t319 * t340;
t331 = 0.1e1 / qJ(2,3);
t233 = (t278 * t287 + t279 * t284) * t331;
t379 = pkin(1) * t233;
t321 = legFrame(2,3);
t310 = sin(t321);
t313 = cos(t321);
t324 = sin(qJ(1,2));
t327 = cos(qJ(1,2));
t280 = t310 * t327 + t313 * t324;
t281 = -t310 * t324 + t313 * t327;
t338 = koppelP(2,2);
t341 = koppelP(2,1);
t285 = -t318 * t341 - t319 * t338;
t288 = -t318 * t338 + t319 * t341;
t333 = 0.1e1 / qJ(2,2);
t234 = (t280 * t288 + t281 * t285) * t333;
t378 = pkin(1) * t234;
t322 = legFrame(1,3);
t311 = sin(t322);
t314 = cos(t322);
t325 = sin(qJ(1,1));
t328 = cos(qJ(1,1));
t282 = t311 * t328 + t314 * t325;
t283 = -t311 * t325 + t314 * t328;
t339 = koppelP(1,2);
t342 = koppelP(1,1);
t286 = -t318 * t342 - t319 * t339;
t289 = -t318 * t339 + t319 * t342;
t335 = 0.1e1 / qJ(2,1);
t235 = (t282 * t289 + t283 * t286) * t335;
t377 = pkin(1) * t235;
t376 = pkin(1) * t278;
t375 = pkin(1) * t279;
t374 = pkin(1) * t280;
t373 = pkin(1) * t281;
t372 = pkin(1) * t282;
t371 = pkin(1) * t283;
t290 = t323 * t340 - t326 * t337;
t291 = t323 * t337 + t326 * t340;
t230 = (t290 * t319 - t318 * t291) * t312 + t309 * (t318 * t290 + t291 * t319);
t370 = t230 * t233;
t292 = t324 * t341 - t327 * t338;
t293 = t324 * t338 + t327 * t341;
t231 = (t292 * t319 - t318 * t293) * t313 + t310 * (t318 * t292 + t293 * t319);
t369 = t231 * t234;
t294 = t325 * t342 - t328 * t339;
t295 = t325 * t339 + t328 * t342;
t232 = (t294 * t319 - t318 * t295) * t314 + t311 * (t318 * t294 + t295 * t319);
t368 = t232 * t235;
t329 = pkin(1) + pkin(2);
t296 = -t326 * qJ(2,3) + t323 * t329;
t299 = t323 * qJ(2,3) + t329 * t326;
t260 = t296 * t312 + t309 * t299;
t367 = (-t260 + 0.2e1 * t376) * t331 ^ 2;
t263 = -t309 * t296 + t299 * t312;
t244 = (-t263 + 0.2e1 * t375) * t331;
t366 = t244 * t331;
t297 = -t327 * qJ(2,2) + t324 * t329;
t300 = t324 * qJ(2,2) + t329 * t327;
t261 = t297 * t313 + t310 * t300;
t365 = (-t261 + 0.2e1 * t374) * t333 ^ 2;
t264 = -t310 * t297 + t300 * t313;
t248 = (-t264 + 0.2e1 * t373) * t333;
t364 = t248 * t333;
t298 = -t328 * qJ(2,1) + t325 * t329;
t301 = t325 * qJ(2,1) + t329 * t328;
t262 = t298 * t314 + t311 * t301;
t363 = (-t262 + 0.2e1 * t372) * t335 ^ 2;
t265 = -t311 * t298 + t301 * t314;
t252 = (-t265 + 0.2e1 * t371) * t335;
t362 = t252 * t335;
t361 = t278 * t331;
t346 = qJ(2,3) ^ 2;
t332 = 0.1e1 / t346;
t360 = t278 * t332;
t359 = t279 * t331;
t358 = t279 * t332;
t357 = t280 * t333;
t345 = qJ(2,2) ^ 2;
t334 = 0.1e1 / t345;
t356 = t280 * t334;
t355 = t281 * t333;
t354 = t281 * t334;
t353 = t282 * t335;
t344 = qJ(2,1) ^ 2;
t336 = 0.1e1 / t344;
t352 = t282 * t336;
t351 = t283 * t335;
t350 = t283 * t336;
t349 = (t260 * t287 + t263 * t284) * t331;
t348 = (t261 * t288 + t264 * t285) * t333;
t347 = (t262 * t289 + t265 * t286) * t335;
t343 = pkin(1) ^ 2;
t317 = t343 + t344;
t316 = t343 + t345;
t315 = t343 + t346;
t308 = qJ(2,1) * t342 + t329 * t339;
t307 = qJ(2,2) * t341 + t329 * t338;
t306 = qJ(2,3) * t340 + t329 * t337;
t305 = -qJ(2,1) * t339 + t329 * t342;
t304 = -qJ(2,2) * t338 + t329 * t341;
t303 = -qJ(2,3) * t337 + t329 * t340;
t302 = t318 ^ 2 + t319 ^ 2;
t277 = t283 ^ 2;
t276 = t282 ^ 2;
t275 = t281 ^ 2;
t274 = t280 ^ 2;
t273 = t279 ^ 2;
t272 = t278 ^ 2;
t271 = t325 * t305 - t308 * t328;
t270 = t324 * t304 - t307 * t327;
t269 = t323 * t303 - t306 * t326;
t268 = t305 * t328 + t308 * t325;
t267 = t304 * t327 + t307 * t324;
t266 = t303 * t326 + t306 * t323;
t253 = (t265 - t371) * t335;
t251 = (t262 - t372) * t335;
t249 = (t264 - t373) * t333;
t247 = (t261 - t374) * t333;
t245 = (t263 - t375) * t331;
t243 = (t260 - t376) * t331;
t241 = (-pkin(1) * t265 + t283 * t317) * t335;
t240 = (-pkin(1) * t262 + t282 * t317) * t335;
t239 = (-pkin(1) * t264 + t281 * t316) * t333;
t238 = (-pkin(1) * t261 + t280 * t316) * t333;
t237 = (-pkin(1) * t263 + t279 * t315) * t331;
t236 = (-pkin(1) * t260 + t278 * t315) * t331;
t229 = (-t318 * t268 + t271 * t319) * t314 + (t268 * t319 + t271 * t318) * t311;
t228 = (-t318 * t267 + t270 * t319) * t313 + (t267 * t319 + t270 * t318) * t310;
t227 = (-t318 * t266 + t269 * t319) * t312 + (t266 * t319 + t269 * t318) * t309;
t1 = [t332 * t273 + t334 * t275 + t336 * t277, 0, 0, (-t265 * t336 + t362) * t283 + (-t264 * t334 + t364) * t281 + (-t263 * t332 + t366) * t279, 0.2e1 * t273 * t331 + 0.2e1 * t275 * t333 + 0.2e1 * t277 * t335, (t241 * t283 + t253 * t265) * t335 + (t239 * t281 + t249 * t264) * t333 + (t237 * t279 + t245 * t263) * t331, 0, 0, 0, t302; t278 * t358 + t280 * t354 + t282 * t350, 0, 0, t244 * t361 + t248 * t357 + t252 * t353 - t260 * t358 - t261 * t354 - t262 * t350, 0.2e1 * t278 * t359 + 0.2e1 * t280 * t355 + 0.2e1 * t282 * t351, (t241 * t282 + t253 * t262) * t335 + (t239 * t280 + t249 * t261) * t333 + (t237 * t278 + t245 * t260) * t331, 0, 0, 0, 0; t332 * t272 + t334 * t274 + t336 * t276, 0, 0, (-t262 * t336 + t363) * t282 + (-t261 * t334 + t365) * t280 + (-t260 * t332 + t367) * t278, 0.2e1 * t272 * t331 + 0.2e1 * t274 * t333 + 0.2e1 * t276 * t335, (t240 * t282 + t251 * t262) * t335 + (t238 * t280 + t247 * t261) * t333 + (t236 * t278 + t243 * t260) * t331, 0, 0, 0, t302; t230 * t358 + t231 * t354 + t232 * t350, 0, 0, -t227 * t358 - t228 * t354 - t229 * t350 + t230 * t366 + t231 * t364 + t232 * t362, 0.2e1 * t230 * t359 + 0.2e1 * t231 * t355 + 0.2e1 * t232 * t351, (t229 * t253 + t232 * t241) * t335 + (t228 * t249 + t231 * t239) * t333 + (t227 * t245 + t230 * t237) * t331, 0, -t318, -t319, 0; t230 * t360 + t231 * t356 + t232 * t352, 0, 0, -t227 * t360 - t228 * t356 - t229 * t352 + t230 * t367 + t231 * t365 + t232 * t363, 0.2e1 * t230 * t361 + 0.2e1 * t231 * t357 + 0.2e1 * t232 * t353, (t229 * t251 + t232 * t240) * t335 + (t228 * t247 + t231 * t238) * t333 + (t227 * t243 + t230 * t236) * t331, 0, t319, -t318, 0; t331 * t370 + t333 * t369 + t335 * t368, 0, 0, (t232 * (-t347 + 0.2e1 * t377) - t229 * t235) * t335 + (t231 * (-t348 + 0.2e1 * t378) - t228 * t234) * t333 + (t230 * (-t349 + 0.2e1 * t379) - t227 * t233) * t331, 0.2e1 * t368 + 0.2e1 * t369 + 0.2e1 * t370, (t232 * (-pkin(1) * t347 + t317 * t235) + t229 * (t347 - t377)) * t335 + (t231 * (-pkin(1) * t348 + t316 * t234) + t228 * (t348 - t378)) * t333 + (t230 * (-pkin(1) * t349 + t315 * t233) + t227 * (t349 - t379)) * t331, 1, 0, 0, 0;];
tau_reg  = t1;
