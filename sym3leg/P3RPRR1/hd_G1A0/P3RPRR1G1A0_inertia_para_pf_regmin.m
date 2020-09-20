% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RPRR1G1P1A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3*3x8]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:23
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RPRR1G1P1A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G1P1A0_inertia_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G1P1A0_inertia_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G1P1A0_inertia_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G1P1A0_inertia_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G1P1A0_inertia_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:23:18
% EndTime: 2020-03-09 21:23:19
% DurationCPUTime: 0.40s
% Computational Cost: add. (1708->124), mult. (1092->222), div. (276->7), fcn. (984->29), ass. (0->131)
t308 = sin(pkin(7));
t368 = pkin(1) * t308;
t305 = legFrame(3,3) + qJ(1,3);
t301 = pkin(7) + t305;
t298 = qJ(3,3) + t301;
t292 = sin(t298);
t256 = -pkin(1) * sin(t305) - pkin(2) * sin(t301) - pkin(3) * t292;
t313 = sin(qJ(3,3));
t283 = pkin(1) * sin(pkin(7) + qJ(3,3)) + t313 * pkin(2);
t277 = 0.1e1 / t283;
t367 = t256 * t277;
t306 = legFrame(2,3) + qJ(1,2);
t302 = pkin(7) + t306;
t299 = qJ(3,2) + t302;
t293 = sin(t299);
t257 = -pkin(1) * sin(t306) - pkin(2) * sin(t302) - pkin(3) * t293;
t314 = sin(qJ(3,2));
t284 = pkin(1) * sin(pkin(7) + qJ(3,2)) + t314 * pkin(2);
t279 = 0.1e1 / t284;
t366 = t257 * t279;
t307 = legFrame(1,3) + qJ(1,1);
t303 = pkin(7) + t307;
t300 = qJ(3,1) + t303;
t294 = sin(t300);
t258 = -pkin(1) * sin(t307) - pkin(2) * sin(t303) - pkin(3) * t294;
t315 = sin(qJ(3,1));
t285 = pkin(1) * sin(pkin(7) + qJ(3,1)) + t315 * pkin(2);
t281 = 0.1e1 / t285;
t365 = t258 * t281;
t295 = cos(t298);
t259 = -pkin(1) * cos(t305) - pkin(2) * cos(t301) - pkin(3) * t295;
t364 = t259 * t277;
t296 = cos(t299);
t260 = -pkin(1) * cos(t306) - pkin(2) * cos(t302) - pkin(3) * t296;
t363 = t260 * t279;
t297 = cos(t300);
t261 = -pkin(1) * cos(t307) - pkin(2) * cos(t303) - pkin(3) * t297;
t362 = t261 * t281;
t319 = 0.1e1 / pkin(3);
t361 = t277 * t319;
t278 = 0.1e1 / t283 ^ 2;
t360 = t278 * t292;
t359 = t278 * t295;
t358 = t279 * t319;
t280 = 0.1e1 / t284 ^ 2;
t357 = t280 * t293;
t356 = t280 * t296;
t355 = t281 * t319;
t282 = 0.1e1 / t285 ^ 2;
t354 = t282 * t294;
t353 = t282 * t297;
t262 = t292 * t277;
t263 = t293 * t279;
t264 = t294 * t281;
t265 = t295 * t277;
t266 = t296 * t279;
t267 = t297 * t281;
t352 = t308 * t313;
t351 = t308 * t314;
t350 = t308 * t315;
t309 = cos(pkin(7));
t304 = t309 * pkin(1) + pkin(2);
t316 = cos(qJ(3,3));
t274 = t313 * t304 + t316 * t368;
t349 = -0.2e1 * t274 * t277;
t317 = cos(qJ(3,2));
t275 = t314 * t304 + t317 * t368;
t348 = -0.2e1 * t275 * t279;
t318 = cos(qJ(3,1));
t276 = t315 * t304 + t318 * t368;
t347 = -0.2e1 * t276 * t281;
t268 = -t316 * pkin(2) + (-t309 * t316 + t352) * pkin(1);
t346 = t268 * t262;
t345 = t268 * t265;
t269 = -t317 * pkin(2) + (-t309 * t317 + t351) * pkin(1);
t344 = t269 * t263;
t343 = t269 * t266;
t270 = -t318 * pkin(2) + (-t309 * t318 + t350) * pkin(1);
t342 = t270 * t264;
t341 = t270 * t267;
t271 = -pkin(1) * t352 + t316 * t304;
t340 = t271 * t360;
t339 = t271 * t359;
t272 = -pkin(1) * t351 + t317 * t304;
t338 = t272 * t357;
t337 = t272 * t356;
t273 = -pkin(1) * t350 + t318 * t304;
t336 = t273 * t354;
t335 = t273 * t353;
t334 = t274 * t360;
t333 = t274 * t359;
t332 = t275 * t357;
t331 = t275 * t356;
t330 = t276 * t354;
t329 = t276 * t353;
t328 = t292 * t349;
t327 = t295 * t349;
t326 = t293 * t348;
t325 = t296 * t348;
t324 = t294 * t347;
t323 = t297 * t347;
t322 = t278 * t292 ^ 2 + t280 * t293 ^ 2 + t282 * t294 ^ 2;
t321 = t278 * t295 ^ 2 + t280 * t296 ^ 2 + t282 * t297 ^ 2;
t231 = t292 * t359 + t293 * t356 + t294 * t353;
t320 = pkin(1) ^ 2;
t255 = t261 * t355;
t254 = t260 * t358;
t253 = t259 * t361;
t252 = t258 * t355;
t251 = t257 * t358;
t250 = t256 * t361;
t249 = t267 + t255;
t248 = t266 + t254;
t247 = t265 + t253;
t246 = t264 + t252;
t245 = t263 + t251;
t244 = t262 + t250;
t243 = t255 + 0.2e1 * t267;
t242 = t254 + 0.2e1 * t266;
t241 = t253 + 0.2e1 * t265;
t240 = t252 + 0.2e1 * t264;
t239 = t251 + 0.2e1 * t263;
t238 = t250 + 0.2e1 * t262;
t237 = t267 + t255 / 0.2e1;
t236 = t266 + t254 / 0.2e1;
t235 = t265 + t253 / 0.2e1;
t234 = t264 + t252 / 0.2e1;
t233 = t263 + t251 / 0.2e1;
t232 = t262 + t250 / 0.2e1;
t230 = t231 * t320;
t1 = [t321, 0, 0, t321 * t320, t247 * t265 + t248 * t266 + t249 * t267 + (t247 * t364 + t248 * t363 + t249 * t362) * t319, -t241 * t345 - t242 * t343 - t243 * t341 + (t259 * t339 + t260 * t337 + t261 * t335) * t319, t235 * t327 + t236 * t325 + t237 * t323 + (-t259 * t333 - t260 * t331 - t261 * t329) * t319, 1; t231, 0, 0, t230, t244 * t265 + t245 * t266 + t246 * t267 + (t244 * t364 + t245 * t363 + t246 * t362) * t319, -t238 * t345 - t239 * t343 - t240 * t341 + (t259 * t340 + t260 * t338 + t261 * t336) * t319, t232 * t327 + t233 * t325 + t234 * t323 + (-t259 * t334 - t260 * t332 - t261 * t330) * t319, 0; 0, 0, 0, 0, 0, 0, 0, 0; t231, 0, 0, t230, t247 * t262 + t248 * t263 + t249 * t264 + (t247 * t367 + t248 * t366 + t249 * t365) * t319, -t241 * t346 - t242 * t344 - t243 * t342 + (t256 * t339 + t257 * t337 + t258 * t335) * t319, t235 * t328 + t236 * t326 + t237 * t324 + (-t256 * t333 - t257 * t331 - t258 * t329) * t319, 0; t322, 0, 0, t322 * t320, t244 * t262 + t245 * t263 + t246 * t264 + (t244 * t367 + t245 * t366 + t246 * t365) * t319, -t238 * t346 - t239 * t344 - t240 * t342 + (t256 * t340 + t257 * t338 + t258 * t336) * t319, t232 * t328 + t233 * t326 + t234 * t324 + (-t256 * t334 - t257 * t332 - t258 * t330) * t319, 1; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 3, 0, 0, 0, 1;];
tau_reg  = t1;
