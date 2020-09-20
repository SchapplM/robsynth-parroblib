% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3PRRR1G1A0
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2020-03-09 21:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRR1G1A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G1A0_inertia_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G1A0_inertia_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G1A0_inertia_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G1A0_inertia_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G1A0_inertia_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:14:54
% EndTime: 2020-03-09 21:14:54
% DurationCPUTime: 0.55s
% Computational Cost: add. (3044->105), mult. (1452->219), div. (456->6), fcn. (2112->24), ass. (0->123)
t284 = pkin(7) + qJ(2,1);
t312 = qJ(3,1) + t284;
t272 = cos(t312);
t364 = legFrame(1,3);
t281 = cos(t364);
t297 = sin(t312);
t315 = sin(t364);
t265 = t315 * t272 + t281 * t297;
t275 = sin(t284);
t278 = cos(t284);
t365 = 0.1e1 / (-t275 * t272 + t297 * t278);
t341 = t365 * t265;
t283 = pkin(7) + qJ(2,2);
t311 = qJ(3,2) + t283;
t271 = cos(t311);
t363 = legFrame(2,3);
t280 = cos(t363);
t296 = sin(t311);
t314 = sin(t363);
t263 = t314 * t271 + t280 * t296;
t274 = sin(t283);
t277 = cos(t283);
t366 = 0.1e1 / (-t274 * t271 + t296 * t277);
t345 = t366 * t263;
t282 = pkin(7) + qJ(2,3);
t310 = qJ(3,3) + t282;
t270 = cos(t310);
t362 = legFrame(3,3);
t279 = cos(t362);
t295 = sin(t310);
t313 = sin(t362);
t261 = t313 * t270 + t279 * t295;
t273 = sin(t282);
t276 = cos(t282);
t367 = 0.1e1 / (-t273 * t270 + t295 * t276);
t349 = t367 * t261;
t246 = pkin(2) * (t273 * t279 + t313 * t276) + t261 * pkin(3);
t361 = t246 * t367;
t247 = pkin(2) * (t274 * t280 + t314 * t277) + t263 * pkin(3);
t360 = t247 * t366;
t248 = pkin(2) * (t275 * t281 + t315 * t278) + t265 * pkin(3);
t359 = t248 * t365;
t262 = t279 * t270 - t295 * t313;
t249 = -pkin(2) * (t273 * t313 - t279 * t276) + t262 * pkin(3);
t358 = t249 * t367;
t264 = t280 * t271 - t296 * t314;
t250 = -pkin(2) * (t274 * t314 - t280 * t277) + t264 * pkin(3);
t357 = t250 * t366;
t266 = t281 * t272 - t297 * t315;
t251 = -pkin(2) * (t275 * t315 - t281 * t278) + t266 * pkin(3);
t356 = t251 * t365;
t355 = t367 ^ 2;
t292 = 0.1e1 / pkin(2);
t354 = t367 * t292;
t353 = t366 ^ 2;
t352 = t366 * t292;
t351 = t365 ^ 2;
t350 = t365 * t292;
t348 = t367 * t262;
t347 = t367 * sin(qJ(3,3));
t346 = t367 * cos(qJ(3,3));
t344 = t366 * t264;
t343 = t366 * sin(qJ(3,2));
t342 = t366 * cos(qJ(3,2));
t340 = t365 * t266;
t339 = t365 * sin(qJ(3,1));
t338 = t365 * cos(qJ(3,1));
t291 = 0.1e1 / pkin(3);
t337 = t291 * t292;
t336 = t261 * t354;
t335 = t262 * t354;
t334 = t263 * t352;
t333 = t264 * t352;
t332 = t265 * t350;
t331 = t266 * t350;
t330 = t261 * t347;
t329 = t261 * t346;
t328 = t262 * t347;
t327 = t262 * t346;
t326 = t367 * t337;
t325 = t263 * t343;
t324 = t263 * t342;
t323 = t264 * t343;
t322 = t264 * t342;
t321 = t366 * t337;
t320 = t265 * t339;
t319 = t265 * t338;
t318 = t266 * t339;
t317 = t266 * t338;
t316 = t365 * t337;
t309 = t367 * t330;
t308 = t367 * t329;
t307 = t367 * t328;
t306 = t367 * t327;
t305 = t366 * t325;
t304 = t366 * t324;
t303 = t366 * t323;
t302 = t366 * t322;
t301 = t365 * t320;
t300 = t365 * t319;
t299 = t365 * t318;
t298 = t365 * t317;
t293 = 0.1e1 / pkin(2) ^ 2;
t294 = (t340 * t341 + t344 * t345 + t348 * t349) * t293;
t245 = t251 * t316;
t244 = t250 * t321;
t243 = t249 * t326;
t242 = t248 * t316;
t241 = t247 * t321;
t240 = t246 * t326;
t239 = -t245 + t331;
t238 = -t242 + t332;
t237 = -t244 + t333;
t236 = -t241 + t334;
t235 = -t243 + t335;
t234 = -t240 + t336;
t233 = -t245 + 0.2e1 * t331;
t232 = -t244 + 0.2e1 * t333;
t231 = -t243 + 0.2e1 * t335;
t230 = -t242 + 0.2e1 * t332;
t229 = -t241 + 0.2e1 * t334;
t228 = -t240 + 0.2e1 * t336;
t1 = [0, (t262 ^ 2 * t355 + t264 ^ 2 * t353 + t266 ^ 2 * t351) * t293, 0, 0, (t235 * t348 + t237 * t344 + t239 * t340 + (-t235 * t358 - t237 * t357 - t239 * t356) * t291) * t292, t231 * t327 + t232 * t322 + t233 * t317 + (-t249 * t306 - t250 * t302 - t251 * t298) * t337, -t231 * t328 - t232 * t323 - t233 * t318 + (t249 * t307 + t250 * t303 + t251 * t299) * t337, 1; 0, t294, 0, 0, (t234 * t348 + t236 * t344 + t238 * t340 + (-t234 * t358 - t236 * t357 - t238 * t356) * t291) * t292, t228 * t327 + t229 * t322 + t230 * t317 + (-t249 * t308 - t250 * t304 - t251 * t300) * t337, -t228 * t328 - t229 * t323 - t230 * t318 + (t249 * t309 + t250 * t305 + t251 * t301) * t337, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, t294, 0, 0, (t235 * t349 + t237 * t345 + t239 * t341 + (-t235 * t361 - t237 * t360 - t239 * t359) * t291) * t292, t231 * t329 + t232 * t324 + t233 * t319 + (-t246 * t306 - t247 * t302 - t248 * t298) * t337, -t231 * t330 - t232 * t325 - t233 * t320 + (t246 * t307 + t247 * t303 + t248 * t299) * t337, 0; 0, (t261 ^ 2 * t355 + t263 ^ 2 * t353 + t265 ^ 2 * t351) * t293, 0, 0, (t234 * t349 + t236 * t345 + t238 * t341 + (-t234 * t361 - t236 * t360 - t238 * t359) * t291) * t292, t228 * t329 + t229 * t324 + t230 * t319 + (-t246 * t308 - t247 * t304 - t248 * t300) * t337, -t228 * t330 - t229 * t325 - t230 * t320 + (t246 * t309 + t247 * t305 + t248 * t301) * t337, 1; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 3, 0, 0, 0, 0, 0, 0, 1;];
tau_reg  = t1;
